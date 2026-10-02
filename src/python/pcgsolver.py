from __future__ import annotations
import python.cuda_functions as cufunc

import faulthandler
import sys
import time

import numpy as np

import cupy as cp
from cupy import cublas
from cupy._core import _dtype
from cupy.cuda import device
from cupy_backends.cuda.libs import cublas as _cublas
from cupyx.scipy.sparse import _csr
from cupyx.scipy.sparse.linalg import _interface, spilu
from cupy_backends.cuda.libs import cusparse as _cusparse
from cupyx import cusparse
from cupyx.scipy.sparse import csr_matrix

from iop.cudaKernelWrapper import CudaDeviceConfig
from line_profiler import LineProfiler

checkEveryNIterations = 1


# ---------------------------------------------------------------------------
# Fused CUDA kernels for the multi-GPU PCG loop.
#
# All reductions are two-stage (per-block partials, then one block), so results
# are deterministic run to run. Scalars (rho, beta, alpha, rr) live on the GPU
# in a small array `sc`; partial dot products of all GPUs are summed in the same
# fixed order on every GPU, so every GPU computes bit-identical scalars.
# ---------------------------------------------------------------------------
_PCG_KERNELS = r'''
template<int NV>
__device__ __forceinline__ void block_reduce(double (&v)[NV], double* __restrict__ out)
{
    __shared__ double sh[NV][32];
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
    #pragma unroll
    for (int k = 0; k < NV; ++k) {
        double s = v[k];
        for (int o = 16; o > 0; o >>= 1) s += __shfl_down_sync(0xffffffffu, s, o);
        if (lane == 0) sh[k][warp] = s;
    }
    __syncthreads();
    if (warp == 0) {
        const int nw = blockDim.x >> 5;
        #pragma unroll
        for (int k = 0; k < NV; ++k) {
            double s = (lane < nw) ? sh[k][lane] : 0.0;
            for (int o = 16; o > 0; o >>= 1) s += __shfl_down_sync(0xffffffffu, s, o);
            if (lane == 0) out[k] = s;
        }
    }
}

#define GRID_STRIDE(i, n) \
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < (n); i += gridDim.x * blockDim.x)

// y = A x for a CSR row block, T threads per row (T = power of two, 2..32).
// Every thread of a warp reaches the shuffles (rows past the end contribute 0),
// so full-warp masks are valid. The reduction order is fixed -> deterministic.
template<int T>
__device__ __forceinline__ void csr_spmv_impl(const int* __restrict__ rowptr,
                                              const int* __restrict__ col,
                                              const double* __restrict__ val,
                                              const double* __restrict__ x,
                                              double* __restrict__ y, int nrows)
{
    const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    const int row  = tid / T;
    const int lane = threadIdx.x & (T - 1);
    double s = 0.0;
    if (row < nrows) {
        const int end = rowptr[row + 1];
        for (int k = rowptr[row] + lane; k < end; k += T)
            s += val[k] * x[col[k]];
    }
    #pragma unroll
    for (int o = T / 2; o > 0; o >>= 1)
        s += __shfl_down_sync(0xffffffffu, s, o, T);
    if (row < nrows && lane == 0) y[row] = s;
}

// Slots in the per-GPU scalar array sc (keep in sync with SC_* in Python)
enum { RHO = 0, BETA = 1, ALPHA = 2, RR = 3, ATOL2 = 4, CONV = 5, ITERS = 6 };

extern "C" {

#define CSR_SPMV(T) \
__global__ void csr_spmv_##T(const int* __restrict__ rowptr, const int* __restrict__ col, \
                             const double* __restrict__ val, const double* __restrict__ x, \
                             double* __restrict__ y, int nrows) \
{ csr_spmv_impl<T>(rowptr, col, val, x, y, nrows); }
CSR_SPMV(2)
CSR_SPMV(4)
CSR_SPMV(8)
CSR_SPMV(16)
CSR_SPMV(32)

// p = z + beta * p
__global__ void update_p(const double* __restrict__ z, double* __restrict__ p,
                         const double* __restrict__ sc, int n)
{
    const double beta = sc[BETA];
    GRID_STRIDE(i, n) p[i] = z[i] + beta * p[i];
}

// per-block partials of a·b
__global__ void dot_blocks(const double* __restrict__ a, const double* __restrict__ b,
                           double* __restrict__ blk, int n)
{
    double v[1] = {0.0};
    GRID_STRIDE(i, n) v[0] += a[i] * b[i];
    block_reduce<1>(v, blk + blockIdx.x);
}

// x += alpha p;  r -= alpha q;  z = invD r;  per-block partials of (r·r, r·z)
__global__ void update_xrz(double* __restrict__ x, double* __restrict__ r, double* __restrict__ z,
                           const double* __restrict__ p, const double* __restrict__ q,
                           const double* __restrict__ invD, const double* __restrict__ sc,
                           double* __restrict__ blk, int n)
{
    const double alpha = sc[ALPHA];
    double v[2] = {0.0, 0.0};
    GRID_STRIDE(i, n) {
        x[i] += alpha * p[i];
        const double ri = r[i] - alpha * q[i];
        const double zi = invD[i] * ri;
        r[i] = ri;  z[i] = zi;
        v[0] += ri * ri;  v[1] += ri * zi;
    }
    block_reduce<2>(v, blk + 2 * blockIdx.x);
}

// single block: sum nblk block partials into out[0]
__global__ void finish1(const double* __restrict__ blk, int nblk, double* __restrict__ out)
{
    double v[1] = {0.0};
    for (int b = threadIdx.x; b < nblk; b += blockDim.x) v[0] += blk[b];
    block_reduce<1>(v, out);
}

// single block: sum nblk pairs of block partials into out[0], out[1]
__global__ void finish2(const double* __restrict__ blk, int nblk, double* __restrict__ out)
{
    double v[2] = {0.0, 0.0};
    for (int b = threadIdx.x; b < nblk; b += blockDim.x) { v[0] += blk[2*b]; v[1] += blk[2*b+1]; }
    block_reduce<2>(v, out);
}

// alpha = rho / sum_j pq_j   (fixed summation order on every GPU)
// Once converged (CONV set), alpha = 0 so any remaining iterations of a graph chunk
// leave x and r untouched.
__global__ void scalar_alpha(const double* __restrict__ partPQ, int P, double* __restrict__ sc)
{
    if (sc[CONV] != 0.0) { sc[ALPHA] = 0.0; return; }
    double pq = 0.0;
    for (int j = 0; j < P; ++j) pq += partPQ[j];
    sc[ALPHA] = sc[RHO] / pq;
}

// beta = rz_new / rho_old (0 on the first pass);  rho = rz_new;  rr = sum_j rr_j.
// Counts iterations and sets CONV when rr <= atol^2. Once converged, beta = 0 and
// nothing else changes (frozen).
__global__ void scalar_beta(const double* __restrict__ partRZ, int P, double* __restrict__ sc, int first)
{
    if (sc[CONV] != 0.0) { sc[BETA] = 0.0; return; }
    double rr = 0.0, rz = 0.0;
    for (int j = 0; j < P; ++j) { rr += partRZ[2*j]; rz += partRZ[2*j+1]; }
    if (!isfinite(rr) || !isfinite(rz)) {          // NaN/Inf: flag error (-1) and freeze
        sc[BETA] = 0.0;  sc[RR] = rr;  sc[CONV] = -1.0;
        return;
    }
    sc[BETA] = (first || sc[RHO] == 0.0) ? 0.0 : rz / sc[RHO];
    sc[RHO]  = rz;
    sc[RR]   = rr;
    if (!first) sc[ITERS] += 1.0;
    if (rr <= sc[ATOL2]) sc[CONV] = 1.0;
}

}
'''




class PCGSolver:
    def __init__(self, N, NNZ):
        """Pre-allocate everything for maximum performance."""


        # # Create independent profiler for solver
        # self.solver_profiler = LineProfiler()
        # self.solveprof = self.solver_profiler(self.solve2)
        # self.solver_profiler.add_function(self.solve) 
        self.n = N
        
        self.data     = cp.zeros(NNZ, dtype=cp.float64)
        self.b        = cp.zeros(N,  dtype=cp.float64)


        self.indices= cp.zeros(NNZ,    dtype=cp.int32)
        self.indptr = cp.zeros(N + 1, dtype=cp.int32)
        self.indptr[N] = NNZ

        self.A   = csr_matrix((self.data  , self.indices, self.indptr), shape=(N, N))

        
        # Pre-allocate all work vectors (no allocations in solve loop)
        self.diag = cp.empty(self.n, dtype=cp.float64)
        self.r    = cp.empty(self.n, dtype=cp.float64)
        self.z    = cp.empty(self.n, dtype=cp.float64) 
        self.p    = cp.empty(self.n, dtype=cp.float64)
        self.q    = cp.empty(self.n, dtype=cp.float64)
        self.x    = cp.empty(self.n, dtype=cp.float64)
        self.x0   = cp.empty(self.n, dtype=cp.float64)
        self.x1   = cp.empty(self.n, dtype=cp.float64)
        self.x2   = cp.empty(self.n, dtype=cp.float64)

        # Set up scalars
        self.one  = np.array(1.0, np.float64)
        self.zero = np.array(0.0, np.float64)        
        
        self.cuda_dtype = _dtype.to_cuda_dtype(self.zero.dtype)

        # Pre-allocate cuSPARSE descriptors
        self.handle = device.get_cusparse_handle()
        self.desc_A = cusparse.SpMatDescriptor.create(self.A)
        self.desc_p = cusparse.DnVecDescriptor.create(self.p)
        self.desc_q = cusparse.DnVecDescriptor.create(self.q)
        self.desc_r = cusparse.DnVecDescriptor.create(self.r)
        self.desc_x = cusparse.DnVecDescriptor.create(self.x)
        
        # Pre-allocate buffer
        alg = _cusparse.CUSPARSE_MV_ALG_DEFAULT
        self.buff_size = _cusparse.spMV_bufferSize(
            self.handle, _cusparse.CUSPARSE_OPERATION_NON_TRANSPOSE,
            self.one.ctypes.data, self.desc_A.desc, self.desc_p.desc,
            self.zero.ctypes.data, self.desc_q.desc, self.cuda_dtype, alg)
        self.buff = cp.empty(self.buff_size, cp.int8)


    def mulA(self, p, q):
        """Compute q = A @ p via cuSPARSE SpMV."""
        # Reuse preallocated descriptors when callers pass self.p / self.q
        # (the CG hot path). Otherwise bind descriptors to the given vectors.
        desc_p = self.desc_p if p is self.p else cusparse.DnVecDescriptor.create(p)
        desc_q = self.desc_q if q is self.q else cusparse.DnVecDescriptor.create(q)
        _cusparse.spMV(
            self.handle, _cusparse.CUSPARSE_OPERATION_NON_TRANSPOSE,
            self.one.ctypes.data, self.desc_A.desc, desc_p.desc,
            self.zero.ctypes.data, desc_q.desc, self.cuda_dtype,
            _cusparse.CUSPARSE_MV_ALG_DEFAULT, self.buff.data.ptr)


    def solve(self, rtol=1e-5, maxiter = None, ensureSymmetric = False):
        """Solve with pre-allocated vectors - zero allocations in hot loop."""

        
        cufunc.setupCGSolver(self.data, self.b, self.indices, self.indptr)
        
        self.A = csr_matrix((self.data  , self.indices, self.indptr), shape=(self.n, self.n))
        self.A.sort_indices()

        if ensureSymmetric:
            # Make the matrix symmetric (removes small off-diagonal elements)
            self.A.data[:]   = 0.5*(self.A   + self.A.T).data  


        # Recreate descriptors after matrix data is populated
        self.desc_A = cusparse.SpMatDescriptor.create(self.A)

        # Uncomment this to use ILU preconditioner
        # self.ilu = spilu(self.A.astype(cp.float32), drop_tol=1e-6)       
        # store the inverse of the diagonal of the matrix
        self.diag[:] = 1.0/self.A.diagonal()   # XXXXX Can it be made faster?

        resid = 1e100

        if maxiter is None:
            maxiter = self.n * 10

        old_iters = 0

 
        # Initial residual: r = b - A*x
        self.r[:] = self.b - self.A @ self.x 
        
        # Tolerances are based on the first residual
        abs_r  = cp.linalg.norm(self.r)
        atol   = rtol * abs_r



        # Initialize solution using parabolic extrapolation from previous solutions
        self.x[:] = 3.0*self.x0 - 3.0*self.x1 + self.x2
        self.r[:] = self.b - self.A @ self.x
        self.p[:] = 0.0



        # CG iterations with pre-allocated vectors
        rho = 0.0
        for iters in range(old_iters, maxiter):
            # Apply the preconditioner
            # z = diag * r (in-place)
            # Comment this and uncomment the ILU preconditioner to use it
            cp.multiply(self.r, self.diag, out=self.z[:])
            # self.z[:] = self.ilu.solve(self.r)  # ILU preconditioner

            rho_new = cublas.dotc(self.r, self.z)
            beta = rho_new / rho if rho != 0.0 else 0.0
            
            # p = z + beta * p (in-place)
            self.p[:] *= beta
            cp.add(self.p, self.z, out=self.p[:])   # self.p[:] += self.z
            

            
            # q = A * p 
            # self.q[:] = self.A @ self.p
            self.mulA(self.p, self.q)

          
           
            alpha = rho_new / cublas.dotc(self.p, self.q)
            
            # x += alpha * p, r -= alpha * q (in-place)
            cublas.axpy( alpha, self.p, self.x)
            cublas.axpy(-alpha, self.q, self.r)
            
            rho = rho_new
            
            # Check convergence
            if iters % checkEveryNIterations == 0:  # Avoids checking every iteration because it's expensive
                # resid = cublas.nrm2(self.r)
                resid = cublas.dotc(self.r, self.r) # nrm2 is more stable but more expensive. At this point there is no risk of overflow, so we use dotc.

                # print(f"iters2: {iters}, residual: {resid}, residual2: {cp.linalg.norm(self.b - self.A @ self.x)**2}")

                if resid <= atol*atol:
                    break


        if iters >= maxiter - 1 and not (resid <= atol*atol):
            raise ValueError("PCG solver did not converge")


        cufunc.solverCopyResult(self.x)

        # Roll the previous solutions (for parabolic extrapolation)
        self.x2[:] = self.x1
        self.x1[:] = self.x0
        self.x0[:] = self.x
        
        return iters



    def csr_data(self, m=0):    return self.data
    def csr_indices(self, m=0): return self.indices
    def csr_indptr(self, m=0):  return self.indptr


    def print_stats(self):
        self.solver_profiler.print_stats()

class PCGSolverMultiGPU:
    TPB = 256   # threads per block for the fused kernels

    # Slots in the device-side scalar array sc (must match the enum in _PCG_KERNELS)
    SC_RHO, SC_BETA, SC_ALPHA, SC_RR, SC_ATOL2, SC_CONV, SC_ITERS = range(7)
    SC_SIZE = 8

    useGraph = True   # capture the CG loop in a CUDA Graph (falls back to eager launches on failure)
    waitTimeout = 60.0       # s: max wait for one chunk of GPU work before raising TimeoutError
    watchdogSeconds = 120.0  # s: if one solve takes longer, dump all thread stacks to stderr (None = off)

    def __init__(self, configCuda: CudaDeviceConfig, N: list[int], NNZ: list[int]):
        """Pre-allocate row-block matrices A_m and per-partition vectors.

        Naming:
          m            : partition / device index
          A_m[m]       : CSR row-block on device m, shape (N_m[m], N_total), nnz NNZ_m[m]
          *_m[m]       : vectors owned by partition m (length N_m[m])
          p[m], x[m]   : full packed vectors on device m (length N_total)
          p_m[m]       : owned segment of p, a VIEW into p[m] (no self-copy when gathering)
          x_m[m]       : owned segment of x (separate array)
          q_m[m]       : local matvec result  q_m = A_m @ p  (or A_m @ x)

        Column packing for A_m / p / x is [partition 0 | partition 1 | ...] via nodeOffset.
        """

        numParts = len(N)
        if len(NNZ) != numParts:
            raise ValueError(f"NNZ must have length numParts={numParts}, got {len(NNZ)}")

        if len(configCuda.devices) != numParts:
            raise ValueError(f"Number of devices ({len(configCuda.devices)}) must match numParts ({numParts})")

        self.configCuda = configCuda
        self.numParts = numParts
        self.N_m = list(N)            # N_m[m] = rows owned by partition m
        self.NNZ_m = list(NNZ)        # NNZ_m[m] = nonzeros in row-block A_m
        self.N_total = int(sum(self.N_m))

        # Packed column offsets for p / CSR column indices
        self.nodeOffset = [0] * (numParts + 1)
        for m in range(numParts):
            self.nodeOffset[m + 1] = self.nodeOffset[m] + self.N_m[m]

        # Row-block CSR storage (one A_m per device)
        self.A_m       = [None] * numParts
        self.data_m    = [None] * numParts
        self.indices_m = [None] * numParts
        self.indptr_m  = [None] * numParts

        # Local (partition-owned) vectors, length N_m[m]
        self.b_m    = [None] * numParts
        self.invD_m = [None] * numParts
        self.r_m    = [None] * numParts
        self.z_m    = [None] * numParts
        self.p_m    = [None] * numParts   # owned segment of CG search direction (view into p[m])
        self.q_m    = [None] * numParts   # q_m = A_m @ p
        self.x_m    = [None] * numParts
        self.x0_m   = [None] * numParts
        self.x1_m   = [None] * numParts
        self.x2_m   = [None] * numParts

        # Full packed vectors on each device (SpMV inputs)
        self.p = [None] * numParts
        self.x = [None] * numParts

        # cuSPARSE state (one SpMV of A_m per device)
        self.handle   = [None] * numParts
        self.desc_A_m = [None] * numParts
        self.desc_p   = [None] * numParts
        self.desc_x   = [None] * numParts
        self.desc_q_m = [None] * numParts
        self.buff_m   = [None] * numParts

        self.one  = np.array(1.0, np.float64)
        self.zero = np.array(0.0, np.float64)
        self.cuda_dtype = _dtype.to_cuda_dtype(self.zero.dtype)
        alg = _cusparse.CUSPARSE_MV_ALG_DEFAULT

        for m, device_id in enumerate(configCuda.devices):
            with cp.cuda.Device(device_id):
                Nm = self.N_m[m]
                NNZm = self.NNZ_m[m]

                self.data_m[m] = cp.zeros(NNZm, dtype=cp.float64)
                self.indices_m[m] = cp.zeros(NNZm, dtype=cp.int32)
                self.indptr_m[m] = cp.zeros(Nm + 1, dtype=cp.int32)
                self.indptr_m[m][Nm] = NNZm

                # A_m: Nm rows (local), N_total columns (packed partitions)
                self.A_m[m] = csr_matrix(
                    (self.data_m[m], self.indices_m[m], self.indptr_m[m]),
                    shape=(Nm, self.N_total),
                )

                self.b_m[m]    = cp.zeros(Nm, dtype=cp.float64)
                self.invD_m[m] = cp.zeros(Nm, dtype=cp.float64)
                self.r_m[m]    = cp.zeros(Nm, dtype=cp.float64)
                self.z_m[m]    = cp.zeros(Nm, dtype=cp.float64)
                self.q_m[m]    = cp.zeros(Nm, dtype=cp.float64)
                self.x_m[m]    = cp.zeros(Nm, dtype=cp.float64)
                self.x0_m[m]   = cp.zeros(Nm, dtype=cp.float64)
                self.x1_m[m]   = cp.zeros(Nm, dtype=cp.float64)
                self.x2_m[m]   = cp.zeros(Nm, dtype=cp.float64)

                self.p[m] = cp.zeros(self.N_total, dtype=cp.float64)
                self.x[m] = cp.zeros(self.N_total, dtype=cp.float64)

                # Owned segment of p is a view into the full p: no self-copy in the gather
                self.p_m[m] = self.p[m][self.nodeOffset[m]:self.nodeOffset[m + 1]]

                self.handle[m] = device.get_cusparse_handle()
                self.desc_A_m[m] = cusparse.SpMatDescriptor.create(self.A_m[m])
                self.desc_p[m]   = cusparse.DnVecDescriptor.create(self.p[m])
                self.desc_x[m]   = cusparse.DnVecDescriptor.create(self.x[m])
                self.desc_q_m[m] = cusparse.DnVecDescriptor.create(self.q_m[m])

                buff_size = _cusparse.spMV_bufferSize(
                    self.handle[m],
                    _cusparse.CUSPARSE_OPERATION_NON_TRANSPOSE,
                    self.one.ctypes.data,
                    self.desc_A_m[m].desc,
                    self.desc_p[m].desc,
                    self.zero.ctypes.data,
                    self.desc_q_m[m].desc,
                    self.cuda_dtype,
                    alg,
                )
                self.buff_m[m] = cp.empty(buff_size, dtype=cp.int8)

        self._setupFused()


    # ------------------------------------------------------------------
    # Fused-kernel / stream / event setup
    # ------------------------------------------------------------------

    def _setupFused(self):
        """Solver streams, events, compiled kernels and device-side scalar buffers."""
        P = self.numParts
        names = ('update_p', 'dot_blocks', 'update_xrz', 'finish1', 'finish2',
                 'scalar_alpha', 'scalar_beta',
                 'csr_spmv_2', 'csr_spmv_4', 'csr_spmv_8', 'csr_spmv_16', 'csr_spmv_32')
        self.stream = [None] * P
        self.evP    = [None] * P     # "p_m is ready" (before the gather)
        self.evRed  = [None] * P     # "my partial dot is ready" (before an exchange)
        self.kern   = [None] * P
        self.nblk   = [0] * P
        self.nblk32 = [None] * P
        self.blk    = [None] * P     # per-block partial sums
        self.sc     = [None] * P     # [rho, beta, alpha, rr, atol2, conv, iters] on every GPU
        self.partPQ = [None] * P     # p·q partial of every GPU (P values)
        self.partRZ = [None] * P     # (r·r, r·z) partials of every GPU (2P values)
        self.n32    = [np.int32(n) for n in self.N_m]

        # CUDA Graph state
        self.evJoin  = [None] * P    # graph capture: other GPUs join back into stream[0]
        self.evPre   = [None] * P    # graph launch: stream[0] waits for prior work on other GPUs
        self.evFork  = None          # graph capture: other GPUs join the capture
        self._graph  = None
        self._graphK = None

        for m, dev in enumerate(self.configCuda.devices):
            with cp.cuda.Device(dev):
                self.stream[m] = cp.cuda.Stream(non_blocking=True)
                self.evP[m]    = cp.cuda.Event(disable_timing=True)
                self.evRed[m]  = cp.cuda.Event(disable_timing=True)
                self.evJoin[m] = cp.cuda.Event(disable_timing=True)
                self.evPre[m]  = cp.cuda.Event(disable_timing=True)
                if m == 0:
                    self.evFork = cp.cuda.Event(disable_timing=True)
                    self.evDone = cp.cuda.Event(disable_timing=True)   # chunk completion (polled)
                mod = cp.RawModule(code=_PCG_KERNELS, options=('-std=c++17',))
                self.kern[m] = {nm: mod.get_function(nm) for nm in names}
                sms = cp.cuda.Device(dev).attributes['MultiProcessorCount']
                self.nblk[m]   = min(4 * sms, 1024)
                self.nblk32[m] = np.int32(self.nblk[m])
                self.blk[m]    = cp.zeros(2 * self.nblk[m], cp.float64)
                self.sc[m]     = cp.zeros(self.SC_SIZE, cp.float64)
                self.partPQ[m] = cp.zeros(P, cp.float64)
                self.partRZ[m] = cp.zeros(2 * P, cp.float64)

        # Custom CSR SpMV for the solver loop (a plain kernel: capturable and replayable).
        # Threads per row: smallest power of two >= average nonzeros per row (2..32).
        self.spmvLanes = [0] * P
        self.spmvGrid  = [0] * P
        for m in range(P):
            avg = self.NNZ_m[m] / max(self.N_m[m], 1)
            T = 2
            while T < avg and T < 32:
                T *= 2
            self.spmvLanes[m] = T
            self.spmvGrid[m]  = max(1, -(-self.N_m[m] * T // self.TPB))   # ceil(N_m * T / TPB)
        self._useCustomSpMV = True     # set False if the self-check below fails
        self._spmvChecked   = False


    def _spmvCustom(self, m, x_full):
        """q_m = A_m @ x_full with the custom CSR kernel, on the current stream of device m."""
        self.kern[m][f'csr_spmv_{self.spmvLanes[m]}'](
            (self.spmvGrid[m],), (self.TPB,),
            (self.indptr_m[m], self.indices_m[m], self.data_m[m], x_full, self.q_m[m], self.n32[m]))


    def _checkCustomSpMV(self):
        """Compare the custom SpMV with cuSPARSE on the real matrix (first solve only).

        Uses the current x (already gathered by computeResiduals). The two differ only
        in summation order, so a tolerance is used. On a mismatch or error, falls back
        to cuSPARSE without CUDA Graphs (the path that is known to work).
        """
        self._spmvChecked = True
        try:
            ok = True
            for m, dev in enumerate(self.configCuda.devices):
                with cp.cuda.Device(dev):
                    self._spmvLocal(m, self.desc_x)            # cuSPARSE reference
                    ref = self.q_m[m].copy()
                    self._spmvCustom(m, self.x[m])             # same (default) stream
                    scale = float(cp.max(cp.abs(ref))) if ref.size else 0.0
                    same = bool(cp.allclose(self.q_m[m], ref, rtol=1e-10, atol=1e-10 * scale))
                    if not same:
                        err = float(cp.max(cp.abs(self.q_m[m] - ref)))
                        print(f"PCGSolverMultiGPU: custom SpMV mismatch on partition {m} "
                              f"(max abs diff {err:.3e}, max |ref| {scale:.3e}).", flush=True)
                    ok = ok and same
        except Exception as e:
            print(f"PCGSolverMultiGPU: custom SpMV failed ({e!r}).", flush=True)
            ok = False
        if ok:
            print(f"PCGSolverMultiGPU: custom SpMV verified against cuSPARSE "
                  f"(threads per row: {self.spmvLanes}).", flush=True)
        else:
            print("PCGSolverMultiGPU: using cuSPARSE without CUDA Graphs.", flush=True)
            self._useCustomSpMV = False
            self.useGraph = False
            self.resetGraph()
        self._syncAll()


    def _syncAll(self):
        """Wait for all work on all devices (all streams)."""
        for dev in self.configCuda.devices:
            cp.cuda.Device(dev).synchronize()


    # ------------------------------------------------------------------
    # Setup-phase helpers (default stream; host syncs allowed)
    # ------------------------------------------------------------------

    def gatherVector(self, full, local):
        """Gather full[dest] on every GPU from local[src] pieces (packed by partition)."""
        # Make sure all source partitions finished writing before peer copies.
        for devId in self.configCuda.devices:
            with cp.cuda.Device(devId):
                cp.cuda.Stream.null.synchronize()

        for dest_m, dest_id in enumerate(self.configCuda.devices):
            with cp.cuda.Device(dest_id):
                for src_m, src_id in enumerate(self.configCuda.devices):
                    lo = self.nodeOffset[src_m]
                    hi = self.nodeOffset[src_m + 1]
                    dst = full[dest_m][lo:hi]
                    src = local[src_m]
                    if dest_id == src_id:
                        if dst.data.ptr != src.data.ptr:   # p_m is a view of p: nothing to copy
                            dst[:] = src
                    else:
                        cp.cuda.runtime.memcpyPeerAsync(
                            dst.data.ptr,
                            dest_id,
                            src.data.ptr,
                            src_id,
                            dst.nbytes,
                            cp.cuda.get_current_stream().ptr,
                        )


        # Synchronize after all peer copies are enqueued
        for devId in self.configCuda.devices:
            with cp.cuda.Device(devId):
                cp.cuda.Stream.null.synchronize()


    def gatherP(self):
        """Gather full p from local p_m pieces (setup/debug path; the solver loop uses _gatherPAsync)."""
        self.gatherVector(self.p, self.p_m)


    def gather_x(self):
        """Gather full x from local x_m pieces."""
        self.gatherVector(self.x, self.x_m)


    def spmv(self, desc_in):
        """On each device m: q_m[m] = A_m[m] @ in[m], on the current stream."""
        for m, device_id in enumerate(self.configCuda.devices):
            with cp.cuda.Device(device_id):
                self._spmvLocal(m, desc_in)


    def _spmvLocal(self, m, desc_in):
        """q_m = A_m @ in on the *current* stream of the current device.

        CuPy's spMV wrapper binds the handle to the current stream itself.
        """
        _cusparse.spMV(
            self.handle[m],
            _cusparse.CUSPARSE_OPERATION_NON_TRANSPOSE,
            self.one.ctypes.data,
            self.desc_A_m[m].desc,
            desc_in[m].desc,
            self.zero.ctypes.data,
            self.desc_q_m[m].desc,
            self.cuda_dtype,
            _cusparse.CUSPARSE_MV_ALG_DEFAULT,
            self.buff_m[m].data.ptr,
        )


    def mulA(self):
        """q_m = A_m @ p. Call gatherP() first."""
        self.spmv(self.desc_p)


    def extrapolateGuess(self):
        """Parabolic initial guess from previous solutions; zero local search direction."""
        for m, device_id in enumerate(self.configCuda.devices):
            with cp.cuda.Device(device_id):
                self.x_m[m][:] = 3.0 * self.x0_m[m] - 3.0 * self.x1_m[m] + self.x2_m[m]
                self.p_m[m][:] = 0.0


    def computeResiduals(self):
        """r_m = b_m - A_m @ x  (gathers full x first)."""
        self.gather_x()
        self.spmv(self.desc_x)
        for m, device_id in enumerate(self.configCuda.devices):
            with cp.cuda.Device(device_id):
                cp.subtract(self.b_m[m], self.q_m[m], out=self.r_m[m])


    def dotSums(self, pairs):
        """For each (a_m, b_m) in pairs, return sum_m a_m[m] · b_m[m] (host float64).

        Setup-phase helper: all local dots are enqueued on all GPUs first, then read back.
        """
        k = len(pairs)
        devs = self.configCuda.devices
        partial = [[None] * k for _ in devs]
        for m, dev in enumerate(devs):
            with cp.cuda.Device(dev):
                for j, (a_m, b_m) in enumerate(pairs):
                    partial[m][j] = cp.dot(a_m[m], b_m[m])
        totals = np.zeros(k, dtype=np.float64)
        for m, dev in enumerate(devs):
            with cp.cuda.Device(dev):
                for j in range(k):
                    totals[j] += float(partial[m][j])
        return totals


    def dotSum(self, a_m, b_m):
        """sum_m a_m · b_m (host float64)."""
        return self.dotSums(((a_m, b_m),))[0]


    def norm2Sum(self, a_m):
        """||a|| = sqrt(sum_m a_m · a_m)."""
        return np.sqrt(self.dotSum(a_m, a_m))


    def computeAtol(self, rtol):
        """atol = rtol * ||b|| (stable across timesteps; matches checkSolution)."""
        return rtol * self.norm2Sum(self.b_m)


    def computeInvDiag(self):
        """invD_m = 1.0 / A_m.diagonal()."""
        for m, device_id in enumerate(self.configCuda.devices):
            with cp.cuda.Device(device_id):
                k = self.nodeOffset[m]
                self.invD_m[m][:] = 1.0 / self.A_m[m].diagonal(k)


    def rollHistory(self):
        """x2_m <- x1_m <- x0_m <- x_m (parabolic extrapolation history)."""
        for m, device_id in enumerate(self.configCuda.devices):
            with cp.cuda.Device(device_id):
                self.x2_m[m][:] = self.x1_m[m]
                self.x1_m[m][:] = self.x0_m[m]
                self.x0_m[m][:] = self.x_m[m]


    def gatherBFromMesh(self):
        """b_m <- qel on each partition."""
        for m, device_id in enumerate(self.configCuda.devices):
            with cp.cuda.Device(device_id):
                cufunc.gatherBFromMesh(self.b_m[m])


    def scatterXToMesh(self):
        """mX <- x_m on each partition."""
        for m, device_id in enumerate(self.configCuda.devices):
            with cp.cuda.Device(device_id):
                cufunc.scatterXToMesh(self.x_m[m])


    # ------------------------------------------------------------------
    # Solver-loop helpers (solver streams; no host synchronization)
    #
    # Ordering between GPUs uses events only. Buffers read by another GPU are
    # never overwritten early:
    #   - GPU A rewrites its p segment in the next update_p, after its stream
    #     waited for GPU B's p·q partial, which B computed after pulling A's p.
    #   - A p·q partial slot is rewritten one iteration later, after its owner
    #     waited for the other GPU's (r·r, r·z) partial, computed after the pull.
    #   - Same for the (r·r, r·z) slot, with the exchanges swapped.
    # ------------------------------------------------------------------

    def _gatherPAsync(self):
        """Full p on every GPU; ordered with events, no host sync."""
        devs = self.configCuda.devices
        for m, dev in enumerate(devs):
            with cp.cuda.Device(dev):
                self.evP[m].record(self.stream[m])
        for dst, ddev in enumerate(devs):
            with cp.cuda.Device(ddev):
                s = self.stream[dst]
                for src, sdev in enumerate(devs):
                    if src == dst:
                        continue
                    s.wait_event(self.evP[src])
                    cp.cuda.runtime.memcpyPeerAsync(
                        self.p[dst].data.ptr + self.nodeOffset[src] * 8, ddev,
                        self.p_m[src].data.ptr, sdev,
                        self.N_m[src] * 8, s.ptr)


    def _exchange(self, part, nv):
        """Every GPU pulls the nv partial values of every other GPU into its own part array."""
        devs = self.configCuda.devices
        nbytes = nv * 8
        for m, dev in enumerate(devs):
            with cp.cuda.Device(dev):
                self.evRed[m].record(self.stream[m])
        for dst, ddev in enumerate(devs):
            with cp.cuda.Device(ddev):
                s = self.stream[dst]
                for src, sdev in enumerate(devs):
                    if src == dst:
                        continue
                    s.wait_event(self.evRed[src])
                    cp.cuda.runtime.memcpyPeerAsync(
                        part[dst].data.ptr + src * nbytes, ddev,
                        part[src].data.ptr + src * nbytes, sdev,
                        nbytes, s.ptr)


    def _updateXRZ(self, computeAlpha):
        """[alpha = rho / p·q]; x += alpha p; r -= alpha q; z = invD r; local (r·r, r·z)."""
        P32 = np.int32(self.numParts)
        for m, dev in enumerate(self.configCuda.devices):
            with cp.cuda.Device(dev), self.stream[m]:
                k = self.kern[m]
                if computeAlpha:
                    k['scalar_alpha']((1,), (1,), (self.partPQ[m], P32, self.sc[m]))
                k['update_xrz']((self.nblk[m],), (self.TPB,),
                    (self.x_m[m], self.r_m[m], self.z_m[m], self.p_m[m], self.q_m[m],
                     self.invD_m[m], self.sc[m], self.blk[m], self.n32[m]))
                k['finish2']((1,), (self.TPB,),
                    (self.blk[m], self.nblk32[m], self.partRZ[m][2 * m:2 * m + 2]))


    def _scalarBeta(self, first):
        """beta = rz_new / rho_old (0 if first); rho = rz_new; rr stored for the convergence check."""
        P32, f = np.int32(self.numParts), np.int32(first)
        for m, dev in enumerate(self.configCuda.devices):
            with cp.cuda.Device(dev), self.stream[m]:
                self.kern[m]['scalar_beta']((1,), (1,), (self.partRZ[m], P32, self.sc[m], f))


    def _iteration(self):
        """One PCG iteration on all GPUs (solver streams only; capturable in a CUDA Graph)."""
        devs = self.configCuda.devices
        for m, dev in enumerate(devs):                        # p = z + beta p
            with cp.cuda.Device(dev), self.stream[m]:
                self.kern[m]['update_p']((self.nblk[m],), (self.TPB,),
                    (self.z_m[m], self.p_m[m], self.sc[m], self.n32[m]))

        self._gatherPAsync()                                   # full p everywhere

        for m, dev in enumerate(devs):                        # q = A p; local p·q
            with cp.cuda.Device(dev), self.stream[m]:
                k = self.kern[m]
                if self._useCustomSpMV:
                    self._spmvCustom(m, self.p[m])
                else:
                    self._spmvLocal(m, self.desc_p)                # cuSPARSE fallback (not capturable)
                k['dot_blocks']((self.nblk[m],), (self.TPB,),
                    (self.p_m[m], self.q_m[m], self.blk[m], self.n32[m]))
                k['finish1']((1,), (self.TPB,),
                    (self.blk[m], self.nblk32[m], self.partPQ[m][m:m + 1]))

        self._exchange(self.partPQ, 1)
        self._updateXRZ(computeAlpha=True)                     # alpha, x, r, z, local (r·r, r·z)
        self._exchange(self.partRZ, 2)
        self._scalarBeta(first=0)                              # beta, rho, rr, iters, conv


    # ------------------------------------------------------------------
    # CUDA Graph: capture K iterations once, replay them for every chunk
    #
    # Capture starts on stream[0]; the other GPUs' streams join it by waiting
    # on evFork and are joined back into stream[0] at the end, so one graph
    # holds the work of all GPUs. The graph is reused across solves: every
    # buffer has a fixed address and all scalars live on the GPU. If buffers
    # are ever reallocated (not just overwritten), call resetGraph().
    # ------------------------------------------------------------------

    def resetGraph(self):
        """Discard the captured graph; it is recaptured on the next solve."""
        self._graph = None
        self._graphK = None


    def _captureGraph(self, K):
        devs = self.configCuda.devices
        P = self.numParts
        s0 = self.stream[0]
        with cp.cuda.Device(devs[0]):
            s0.begin_capture()
        try:
            with cp.cuda.Device(devs[0]):
                self.evFork.record(s0)
            for m in range(1, P):                              # fork: other GPUs join the capture
                with cp.cuda.Device(devs[m]):
                    self.stream[m].wait_event(self.evFork)
            for _ in range(K):
                self._iteration()
            for m in range(1, P):                              # join back into stream[0]
                with cp.cuda.Device(devs[m]):
                    self.evJoin[m].record(self.stream[m])
            with cp.cuda.Device(devs[0]):
                for m in range(1, P):
                    s0.wait_event(self.evJoin[m])
        except Exception:
            with cp.cuda.Device(devs[0]):
                try:
                    s0.end_capture()
                except Exception:
                    pass
            raise
        with cp.cuda.Device(devs[0]):
            self._graph = s0.end_capture()
        self._graphK = K


    def _launchGraph(self):
        """Replay K iterations. The launch is ordered after all prior work on every solver stream."""
        devs = self.configCuda.devices
        s0 = self.stream[0]
        for m in range(1, self.numParts):
            with cp.cuda.Device(devs[m]):
                self.evPre[m].record(self.stream[m])
        with cp.cuda.Device(devs[0]):
            for m in range(1, self.numParts):
                s0.wait_event(self.evPre[m])
            self._graph.launch(s0)


    def _cgLoop(self, maxiter, K):
        """Run PCG in chunks of K iterations until converged. Returns the iteration count.

        After convergence the GPUs freeze (alpha = beta = 0), so the rest of a chunk
        changes nothing and the count read from the GPU is exact.
        """
        devs = self.configCuda.devices
        C, I = self.SC_CONV, self.SC_ITERS
        done = 0
        iters = 0
        chunk = 0
        while done < maxiter:
            chunk += 1
            ranGraph = self.useGraph and self._graph is not None and self._graphK == K
            if ranGraph:
                self._launchGraph()
            else:
                for _ in range(K):                             # eager chunk (also the warm-up)
                    self._iteration()
                if self.useGraph:
                    try:
                        self._captureGraph(K)
                        print(f"PCGSolverMultiGPU: CUDA Graph captured ({K} iterations).", flush=True)
                    except Exception as e:
                        print(f"PCGSolverMultiGPU: CUDA Graph capture failed ({e!r}); "
                              f"continuing without graphs.")
                        self.useGraph = False
                        self.resetGraph()
            done += K

            # The only host read: converged flag and iteration count (stream[0]
            # holds the joined work of all GPUs). Wait with a time limit, so a
            # stuck chunk raises an error instead of hanging the job.
            what = (f"chunk {chunk}, " +
                    ("replayed from the CUDA Graph" if ranGraph
                     else "eager launches" + (" + graph capture" if self._graph is not None else "")))
            conv, it = self._readStatus(C, I, what)
            iters = int(it)
            if conv < 0.0:
                raise RuntimeError(f"PCGSolverMultiGPU: non-finite residual (NaN/Inf) "
                                   f"after {iters} iterations.")
            if conv != 0.0:
                return iters
        print(f"PCGSolverMultiGPU: WARNING: not converged after {iters} iterations "
              f"(maxiter={maxiter}).", flush=True)
        return iters


    def _readStatus(self, C, I, what):
        """Wait (bounded by waitTimeout) for stream[0] to finish, then read sc[C..I]."""
        devs = self.configCuda.devices
        with cp.cuda.Device(devs[0]):
            self.evDone.record(self.stream[0])
            t0 = time.perf_counter()
            while not self.evDone.done:              # cudaEventQuery: never blocks
                if time.perf_counter() - t0 > self.waitTimeout:
                    raise TimeoutError(f"PCGSolverMultiGPU: GPU work did not finish within "
                                       f"{self.waitTimeout:.0f} s ({what}).")
            vals = self.sc[0][C:I + 1].get(stream=self.stream[0])
        return float(vals[0]), float(vals[1])


    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def solve(self, rtol=1e-5, maxiter=None, ensureSymmetric=False, checkEvery=10):
        """Jacobi-preconditioned CG across all devices.

        Iterations run in chunks of `checkEvery` on per-device solver streams with
        device-side scalars, replayed from a CUDA Graph after the first chunk. The
        host reads the convergence flag once per chunk. The returned iteration
        count is exact (iterations after convergence are frozen no-ops).
        """
        if maxiter is None:
            maxiter = self.N_total * 10

        # Watchdog: if this solve hangs (even inside a C/CUDA call), print every
        # thread's stack to stderr after watchdogSeconds.
        if self.watchdogSeconds:
            faulthandler.dump_traceback_later(self.watchdogSeconds, repeat=True, file=sys.stderr)
        try:
            return self._solve(rtol, maxiter, checkEvery)
        finally:
            if self.watchdogSeconds:
                faulthandler.cancel_dump_traceback_later()


    def _solve(self, rtol, maxiter, checkEvery):
        # Per-solve setup on the default stream (host syncs are fine here)
        self.gatherBFromMesh()
        self.computeInvDiag()
        atol = self.computeAtol(rtol)
        self.extrapolateGuess()
        self.computeResiduals()                    # r_m = b_m - A_m @ x_guess (gathers x)
        if not self._spmvChecked:
            self._checkCustomSpMV()                # first solve only; uses the gathered x
        for m, dev in enumerate(self.configCuda.devices):
            with cp.cuda.Device(dev):
                self.sc[m].fill(0.0)               # alpha = 0, conv = 0, iters = 0
                self.sc[m][self.SC_ATOL2] = atol * atol
        self._syncAll()                            # hand over to the solver streams

        # Init: z = invD r, rho = r·z, rr = r·r  (update_xrz with alpha = 0)
        self._updateXRZ(computeAlpha=False)
        self._exchange(self.partRZ, 2)
        self._scalarBeta(first=1)                  # also sets conv if already below atol

        iters = self._cgLoop(maxiter, checkEvery)

        self._syncAll()                            # back to the default stream
        self.scatterXToMesh()
        self.rollHistory()
        return iters


    def checkSolution(self, label=""):
        """True residual after solve: ||b - A x|| / ||b|| (and absolute norms)."""
        self.gather_x()
        self.spmv(self.desc_x)  # q_m = A_m @ x
        rnorm2 = np.float64(0.0)
        bnorm2 = np.float64(0.0)
        for m, device_id in enumerate(self.configCuda.devices):
            with cp.cuda.Device(device_id):
                r = self.b_m[m] - self.q_m[m]
                rnorm2 += np.float64(cublas.dotc(r, r))
                bnorm2 += np.float64(cublas.dotc(self.b_m[m], self.b_m[m]))
        rnorm = np.sqrt(rnorm2)
        bnorm = np.sqrt(bnorm2)
        rel = rnorm / bnorm if bnorm > 0 else rnorm
        print(f"[PCG check{(' ' + label) if label else ''}] "
            f"||r||={rnorm:.6e}  ||b||={bnorm:.6e}  ||r||/||b||={rel:.6e}")
        return rel


    def csr_data(self, m):      return self.data_m[m]
    def csr_indices(self, m):   return self.indices_m[m]
    def csr_indptr(self, m):    return self.indptr_m[m]


    def print_stats(self):
        self.solver_profiler.print_stats()