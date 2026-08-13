from __future__ import annotations
import python.cuda_functions as cufunc

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
    def __init__(self, configCuda: CudaDeviceConfig, N: list[int], NNZ: list[int]):
        """Pre-allocate row-block matrices A_m and per-partition vectors.

        Naming:
          m            : partition / device index
          A_m[m]       : CSR row-block on device m, shape (N_m[m], N_total), nnz NNZ_m[m]
          *_m[m]       : vectors owned by partition m (length N_m[m])
          p[m], x[m]   : full packed vectors on device m (length N_total)
          p_m[m], x_m[m] : owned segments (length N_m[m])
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
        self.p_m    = [None] * numParts   # owned segment of CG search direction
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
                self.p_m[m]    = cp.zeros(Nm, dtype=cp.float64)
                self.q_m[m]    = cp.zeros(Nm, dtype=cp.float64)
                self.x_m[m]    = cp.zeros(Nm, dtype=cp.float64)
                self.x0_m[m]   = cp.zeros(Nm, dtype=cp.float64)
                self.x1_m[m]   = cp.zeros(Nm, dtype=cp.float64)
                self.x2_m[m]   = cp.zeros(Nm, dtype=cp.float64)

                self.p[m] = cp.zeros(self.N_total, dtype=cp.float64)
                self.x[m] = cp.zeros(self.N_total, dtype=cp.float64)

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


    def gatherVector(self, full, local):
        """Gather full[dest] on every GPU from local[src] pieces (packed by partition)."""
        # Make sure all source partitions finished writing before peer copies.
        for src_id in self.configCuda.devices:
            with cp.cuda.Device(src_id):
                cp.cuda.Stream.null.synchronize()

        for dest_m, dest_id in enumerate(self.configCuda.devices):
            with cp.cuda.Device(dest_id):
                for src_m, src_id in enumerate(self.configCuda.devices):
                    lo = self.nodeOffset[src_m]
                    hi = self.nodeOffset[src_m + 1]
                    dst = full[dest_m][lo:hi]
                    src = local[src_m]
                    if dest_id == src_id:
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
                cp.cuda.Stream.null.synchronize()


    def gatherP(self):
        """Gather full p from local p_m pieces."""
        self.gatherVector(self.p, self.p_m)


    def gather_x(self):
        """Gather full x from local x_m pieces."""
        self.gatherVector(self.x, self.x_m)


    def spmv(self, desc_in):
        """On each device m: q_m[m] = A_m[m] @ in[m]."""
        for m, device_id in enumerate(self.configCuda.devices):
            with cp.cuda.Device(device_id):
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
        """Parabolic initial guess from previous solutions; zero local search direction.

        Call after computeAtol() so atol is based on the previous x, not this guess.
        """
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


    def norm2Sum(self, a_m):
        """||a|| = sqrt(sum_m a_m · a_m)."""
        return np.sqrt(self.dotSum(a_m, a_m))


    def computeAtol(self, rtol):
        """atol = rtol * ||b|| (stable across timesteps; matches checkSolution)."""
        return rtol * self.norm2Sum(self.b_m)


    def computeInvDiag(self):
        """invD_m = 1.0 / A_m.diagonal() (gathers full A first)."""

        for m, device_id in enumerate(self.configCuda.devices):
            with cp.cuda.Device(device_id):
                k = self.nodeOffset[m]
                self.invD_m[m][:] = 1.0 / self.A_m[m].diagonal(k)

    def applyPreconditioner(self):
        """z_m = invD_m * r_m (Jacobi)."""
        for m, device_id in enumerate(self.configCuda.devices):
            with cp.cuda.Device(device_id):
                cp.multiply(self.r_m[m], self.invD_m[m], out=self.z_m[m])


    def dotSum(self, a_m, b_m):
        """Sum of local dots across partitions: sum_m a_m · b_m (host float64)."""
        total = np.float64(0.0)
        for m, device_id in enumerate(self.configCuda.devices):
            with cp.cuda.Device(device_id):
                total += np.float64(cublas.dotc(a_m[m], b_m[m]))
        return total


    def dotRZ(self):
        """rho = sum_m r_m · z_m."""
        return self.dotSum(self.r_m, self.z_m)


    def dotPQ(self):
        """sum_m p_m · q_m."""
        return self.dotSum(self.p_m, self.q_m)


    def updateP(self, rhoNew, rho):
        """p_m = z_m + beta * p_m, with beta = rhoNew/rho (0 if rho == 0)."""
        beta = rhoNew / rho if rho != 0.0 else np.float64(0.0)
        for m, device_id in enumerate(self.configCuda.devices):
            with cp.cuda.Device(device_id):
                self.p_m[m][:] *= beta
                cp.add(self.p_m[m], self.z_m[m], out=self.p_m[m])


    def axpyUpdate(self, alpha):
        """x_m += alpha * p_m, r_m -= alpha * q_m."""
        for m, device_id in enumerate(self.configCuda.devices):
            with cp.cuda.Device(device_id):
                cublas.axpy( alpha, self.p_m[m], self.x_m[m])
                cublas.axpy(-alpha, self.q_m[m], self.r_m[m])


    def rollHistory(self):
        """x2_m <- x1_m <- x0_m <- x_m (parabolic extrapolation history)."""
        for m, device_id in enumerate(self.configCuda.devices):
            with cp.cuda.Device(device_id):
                self.x2_m[m][:] = self.x1_m[m]
                self.x1_m[m][:] = self.x0_m[m]
                self.x0_m[m][:] = self.x_m[m]


    def residualBelow(self, atol):
        """True if sum_m r_m · r_m <= atol^2."""
        return self.dotSum(self.r_m, self.r_m) <= atol * atol


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


    def solve(self, rtol=1e-5, maxiter=None, ensureSymmetric=False):

        if maxiter is None:
            maxiter = self.N_total * 10

        self.gatherBFromMesh()
        self.computeInvDiag()
        atol = self.computeAtol(rtol)

        print("N_m", self.N_m,
            "N_total", self.N_total,
            "sum(N_m)", sum(self.N_m),
            "indices max", [int(self.indices_m[m].max()) for m in range(self.numParts)])
        for m, device_id in enumerate(self.configCuda.devices):
            with cp.cuda.Device(device_id):
                d = self.A_m[m].diagonal(self.nodeOffset[m])
                print(m, "||b||", float(cp.linalg.norm(self.b_m[m])),
                    "nnz(data)", int(cp.count_nonzero(self.data_m[m])),
                    "min|diag|", float(cp.min(cp.abs(d))),
                    "num zero diag", int(cp.sum(d == 0)),
                    "indices min/max", int(self.indices_m[m].min()), int(self.indices_m[m].max()))

        self.extrapolateGuess()
        self.computeResiduals()          # r_m = b_m - A_m @ x_guess
        r0 = self.norm2Sum(self.r_m)
        print(f"PCG: ||r0||={r0:.6e}  atol={atol:.6e}  ||b||={self.norm2Sum(self.b_m):.6e}")

        rho = 0.0
        iters = 0
        for iters in range(maxiter):

            self.applyPreconditioner()   # z_m = invD_m * r_m
            rhoNew = self.dotRZ()        # sum_m r_m · z_m

            self.updateP(rhoNew, rho)    # p_m = z_m + beta p_m
            self.gatherP()
            self.mulA()                  # q_m = A_m @ p
            alpha = rhoNew / self.dotPQ()

            self.axpyUpdate(alpha)       # x_m, r_m
            rho = rhoNew
            if self.residualBelow(atol):
                break

        self.scatterXToMesh()
        self.rollHistory()
        print(f"PCG: iters={iters}  ||r0||={r0:.6e}  atol={atol:.6e}")
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




