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
        # XXXXX It is not using p or q
        # Computes q = A * p 
        _cusparse.spMV(self.handle, _cusparse.CUSPARSE_OPERATION_NON_TRANSPOSE,
                           self.one.ctypes.data, self.desc_A.desc, self.desc_p.desc,
                           self.zero.ctypes.data, self.desc_q.desc, self.cuda_dtype, _cusparse.CUSPARSE_MV_ALG_DEFAULT, 
                           self.buff.data.ptr)

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

    def print_stats(self):
        self.solver_profiler.print_stats()

class PCGSolverMultiGPU:
    def __init__(self, configCuda: CudaDeviceConfig, N: list[int], NNZ: list[int]):
        """Pre-allocate everything for maximum performance."""

        self.configCuda = configCuda
        self.deviceCount = configCuda.deviceCount
        self.N = N
        self.NNZ = NNZ

        self.A       = [[None]*self.deviceCount for _ in range(self.deviceCount)]
        self.b       = [[None]*self.deviceCount]
        self.indices = [[None]*self.deviceCount]
        self.indptr  = [[None]*self.deviceCount]
        self.diag    = [[None]*self.deviceCount]
        self.r       = [[None]*self.deviceCount]
        self.z       = [[None]*self.deviceCount]
        self.p       = [[None]*self.deviceCount]
        self.q       = [[None]*self.deviceCount]

        for i, device_id in enumerate(configCuda.devices):
            
            with cp.cuda.Device(device_id):
                for j in range(self.deviceCount):

                    self.data     = cp.zeros(self.NNZ[i], dtype=cp.float64)
                    self.b        = cp.zeros(self.N[i],  dtype=cp.float64)


                    self.indices= cp.zeros(self.NNZ[i],    dtype=cp.int32)
                    self.indptr = cp.zeros(self.N[i] + 1, dtype=cp.int32)
                    self.indptr[self.N[i]] = self.NNZ[i]

                    self.A[i][j]   = csr_matrix((self.data  , self.indices, self.indptr), shape=(self.N[i], self.N[j]))

                    
                    # Pre-allocate all work vectors (no allocations in solve loop)
                    self.diag = cp.empty(self.N[i], dtype=cp.float64)
                    self.r    = cp.empty(self.N[i], dtype=cp.float64)
                    self.z    = cp.empty(self.N[i], dtype=cp.float64) 
                    self.p    = cp.empty(self.N[i], dtype=cp.float64)
                    self.q    = cp.empty(self.N[i], dtype=cp.float64)
                    self.x    = cp.empty(self.N[i], dtype=cp.float64)
                    self.x0   = cp.empty(self.N[i], dtype=cp.float64)
                    self.x1   = cp.empty(self.N[i], dtype=cp.float64)
                    self.x2   = cp.empty(self.N[i], dtype=cp.float64)

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
        # XXXXX It is not using p or q
        # Computes q = A * p 
        _cusparse.spMV(self.handle, _cusparse.CUSPARSE_OPERATION_NON_TRANSPOSE,
                           self.one.ctypes.data, self.desc_A.desc, self.desc_p.desc,
                           self.zero.ctypes.data, self.desc_q.desc, self.cuda_dtype, _cusparse.CUSPARSE_MV_ALG_DEFAULT, 
                           self.buff.data.ptr)

 
    def solveMultiGPU(self, rtol=1e-5, maxiter = None, ensureSymmetric = False):

        # Funtions below are not implemented yet. They will operate on class managed
        # data structures.
        self.packIntoCacheFriendlyFormat()

        self.extrapolateGuess()

        self.computeResiduals()          # ri = bi - Aii xi - sum_j Aij xj

        atol = self.computeAtol(rtol)    # all-reduce of norms
        rho = 0.0
        for iters in range(maxiter):
            self.applyPreconditioner()   # zi = Diinv * ri
            rhoNew = self.dotRZ()        # sum_i ri·zi
            self.updateP(rhoNew, rho)    # pi = zi + beta pi
            self.matvecQ()               # qi = Aii pi + sum_j Aij pj
            alpha = rhoNew / self.dotPQ()
            self.axpyUpdate(alpha)       # xi, ri
            rho = rhoNew
            if self.residualBelow(atol):
                break

        self.scatterXToMesh()
        self.rollHistory()
        return iters


    def print_stats(self):
        self.solver_profiler.print_stats()




