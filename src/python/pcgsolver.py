from __future__ import annotations
import python.cuda_functions as cufunc

import numpy as np

import cupy as cp
from cupy import cublas
from cupy._core import _dtype
from cupy.cuda import device
from cupy_backends.cuda.libs import cublas as _cublas
from cupyx.scipy.sparse import _csr
from cupyx.scipy.sparse.linalg import _interface
from cupy_backends.cuda.libs import cusparse as _cusparse
from cupyx import cusparse
from cupyx.scipy.sparse import csr_matrix


from line_profiler import LineProfiler

class PCGSolver:
    def __init__(self):
        """Pre-allocate everything for maximum performance."""

        N = 1457889
        NN = 10141117

        # # Create independent profiler for solver
        # self.solver_profiler = LineProfiler()
        # self.solveprof = self.solver_profiler(self.solve2)
        # self.solver_profiler.add_function(self.solve) 

        self.data     = cp.zeros(NN, dtype=cp.float64)
        self.data_f   = cp.zeros(NN, dtype=cp.float32)
        self.b        = cp.zeros(N,  dtype=cp.float64)
        self.b_f      = cp.zeros(N,  dtype=cp.float32)
        self.diag     = cp.zeros(N,  dtype=cp.float64)
        self.diag_f   = cp.zeros(N,  dtype=cp.float32)

        self.indices= cp.zeros(NN,    dtype=cp.int32)
        self.indptr = cp.zeros(N + 1, dtype=cp.int32)
        self.indptr[N] = NN

        self.A   = csr_matrix((self.data  , self.indices, self.indptr), shape=(N, N))
        self.A_f = csr_matrix((self.data_f, self.indices, self.indptr), shape=(N, N))

        self.n = N

        
        # Pre-allocate all work vectors (no allocations in solve loop)
        self.r    = cp.empty(self.n, dtype=cp.float64)
        self.z    = cp.empty(self.n, dtype=cp.float64) 
        self.p    = cp.empty(self.n, dtype=cp.float64)
        self.q    = cp.empty(self.n, dtype=cp.float64)
        self.x    = cp.empty(self.n, dtype=cp.float64)

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
        buff_size = _cusparse.spMV_bufferSize(
            self.handle, _cusparse.CUSPARSE_OPERATION_NON_TRANSPOSE,
            self.one.ctypes.data, self.desc_A.desc, self.desc_p.desc,
            self.zero.ctypes.data, self.desc_q.desc, self.cuda_dtype, alg)
        self.buff = cp.empty(buff_size, cp.int8)


        # Repeat for float
        self.r_f    = cp.zeros(self.n, dtype=cp.float32)
        self.z_f    = cp.zeros(self.n, dtype=cp.float32) 
        self.p_f    = cp.zeros(self.n, dtype=cp.float32)
        self.q_f    = cp.zeros(self.n, dtype=cp.float32)
        self.x_f    = cp.zeros(self.n, dtype=cp.float32)

        # Previous solutions
        self.x1_f   = cp.zeros(self.n, dtype=cp.float32)
        self.x2_f   = cp.zeros(self.n, dtype=cp.float32)

        # Set up scalars
        self.one_f  = np.array(1.0, np.float32)
        self.zero_f = np.array(0.0, np.float32)        
        
        self.cuda_dtype_f = _dtype.to_cuda_dtype(self.zero_f.dtype)

        # Pre-allocate cuSPARSE descriptors
        self.handle = device.get_cusparse_handle()
        self.desc_A_f = cusparse.SpMatDescriptor.create(self.A_f)
        self.desc_p_f = cusparse.DnVecDescriptor.create(self.p_f)
        self.desc_q_f = cusparse.DnVecDescriptor.create(self.q_f)
        self.desc_r_f = cusparse.DnVecDescriptor.create(self.r_f)
        self.desc_x_f = cusparse.DnVecDescriptor.create(self.x_f)
        
        # Pre-allocate buffer
        buff_size_f = _cusparse.spMV_bufferSize(
            self.handle, _cusparse.CUSPARSE_OPERATION_NON_TRANSPOSE,
            self.one_f.ctypes.data, self.desc_A_f.desc, self.desc_p_f.desc,
            self.zero_f.ctypes.data, self.desc_q_f.desc, self.cuda_dtype_f, alg)

        self.buff_f = cp.empty(buff_size_f, cp.int8)
        
               

    # def solve(self, rtol=1e-5, maxiter = None):
    #     self.solveprof(rtol, maxiter)

    def solve(self, rtol=1e-5, maxiter = None):
        """Solve with pre-allocated vectors - zero allocations in hot loop."""
 
        
        cufunc.setupCGSolvers(self.diag, self.data, self.b, self.diag_f, self.data_f, self.b_f, self.indices, self.indptr)

        # Initialize solution using parabolic extrapolation from previous solutions
        self.x_f[:] = 2.0*self.x1_f - self.x2_f

        if maxiter is None:
            maxiter = self.n * 10

        old_iters = 0

        atol_f = 1e-6 * cp.linalg.norm(self.b_f)
        atol   = rtol * cp.linalg.norm(self.b)

        # Starts solving on single precision
        # ----------------------------------

        # Initial residual: r = b - A*x
        # r = A * x (in-place)
        _cusparse.spMV(self.handle, _cusparse.CUSPARSE_OPERATION_NON_TRANSPOSE,
                       self.one_f.ctypes.data,  self.desc_A_f.desc, self.desc_x_f.desc,
                       self.zero_f.ctypes.data, self.desc_q_f.desc, self.cuda_dtype_f, _cusparse.CUSPARSE_MV_ALG_DEFAULT, 
                       self.buff_f.data.ptr)        
        self.r_f[:] = self.b_f 
        cublas.axpy(-1.0, self.q_f, self.r_f)  
        # cp.subtract(b_f, self.q_f, out=self.r_f) # r = b - A*x
        
        # CG iterations with pre-allocated vectors
        rho = 0
        for iters in range(maxiter):
            # z = diag * r (in-place)
            cp.multiply(self.r_f, self.diag_f, out=self.z_f)
            
            rho_new = cublas.dotc(self.r_f, self.z_f)
            beta = rho_new / rho if rho != 0 else 0
            
            # p = z + beta * p (in-place)
            self.p_f *= beta
            self.p_f += self.z_f
            
            # q = A * p (in-place)
            _cusparse.spMV(self.handle, _cusparse.CUSPARSE_OPERATION_NON_TRANSPOSE,
                           self.one_f.ctypes.data, self.desc_A_f.desc, self.desc_p_f.desc,
                           self.zero_f.ctypes.data, self.desc_q_f.desc, self.cuda_dtype_f, _cusparse.CUSPARSE_MV_ALG_DEFAULT, 
                           self.buff_f.data.ptr)
        
            
            alpha = rho_new / cublas.dotc(self.p_f, self.q_f)
            
            # x += alpha * p, r -= alpha * q (in-place)
            cublas.axpy( alpha, self.p_f, self.x_f)
            cublas.axpy(-alpha, self.q_f, self.r_f)
            
            rho = rho_new
            
            # Check convergence
            if iters % 10 == 0:  # Avoids checking every iteration because it's expensive
                resid = cublas.nrm2(self.r_f)
                if resid <= atol_f:
                    old_iters = iters + 1
                    break

        if iters >= maxiter and not (resid <= atol_f):
            raise ValueError("PCG solver did not converge")

        # Uses the solution from single precision as the initial guess for double precision
        self.x[:] = self.x_f

        # Initial residual: r = b - A*x
        # r = A * x (in-place)
        _cusparse.spMV(self.handle, _cusparse.CUSPARSE_OPERATION_NON_TRANSPOSE,
                       self.one.ctypes.data,  self.desc_A.desc, self.desc_x.desc,
                       self.zero.ctypes.data, self.desc_q.desc, self.cuda_dtype, _cusparse.CUSPARSE_MV_ALG_DEFAULT, 
                       self.buff.data.ptr)        
        self.r[:] = self.b 
        cublas.axpy(-1.0, self.q, self.r)  
        # cp.subtract(b, self.q, out=self.r) # r = b - A*x
        
        # CG iterations with pre-allocated vectors
        rho = 0
        for iters in range(old_iters, maxiter):
            # z = diag * r (in-place)
            cp.multiply(self.r, self.diag, out=self.z)
            
            rho_new = cublas.dotc(self.r, self.z)
            beta = rho_new / rho if rho != 0 else 0
            
            # p = z + beta * p (in-place)
            self.p *= beta
            self.p += self.z
            
            # q = A * p (in-place)
            _cusparse.spMV(self.handle, _cusparse.CUSPARSE_OPERATION_NON_TRANSPOSE,
                           self.one.ctypes.data, self.desc_A.desc, self.desc_p.desc,
                           self.zero.ctypes.data, self.desc_q.desc, self.cuda_dtype, _cusparse.CUSPARSE_MV_ALG_DEFAULT, 
                           self.buff.data.ptr)
        
            
            alpha = rho_new / cublas.dotc(self.p, self.q)
            
            # x += alpha * p, r -= alpha * q (in-place)
            cublas.axpy( alpha, self.p, self.x)
            cublas.axpy(-alpha, self.q, self.r)
            
            rho = rho_new
            
            # Check convergence
            if iters % 10 == 0:  # Avoids checking every iteration because it's expensive
                # resid = cublas.nrm2(self.r)
                resid = cublas.dotc(self.r, self.r) # nrm2 is more stable but more expensive. At this point there is no risk of overflow, so we use dotc.
                if resid <= atol*atol:
                    break

        # store the previous solutions
        self.x2_f[:] = self.x1_f
        self.x1_f[:] = self.x_f
        self.x_f [:] = self.x

        cufunc.solverCopyResult(self.x)
        
        return iters


    def print_stats(self):
        self.solver_profiler.print_stats()






