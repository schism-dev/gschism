from __future__ import annotations

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



class PCGSolver:
    def __init__(self, A):
        """Pre-allocate everything for maximum performance."""
        self.n = A.shape[0]
        
        # Pre-allocate all work vectors (no allocations in solve loop)
        self.r    = cp.empty(self.n, dtype=cp.float64)
        self.z    = cp.empty(self.n, dtype=cp.float64) 
        self.p    = cp.empty(self.n, dtype=cp.float64)
        self.q    = cp.empty(self.n, dtype=cp.float64)
        self.x    = cp.empty(self.n, dtype=cp.float64)
        self.diag = cp.empty(self.n, dtype=cp.float64)

        # Set up scalars
        self.one  = np.array(1.0, np.float64)
        self.zero = np.array(0.0, np.float64)        
        
        self.cuda_dtype = _dtype.to_cuda_dtype(self.zero.dtype)

        if not _csr.isspmatrix_csr(A):
            raise ValueError("A must be a sparse matrix")

        # Pre-allocate cuSPARSE descriptors
        self.handle = device.get_cusparse_handle()
        self.desc_A = cusparse.SpMatDescriptor.create(A)
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
       
        
    
    def solve(self, b, x0=None, rtol=1e-5, maxiter=None):
        """Solve with pre-allocated vectors - zero allocations in hot loop."""
        # Update only the values of A (same sparsity)
        # self.desc_A.values[:] = A_values
        # self.diag[:] = diag

        # Initialize solution
        if x0 is None:
            self.x.fill(0)
        else:
            self.x[:] = x0
        
        # Initial residual: r = b - A*x
        # r = A * x (in-place)
        _cusparse.spMV(self.handle, _cusparse.CUSPARSE_OPERATION_NON_TRANSPOSE,
                       self.one.ctypes.data,  self.desc_A.desc, self.desc_x.desc,
                       self.zero.ctypes.data, self.desc_q.desc, self.cuda_dtype, _cusparse.CUSPARSE_MV_ALG_DEFAULT, 
                       self.buff.data.ptr)        
        self.r[:] = b 
        # cublas.axpy(-1.0, self.q, self.r)  
        cp.subtract(b, self.q, out=self.r) # r = b - A*x
        
        # CG iterations with pre-allocated vectors
        rho = 0
        for iters in range(maxiter):
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
            resid = cublas.nrm2(self.r)
            if resid <= rtol * cp.linalg.norm(b):
                break
        
        return self.x, iters












# def cg2(A, b, x0=None, rtol=1e-5, maxiter=None, diag=None):
#     """Uses Conjugate Gradient iteration to solve ``Ax = b``.

#     Args:
#         A (ndarray, spmatrix or LinearOperator): The real or complex matrix of
#             the linear system with shape ``(n, n)``. ``A`` must be a hermitian,
#             positive definitive matrix with type of :class:`cupy.ndarray`,
#             :class:`cupyx.scipy.sparse.spmatrix` or
#             :class:`cupyx.scipy.sparse.linalg.LinearOperator`.
#         b (cupy.ndarray): Right hand side of the linear system with shape
#             ``(n,)`` or ``(n, 1)``.
#         x0 (cupy.ndarray): Starting guess for the solution.
#         rtol (float): Tolerance for convergence.
#         maxiter (int): Maximum number of iterations.
#         diag (cupy.ndarray): Inverse of the diagonal matrix.

#     Returns:
#         tuple:
#             It returns ``x`` (cupy.ndarray) and ``info`` (int) where ``x`` is
#             the converged solution and ``info`` provides convergence
#             information.

#     .. seealso:: :func:`scipy.sparse.linalg.cg`
#     """
#     A, x, b = _make_system(A, x0, b)
#     matvec = A.matvec

#     b_norm = cupy.linalg.norm(b)
#     if b_norm == 0:
#         return b, 0
#     atol = rtol * float(b_norm)

#     r = b - matvec(x)
#     iters = 0
#     rho = 0
#     rho1 = 1
#     beta = 0.0
#     resid = 1e10
#     p = cupy.zeros_like(x)
#     z = cupy.zeros_like(x)
#     while iters < maxiter and resid > atol:
#         # z = diag * r 
        
#         cupy.multiply(r, diag, out=z)
        
#         rho = cublas.dotc(r, z)

#         beta = rho / rho1
#         rho1 = rho
#         # p = z + beta * p
#         cublas.scal(p, beta)
#         cublas.axpy(1.0, z, p)
#         # cublas.axpy(beta, z, p)

#         q = matvec(p)
#         alpha = rho / cublas.dotc(p, q)
#         cublas.axpy( alpha, p, x)
#         cublas.axpy(-alpha, q, r)

#         iters += 1

#         resid = cublas.nrm2(r)

#     info = 0
#     if iters == maxiter and not (resid <= atol):
#         info = iters

#     return x, info





# def _make_system(A, x0, b):
#     """Make a linear system Ax = b

#     Args:
#         A (cupy.ndarray or cupyx.scipy.sparse.spmatrix or
#             cupyx.scipy.sparse.LinearOperator): sparse or dense matrix.
#         x0 (cupy.ndarray): initial guess to iterative method.
#         b (cupy.ndarray): right hand side.

#     Returns:
#         tuple:
#             It returns (A, M, x, b).
#             A (LinaerOperator): matrix of linear system
#             M (LinearOperator): preconditioner
#             x (cupy.ndarray): initial guess
#             b (cupy.ndarray): right hand side.
#     """
#     fast_matvec = _make_fast_matvec(A)
#     A = _interface.aslinearoperator(A)
#     if fast_matvec is not None:
#         A = _interface.LinearOperator(A.shape, matvec=fast_matvec,
#                                       rmatvec=A.rmatvec, dtype=A.dtype)
#     if A.shape[0] != A.shape[1]:
#         raise ValueError('expected square matrix (shape: {})'.format(A.shape))
#     if A.dtype.char not in 'fdFD':
#         raise TypeError('unsupprted dtype (actual: {})'.format(A.dtype))
#     n = A.shape[0]
#     if not (b.shape == (n,) or b.shape == (n, 1)):
#         raise ValueError('b has incompatible dimensions')
#     b = b.astype(A.dtype).ravel()
#     if x0 is None:
#         x = cupy.zeros((n,), dtype=A.dtype)
#     else:
#         if not (x0.shape == (n,) or x0.shape == (n, 1)):
#             raise ValueError('x0 has incompatible dimensions')
#         x = x0.astype(A.dtype).ravel()
 
#     return A, x, b


# def _make_fast_matvec(A):
#     from cupy_backends.cuda.libs import cusparse as _cusparse
#     from cupyx import cusparse

#     if _csr.isspmatrix_csr(A) and cusparse.check_availability('spmv'):
#         handle = device.get_cusparse_handle()
#         op_a = _cusparse.CUSPARSE_OPERATION_NON_TRANSPOSE
#         alpha = numpy.array(1.0, A.dtype)
#         beta = numpy.array(0.0, A.dtype)
#         cuda_dtype = _dtype.to_cuda_dtype(A.dtype)
#         alg = _cusparse.CUSPARSE_MV_ALG_DEFAULT
#         x = cupy.empty((A.shape[0],), dtype=A.dtype)
#         y = cupy.empty((A.shape[0],), dtype=A.dtype)
#         desc_A = cusparse.SpMatDescriptor.create(A)
#         desc_x = cusparse.DnVecDescriptor.create(x)
#         desc_y = cusparse.DnVecDescriptor.create(y)
#         buff_size = _cusparse.spMV_bufferSize(
#             handle, op_a, alpha.ctypes.data, desc_A.desc, desc_x.desc,
#             beta.ctypes.data, desc_y.desc, cuda_dtype, alg)
#         buff = cupy.empty(buff_size, cupy.int8)
#         del x, desc_x, y, desc_y

#         def matvec(x):
#             y = cupy.empty_like(x)
#             desc_x = cusparse.DnVecDescriptor.create(x)
#             desc_y = cusparse.DnVecDescriptor.create(y)
#             _cusparse.spMV(
#                 handle, op_a, alpha.ctypes.data, desc_A.desc, desc_x.desc,
#                 beta.ctypes.data, desc_y.desc, cuda_dtype, alg, buff.data.ptr)
#             return y

#         return matvec
#     return None


