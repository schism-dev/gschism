import numpy as np
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
from cupyx.scipy.sparse import csr_matrix

class CGSolver:
    def __init__(self, mesh, block_size=256):
        self.mesh = mesh
        self.block_size = block_size
        self.num_blocks = (mesh.np + block_size - 1) // block_size
        
        # Compile CUDA kernels
        with open('src/kernels/meshdata.cu', 'r') as f:
            cuda_source = f.read()
        self.mod = SourceModule(cuda_source)
        
        # Get kernel functions
        self.cg_iteration = self.mod.get_function("cgIteration")
        self.reduce_block_results = self.mod.get_function("reduceBlockResults")
        
    def solve_wave_equation(self, x, b, max_iter=1000, tol=1e-6):
        """
        Solve the wave equation using Conjugate Gradient method
        
        Parameters:
        -----------
        x : np.ndarray
            Initial guess vector
        b : np.ndarray
            Right-hand side vector
        max_iter : int
            Maximum number of iterations
        tol : float
            Convergence tolerance
        """
        # Allocate device memory
        d_x = cuda.mem_alloc(x.nbytes)
        d_r = cuda.mem_alloc(x.nbytes)
        d_p = cuda.mem_alloc(x.nbytes)
        d_z = cuda.mem_alloc(x.nbytes)
        d_sp = cuda.mem_alloc(x.nbytes)
        d_block_results = cuda.mem_alloc(3 * self.num_blocks * np.float64().nbytes)
        d_final_results = cuda.mem_alloc(3 * np.float64().nbytes)
        
        # Copy initial data to device
        cuda.memcpy_htod(d_x, x)
        cuda.memcpy_htod(d_r, b)  # Initial residual r = b (since x = 0)
        
        # Set up grid dimensions
        grid_dim = (self.num_blocks, 1, 1)
        block_dim = (self.block_size, 1, 1)
        
        # For reduction kernel
        red_grid_dim = (1, 3, 1)  # One block for each reduction
        red_block_dim = (self.block_size, 1, 1)
        
        # Initialize variables
        iter_count = 0
        rdotr = np.float64(1.0)
        rdotr0 = np.float64(1.0)
        results = np.zeros(3, dtype=np.float64)
        
        # Main CG loop
        while iter_count < max_iter and rdotr > tol * tol * rdotr0:
            # Launch main CG iteration kernel
            shared_mem_size = 3 * self.block_size * np.float64().nbytes
            self.cg_iteration(d_x, d_r, d_p, d_sp, d_z, d_block_results,
                            grid=grid_dim, block=block_dim,
                            shared=shared_mem_size)
            
            # Launch reduction kernel
            self.reduce_block_results(d_block_results,
                                    np.int32(self.num_blocks),
                                    d_final_results,
                                    grid=red_grid_dim,
                                    block=red_block_dim,
                                    shared=self.block_size * np.float64().nbytes)
            
            # Copy results back to host
            cuda.memcpy_dtoh(results, d_final_results)
            
            # Update solution
            rdotr = results[0]
            rdotz = results[1]
            pdotsp = results[2]
            
            if iter_count == 0:
                rdotr0 = rdotr
            
            # Compute alpha and update x, r
            alpha = rdotz / pdotsp
            
            # Update solution and residual on device
            # These could be separate CUDA kernels or part of cgIteration
            self.update_solution(d_x, d_p, alpha, grid=grid_dim, block=block_dim)
            self.update_residual(d_r, d_sp, alpha, grid=grid_dim, block=block_dim)
            
            # Apply preconditioner and update search direction
            beta = results[1] / rdotz  # new_rdotz / old_rdotz
            self.update_search_direction(d_z, d_p, beta, grid=grid_dim, block=block_dim)
            
            iter_count += 1
        
        # Copy final solution back to host
        cuda.memcpy_dtoh(x, d_x)
        
        # Free device memory
        for arr in [d_x, d_r, d_p, d_z, d_sp, d_block_results, d_final_results]:
            arr.free()
        
        return iter_count, rdotr
    
    def update_solution(self, d_x, d_p, alpha, grid, block):
        """Update solution: x = x + alpha*p"""
        self.mod.get_function("vectorAdd")(d_x, d_p, np.float64(alpha), d_x,
                                         grid=grid, block=block)
    
    def update_residual(self, d_r, d_sp, alpha, grid, block):
        """Update residual: r = r - alpha*sp"""
        self.mod.get_function("vectorAdd")(d_r, d_sp, np.float64(-alpha), d_r,
                                         grid=grid, block=block)
    
    def update_search_direction(self, d_z, d_p, beta, grid, block):
        """Update search direction: p = z + beta*p"""
        self.mod.get_function("vectorAdd")(d_z, d_p, np.float64(beta), d_p,
                                         grid=grid, block=block) 