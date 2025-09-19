import cupy as cp
from functools import partial

THREADS_PER_BLOCK = 128  # XXXXXXX

class CudaKernelWrapper:
    def __init__(self, kernel_func, sizes, kernel_type='global', n_gpus=None):
        """
        Initialize wrapper for a specific kernel function
        
        Args:
            kernel_func: CuPy kernel function from get_function()
            kernel_type: Type of kernel ('global', 'nodes', 'elements', 'sides')
            n_gpus: Number of GPUs to use (defaults to all available)
        """
        self.n_gpus = n_gpus or cp.cuda.runtime.getDeviceCount()
        self.kernel_func = kernel_func
        self.kernel_type = kernel_type
        
        # Set number of items and threads based on kernel type

        self.num_items = sizes[kernel_type]
        if kernel_type != 'global':
            self.threads_per_block = THREADS_PER_BLOCK
        else:
            self.threads_per_block = 1
        

    def __call__(self, *args, async_launch=False, **kwargs):
        """Execute the kernel on all devices"""

        if self.n_gpus == 1:
            # Launch on one device only
            num_blocks = (self.num_items + self.threads_per_block - 1) // self.threads_per_block
            self.kernel_func((num_blocks,), (self.threads_per_block,), args)
            if not async_launch:
                cp.cuda.Stream.null.synchronize()
        else:
            # Launch on all devices
            for device_id in range(self.n_gpus):
                with cp.cuda.Device(device_id):
                    num_blocks = (self.num_items + self.threads_per_block - 1) // self.threads_per_block
                    self.kernel_func((num_blocks,), (self.threads_per_block,), args)
                    if not async_launch:
                        cp.cuda.Stream.null.synchronize()

    def async_(self, *args, **kwargs):
        """Async version of the kernel execution"""
        return self(*args, async_launch=True, **kwargs)

