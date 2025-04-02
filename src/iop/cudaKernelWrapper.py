import cupy as cp
from functools import partial

class CudaKernelWrapper:
    def __init__(self, n_gpus=None):
        self.n_gpus = n_gpus or cp.cuda.runtime.getDeviceCount()
        
    def wrap(self, kernel_func, kernel_type='global'):
        """
        Wraps a CuPy kernel function with multi-GPU launch logic
        
        Args:
            kernel_func: CuPy kernel function from get_function()
            kernel_type: Type of kernel ('global', 'nodes', 'elements', 'sides')
        """
        def get_num_items(device_id):
            with cp.cuda.Device(device_id):
                if kernel_type == 'nodes':
                    return mesh.np
                elif kernel_type == 'elements':
                    return mesh.ne
                elif kernel_type == 'sides':
                    return mesh.ns
                else:  # global
                    return 1

        def wrapped_kernel(*args, threads_per_block=256, async_launch=False, **kwargs):
            for device_id in range(self.n_gpus):
                with cp.cuda.Device(device_id):
                    num_items = get_num_items(device_id)
                    num_blocks = (num_items + threads_per_block - 1) // threads_per_block
                    kernel_func((num_blocks,), (threads_per_block,), *args, **kwargs)
                    if not async_launch:
                        cp.cuda.Stream.null.synchronize()

        # Add async version
        wrapped_kernel.async_ = partial(wrapped_kernel, async_launch=True)
        return wrapped_kernel

