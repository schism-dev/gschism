import cupy as cp
from functools import partial
from schism_consts import VISTA

THREADS_PER_BLOCK = 64  # XXXXXXX


class CudaDeviceConfig:
    """Information about the devices available in the system."""
    def __init__(self, devices=None):
        if devices is None:
            self.devices = list(range(cp.cuda.runtime.getDeviceCount()))
        else:
            self.devices = list(devices)

        self.deviceCount = len(self.devices)


configCuda = None
def initCudaDeviceConfig():
    global configCuda

    # THIS IS FOR TESTING ONLY XXXXXXX
    if VISTA:
        devices = (0, )
    else:
        devices = (0,)


    configCuda = CudaDeviceConfig(devices)

    # Enable P2P both ways (test only) XXXXXXXXX
    try:
        with cp.cuda.Device(0):
            cp.cuda.runtime.deviceEnablePeerAccess(1)
        with cp.cuda.Device(1):
            cp.cuda.runtime.deviceEnablePeerAccess(0)
    except:
        print("Error enabling peer access between devices")

    return configCuda


class CudaKernelWrapper:
    def __init__(self, kernel_func, sizes, kernel_type='global'):
        """
        Initialize wrapper for a specific kernel function
        
        Args:
            kernel_func: CuPy kernel function from get_function()
            kernel_type: Type of kernel ('global', 'nodes', 'elements', 'sides'). Can also
                         include 'onlycurrent' to indicate that the kernel should only be executed on the current device.
        """

        tags = [t.strip() for t in kernel_type.split(',') if t.strip()]

        self.onlyCurrent = 'onlycurrent' in tags
        if self.onlyCurrent:
            tags.remove('onlycurrent')

        if len(tags) > 1 or len(tags) == 0:
            raise ValueError(f"Too many tags in kernel type: {kernel_type}")

        self.n_gpus = configCuda.deviceCount
        self.kernel_func = kernel_func
        self.kernel_type = tags[0]
        
        # Set number of items and threads based on kernel type

        self.num_items = sizes[self.kernel_type]
        if self.kernel_type != 'global':
            self.threads_per_block = THREADS_PER_BLOCK
        else:
            self.threads_per_block = 1
        

    def __call__(self, *args, async_launch=False, **kwargs):
        """Execute the kernel on all devices"""


        numGPUs = configCuda.deviceCount
        if numGPUs == 1 or self.onlyCurrent:
            # Launch on one device only
            num_blocks = (self.num_items + self.threads_per_block - 1) // self.threads_per_block
            self.kernel_func((num_blocks,), (self.threads_per_block,), args)
            if not async_launch:
                cp.cuda.Stream.null.synchronize()
        else:
            print(f">>>>>*** 1 launching {self.kernel_func.__name__} on all devices")
            # Launch on all devices
            for device_id in configCuda.devices:
                with cp.cuda.Device(device_id):
                    num_blocks = (self.num_items + self.threads_per_block - 1) // self.threads_per_block
                    print(f">>>>>*** 2 launching on device {device_id}, num_blocks = {num_blocks}, threads_per_block = {self.threads_per_block}")
                    self.kernel_func((num_blocks,), (self.threads_per_block,), args)


            if not async_launch:
                for device_id in configCuda.devices:
                    with cp.cuda.Device(device_id):
                        cp.cuda.Stream.null.synchronize()

    def async_(self, *args, **kwargs):
        """Async version of the kernel execution"""
        return self(*args, async_launch=True, **kwargs)

