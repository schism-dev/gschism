import cupy as cp
import time

# ---------------------------
# Kernel (your test pattern)
# ---------------------------
kernel_code = r'''

__device__ __forceinline__ double compute(double x, double y, double z) {
    double result = 0.0;
    x = 2*x * cos(x);
    y = 2*y;
    z = 2*z;
#pragma unroll
    for (int k=0; k<10; k++) {
        if (x - y < 0.0) {
            result += x * cos(x);
            if (x > 0.0) {
                result -= x * cos(x);
            } else {
                result += x * cos(x);
            }
            result += cos(z);
            result += cos(x);
            result += sin(z);
        } else {
            if (x > 0.0) {
                result += x * cos(x);
            } else {
                result -= x * cos(x);
            }
            result += x * cos(x);
            result += cos(z);
            result += cos(x);
            result += sin(z);
        }
        result += x * cos(x);
        result += cos(z);
        result += cos(x);
        result += sin(z);
    }
    return result;
}

extern "C" __global__
void trig_kernel(double* array1, double* array2, double* array3, int n) {
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    if (tid < n) {
        int i = (int(array1[tid]) * 23) % n;
        int j = (int(array1[tid]) * 47) % n;
        double x = array1[tid];
        double y = array2[i];
        double z = array3[j];

#pragma unroll 32




    //    for (int k=0; k<100; k++) {
    //        y += x * cos(x);
    //        y += cos(z);
    //        y += cos(x);
    //        y += sin(z);
    //    }

     //   array3[tid] = y;
       // array2[tid] = compute(x, y, z);
        array1[tid] = compute(x, y, z);
    }
}
'''

# ---------------------------
# Build modules
# ---------------------------
def build_module(backend):
    if backend == "nvcc":
        options = ( "-O3", "-Xptxas", "-O3", "-use_fast_math")
    elif backend == "nvrtc":
        options = ("--gpu-architecture=compute_86",)
    else:
        raise ValueError(f"Unsupported backend: {backend}")
    return cp.RawModule(
        code=kernel_code,
        backend=backend,
        options=options,
    )

mod_nvrtc = build_module("nvrtc")
mod_nvcc  = build_module("nvcc")

ker_nvrtc = mod_nvrtc.get_function("trig_kernel")
ker_nvcc  = mod_nvcc.get_function("trig_kernel")

# ---------------------------
# Data
# ---------------------------
n = 100000000
x = cp.random.rand(n, dtype=cp.float64)
y = cp.random.rand(n, dtype=cp.float64)
z = cp.random.rand(n, dtype=cp.float64)

# Make copies so both kernels get identical input
x_nvrtc = x.copy()
x_nvcc  = x.copy()
y_nvrtc = y.copy()
y_nvcc  = y.copy()
z_nvrtc = z.copy()
z_nvcc  = z.copy()

threads = 256
blocks = (n + threads - 1) // threads

# ---------------------------
# Timing utility (GPU events)
# ---------------------------
def time_kernel(kernel, arr1, arr2, arr3, n, n_iter=20):
    start = cp.cuda.Event()
    end = cp.cuda.Event()

    # warmup
    for _ in range(5):
        kernel((blocks,), (threads,), (arr1, arr2, arr3, n))

    cp.cuda.Device().synchronize()

    start.record()
    kernel((blocks,), (threads,), (arr1, arr2, arr3, n))
    end.record()

    end.synchronize()
    return cp.cuda.get_elapsed_time(start, end) / n_iter  # ms per run

# ---------------------------
# Run benchmark
# ---------------------------
t_nvrtc = time_kernel(ker_nvrtc, x_nvrtc, y_nvrtc, z_nvrtc, n)
t_nvcc  = time_kernel(ker_nvcc,  x_nvcc,  y_nvcc,  z_nvcc,  n)

print(f"NVRTC time: {t_nvrtc:.3f} ms")
print(f"NVCC  time: {t_nvcc:.3f} ms")
print(f"Speedup (NVRTC / NVCC): {t_nvrtc / t_nvcc:.2f}x")

# ---------------------------
# Sanity check (optional)
# ---------------------------
diff = cp.max(cp.abs(x_nvrtc - x_nvcc))
print(f"Max difference: {diff}")
for i in range(3):
    print(f"x_nvrtc[{i}] = {x_nvrtc[i]}, x_nvcc[{i}] = {x_nvcc[i]}")
    print(f"y_nvrtc[{i}] = {y_nvrtc[i]}, y_nvcc[{i}] = {y_nvcc[i]}")
    print(f"z_nvrtc[{i}] = {z_nvrtc[i]}, z_nvcc[{i}] = {z_nvcc[i]}")