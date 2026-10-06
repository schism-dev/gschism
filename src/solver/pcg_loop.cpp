// pcg_loop.cpp -- compiled host loop for PCGSolverMultiGPU (step 1a: one GPU).
//
// Runs exactly the same sequence as PCGSolverMultiGPU._cgLoop (plus the
// initialization pass from solve()), but issues the calls from C++ instead of
// Python:
//   - the kernels are the ones CuPy compiled with NVRTC (pcgsolver._PCG_KERNELS);
//     Python passes their CUfunction handles and this code launches them with
//     cuLaunchKernel,
//   - the SpMV is cuSPARSE (same algorithm as the Python path),
//   - all arrays are the CuPy arrays owned by PCGSolverMultiGPU; Python passes
//     their device addresses once (they never change).
// Built and imported by pcgsolver._load_pcg_loop() (plain g++, no nvcc).

#include <pybind11/pybind11.h>

#include <cuda.h>
#include <cusparse.h>

#include <chrono>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;

// ---------------------------------------------------------------------------
// Error checking
// ---------------------------------------------------------------------------
static void cuCheck(CUresult r, const char* what)
{
    if (r != CUDA_SUCCESS) {
        const char* name = nullptr;
        cuGetErrorName(r, &name);
        throw std::runtime_error(std::string("pcg_loop: ") + what + " failed: " +
                                 (name ? name : "unknown CUDA error"));
    }
}

static void spCheck(cusparseStatus_t s, const char* what)
{
    if (s != CUSPARSE_STATUS_SUCCESS)
        throw std::runtime_error(std::string("pcg_loop: ") + what + " failed: " +
                                 cusparseGetErrorString(s));
}

// Makes a context current for the lifetime of the object.
struct CtxGuard {
    explicit CtxGuard(CUcontext c) { cuCheck(cuCtxPushCurrent(c), "cuCtxPushCurrent"); }
    ~CtxGuard() { CUcontext dummy; cuCtxPopCurrent(&dummy); }
};

// Layout of the accumulator array acc[2][4] (must match A_* in _PCG_KERNELS).
enum { ACC_PQ = 0, ACC_RR = 1, ACC_RZ = 2, ACC_STRIDE = 4 };

// ---------------------------------------------------------------------------
// One partition (one GPU)
// ---------------------------------------------------------------------------
struct Partition {
    int device = -1;
    CUcontext ctx = nullptr;
    CUstream stream = nullptr;

    // sizes / launch configuration
    int N = 0, nblk = 0, tpb = 0;
    int64_t Ntot = 0, nnz = 0;

    // device arrays (owned by CuPy)
    CUdeviceptr data = 0, indices = 0, indptr = 0;           // CSR row block
    CUdeviceptr p_full = 0, p_m = 0, q = 0, r = 0, z = 0, x = 0, invD = 0;
    CUdeviceptr acc = 0;                                     // acc[2][4] for the fused kernels

    // fused kernels (compiled by CuPy/NVRTC, see pcgsolver._PCG_KERNELS)
    CUfunction kInit = nullptr, kUpdateP = nullptr, kDot = nullptr, kUpdateXRZ = nullptr;

    // cuSPARSE
    cusparseHandle_t sp = nullptr;
    cusparseSpMatDescr_t A = nullptr;
    cusparseDnVecDescr_t vecP = nullptr, vecQ = nullptr;
    CUdeviceptr spBuffer = 0;

    // pinned host slot for the convergence read
    double* hostRR = nullptr;
};

static CUdeviceptr getPtr(const py::dict& d, const char* key)
{
    if (!d.contains(key))
        throw std::runtime_error(std::string("pcg_loop: missing pointer '") + key + "'");
    return static_cast<CUdeviceptr>(d[key].cast<uint64_t>());
}

static CUfunction getFunc(const py::dict& d, const char* key)
{
    if (!d.contains(key))
        throw std::runtime_error(std::string("pcg_loop: missing kernel '") + key + "'");
    return reinterpret_cast<CUfunction>(static_cast<uintptr_t>(d[key].cast<uint64_t>()));
}

class NativePCG {
public:
    NativePCG() { cuCheck(cuInit(0), "cuInit"); }

    ~NativePCG()
    {
        // Best effort: ignore errors (e.g. at interpreter shutdown).
        for (auto& p : parts_) {
            if (!p.ctx) continue;
            if (cuCtxPushCurrent(p.ctx) != CUDA_SUCCESS) continue;
            if (p.vecQ) cusparseDestroyDnVec(p.vecQ);
            if (p.vecP) cusparseDestroyDnVec(p.vecP);
            if (p.A) cusparseDestroySpMat(p.A);
            if (p.sp) cusparseDestroy(p.sp);
            if (p.spBuffer) cuMemFree(p.spBuffer);
            if (p.hostRR) cuMemFreeHost(p.hostRR);
            if (p.stream) cuStreamDestroy(p.stream);
            CUcontext dummy;
            cuCtxPopCurrent(&dummy);
            CUdevice dev;
            if (cuDeviceGet(&dev, p.device) == CUDA_SUCCESS) cuDevicePrimaryCtxRelease(dev);
        }
    }

    // Register one partition. ptrs/kernels: dicts of device addresses / CUfunction handles.
    void addPartition(int device, int N, int64_t Ntot, int64_t nnz, int nblk, int tpb,
                      py::dict ptrs, py::dict kernels)
    {
        if (!parts_.empty())
            throw std::runtime_error("pcg_loop: step 1a supports a single partition (one GPU)");

        Partition p;
        p.device = device;
        p.N = N; p.Ntot = Ntot; p.nnz = nnz; p.nblk = nblk; p.tpb = tpb;

        p.data    = getPtr(ptrs, "data");
        p.indices = getPtr(ptrs, "indices");
        p.indptr  = getPtr(ptrs, "indptr");
        p.p_full  = getPtr(ptrs, "p_full");
        p.p_m     = getPtr(ptrs, "p_m");
        p.q       = getPtr(ptrs, "q");
        p.r       = getPtr(ptrs, "r");
        p.z       = getPtr(ptrs, "z");
        p.x       = getPtr(ptrs, "x");
        p.invD    = getPtr(ptrs, "invD");
        p.acc     = getPtr(ptrs, "acc");

        p.kInit      = getFunc(kernels, "f_init");
        p.kUpdateP   = getFunc(kernels, "f_update_p");
        p.kDot       = getFunc(kernels, "f_dot");
        p.kUpdateXRZ = getFunc(kernels, "f_update_xrz");

        // Use the device's primary context: the same one CuPy uses.
        CUdevice dev;
        cuCheck(cuDeviceGet(&dev, device), "cuDeviceGet");
        cuCheck(cuDevicePrimaryCtxRetain(&p.ctx, dev), "cuDevicePrimaryCtxRetain");
        CtxGuard g(p.ctx);

        cuCheck(cuStreamCreate(&p.stream, CU_STREAM_NON_BLOCKING), "cuStreamCreate");
        cuCheck(cuMemAllocHost(reinterpret_cast<void**>(&p.hostRR), sizeof(double)),
                "cuMemAllocHost");

        spCheck(cusparseCreate(&p.sp), "cusparseCreate");
        spCheck(cusparseSetStream(p.sp, reinterpret_cast<cudaStream_t>(p.stream)),
                "cusparseSetStream");
        spCheck(cusparseCreateCsr(&p.A, N, Ntot, nnz,
                                  reinterpret_cast<void*>(p.indptr),
                                  reinterpret_cast<void*>(p.indices),
                                  reinterpret_cast<void*>(p.data),
                                  CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                  CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F),
                "cusparseCreateCsr");
        spCheck(cusparseCreateDnVec(&p.vecP, Ntot, reinterpret_cast<void*>(p.p_full), CUDA_R_64F),
                "cusparseCreateDnVec(p)");
        spCheck(cusparseCreateDnVec(&p.vecQ, N, reinterpret_cast<void*>(p.q), CUDA_R_64F),
                "cusparseCreateDnVec(q)");

        size_t bufSize = 0;
        const double one = 1.0, zero = 0.0;
        spCheck(cusparseSpMV_bufferSize(p.sp, CUSPARSE_OPERATION_NON_TRANSPOSE, &one, p.A, p.vecP,
                                        &zero, p.vecQ, CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT,
                                        &bufSize),
                "cusparseSpMV_bufferSize");
        if (bufSize > 0) cuCheck(cuMemAlloc(&p.spBuffer, bufSize), "cuMemAlloc(spmv buffer)");

        parts_.push_back(p);
    }

    // Initialization pass + CG iterations with the fused kernels:
    //   f_update_p, cuSPARSE SpMV, f_dot, f_update_xrz   (per iteration)
    // Returns it+1 at the first check where r.r <= atol2, or maxiter (same
    // convention as the Python loop). Python must have synchronized its
    // default stream before calling (r and x are prepared there).
    int solve(double atol2, int maxiter, int checkEvery, double timeoutSeconds)
    {
        if (parts_.size() != 1)
            throw std::runtime_error("pcg_loop: no partition registered");
        if (checkEvery < 1) checkEvery = 1;

        Partition& p = parts_[0];
        int result = maxiter;
        {
            py::gil_scoped_release noGil;       // pure C++/CUDA from here on
            CtxGuard g(p.ctx);

            // Init: acc = 0;  z = invD r;  r.r and r.z into acc[1]
            cuCheck(cuMemsetD8Async(p.acc, 0, 8 * sizeof(double), p.stream), "cuMemsetD8Async");
            launch(p, p.kInit, p.nblk, p.tpb, {&p.r, &p.z, &p.invD, &p.acc, &p.N});

            for (int it = 0; it < maxiter; ++it) {
                int k = it;
                // p = z + beta p   (beta computed in the kernel)
                launch(p, p.kUpdateP, p.nblk, p.tpb, {&p.z, &p.p_m, &p.acc, &p.N, &k});
                // q = A p   (cuSPARSE)
                spmv(p);
                // p.q   (atomic accumulation)
                launch(p, p.kDot, p.nblk, p.tpb, {&p.p_m, &p.q, &p.acc, &p.N, &k});
                // alpha in the kernel; x, r, z update; r.r and r.z (atomic accumulation)
                launch(p, p.kUpdateXRZ, p.nblk, p.tpb,
                       {&p.x, &p.r, &p.z, &p.p_m, &p.q, &p.invD, &p.acc, &p.N, &k});

                if ((it + 1) % checkEvery == 0 || it == maxiter - 1) {
                    if (readRR(p, it, timeoutSeconds) <= atol2) { result = it + 1; break; }
                }
            }
            waitStream(p, timeoutSeconds);      // leave nothing queued for the caller
        }
        return result;
    }

private:
    std::vector<Partition> parts_;
    int numParts() const { return static_cast<int>(parts_.size()); }

    void launch(Partition& p, CUfunction f, int grid, int block, std::vector<void*> args)
    {
        cuCheck(cuLaunchKernel(f, grid, 1, 1, block, 1, 1, 0, p.stream, args.data(), nullptr),
                "cuLaunchKernel");
    }

    void spmv(Partition& p)
    {
        const double one = 1.0, zero = 0.0;
        spCheck(cusparseSpMV(p.sp, CUSPARSE_OPERATION_NON_TRANSPOSE, &one, p.A, p.vecP, &zero,
                             p.vecQ, CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT,
                             reinterpret_cast<void*>(p.spBuffer)),
                "cusparseSpMV");
    }

    // Wait for the stream with a time limit, so a stuck GPU raises instead of hanging.
    void waitStream(Partition& p, double timeoutSeconds)
    {
        auto t0 = std::chrono::steady_clock::now();
        for (;;) {
            CUresult r = cuStreamQuery(p.stream);
            if (r == CUDA_SUCCESS) return;
            if (r != CUDA_ERROR_NOT_READY) cuCheck(r, "cuStreamQuery");
            double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            if (dt > timeoutSeconds)
                throw std::runtime_error("pcg_loop: GPU work did not finish within the time limit");
        }
    }

    // r.r after iteration it lives in acc[it & 1][A_RR]
    double readRR(Partition& p, int it, double timeoutSeconds)
    {
        const CUdeviceptr src = p.acc + ((it & 1) * ACC_STRIDE + ACC_RR) * sizeof(double);
        cuCheck(cuMemcpyDtoHAsync(p.hostRR, src, sizeof(double), p.stream),
                "cuMemcpyDtoHAsync");
        waitStream(p, timeoutSeconds);
        return *p.hostRR;
    }
};

PYBIND11_MODULE(pcg_loop, m)
{
    m.doc() = "Compiled host loop for PCGSolverMultiGPU (one GPU, fused kernels)";
    py::class_<NativePCG>(m, "NativePCG")
        .def(py::init<>())
        .def("add_partition", &NativePCG::addPartition,
             py::arg("device"), py::arg("N"), py::arg("Ntot"), py::arg("nnz"),
             py::arg("nblk"), py::arg("tpb"), py::arg("ptrs"), py::arg("kernels"))
        .def("solve", &NativePCG::solve,
             py::arg("atol2"), py::arg("maxiter"), py::arg("check_every"),
             py::arg("timeout_seconds") = 60.0);
}