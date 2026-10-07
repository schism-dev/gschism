// pcg_loop.cpp -- compiled host loop for PCGSolverMultiGPU (1 or more GPUs).
//
// One host thread per GPU; each issues only its own GPU's work on its own stream.
// Per iteration on GPU m (fused kernels from pcgsolver._PCG_KERNELS, launched by
// their CUfunction handles; SpMV by cuSPARSE):
//
//   barrier; wait evXRZ of every other GPU      (their f_update_xrz(k-1) is done)
//   f_update_p        beta from all GPUs' r.z partials (read over NVLink)
//   record evP; barrier; wait evP of every other GPU
//                     (direct exchange: f_update_p already stored the values other GPUs
//                      use into their full p; otherwise: copy the others' p segments here)
//   cuSPARSE SpMV     q = A p
//   f_dot             own p.q partial
//   record evDot; barrier; wait evDot of every other GPU
//   f_update_xrz      alpha from all GPUs' p.q partials; own r.r and r.z partials
//   record evXRZ
//
// The host barriers guarantee a thread never waits on another thread's event
// before that thread has recorded it for the current phase. Every 10 iterations
// each thread reads its own r.r partial; after a barrier all threads add the
// partials in the same order and stop at the same iteration.
// With one GPU there are no barriers, waits or copies.
//
// Built and imported by pcgsolver._load_pcg_loop() (plain g++, no nvcc).

#include <pybind11/pybind11.h>

#include <cuda.h>
#include <cusparse.h>

#include <atomic>
#include <chrono>
#include <cstdint>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
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
// Spin barrier for the host threads; abortable so an error in one thread
// releases the others instead of leaving them waiting forever.
// ---------------------------------------------------------------------------
class SpinBarrier {
public:
    SpinBarrier(int n, const std::atomic<bool>& abort) : n_(n), abort_(abort) {}
    void wait()
    {
        if (n_ <= 1) return;
        const int gen = gen_.load(std::memory_order_acquire);
        if (count_.fetch_add(1, std::memory_order_acq_rel) == n_ - 1) {
            count_.store(0, std::memory_order_relaxed);
            gen_.fetch_add(1, std::memory_order_release);
        } else {
            // Spin briefly (lowest latency when every thread has its own core), then
            // yield the core while waiting, so the threads still progress if the job
            // has fewer cores than GPU threads.
            for (int spins = 0; gen_.load(std::memory_order_acquire) == gen; ++spins) {
                if (abort_.load(std::memory_order_relaxed))
                    throw std::runtime_error("pcg_loop: aborted because another GPU's thread failed");
                if (spins > 2000) std::this_thread::yield();
            }
        }
    }
private:
    int n_;
    const std::atomic<bool>& abort_;
    std::atomic<int> count_{0};
    std::atomic<int> gen_{0};
};

// ---------------------------------------------------------------------------
// One partition (one GPU)
// ---------------------------------------------------------------------------
struct Partition {
    int device = -1;
    CUcontext ctx = nullptr;
    CUstream stream = nullptr;
    CUevent evP = nullptr, evDot = nullptr, evXRZ = nullptr;

    // sizes / launch configuration
    int N = 0, nblk = 0, tpb = 0;
    int64_t Ntot = 0, nnz = 0, offset = 0;      // offset: first packed column/row of this GPU

    // device arrays (owned by CuPy)
    CUdeviceptr data = 0, indices = 0, indptr = 0;            // CSR row block (packed columns)
    CUdeviceptr p_full = 0, p_m = 0, q = 0, r = 0, z = 0, x = 0, invD = 0;
    CUdeviceptr acc = 0;                                      // acc[2][4] of this GPU
    CUdeviceptr accAll = 0;                                   // device array: acc pointer of every GPU
    CUdeviceptr sendMask = 0;                                 // one byte per row: bit q = GPU q uses it
    CUdeviceptr pAll = 0;                                     // device array: full-p pointer of every GPU

    // fused kernels (compiled by CuPy/NVRTC on this GPU)
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
            if (p.accAll) cuMemFree(p.accAll);
            if (p.pAll) cuMemFree(p.pAll);
            if (p.hostRR) cuMemFreeHost(p.hostRR);
            if (p.evP) cuEventDestroy(p.evP);
            if (p.evDot) cuEventDestroy(p.evDot);
            if (p.evXRZ) cuEventDestroy(p.evXRZ);
            if (p.stream) cuStreamDestroy(p.stream);
            CUcontext dummy;
            cuCtxPopCurrent(&dummy);
            CUdevice dev;
            if (cuDeviceGet(&dev, p.device) == CUDA_SUCCESS) cuDevicePrimaryCtxRelease(dev);
        }
    }

    // Register partition m (in order m = 0, 1, ...).
    void addPartition(int device, int N, int64_t Ntot, int64_t nnz, int64_t offset,
                      int nblk, int tpb, py::dict ptrs, py::dict kernels)
    {
        if (finalized_)
            throw std::runtime_error("pcg_loop: partitions must be added before the first solve");

        Partition p;
        p.device = device;
        p.N = N; p.Ntot = Ntot; p.nnz = nnz; p.offset = offset; p.nblk = nblk; p.tpb = tpb;

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
        p.sendMask = ptrs.contains("send_mask")
                         ? static_cast<CUdeviceptr>(ptrs["send_mask"].cast<uint64_t>()) : 0;

        p.kInit      = getFunc(kernels, "f_init");
        p.kUpdateP   = getFunc(kernels, "f_update_p");
        p.kDot       = getFunc(kernels, "f_dot");
        p.kUpdateXRZ = getFunc(kernels, "f_update_xrz");

        // The device's primary context: the same one CuPy uses.
        CUdevice dev;
        cuCheck(cuDeviceGet(&dev, device), "cuDeviceGet");
        cuCheck(cuDevicePrimaryCtxRetain(&p.ctx, dev), "cuDevicePrimaryCtxRetain");
        CtxGuard g(p.ctx);

        cuCheck(cuStreamCreate(&p.stream, CU_STREAM_NON_BLOCKING), "cuStreamCreate");
        cuCheck(cuEventCreate(&p.evP, CU_EVENT_DISABLE_TIMING), "cuEventCreate");
        cuCheck(cuEventCreate(&p.evDot, CU_EVENT_DISABLE_TIMING), "cuEventCreate");
        cuCheck(cuEventCreate(&p.evXRZ, CU_EVENT_DISABLE_TIMING), "cuEventCreate");
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

    // Initialization pass + CG iterations on all registered GPUs.
    // Returns it+1 at the first check where the total r.r <= atol2, or maxiter.
    // Python must have synchronized its default streams before calling.
    int solve(double atol2, int maxiter, int checkEvery, double timeoutSeconds, bool preprocess,
              bool direct)
    {
        if (parts_.empty())
            throw std::runtime_error("pcg_loop: no partition registered");
        if (checkEvery < 1) checkEvery = 1;
        finalize();

        const int P = static_cast<int>(parts_.size());
        std::atomic<bool> abort{false};
        SpinBarrier barrier(P, abort);
        std::vector<int> result(P, maxiter);
        std::mutex errMutex;
        std::string firstError;

        auto worker = [&](int m) {
            try {
                result[m] = runPartition(m, atol2, maxiter, checkEvery, timeoutSeconds,
                                         preprocess, direct, barrier);
            } catch (const std::exception& e) {
                std::lock_guard<std::mutex> lock(errMutex);
                if (firstError.empty()) firstError = e.what();
                abort.store(true);
            }
        };

        {
            py::gil_scoped_release noGil;       // pure C++/CUDA from here on
            if (P == 1) {
                worker(0);
            } else {
                std::vector<std::thread> threads;
                for (int m = 0; m < P; ++m) threads.emplace_back(worker, m);
                for (auto& t : threads) t.join();
            }
        }
        if (!firstError.empty()) throw std::runtime_error(firstError);
        return result[0];
    }

private:
    std::vector<Partition> parts_;
    bool finalized_ = false;

    // Once, before the first solve: peer access between all GPUs, and on each GPU
    // the device array with the acc pointers of all GPUs.
    void finalize()
    {
        if (finalized_) return;
        const int P = static_cast<int>(parts_.size());
        std::vector<uint64_t> accPtrs(P), pPtrs(P);
        for (int j = 0; j < P; ++j) { accPtrs[j] = parts_[j].acc; pPtrs[j] = parts_[j].p_full; }

        for (int a = 0; a < P; ++a) {
            CtxGuard g(parts_[a].ctx);
            for (int b = 0; b < P; ++b) {
                if (a == b) continue;
                CUdevice da, db;
                cuCheck(cuDeviceGet(&da, parts_[a].device), "cuDeviceGet");
                cuCheck(cuDeviceGet(&db, parts_[b].device), "cuDeviceGet");
                int can = 0;
                cuCheck(cuDeviceCanAccessPeer(&can, da, db), "cuDeviceCanAccessPeer");
                if (!can)
                    throw std::runtime_error("pcg_loop: GPU " + std::to_string(parts_[a].device) +
                                             " cannot access GPU " +
                                             std::to_string(parts_[b].device) + " (peer access)");
                CUresult r = cuCtxEnablePeerAccess(parts_[b].ctx, 0);
                if (r != CUDA_SUCCESS && r != CUDA_ERROR_PEER_ACCESS_ALREADY_ENABLED)
                    cuCheck(r, "cuCtxEnablePeerAccess");
            }
            cuCheck(cuMemAlloc(&parts_[a].accAll, P * sizeof(uint64_t)), "cuMemAlloc(accAll)");
            cuCheck(cuMemcpyHtoD(parts_[a].accAll, accPtrs.data(), P * sizeof(uint64_t)),
                    "cuMemcpyHtoD(accAll)");
            cuCheck(cuMemAlloc(&parts_[a].pAll, P * sizeof(uint64_t)), "cuMemAlloc(pAll)");
            cuCheck(cuMemcpyHtoD(parts_[a].pAll, pPtrs.data(), P * sizeof(uint64_t)),
                    "cuMemcpyHtoD(pAll)");
        }
        finalized_ = true;
    }

    // Everything one GPU's thread does during a solve.
    int runPartition(int m, double atol2, int maxiter, int checkEvery, double timeoutSeconds,
                     bool preprocess, bool direct, SpinBarrier& barrier)
    {
        Partition& p = parts_[m];
        const int P = static_cast<int>(parts_.size());
        // Direct exchange: f_update_p stores the values other GPUs use into their full p.
        direct = direct && P > 1;
        for (const auto& o : parts_)
            if (direct && o.sendMask == 0)
                throw std::runtime_error("pcg_loop: direct exchange requested without send masks");
        CUdeviceptr maskArg = direct ? p.sendMask : 0;       // 0 = no stores to other GPUs
        int offsetArg = static_cast<int>(p.offset);
        CtxGuard g(p.ctx);
        int result = maxiter;

        // Optional: cuSPARSE preprocessing for this matrix, once per solve (the matrix
        // values change every time step). Every SpMV below uses the same parameters
        // and the same, otherwise untouched buffer, as cuSPARSE requires.
        if (preprocess) {
            const double one = 1.0, zero = 0.0;
            spCheck(cusparseSpMV_preprocess(p.sp, CUSPARSE_OPERATION_NON_TRANSPOSE, &one, p.A,
                                            p.vecP, &zero, p.vecQ, CUDA_R_64F,
                                            CUSPARSE_SPMV_ALG_DEFAULT,
                                            reinterpret_cast<void*>(p.spBuffer)),
                    "cusparseSpMV_preprocess");
        }

        // Init: acc = 0;  z = invD r;  own r.r and r.z partials into acc[1]
        cuCheck(cuMemsetD8Async(p.acc, 0, 8 * sizeof(double), p.stream), "cuMemsetD8Async");
        launch(p, p.kInit, p.nblk, {&p.r, &p.z, &p.invD, &p.acc, &p.N});
        record(p.evXRZ, p.stream);

        int Pk = P;
        for (int it = 0; it < maxiter; ++it) {
            int k = it;

            // every GPU's f_update_xrz(k-1) (or f_init) is done -> beta, p = z + beta p
            barrier.wait();
            waitOthers(m, &Partition::evXRZ);
            launch(p, p.kUpdateP, p.nblk,
                   {&p.z, &p.p_m, &p.acc, &p.accAll, &Pk, &p.N, &k, &maskArg, &p.pAll, &offsetArg});
            record(p.evP, p.stream);

            // every GPU's f_update_p is done. Direct exchange: their stores into this GPU's
            // full p are complete, nothing to copy. Otherwise: copy the others' segments.
            barrier.wait();
            if (direct) {
                waitOthers(m, &Partition::evP);
            } else {
                for (int j = 0; j < P; ++j) {
                    if (j == m) continue;
                    Partition& o = parts_[j];
                    cuCheck(cuStreamWaitEvent(p.stream, o.evP, 0), "cuStreamWaitEvent");
                    cuCheck(cuMemcpyPeerAsync(p.p_full + o.offset * sizeof(double), p.ctx,
                                              o.p_m, o.ctx, o.N * sizeof(double), p.stream),
                            "cuMemcpyPeerAsync");
                }
            }

            // q = A p;  own p.q partial
            spmv(p);
            launch(p, p.kDot, p.nblk, {&p.p_m, &p.q, &p.acc, &p.N, &k});
            record(p.evDot, p.stream);

            // every GPU's p.q partial is ready -> alpha; x, r, z; own r.r and r.z partials
            barrier.wait();
            waitOthers(m, &Partition::evDot);
            launch(p, p.kUpdateXRZ, p.nblk,
                   {&p.x, &p.r, &p.z, &p.p_m, &p.q, &p.invD, &p.acc, &p.accAll, &Pk, &p.N, &k});
            record(p.evXRZ, p.stream);

            if ((it + 1) % checkEvery == 0 || it == maxiter - 1) {
                readRR(p, it, timeoutSeconds);        // own partial -> p.hostRR
                barrier.wait();                       // every thread has its partial
                double rr = 0.0;
                for (int j = 0; j < P; ++j) rr += *parts_[j].hostRR;   // same order everywhere
                barrier.wait();                       // all have read before anyone overwrites
                if (rr <= atol2) { result = it + 1; break; }
            }
        }
        waitStream(p, timeoutSeconds);                // leave nothing queued for the caller
        return result;
    }

    void waitOthers(int m, CUevent Partition::*ev)
    {
        Partition& p = parts_[m];
        for (int j = 0; j < static_cast<int>(parts_.size()); ++j)
            if (j != m)
                cuCheck(cuStreamWaitEvent(p.stream, parts_[j].*ev, 0), "cuStreamWaitEvent");
    }

    static void record(CUevent ev, CUstream s) { cuCheck(cuEventRecord(ev, s), "cuEventRecord"); }

    void launch(Partition& p, CUfunction f, int grid, std::vector<void*> args)
    {
        cuCheck(cuLaunchKernel(f, grid, 1, 1, p.tpb, 1, 1, 0, p.stream, args.data(), nullptr),
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
                throw std::runtime_error("pcg_loop: GPU " + std::to_string(p.device) +
                                         " did not finish within the time limit");
        }
    }

    // own r.r partial after iteration it lives in acc[it & 1][ACC_RR]
    void readRR(Partition& p, int it, double timeoutSeconds)
    {
        const CUdeviceptr src = p.acc + ((it & 1) * ACC_STRIDE + ACC_RR) * sizeof(double);
        cuCheck(cuMemcpyDtoHAsync(p.hostRR, src, sizeof(double), p.stream), "cuMemcpyDtoHAsync");
        waitStream(p, timeoutSeconds);
    }
};

PYBIND11_MODULE(pcg_loop, m)
{
    m.doc() = "Compiled host loop for PCGSolverMultiGPU (1 or more GPUs, fused kernels)";
    py::class_<NativePCG>(m, "NativePCG")
        .def(py::init<>())
        .def("add_partition", &NativePCG::addPartition,
             py::arg("device"), py::arg("N"), py::arg("Ntot"), py::arg("nnz"), py::arg("offset"),
             py::arg("nblk"), py::arg("tpb"), py::arg("ptrs"), py::arg("kernels"))
        .def("solve", &NativePCG::solve,
             py::arg("atol2"), py::arg("maxiter"), py::arg("check_every"),
             py::arg("timeout_seconds") = 60.0, py::arg("preprocess") = true,
             py::arg("direct") = false);
}