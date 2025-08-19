#pragma once
#include "kernel.h"

// Simple RAII wrapper for CUDA device memory
// Allocates on construction and frees on destruction.
template <typename T>
class CudaMemory {
    T* ptr_ = nullptr;

   public:
    CudaMemory() = default;
    explicit CudaMemory(size_t count) { allocate(count); }
    CudaMemory(const CudaMemory&) = delete;
    CudaMemory& operator=(const CudaMemory&) = delete;
    CudaMemory(CudaMemory&& other) noexcept : ptr_(other.ptr_) { other.ptr_ = nullptr; }
    CudaMemory& operator=(CudaMemory&& other) noexcept {
        if (this != &other) {
            reset_noexcept();
            ptr_ = other.ptr_;
            other.ptr_ = nullptr;
        }
        return *this;
    }
    ~CudaMemory() noexcept { reset_noexcept(); }

    void allocate(size_t count) {
        reset_noexcept();
        checkCUDA(cudaMalloc(&ptr_, count * sizeof(T)));
    }

    T* get() const { return ptr_; }
    operator T*() const { return ptr_; }

    void reset() { reset_noexcept(); }

   private:
    void reset_noexcept() noexcept {
        if (ptr_) {
            cudaFree(ptr_); // ignore error in noexcept context
            ptr_ = nullptr;
        }
    }
};
