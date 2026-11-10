#ifndef DES_MPI_GPU_BUFFER_HPP
#define DES_MPI_GPU_BUFFER_HPP
#ifdef ACC
#include <cuda_runtime.h>
#include <cstddef>
#include <utility>
#include "utils.hpp"

inline void mpi_cuda_check(cudaError_t status) {
    if(status!=cudaSuccess)die(EXIT_RUNTIME_RESOURCE,cudaGetErrorString(status));
}

// GPU pack/unpack storage paired with pinned HOST storage for ordinary MPI.
// The MPI transport never receives a managed or device pointer in this mode.
// Capacity persists across timesteps; ownership moves when the peer vector grows.
class MpiGpuBuffer {
public:
    MpiGpuBuffer() = default;
    ~MpiGpuBuffer() {release();}
    MpiGpuBuffer(const MpiGpuBuffer&) = delete;
    MpiGpuBuffer& operator=(const MpiGpuBuffer&) = delete;
    MpiGpuBuffer(MpiGpuBuffer&& other) noexcept {swap(other);}
    MpiGpuBuffer& operator=(MpiGpuBuffer&& other) noexcept {
        if(this!=&other){release();swap(other);}return *this;
    }
    void resize(std::size_t n) {
        if(n>capacity_) {
            release();
            mpi_cuda_check(cudaMalloc(reinterpret_cast<void**>(&device_),n*sizeof(double)));
            mpi_cuda_check(cudaMallocHost(reinterpret_cast<void**>(&host_),n*sizeof(double)));
            capacity_=n;
        }
        size_=n;
    }
    double* device() {return device_;}
    double* host() {return host_;}
    std::size_t size() const {return size_;}
    void download() {if(size_)mpi_cuda_check(cudaMemcpy(host_,device_,size_*sizeof(double),cudaMemcpyDeviceToHost));}
    void upload() {if(size_)mpi_cuda_check(cudaMemcpy(device_,host_,size_*sizeof(double),cudaMemcpyHostToDevice));}
private:
    double *device_=nullptr,*host_=nullptr;
    std::size_t size_=0,capacity_=0;
    void release() noexcept {
        if(device_)cudaFree(device_);
        if(host_)cudaFreeHost(host_);
        device_=host_=nullptr;size_=capacity_=0;
    }
    void swap(MpiGpuBuffer& other) noexcept {
        std::swap(device_,other.device_);std::swap(host_,other.host_);
        std::swap(size_,other.size_);std::swap(capacity_,other.capacity_);
    }
};
#endif
#endif
