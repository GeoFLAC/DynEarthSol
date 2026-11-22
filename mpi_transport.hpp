#ifndef DES_MPI_TRANSPORT_HPP
#define DES_MPI_TRANSPORT_HPP
#include <mpi.h>
#include <cstddef>
#include <new>
#include <vector>
#ifdef ACC
#include <cuda_runtime.h>
#endif

// Vectors that MPI reads or writes. A GPU build (-gpu=mem:managed) puts every heap
// allocation in CUDA managed memory, and a CUDA-aware transport (UCX) corrupted managed
// payloads on two GPUs. It treats pinned host memory as host memory, so a GPU build pins
// these vectors, as it pins the halo staging buffers.
#ifdef ACC
template<class T> struct MpiHostAllocator {
    using value_type = T;
    MpiHostAllocator() = default;
    template<class U> MpiHostAllocator(const MpiHostAllocator<U>&) {}
    T* allocate(std::size_t n) {
        void* pointer = nullptr;
        if (cudaMallocHost(&pointer, n * sizeof(T)) != cudaSuccess) {
            cudaGetLastError();
            throw std::bad_alloc();
        }
        return static_cast<T*>(pointer);
    }
    void deallocate(T* pointer, std::size_t) { cudaFreeHost(pointer); }
};
template<class T, class U> bool operator==(const MpiHostAllocator<T>&, const MpiHostAllocator<U>&) { return true; }
template<class T, class U> bool operator!=(const MpiHostAllocator<T>&, const MpiHostAllocator<U>&) { return false; }
template<class T> using mpi_host_vector = std::vector<T, MpiHostAllocator<T>>;
#else
template<class T> using mpi_host_vector = std::vector<T>;
#endif

// Contiguous buffers only. The transport never dereferences pointers: a GPU
// backend can pass device buffers when MPI is CUDA-aware. Producers must finish
// before this call; completion makes receives safe for consumers and buffers reusable.
struct MpiPeerBuffer {
    int rank;
    const double* send;
    double* receive;
    int send_count, receive_count;
};
inline int mpi_post_buffers(const MpiPeerBuffer* peers,int n,int tag,MPI_Request* requests) {
    int pending=0;
    for(int i=0;i<n;++i) if(peers[i].receive_count)
        MPI_Irecv(peers[i].receive,peers[i].receive_count,MPI_DOUBLE,peers[i].rank,tag,MPI_COMM_WORLD,&requests[pending++]);
    for(int i=0;i<n;++i) if(peers[i].send_count)
        MPI_Isend(peers[i].send,peers[i].send_count,MPI_DOUBLE,peers[i].rank,tag,MPI_COMM_WORLD,&requests[pending++]);
    return pending;
}
inline void mpi_finish_buffers(int pending,MPI_Request* requests) {
    if(pending) MPI_Waitall(pending,requests,MPI_STATUSES_IGNORE);
}
#endif
