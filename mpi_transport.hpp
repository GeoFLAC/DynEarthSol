#ifndef DES_MPI_TRANSPORT_HPP
#define DES_MPI_TRANSPORT_HPP
#include <mpi.h>

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
    MPI_Waitall(pending,requests,MPI_STATUSES_IGNORE);
}
inline void mpi_exchange_buffers(const MpiPeerBuffer* peers,int n,int tag,MPI_Request* requests) {
    mpi_finish_buffers(mpi_post_buffers(peers,n,tag,requests),requests);
}
#endif
