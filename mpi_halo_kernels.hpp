#ifndef DES_MPI_HALO_KERNELS_HPP
#define DES_MPI_HALO_KERNELS_HPP
#include <cstddef>

// Host management, GPU halos and the accelerator harness share this indexing.
// With OnDevice=true, all three pointers must be device addresses.
// Kernels complete before return, so producers are finished before MPI reads.
template<bool OnDevice=false>
inline void mpi_pack_columns(const double* field,int stride,int columns,
                             const int* rows,int nrows,double* packed) {
    #pragma acc parallel loop collapse(2) if(OnDevice) deviceptr(field,rows,packed)
    for(int c=0;c<columns;++c)for(int i=0;i<nrows;++i)
        packed[std::size_t(c)*nrows+i]=field[std::size_t(c)*stride+rows[i]];
}
template<bool OnDevice=false>
inline void mpi_unpack_columns(double* field,int stride,int columns,
                               const int* rows,int nrows,const double* packed) {
    #pragma acc parallel loop collapse(2) if(OnDevice) deviceptr(field,rows,packed)
    for(int c=0;c<columns;++c)for(int i=0;i<nrows;++i)
        field[std::size_t(c)*stride+rows[i]]=packed[std::size_t(c)*nrows+i];
}

#endif
