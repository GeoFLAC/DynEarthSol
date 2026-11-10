#ifndef DES_MPI_HALO_KERNELS_HPP
#define DES_MPI_HALO_KERNELS_HPP
#include <cstddef>

// The production CPU packers and the accelerator harness share this indexing.
// With DES_MPI_DEVICE_KERNELS, all three pointers must be device addresses.
// Kernels complete before return, so producers are finished before MPI reads.
inline void mpi_pack_columns(const double* field,int stride,int columns,
                             const int* rows,int nrows,double* packed) {
#ifdef DES_MPI_DEVICE_KERNELS
    #pragma acc parallel loop collapse(2) deviceptr(field,rows,packed)
#endif
    for(int c=0;c<columns;++c)for(int i=0;i<nrows;++i)
        packed[std::size_t(c)*nrows+i]=field[std::size_t(c)*stride+rows[i]];
}
inline void mpi_unpack_columns(double* field,int stride,int columns,
                               const int* rows,int nrows,const double* packed) {
#ifdef DES_MPI_DEVICE_KERNELS
    #pragma acc parallel loop collapse(2) deviceptr(field,rows,packed)
#endif
    for(int c=0;c<columns;++c)for(int i=0;i<nrows;++i)
        field[std::size_t(c)*stride+rows[i]]=packed[std::size_t(c)*nrows+i];
}

// Explicit device variants let production management keep using host packers.
// Field/row pointers may be CUDA managed addresses; packed is a device allocation.
inline void mpi_pack_device_columns(const double* field,int stride,int columns,
                                    const int* rows,int nrows,double* packed) {
    #pragma acc parallel loop collapse(2) deviceptr(field,rows,packed)
    for(int c=0;c<columns;++c)for(int i=0;i<nrows;++i)
        packed[std::size_t(c)*nrows+i]=field[std::size_t(c)*stride+rows[i]];
}
inline void mpi_unpack_device_columns(double* field,int stride,int columns,
                                      const int* rows,int nrows,const double* packed) {
    #pragma acc parallel loop collapse(2) deviceptr(field,rows,packed)
    for(int c=0;c<columns;++c)for(int i=0;i<nrows;++i)
        field[std::size_t(c)*stride+rows[i]]=packed[std::size_t(c)*nrows+i];
}
#endif
