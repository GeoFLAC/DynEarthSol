#ifndef DYNEARTHSOL3D_BINARYIO_HPP
#define DYNEARTHSOL3D_BINARYIO_HPP

#include <map>
#ifdef HDF5
#include <H5Cpp.h>
#endif
#include "array2d.hpp"

#ifndef HDF5

class BinaryOutput
{
private:
    long eof_pos;
    char *header;
    char *hd_pos;
    std::FILE* f;

    void write_header(const char *name);

public:
    BinaryOutput(const char *filename);
    ~BinaryOutput();

    void close();

    template <typename T>
    void write_array(const std::vector<T>& A, const char *name, std::size_t size);

    void write_array(const std::vector<uint>& A, const char *name, std::size_t size);

    template <typename T, int N>
    void write_array(const Array2D<T,N>& A, const char *name, std::size_t size);
};


class BinaryInput
{
private:
    std::FILE* f;
    std::map<std::string, std::size_t> offset;

    void read_header();
    void seek_to_array(const char *name);

public:
    BinaryInput(const char *filename);
    ~BinaryInput();

    template <typename T>
    void read_array(std::vector<T>& A, const char *name);

    template <typename T, int N>
    void read_array(Array2D<T,N>& A, const char *name);
};

#else

using namespace H5;

class HDF5Output
{
private:
    H5File h5_file;
    const int compression_level;

    void write_header();

public:
    HDF5Output(const char *filename, const int hdf5_compression_level);
    ~HDF5Output();

    template <typename T>
    void write_array(const std::vector<T>& A, const char *name, hsize_t len);

    template <typename T, int N>
    void write_array(const Array2D<T,N>& A, const char *name, hsize_t len);
};

class HDF5Input
{
private:
    H5File h5_file;

    void read_header();

public:
    HDF5Input(const char *filename);
    ~HDF5Input();

    template <typename T>
    void read_array(std::vector<T>& A, const char *name);

    template <typename T, int N>
    void read_array(Array2D<T,N>& A, const char *name);
};

#endif // HDF5

#endif
