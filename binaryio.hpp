#ifndef DYNEARTHSOL3D_BINARYIO_HPP
#define DYNEARTHSOL3D_BINARYIO_HPP

#include <map>
#ifdef NETCDF
#include <netcdf>
#endif
#include "array2d.hpp"


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

#ifdef NETCDF

using namespace netCDF;

class NetCDFOutput
{
public:
    NetCDFOutput(const char *filename);
    ~NetCDFOutput();

    void write_header();

    template <typename T>
    void write_array(const std::vector<T>& A, const char *name, std::size_t size);

    template <typename T, int N>
    void write_array(const Array2D<T,N>& A, const char *name, std::size_t size);

private:
    NcFile nc_file;
};

#endif // NETCDF

#endif
