#ifndef DYNEARTHSOL3D_BINARYIO_HPP
#define DYNEARTHSOL3D_BINARYIO_HPP

#include <map>
#ifdef HDF5
#include "H5Lpublic.h"
#include "H5Gpublic.h"
#include "H5Ppublic.h"
#include "H5Tpublic.h"
#include "H5Fpublic.h"
#include "H5Apublic.h"
#include "H5Spublic.h"
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
    void read_array(std::vector<T>& A, const char *name, std::size_t size = 0);

    template <typename T, int N>
    void read_array(Array2D<T,N>& A, const char *name, std::size_t size = 0);
};

#else

class HDF5Output
{
private:
    hid_t file_id = -1;
    const int compression_level;
    long nnode = 0, nelem = 0, nseg = 0, etop = 0, nnode_cell = 0;
    bool has_metadata = false;
    const bool is_checkpoint;
    std::string kind, block_base;

    void write_header();

public:
    HDF5Output(const char *filename, const int hdf5_compression_level, const bool is_chkpt=false);
    ~HDF5Output();

    template<typename T>
    void write_fieldData(const T& A, const std::string& name);

    template <typename T>
    void write_scalar(const T& A, const std::string& name);

    template <typename T>
    void write_array(const std::vector<T>& A, const char *name, hsize_t len);

    template <typename T, int N>
    void write_array(const Array2D<T,N>& A, const char *name, hsize_t len);

    template <typename T>
    void write_attribute(const T& A, const std::string& name, hid_t& vtkgrpBlock_id);
    template <typename T>
    void write_attribute(const std::vector<T>& A, const std::string& name, hsize_t len, hid_t& vtkgrpBlock_id);

    void write_block_metadata(const Variables& var, const std::string& base, MarkerSet* ms = nullptr);

    void create_virtual_dataset(const std::string& src_name, const std::string& dest_name, hid_t& src_space_id, hid_t& dtype_id);
    void create_virtual_dataset(const std::string& src_name, const std::string& dest_name, hid_t& space_id, hid_t& dtype_id, hsize_t len);
    void create_virtual_dataset(const std::string& src_name, const std::string& dest_name, hid_t& space_id, hid_t& dtype_id, hsize_t len, int N);

    void add_soft_link(const std::string& assemblyNodePath, const std::string& linkName,
                        const std::string& targetAbsPath);
    hid_t create_group_with_order(const std::string& path);
};

class HDF5Input
{
private:
    hid_t file_id = -1;

    void read_header();

public:
    HDF5Input(const char *filename);
    ~HDF5Input();

    template <typename T>
    void read_scaler(T& A, const std::string& name);

    template <typename T>
    void read_array(std::vector<T>& A, const char *name, std::size_t size = 0);

    template <typename T, int N>
    void read_array(Array2D<T,N>& A, const char *name, std::size_t size = 0);
};

#endif // HDF5

#endif
