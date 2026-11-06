#ifndef DYNEARTHSOL3D_BINARYIO_HPP
#define DYNEARTHSOL3D_BINARYIO_HPP

#include "H5Lpublic.h"
#include "H5Gpublic.h"
#include "H5Ppublic.h"
#include "H5Tpublic.h"
#include "H5Fpublic.h"
#include "H5Apublic.h"
#include "H5Spublic.h"
#include "array2d.hpp"

// Layout revision of the frame and checkpoint files: each file's "revision" attribute, which
// a restart checks. Bump it whenever the layout changes.
#define FILE_REVISION 4

struct BuildInfo;
struct CpuInfo;
struct DeviceInfo;

void rename_to_old_backup(const char *filename);


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
    HDF5Output(const char *filename, const int hdf5_compression_level,
               const bool is_chkpt=false, const bool rename_if_exists=false);
    ~HDF5Output();

    template<typename T>
    void write_fieldData(const T& A, const std::string& name);

    template <typename T>
    void write_scalar(const T& A, const std::string& name);

    template <typename T>
    void write_array(const std::vector<T>& A, const char *name, hsize_t len);

    template <typename T>
    void write_aux_array(const std::vector<T>& A, const char *name, hsize_t len);

    template <typename T, int N>
    void write_array(const Array2D<T,N>& A, const char *name, hsize_t len, int dest_N = -1);

    void write_nodal_vec_array(const Array2D<double,NDIMS>& A, const char *name, hsize_t len);

    void write_attribute(const std::string& A, const std::string& name, hid_t& vtkgrpBlock_id);
    template <typename T>
    void write_attribute(const T& A, const std::string& name, hid_t& vtkgrpBlock_id);
    template <typename T>
    void write_attribute(const std::vector<T>& A, const std::string& name, hsize_t len, hid_t& vtkgrpBlock_id);

    void write_block_metadata(const Variables& var, const std::string& base, MarkerSet* ms = nullptr);

    // peak_rss_gib: the caller's running peak, which this write's rss sample may raise.
    void write_run_provenance(const BuildInfo& build, const CpuInfo& cpu,
                              const DeviceInfo& dev, const std::string& restart_from,
                              double& peak_rss_gib);

    void create_virtual_dataset(const std::string& src_name, const std::string& dest_name, hid_t& src_space_id, hid_t& dtype_id);
    void create_virtual_dataset(const std::string& src_name, const std::string& dest_name, hid_t& space_id, hid_t& dtype_id, hsize_t len);
    void create_virtual_dataset(const std::string& src_name, const std::string& dest_name, hid_t& space_id, hid_t& dtype_id, hsize_t len, int N, int dest_N = -1);

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

    bool has_array(const char *name) const;

    template <typename T>
    void read_scalar(T& A, const std::string& name);

    template <typename T>
    void read_array(std::vector<T>& A, const char *name, std::size_t size = 0);

    template <typename T, int N>
    void read_array(Array2D<T,N>& A, const char *name, std::size_t size = 0);
};


#endif
