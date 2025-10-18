#include <cstdio>
#include <cstring>
#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "binaryio.hpp"
#include "markerset.hpp"

#ifdef WIN32
#ifdef _MSC_VER
#define snprintf _snprintf
#endif // _MSC_VER
namespace std { using ::snprintf; }
#endif // WIN32

/*****************************************************************************
 * The format of the binary file:
 * 1  The first 'headerlen' bytes are ASCII text.
 *   1.1  The 1st line in the header is the revision string. Starting with
 *        "# DynEarthSol ndims=%1 revision=%2", with %1 equal to 2 or 3
 *        (indicating 2D or 3D simulation) and %2 an integer.
 *   1.2  The following lines are 'name', 'position' pairs, separated by a
 *        TAB character. This line tells the name of the data and the
 *        starting position (in bytes) of the data in this file.
 * 2  The rests are binary data.
 ****************************************************************************/

namespace {
    const std::size_t headerlen = 4096;
    const char revision_str[] = "# DynEarthSol ndims="
#ifdef THREED
        "3"
#else
        "2"
#endif
        " revision=3\n";
}


/* Not using C++ stream IO for bulk file io since it can be much slower than C stdio. */

#ifndef HDF5

BinaryOutput::BinaryOutput(const char *filename)
{
    f = std::fopen(filename, "wb");
    if (f == NULL) {
        std::cerr << "Error: cannot open file: " << filename << '\n';
        std::exit(2);
    }

    header = new char[headerlen]();
    hd_pos = std::strcat(header, revision_str);
    eof_pos = headerlen;

    std::fseek(f, eof_pos, SEEK_SET);
}


BinaryOutput::~BinaryOutput()
{
    close();
}


void BinaryOutput::close()
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    if (f) {
        /* write header buffer to the beginning of file */
        std::fseek(f, 0, SEEK_SET);
        std::fwrite(header, sizeof(char), headerlen, f);
        std::fclose(f);
        f = NULL;
    }

    delete [] header;
    header = NULL;
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void BinaryOutput::write_header(const char *name)
{
    /* write to header buffer */
    const std::size_t bsize = 256;
    char buffer[bsize];
    std::size_t len = std::snprintf(buffer, bsize, "%s\t%ld\n", name, eof_pos);
    if (len >= bsize) {
        std::cerr << "Error: exceeding buffer length at Output::write_array, name=" << name
                  << " eof_position=" << eof_pos << '\n';
        std::exit(12);
    }
    if (len >= headerlen - (hd_pos - header)*sizeof(char)) {
        std::cerr << "Error: exceeding header length at Output::write_array, name=" << name
                  << " eof_position=" << eof_pos << '\n';
        std::exit(12);
    }
    hd_pos = std::strncat(hd_pos, buffer, len);
}

// XXX: when A is *var.bcflag, i.e. T is uint, g++ cannot instantiate the template
template <typename T>
void BinaryOutput::write_array(const std::vector<T>& A, const char *name, std::size_t size)
{
    write_header(name);
    std::size_t n = std::fwrite(A.data(), sizeof(T), size, f);
    eof_pos += n * sizeof(T);
}


// specialize for uint
void BinaryOutput::write_array(const std::vector<uint>& A, const char *name, std::size_t size)
{
    write_header(name);
    std::size_t n = std::fwrite(A.data(), sizeof(uint), size, f);
    eof_pos += n * sizeof(uint);
}


template <typename T, int N>
void BinaryOutput::write_array(const Array2D<T,N>& A, const char *name, std::size_t size)
{
    write_header(name);
    std::size_t n = std::fwrite(A.data(), sizeof(T), size*N, f);
    eof_pos += n * sizeof(T);
}


// explicit instantiation
template
void BinaryOutput::write_array<int>(const int_vec& A, const char *name, std::size_t);
template
void BinaryOutput::write_array<double>(const double_vec& A, const char *name, std::size_t);

template
void BinaryOutput::write_array<double,NDIMS>(const Array2D<double,NDIMS>& A, const char *name, std::size_t);
template
void BinaryOutput::write_array<double,NSTR>(const Array2D<double,NSTR>& A, const char *name, std::size_t);
#ifdef THREED // when 2d, NSTR == NODES_PER_ELEM == 3
template
void BinaryOutput::write_array<double,NODES_PER_ELEM>(const Array2D<double,NODES_PER_ELEM>& A, const char *name, std::size_t);
#endif
template
void BinaryOutput::write_array<double,1>(const Array2D<double,1>& A, const char *name, std::size_t);
template
void BinaryOutput::write_array<int,NODES_PER_ELEM>(const Array2D<int,NODES_PER_ELEM>& A, const char *name, std::size_t);
template
void BinaryOutput::write_array<int,NDIMS>(const Array2D<int,NDIMS>& A, const char *name, std::size_t);
template
void BinaryOutput::write_array<int,1>(const Array2D<int,1>& A, const char *name, std::size_t);

//////////////////////////////////////////////////////////////////////////////

BinaryInput::BinaryInput(const char *filename)
{
    f = std::fopen(filename, "r");
    if (f == NULL) {
        std::cerr << "Error: cannot open file: " << filename << '\n';
        std::exit(2);
    }
    read_header();
}


BinaryInput::~BinaryInput()
{
    std::fclose(f);
}


void BinaryInput::read_header()
{
    /* Read into header buffer */
    std::fseek(f, 0, SEEK_SET);
    char *header = new char[headerlen]();
    std::size_t n = std::fread(header, sizeof(char), headerlen, f);
    if (n != headerlen) {
        std::cerr << "Error: error reading file header\n";
        std::exit(2);
    }

    /* Parse the content of header buffer */
    char *line = header;

    // Compare revision string (excluding the trailing new line)
    line = std::strtok(header, "\n");
    if (strncmp(line, revision_str, strlen(revision_str)-1) != 0) {
        std::cerr << "Error: mismatching revision string in header\n"
                  << "  Expect: " << revision_str
                  << "  Got: "<< line << '\n';
        std::exit(1);
    }

    line = std::strtok(NULL, "\n");
    while (line != NULL) {
        /* Each line is a string (might contain space), a tab, and an integer */
        char *tab = std::strchr(line, '\t');
        if (tab == NULL) {
            std::cerr << "Error: error parsing file header\n"
                      << " Line is:" << line << '\n';
            std::exit(1);
        }
        std::string name(line, tab-line);
        std::size_t loc;
        std::sscanf(tab, "%zu", &loc);

        offset[name] = loc;
        line = std::strtok(NULL, "\n");
    }

    delete [] header;
}


void BinaryInput::seek_to_array(const char *name)
{
    std::string name2(name);
    auto it = offset.find(name);
    if (it == offset.end()) {
        std::cerr << "Error: no array with a name: " << name << '\n';
        std::exit(1);
    }
    std::size_t loc = it->second;
    //std::cout << name << ' ' << loc << '\n';
    std::fseek(f, loc, SEEK_SET);
}


template <typename T>
void BinaryInput::read_array(std::vector<T>& A, const char *name)
{
    /* The caller must ensure A is of right size to hold the array */

    int size = A.size();
    if (A.size() == 0) {
        std::cerr << "Error: array size is 0: " << name << '\n';
        std::exit(1);
    }

    seek_to_array(name);
    int n = std::fread(A.data(), sizeof(T), size, f);

    if (n != size) {
        std::cerr << "Error: cannot read array: " << name << '\n';
        std::exit(1);
    }
}


template <typename T, int N>
void BinaryInput::read_array(Array2D<T,N>& A, const char *name)
{
    /* The caller must ensure A is of right size to hold the array */

    int size = A.size();
    if (A.size() == 0) {
        std::cerr << "Error: array size is 0: " << name << '\n';
        std::exit(1);
    }

    seek_to_array(name);
    int n = std::fread(A.data(), sizeof(T), size*N, f);

    if (n != N*size) {
        std::cerr << "Error: cannot read array: " << name << '\n';
        std::exit(1);
    }
}


// explicit instantiation
template
void BinaryInput::read_array<double>(double_vec& A, const char *name);
template
void BinaryInput::read_array<int>(int_vec& A, const char *name);
template
void BinaryInput::read_array<double,NDIMS>(Array2D<double,NDIMS>& A, const char *name);
template
void BinaryInput::read_array<double,NSTR>(Array2D<double,NSTR>& A, const char *name);
#ifdef THREED // when 2d, NSTR == NODES_PER_ELEM == 3
template
void BinaryInput::read_array<double,NODES_PER_ELEM>(Array2D<double,NODES_PER_ELEM>& A, const char *name);
#endif
template
void BinaryInput::read_array<double,1>(Array2D<double,1>& A, const char *name);
template
void BinaryInput::read_array<int,NDIMS>(Array2D<int,NDIMS>& A, const char *name);
template
void BinaryInput::read_array<int,NODES_PER_ELEM>(Array2D<int,NODES_PER_ELEM>& A, const char *name);
template
void BinaryInput::read_array<int,1>(Array2D<int,1>& A, const char *name);

#else

HDF5Output::HDF5Output(const char *filename, const int hdf5_compression_level, const bool is_chkpt)
    : h5_file(filename, H5F_ACC_TRUNC), compression_level(hdf5_compression_level), is_chechkpoint(is_chkpt)
{
    write_header();
}

HDF5Output::~HDF5Output()
{
    // nothing to do
}

// Create a group with link creation order tracking (required by VTKHDF Assembly trees)
static Group create_group_with_order(H5File& file, const std::string& path) {
    hid_t fid  = file.getId();
    hid_t gcpl = H5Pcreate(H5P_GROUP_CREATE);
    H5Pset_link_creation_order(gcpl, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);

    hid_t gid  = H5Gcreate2(fid, path.c_str(), H5P_DEFAULT, gcpl, H5P_DEFAULT);
    H5Pclose(gcpl);

    if (gid < 0) throw std::runtime_error(std::string("H5Gcreate2 failed for ") + path);
    H5::Group g(gid);
    H5Gclose(gid);
    return g;
}

// Add a soft link named `linkName` under `assemblyNodePath` that points to `targetAbsPath`
static void add_soft_link(H5File& file,
                          const std::string& assemblyNodePath,
                          const std::string& linkName,
                          const std::string& targetAbsPath) {
    hid_t gid = H5Gopen2(file.getId(), assemblyNodePath.c_str(), H5P_DEFAULT);
    if (gid < 0) throw std::runtime_error(std::string("H5Gopen2 failed for ") + assemblyNodePath);

    herr_t status = H5Lcreate_soft(targetAbsPath.c_str(),
                                   gid,
                                   linkName.c_str(),
                                   H5P_DEFAULT,
                                   H5P_DEFAULT);
    H5Gclose(gid);
    if (status < 0) throw std::runtime_error(std::string("H5Lcreate_soft failed for link '") + linkName + "'");
}

void HDF5Output::write_header()
{
    write_attribute(NDIMS, "ndims", h5_file);
    write_attribute(3, "revision", h5_file);

    create_group_with_order(h5_file, "VTKHDF");
    Group vtkgrp = h5_file.openGroup("/VTKHDF");

    std::string vtkhdf_type = "PartitionedDataSetCollection";
    write_attribute(vtkhdf_type, "Type", vtkgrp);

    int_vec version = {2, 1};
    write_attribute(version, "Version", 2, vtkgrp);

    create_group_with_order(h5_file, "VTKHDF/Assembly");
}

void HDF5Output::write_block_metadata(const Variables& var, const std::string& base, MarkerSet* ms)
{
    int cell_type, link_idx;

    has_metadata = false;
    block_base = base;
    std::string block_path = "/VTKHDF/" + base;

    create_group_with_order(h5_file, block_path);
    Group vtkgrpBlock = h5_file.openGroup(block_path);

    std::string vtkhdf_type = "UnstructuredGrid";
    write_attribute(vtkhdf_type, "Type", vtkgrpBlock);

    int_vec version = {2, 1};
    write_attribute(version, "Version", 2, vtkgrpBlock);
    
    create_group_with_order(h5_file, "/VTKHDF/"+base+"/PointData");
    create_group_with_order(h5_file, "/VTKHDF/"+base+"/CellData");
    create_group_with_order(h5_file, "/VTKHDF/Assembly/"+base);
    add_soft_link(h5_file, "/VTKHDF/Assembly/"+base, base, block_path);

    if (base == "grid") {
        kind = "grid";
        link_idx = 0;
        cell_type = NDIMS == 3 ? 10 : 5; // VTK_TETRA=10, VTK_TRIANGLE=5
        nnode_cell = NODES_PER_ELEM;

        nnode = var.nnode;
        nelem = var.nelem;

        if (!is_chechkpoint) {
            write_array(*var.coord, "Points",  nnode);

            int* conn_ptr = var.connectivity->data();
            int_vec int_tmp(conn_ptr, conn_ptr + nelem*nnode_cell);
            write_array(int_tmp, "Connectivity",  nelem*nnode_cell);
        } else {
            nseg = var.segment->size();
            etop = var.surfinfo.etop;
        }
    } else {
        kind = "marker";
        link_idx = 1;
        cell_type = 1; // VTK_VERTEX=1
        nnode_cell = 1;

        nnode = ms->get_nmarkers();
        nelem = nnode;

        if (!is_chechkpoint) {
            array_t *mcoord = ms->calculate_marker_coord(var); // coordinate of markers
            write_array(*mcoord, "Points",  nnode);
            delete mcoord;

            int_vec int_tmp(nelem);
            for (int i=0; i<nelem; i++) int_tmp[i] = i;
            write_array(int_tmp, "Connectivity",  nelem);
        } else {
            nseg = 0;
            etop = 0;
        }
    }

    if (!is_chechkpoint) {
        int_vec offset(nelem+1);
        for (int i=0; i<nelem+1; ++i) offset[i] = nnode_cell*i;
        write_array(offset, "Offsets",  nelem+1);

        uchar_vec types(nelem, cell_type);
        write_array(types, "Types",  nelem);

        write_scaler(nnode, "NumberOfPoints");
        write_scaler(nelem, "NumberOfCells");
        write_scaler(nelem * nnode_cell, "NumberOfConnectivityIds");
    }
    write_attribute(link_idx, "Index", vtkgrpBlock);
    has_metadata = true;
}

template<typename T>
struct H5TypeMap;

template<>
struct H5TypeMap<int> {
    static PredType type() { return PredType::NATIVE_INT; }
};

template<>
struct H5TypeMap<unsigned int> {
    static PredType type() { return PredType::NATIVE_UINT; }
};

template<>
struct H5TypeMap<long> {
    static PredType type() { return PredType::NATIVE_LONG; }
};

template<>
struct H5TypeMap<float> {
    static PredType type() { return PredType::NATIVE_FLOAT; }
};

template<>
struct H5TypeMap<double> {
    static PredType type() { return PredType::NATIVE_DOUBLE; }
};

template<>
struct H5TypeMap<unsigned char> {
    static PredType type() { return PredType::NATIVE_UCHAR; }
};

template<>
struct H5TypeMap<std::string> {
    static StrType type() { return StrType(PredType::C_S1, H5T_VARIABLE); }
};

// scalear
template<typename T>
void HDF5Output::write_attribute(const T& A, const std::string& name, Group& vtkgrpBlock)
{
    auto dtype = H5TypeMap<T>::type();
    DataSpace attrSpace = DataSpace(H5S_SCALAR);
    Attribute attr = vtkgrpBlock.createAttribute(name, dtype, attrSpace);
    attr.write(dtype, &A);
}

// 1D array
template<typename T>
void HDF5Output::write_attribute(const std::vector<T>& A, const std::string& name, hsize_t len, Group& vtkgrpBlock)
{
    PredType dtype = H5TypeMap<T>::type();
    hsize_t dims[1] = {len};
    DataSpace attrSpace = DataSpace(1, dims);
    Attribute attr = vtkgrpBlock.createAttribute(name, dtype, attrSpace);
    attr.write(dtype, A.data());
}

// explicit instantiation
template
void HDF5Output::write_attribute<int>(const int& A, const std::string& name, Group& vtkgrpBlock);
template
void HDF5Output::write_attribute<double>(const double& A, const std::string& name, Group& vtkgrpBlock);
template
void HDF5Output::write_attribute<uint>(const uint& A, const std::string& name, Group& vtkgrpBlock);
template
void HDF5Output::write_attribute<std::string>(const std::string& A, const std::string& name, Group& vtkgrpBlock);
template
void HDF5Output::write_attribute<int>(const int_vec& A, const std::string& name, hsize_t len, Group& vtkgrpBlock);

// 1D array
template<typename T>
void HDF5Output::write_scaler(const T &A, const std::string& name)
{
    std::string full_name = "/VTKHDF/" + block_base + "/" + name;
    PredType dtype = H5TypeMap<T>::type();

    hsize_t one = 1;
    DataSpace ps = DataSpace(1, &one);
    DataSet ds_np = h5_file.createDataSet(full_name, dtype, ps);

    ds_np.write(&A, dtype);

    if (name == "NumberOfConnectivityIds") return;

    std::string vis_name = name;
    if (kind == "marker") {
        if (name == "NumberOfPoints") {
            vis_name = block_base + ".nmarkers";
        } else if (name == "NumberOfCells") {
            return;
        }
    } else if (kind == "grid") {
        if (name == "NumberOfPoints") {
            vis_name = "nnode";
        } else if (name == "NumberOfCells") {
            vis_name = "nelem";
        }
    }
    create_virtual_dataset(full_name, vis_name, ps, dtype);
}

template
void HDF5Output::write_scaler<int>(const int& A, const std::string& name);
template
void HDF5Output::write_scaler<double>(const double& A, const std::string& name);

// 1D array
template<typename T>
void HDF5Output::write_array(const std::vector<T> &A, const char *name, hsize_t len)
{
    std::string mid;
    if (has_metadata) {
        if (len == nnode) {
            mid = "PointData/";
        } else if (len == nelem) {
            mid = "CellData/";
        } else if (len == nseg || len == etop) {
        } else {
            printf("name = %s\n", name);
            std::exit(13);
        }
    }
    std::string full_name = "/VTKHDF/" + block_base + "/" + mid + name;

    hsize_t dims[1] = { len };
    DataSpace dataspace(1, dims);

    PredType dtype = H5TypeMap<T>::type();

    DSetCreatPropList prop;
    hsize_t chunk_dim = (len < 1024 ? len : 1024);
    prop.setChunk(1, &chunk_dim);
    prop.setDeflate(compression_level);
    prop.setShuffle();

    DataSet dataset = h5_file.createDataSet(full_name, dtype, dataspace, prop);
    dataset.write(A.data(), dtype);

    if (std::string(name) == "Offsets" || std::string(name) == "Types") return;
    if (kind == "marker" && std::string(name) == "Connectivity") return;

    if (std::string(name) == "Connectivity") {
        int len2D = len / nnode_cell;
        create_virtual_dataset(full_name, "connectivity", dataspace, dtype, len2D, nnode_cell);
    } else {
        create_virtual_dataset(full_name, name, dataspace, dtype, len);
    }
}

// 2D array
template<typename T, int N>
void HDF5Output::write_array(const Array2D<T, N>& A, const char *name, hsize_t len)
{
    std::string mid;
    if (has_metadata) {
        if (len == nnode) {
            mid = "PointData/";
        } else if (len == nelem) {
            mid = "CellData/";
        } else if (len == nseg || len == etop) {
        } else {
            printf("name = %s\n", name);
            std::exit(13);
        }
    }

    std::string full_name = "/VTKHDF/" + block_base + "/" + mid + name;

    hsize_t dims[2] = { len, (hsize_t)N };
    DataSpace dataspace(2, dims);

    PredType dtype = H5TypeMap<T>::type();

    DSetCreatPropList prop;
    hsize_t chunk_dims[2];
    chunk_dims[0] = (len < 128 ? len : 128);
    chunk_dims[1] = N;
    prop.setChunk(2, chunk_dims);
    prop.setDeflate(compression_level);
    prop.setShuffle();

    DataSet dataset = h5_file.createDataSet(full_name, dtype, dataspace, prop);
    dataset.write(A.data(), dtype);

    std::string vis_name = name;
    if (std::string(name) == "Points") {
        if (kind == "grid") {
            vis_name = "coordinate";
        } else {
            vis_name = block_base + ".coord";
        }
    }
    create_virtual_dataset(full_name, vis_name, dataspace, dtype, len, N);
}

// explicit instantiation
template
void HDF5Output::write_array<int>(const int_vec& A, const char *name, hsize_t);
template
void HDF5Output::write_array<double>(const double_vec& A, const char *name, hsize_t);
template
void HDF5Output::write_array<uint>(const std::vector<uint>& A, const char *name, hsize_t);
template
void HDF5Output::write_array<unsigned char>(const std::vector<unsigned char>& A, const char *name, hsize_t);
template
void HDF5Output::write_array<double,NDIMS>(const Array2D<double,NDIMS>& A, const char *name, hsize_t);
template
void HDF5Output::write_array<double,NSTR>(const Array2D<double,NSTR>& A, const char *name, hsize_t);
#ifdef THREED // when 2d, NSTR == NODES_PER_ELEM == 3
template
void HDF5Output::write_array<double,NODES_PER_ELEM>(const Array2D<double,NODES_PER_ELEM>& A, const char *name, hsize_t);
#endif
template
void HDF5Output::write_array<double,1>(const Array2D<double,1>& A, const char *name, hsize_t);
template
void HDF5Output::write_array<int,NODES_PER_ELEM>(const Array2D<int,NODES_PER_ELEM>& A, const char *name, hsize_t);
template
void HDF5Output::write_array<int,NDIMS>(const Array2D<int,NDIMS>& A, const char *name, hsize_t);
template
void HDF5Output::write_array<int,1>(const Array2D<int,1>& A, const char *name, hsize_t);

// scaler
void HDF5Output::create_virtual_dataset(const std::string& src_name, const std::string& dest_name, DataSpace& dataspace, PredType& dtype)
{
    hsize_t one = 1;
    hid_t vds_space_id = H5Screate_simple(1, &one, nullptr);

    hid_t vds_dcpl = H5Pcreate(H5P_DATASET_CREATE);
    hid_t src_space_id = H5Scopy(dataspace.getId());
    hid_t file_id = h5_file.getId();
    hid_t dtype_id = dtype.getId();

    herr_t status = H5Pset_virtual(vds_dcpl, vds_space_id, ".", src_name.c_str(), src_space_id);

    hid_t vds_dset_id = H5Dcreate2(file_id, dest_name.c_str(), dtype_id, vds_space_id, H5P_DEFAULT, vds_dcpl, H5P_DEFAULT);

    H5Dclose(vds_dset_id);
    H5Pclose(vds_dcpl);
    H5Sclose(vds_space_id);
    H5Sclose(src_space_id);
}

// 1D array
void HDF5Output::create_virtual_dataset(const std::string& src_name, const std::string& dest_name, DataSpace& dataspace, PredType& dtype, hsize_t len)
{
    hsize_t vds_dims[1] = { len };
    hid_t vds_space_id = H5Screate_simple(1, vds_dims, nullptr);

    hid_t vds_dcpl = H5Pcreate(H5P_DATASET_CREATE);
    hid_t src_space_id = H5Scopy(dataspace.getId());
    hid_t file_id = h5_file.getId();
    hid_t dtype_id = dtype.getId();

    herr_t status = H5Pset_virtual(vds_dcpl, vds_space_id, ".", src_name.c_str(), src_space_id);

    hid_t vds_dset_id = H5Dcreate2(file_id, dest_name.c_str(), dtype_id, vds_space_id, H5P_DEFAULT, vds_dcpl, H5P_DEFAULT);

    H5Dclose(vds_dset_id);
    H5Pclose(vds_dcpl);
    H5Sclose(vds_space_id);
    H5Sclose(src_space_id);
}

// 2D array
void HDF5Output::create_virtual_dataset(const std::string& src_name, const std::string& dest_name, DataSpace& dataspace, PredType& dtype, hsize_t len, int N)
{
    hsize_t vds_dims[2] = { len, static_cast<hsize_t>(N) };
    hid_t vds_space_id = H5Screate_simple(2, vds_dims, nullptr);

    hid_t vds_dcpl = H5Pcreate(H5P_DATASET_CREATE);
    hid_t src_space_id = H5Scopy(dataspace.getId());
    hid_t file_id = h5_file.getId();
    hid_t dtype_id = dtype.getId();

    herr_t status = H5Pset_virtual(vds_dcpl, vds_space_id, ".", src_name.c_str(), src_space_id);

    hid_t vds_dset_id = H5Dcreate2(file_id, dest_name.c_str(), dtype_id, vds_space_id, H5P_DEFAULT, vds_dcpl, H5P_DEFAULT);

    H5Dclose(vds_dset_id);
    H5Pclose(vds_dcpl);
    H5Sclose(vds_space_id);
    H5Sclose(src_space_id);
}

HDF5Input::HDF5Input(const char *filename)
    : h5_file(filename, H5F_ACC_RDONLY)
{
    read_header();
}

void HDF5Input::read_header()
{
    if (!h5_file.attrExists("ndims")) {
        std::cerr << "Error: missing attribute ndims in HDF5 file\n";
        std::exit(1);
    }
    
    Attribute ndims_attr = h5_file.openAttribute("ndims");
    DataType ndims_type = ndims_attr.getDataType();

    int ndims = -1;
    ndims_attr.read(ndims_type, &ndims);
    // std::cout << "ndims = " << ndims << std::endl;

    if (!h5_file.attrExists("revision")) {
        std::cerr << "Error: missing attribute revision in HDF5 file\n";
        std::exit(1);
    }

    Attribute rev_attr = h5_file.openAttribute("revision");
    DataType rev_type = rev_attr.getDataType();

    int rev = -1;
    rev_attr.read(rev_type, &rev);
    // std::cout << "revision = " << rev << std::endl;
}

HDF5Input::~HDF5Input()
{
    // nothing to do
}

template <typename T>
void HDF5Input::read_scaler(T& A, const std::string& name)
{
    DataSet dataset;
    try {
        dataset = h5_file.openDataSet(name);
    }
    catch (FileIException &e) {
        std::cerr << "Error: cannot open dataset (file error): " << name << '\n';
        e.printErrorStack();
        std::exit(1);
    }
    catch (DataSetIException &e) {
        std::cerr << "Error: cannot open dataset (dataset error): " << name << '\n';
        e.printErrorStack();
        std::exit(1);
    }

    DataSpace filespace = dataset.getSpace();
    int rank = filespace.getSimpleExtentNdims();
    if (rank != 1) {
        std::cerr << "Error: dataset rank mismatch for " << name
                  << ", expected 1 dim, got " << rank << "\n";
        std::exit(1);
    }

    hsize_t one = 1;
    filespace.getSimpleExtentDims(&one, nullptr);

    DataSpace memspace(1, &one);
    PredType dtype = H5TypeMap<T>::type();

    try {
        dataset.read(&A, dtype, memspace, filespace);
    }
    catch (DataSetIException &e) {
        std::cerr << "Error: failed to read dataset: " << name << "\n";
        e.printErrorStack();
        std::exit(1);
    }
    catch (DataSpaceIException &e) {
        std::cerr << "Error: dataspace error reading dataset: " << name << "\n";
        e.printErrorStack();
        std::exit(1);
    }
}

template
void HDF5Input::read_scaler<int>(int& A, const std::string& name);
template
void HDF5Input::read_scaler<double>(double& A, const std::string& name);

template <typename T>
void HDF5Input::read_array(std::vector<T>& A, const char *name)
{
    std::size_t size = A.size();
    if (size == 0) {
        std::cerr << "Error: array size is 0: " << name << '\n';
        std::exit(1);
    }

    DataSet dataset;
    try {
        dataset = h5_file.openDataSet(name);
    }
    catch (FileIException &e) {
        std::cerr << "Error: cannot open dataset (file error): " << name << '\n';
        e.printErrorStack();
        std::exit(1);
    }
    catch (DataSetIException &e) {
        std::cerr << "Error: cannot open dataset (dataset error): " << name << '\n';
        e.printErrorStack();
        std::exit(1);
    }

    DataSpace filespace = dataset.getSpace();
    int rank = filespace.getSimpleExtentNdims();
    if (rank != 1) {
        std::cerr << "Error: dataset rank mismatch for " << name
                  << ", expected 1 dim, got " << rank << "\n";
        std::exit(1);
    }
    hsize_t dims[1];
    filespace.getSimpleExtentDims(dims, nullptr);
    if (dims[0] != size) {
        std::cerr << "Error: array size is not matched: " << name
                  << " (file dim = " << dims[0] << ", expected = " << size << ")\n";
        std::exit(1);
    }

    DataSpace memspace(1, dims);
    PredType dtype = H5TypeMap<T>::type();

    try {
        dataset.read(A.data(), dtype, memspace, filespace);
    }
    catch (DataSetIException &e) {
        std::cerr << "Error: failed to read dataset: " << name << "\n";
        e.printErrorStack();
        std::exit(1);
    }
    catch (DataSpaceIException &e) {
        std::cerr << "Error: dataspace error reading dataset: " << name << "\n";
        e.printErrorStack();
        std::exit(1);
    }
}


template <typename T, int N>
void HDF5Input::read_array(Array2D<T,N>& A, const char *name)
{
    /* The caller must ensure A is of right size to hold the array */

    int size = A.size();
    if (A.size() == 0) {
        std::cerr << "Error: array size is 0: " << name << '\n';
        std::exit(1);
    }

    DataSet dataset;
    try {
        dataset = h5_file.openDataSet(name);
    } catch (FileIException &e) {
        std::cerr << "Error: dataset not found: " << name << "\n";
        e.printErrorStack();
        std::exit(1);
    } catch (DataSetIException &e) {
        std::cerr << "Error: cannot open dataset: " << name << "\n";
        e.printErrorStack();
        std::exit(1);
    }

    DataSpace filespace = dataset.getSpace();
    int rank = filespace.getSimpleExtentNdims();
    if (rank != 2) {
        std::cerr << "Error: dataset rank mismatch for " << name 
                  << ", expected 2 dims, got " << rank << '\n';
        std::exit(1);
    }
    
    hsize_t dims[2];
    filespace.getSimpleExtentDims(dims, nullptr);
    if (dims[0] != size || dims[1] != static_cast<hsize_t>(N)) {
        std::cerr << "Error: dataset dimensions mismatch for " << name
                  << ": file dims = (" << dims[0] << ", " << dims[1]
                  << "), expected (" << size << ", " << N << ")\n";
        std::exit(1);
    }

    DataSpace memspace(2, dims);
    PredType dtype = H5TypeMap<T>::type();

    try {
        dataset.read(A.data(), dtype, memspace, filespace);
    } catch (DataSetIException &e) {
        std::cerr << "Error: failed to read dataset: " << name << "\n";
        e.printErrorStack();
        std::exit(1);
    } catch (DataSpaceIException &e) {
        std::cerr << "Error: dataspace error reading dataset: " << name << "\n";
        e.printErrorStack();
        std::exit(1);
    }
}

// explicit instantiation
template
void HDF5Input::read_array<double>(double_vec& A, const char *name);
template
void HDF5Input::read_array<int>(int_vec& A, const char *name);
template
void HDF5Input::read_array<double,NDIMS>(Array2D<double,NDIMS>& A, const char *name);
template
void HDF5Input::read_array<double,NSTR>(Array2D<double,NSTR>& A, const char *name);
#ifdef THREED // when 2d, NSTR == NODES_PER_ELEM == 3
template
void HDF5Input::read_array<double,NODES_PER_ELEM>(Array2D<double,NODES_PER_ELEM>& A, const char *name);
#endif
template
void HDF5Input::read_array<double,1>(Array2D<double,1>& A, const char *name);
template
void HDF5Input::read_array<int,NDIMS>(Array2D<int,NDIMS>& A, const char *name);
template
void HDF5Input::read_array<int,NODES_PER_ELEM>(Array2D<int,NODES_PER_ELEM>& A, const char *name);
template
void HDF5Input::read_array<int,1>(Array2D<int,1>& A, const char *name);

#endif