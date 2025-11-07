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
void BinaryInput::read_array(std::vector<T>& A, const char *name, std::size_t size)
{
    /* The caller must ensure A is of right size to hold the array */

    size = size > 0 ? size : A.size();
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
void BinaryInput::read_array(Array2D<T,N>& A, const char *name, std::size_t size)
{
    /* The caller must ensure A is of right size to hold the array */

    size = size > 0 ? size : A.size();
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
void BinaryInput::read_array<double>(double_vec& A, const char *name, std::size_t size);
template
void BinaryInput::read_array<int>(int_vec& A, const char *name, std::size_t size);
template
void BinaryInput::read_array<double,NDIMS>(Array2D<double,NDIMS>& A, const char *name, std::size_t size);
template
void BinaryInput::read_array<double,NSTR>(Array2D<double,NSTR>& A, const char *name, std::size_t size);
#ifdef THREED // when 2d, NSTR == NODES_PER_ELEM == 3
template
void BinaryInput::read_array<double,NODES_PER_ELEM>(Array2D<double,NODES_PER_ELEM>& A, const char *name, std::size_t size);
#endif
template
void BinaryInput::read_array<double,1>(Array2D<double,1>& A, const char *name, std::size_t size);
template
void BinaryInput::read_array<int,NDIMS>(Array2D<int,NDIMS>& A, const char *name, std::size_t size);
template
void BinaryInput::read_array<int,NODES_PER_ELEM>(Array2D<int,NODES_PER_ELEM>& A, const char *name, std::size_t size);
template
void BinaryInput::read_array<int,1>(Array2D<int,1>& A, const char *name, std::size_t size);

#else

HDF5Output::HDF5Output(const char *filename, const int hdf5_compression_level, const bool is_chkpt)
    : compression_level(hdf5_compression_level), is_checkpoint(is_chkpt)
{
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    H5Pclose(fapl_id);

    if (file_id < 0) {
        throw std::runtime_error(std::string("H5Fcreate failed: ") + filename);
    }

    write_header();
}

HDF5Output::~HDF5Output()
{
    if (file_id >= 0) {
        H5Fflush(file_id, H5F_SCOPE_GLOBAL);
        H5Fclose(file_id);
        file_id = -1;
    }
}

// Create a group with link creation order tracking (required by VTKHDF Assembly trees)
hid_t HDF5Output::create_group_with_order(const std::string& path) {
    hid_t gcpl_id = H5Pcreate(H5P_GROUP_CREATE);
    H5Pset_link_creation_order(gcpl_id, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);

    hid_t gid  = H5Gcreate2(file_id, path.c_str(), H5P_DEFAULT, gcpl_id, H5P_DEFAULT);
    H5Pclose(gcpl_id);

    if (gid < 0) throw std::runtime_error(std::string("H5Gcreate2 failed for ") + path);
    return gid;
}

// Add a soft link named `linkName` under `assemblyNodePath` that points to `targetAbsPath`
void HDF5Output::add_soft_link(const std::string& assemblyNodePath,
                          const std::string& linkName,
                          const std::string& targetAbsPath) {
    hid_t gid = H5Gopen2(file_id, assemblyNodePath.c_str(), H5P_DEFAULT);
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
    write_attribute(NDIMS, "ndims", file_id);
    write_attribute(3, "revision", file_id);

    hid_t gid = create_group_with_order("/VTKHDF");

    std::string vtkhdf_type = "PartitionedDataSetCollection";
    write_attribute(vtkhdf_type, "Type", gid);

    int_vec version = {2, 1};
    write_attribute(version, "Version", 2, gid);
    H5Gclose(gid);

    gid = create_group_with_order("/VTKHDF/Assembly");
    H5Gclose(gid);
}

void HDF5Output::write_block_metadata(const Variables& var, const std::string& base, MarkerSet* ms)
{
    int cell_type, link_idx;

    has_metadata = false;
    block_base = base;
    std::string block_path = "/VTKHDF/" + base;

    hid_t gid_block = create_group_with_order(block_path);

    std::string vtkhdf_type = "UnstructuredGrid";
    write_attribute(vtkhdf_type, "Type", gid_block);

    int_vec version = {2, 1};
    write_attribute(version, "Version", 2, gid_block);
    
    hid_t gid = create_group_with_order("/VTKHDF/"+base+"/PointData");
    H5Gclose(gid);
    gid = create_group_with_order("/VTKHDF/"+base+"/CellData");
    H5Gclose(gid);
    gid = create_group_with_order("/VTKHDF/Assembly/"+base);
    H5Gclose(gid);
    add_soft_link("/VTKHDF/Assembly/"+base, base, block_path);

    if (base == "grid") {
        gid = create_group_with_order("/VTKHDF/"+base+"/FieldData");
        H5Gclose(gid);

        kind = "grid";
        link_idx = 0;
        cell_type = NDIMS == 3 ? 10 : 5; // VTK_TETRA=10, VTK_TRIANGLE=5
        nnode_cell = NODES_PER_ELEM;

        nnode = var.nnode;
        nelem = var.nelem;

        if (!is_checkpoint) {
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

        if (!is_checkpoint) {
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

    if (!is_checkpoint) {
        int_vec offset(nelem+1);
        for (int i=0; i<nelem+1; ++i) offset[i] = nnode_cell*i;
        write_array(offset, "Offsets",  nelem+1);

        uchar_vec types(nelem, cell_type);
        write_array(types, "Types",  nelem);

        write_scalar(nnode, "NumberOfPoints");
        write_scalar(nelem, "NumberOfCells");
        write_scalar(nelem * nnode_cell, "NumberOfConnectivityIds");
    }
    write_attribute(link_idx, "Index", gid_block);
    H5Gclose(gid_block);
    has_metadata = true;
}

template<typename T> struct H5Native;
template<> struct H5Native<int>            { static hid_t id() { return H5T_NATIVE_INT; } };
template<> struct H5Native<unsigned int>   { static hid_t id() { return H5T_NATIVE_UINT; } };
template<> struct H5Native<long>           { static hid_t id() { return H5T_NATIVE_LONG; } };
template<> struct H5Native<float>          { static hid_t id() { return H5T_NATIVE_FLOAT; } };
template<> struct H5Native<double>         { static hid_t id() { return H5T_NATIVE_DOUBLE; } };
template<> struct H5Native<unsigned char>  { static hid_t id() { return H5T_NATIVE_UCHAR; } };
template<> struct H5Native<std::string>    { static hid_t id() { return H5T_C_S1; } };

template<typename T>
void HDF5Output::write_fieldData(const T& A, const std::string& name)
{
    std::string full_name = "/VTKHDF/" + block_base + "/FieldData/" + name;
    hid_t dtype_id = H5Native<T>::id();

    hsize_t one = 1;
    hid_t space_id = H5Screate_simple(1, &one, nullptr);
    hid_t dset_id = H5Dcreate2(file_id, full_name.c_str(), dtype_id, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &A);

    create_virtual_dataset(full_name, name, space_id, dtype_id);

    H5Dclose(dset_id);
    H5Sclose(space_id);
}

template
void HDF5Output::write_fieldData<int>(const int& A, const std::string& name);
template
void HDF5Output::write_fieldData<double>(const double& A, const std::string& name);
template
void HDF5Output::write_fieldData<long>(const long& A, const std::string& name);


// scalear
template<typename T>
void HDF5Output::write_attribute(const T& A, const std::string& name, hid_t& vtkgrpBlock_id)
{
    hid_t dtype_id = H5Native<T>::id();
    hid_t space_id = H5Screate(H5S_SCALAR);
    hid_t attr_id = H5Acreate2(vtkgrpBlock_id, name.c_str(), dtype_id, space_id, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(attr_id, dtype_id, &A);

	H5Aclose(attr_id);
	H5Sclose(space_id);
}

void HDF5Output::write_attribute(const std::string& A, const std::string& name, hid_t& vtkgrpBlock_id)
{
    // for string type, create a copy to avoid closing H5T_C_S1 later
    hid_t str_t = H5Tcopy(H5T_C_S1);
    H5Tset_size(str_t, H5T_VARIABLE);

    hid_t space_id = H5Screate(H5S_SCALAR);
    hid_t attr_id = H5Acreate2(vtkgrpBlock_id, name.c_str(), str_t, space_id, H5P_DEFAULT, H5P_DEFAULT);
    
    // For variable-length strings, HDF5 expects a pointer to a C string (const char*)
    const char* c_str = A.c_str();
    H5Awrite(attr_id, str_t, &c_str);
    
    H5Aclose(attr_id);
    H5Sclose(space_id);
    H5Tclose(str_t);
}

// 1D array
template<typename T>
void HDF5Output::write_attribute(const std::vector<T>& A, const std::string& name, hsize_t len, hid_t& vtkgrpBlock_id)
{
    hid_t dtype_id = H5Native<T>::id();
    hid_t space_id = H5Screate_simple(1, &len, nullptr);
    hid_t attr_id = H5Acreate2(vtkgrpBlock_id, name.c_str(), dtype_id, space_id, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(attr_id, dtype_id, A.data());

	H5Aclose(attr_id);
	H5Sclose(space_id);
}

// explicit instantiation
template void HDF5Output::write_attribute<int>(const int& A, const std::string& name, hid_t& vtkgrpBlock_id);
template void HDF5Output::write_attribute<double>(const double& A, const std::string& name, hid_t& vtkgrpBlock_id);
template void HDF5Output::write_attribute<uint>(const uint& A, const std::string& name, hid_t& vtkgrpBlock_id);
template void HDF5Output::write_attribute<int>(const int_vec& A, const std::string& name, hsize_t len, hid_t& vtkgrpBlock_id);

// 1D array
template<typename T>
void HDF5Output::write_scalar(const T &A, const std::string& name)
{
    std::string full_name = "/VTKHDF/" + block_base + "/" + name;
    hid_t dtype_id = H5Native<T>::id();

    hsize_t one = 1;
    hid_t space_id = H5Screate_simple(1, &one, nullptr);
    hid_t dset_id = H5Dcreate2(file_id, full_name.c_str(), dtype_id, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &A);

    if (name == "NumberOfConnectivityIds") return;

    std::string vis_name = name;
    bool is_field = false;
    if (kind == "marker") {
        if (name == "NumberOfPoints") {
            vis_name = block_base + ".nmarkers";
            is_field = true;
        } else if (name == "NumberOfCells") {
            return;
        }
    } else if (kind == "grid") {
        if (name == "NumberOfPoints") {
            vis_name = "nnode";
            is_field = true;
        } else if (name == "NumberOfCells") {
            vis_name = "nelem";
            is_field = true;
        }
    }
    create_virtual_dataset(full_name, vis_name, space_id, dtype_id);
    if (is_field) {
        vis_name = "/VTKHDF/grid/FieldData/" + vis_name;
        create_virtual_dataset(full_name, vis_name, space_id, dtype_id);
    }
    H5Dclose(dset_id);
    H5Sclose(space_id);
}

template
void HDF5Output::write_scalar<int>(const int& A, const std::string& name);
template
void HDF5Output::write_scalar<double>(const double& A, const std::string& name);

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

    hid_t space_id = H5Screate_simple(1, &len, nullptr);

    hid_t dtype_id = H5Native<T>::id();

    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk_dim = (len < 1024 ? len : 1024);
    H5Pset_chunk(dcpl_id, 1, &chunk_dim);
    H5Pset_shuffle(dcpl_id);
    H5Pset_deflate(dcpl_id, compression_level);

    hid_t dset_id = H5Dcreate2(file_id, full_name.c_str(), dtype_id, space_id,
                        H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, A.data());

    if (std::string(name) == "Offsets" || std::string(name) == "Types") return;
    if (kind == "marker" && std::string(name) == "Connectivity") return;

    if (std::string(name) == "Connectivity") {
        int len2D = len / nnode_cell;
        create_virtual_dataset(full_name, "connectivity", space_id, dtype_id, len2D, nnode_cell);
    } else {
        create_virtual_dataset(full_name, name, space_id, dtype_id, len);
    }
    H5Dclose(dset_id);
    H5Pclose(dcpl_id);
    H5Sclose(space_id);
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
    hid_t space_id = H5Screate_simple(2, dims, nullptr);

    hid_t dtype_id = H5Native<T>::id();

    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk_dims[2];
    chunk_dims[0] = (len < 128 ? len : 128);
    chunk_dims[1] = N;
    H5Pset_chunk(dcpl_id, 2, chunk_dims);
    H5Pset_shuffle(dcpl_id);
    H5Pset_deflate(dcpl_id, compression_level);

    hid_t dset_id = H5Dcreate2(file_id, full_name.c_str(), dtype_id, space_id,
                        H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, A.data());

    std::string vis_name = name;
    if (std::string(name) == "Points") {
        if (kind == "grid") {
            vis_name = "coordinate";
        } else {
            vis_name = block_base + ".coord";
        }
    }
    create_virtual_dataset(full_name, vis_name, space_id, dtype_id, len, N);

    H5Dclose(dset_id);
    H5Pclose(dcpl_id);
    H5Sclose(space_id);
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
void HDF5Output::create_virtual_dataset(const std::string& src_name, const std::string& dest_name, hid_t& src_space_id, hid_t& dtype_id)
{
    hsize_t one = 1;
    hid_t vds_space_id = H5Screate_simple(1, &one, nullptr);
    hid_t vds_dcpl = H5Pcreate(H5P_DATASET_CREATE);
    // hid_t file_id = h5_file.getId();

    herr_t status = H5Pset_virtual(vds_dcpl, vds_space_id, ".", src_name.c_str(), src_space_id);

    hid_t vds_dset_id = H5Dcreate2(file_id, dest_name.c_str(), dtype_id, vds_space_id, H5P_DEFAULT, vds_dcpl, H5P_DEFAULT);

    H5Dclose(vds_dset_id);
    H5Pclose(vds_dcpl);
    H5Sclose(vds_space_id);
}

// 1D array
void HDF5Output::create_virtual_dataset(const std::string& src_name, const std::string& dest_name, hid_t& space_id, hid_t& dtype_id, hsize_t len)
{
    hsize_t vds_dims[1] = { len };
    hid_t vds_space_id = H5Screate_simple(1, vds_dims, nullptr);

    hid_t vds_dcpl = H5Pcreate(H5P_DATASET_CREATE);
    hid_t src_space_id = H5Scopy(space_id);

    herr_t status = H5Pset_virtual(vds_dcpl, vds_space_id, ".", src_name.c_str(), src_space_id);

    hid_t vds_dset_id = H5Dcreate2(file_id, dest_name.c_str(), dtype_id, vds_space_id, H5P_DEFAULT, vds_dcpl, H5P_DEFAULT);

    H5Dclose(vds_dset_id);
    H5Pclose(vds_dcpl);
    H5Sclose(vds_space_id);
    H5Sclose(src_space_id);
}

// 2D array
void HDF5Output::create_virtual_dataset(const std::string& src_name, const std::string& dest_name, hid_t& space_id, hid_t& dtype_id, hsize_t len, int N)
{
    hsize_t vds_dims[2] = { len, static_cast<hsize_t>(N) };
    hid_t vds_space_id = H5Screate_simple(2, vds_dims, nullptr);

    hid_t vds_dcpl = H5Pcreate(H5P_DATASET_CREATE);
    hid_t src_space_id = H5Scopy(space_id);

    herr_t status = H5Pset_virtual(vds_dcpl, vds_space_id, ".", src_name.c_str(), src_space_id);

    hid_t vds_dset_id = H5Dcreate2(file_id, dest_name.c_str(), dtype_id, vds_space_id, H5P_DEFAULT, vds_dcpl, H5P_DEFAULT);

    H5Dclose(vds_dset_id);
    H5Pclose(vds_dcpl);
    H5Sclose(vds_space_id);
    H5Sclose(src_space_id);
}

HDF5Input::HDF5Input(const char *filename)
{
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        std::cerr << "Error: cannot open HDF5 file for reading: " << filename << "\n";
        std::exit(1);
    }

    read_header();
}

void HDF5Input::read_header()
{
    if (H5Aexists(file_id, "ndims") <= 0) {
        std::cerr << "Error: missing attribute ndims in HDF5 file\n";
        std::exit(1);
    }

    hid_t attr = H5Aopen(file_id, "ndims", H5P_DEFAULT);
    hid_t atype = H5Aget_type(attr);

    int ndims = -1;
    H5Aread(attr, atype, &ndims);
    H5Aclose(attr);
    H5Tclose(atype);

    if (H5Aexists(file_id, "revision") <= 0) {
        std::cerr << "Error: missing attribute revision in HDF5 file\n";
        std::exit(1);
    }

    attr = H5Aopen(file_id, "revision", H5P_DEFAULT);
    atype = H5Aget_type(attr);

    int revision = -1;
    H5Aread(attr, atype, &revision);
    H5Aclose(attr);
    H5Tclose(atype);
}

HDF5Input::~HDF5Input()
{
    if (file_id >= 0) {
        H5Fclose(file_id);
        file_id = -1;
    }
}

template <typename T>
void HDF5Input::read_scaler(T& A, const std::string& name)
{
    hid_t dset_id = H5Dopen2(file_id, name.c_str(), H5P_DEFAULT);
    if (dset_id < 0) {
        std::cerr << "Error: cannot open dataset: " << name << "\n";
        std::exit(1);
    }
    hid_t space_id = H5Dget_space(dset_id);
    if (space_id < 0) {
        H5Dclose(dset_id);
        std::cerr << "Error: cannot get dataspace for " << name << "\n";
        std::exit(1);
    }
    int rank = H5Sget_simple_extent_ndims(space_id);
    if (rank < 0) {
        H5Sclose(space_id); H5Dclose(dset_id);
        std::cerr << "Error: cannot get rank for " << name << "\n";
        std::exit(1);
    }
    if (rank == 0 || rank > 1) {
        H5Sclose(space_id); H5Dclose(dset_id);
        std::cerr << "Error: dataset rank mismatch for " << name
                  << ", expected rank 1, got " << rank << "\n";
        std::exit(1);
    }

    hid_t mspace_id =  H5Screate(H5S_SCALAR);

    if (mspace_id < 0) {
        H5Sclose(space_id); H5Dclose(dset_id);
        std::cerr << "Error: cannot create memspace for " << name << "\n";
        std::exit(1);
    }

    hid_t dtype_id = H5Native<T>::id();

    if (H5Dread(dset_id, dtype_id, mspace_id, space_id, H5P_DEFAULT, &A) < 0) {
        H5Sclose(mspace_id); H5Sclose(space_id); H5Dclose(dset_id);
        std::cerr << "Error: failed to read dataset: " << name << "\n";
        std::exit(1);
    }

    H5Sclose(mspace_id);
    H5Sclose(space_id);
    H5Dclose(dset_id);
}

template
void HDF5Input::read_scaler<int>(int& A, const std::string& name);
template
void HDF5Input::read_scaler<double>(double& A, const std::string& name);

template <typename T>
void HDF5Input::read_array(std::vector<T>& A, const char *name, std::size_t size)
{
    size = size > 0 ? size : A.size();
    if (size == 0) {
        std::cerr << "Error: array size is 0: " << name << '\n';
        std::exit(1);
    }

    hid_t dset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
    if (dset_id < 0) {
        std::cerr << "Error: cannot open dataset: " << name << "\n";
        std::exit(1);
    }
    hid_t space_id = H5Dget_space(dset_id);
    if (space_id < 0) {
        H5Dclose(dset_id);
        std::cerr << "Error: cannot get dataspace for " << name << "\n";
        std::exit(1);
    }
    int rank = H5Sget_simple_extent_ndims(space_id);
    if (rank != 1) {
        H5Sclose(space_id); H5Dclose(dset_id);
        std::cerr << "Error: dataset rank mismatch for " << name
                  << ", expected rank 0 or 1, got " << rank << "\n";
        std::exit(1);
    }
    hsize_t dims[1];
    H5Sget_simple_extent_dims(space_id, dims, nullptr);
    if (dims[0] != size) {
        std::cerr << "Error: array size is not matched: " << name
                  << " (file dim = " << dims[0] << ", expected = " << size << ")\n";
        std::exit(1);
    }

    hid_t mspace_id = H5Screate_simple(1, dims, nullptr);
    hid_t dtype_id = H5Native<T>::id();

    if (H5Dread(dset_id, dtype_id, mspace_id, space_id, H5P_DEFAULT, A.data()) < 0) {
        H5Sclose(mspace_id); H5Sclose(space_id); H5Dclose(dset_id);
        std::cerr << "Error: failed to read dataset: " << name << "\n";
        std::exit(1);
    }

    H5Sclose(mspace_id);
    H5Sclose(space_id);
    H5Dclose(dset_id);
}


template <typename T, int N>
void HDF5Input::read_array(Array2D<T,N>& A, const char *name, std::size_t size)
{
    /* The caller must ensure A is of right size to hold the array */

    size = size > 0 ? size : A.size();
    if (A.size() == 0) {
        std::cerr << "Error: array size is 0: " << name << '\n';
        std::exit(1);
    }

    hid_t dset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
    if (dset_id < 0) {
        std::cerr << "Error: cannot open dataset: " << name << "\n";
        std::exit(1);
    }
    hid_t space_id = H5Dget_space(dset_id);
    if (space_id < 0) {
        H5Dclose(dset_id);
        std::cerr << "Error: cannot get dataspace for " << name << "\n";
        std::exit(1);
    }
    int rank = H5Sget_simple_extent_ndims(space_id);
    if (rank != 2) {
        std::cerr << "Error: dataset rank mismatch for " << name 
                  << ", expected 2 dims, got " << rank << '\n';
        std::exit(1);
    }

    hsize_t dims[2];
    H5Sget_simple_extent_dims(space_id, dims, nullptr);
    if (dims[0] != size || dims[1] != static_cast<hsize_t>(N)) {
        std::cerr << "Error: dataset dimensions mismatch for " << name
                  << ": file dims = (" << dims[0] << ", " << dims[1]
                  << "), expected (" << size << ", " << N << ")\n";
        std::exit(1);
    }
    hid_t mspace_id = H5Screate_simple(2, dims, nullptr);
    hid_t dtype_id = H5Native<T>::id();

    if (H5Dread(dset_id, dtype_id, mspace_id, space_id, H5P_DEFAULT, A.data()) < 0) {
        H5Sclose(mspace_id); H5Sclose(space_id); H5Dclose(dset_id);
        std::cerr << "Error: failed to read dataset: " << name << "\n";
        std::exit(1);
    }

    H5Sclose(mspace_id);
    H5Sclose(space_id);
    H5Dclose(dset_id);
}

// explicit instantiation
template
void HDF5Input::read_array<double>(double_vec& A, const char *name, std::size_t size);
template
void HDF5Input::read_array<int>(int_vec& A, const char *name, std::size_t size);
template
void HDF5Input::read_array<double,NDIMS>(Array2D<double,NDIMS>& A, const char *name, std::size_t size);
template
void HDF5Input::read_array<double,NSTR>(Array2D<double,NSTR>& A, const char *name, std::size_t size);
#ifdef THREED // when 2d, NSTR == NODES_PER_ELEM == 3
template
void HDF5Input::read_array<double,NODES_PER_ELEM>(Array2D<double,NODES_PER_ELEM>& A, const char *name, std::size_t size);
#endif
template
void HDF5Input::read_array<double,1>(Array2D<double,1>& A, const char *name, std::size_t size);
template
void HDF5Input::read_array<int,NDIMS>(Array2D<int,NDIMS>& A, const char *name, std::size_t size);
template
void HDF5Input::read_array<int,NODES_PER_ELEM>(Array2D<int,NODES_PER_ELEM>& A, const char *name, std::size_t size);
template
void HDF5Input::read_array<int,1>(Array2D<int,1>& A, const char *name, std::size_t size);

#endif