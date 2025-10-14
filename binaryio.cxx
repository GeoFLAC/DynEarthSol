#include <cstdio>
#include <cstring>
#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "binaryio.hpp"

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

HDF5Output::HDF5Output(const char *filename, const int hdf5_compression_level)
    : h5_file(filename, H5F_ACC_TRUNC), compression_level(hdf5_compression_level)
{
    write_header();
}

HDF5Output::~HDF5Output()
{
    // nothing to do
}

void HDF5Output::write_header()
{
    hsize_t dims = 1;
    DataSpace attr_space = DataSpace(H5S_SCALAR);
    Attribute attr_nd = h5_file.createAttribute("ndims", PredType::NATIVE_INT, attr_space);
    int nd = NDIMS;
    attr_nd.write(PredType::NATIVE_INT, &nd);

    Attribute attr_rev = h5_file.createAttribute("revision", PredType::NATIVE_INT, attr_space);
    int rev = 3;
    attr_rev.write(PredType::NATIVE_INT, &rev);
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

// 1D array
template<typename T>
void HDF5Output::write_array(const std::vector<T> &A, const char *name, hsize_t len)
{
    hsize_t dims[1] = { len };
    DataSpace dataspace(1, dims);

    PredType dtype = H5TypeMap<T>::type();

    DSetCreatPropList prop;
    hsize_t chunk_dim = (len < 1024 ? len : 1024);
    prop.setChunk(1, &chunk_dim);
    prop.setDeflate(compression_level);
    prop.setShuffle();

    DataSet dataset = h5_file.createDataSet(name, dtype, dataspace, prop);
    dataset.write(A.data(), dtype);
}

// 2D array
template<typename T, int N>
void HDF5Output::write_array(const Array2D<T, N>& A, const char *name, hsize_t len)
{
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

    DataSet dataset = h5_file.createDataSet(name, dtype, dataspace);
    dataset.write(A.data(), dtype);
}

// explicit instantiation
template
void HDF5Output::write_array<int>(const int_vec& A, const char *name, hsize_t);
template
void HDF5Output::write_array<double>(const double_vec& A, const char *name, hsize_t);
template
void HDF5Output::write_array<uint>(const std::vector<uint>& A, const char *name, hsize_t);
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