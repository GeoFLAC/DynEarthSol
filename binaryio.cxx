#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "binaryio.hpp"
#include "utils.hpp"
#include "markerset.hpp"
#include "runtime_info.hpp"

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
    const std::string revision_str = "# DynEarthSol ndims=" + std::to_string(NDIMS)
                                   + " revision=" + std::to_string(BINARY_FILE_REVISION) + "\n";
}


/* Not using C++ stream IO for bulk file io since it can be much slower than C stdio. */

void rename_to_old_backup(const char *filename) {
    // Find the highest-numbered existing backup (.old, .old2, .old3, ...) and
    // rename filename to max+1, so new backups always append after the largest
    // rather than filling gaps left by manual deletions.
    // fopen on a nonexistent file returns NULL immediately (ENOENT, no disk IO),
    // so probing 200 past the last found is effectively free.
    std::string fpath(filename);
    int max_n = 0;
    for (int n = 1; n <= max_n + 200; ++n) {
        std::string candidate = fpath + ".old" + (n == 1 ? "" : std::to_string(n));
        std::FILE *f = std::fopen(candidate.c_str(), "r");
        if (f) { std::fclose(f); max_n = n; }
    }
    int next_n = max_n + 1;
    std::string backup = fpath + ".old" + (next_n == 1 ? "" : std::to_string(next_n));
    if (std::rename(filename, backup.c_str()) == 0)
        std::cerr << "[Runtime][IO] Renamed '" << filename << "' -> '" << backup
                  << "' (preserving previous output)\n";
}

#ifndef HDF5

BinaryOutput::BinaryOutput(const char *filename, const bool rename_if_exists)
{
    if (rename_if_exists) rename_to_old_backup(filename);
    f = std::fopen(filename, "wb");
    if (f == NULL) {
        std::cerr << "Error: cannot open file: " << filename << '\n';
        die(EXIT_IO_OPEN);
    }

    header = new char[headerlen]();
    hd_pos = std::strcat(header, revision_str.c_str());
    eof_pos = headerlen;

    std::fseek(f, eof_pos, SEEK_SET);
}


BinaryOutput::~BinaryOutput()
{
#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
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
#ifdef NPROF
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
        die(EXIT_INTERNAL_ASSERT);
    }
    if (len >= headerlen - (hd_pos - header)*sizeof(char)) {
        std::cerr << "Error: exceeding header length at Output::write_array, name=" << name
                  << " eof_position=" << eof_pos << '\n';
        die(EXIT_INTERNAL_ASSERT);
    }
    hd_pos = std::strncat(hd_pos, buffer, len);
}

template <typename T>
void BinaryOutput::write_scalar(const T& A, const std::string& name)
{
    write_header(name.c_str());
    std::size_t n = std::fwrite(&A, sizeof(T), 1, f);
    eof_pos += n * sizeof(T);
}


// explicit instantiation
template
void BinaryOutput::write_scalar<int>(const int& A, const std::string& name);
template
void BinaryOutput::write_scalar<double>(const double& A, const std::string& name);

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

    std::size_t written_elements = 0;
    static std::vector<T> buffer; 
    A.pack_to(buffer, size);

    written_elements = std::fwrite(buffer.data(), sizeof(T), size * N, f);

    if (written_elements != size * N) {
        std::cerr << "Error: cannot write array: " << name << '\n';
    }

    eof_pos += written_elements * sizeof(T);
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

namespace {

void kv(std::string& p, const char* name, const std::string& value)
{
    p += name; p += '='; p += value; p += '\n';
}

void kv(std::string& p, const char* name, int value)
{
    char b[32];
    std::snprintf(b, sizeof(b), "%d", value);
    kv(p, name, std::string(b));
}

// %.9g: 6 digits quantize mem_rss_gib to ~107 KiB at 10 GiB, coarser than the pages a
// footprint is counted in. %g strips trailing zeros at any precision, so "64" for a whole
// number is not precision loss.
void kv(std::string& p, const char* name, double value)
{
    char b[32];
    std::snprintf(b, sizeof(b), "%.9g", value);
    kv(p, name, std::string(b));
}

} // anonymous namespace

// The identity an HDF5 frame carries in its /provenance group, as one text record so a
// des-binary frame explains itself too -- readable with strings(1). Keep the field set
// in step with HDF5Output::write_run_provenance.
// No revision bump: records are looked up by name, so one more is ignored by a
// reader that does not know it, and bumping would refuse every existing checkpoint.
void BinaryOutput::write_run_provenance(const BuildInfo& build, const CpuInfo& cpu,
                                        const DeviceInfo& dev, const std::string& restart_from,
                                        double rss_gib, double peak_rss_gib)
{
    const std::string unknown("unknown");
    std::string p;

    kv(p, "code_rev", build.rev);
    kv(p, "code_branch", build.branch);
    kv(p, "code_dirty", build.dirty);
    kv(p, "code_origin", build.origin);
    kv(p, "code_state_utc", build.state_utc);
    kv(p, "build_os", build.build_os);
    kv(p, "builder", build.builder);
    kv(p, "exe_mtime_utc", build.exe_mtime_utc);

    kv(p, "os", cpu.os);
    kv(p, "runner", cpu.runner);
    kv(p, "cpu_model", cpu.model);
    // Machine context, not a denominator: the cpu_time ratio below divides by
    // omp_threads, and logical_cores is blind to affinity and cgroup limits.
    kv(p, "logical_cores", cpu.logical_cores);
    kv(p, "mem_total_gib", cpu.mem_total_gib);
    // The parallel width, not the machine's core counts: the PT loop exits on an
    // OpenMP reduction whose summation order depends on the team size. Sampled here
    // rather than at startup, since a dynamic team drifts between frames.
    kv(p, "omp_threads", omp_team_size_now());

    kv(p, "kernel", dev.kernel);
    // The filter reads "0/1" on a filtered 8-GPU node, so gpu_device alone cannot name the
    // physical card. Guarded on acc_build rather than using_gpu, which an ACC build that
    // reaches a frame now implies -- a build that could fall back to the host would not.
    if (dev.acc_build)
        kv(p, "gpu_visible_devices", dev.visible_devices);
    if (dev.using_gpu) {
        kv(p, "gpu_model", dev.name.empty() ? unknown : dev.name);
        kv(p, "gpu_device", std::to_string(dev.active_dev) + "/"
                          + std::to_string(dev.num_devices));
        kv(p, "gpu_cuda_driver", dev.cuda_driver.empty() ? unknown : dev.cuda_driver);
        kv(p, "gpu_mem_total_gib", dev.mem_total_gib);
        kv(p, "gpu_mem_free_at_start_gib", dev.mem_free_gib);
        kv(p, "gpu_mem_used_dev_gib", device_mem_used_dev_gib(dev));
    }

    // Without this a lone frame from a restarted run credits the whole state to the
    // last leg's binary and host -- a wrong claim, not a missing one.
    kv(p, "restart_from", restart_from);

    // Machine health and resource use AT THIS FRAME. Sampled every frame, so a run's
    // frames form a time series: rss against frame number shows a leak, the
    // cpu_time/(walltime x omp_threads) ratio a stalling region, and free/load
    // whether the machine itself went bad.
    kv(p, "mem_rss_gib", rss_gib);
    kv(p, "mem_peak_rss_gib", peak_rss_gib);
    kv(p, "mem_avail_gib", host_mem_avail_gib_now());
    kv(p, "cpu_time_sec", process_cpu_time_sec());
    kv(p, "load_avg_1m", host_load_avg_1m());
    kv(p, "write_utc", utc_now());

    write_header("provenance");
    eof_pos += std::fwrite(p.data(), sizeof(char), p.size(), f) * sizeof(char);
}


void BinaryOutput::write_nodal_vec_array(const Array2D<double,NDIMS>& A, const char *name, std::size_t len)
{
    write_array(A, name, len);
}

//////////////////////////////////////////////////////////////////////////////

BinaryInput::BinaryInput(const char *filename)
{
    f = std::fopen(filename, "r");
    if (f == NULL) {
        std::cerr << "Error: cannot open file: " << filename << '\n';
        die(EXIT_IO_OPEN);
    }
    read_header();
}


BinaryInput::~BinaryInput()
{
    std::fclose(f);
}


bool BinaryInput::has_array(const char *name) const
{
    return offset.find(name) != offset.end();
}


void BinaryInput::read_header()
{
    /* Read into header buffer */
    std::fseek(f, 0, SEEK_SET);
    char *header = new char[headerlen]();
    std::size_t n = std::fread(header, sizeof(char), headerlen, f);
    if (n != headerlen) {
        die(EXIT_IO_RW, "error reading file header");
    }

    /* Parse the content of header buffer */
    char *line = header;

    // Compare revision string (excluding the trailing new line)
    line = std::strtok(header, "\n");
    if (strncmp(line, revision_str.c_str(), revision_str.size()-1) != 0) {
        std::cerr << "Error: mismatching revision string in header\n"
                  << "  Expect: " << revision_str
                  << "  Got: "<< line << '\n';
        die(EXIT_IO_RESTART);
    }

    line = std::strtok(NULL, "\n");
    while (line != NULL) {
        /* Each line is a string (might contain space), a tab, and an integer */
        char *tab = std::strchr(line, '\t');
        if (tab == NULL) {
            std::cerr << "Error: error parsing file header\n"
                      << " Line is:" << line << '\n';
            die(EXIT_IO_RW);
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
        die(EXIT_IO_RESTART);
    }
    std::size_t loc = it->second;
    //std::cout << name << ' ' << loc << '\n';
    std::fseek(f, loc, SEEK_SET);
}


template <typename T>
void BinaryInput::read_scalar(T& A, const std::string& name)
{
    seek_to_array(name.c_str());
    std::size_t n = std::fread(&A, sizeof(T), 1, f);
    if (n != 1) {
        std::cerr << "Error: cannot read scalar: " << name << '\n';
        die(EXIT_IO_RW);
    }
}


// explicit instantiation
template
void BinaryInput::read_scalar<int>(int& A, const std::string& name);
template
void BinaryInput::read_scalar<double>(double& A, const std::string& name);


template <typename T>
void BinaryInput::read_array(std::vector<T>& A, const char *name, std::size_t size)
{
    /* The caller must ensure A is of right size to hold the array */

    size = size > 0 ? size : A.size();
    if (A.size() == 0) {
        std::cerr << "Error: array size is 0: " << name << '\n';
        die(EXIT_IO_RW);
    }

    seek_to_array(name);
    const std::size_t n = std::fread(A.data(), sizeof(T), size, f);

    if (n != size) {
        std::cerr << "Error: cannot read array: " << name << '\n';
        die(EXIT_IO_RW);
    }
}


template <typename T, int N>
void BinaryInput::read_array(Array2D<T,N>& A, const char *name, std::size_t size)
{
    /* The caller must ensure A is of right size to hold the array */

    size = size > 0 ? size : A.size();
    if (A.size() == 0) {
        std::cerr << "Error: array size is 0: " << name << '\n';
        die(EXIT_IO_RW);
    }

    seek_to_array(name);

    static std::vector<T> buffer;
    std::size_t total_elements = size * N;
    if (buffer.size() < total_elements)
        buffer.resize(total_elements);

    const std::size_t n = std::fread(buffer.data(), sizeof(T), size * N, f);
    if (n != N * size) {
        std::cerr << "Error: cannot read array (buffered path): " << name << '\n';
        die(EXIT_IO_RW);
    }

    A.load_from_buffer(buffer.data(), size, false);
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

HDF5Output::HDF5Output(const char *filename, const int hdf5_compression_level,
                       const bool is_chkpt, const bool rename_if_exists)
    : compression_level(hdf5_compression_level), is_checkpoint(is_chkpt)
{
    if (rename_if_exists) rename_to_old_backup(filename);

    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    // Locking off (use_file_locking = false) so ParaView can read a frame while the run
    // still holds it; ignore_when_disabled = true tolerates builds that already disabled
    // it. HDF5 then no longer refuses a second writer, so two runs sharing a modelname
    // corrupt the file silently. H5Pset_file_locking needs HDF5 >= 1.10.7.
#if H5_VERSION_GE(1, 10, 7)
    H5Pset_file_locking(fapl_id, false, true);
#endif
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    H5Pclose(fapl_id);

    if (file_id < 0) {
        throw std::runtime_error(std::string("H5Fcreate failed: ") + filename);
    }

    write_header();
}

HDF5Output::~HDF5Output()
{
#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    if (file_id >= 0) {
        H5Fflush(file_id, H5F_SCOPE_GLOBAL);
        H5Fclose(file_id);
        file_id = -1;
    }
#ifdef NPROF
    nvtxRangePop();
#endif
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
    write_attribute(BINARY_FILE_REVISION, "revision", file_id);

    hid_t gid = create_group_with_order("/VTKHDF");

    std::string vtkhdf_type = "PartitionedDataSetCollection";
    write_attribute(vtkhdf_type, "Type", gid);

    int_vec version = {2, 1};
    write_attribute(version, "Version", 2, gid);
    H5Gclose(gid);

    gid = create_group_with_order("/VTKHDF/Assembly");
    H5Gclose(gid);
}

// Host, build and device provenance in a /provenance group, outside /VTKHDF where it
// cannot collide with the schema. Per frame on purpose, so a lone frame explains itself
// without the .manifest. Integer counts use -1 = unknown and the _gib doubles use 0, so
// consumers must reject non-positives rather than computing with them.
void HDF5Output::write_run_provenance(const BuildInfo& build, const CpuInfo& cpu,
                                      const DeviceInfo& dev, const std::string& restart_from,
                                      double rss_gib, double peak_rss_gib)
{
    // Its own group, not the file root: the root is the container's namespace, and a
    // frame's identity is not part of it. Same name as the des-binary record.
    hid_t g = H5Gcreate2(file_id, "/provenance", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (g < 0) return;   // metadata must never abort a run

    // The build's code and host groups: which source state, built where and when.
    write_attribute(build.rev, "code_rev", g);
    write_attribute(build.branch, "code_branch", g);
    write_attribute(build.dirty, "code_dirty", g);
    write_attribute(build.origin, "code_origin", g);
    write_attribute(build.state_utc, "code_state_utc", g);
    write_attribute(build.build_os, "build_os", g);
    write_attribute(build.builder, "builder", g);
    write_attribute(build.exe_mtime_utc, "exe_mtime_utc", g);

    // The run's host: a g++ and an nvc++ build, and a GPU run and its host
    // fallback, do not agree bit for bit -- and neither do two hosts.
    write_attribute(cpu.os, "os", g);
    write_attribute(cpu.runner, "runner", g);
    write_attribute(cpu.model, "cpu_model", g);
    // Machine context, not a denominator: the cpu_time ratio below divides by
    // omp_threads, and logical_cores is blind to affinity and cgroup limits.
    write_attribute(cpu.logical_cores, "logical_cores", g);
    write_attribute(cpu.mem_total_gib, "mem_total_gib", g);
    // The parallel width, not the machine's core counts: the PT loop exits on an
    // OpenMP reduction whose summation order depends on the team size. Sampled at this
    // write, since a dynamic team drifts.
    write_attribute(omp_team_size_now(), "omp_threads", g);

    write_attribute(dev.kernel, "kernel", g);
    // Guarded on acc_build rather than using_gpu, which an ACC build that reaches a frame
    // now implies -- a build that could fall back to the host would not.
    if (dev.acc_build)
        write_attribute(dev.visible_devices, "gpu_visible_devices", g);
    if (dev.using_gpu) {
        // The contract for a frame is that an unknown string reads "unknown"; the
        // probes leave a field they could not read empty.
        const std::string unknown("unknown");
        write_attribute(dev.name.empty() ? unknown : dev.name, "gpu_model", g);
        // Same N/M form as the manifest's selected=, so a grep for either finds both.
        const std::string selected(std::to_string(dev.active_dev) + "/"
                                   + std::to_string(dev.num_devices));
        write_attribute(selected, "gpu_device", g);
        write_attribute(dev.cuda_driver.empty() ? unknown : dev.cuda_driver,
                        "gpu_cuda_driver", g);
        write_attribute(dev.mem_total_gib, "gpu_mem_total_gib", g);
        write_attribute(dev.mem_free_gib, "gpu_mem_free_at_start_gib", g);
    }

    // Status at write time: how the run was doing when this frame was written.
    // Without this a lone frame from a restarted run credits the whole state to the
    // last leg's binary and host -- a wrong claim, not a missing one.
    write_attribute(restart_from, "restart_from", g);

    // Machine health and resource use AT THIS FRAME. Sampled every frame, so a run's
    // frames form a time series: rss against frame number shows a leak, the
    // cpu_time/(walltime x omp_threads) ratio a stalling region, and free/load
    // whether the machine itself went bad.
    write_attribute(rss_gib, "mem_rss_gib", g);
    write_attribute(peak_rss_gib, "mem_peak_rss_gib", g);
    write_attribute(host_mem_avail_gib_now(), "mem_avail_gib", g);
    write_attribute(process_cpu_time_sec(), "cpu_time_sec", g);
    write_attribute(host_load_avg_1m(), "load_avg_1m", g);
    if (dev.using_gpu)
        write_attribute(device_mem_used_dev_gib(dev), "gpu_mem_used_dev_gib", g);
    write_attribute(utc_now(), "write_utc", g);

    H5Gclose(g);
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
            write_nodal_vec_array(*var.coord, "Points", nnode);

            int_vec buffer;
            var.connectivity->pack_to(buffer);
            int* conn_ptr = buffer.data();
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
            array_t mcoord(nnode);
            ms->calculate_marker_coord(var, mcoord); // coordinate of markers
            write_nodal_vec_array(mcoord, "Points", nnode);

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

    if (name != "NumberOfConnectivityIds") {
        std::string vis_name = name;
        bool is_field = false;
        if (kind == "marker") {
            if (name == "NumberOfPoints") {
                vis_name = block_base + ".nmarkers";
                is_field = true;
            } else if (name == "NumberOfCells") {
                // do nothing
            } else {
                create_virtual_dataset(full_name, vis_name, space_id, dtype_id);
            }
        } else if (kind == "grid") {
            if (name == "NumberOfPoints") {
                vis_name = "nnode";
                is_field = true;
            } else if (name == "NumberOfCells") {
                vis_name = "nelem";
                is_field = true;
            } else {
                create_virtual_dataset(full_name, vis_name, space_id, dtype_id);
            }
        }
        
        if (is_field) {
            create_virtual_dataset(full_name, vis_name, space_id, dtype_id); // Create at root for restart/legacy
            std::string field_name = "/VTKHDF/grid/FieldData/" + vis_name;
            create_virtual_dataset(full_name, field_name, space_id, dtype_id); // Create in FieldData for ParaView
        }
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
            die(EXIT_INTERNAL_ASSERT);
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

    bool skip_virtual = (std::string(name) == "Offsets" || std::string(name) == "Types");
    if (kind == "marker" && std::string(name) == "Connectivity") skip_virtual = true;

    if (!skip_virtual) {
        if (std::string(name) == "Connectivity") {
            int len2D = len / nnode_cell;
            create_virtual_dataset(full_name, "connectivity", space_id, dtype_id, len2D, nnode_cell);
        } else {
            create_virtual_dataset(full_name, name, space_id, dtype_id, len);
        }
    }
    H5Dclose(dset_id);
    H5Pclose(dcpl_id);
    H5Sclose(space_id);
}

// 2D array
template<typename T, int N>
void HDF5Output::write_array(const Array2D<T, N>& A, const char *name, hsize_t len, int dest_N)
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
            die(EXIT_INTERNAL_ASSERT);
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

    static std::vector<T> buffer;
    A.pack_to(buffer, len);

    H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());

    std::string vis_name = name;
    if (std::string(name) == "Points") {
        if (kind == "grid") {
            vis_name = "coordinate";
        } else {
            vis_name = block_base + ".coord";
        }
    }

    create_virtual_dataset(full_name, vis_name, space_id, dtype_id, len, N, dest_N);

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
void HDF5Output::write_array<double,NDIMS>(const Array2D<double,NDIMS>& A, const char *name, hsize_t, int);
template
void HDF5Output::write_array<double,NSTR>(const Array2D<double,NSTR>& A, const char *name, hsize_t, int);
#ifdef THREED // when 2d, NSTR == NODES_PER_ELEM == 3
template
void HDF5Output::write_array<double,NODES_PER_ELEM>(const Array2D<double,NODES_PER_ELEM>& A, const char *name, hsize_t, int);
#endif
template
void HDF5Output::write_array<double,1>(const Array2D<double,1>& A, const char *name, hsize_t, int);
template
void HDF5Output::write_array<int,NODES_PER_ELEM>(const Array2D<int,NODES_PER_ELEM>& A, const char *name, hsize_t, int);
template
void HDF5Output::write_array<int,NDIMS>(const Array2D<int,NDIMS>& A, const char *name, hsize_t, int);
template
void HDF5Output::write_array<int,1>(const Array2D<int,1>& A, const char *name, hsize_t, int);

void HDF5Output::write_nodal_vec_array(const Array2D<double,NDIMS>& A, const char *name, hsize_t len)
{
#ifdef THREED
    write_array(A, name, len);
#else
    // Store as 3-component in PointData for ParaView glyph arrows,
    // but keep the root virtual dataset at NDIMS components so legacy
    // readers (Dynearthsol.py / 2vtk.py) continue to read the correct shape.
    Array2D<double, 3> A3d(len);
    #pragma omp parallel for default(none) shared(len, A3d, A)
    for (hsize_t i=0; i<len; ++i) {
        A3d[i][0] = A[i][0];
        A3d[i][1] = A[i][1];
        A3d[i][2] = 0.0;
    }
    write_array(A3d, name, len, NDIMS);
#endif
}

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
void HDF5Output::create_virtual_dataset(const std::string& src_name, const std::string& dest_name, hid_t& space_id, hid_t& dtype_id, hsize_t len, int N, int dest_N)
{
    if (dest_N == -1) dest_N = N;

    hsize_t vds_dims[2] = { len, static_cast<hsize_t>(dest_N) };
    hid_t vds_space_id = H5Screate_simple(2, vds_dims, nullptr);

    hid_t vds_dcpl = H5Pcreate(H5P_DATASET_CREATE);
    hid_t src_space_id = H5Scopy(space_id);

    // If source and dest dimensions differ (e.g. 3D source -> 2D virtual), select hyperslab
    if (N != dest_N) {
        hsize_t start[2] = {0, 0};
        hsize_t count[2] = {len, static_cast<hsize_t>(dest_N)};

        // Select first dest_N columns from source
        herr_t status = H5Sselect_hyperslab(src_space_id, H5S_SELECT_SET, start, nullptr, count, nullptr);
        if (status < 0) {
           std::cerr << "Error selecting source hyperslab for VDS " << dest_name << "\n";
        }

        // Select all of virtual dataset (which is dest_N wide)
        status = H5Sselect_hyperslab(vds_space_id, H5S_SELECT_SET, start, nullptr, count, nullptr);
        if (status < 0) {
           std::cerr << "Error selecting virtual hyperslab for VDS " << dest_name << "\n";
        }
    }

    herr_t status = H5Pset_virtual(vds_dcpl, vds_space_id, ".", src_name.c_str(), src_space_id);
    if (status < 0) {
        std::cerr << "Error setting virtual dataset mapping for " << dest_name << "\n";
    }

    hid_t vds_dset_id = H5Dcreate2(file_id, dest_name.c_str(), dtype_id, vds_space_id, H5P_DEFAULT, vds_dcpl, H5P_DEFAULT);
    if (vds_dset_id < 0) {
        std::cerr << "Error creating virtual dataset " << dest_name << "\n";
    }

    H5Dclose(vds_dset_id);
    H5Sclose(src_space_id);
    H5Pclose(vds_dcpl);
    H5Sclose(vds_space_id);
}

HDF5Input::HDF5Input(const char *filename)
{
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    // Locking off so a restart can read a file ParaView holds open; the trade-off is
    // spelled out in the HDF5Output constructor above.
#if H5_VERSION_GE(1, 10, 7)
    H5Pset_file_locking(fapl_id, false, true);
#endif
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, fapl_id);
    H5Pclose(fapl_id);
    if (file_id < 0) {
        std::cerr << "Error: cannot open HDF5 file for reading: " << filename << "\n";
        die(EXIT_IO_OPEN);
    }

    read_header();
}

void HDF5Input::read_header()
{
    if (H5Aexists(file_id, "ndims") <= 0) {
        die(EXIT_IO_RESTART, "missing attribute ndims in HDF5 file");
    }

    hid_t attr = H5Aopen(file_id, "ndims", H5P_DEFAULT);
    hid_t atype = H5Aget_type(attr);

    int ndims = -1;
    H5Aread(attr, atype, &ndims);
    H5Tclose(atype);
    H5Aclose(attr);

    if (ndims != NDIMS) {
        std::cerr << "Error: mismatching ndims in HDF5 file\n"
                  << "  Expect: " << NDIMS << "  Got: " << ndims << '\n';
        die(EXIT_IO_RESTART);
    }

    if (H5Aexists(file_id, "revision") <= 0) {
        die(EXIT_IO_RESTART, "missing attribute revision in HDF5 file");
    }

    attr = H5Aopen(file_id, "revision", H5P_DEFAULT);
    atype = H5Aget_type(attr);

    int revision = -1;
    H5Aread(attr, atype, &revision);
    H5Tclose(atype);
    H5Aclose(attr);

    if (revision != BINARY_FILE_REVISION) {
        std::cerr << "Error: mismatching revision in HDF5 file\n"
                  << "  Expect: " << BINARY_FILE_REVISION << "  Got: " << revision << '\n';
        die(EXIT_IO_RESTART);
    }
}

HDF5Input::~HDF5Input()
{
    if (file_id >= 0) {
        H5Fclose(file_id);
        file_id = -1;
    }
}


bool HDF5Input::has_array(const char *name) const
{
    return H5Lexists(file_id, name, H5P_DEFAULT) > 0;
}

template <typename T>
void HDF5Input::read_scalar(T& A, const std::string& name)
{
    hid_t dset_id = H5Dopen2(file_id, name.c_str(), H5P_DEFAULT);
    if (dset_id < 0) {
        std::cerr << "Error: cannot open dataset: " << name << "\n";
        die(EXIT_IO_RW);
    }
    hid_t space_id = H5Dget_space(dset_id);
    if (space_id < 0) {
        H5Dclose(dset_id);
        std::cerr << "Error: cannot get dataspace for " << name << "\n";
        die(EXIT_IO_RW);
    }
    int rank = H5Sget_simple_extent_ndims(space_id);
    if (rank < 0) {
        H5Sclose(space_id); H5Dclose(dset_id);
        std::cerr << "Error: cannot get rank for " << name << "\n";
        die(EXIT_IO_RESTART);
    }
    if (rank == 0 || rank > 1) {
        H5Sclose(space_id); H5Dclose(dset_id);
        std::cerr << "Error: dataset rank mismatch for " << name
                  << ", expected rank 1, got " << rank << "\n";
        die(EXIT_IO_RW);
    }

    hid_t mspace_id =  H5Screate(H5S_SCALAR);

    if (mspace_id < 0) {
        H5Sclose(space_id); H5Dclose(dset_id);
        std::cerr << "Error: cannot create memspace for " << name << "\n";
        die(EXIT_IO_RW);
    }

    hid_t dtype_id = H5Native<T>::id();

    if (H5Dread(dset_id, dtype_id, mspace_id, space_id, H5P_DEFAULT, &A) < 0) {
        H5Sclose(mspace_id); H5Sclose(space_id); H5Dclose(dset_id);
        std::cerr << "Error: failed to read dataset: " << name << "\n";
        die(EXIT_IO_RW);
    }

    H5Sclose(mspace_id);
    H5Sclose(space_id);
    H5Dclose(dset_id);
}

template
void HDF5Input::read_scalar<int>(int& A, const std::string& name);
template
void HDF5Input::read_scalar<double>(double& A, const std::string& name);

template <typename T>
void HDF5Input::read_array(std::vector<T>& A, const char *name, std::size_t size)
{
    size = size > 0 ? size : A.size();
    if (size == 0) {
        std::cerr << "Error: array size is 0: " << name << '\n';
        die(EXIT_IO_RW);
    }

    hid_t dset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
    if (dset_id < 0) {
        std::cerr << "Error: cannot open dataset: " << name << "\n";
        die(EXIT_IO_RW);
    }
    hid_t space_id = H5Dget_space(dset_id);
    if (space_id < 0) {
        H5Dclose(dset_id);
        std::cerr << "Error: cannot get dataspace for " << name << "\n";
        die(EXIT_IO_RW);
    }
    int rank = H5Sget_simple_extent_ndims(space_id);
    if (rank != 1) {
        H5Sclose(space_id); H5Dclose(dset_id);
        std::cerr << "Error: dataset rank mismatch for " << name
                  << ", expected rank 0 or 1, got " << rank << "\n";
        die(EXIT_IO_RESTART);
    }
    hsize_t dims[1];
    H5Sget_simple_extent_dims(space_id, dims, nullptr);
    if (dims[0] != size) {
        std::cerr << "Error: array size is not matched: " << name
                  << " (file dim = " << dims[0] << ", expected = " << size << ")\n";
        die(EXIT_IO_RW);
    }

    hid_t mspace_id = H5Screate_simple(1, dims, nullptr);
    hid_t dtype_id = H5Native<T>::id();

    if (H5Dread(dset_id, dtype_id, mspace_id, space_id, H5P_DEFAULT, A.data()) < 0) {
        H5Sclose(mspace_id); H5Sclose(space_id); H5Dclose(dset_id);
        std::cerr << "Error: failed to read dataset: " << name << "\n";
        die(EXIT_IO_RW);
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
        die(EXIT_IO_RW);
    }

    hid_t dset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
    if (dset_id < 0) {
        std::cerr << "Error: cannot open dataset: " << name << "\n";
        die(EXIT_IO_RW);
    }
    hid_t space_id = H5Dget_space(dset_id);
    if (space_id < 0) {
        H5Dclose(dset_id);
        std::cerr << "Error: cannot get dataspace for " << name << "\n";
        die(EXIT_IO_RW);
    }
    int rank = H5Sget_simple_extent_ndims(space_id);
    if (rank != 2) {
        std::cerr << "Error: dataset rank mismatch for " << name 
                  << ", expected 2 dims, got " << rank << '\n';
        die(EXIT_IO_RESTART);
    }

    hsize_t dims[2];
    H5Sget_simple_extent_dims(space_id, dims, nullptr);
    if (dims[0] != size || dims[1] != static_cast<hsize_t>(N)) {
        std::cerr << "Error: dataset dimensions mismatch for " << name
                  << ": file dims = (" << dims[0] << ", " << dims[1]
                  << "), expected (" << size << ", " << N << ")\n";
        die(EXIT_IO_RW);
    }
    hid_t mspace_id = H5Screate_simple(2, dims, nullptr);
    hid_t dtype_id = H5Native<T>::id();

    static std::vector<T> buffer;
    std::size_t total_elements = size * N;
    if (buffer.size() < total_elements)
        buffer.resize(total_elements);

    if (H5Dread(dset_id, dtype_id, mspace_id, space_id, H5P_DEFAULT, buffer.data()) < 0) {
        H5Sclose(mspace_id); H5Sclose(space_id); H5Dclose(dset_id);
        std::cerr << "Error: failed to read dataset: " << name << "\n";
        die(EXIT_IO_RW);
    }

    A.load_from_buffer(buffer.data(), size, false);

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
