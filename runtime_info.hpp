#ifndef DYNEARTHSOL_RUNTIME_INFO_HPP
#define DYNEARTHSOL_RUNTIME_INFO_HPP

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "parameters.hpp"

// Host CPU and OpenMP thread assignment of this process. Integer counts are -1
// when the platform does not expose them; the _gib fields are <= 0 when unknown.
struct CpuInfo {
    std::string model;     // e.g. "Apple M3 Max"
    std::string os;        // the run host's OS, e.g. "macos-26.6", "ubuntu-22.04"
    std::string runner;    // user@host -- the run-side counterpart of the build's builder=
    int logical_cores;     // hardware threads on the machine
    int physical_cores;
    int perf_cores;        // hybrid CPUs only (Apple performance cores)
    int eff_cores;         // hybrid CPUs only (Apple efficiency cores)
    double mem_total_gib;
    double mem_avail_gib;  // at probe time: what other processes left, not a host property
    bool openmp_enabled;   // false in an openmp=0 build; the omp_* fields are -1
    int omp_team_size;     // threads the runtime actually hands out (measured)
    int omp_max_threads;
    int omp_num_procs;     // processors available to THIS process: affinity- and cgroup-scoped
    bool omp_dynamic;
    std::string wait_policy_src;   // "env", "des-default" (main()'s macOS default) or "runtime"
};

// The embedded build.snapshot's code and host groups as fields, for the frame writers;
// the build configuration stays in the block alone.
struct BuildInfo {
    std::string rev;        // git describe against version tags
    std::string branch;
    std::string dirty;      // file/±line counts ("2f/+13/-12"), or "clean"
    std::string origin;     // credential-stripped remote URL
    std::string state_time; // when the source/configuration state last changed, not a build time
    std::string build_os;   // the build host's OS
    std::string builder;    // user@host that built the binary
    std::string exe_mtime;  // the executable file's mtime: a cp or a touch moves it
};

// The offload device this process got, not a build property: an openacc=1 binary falls
// back to the host when no NVIDIA device is visible. Empty or non-positive = unknown.
struct DeviceInfo {
    bool acc_build;        // compiled with -acc
    bool using_gpu;        // an NVIDIA device was selected AND initialized
    std::string kernel;    // "CPU" or "GPU" -- where this run's arithmetic happened
    std::string name;      // e.g. "NVIDIA A100-SXM4-40GB"
    std::string cuda_driver;   // the CUDA driver API version, not nvidia-smi's driver_version
    // CUDA_VISIBLE_DEVICES as it selected the device: "all" unset, "none" set but empty
    // (which hides every device), else its value.
    std::string visible_devices;
    int num_devices;       // NVIDIA devices visible to this process
    int active_dev;        // index among them, -1 when running on the host
    double mem_total_gib;
    double mem_free_gib;   // free at selection, after acc_init took this run's own context
};

BuildInfo probe_build_info();
// wait_policy_from_des: main() set OMP_WAIT_POLICY itself, which only main() can know.
CpuInfo probe_cpu_info(bool wait_policy_from_des);
// Selects and initializes the offload device, then describes it: call it once, before any
// compute. A non-ACC build only fills in the CPU answer.
DeviceInfo init_offload_device();
// A non-empty OMP_WAIT_POLICY or KMP_BLOCKTIME is set.
bool env_sets_wait_policy();

// One "[section]" of the .manifest, fields in file order; the screen prints the same ones.
struct ManifestSection {
    std::string name;
    std::vector<std::pair<std::string, std::string> > fields;
};
typedef std::vector<ManifestSection> Manifest;

// The [runtime.*] sections, then the embedded block as [build.*] sections: composed once,
// so file and screen agree.
Manifest compose_manifest(const Param& param, const BuildInfo& build,
                          const CpuInfo& cpu, const DeviceInfo& dev);
// One screen line per section, "[group][topic] key=value, ...", the [build] group first.
void report_build_and_runtime_info(const Manifest& manifest);
// keep_existing: append even on a fresh run, for a run about to fail.
void write_manifest(const Param& param, const Manifest& manifest,
                    bool keep_existing = false);
// Appends [runtime.end] when the time loop ends; a record without one did not end normally.
// manifest is the start record, written again first if the file went missing.
void write_manifest_end(const Param& param, const Manifest& manifest, const Variables& var,
                        int steps_this_run, int64_t wall_ns, int64_t init_ns,
                        int64_t compute_ns, double peak_rss_gib);
void report_mesh_info(const Variables& var, const char* tag);

// The team size a region gets now; CpuInfo's is one startup sample, which a dynamic team
// drifts from. -1 without OpenMP. Opens a parallel region: never call it from inside one.
int omp_team_size_now();

// Samples for the frame provenance: each degrades to "unknown" or 0 rather than fail.
std::string local_now();          // "2026-09-26T18:53:59-05:00 CDT": local, offset, zone
// This process's rss, sampled once and folded into the caller's running peak, floored at
// that sample: a record never shows a peak below its own rss.
double sample_rss_gib(double& peak_rss_gib);
double host_mem_avail_gib_now();  // what other processes have left
double process_cpu_time_sec();    // user + system, this process
double host_load_avg_1m();        // whole-machine contention
// Memory in use on the whole device, a neighbour's included; 0 unless on a GPU.
double device_mem_used_dev_gib(const DeviceInfo& dev);

#endif
