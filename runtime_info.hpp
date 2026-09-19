#ifndef DYNEARTHSOL_RUNTIME_INFO_HPP
#define DYNEARTHSOL_RUNTIME_INFO_HPP

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
    double mem_avail_gib;    // at probe time: what other processes left, not a host property
    bool openmp_enabled;   // false in an openmp=0 build; the omp_* fields are -1
    int omp_team_size;     // threads the runtime actually hands out (measured)
    int omp_max_threads;
    int omp_num_procs;     // processors available to THIS process: affinity- and cgroup-scoped
    bool omp_dynamic;
    // Who supplied the wait policy: "env" (the operator set OMP_WAIT_POLICY or
    // KMP_BLOCKTIME), "des-default" (main()'s macOS default), or "runtime" (neither, so
    // the OpenMP runtime's own default governs).
    std::string wait_policy_src;
};

// What this binary was built as: the embedded build.snapshot code/host groups as
// fields, for the consumers that write them one at a time (the HDF5 frame writer).
// The build CONFIGURATION is deliberately absent -- the block's own
// options/toolchain/defines lines are its single source of truth.
struct BuildInfo {
    std::string rev;       // git describe against version tags
    std::string branch;
    std::string dirty;     // file/±line counts ("2f/+13/-12"), or "clean"
    std::string origin;    // credential-stripped remote URL
    std::string state_utc; // when the build's source/configuration state last changed
    std::string build_os;  // the build host's OS
    std::string builder;   // user@host that built the binary
    // The executable FILE's mtime, so a cp without -p or a touch moves it: near enough to
    // the link time to spot a swapped binary, not an authenticated build time.
    std::string exe_mtime_utc;
};

// The offload device this process actually got, which is not a property of the build:
// an openacc=1 binary falls back to the host when no NVIDIA device is visible, and the
// two do not give the same numbers. Empty strings and non-positive numbers = unknown.
struct DeviceInfo {
    bool acc_build;        // compiled with -acc
    bool using_gpu;        // an NVIDIA device was selected AND initialized
    std::string kernel;    // "CPU" or "GPU" -- where this run's arithmetic happened
    std::string name;      // e.g. "NVIDIA A100-SXM4-40GB"
    std::string vendor;
    std::string cuda_driver;   // the CUDA driver API version, not nvidia-smi's driver_version
    // CUDA_VISIBLE_DEVICES as it stood when the device was chosen: "all" unset, "none"
    // set but empty (which is what hides every device), else the value itself.
    std::string visible_devices;
    int num_devices;       // NVIDIA devices visible to this process
    int active_dev;        // index among them, -1 when running on the host
    double mem_total_gib;
    // Genuinely free device memory, unlike the host's mem_avail_gib estimate -- but
    // sampled at selection time, AFTER acc_init, so our own context (~230 MiB) is
    // already allocated and this is not a clean "before this run" baseline.
    double mem_free_gib;
};

BuildInfo probe_build_info();
// wait_policy_from_des: main() set OMP_WAIT_POLICY itself. Only main() can know -- by the
// time this runs, getenv cannot tell its own setenv from the operator's export.
CpuInfo probe_cpu_info(bool wait_policy_from_des);
// Selects and initializes the OpenACC device, then describes what it selected. The
// acc_init() is the job, not a side effect, so call it exactly once and before any
// compute; a non-ACC build only fills in the CPU answer.
DeviceInfo init_offload_device();
// One "[section]" of the .manifest: key/value pairs in file order. The screen
// report prints these same sections, so file and screen cannot drift.
struct ManifestSection {
    std::string name;
    std::vector<std::pair<std::string, std::string> > fields;
};
typedef std::vector<ManifestSection> Manifest;

// Everything the .manifest and the start-up screen report say about this run: the
// [runtime.*] sections measured here, then the executable's build.snapshot block as
// [build.*] sections. Composed once, so file and screen carry the same start_utc.
Manifest compose_manifest(const Param& param, const BuildInfo& build,
                          const CpuInfo& cpu, const DeviceInfo& dev);
// One screen line per section, "[group][topic] key=value, ...", the [build] group
// first and then [runtime]. Gated by the caller on sim.has_runtime_info_display.
void report_build_and_runtime_info(const Manifest& manifest);
// keep_existing: append even on a fresh run, for a run that is about to fail -- truncating
// would replace a completed run's record with the record of one that never started.
void write_manifest(const Param& param, const Manifest& manifest,
                    bool keep_existing = false);
void report_mesh_info(const Variables& var, const char* tag);

// The team size a region gets NOW: CpuInfo's is one startup sample, and libgomp's dynamic
// adjustment moves the team between regions (measured 55 -> 54 under load). nvomp sets
// OMP_DYNAMIC on by default but was measured never to shrink a team, so there it changes
// nothing. -1 in an openmp=0 build. Opens a parallel region, so never call it from inside
// one: it would measure a nested team of 1.
int omp_team_size_now();

// Status-at-write samples for the per-frame provenance attributes. All degrade to
// "unknown" / 0 rather than fail -- metadata must never abort a run.
std::string utc_now();
double host_mem_rss_gib();                          // this process's resident set
// A LOWER BOUND on the peak right now: Linux ru_maxrss trails the live rss, so this is
// max(the two) and can fall when the rss does. Callers must keep the running maximum.
double host_peak_rss_gib();
double host_mem_avail_gib_now();                     // what other processes have left
double process_cpu_time_sec();                     // user + system, this process
double host_load_avg_1m();                         // whole-machine contention
// Memory in use on the whole DEVICE, not by this process: the OpenACC runtime reports
// only free memory, so a neighbour's allocation counts here too. 0 unless on a GPU.
double device_mem_used_dev_gib(const DeviceInfo& dev);

#endif
