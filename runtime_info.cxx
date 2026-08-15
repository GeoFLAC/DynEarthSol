#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>

#include "runtime_info.hpp"

#ifdef ACC
#include <openacc.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __APPLE__
#include <sys/sysctl.h>
#include <sys/types.h>
#endif

namespace {

#ifdef ACC
const char* openacc_device_name(acc_device_t t)
{
    switch (t) {
    case acc_device_nvidia: return "nvidia";
    case acc_device_host: return "host";
    case acc_device_none: return "none";
    case acc_device_default: return "default";
    case acc_device_not_host: return "not_host";
    default: return "unknown";
    }
}

acc_device_t parse_openacc_device_env(const char* env)
{
    if (!env) return acc_device_default;
    std::string s(env);
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (s == "nvidia" || s == "gpu") return acc_device_nvidia;
    if (s == "host" || s == "cpu") return acc_device_host;
    if (s == "none") return acc_device_none;
    if (s == "default") return acc_device_default;
    if (s == "not_host") return acc_device_not_host;
    return acc_device_default;
}

bool env_forces_device(const char* env, acc_device_t& device)
{
    if (!env) return false;
    device = parse_openacc_device_env(env);
    return (device != acc_device_default);
}

bool using_offload_device()
{
    return (acc_get_device_type() == acc_device_nvidia &&
            acc_get_num_devices(acc_device_nvidia) > 0);
}

const char* env_or_unset(const char* name)
{
    const char* v = std::getenv(name);
    return v ? v : "(unset)";
}
#endif

std::string read_cpu_model()
{
#ifdef __APPLE__
    // No /proc on macOS; the brand string lives in sysctl.
    std::size_t len = 0;
    if (sysctlbyname("machdep.cpu.brand_string", NULL, &len, NULL, 0) != 0 || len == 0)
        return "unknown";
    std::string model(len, '\0');
    if (sysctlbyname("machdep.cpu.brand_string", &model[0], &len, NULL, 0) != 0)
        return "unknown";
    // len counts the trailing NUL, which does not belong in a std::string.
    model.resize(len > 0 ? len - 1 : 0);
    return model.empty() ? "unknown" : model;
#else
    std::ifstream cpuinfo("/proc/cpuinfo");
    if (!cpuinfo) return "unknown";

    std::string line;
    while (std::getline(cpuinfo, line)) {
        const std::size_t pos = line.find("model name");
        if (pos == std::string::npos) continue;
        const std::size_t colon = line.find(':', pos);
        if (colon == std::string::npos) continue;
        std::string model = line.substr(colon + 1);
        while (!model.empty() && std::isspace(static_cast<unsigned char>(model.front())))
            model.erase(model.begin());
        if (!model.empty()) return model;
    }
    return "unknown";
#endif
}

} // namespace

void report_host_runtime_status()
{
    const std::string cpu_model = read_cpu_model();
    const unsigned int hw_threads = std::thread::hardware_concurrency();
    std::cout << "[Runtime][Host] name=" << cpu_model
              << ", hw_threads="
              << (hw_threads ? std::to_string(hw_threads) : "unknown");
#ifdef _OPENMP
    std::cout << ", openmp_max_threads=" << omp_get_max_threads()
              << ", openmp_num_procs=" << omp_get_num_procs()
              << ", openmp_dynamic=" << (omp_get_dynamic() ? "on" : "off");
#else
    std::cout << ", openmp=disabled";
#endif
    std::cout << '\n';
}

void init_offload_device()
{
#ifdef ACC
    acc_device_t forced = acc_device_default;
    const bool forced_by_env =
        env_forces_device(std::getenv("ACC_DEVICE_TYPE"), forced) ||
        env_forces_device(std::getenv("NVCOMPILER_ACC_DEVICE_TYPE"), forced);

    if (forced_by_env) {
        acc_set_device_type(forced);
        acc_init(forced);
    }
    else if (acc_get_num_devices(acc_device_nvidia) > 0) {
        acc_set_device_type(acc_device_nvidia);
        acc_init(acc_device_nvidia);
    }
    else {
        acc_set_device_type(acc_device_host);
        acc_init(acc_device_host);
    }
#endif
}

void report_device_runtime_status()
{
#ifndef ACC
    std::cout << "[Runtime][Device] build=non-ACC, kernel=CPU\n";
#else
    const acc_device_t active_type = acc_get_device_type();
    const int active_dev = acc_get_device_num(active_type);
    const bool using_device = using_offload_device();

    // The selection picture: the steering environment (these two can force the
    // device type), what was visible, and what was chosen.
    std::cout << "[Runtime][Device] build=ACC"
              << ", env.ACC_DEVICE_TYPE=" << env_or_unset("ACC_DEVICE_TYPE")
              << ", env.NVCOMPILER_ACC_DEVICE_TYPE="
              << env_or_unset("NVCOMPILER_ACC_DEVICE_TYPE");
    // When set it renumbers devices from 0, so num_nvidia and active_dev below
    // count within the remapped set, not the physical machine.
    const char* visible = std::getenv("CUDA_VISIBLE_DEVICES");
    if (visible)
        std::cout << ", env.CUDA_VISIBLE_DEVICES=" << visible;
    std::cout << ", num_nvidia=" << acc_get_num_devices(acc_device_nvidia)
              << ", num_host=" << acc_get_num_devices(acc_device_host)
              << ", active_type=" << openacc_device_name(active_type)
              << ", active_dev=" << active_dev
              << '\n';

    const char* name = acc_get_property_string(active_dev, active_type, acc_property_name);
    const char* vendor = acc_get_property_string(active_dev, active_type, acc_property_vendor);
    const char* driver = acc_get_property_string(active_dev, active_type, acc_property_driver);
    const size_t mem_total = acc_get_property(active_dev, active_type, acc_property_memory);
    // A snapshot, not a device property: whatever other processes hold on the device
    // right now, before DES allocates anything. Not comparable across runs.
    const size_t mem_free = acc_get_property(active_dev, active_type, acc_property_free_memory);

    if (name || vendor || driver || mem_total || mem_free) {
        std::cout << "[Runtime][Device] name=" << (name ? name : "(unknown)")
                  << ", vendor=" << (vendor ? vendor : "(unknown)")
                  << ", driver=" << (driver ? driver : "(unknown)");
        char mem[32];
        if (mem_total) {
            std::snprintf(mem, sizeof(mem), "%.1f",
                          static_cast<double>(mem_total) / (1024.0 * 1024.0 * 1024.0));
            std::cout << ", mem_total_gb=" << mem;
        }
        if (mem_free) {
            std::snprintf(mem, sizeof(mem), "%.1f",
                          static_cast<double>(mem_free) / (1024.0 * 1024.0 * 1024.0));
            std::cout << ", mem_free_at_start_gb=" << mem;
        }
        std::cout << '\n';
    }

    std::cout << "[Runtime][Device] kernel=" << (using_device ? "GPU" : "CPU")
              << ", status=" << (using_device ? "gpu-offload" : "cpu-fallback") << '\n';
#endif
}

void report_mesh_info(const Variables& var, const char* tag)
{
    std::cout << "[Runtime][Mesh] " << tag
              << " nodes=" << var.nnode
              << ", elements=" << var.nelem
              << ", segments=" << var.nseg
              << '\n';
}
