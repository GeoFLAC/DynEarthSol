#include <algorithm>
#include <cctype>
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
#endif

std::string read_cpu_model()
{
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
}

} // namespace

void report_cpu_runtime_status()
{
    const std::string cpu_model = read_cpu_model();
    const unsigned int hw_threads = std::thread::hardware_concurrency();
    std::cout << "[Runtime][CPU] model=" << cpu_model
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

void report_openacc_runtime_status()
{
#ifdef ACC
    const char* env_acc_device_type = std::getenv("ACC_DEVICE_TYPE");
    const char* env_nv_acc_device_type = std::getenv("NVCOMPILER_ACC_DEVICE_TYPE");
    const int n_nvidia = acc_get_num_devices(acc_device_nvidia);
    const int n_host = acc_get_num_devices(acc_device_host);

    acc_device_t forced = acc_device_default;
    const bool forced_by_env = env_forces_device(env_acc_device_type, forced) ||
                               env_forces_device(env_nv_acc_device_type, forced);

    if (forced_by_env) {
        acc_set_device_type(forced);
        acc_init(forced);
    }
    else if (n_nvidia > 0) {
        acc_set_device_type(acc_device_nvidia);
        acc_init(acc_device_nvidia);
    }
    else {
        acc_set_device_type(acc_device_host);
        acc_init(acc_device_host);
    }

    const acc_device_t active_type = acc_get_device_type();
    const int active_dev = acc_get_device_num(active_type);
    const bool using_gpu = (active_type == acc_device_nvidia && n_nvidia > 0);

    std::cout << "[Runtime][OpenACC] build=ACC"
              << ", env.ACC_DEVICE_TYPE="
              << (env_acc_device_type ? env_acc_device_type : "(unset)")
              << ", env.NVCOMPILER_ACC_DEVICE_TYPE="
              << (env_nv_acc_device_type ? env_nv_acc_device_type : "(unset)")
              << ", num_nvidia=" << n_nvidia
              << ", num_host=" << n_host
              << ", active_type=" << openacc_device_name(active_type)
              << ", active_dev=" << active_dev
              << '\n';

    const char* dev_name =
        acc_get_property_string(active_dev, active_type, acc_property_name);
    const char* dev_vendor =
        acc_get_property_string(active_dev, active_type, acc_property_vendor);
    const char* dev_driver =
        acc_get_property_string(active_dev, active_type, acc_property_driver);
    const size_t mem_total =
        acc_get_property(active_dev, active_type, acc_property_memory);
    const size_t mem_free =
        acc_get_property(active_dev, active_type, acc_property_free_memory);

    if (dev_name || dev_vendor || dev_driver || mem_total || mem_free) {
        std::cout << "[Runtime][OpenACC] device"
                  << " name=" << (dev_name ? dev_name : "(unknown)")
                  << ", vendor=" << (dev_vendor ? dev_vendor : "(unknown)")
                  << ", driver=" << (dev_driver ? dev_driver : "(unknown)");
        if (mem_total) {
            const double gb = static_cast<double>(mem_total) /
                              (1024.0 * 1024.0 * 1024.0);
            std::cout << ", mem_total_gb=" << gb;
        }
        if (mem_free) {
            const double gb = static_cast<double>(mem_free) /
                              (1024.0 * 1024.0 * 1024.0);
            std::cout << ", mem_free_gb=" << gb;
        }
        std::cout << '\n';
    }

    std::cout << "[Runtime][OpenACC] status="
              << (using_gpu ? "gpu-offload" : "cpu-fallback") << '\n';
#else
    std::cout << "[Runtime][OpenACC] build=non-ACC\n";
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
