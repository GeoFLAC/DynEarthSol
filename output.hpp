#ifndef DYNEARTHSOL3D_OUTPUT_HPP
#define DYNEARTHSOL3D_OUTPUT_HPP

#include "array2d.hpp"
#include "runtime_info.hpp"

class Output
{
private:
    const std::string &modelname;
    // "no", or "<model>:<frame>" -- the frame's own answer to "was this a restart?"
    const std::string restart_from;
    // Copied, unlike modelname: it spares callers from keeping the originals alive.
    const BuildInfo build;
    const CpuInfo cpu;
    const DeviceInfo dev;
    // The run's high-water resident set, remembered across writes: each sample is only a
    // lower bound (see host_peak_rss_gib), so the maximum has to live somewhere. Per
    // PROCESS on purpose -- a restart leg starts its own, since carrying a predecessor's
    // value forward would report a peak this address space never reached. The legs chain
    // through restart_from instead.
    double peak_rss_gib;
    const int64_t start_time;
    const bool is_averaged;
    const int average_interval;
    const bool has_marker_output;
    const int hdf5_compression_level;
    const bool may_overwrite_;
    const int start_frame_;
    int frame;
    int64_t run_time_ns;

    // stuffs for averging fields
    double time0;
    array_t coord0;
    tensor_t strain0;
    tensor_t stress_avg;
    double_vec delta_plstrain_avg;

    void write_info(const Variables& var, double dt);
    void _write(const Variables& var, bool disable_averaging=false);

public:
    Output(const Param& param, const BuildInfo& build, const CpuInfo& cpu,
           const DeviceInfo& dev, int64_t start_time, int start_frame);
    ~Output();
    void write(Variables& var);
    void write_exact(Variables& var);
    void write_exact_error(const Variables& var);
    void write_checkpoint(const Param& param, const Variables& var);
    void average_fields(Variables& var);

};


#endif
