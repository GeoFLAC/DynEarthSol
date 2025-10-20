#ifndef DYNEARTHSOL3D_OUTPUT_HPP
#define DYNEARTHSOL3D_OUTPUT_HPP

#include "array2d.hpp"

class Output
{
private:
    const std::string &modelname;
    const int64_t start_time;
    const bool is_averaged;
    const int average_interval;
    const bool has_marker_output;
    const int hdf5_compression_level;
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
    Output(const Param& param, int64_t start_time, int start_frame);
    ~Output();
    void write(Variables& var);
    void write_exact(Variables& var);
    void write_checkpoint(const Param& param, const Variables& var);
    void average_fields(Variables& var);

};


#endif
