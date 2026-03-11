#ifndef DYNEARTHSOL3D_ATS_OUTPUT_SCHEDULER_HPP
#define DYNEARTHSOL3D_ATS_OUTPUT_SCHEDULER_HPP

#include "earthquake_state.hpp"

class Output;

void handle_ats_output(const Param& param,
                       Variables& var,
                       Output& output,
                       EarthquakeState& state,
                       double starting_time,
                       double starting_step,
                       int& next_regular_frame);

#endif
