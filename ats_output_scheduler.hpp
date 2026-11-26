#ifndef DYNEARTHSOL3D_ATS_OUTPUT_SCHEDULER_HPP
#define DYNEARTHSOL3D_ATS_OUTPUT_SCHEDULER_HPP

#include "earthquake_state.hpp"

class Output;

// Averaged fields can be saved only at a mesh-quality step.
bool output_allowed(const Param& param, const Variables& var);
// catch_up (GVS): due once the next frame's step or time is reached or passed.
// Otherwise due at the exact step or at the first time past it.
bool regular_output_is_due(const Param& param, const Variables& var, int next_regular_frame,
                           bool catch_up);

void handle_ats_output(const Param& param,
                       Variables& var,
                       Output& output,
                       EarthquakeState& state,
                       int& next_regular_frame);

#endif
