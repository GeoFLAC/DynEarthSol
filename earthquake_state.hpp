#ifndef DYNEARTHSOL3D_EARTHQUAKE_STATE_HPP
#define DYNEARTHSOL3D_EARTHQUAKE_STATE_HPP

#include <vector>

#include "parameters.hpp"

struct EarthquakeState {
    bool in_earthquake_mode = false;
    bool allow_earthquake_output = false;
    int last_output_step = 0;
    std::vector<double> cumulative_moment_by_mat;
};

void init_earthquake_state(const Param& param, EarthquakeState& state);

void update_earthquake_tracking(const Param& param,
                                const Variables& var,
                                EarthquakeState& state);

#endif
