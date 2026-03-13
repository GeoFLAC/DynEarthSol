#include <algorithm>
#include <cmath>
#include <limits>

#include "ats_output_scheduler.hpp"

#include "constants.hpp"
#include "output.hpp"
#include "utils.hpp"

void handle_ats_output(const Param& param,
                       Variables& var,
                       Output& output,
                       EarthquakeState& state,
                       double starting_time,
                       double starting_step,
                       int& next_regular_frame)
{
    const bool output_allowed =
        (!param.sim.is_outputting_averaged_fields ||
         (var.steps % param.mesh.quality_check_step_interval == 0));

    const bool regular_output_due =
        (((param.sim.output_step_interval != std::numeric_limits<int>::max() &&
           (var.steps - starting_step) >= next_regular_frame * param.sim.output_step_interval) ||
          (param.sim.output_time_interval_in_yr != std::numeric_limits<double>::max() &&
           (var.time - starting_time) >= next_regular_frame * param.sim.output_time_interval_in_yr * YEAR2SEC)) &&
         output_allowed);

    bool wrote_earthquake_output = false;

    // (1) Earthquake-triggered output first.
    if (state.in_earthquake_mode &&
        state.allow_earthquake_output &&
        output_allowed) {

        if (next_regular_frame % param.sim.checkpoint_frame_interval == 0)
            output.write_checkpoint(param, var);

        const int64_t t0 = get_nanoseconds();
        output.write(var);
        var.func_time.output_time += get_nanoseconds() - t0;

        state.last_output_step = var.steps;
        wrote_earthquake_output = true;
    }

    // (2) Regular (interseismic) output.
    if (regular_output_due) {
        const bool regular_frame_satisfied_by_earthquake_output =
            wrote_earthquake_output && state.last_output_step == var.steps;

        if (!regular_frame_satisfied_by_earthquake_output) {
            if (next_regular_frame % param.sim.checkpoint_frame_interval == 0)
                output.write_checkpoint(param, var);

            const int64_t t0 = get_nanoseconds();
            output.write(var);
            var.func_time.output_time += get_nanoseconds() - t0;
        }

        // Catch-up logic for mixed step/time schedules.
        int frames_due_step = 0;
        if (param.sim.output_step_interval != std::numeric_limits<int>::max()) {
            const int64_t steps_since = var.steps - starting_step;
            frames_due_step = static_cast<int>(steps_since / param.sim.output_step_interval);
        }

        int frames_due_time = 0;
        if (param.sim.output_time_interval_in_yr != std::numeric_limits<double>::max()) {
            const double elapsed_years = (var.time - starting_time) / YEAR2SEC;
            frames_due_time = static_cast<int>(std::floor(elapsed_years / param.sim.output_time_interval_in_yr));
        }

        const int frames_due = std::max(frames_due_step, frames_due_time);
        next_regular_frame = frames_due + 1;
    }
}
