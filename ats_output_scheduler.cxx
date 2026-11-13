#include <algorithm>
#include <cmath>
#include <limits>

#include "ats_output_scheduler.hpp"

#include "constants.hpp"
#include "output.hpp"

bool output_allowed(const Param& param, const Variables& var)
{
    return !param.sim.is_outputting_averaged_fields ||
        var.steps % param.mesh.quality_check_step_interval == 0;
}

bool regular_output_is_due(const Param& param, const Variables& var, int next_regular_frame,
                           bool catch_up)
{
    if (!output_allowed(param, var)) return false;

    const int64_t step_target = static_cast<int64_t>(next_regular_frame) * param.sim.output_step_interval;
    const int64_t steps = var.steps - var.output_start_step;
    const bool step_due = param.sim.output_step_interval != std::numeric_limits<int>::max() &&
        (catch_up ? steps >= step_target : steps == step_target);
    if (step_due) return true;
    if (param.sim.output_time_interval_in_yr == std::numeric_limits<double>::max()) return false;
    const double time_target = next_regular_frame * param.sim.output_time_interval_in_yr * YEAR2SEC;
    const double elapsed = var.time - var.output_start_time;
    return catch_up ? elapsed >= time_target : elapsed > time_target;
}

void handle_ats_output(const Param& param,
                       Variables& var,
                       Output& output,
                       EarthquakeState& state,
                       int& next_regular_frame)
{
    const bool regular_output_due = regular_output_is_due(param, var, next_regular_frame, true);

    int following_regular_frame = next_regular_frame;
    if (regular_output_due) {
        // Catch-up logic for mixed step/time schedules.
        int frames_due_step = 0;
        if (param.sim.output_step_interval != std::numeric_limits<int>::max()) {
            const int64_t steps_since = var.steps - var.output_start_step;
            frames_due_step = static_cast<int>(steps_since / param.sim.output_step_interval);
        }

        int frames_due_time = 0;
        if (param.sim.output_time_interval_in_yr != std::numeric_limits<double>::max()) {
            const double elapsed_years = (var.time - var.output_start_time) / YEAR2SEC;
            frames_due_time = static_cast<int>(std::floor(elapsed_years / param.sim.output_time_interval_in_yr));
        }

        // The due test compares seconds, the floor above years; rounding can
        // leave the floor one short, so the frame written now always counts.
        const int frames_due = std::max({frames_due_step, frames_due_time, next_regular_frame});
        following_regular_frame = frames_due + 1;
    }
    // Checkpoints below must carry the index for the next physical step.
    var.output_next_regular_frame = following_regular_frame;

    bool wrote_earthquake_output = false;

    // (1) Earthquake-triggered output first.
    if (state.in_earthquake_mode &&
        state.allow_earthquake_output &&
        output_allowed(param, var)) {

        if (next_regular_frame % param.sim.checkpoint_frame_interval == 0)
            output.write_checkpoint(param, var);

        output.write(var);

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

            output.write(var);
        }

        next_regular_frame = following_regular_frame;
    }
}
