#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include "earthquake_state.hpp"

#include "constants.hpp"
#include "matprops.hpp"

namespace {

// Prevent aseismic creep from being misclassified as coseismic when plate rates are tiny.
constexpr double kMinEarthquakeSpeed = 1e-6; // m/s

double compute_max_plastic_strain_rate(const Variables& var)
{
    const double inv_dt = 1.0 / std::max(var.dt, 1e-30);
    double max_depls = 0.0;

#ifdef ACC
    #pragma acc parallel loop reduction(max:max_depls)
#else
    #pragma omp parallel for reduction(max:max_depls) default(none) shared(var)
#endif
    for (int e = 0; e < var.nelem; ++e) {
        const double depls = std::fabs((*var.delta_plstrain)[e]);
        max_depls = std::max(max_depls, depls);
    }
    return max_depls * inv_dt;
}

double compute_max_global_velocity_magnitude(const Variables& var)
{
    double vmax = 0.0;

#ifdef ACC
    #pragma acc parallel loop reduction(max:vmax)
#else
    #pragma omp parallel for reduction(max:vmax) default(none) shared(var)
#endif
    for (int n = 0; n < var.nnode; ++n) {
#ifdef THREED
        const double vx = (*var.vel)[n][0];
        const double vy = (*var.vel)[n][1];
        const double vz = (*var.vel)[n][2];
#else
        const double vx = (*var.vel)[n][0];
        const double vy = (*var.vel)[n][1];
        const double vz = 0.0;
#endif
        const double vmag = std::sqrt(vx * vx + vy * vy + vz * vz);
        vmax = std::max(vmax, vmag);
    }
    return vmax;
}

std::vector<double> compute_seismic_moment_rate_by_material(const Variables& var)
{
    const int nmat = var.mat->nmat;
    std::vector<double> moment_rate_by_mat(nmat, 0.0);
    if (nmat == 0) return moment_rate_by_mat;

    int_vec dominant_mat(var.nelem, -1);

#pragma omp parallel for default(none) shared(var, dominant_mat)
    for (int e = 0; e < var.nelem; ++e) {
        int_vec& a = (*var.elemmarkers)[e];
        dominant_mat[e] = std::distance(a.begin(), std::max_element(a.begin(), a.end()));
    }

    for (int e = 0; e < var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];

        double vx = 0.0, vy = 0.0, vz = 0.0;
        for (int j = 0; j < NODES_PER_ELEM; ++j) {
            vx += (*var.vel)[conn[j]][0];
            vy += (*var.vel)[conn[j]][1];
#ifdef THREED
            vz += (*var.vel)[conn[j]][2];
#endif
        }

        vx /= NODES_PER_ELEM;
        vy /= NODES_PER_ELEM;
#ifdef THREED
        vz /= NODES_PER_ELEM;
#else
        vz = 0.0;
#endif

        const double vmag = std::sqrt(vx * vx + vy * vy + vz * vz);
        const double volume = (*var.volume)[e];
        const double shearm = var.mat->shearm(e);

        if (dominant_mat[e] >= 0 && dominant_mat[e] < nmat) {
            moment_rate_by_mat[dominant_mat[e]] += shearm * volume * vmag;
        }
    }

    return moment_rate_by_mat;
}

} // namespace

void init_earthquake_state(const Param& param, EarthquakeState& state)
{
    state.in_earthquake_mode = false;
    state.allow_earthquake_output = false;
    state.last_output_step = 0;
    state.cumulative_moment_by_mat.assign(param.mat.nmat, 0.0);
}

void update_earthquake_tracking(const Param& param,
                                const Variables& var,
                                EarthquakeState& state)
{
    const bool has_rsf = (param.mat.rheol_type & MatProps::rh_rsf) != 0;
    if (!has_rsf) {
        state.in_earthquake_mode = false;
        state.allow_earthquake_output = false;
        return;
    }

    const double max_global_vel_mag = compute_max_global_velocity_magnitude(var);
    const double max_plastic_strain_rate = compute_max_plastic_strain_rate(var);
    const bool plastic_active = (max_plastic_strain_rate > 0.0);

    const double start_velocity_threshold =
        std::max(param.sim.earthquake_start_factor * var.max_vbc_val, kMinEarthquakeSpeed);
    const double end_velocity_threshold =
        std::max(param.sim.earthquake_end_factor * var.max_vbc_val, 0.5 * kMinEarthquakeSpeed);

    // Require both dynamic velocity and plastic activity to enter earthquake mode.
    // Exit when dynamic velocity relaxes OR plastic activity vanishes.
    const bool earthquake_now =
        (max_global_vel_mag > start_velocity_threshold) && plastic_active;
    const bool earthquake_end =
        (max_global_vel_mag < end_velocity_threshold) || !plastic_active;

    if (!state.in_earthquake_mode && earthquake_now) {
        state.in_earthquake_mode = true;
        // Do not output immediately on mode transition.
        state.last_output_step = var.steps;
        if (param.sim.seismic_moment_calculate_output) {
            std::fill(state.cumulative_moment_by_mat.begin(),
                      state.cumulative_moment_by_mat.end(),
                      0.0);
            std::ofstream ofs("seismic_moment_magnitude.txt", std::ios_base::app);
            if (ofs) ofs << "Earthquake event started at time: " << var.time << " s\n";
            else std::cerr << "Cannot open seismic_moment_magnitude.txt (start)\n";
        }
    } else if (state.in_earthquake_mode && earthquake_end) {
        state.in_earthquake_mode = false;
        if (param.sim.seismic_moment_calculate_output) {
            double total_moment = 0.0;
            for (double m : state.cumulative_moment_by_mat) total_moment += m;
            if (total_moment > 0.0) {
                const double mw = (2.0 / 3.0) * (std::log10(total_moment) - 9.1);
                std::cout << "M0=" << total_moment << ", Mw=" << mw << "\n";
            } else {
                std::cout << "Total seismic moment (M0) <= 0\n";
            }
            std::ofstream ofs("seismic_moment_magnitude.txt", std::ios_base::app);
            if (ofs) {
                ofs << "Earthquake event ended at time: " << var.time << " s\n";
                ofs << "Total seismic moment (M0): " << total_moment << "\n";
                if (total_moment > 0.0)
                    ofs << "Moment magnitude (Mw): " << (2.0 / 3.0) * (std::log10(total_moment) - 9.1) << "\n";
                ofs << "----------------------------------------\n";
            } else std::cerr << "Cannot open seismic_moment_magnitude.txt (end)\n";
        }
    }

    if (state.in_earthquake_mode && param.sim.seismic_moment_calculate_output) {
        std::vector<double> rate = compute_seismic_moment_rate_by_material(var);
        const size_t n = std::min(rate.size(), state.cumulative_moment_by_mat.size());
        for (size_t i = 0; i < n; ++i)
            state.cumulative_moment_by_mat[i] += rate[i] * var.dt;
    }

    state.allow_earthquake_output =
        (var.steps - state.last_output_step >= param.sim.earthquake_output_step_interval);
}
