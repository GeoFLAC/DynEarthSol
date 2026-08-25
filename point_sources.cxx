#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

#include "point_sources.hpp"
#include "utils.hpp"

namespace {

#ifdef THREED
double det3(const double a[3], const double b[3], const double c[3])
{
    return a[0] * (b[1] * c[2] - b[2] * c[1]) -
           a[1] * (b[0] * c[2] - b[2] * c[0]) +
           a[2] * (b[0] * c[1] - b[1] * c[0]);
}
#endif

bool compute_point_shape_weights(const Variables& var, int e,
                                 const double q[NDIMS],
                                 double weights[NODES_PER_ELEM])
{
    ConstConnAccessor conn = (*var.connectivity)[e];
#ifdef THREED
    ConstArrayAccessor a = (*var.coord)[conn[0]];
    ConstArrayAccessor b = (*var.coord)[conn[1]];
    ConstArrayAccessor c = (*var.coord)[conn[2]];
    ConstArrayAccessor d = (*var.coord)[conn[3]];
    double ba[3], ca[3], da[3], qa[3];
    for (int i = 0; i < 3; ++i) {
        ba[i] = b[i] - a[i];
        ca[i] = c[i] - a[i];
        da[i] = d[i] - a[i];
        qa[i] = q[i] - a[i];
    }

    const double det = det3(ba, ca, da);
    if (std::abs(det) <= 0.0) return false;
    weights[1] = det3(qa, ca, da) / det;
    weights[2] = det3(ba, qa, da) / det;
    weights[3] = det3(ba, ca, qa) / det;
    weights[0] = 1.0 - weights[1] - weights[2] - weights[3];
    const double tolerance = 5e-11;
#else
    ConstArrayAccessor a = (*var.coord)[conn[0]];
    ConstArrayAccessor b = (*var.coord)[conn[1]];
    ConstArrayAccessor c = (*var.coord)[conn[2]];
    const double det = (b[1] - c[1]) * (a[0] - c[0]) +
                       (c[0] - b[0]) * (a[1] - c[1]);
    if (std::abs(det) <= 0.0) return false;
    weights[0] = ((b[1] - c[1]) * (q[0] - c[0]) +
                  (c[0] - b[0]) * (q[1] - c[1])) / det;
    weights[1] = ((c[1] - a[1]) * (q[0] - c[0]) +
                  (a[0] - c[0]) * (q[1] - c[1])) / det;
    weights[2] = 1.0 - weights[0] - weights[1];
    const double tolerance = 1e-12;
#endif

    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        if (weights[i] < -tolerance || weights[i] > 1.0 + tolerance)
            return false;
    }

    double sum = 0.0;
    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        weights[i] = std::min(1.0, std::max(0.0, weights[i]));
        sum += weights[i];
    }
    if (!(sum > 0.0)) return false;
    for (int i = 0; i < NODES_PER_ELEM; ++i)
        weights[i] /= sum;

    return true;
}

double injection_rate_for_step(const Injection& injection,
                               const Variables& var, int point_id)
{
    const double step_start = var.time - var.dt;
    const double step_end = var.time;
    const double source_start = injection.start_time[point_id];
    const double source_end = injection.end_time[point_id];
    const double overlap =
        std::min(step_end, source_end) - std::max(step_start, source_start);
    if (!(overlap > 0.0)) return 0.0;

    double base_rate = injection.rate[point_id];
    if (injection.rate_model == injection_rate_total_amount) {
        base_rate = injection.total_amount[point_id] /
                    (source_end - source_start);
    }
    return base_rate * overlap / var.dt;
}

bool source_totals_match(long double actual, long double expected,
                         long double absolute_rate_total,
                         std::size_t accumulation_terms)
{
    const long double scale =
        std::max(absolute_rate_total,
                 std::max(std::abs(actual), std::abs(expected)));
    if (scale == 0.0L) return actual == expected;

    const std::size_t terms = std::max<std::size_t>(1, accumulation_terms);
    const long double tolerance =
        32.0L * terms * std::numeric_limits<double>::epsilon() * scale;
    return std::abs(actual - expected) <= tolerance;
}

[[noreturn]] void point_source_error(ExitCode code, int point_id,
                                     const double q[NDIMS],
                                     const Variables& var, const char* reason)
{
    std::ostringstream message;
    message << "Fluid source point " << point_id << " at (";
    for (int d = 0; d < NDIMS; ++d) {
        if (d != 0) message << ", ";
        message << q[d];
    }
    message << ") " << reason << " at step " << var.steps
            << " and time " << var.time << " s.";
    const std::string text = message.str();
    die(code, text.c_str());
}

} // namespace

void assemble_fluid_point_sources(const Param& param, const Variables& var,
                                  double_vec& fluid_source)
{
    std::fill(fluid_source.begin(), fluid_source.end(), 0.0);
    if (!param.injection.enabled) return;
    if (!std::isfinite(var.dt) || !(var.dt > 0.0)) {
        die(EXIT_RUNTIME_NAN,
            "Fluid point-source hydraulic integration interval must be finite and positive.");
    }
    if (fluid_source.size() != static_cast<std::size_t>(var.nnode)) {
        die(EXIT_INTERNAL_ASSERT,
            "Fluid point-source vector size does not match the current mesh node count.");
    }

    // Coordinates and connectivity can be produced by asynchronous device
    // work. The point lookup below is deliberately host-side and authoritative
    // for the current mesh, so complete that work before reading either array.
    #pragma acc wait

    long double requested_total = 0.0L;
    long double absolute_rate_total = 0.0L;
    std::size_t active_contributions = 0;

    for (int p = 0; p < param.injection.num_points; ++p) {
        double q[NDIMS];
        q[0] = param.injection.points_x[p];
#ifdef THREED
        q[1] = param.injection.points_y[p];
        q[2] = param.injection.points_z[p];
#else
        q[1] = param.injection.points_z[p];
#endif

        for (int d = 0; d < NDIMS; ++d) {
            if (!std::isfinite(q[d])) {
                point_source_error(EXIT_RUNTIME_NAN, p, q, var,
                                   "has a non-finite coordinate");
            }
        }

        int mapped_element = -1;
        double mapped_weights[NODES_PER_ELEM];
        for (int e = 0; e < var.nelem; ++e) {
            double weights[NODES_PER_ELEM];
            if (!compute_point_shape_weights(var, e, q, weights)) continue;

            // Shared faces can have several candidates. Element order is
            // stable, so selecting the first maps the requested rate once and
            // makes ownership deterministic across thread counts.
            mapped_element = e;
            for (int i = 0; i < NODES_PER_ELEM; ++i)
                mapped_weights[i] = weights[i];
            break;
        }
        if (mapped_element < 0) {
            point_source_error(EXIT_RUNTIME_LOOKUP, p, q, var,
                               "could not be mapped into the current mesh");
        }

        const double source_rate =
            injection_rate_for_step(param.injection, var, p);
        if (!std::isfinite(source_rate)) {
            point_source_error(EXIT_RUNTIME_NAN, p, q, var,
                               "produced a non-finite rate");
        }
        if (source_rate == 0.0) continue;

        ConstConnAccessor conn = (*var.connectivity)[mapped_element];
        long double distributed_rate = 0.0L;
        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            const double contribution = source_rate * mapped_weights[i];
            fluid_source[conn[i]] += contribution;
            if (!std::isfinite(fluid_source[conn[i]])) {
                point_source_error(EXIT_RUNTIME_NAN, p, q, var,
                                   "overflowed a nodal source value");
            }
            distributed_rate += static_cast<long double>(contribution);
        }

        const long double point_rate = static_cast<long double>(source_rate);
        if (!source_totals_match(distributed_rate, point_rate,
                                 std::abs(point_rate), NODES_PER_ELEM)) {
            point_source_error(EXIT_INTERNAL_ASSERT, p, q, var,
                               "did not conserve its requested nodal source total");
        }
        requested_total += point_rate;
        absolute_rate_total += std::abs(point_rate);
        active_contributions += NODES_PER_ELEM;
    }

    long double assembled_total = 0.0L;
    for (double nodal_source : fluid_source) {
        if (!std::isfinite(nodal_source)) {
            die(EXIT_RUNTIME_NAN,
                "Fluid point-source assembly produced a non-finite nodal value.");
        }
        assembled_total += static_cast<long double>(nodal_source);
    }
    if (!source_totals_match(assembled_total, requested_total,
                             absolute_rate_total, active_contributions)) {
        die(EXIT_INTERNAL_ASSERT,
            "Fluid point-source nodal total does not equal the requested active rate.");
    }
}
