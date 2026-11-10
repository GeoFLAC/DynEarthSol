#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

#include "point_sources.hpp"
#include "fields.hpp"
#include "utils.hpp"

namespace {

bool compute_point_shape_weights(const Variables &var, int e,
                                 const double q[NDIMS],
                                 double weights[NODES_PER_ELEM])
{
    if (!((*var.volume)[e] > 0 && std::isfinite((*var.volume)[e]))) return false;
    ConstConnAccessor conn = (*var.connectivity)[e];
    ConstArrayAccessor origin = (*var.coord)[conn[0]];
#ifdef THREED
    double dx[NODES_PER_ELEM], dy[NODES_PER_ELEM], dz[NODES_PER_ELEM];
    get_local_shape_fn(var, e, dx, dy, dz);
    const double *gradient[NDIMS] = {dx, dy, dz};
    const double tolerance = 5e-11;
#else
    double dx[NODES_PER_ELEM], dz[NODES_PER_ELEM];
    get_local_shape_fn(var, e, dx, dz);
    const double *gradient[NDIMS] = {dx, dz};
    const double tolerance = 1e-12;
#endif
    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        weights[i] = (i == 0);
        for (int d = 0; d < NDIMS; ++d)
            weights[i] += gradient[d][i] * (q[d] - origin[d]);
    }

    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        if (!std::isfinite(weights[i]) || weights[i] < -tolerance || weights[i] > 1.0 + tolerance)
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
                               double step_start, double step_end, int point_id, bool peak_rate)
{
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
    return peak_rate ? base_rate : base_rate * overlap / (step_end - step_start);
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
    message << "Source point " << point_id << " at (";
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

void deposit_point_sources(const Injection &injection, const Variables &var,
                           double start, double stop, PointSourceSink sink, void *context,
                           const double_vec *scale, bool peak_rate)
{
    // The previous substep can still be reading this derived source array.
    // Finish device work before any host rewrite or current-mesh lookup.
    #pragma acc wait
    if (!injection.enabled) return;
    if (!std::isfinite(stop - start) || !(stop > start)) {
        die(EXIT_RUNTIME_NAN,
            "Point-source integration interval must be finite and positive.");
    }

    for (int p = 0; p < injection.num_points; ++p) {
        double q[NDIMS];
        q[0] = injection.points_x[p];
#ifdef THREED
        q[1] = injection.points_y[p];
#endif
        q[NDIMS - 1] = injection.points_z[p];

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

        const double source_rate = (scale ? (*scale)[p] : 1.0) *
            injection_rate_for_step(injection, start, stop, p, peak_rate);
        if (!std::isfinite(source_rate)) {
            point_source_error(EXIT_RUNTIME_NAN, p, q, var,
                               "produced a non-finite rate");
        }
        if (source_rate == 0.0) continue;

        ConstConnAccessor conn = (*var.connectivity)[mapped_element];
        long double distributed_rate = 0.0L;
        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            const double contribution = source_rate * mapped_weights[i];
            sink(p,conn[i],contribution,context);
            distributed_rate += static_cast<long double>(contribution);
        }

        const long double point_rate = static_cast<long double>(source_rate);
        if (!source_totals_match(distributed_rate, point_rate,
                                 std::abs(point_rate), NODES_PER_ELEM)) {
            point_source_error(EXIT_INTERNAL_ASSERT, p, q, var,
                               "did not conserve its requested nodal source total");
        }
    }

}

namespace {
void add_point_source(int, int node, double rate, void *context)
{
    auto &source = *static_cast<double_vec *>(context);
    source[node] += rate;
    if (!std::isfinite(source[node]))
        die(EXIT_RUNTIME_NAN,"Point source overflowed a nodal source value.");
}
} // namespace

void assemble_point_sources(const Injection &injection, const Variables &var,
                           double start, double stop, double_vec &source,
                           const double_vec *scale, bool peak_rate)
{
    #pragma acc wait
    std::fill(source.begin(),source.end(),0.0);
    if (!injection.enabled) return;
    if (source.size() != static_cast<std::size_t>(var.nnode))
        die(EXIT_INTERNAL_ASSERT,"Point-source vector size does not match the current mesh node count.");
    deposit_point_sources(injection,var,start,stop,add_point_source,&source,scale,peak_rate);
    long double assembled_total = 0, requested_total = 0, absolute_rate_total = 0;
    std::size_t active_contributions = 0;
    for (int p = 0; p < injection.num_points; ++p) {
        const double rate = (scale ? (*scale)[p] : 1.0)*
            injection_rate_for_step(injection,start,stop,p,peak_rate);
        requested_total += rate;
        absolute_rate_total += std::abs(rate);
        if (rate != 0) active_contributions += NODES_PER_ELEM;
    }
    for (double nodal_source : source) {
        if (!std::isfinite(nodal_source)) {
            die(EXIT_RUNTIME_NAN,
                "Point-source assembly produced a non-finite nodal value.");
        }
        assembled_total += static_cast<long double>(nodal_source);
    }
    if (!source_totals_match(assembled_total, requested_total,
                             absolute_rate_total, active_contributions)) {
        die(EXIT_INTERNAL_ASSERT,
            "Point-source nodal total does not equal the requested active rate.");
    }
}
