#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

#include "point_sources.hpp"
#include "hydraulics.hpp"
#include "mpi_simulation.hpp"
#include "fields.hpp"
#include "geometry.hpp"
#include "utils.hpp"

namespace {

bool compute_point_shape_weights(const Variables &var, int e,
                                 const double q[NDIMS],
                                 double weights[NODES_PER_ELEM])
{
    if (!((*var.volume)[e] > 0 && std::isfinite((*var.volume)[e]))) return false;
#ifdef THREED
    const double tolerance = 5e-11;
#else
    const double tolerance = 1e-12;
#endif
    ConstArrayAccessor origin = (*var.coord)[(*var.connectivity)[e][0]];
    double gradient[NDIMS][NODES_PER_ELEM];
    get_shape_gradients(var, e, gradient);
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

void point_position(const Injection &injection, int p, double q[NDIMS])
{
    q[0] = injection.points_x[p];
#ifdef THREED
    q[1] = injection.points_y[p];
#endif
    q[NDIMS - 1] = injection.points_z[p];
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

// Current positions of the material points that owners hold.
void owner_positions(const Variables &var, const PointSourceOwners &owners, double_vec &positions)
{
    positions.resize(owners.element.size()*NDIMS);
    for (std::size_t p = 0; p < owners.element.size(); ++p)
        barycentric_point(var, owners.element[p], &owners.weight[p*NODES_PER_ELEM], &positions[p*NDIMS]);
}

} // namespace

void place_point_sources(const Injection &injection, const Variables &var,
                         PointSourceOwners &owners, const double_vec *positions)
{
    owners.element.clear();
    owners.weight.clear();
    if (!injection.enabled) return;
    owners.element.assign(injection.num_points, -1);
    owners.weight.assign(std::size_t(injection.num_points)*NODES_PER_ELEM, 0.0);
    for (int p = 0; p < injection.num_points; ++p) {
        double q[NDIMS];
        if (positions) std::copy_n(&(*positions)[p*NDIMS], NDIMS, q);
        else point_position(injection, p, q);
        for (int d = 0; d < NDIMS; ++d) {
            if (!std::isfinite(q[d]))
                point_source_error(EXIT_RUNTIME_NAN, p, q, var, "has a non-finite coordinate");
        }
        for (int e = 0; e < var.nelem; ++e) {
            double weights[NODES_PER_ELEM];
            if (!compute_point_shape_weights(var, e, q, weights)) continue;
            // Shared faces can have several candidates. Element order is
            // stable, so selecting the first maps the requested rate once.
            owners.element[p] = e;
            std::copy_n(weights, NODES_PER_ELEM, &owners.weight[p*NODES_PER_ELEM]);
            break;
        }
        if (owners.element[p] < 0)
            point_source_error(EXIT_RUNTIME_LOOKUP, p, q, var,
                               "could not be mapped into the current mesh");
    }
}

SourcePositions source_positions(const Variables &var)
{
    #pragma acc wait
    SourcePositions at;
    owner_positions(var, var.fluid_points, at.fluid);
    owner_positions(var, var.heat_points, at.heat);
    well_screen_positions(var, at.wells);
    return at;
}

void place_sources(const Param &param, Variables &var, const SourcePositions *at)
{
    #pragma acc wait
    place_point_sources(param.injection, var, var.fluid_points, at ? &at->fluid : nullptr);
    place_point_sources(param.heat_injection, var, var.heat_points, at ? &at->heat : nullptr);
    place_well_screens(param, var, at ? &at->wells : nullptr);
    assemble_well_source(param, var);
}

void deposit_point_sources(const Injection &injection, const PointSourceOwners &owners,
                           const Variables &var, double start, double stop,
                           PointSourceSink sink, void *context,
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
    if (owners.element.size() != std::size_t(injection.num_points))
        die(EXIT_INTERNAL_ASSERT, "Point sources were not placed in the mesh.");

    for (int p = 0; p < injection.num_points; ++p) {
        // The owner is an element of the whole mesh; a compute rank holds it
        // as an owned or ghost element, or not at all.
        const int mapped_element = mpi_local_element(var, owners.element[p]);
        const double *mapped_weights = &owners.weight[p*NODES_PER_ELEM];
        if (mapped_element < 0) continue;  // on another rank
        double q[NDIMS];
        point_position(injection, p, q);

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
struct PointSourceAssembly {
    double_vec &source;
    double_vec *weighted;
    const double_vec *point_weight;
};

void add_point_source(int point, int node, double rate, void *context)
{
    auto &assembly = *static_cast<PointSourceAssembly *>(context);
    assembly.source[node] += rate;
    if (!std::isfinite(assembly.source[node]))
        die(EXIT_RUNTIME_NAN,"Point source overflowed a nodal source value.");
    if (assembly.weighted) {
        auto &value = (*assembly.weighted)[node];
        value += rate * (*assembly.point_weight)[point];
        if (!std::isfinite(value))
            die(EXIT_RUNTIME_NAN,"Point source overflowed a weighted nodal source value.");
    }
}
} // namespace

void assemble_point_sources(const Injection &injection, const PointSourceOwners &owners,
                            const Variables &var, double start, double stop, double_vec &source,
                            const double_vec *scale, bool peak_rate,
                            double_vec *weighted, const double_vec *point_weight)
{
    #pragma acc wait
    if ((weighted == nullptr) != (point_weight == nullptr) || weighted == &source ||
        (weighted && (weighted->size() != source.size() ||
                      point_weight->size() != static_cast<std::size_t>(injection.num_points))))
        die(EXIT_INTERNAL_ASSERT,"Weighted point source requires distinct matching output and point-weight arrays.");
    std::fill(source.begin(),source.end(),0.0);
    if (weighted) std::fill(weighted->begin(),weighted->end(),0.0);
    if (!injection.enabled) return;
    if (source.size() != static_cast<std::size_t>(var.nnode))
        die(EXIT_INTERNAL_ASSERT,"Point-source vector size does not match the current mesh node count.");
    PointSourceAssembly assembly{source,weighted,point_weight};
    deposit_point_sources(injection,owners,var,start,stop,add_point_source,&assembly,scale,peak_rate);
    for (int field = 0; field < (weighted ? 2 : 1); ++field) {
        long double assembled_total = 0, requested_total = 0, absolute_rate_total = 0;
        std::size_t active_contributions = 0;
        for (int p = 0; p < injection.num_points; ++p) {
            double rate = (scale ? (*scale)[p] : 1.0)*
                injection_rate_for_step(injection,start,stop,p,peak_rate);
            if (field) rate *= (*point_weight)[p];
            requested_total += rate;
            absolute_rate_total += std::abs(rate);
            if (rate != 0) active_contributions += NODES_PER_ELEM;
        }
        // Owned nodes only: each node is counted once over the ranks.
        const double_vec &nodal = field ? *weighted : source;
        for (int n = 0; n < mpi_node_end(var); ++n) {
            if (!std::isfinite(nodal[n])) {
                die(EXIT_RUNTIME_NAN,
                    "Point-source assembly produced a non-finite nodal value.");
            }
            assembled_total += static_cast<long double>(nodal[n]);
        }
#ifdef USEMPI
        if (var.mpi && var.mpi->in_kernel())
            assembled_total = var.mpi->sum(static_cast<double>(assembled_total));
#endif
        if (!source_totals_match(assembled_total, requested_total,
                                 absolute_rate_total, active_contributions)) {
            die_together(EXIT_INTERNAL_ASSERT,
                "Point-source nodal total does not equal the requested active rate.");
        }
    }
}
