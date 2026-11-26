#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <sstream>

#include "point_sources.hpp"
#include "mpi_simulation.hpp"
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

} // namespace

std::uint64_t mesh_fingerprint(const Variables &var)
{
    // FNV-1a over 64-bit words: the coordinates, the mesh size and the remesh count.
    std::uint64_t h = 1469598103934665603ull;
    const auto mix = [&h](std::uint64_t word) { h ^= word; h *= 1099511628211ull; };
    mix(std::uint64_t(var.nnode)); mix(std::uint64_t(var.nelem)); mix(std::uint64_t(var.nremesh));
    const double *x = var.coord->data();
    for (std::size_t i = 0; i < std::size_t(var.nnode) * NDIMS; ++i) {
        std::uint64_t word;
        std::memcpy(&word, x + i, sizeof word);
        mix(word);
    }
    return h | 1;
}

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

    // The owners belong to this mesh and these point positions.
    std::uint64_t mesh = mesh_fingerprint(var);
    for (const double_vec *v : {&injection.points_x, &injection.points_y, &injection.points_z})
        for (double x : *v) {
            std::uint64_t word;
            std::memcpy(&word, &x, sizeof word);
            mesh = (mesh ^ word) * 1099511628211ull;
        }
    mesh |= 1;
    PointSourceOwners &owners = var.point_owners[&injection];
    // The owner changes only with the mesh: find it once per mesh. All ranks
    // search again when any rank's mesh changed, for the check over the ranks.
    if (mpi_any(var, owners.mesh != mesh)) {
        owners.element.assign(injection.num_points, -1);
        owners.weight.assign(std::size_t(injection.num_points)*NODES_PER_ELEM, 0.0);
        int_vec placed(injection.num_points, 0);
        for (int p = 0; p < injection.num_points; ++p) {
            double q[NDIMS];
            point_position(injection, p, q);
            for (int d = 0; d < NDIMS; ++d) {
                if (!std::isfinite(q[d])) {
                    point_source_error(EXIT_RUNTIME_NAN, p, q, var,
                                       "has a non-finite coordinate");
                }
            }
            for (int e = 0; e < var.nelem; ++e) {
                double weights[NODES_PER_ELEM];
                if (!compute_point_shape_weights(var, e, q, weights)) continue;

                // Shared faces can have several candidates. Element order is
                // stable, so selecting the first maps the requested rate once and
                // makes ownership deterministic across thread counts.
                owners.element[p] = e;
                for (int i = 0; i < NODES_PER_ELEM; ++i)
                    owners.weight[p*NODES_PER_ELEM+i] = weights[i];
                placed[p] = 1;
                break;
            }
        }
#ifdef USEMPI
        // A compute rank holds part of the mesh: a point outside it can lie on
        // another rank. Rates that cancel hide a lost point from the totals
        // that assemble_point_sources checks, so each point is checked here.
        if (var.mpi && var.mpi->in_kernel())
            var.mpi->maximum(placed.data(), injection.num_points);
#endif
        for (int p = 0; p < injection.num_points; ++p) {
            if (placed[p]) continue;
            double q[NDIMS];
            point_position(injection, p, q);
            if (mpi_reports(var))
                point_source_error(EXIT_RUNTIME_LOOKUP, p, q, var,
                                   "could not be mapped into the current mesh");
            die_together(EXIT_RUNTIME_LOOKUP);
        }
        owners.mesh = mesh;
    }

    for (int p = 0; p < injection.num_points; ++p) {
        const int mapped_element = owners.element[p];
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

void assemble_point_sources(const Injection &injection, const Variables &var,
                           double start, double stop, double_vec &source,
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
    deposit_point_sources(injection,var,start,stop,add_point_source,&assembly,scale,peak_rate);
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
