#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

#include "constants.hpp"
#include "parameters.hpp"
#include "fields.hpp"
#include "geometry.hpp"
#include "hydraulics.hpp"
#include "matprops.hpp"
#include "utils.hpp"
#include "point_sources.hpp"
#include "mpi_simulation.hpp"

namespace {

// Calls apply for each flux-boundary facet with the element that owns it.
template <class Apply>
void flux_facets(const Variables &var, Apply apply)
{
    for (int side = 0; side < nbdrytypes_hydro; ++side) {
        if (var.hbc_types[side] != 3) continue;
        for (const auto &face : *var.bfacets[side]) {
            const int e = face.first, f = face.second;
            int nodes[NODES_PER_FACET];
            for (int j = 0; j < NODES_PER_FACET; ++j)
                nodes[j] = (*var.connectivity)[e][NODE_OF_FACET[f][j]];
            apply(side,e,f,compute_facet_measure(var.coord->view_const(nodes)));
        }
    }
}

// Dirichlet side, or -1 for a free pressure node.
#pragma acc routine seq
inline int pressure_boundary(const int n, const Variables &var)
{
    const int types[nbdrytypes_hydro] = {BOUNDX0, BOUNDX1, BOUNDY0, BOUNDY1, BOUNDZ0, BOUNDZ1};
    int fixed = -1;
    for (int i = 0; i < nbdrytypes_hydro; ++i) {
        if (!((*var.bcflag)[n] & types[i])) continue;
        if (var.hbc_types[i] == 2) return i;
        if (var.hbc_types[i] == 1) fixed = i;
    }
    return fixed;
}

// Pressure of a prescribed-pressure (type 2) side at schedule ratio; false for
// a free node (side -1) or a type-1 side.
#pragma acc routine seq
inline bool prescribed_pressure(const Param &param, const Variables &var, int side, double ratio,
                                double &value)
{
    if (side < 0 || var.hbc_types[side] != 2) return false;
    value = param.bc.hbc_pressure[side] * ratio;
    return true;
}

// Exact mean of the well rate ratio over [t0, t1] (in s). The ratio is linear
// between schedule times, a repeated time being a step, and held outside them.
double schedule_mean(const double_vec &t, const double_vec &r, double t0, double t1)
{
    const double a = t0 / YEAR2SEC, b = t1 / YEAR2SEC;
    const int n = static_cast<int>(t.size());
    if (a == b) {
        if (a < t.front()) return r.front();
        const int i = int(std::upper_bound(t.begin(), t.end(), a) - t.begin()) - 1;
        return i + 1 == n ? r.back() : r[i] + (r[i+1]-r[i])*(a-t[i])/(t[i+1]-t[i]);
    }
    // Integral of y0 + slope*(s - s0) over the part of [s0, s1] inside [a, b].
    auto part = [a, b](double s0, double s1, double y0, double slope) {
        const double lo = std::max(a, s0), hi = std::min(b, s1);
        return hi > lo ? (hi - lo) * (y0 + slope * (0.5 * (lo + hi) - s0)) : 0.0;
    };
    double sum = part(std::min(a, t[0]), t[0], r[0], 0.0) +
                 part(t[n-1], std::max(b, t[n-1]), r[n-1], 0.0);
    for (int i = 0; i + 1 < n; ++i)
        if (t[i+1] > t[i])
            sum += part(t[i], t[i+1], r[i], (r[i+1] - r[i]) / (t[i+1] - t[i]));
    return sum / (b - a);
}

// Well and flux-boundary ratios over [start, stop] and the prescribed-pressure
// ratio at stop.
void schedule_ratios(const BC &bc, double start, double stop,
                     double &well, double &boundary, double &pressure)
{
    const bool scheduled = !bc.hbc_schedule_ratio.empty();
    well = bc.num_wells > 0 ?
        schedule_mean(bc.well_schedule_time_in_yr, bc.well_schedule_ratio, start, stop) : 0.0;
    boundary = scheduled ?
        schedule_mean(bc.hbc_schedule_time_in_yr, bc.hbc_schedule_ratio, start, stop) : 1.0;
    pressure = scheduled ?
        schedule_mean(bc.hbc_schedule_time_in_yr, bc.hbc_schedule_ratio, stop, stop) : 1.0;
}

// Darcy mobility k/mu.
#pragma acc routine seq
inline double hydraulic_mobility(const MatProps &mat, int e)
{
    return mat.perm(e) / mat.mu_fluid(e);
}

// Fixed-stress lag L = alpha^2*compliance, in 1/Pa.
#pragma acc routine seq
inline double hydraulic_fixed_stress_lag(const MatProps &mat, int e)
{
    const double alpha = mat.alpha_biot(e);
    return alpha*alpha*hydraulic_mean_stress_compliance(mat,e);
}

// Mean p - p_old of the previous physical step, after the nodal residual slots.
const int pressure_d_slot = NODES_PER_ELEM;

// A plastic return breaks the elastic relation between mean stress and
// volume change. The source then reads the accepted projected strain and keeps
// the fixed-stress pressure buffer with its pending raw increment.
inline bool has_plastic_source(const Param &param)
{
    return param.control.has_poroelastic_pressure_feedback &&
           !hydraulic_strain_form(param) && (param.mat.rheol_type & MatProps::rh_plastic);
}

void pressure_increments(const Variables &var, const double_vec &ppressure, const double_vec &dppressure,
                         elem_cache &work)
{
    const int owned = mpi_element_end(var);
#ifndef ACC
    #pragma omp parallel for default(none) shared(var, ppressure, dppressure, work) firstprivate(owned)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector firstprivate(owned)
    for (int e = 0; e < owned; ++e) {
        const auto value = element_pressure<true>(var,e,ppressure,dppressure,true,false);
        work[e][pressure_d_slot] = value.raw;
    }
}

} // namespace

void allocate_hydraulic_state(Variables &var)
{
    var.hmass = new double_vec(var.nnode);
    var.fluid_source = new double_vec(var.nnode);
    var.fluid_balance = new double_vec(5, 0);
    var.old_mean_stress = new double_vec(var.nelem, 0);
    var.ppressure = new double_vec(var.nnode);
    var.dppressure = new double_vec(var.nnode);
}

void reset_hydraulic_assembly(const Param &param, Variables &var)
{
    delete var.hmass;
    var.hmass = new double_vec(var.nnode);
    delete var.fluid_source;
    var.fluid_source = new double_vec(var.nnode);
    if (param.control.has_hydraulic_diffusion) {
        delete var.old_mean_stress;
        var.old_mean_stress = new double_vec(var.nelem);
    }
    // Pressure and its pending increment were already transferred by the
    // nodal remapper. Reallocating them here would discard the remapped state.
}

void refresh_hydraulic_diffusivity(const Param &param, const Variables &var)
{
    if (!param.control.has_hydraulic_diffusion) return;
    double diffusivity = 0;
    int invalid = 0;
    const int owned = mpi_element_end(var);
#ifndef ACC
    #pragma omp parallel for default(none) shared(param, var) firstprivate(owned) reduction(max:diffusivity) reduction(|:invalid)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector firstprivate(owned) reduction(max:diffusivity) reduction(|:invalid)
    for (int e = 0; e < owned; ++e) {
        const double diff = hydraulic_diffusivity(param, *var.mat, e);
        if (!std::isfinite(diff) || diff < 0) invalid = 1;
        diffusivity = std::max(diffusivity, diff);
    }
    if (mpi_any(var, invalid))
        die_together(EXIT_RUNTIME_NAN, "Invalid current hydraulic properties.");
#ifdef USEMPI
    if (var.mpi && var.mpi->in_kernel()) diffusivity = var.mpi->maximum(diffusivity);
#endif
    var.mat->hydro_diff_max = diffusivity;
    var.mat->hydro_step_limit = 0;
}

namespace {

constexpr double screen_tolerance = 1e-10;

struct ElementBounds {
    double lower[NDIMS], upper[NDIMS];
};

// The pieces of the segment a + t d, t0 <= t <= t1, in the elements of the whole
// mesh, appended to pieces in order of t, each part counted once. Along a piece
// the barycentric coordinates are affine in t, N_i = l0 + dl t; scale converts a
// span of t to the length a piece records. Returns the length covered.
double clip_segment(const Variables &var, const std::vector<ElementBounds> &bounds,
                    const double a[NDIMS], const double d[NDIMS],
                    double t0, double t1, double scale, int well, WellPieces &pieces)
{
    struct Piece { int e; double t0, t1, l0[NODES_PER_ELEM], dl[NODES_PER_ELEM]; };
    double covered = 0;
    auto append = [&](const Piece &pc, double lo, double hi) {
        const double mid = 0.5 * (lo + hi);
        pieces.well.push_back(well);
        pieces.element.push_back(pc.e);
        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            pieces.start.push_back(pc.l0[i] + pc.dl[i] * lo);
            pieces.end.push_back(pc.l0[i] + pc.dl[i] * hi);
            pieces.mean.push_back(pc.l0[i] + pc.dl[i] * mid);
        }
        pieces.length.push_back((hi - lo) * scale);
        covered += hi - lo;
    };
    std::vector<Piece> found;
    double lower[NDIMS], upper[NDIMS];
    for (int i = 0; i < NDIMS; ++i) {
        lower[i] = std::min(a[i] + d[i] * t0, a[i] + d[i] * t1);
        upper[i] = std::max(a[i] + d[i] * t0, a[i] + d[i] * t1);
    }
    for (int e = 0; e < var.nelem; ++e) {
        const ElementBounds &box = bounds[e];
        bool apart = false;
        for (int i = 0; i < NDIMS; ++i)
            apart = apart || box.upper[i] < lower[i] || upper[i] < box.lower[i];
        if (apart) continue;
        ConstConnAccessor conn = (*var.connectivity)[e];
        ConstArrayAccessor v0 = (*var.coord)[conn[0]];
        double grad[NDIMS][NODES_PER_ELEM];
        get_shape_gradients(var, e, grad);
        Piece pc = {e, t0, t1, {}, {}};
        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            double l0 = (i == 0) + grad[0][i] * (a[0] - v0[0]) + grad[NDIMS-1][i] * (a[NDIMS-1] - v0[NDIMS-1]);
            double dl = grad[0][i] * d[0] + grad[NDIMS-1][i] * d[NDIMS-1];
#ifdef THREED
            l0 += grad[1][i] * (a[1] - v0[1]);
            dl += grad[1][i] * d[1];
#endif
            pc.l0[i] = l0;
            pc.dl[i] = dl;
            // Allow boundary roundoff for a contained segment. Otherwise
            // split at the face, where the two elements share the weights.
            if (l0 + dl*t0 >= -screen_tolerance && l0 + dl*t1 >= -screen_tolerance) continue;
            if (dl > 0) pc.t0 = std::max(pc.t0, -l0 / dl);
            else if (dl < 0) pc.t1 = std::min(pc.t1, -l0 / dl);
            else pc.t1 = pc.t0 - 1;
        }
        if (pc.t0 == t0 && pc.t1 == t1) {
            append(pc, t0, t1);
            return covered * scale;
        }
        if (pc.t1 > pc.t0) found.push_back(pc);
    }

    // A segment along a shared face or edge lies in several elements, whose
    // nodal weights agree there: count each part of it once.
    std::vector<double> cuts = {t0, t1};
    for (const Piece &pc : found) {
        cuts.push_back(pc.t0);
        cuts.push_back(pc.t1);
    }
    std::sort(cuts.begin(), cuts.end());
    const Piece *owner = nullptr;
    double start = t0;
    for (std::size_t k = 0; k + 1 < cuts.size(); ++k) {
        const double lo = cuts[k], hi = cuts[k+1], mid = 0.5 * (lo + hi);
        if (hi <= lo) continue;
        const auto it = std::find_if(found.begin(), found.end(),
                                     [mid](const Piece &pc) { return pc.t0 <= mid && mid <= pc.t1; });
        const Piece *next = it == found.end() ? nullptr : &*it;
        if (next == owner) continue;
        if (owner) append(*owner, start, lo);
        owner = next;
        start = lo;
    }
    if (owner) append(*owner, start, t1);
    return covered * scale;
}

// A single segment can represent adjacent pieces only if their join is at the
// same fraction of the straight line and of the configured screen length,
// within the placement tolerance in barycentric coordinates.
void merge_pieces(WellPieces &pieces)
{
    WellPieces merged;
    for (std::size_t k = 0; k < pieces.element.size(); ++k) {
        const std::size_t n = merged.element.size();
        const double *start = &pieces.start[k*NODES_PER_ELEM];
        const double *end = &pieces.end[k*NODES_PER_ELEM];
        const double *mean = &pieces.mean[k*NODES_PER_ELEM];
        const double length = pieces.length[k];
        if (n && merged.well[n-1] == pieces.well[k] && merged.element[n-1] == pieces.element[k]) {
            const double previous_length = merged.length[n-1];
            const double total = previous_length + length;
            const double fraction = previous_length / total;
            const std::size_t offset = (n-1)*NODES_PER_ELEM;
            bool straight = true;
            for (int i = 0; i < NODES_PER_ELEM; ++i) {
                const double first = merged.start[offset+i];
                const double join = first + fraction * (end[i] - first);
                straight &= std::fabs(start[i] - merged.end[offset+i]) <= screen_tolerance &&
                            std::fabs(start[i] - join) <= screen_tolerance;
            }
            if (straight) {
                for (int i = 0; i < NODES_PER_ELEM; ++i) {
                    double &m = merged.mean[offset+i];
                    m = (m * previous_length + mean[i] * length) / total;
                    merged.end[offset+i] = end[i];
                }
                merged.length[n-1] = total;
                continue;
            }
        }
        merged.well.push_back(pieces.well[k]);
        merged.element.push_back(pieces.element[k]);
        merged.start.insert(merged.start.end(), start, start+NODES_PER_ELEM);
        merged.end.insert(merged.end.end(), end, end+NODES_PER_ELEM);
        merged.mean.insert(merged.mean.end(), mean, mean+NODES_PER_ELEM);
        merged.length.push_back(length);
    }
    pieces = std::move(merged);
}

} // namespace

void well_screen_positions(const Variables &var, double_vec &ends)
{
    const WellPieces &pieces = var.well_pieces;
    ends.resize(pieces.element.size()*2*NDIMS);
    for (std::size_t k = 0; k < pieces.element.size(); ++k) {
        barycentric_point(var, pieces.element[k], &pieces.start[k*NODES_PER_ELEM], &ends[(2*k)*NDIMS]);
        barycentric_point(var, pieces.element[k], &pieces.end[k*NODES_PER_ELEM], &ends[(2*k+1)*NDIMS]);
    }
}

void place_well_screens(const Param &param, Variables &var, const double_vec *ends)
{
    const BC &bc = param.bc;
    if (!param.control.has_hydraulic_diffusion || bc.num_wells == 0) {
        var.well_pieces = {};
        return;
    }
    // Pad the element boxes to include points admitted by the barycentric
    // tolerance in clip_segment.
    std::vector<ElementBounds> bounds(var.nelem);
    for (int e = 0; e < var.nelem; ++e) {
        ConstConnAccessor conn = (*var.connectivity)[e];
        ElementBounds &box = bounds[e];
        double size = 0;
        for (int d = 0; d < NDIMS; ++d) {
            box.lower[d] = box.upper[d] = (*var.coord)[conn[0]][d];
            for (int j = 1; j < NODES_PER_ELEM; ++j) {
                const double x = (*var.coord)[conn[j]][d];
                box.lower[d] = std::min(box.lower[d], x);
                box.upper[d] = std::max(box.upper[d], x);
            }
            size = std::max(size, box.upper[d] - box.lower[d]);
        }
        for (int d = 0; d < NDIMS; ++d) {
            box.lower[d] -= 1e-6 * size;
            box.upper[d] += 1e-6 * size;
        }
    }
    WellPieces placed;
    const WellPieces &previous = var.well_pieces;
    for (int w = 0; w < bc.num_wells; ++w) {
        const double zb = bc.well_z_bottom[w], zt = bc.well_z_top[w];
        double covered = 0;
        if (!ends) {
            // The configured vertical screen, with z for t.
            double a[NDIMS] = {bc.well_x[w]}, d[NDIMS] = {0};
#ifdef THREED
            a[1] = bc.well_y[w];
#endif
            d[NDIMS-1] = 1;
            covered = clip_segment(var, bounds, a, d, zb, zt, 1, w, placed);
        } else {
            // Each piece of the previous mesh, straight between its two ends,
            // keeps its length in the material.
            for (std::size_t k = 0; k < previous.element.size(); ++k) {
                if (previous.well[k] != w) continue;
                const double *a = &(*ends)[(2*k)*NDIMS], *b = &(*ends)[(2*k+1)*NDIMS];
                double d[NDIMS];
                for (int i = 0; i < NDIMS; ++i) d[i] = b[i] - a[i];
                covered += clip_segment(var, bounds, a, d, 0, 1, previous.length[k], w, placed);
            }
        }
        if (covered < (1 - 1e-9) * (zt - zb)) {
            std::cerr << "Error: the screen of well " << w << " is not entirely inside the mesh.\n";
            die(EXIT_CONFIG_VALUE);
        }
    }
    if (ends) merge_pieces(placed);
    var.well_pieces = std::move(placed);
}

void assemble_well_source(const Param &param, const Variables &var)
{
    if (!param.control.has_hydraulic_diffusion || param.bc.num_wells == 0) return;
    #pragma acc wait
    const BC &bc = param.bc;
    const WellPieces &pieces = var.well_pieces;
    double_vec &source = *var.fluid_source;
    std::fill(source.begin(), source.end(), 0.0);
    // The screen's rate, spread uniformly over its configured length.
    for (std::size_t k = 0; k < pieces.element.size(); ++k) {
        const int w = pieces.well[k];
        const double q = bc.well_rate[w] * pieces.length[k] / (bc.well_z_top[w] - bc.well_z_bottom[w]);
        ConstConnAccessor conn = (*var.connectivity)[pieces.element[k]];
        for (int i = 0; i < NODES_PER_ELEM; ++i)
            source[conn[i]] += q * pieces.mean[k*NODES_PER_ELEM+i];
    }
}

void validate_hydraulic_boundaries(const Param &param, const Variables &var)
{
    if (!param.control.has_hydraulic_diffusion) return;
    // Imported boundary tags need not be box faces. Only actual shared nodes
    // constrain two prescribed values; type 2 takes precedence over type 1.
    const uint tags[nbdrytypes_hydro] = {BOUNDX0, BOUNDX1, BOUNDY0, BOUNDY1, BOUNDZ0, BOUNDZ1};
    for (int i = 0; i < nbdrytypes_hydro; ++i) {
        if (var.hbc_types[i] != 2) continue;
        for (int j = 0; j < i; ++j) {
            if (var.hbc_types[j] != 2 || param.bc.hbc_pressure[i] == param.bc.hbc_pressure[j]) continue;
            const uint both = tags[i] | tags[j];
            for (int n = 0; n < var.nnode; ++n) {
                if (((*var.bcflag)[n] & both) != both) continue;
                std::cerr << "Node " << n << ": bc.hbc_pressure[" << j << "] and [" << i << "] disagree.\n";
                die(EXIT_CONFIG_VALUE, "Intersecting prescribed-pressure tags must agree at their shared nodes.");
            }
        }
    }
}


void initial_hydrostatic_state(const Param &param, const Variables &var,
                               double_vec &ppressure, double_vec &dppressure)
{
    // Prescribed pressures at the start, scaled by their schedule as the pressure
    // step scales them; the linear profile runs between the scaled end values.
    const double ratio = param.bc.hbc_schedule_ratio.empty() ? 1 :
        schedule_mean(param.bc.hbc_schedule_time_in_yr, param.bc.hbc_schedule_ratio,
                      var.time, var.time);
    if (param.ic.pore_pressure_option != 0) {
        const double a = param.bc.hbc_pressure[0] * ratio, b = param.bc.hbc_pressure[1] * ratio;
        for (int n = 0; n < var.nnode; ++n) {
            double p = param.ic.pore_pressure_option == 1 ?
                a + (b-a)*(*var.coord)[n][0]/param.mesh.xlength : param.ic.pore_pressure_value;
            prescribed_pressure(param, var, pressure_boundary(n, var), ratio, p);
            ppressure[n] = p; dppressure[n] = 0;
        }
        return;
    }
    for (int n = 0; n < var.nnode; ++n)
        if (prescribed_pressure(param, var, pressure_boundary(n, var), ratio, ppressure[n]))
            dppressure[n] = 0;
    if (param.control.gravity == 0) return;

    // Hydrostatic condition for pressure calculation
    const double rho_fluid = param.mat.fluid_rho0[param.mat.mattype_ref];
    double gravity = param.control.gravity;
    double loading = param.ic.excess_pore_pressure;
    const double bulk = var.mat->bulkm(0);
    const double shear = var.mat->shearm(0);
    const double alpha = var.mat->alpha_biot(0);
    const double loading_efficiency = alpha / ((bulk + 4.0 * shear / 3.0) * hydraulic_inverse_biot_modulus(*var.mat, 0) + alpha * alpha);

    for (int i = 0; i < var.nnode; ++i) {
        ConstArrayAccessor coord = (*var.coord)[i];
        double z = coord[NDIMS - 1];

        const int side = pressure_boundary(i, var);
        if (side >= 0 && var.hbc_types[side] == 2) continue;
        // Pressure is compression-positive; z increases upward.
        ppressure[i] = -rho_fluid * gravity * z;

        // Add excess pore pressure for non-boundary nodes
        if (side < 0) ppressure[i] += loading_efficiency * loading;
        dppressure[i] = 0.0;
    }
}

void update_old_mean_stress(const Param& param, const Variables& var, tensor_t& stress,
                   double_vec& old_mean_stress)
{
#ifdef USEMPI
    if (var.mpi && var.mpi->dispatch(MpiOp::OldMeanStress, var, 0.0)) return;
#endif
    const int mpi_ee = mpi_element_end(var);
#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif

#ifndef ACC
    #pragma omp parallel for default(none)                           \
        shared(param, var, stress, old_mean_stress) firstprivate(mpi_ee)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector async firstprivate(mpi_ee)
    for (int e=0; e<mpi_ee; ++e) {
        ConstTensorAccessor s = stress[e];
        old_mean_stress[e] =trace(s)/NDIMS;
    }
#ifdef NPROF
    nvtxRangePop();
#endif
}


static void pressure_step(const Param &param, const Variables &var,
                          double_vec &ppressure, double_vec &dppressure, double_vec &tdot, elem_cache &tmp_result, tensor_t& stress, double_vec& old_mean_stress, double dt, double time)
{
#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif

    double diff_max_local = 0, source = 0;
    const bool plastic_source = has_plastic_source(param);
    const bool strain_source = hydraulic_strain_form(param) || plastic_source;
    const int mpi_ee = mpi_element_end(var), mpi_ne = mpi_node_end(var);

#ifndef ACC
    #pragma omp parallel for default(none) shared(var, ppressure, tmp_result, stress, old_mean_stress, param) \
        firstprivate(dt, plastic_source, strain_source, mpi_ee) reduction(max:diff_max_local) reduction(+:source)
#endif
        #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector firstprivate(dt, plastic_source, strain_source, mpi_ee) reduction(max:diff_max_local) reduction(+:source) async
    for (int e = 0; e < mpi_ee; e++) {
        ElemCacheAccessor tr = tmp_result[e];
        double current_mean_stress = trace(stress[e])/NDIMS;
        double mean_stress_change = current_mean_stress - old_mean_stress[e];

        double alpha_b = var.mat->alpha_biot(e);         // Biot coefficient
        double bulk_comp = hydraulic_mean_stress_compliance(*var.mat, e);


        // Darcy mobility k/mu in pressure form, well-defined when gravity is zero;
        // it is algebraically equivalent to the head form
        //   (k * gamma_w / mu) * grad(p / gamma_w + z).
        const double kv_pressure = hydraulic_mobility(*var.mat, e) * (*var.volume)[e];

        const double diff_e = hydraulic_diffusivity(param, *var.mat, e);
        diff_max_local = std::max(diff_max_local, diff_e);

        // volume term (poroelastic effect)
        double pe = 0;
        if (param.control.has_poroelastic_pressure_feedback && strain_source)
            // Match the volumetric rate used by the anti-locking constitutive update.
            pe = alpha_b * (*var.edvoldt)[e] * (*var.volume)[e] / NODES_PER_ELEM;
        else if (param.control.has_poroelastic_pressure_feedback)
            pe = alpha_b * mean_stress_change * bulk_comp * (*var.volume)[e] / NODES_PER_ELEM / var.dt;
        if (plastic_source) {
            // The strain source already includes the plastic volume change;
            // only the previous raw pressure increment belongs to the buffer.
            pe -= hydraulic_fixed_stress_lag(*var.mat, e) * tr[pressure_d_slot] * (*var.volume)[e] /
                  (NODES_PER_ELEM * var.dt);
        }
        source += pe * NODES_PER_ELEM * dt;
        double grad[NDIMS][NODES_PER_ELEM];
        get_shape_gradients(var, e, grad);
        // One potential gradient and flux per element, shared by its nodes.
        double drive[NDIMS], flux[NDIMS];
        hydraulic_drive(param, var, ppressure, e, grad, drive);
        for (int d = 0; d < NDIMS; ++d) flux[d] = kv_pressure * drive[d];
        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            tr[i] = pe;
            for (int d = 0; d < NDIMS; ++d) tr[i] += grad[d][i] * flux[d];
        }
    }

    // Ghost elements take their rows from the ranks that own them; the boundary
    // flux below is added afterwards to every local facet, so each owned node
    // receives each facet once.
    mpi_sync(var, tmp_result, MpiEntity::Element, NODES_PER_ELEM);

    // Wells extract their exact schedule volume over this step [time - dt, time].
    const bool has_wells = param.bc.num_wells > 0;
    const bool has_points = param.injection.enabled;
    if (has_points) {
        assemble_point_sources(param.injection,var.fluid_points,var,time-dt,time,*var.fluid_point_source);
        for (int n=0;n<mpi_ne;++n)
            if ((*var.fluid_point_source)[n] != 0 &&
                (pressure_boundary(n,var) >= 0 || !((*var.hmass)[n]>0 && std::isfinite((*var.hmass)[n]))))
                die(EXIT_RUNTIME_LOOKUP, "A fluid point source reaches a fixed-pressure node or nonpositive hydraulic storage.");
    }
    double well_ratio, boundary_ratio, pressure_ratio;
    schedule_ratios(param.bc, time - dt, time, well_ratio, boundary_ratio, pressure_ratio);
    #pragma acc wait
    double prescribed_flux = 0;
    flux_facets(var,[&](int side, int e, int f, double area) {
        const double q = boundary_ratio * param.bc.hbc_flux[side] * area;
        if (e < mpi_ee) prescribed_flux += q * dt;
        for (int k = 0; k < NODES_PER_FACET; ++k)
            tmp_result[e][NODE_OF_FACET[f][k]] += q / NODES_PER_FACET;
    });
    double stored = 0, wells = 0, boundary = -prescribed_flux;

#ifndef ACC
    #pragma omp parallel for default(none) shared(param, var, tdot, ppressure, dppressure, tmp_result) \
        firstprivate(has_wells, has_points, well_ratio, pressure_ratio, dt, mpi_ne) \
        reduction(+:stored, wells, boundary)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector reduction(+:stored, wells, boundary) firstprivate(has_points, mpi_ne)
    for (int n = 0; n < mpi_ne; n++) {
        tdot[n] = 0.0;
        {
            const int npatch = var.support.size(n);
            const int* patch = var.support.patch(n);
            const int* lpatch = var.support.local(n);
            for (int k = 0; k < npatch; ++k)
                tdot[n] += tmp_result[patch[k]][lpatch[k]];
        }
        if (has_wells)
            tdot[n] += well_ratio * (*var.fluid_source)[n];
        if (has_points) tdot[n] -= (*var.fluid_point_source)[n];

        const double p0 = ppressure[n], mass = (*var.hmass)[n];
        const int side = pressure_boundary(n, var);
        double increment = 0;
        if (side < 0 && mass > 0)
            increment = tdot[n] * dt / mass;
        else if (side >= 0 && var.hbc_types[side] == 2)
            increment = p0 - param.bc.hbc_pressure[side] * pressure_ratio;
        const double change = -mass * increment;
        ppressure[n] = p0 - increment;
        dppressure[n] += increment;
        stored += change;
        if (has_wells) wells += well_ratio * (*var.fluid_source)[n] * dt;
        if (has_points) wells -= (*var.fluid_point_source)[n]*dt;
        if (side >= 0) boundary += change + tdot[n]*dt;

    }

    #pragma acc wait
#ifdef USEMPI
    if (var.mpi && var.mpi->in_kernel()) {
        double sums[] = {stored, source, wells, boundary};
        var.mpi->sum(sums, 4);
        stored = sums[0]; source = sums[1]; wells = sums[2]; boundary = sums[3];
        diff_max_local = var.mpi->maximum(diff_max_local);
    }
#endif
    if (var.fluid_balance) {
        const double value[] = {stored, source, wells, boundary, stored + source + wells - boundary};
        for (int i = 0; i < 5; ++i) (*var.fluid_balance)[i] += value[i];
    }
    var.mat->hydro_diff_max = diff_max_local;
#ifdef NPROF
    nvtxRangePop();
#endif
}


static double pressure_step_limit(const Variables &var, elem_cache &work)
{
    // Gershgorin bound on the lumped operator.
    // Reuse assembly scratch; no pressure or material history is copied.
    // These kernels run off the async queue: let it finish with work and the
    // masses first.
    #pragma acc wait
    const int mpi_ee = mpi_element_end(var);
#ifndef ACC
    #pragma omp parallel for default(none) shared(var, work) firstprivate(mpi_ee)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector firstprivate(mpi_ee)
    for (int e = 0; e < mpi_ee; ++e) {
        double g[NDIMS][NODES_PER_ELEM];
        get_shape_gradients(var, e, g);
        const double kv = hydraulic_mobility(*var.mat, e) * (*var.volume)[e];
        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            double row = 0;
            for (int j = 0; j < NODES_PER_ELEM; ++j) {
                double a = 0;
                for (int d = 0; d < NDIMS; ++d) a += g[d][i] * g[d][j];
                row += std::abs(a);
            }
            work[e][i] = kv * row;
        }
    }
    mpi_sync(var, work, MpiEntity::Element, NODES_PER_ELEM);
    return nodal_step_limit(var, work, *var.hmass, std::numeric_limits<double>::max());
}

void update_pore_pressure(const Param &param, const Variables &var,
                         double_vec &p, double_vec &dp, double_vec &rate,
                         elem_cache &work, tensor_t &stress, double_vec &old_stress)
{
#ifdef USEMPI
    if (var.mpi && var.mpi->dispatch(MpiOp::PorePressure, var, 0.0)) return;
#endif
    if (var.mat->hydro_step_limit == 0)
        var.mat->hydro_step_limit = pressure_step_limit(var, work);
    // Before the pending increment is cleared below.
    if (has_plastic_source(param))
        pressure_increments(var, p, dp, work);
    const double limit = var.mat->hydro_step_limit;
    #pragma acc wait
    std::fill(dp.begin(), dp.end(), 0.0);
    const int steps = explicit_substeps(var.dt, limit,
        "Hydraulic timestep needs too many explicit substeps; reduce fixed_dt.");
    const double dt = var.dt / steps;
    if (var.fluid_balance) std::fill(var.fluid_balance->begin(), var.fluid_balance->end(), 0.0);
    for (int i = 0; i < steps; ++i) {
        // Each substep reads its neighbours' ghost pressures; the dispatcher exchanges the last.
        if (i > 0) mpi_sync(var, p, MpiEntity::Node);
        pressure_step(param, var, p, dp, rate, work, stress, old_stress,
                      dt, i+1 == steps ? var.time : var.time - var.dt + (i+1)*dt);
    }
}


namespace {

// Pressure storage and fixed-stress buffer coefficient of element e in the
// remapped content. With mechanical feedback off, storage is the reduced flow
// model's capacity.
void remap_coefficients(const Param &param, const Variables &var, int e, double &storage, double &lag)
{
    const bool feedback = param.control.has_poroelastic_pressure_feedback;
    storage = feedback ? hydraulic_inverse_biot_modulus(*var.mat,e) : hydraulic_storage(param,*var.mat,e);
    lag = feedback && !hydraulic_strain_form(param) ? hydraulic_fixed_stress_lag(*var.mat,e) : 0;
}

// Mobile physical content and numerical fixed-stress buffer.
void remap_content(const Param &param, const Variables &var, int e, int j, double value[2])
{
    const double alpha = var.mat->alpha_biot(e);
    const bool feedback = param.control.has_poroelastic_pressure_feedback;
    double storage, lag;
    remap_coefficients(param,var,e,storage,lag);
    const int n = (*var.connectivity)[e][j];
    value[0] = storage*(*var.ppressure)[n];
    if (feedback) value[0] += alpha*trace((*var.strain)[e]);
    value[1] = -lag*(*var.dppressure)[n];
}

} // namespace

ContentRemap::ContentRemap(const Param &param, const Variables &var) :
    linear(param.control.has_hydraulic_diffusion &&
           (!param.control.has_poroelastic_pressure_feedback ||
            (param.mat.rheol_type & MatProps::rh_elastic)))
{
    if (!linear) return;
    #pragma acc wait
    old_content.resize(var.nelem);
    for (int e = 0; e < var.nelem; ++e)
        for (int j = 0; j < NODES_PER_ELEM; ++j) {
            double value[2];
            remap_content(param,var,e,j,value);
            for (int k = 0; k < 2; ++k) {
                old_content[e][k*NODES_PER_ELEM+j] = value[k];
                total[k] += (*var.volume)[e]*value[k]/NODES_PER_ELEM;
                magnitude[k] += (*var.volume)[e]*std::fabs(value[k])/NODES_PER_ELEM;
            }
        }
}

void ContentRemap::transfer(const Variables &var, const array_t &old_coordinate,
                         const conn_t &old_connectivity)
{
    if (!linear) return;
    #pragma acc wait
    // Centroid of a simplex and the radius about it that encloses the simplex.
    auto enclose = [](const auto &coordinate, double center[NDIMS]) {
        double radius = 0;
        for (int d = 0; d < NDIMS; ++d) {
            center[d] = 0;
            for (int j = 0; j < NODES_PER_ELEM; ++j) center[d] += coordinate[j][d]/NODES_PER_ELEM;
        }
        for (int j = 0; j < NODES_PER_ELEM; ++j) {
            double distance2 = 0;
            for (int d = 0; d < NDIMS; ++d) distance2 += (coordinate[j][d]-center[d])*(coordinate[j][d]-center[d]);
            radius = std::max(radius,std::sqrt(distance2));
        }
        return radius;
    };
    const int count = old_connectivity.size();
    array_t centers(count);
    double_vec radii(count, 0.0);
    double radius = 0;
    for (int e = 0; e < count; ++e) {
        double center[NDIMS];
        radii[e] = enclose(old_coordinate.view_const(old_connectivity[e]),center);
        for (int d = 0; d < NDIMS; ++d) centers[e][d] = center[d];
        radius = std::max(radius,radii[e]);
    }
    PointCloud cloud(centers);
    NANOKDTree tree(NDIMS,cloud);
    std::vector<nanoflann::ResultItem<std::uint32_t,double>> candidates;
    candidates.reserve(64);
    outside.resize(var.nelem,0.0);
    for (int e = 0; e < var.nelem; ++e) {
        const auto coordinate = var.coord->view_const((*var.connectivity)[e]);
        const double volume = compute_volume(coordinate);
        double center[NDIMS];
        const double local_radius = enclose(coordinate,center);
        for (int j = 0; j < NODES_PER_ELEM; ++j) outside[e][j] = volume/NODES_PER_ELEM;
        const double search_radius = (radius+local_radius)*(1+1e-12);
        tree.radiusSearch(center,search_radius*search_radius,candidates);
        for (const auto &candidate : candidates) {
            const int old = candidate.first;
            // The search radius covers the largest donor. Simplices whose
            // enclosing balls about their centroids are apart cannot overlap.
            const double reach = (radii[old]+local_radius)*(1+1e-12);
            if (candidate.second > reach*reach) continue;
            const SimplexOverlap overlap = simplex_overlap(old_coordinate.view_const(old_connectivity[old]),coordinate);
            if (overlap.volume == 0) continue;
            for (int i = 0; i < NODES_PER_ELEM; ++i)
                for (int j = 0; j < NODES_PER_ELEM; ++j) {
                    outside[e][i] -= overlap.product[i][j];
                    for (int k = 0; k < 2; ++k)
                        mapped[k] += overlap.product[i][j]*old_content[old][k*NODES_PER_ELEM+j];
                }
        }
        for (int j = 0; j < NODES_PER_ELEM; ++j) {
            if (outside[e][j] < -1e-9*volume) {
                std::fprintf(stderr,"Error: remap element %d, node %d has excess donor weight %.17g of volume %.17g.\n",
                             e,j,-outside[e][j],volume);
                for (int i = 0; i < NODES_PER_ELEM; ++i) {
                    for (int d = 0; d < NDIMS; ++d) std::fprintf(stderr," %.17g",coordinate[i][d]);
                    std::fputc('\n',stderr);
                }
                die(EXIT_MESH_QUALITY,"Overlapping donor cells in fluid remapping.");
            }
            if (std::fabs(outside[e][j]) < 1e-12*volume) outside[e][j] = 0;
        }
    }
}

void ContentRemap::restore(const Param &param, Variables &var) const
{
    if (!linear) return;
    #pragma acc wait
    double current[2] = {}, current_magnitude[2] = {}, imported[2] = {}, capacity[2] = {};
    for (int e = 0; e < var.nelem; ++e) {
        const double volume = (*var.volume)[e];
        double storage, lag;
        remap_coefficients(param,var,e,storage,lag);
        for (int j = 0; j < NODES_PER_ELEM; ++j) {
            double value[2]; remap_content(param,var,e,j,value);
            for (int k = 0; k < 2; ++k) {
                current[k] += volume*value[k]/NODES_PER_ELEM;
                current_magnitude[k] += volume*std::fabs(value[k])/NODES_PER_ELEM;
                imported[k] += outside[e][j]*value[k];
            }
            if (pressure_boundary((*var.connectivity)[e][j],var) < 0) {
                capacity[0] += volume*storage/NODES_PER_ELEM;
                capacity[1] += volume*lag/NODES_PER_ELEM;
            }
        }
    }
    double correction[2] = {}, reaction[2] = {};
    for (int k = 0; k < 2; ++k) {
        const double difference = mapped[k]+imported[k]-current[k];
        const double tolerance = 128*std::numeric_limits<double>::epsilon()*
                                 std::max(magnitude[k]+current_magnitude[k],1e-30);
        if (std::fabs(difference) <= tolerance) continue;
        if (capacity[k] > 0) correction[k] = difference/capacity[k];
        else reaction[k] = -difference;
        if (!std::isfinite(correction[k])) die(EXIT_RUNTIME_NAN,"Invalid fluid remap correction.");
    }
    // The buffer is -L*dppressure. Reconstruct it independently of real water.
    correction[1] = -correction[1];
    for (int n = 0; n < var.nnode; ++n) {
        if (pressure_boundary(n,var) >= 0) continue;
        (*var.ppressure)[n] += correction[0];
        (*var.dppressure)[n] += correction[1];
    }
    for (int e = 0; e < var.nelem; ++e) {
        int free = 0;
        for (int j = 0; j < NODES_PER_ELEM; ++j)
            free += pressure_boundary((*var.connectivity)[e][j],var) < 0;
        const double mechanical = correction[0]+(hydraulic_strain_form(param) ? 0 : correction[1]);
        const double shift = -var.mat->alpha_biot(e)*mechanical*free/NODES_PER_ELEM;
        for (int d = 0; d < NDIMS; ++d) (*var.stress)[e][d] += shift;
        if (var.mat->is_plane_strain) (*var.stressyy)[e] += shift;
    }
    double remaining[2] = {};
    for (int e = 0; e < var.nelem; ++e)
        for (int j = 0; j < NODES_PER_ELEM; ++j) {
            double value[2]; remap_content(param,var,e,j,value);
            for (int k = 0; k < 2; ++k) remaining[k] += (*var.volume)[e]*value[k]/NODES_PER_ELEM;
        }
    for (int k = 0; k < 2; ++k)
        if (std::fabs(remaining[k]-mapped[k]-imported[k]-reaction[k]) >
            1e-10*std::max(magnitude[k]+current_magnitude[k],1e-30))
            die(EXIT_INTERNAL_ASSERT,"Fluid remapping fails its inventory constraint.");
    std::printf("    Biot content remap: removed/added/Dirichlet [m^%d] = %.17g %.17g %.17g; "
                "dp [Pa] = %.17g; buffer correction [m^%d] = %.17g\n",
                NDIMS,total[0]-mapped[0],imported[0],reaction[0],correction[0],NDIMS,
                capacity[1]*(-correction[1]));
}
