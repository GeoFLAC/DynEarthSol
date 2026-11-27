#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
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

// Calls apply for each flux-boundary facet with the element that owns it:
// density at a shared boundary node can differ between adjoining materials;
// only the facet geometry is common.
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
            apply(side,e,f,nodes,compute_facet_measure(var.coord->view_const(nodes)));
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

// Darcy mobility k/mu, before relative permeability.
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

// Fluid density at node n's pressure and temperature with element e's law.
#pragma acc routine seq
inline double nodal_fluid_density(const Variables &var, int e, int n)
{
    return fluid_density(*var.mat,e,(*var.ppressure)[n],(*var.temperature)[n],var.fluid_reference_temperature);
}

// The Biot model's porosity range: 0 < phi < 1 and phi <= alpha.
#pragma acc routine seq
inline bool biot_porosity_valid(const MatProps &mat, int e, double phi)
{
    return phi > 0 && phi < 1 && phi <= mat.alpha_biot(e);
}

// Pore dilation after its recovery and the accepted plastic increment.
#pragma acc routine seq
inline void update_pore_dilation(const Variables &var, int e, double recovery)
{
    (*var.pore_dilation)[e] = std::max(0.0,(*var.pore_dilation)[e]-recovery+
                                       (*var.pore_plastic_increment)[e]);
}

// Entry (i, j) of g^T K g, the tensor part of an element's Darcy stiffness.
#pragma acc routine seq
inline double permeability_coupling(const MatProps &mat, const double g[NDIMS][NODES_PER_ELEM],
                                    int i, int j)
{
    double a = 0;
    for (int d = 0; d < NDIMS; ++d)
        for (int k = 0; k < NDIMS; ++k)
            a += g[d][i] * mat.perm_tensor[d][k] * g[k][j];
    return a;
}


// Element scratch slot (after the nodal residuals) holding V*phi/(N*markers):
// times a material's marker count it is that material's pore weight at a node.
const int pore_slot = NODES_PER_ELEM;

// Bishop fixed-stress slots: D = mean(p - p_old), E = mean(Se*p - Se_old*p_old)
// of the last physical step. Retain E=D for a wet step through later suction.
const int pressure_d_slot = NODES_PER_ELEM + 1;
const int bishop_e_slot = NODES_PER_ELEM + 2;

// Retention laws in the fixed-stress branch with pressure feedback.
inline bool has_bishop_correction(const Param &param, const Variables &var)
{
    return var.mat->has_retention && param.control.has_poroelastic_pressure_feedback &&
           !hydraulic_strain_form(param);
}

// A plastic return breaks the elastic relation between mean stress and
// volume change. The source then reads the accepted projected strain and keeps
// the fixed-stress pressure buffer with its pending raw increment.
inline bool has_plastic_source(const Param &param)
{
    return param.control.has_poroelastic_pressure_feedback &&
           !hydraulic_strain_form(param) && (param.mat.rheol_type & MatProps::rh_plastic);
}

// Save D and E before clearing the pending increment. The source adds
// L*(chi*E - D)/dt, L = alpha^2*compliance: a lag that vanishes as dt shrinks.
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
        work[e][bishop_e_slot] = value.bishop;
    }
}

// The (element, material) entries of one node for advance_content.
struct PorePatch
{
    const MatProps &mat;
    const int_vec2D &markers;
    const elem_cache &work;
    const int *patch;
    int count;
    #pragma acc routine seq
    int size() const { return count * mat.nmat; }
    #pragma acc routine seq
    const Retention &law(int k) const { return mat.retention(k % mat.nmat); }
    #pragma acc routine seq
    double weight(int k) const
    {
        const int e = patch[k / mat.nmat];
        return work[e][pore_slot] * markers[e][k % mat.nmat];
    }
};

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

namespace {
#pragma acc routine seq
inline void hydraulic_flux(const Param &param, const Variables &var, int e,
                           const double g[NDIMS][NODES_PER_ELEM], double mobility, double density, double *q)
{
    for (int d = 0; d < NDIMS; ++d) q[d] = 0;
    const auto conn = (*var.connectivity)[e];
    double drive[NDIMS];
    if (param.control.hydraulic_content_model == 1) {
        for (int d = 0; d < NDIMS; ++d) {
            drive[d] = d == NDIMS-1 ? density*param.control.gravity : 0;
            for (int j = 0; j < NODES_PER_ELEM; ++j) drive[d] += g[d][j]*(*var.ppressure)[conn[j]];
        }
    } else hydraulic_drive(param, var, *var.ppressure, e, g, drive);
    if (var.mat->has_retention)
        mobility *= element_mean_kr(*var.mat, (*var.elemmarkers)[e], conn, *var.ppressure);
    for (int d = 0; d < NDIMS; ++d)
        for (int k = 0; k < NDIMS; ++k)
            q[d] -= mobility*var.mat->perm_tensor[d][k]*drive[k];
}

} // namespace

#pragma acc routine seq
void hydraulic_element_flux(const Param &param, const Variables &var, int e, double *q,
                            const double (*gradient)[NODES_PER_ELEM])
{
    if (!param.control.has_hydraulic_diffusion) {
        for (int d = 0; d < NDIMS; ++d) q[d] = 0;
        return;
    }
    double g[NDIMS][NODES_PER_ELEM];
    if (!gradient) {
        get_shape_gradients(var, e, g);
        gradient = g;
    }
    double density = 0;
    if (param.control.hydraulic_content_model == 1) {
        const auto conn = (*var.connectivity)[e];
        for (int j = 0; j < NODES_PER_ELEM; ++j)
            density += nodal_fluid_density(var,e,conn[j])/NODES_PER_ELEM;
    }
    hydraulic_flux(param,var,e,gradient,hydraulic_mobility(*var.mat,e),density,q);
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
        if (var.pore_strain && !biot_porosity_valid(*var.mat, e, var.mat->phi(e))) invalid = 1;
        diffusivity = std::max(diffusivity, diff);
    }
    if (mpi_any(var, invalid))
        die_together(EXIT_RUNTIME_NAN, "Invalid current hydraulic properties or porosity state.");
#ifdef USEMPI
    if (var.mpi && var.mpi->in_kernel()) diffusivity = var.mpi->maximum(diffusivity);
#endif
    var.mat->hydro_diff_max = diffusivity;
    var.mat->hydro_step_limit = 0;
}

namespace {

// The pieces of the segment a + t d, t0 <= t <= t1, in the elements of the whole
// mesh, appended to pieces in order of t, each part counted once. Along a piece
// the barycentric coordinates are affine in t, N_i = l0 + dl t; scale converts a
// span of t to the length a piece records. Returns the length covered.
double clip_segment(const Variables &var, const double a[NDIMS], const double d[NDIMS],
                    double t0, double t1, double scale, int well, WellPieces &pieces)
{
    // Barycentric tolerance, so that screens on element faces or edges are kept.
    const double tol = 1e-10;
    struct Piece { int e; double t0, t1, l0[NODES_PER_ELEM], dl[NODES_PER_ELEM]; };
    std::vector<Piece> found;
    // Elements away from the segment's bounding box cannot hold any of it; the
    // margin is far wider than the barycentric tolerance below.
    double lower[NDIMS], upper[NDIMS];
    for (int i = 0; i < NDIMS; ++i) {
        lower[i] = std::min(a[i] + d[i] * t0, a[i] + d[i] * t1);
        upper[i] = std::max(a[i] + d[i] * t0, a[i] + d[i] * t1);
    }
    for (int e = 0; e < var.nelem; ++e) {
        ConstConnAccessor conn = (*var.connectivity)[e];
        double low[NDIMS], high[NDIMS], size = 0;
        for (int i = 0; i < NDIMS; ++i) low[i] = high[i] = (*var.coord)[conn[0]][i];
        for (int k = 1; k < NODES_PER_ELEM; ++k)
            for (int i = 0; i < NDIMS; ++i) {
                low[i] = std::min(low[i], (*var.coord)[conn[k]][i]);
                high[i] = std::max(high[i], (*var.coord)[conn[k]][i]);
            }
        for (int i = 0; i < NDIMS; ++i) size = std::max(size, high[i] - low[i]);
        bool apart = false;
        for (int i = 0; i < NDIMS; ++i)
            apart = apart || high[i] + 1e-6 * size < lower[i] || upper[i] < low[i] - 1e-6 * size;
        if (apart) continue;
        ConstArrayAccessor v0 = (*var.coord)[conn[0]];
        double shpdx[NODES_PER_ELEM], shpdz[NODES_PER_ELEM];
#ifdef THREED
        double shpdy[NODES_PER_ELEM];
        get_local_shape_fn(var, e, shpdx, shpdy, shpdz);
#else
        get_local_shape_fn(var, e, shpdx, shpdz);
#endif
        Piece pc = {e, t0, t1, {}, {}};
        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            double l0 = (i == 0) + shpdx[i] * (a[0] - v0[0]) + shpdz[i] * (a[NDIMS-1] - v0[NDIMS-1]);
            double dl = shpdx[i] * d[0] + shpdz[i] * d[NDIMS-1];
#ifdef THREED
            l0 += shpdy[i] * (a[1] - v0[1]);
            dl += shpdy[i] * d[1];
#endif
            pc.l0[i] = l0;
            pc.dl[i] = dl;
            // Keep t with l0 + dl*t >= -tol.
            if (dl > 0) pc.t0 = std::max(pc.t0, (-tol - l0) / dl);
            else if (dl < 0) pc.t1 = std::min(pc.t1, (-tol - l0) / dl);
            else if (l0 < -tol) pc.t1 = pc.t0 - 1;
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
    double covered = 0;
    for (std::size_t k = 0; k + 1 < cuts.size(); ++k) {
        const double lo = cuts[k], hi = cuts[k+1], mid = 0.5 * (lo + hi);
        if (hi <= lo) continue;
        for (const Piece &pc : found) {
            if (pc.t0 > mid || mid > pc.t1) continue;
            pieces.well.push_back(well);
            pieces.element.push_back(pc.e);
            for (int i = 0; i < NODES_PER_ELEM; ++i) {
                pieces.start.push_back(pc.l0[i] + pc.dl[i] * lo);
                pieces.end.push_back(pc.l0[i] + pc.dl[i] * hi);
                pieces.mean.push_back(pc.l0[i] + pc.dl[i] * mid);
            }
            pieces.length.push_back((hi - lo) * scale);
            covered += hi - lo;
            break;
        }
    }
    return covered * scale;
}

// Consecutive pieces in one element become one, so that the count stays that of
// the elements crossed: the chord between the outer ends stands for the material
// line inside the element. Lengths add, and the mean coordinates, weighted by
// length, keep each node's share.
void merge_pieces(WellPieces &pieces)
{
    WellPieces merged;
    for (std::size_t k = 0; k < pieces.element.size(); ++k) {
        const std::size_t n = merged.element.size();
        const double *mean = &pieces.mean[k*NODES_PER_ELEM], length = pieces.length[k];
        if (n && merged.well[n-1] == pieces.well[k] && merged.element[n-1] == pieces.element[k]) {
            const double total = merged.length[n-1] + length;
            for (int i = 0; i < NODES_PER_ELEM; ++i) {
                double &m = merged.mean[(n-1)*NODES_PER_ELEM+i];
                m = (m * merged.length[n-1] + mean[i] * length) / total;
                merged.end[(n-1)*NODES_PER_ELEM+i] = pieces.end[k*NODES_PER_ELEM+i];
            }
            merged.length[n-1] = total;
            continue;
        }
        merged.well.push_back(pieces.well[k]);
        merged.element.push_back(pieces.element[k]);
        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            merged.start.push_back(pieces.start[k*NODES_PER_ELEM+i]);
            merged.end.push_back(pieces.end[k*NODES_PER_ELEM+i]);
            merged.mean.push_back(mean[i]);
        }
        merged.length.push_back(length);
    }
    pieces = merged;
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
    WellPieces placed;
    if (!param.control.has_hydraulic_diffusion || bc.num_wells == 0) {
        var.well_pieces = placed;
        return;
    }
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
            covered = clip_segment(var, a, d, zb, zt, 1, w, placed);
        } else {
            // Each piece of the previous mesh, straight between its two ends,
            // keeps its length in the material.
            for (std::size_t k = 0; k < previous.element.size(); ++k) {
                if (previous.well[k] != w) continue;
                const double *a = &(*ends)[(2*k)*NDIMS], *b = &(*ends)[(2*k+1)*NDIMS];
                double d[NDIMS];
                for (int i = 0; i < NDIMS; ++i) d[i] = b[i] - a[i];
                covered += clip_segment(var, a, d, 0, 1, previous.length[k], w, placed);
            }
        }
        if (covered < (1 - 1e-9) * (zt - zb)) {
            std::cerr << "Error: the screen of well " << w << " is not entirely inside the mesh.\n";
            die(EXIT_CONFIG_VALUE);
        }
    }
    if (ends) merge_pieces(placed);
    var.well_pieces = placed;
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

// Integral of T - T0 from the lowest node level up to each level, from the
// horizontal mean temperature of each level. One extra level at each end
// extends the end temperatures beyond the mesh, e.g. to a water table above it.
static void buoyancy_profile(const Param &param, const Variables &var, double alpha_f,
                             double_vec &zs, double_vec &area)
{
    const double t0 = param.bc.surface_temperature;
    std::vector<std::pair<double,double>> nodes(var.nnode);
    for (int n = 0; n < var.nnode; ++n) {
        const double t = (*var.temperature)[n];
        if (!(1 - alpha_f * (t - t0) > 0))
            die(EXIT_CONFIG_VALUE, "Buoyant fluid density is not positive at the initial temperature.");
        nodes[n] = std::make_pair((*var.coord)[n][NDIMS-1], t);
    }
    std::sort(nodes.begin(), nodes.end());
    const double tol = 1e-9 * (1 + nodes.back().first - nodes.front().first);
    double_vec ts;
    int_vec count;
    for (int n = 0; n < var.nnode; ++n) {
        if (zs.empty() || nodes[n].first - zs.back() > tol) {
            zs.push_back(nodes[n].first);
            ts.push_back(0);
            count.push_back(0);
        }
        ts.back() += nodes[n].second;
        ++count.back();
    }
    const int nlevel = zs.size();
    area.assign(nlevel, 0);
    for (int k = 0; k < nlevel; ++k) ts[k] = ts[k] / count[k] - t0;
    for (int k = 1; k < nlevel; ++k)
        area[k] = area[k-1] + 0.5 * (ts[k-1] + ts[k]) * (zs[k] - zs[k-1]);
    zs.insert(zs.begin(), zs.front() - 1);
    area.insert(area.begin(), -ts.front());
    zs.push_back(zs.back() + 1);
    area.push_back(area.back() + ts.back());
}

static double profile_integral(const double_vec &zs, const double_vec &area, double z)
{
    if (zs.empty()) return 0;
    const int k = std::upper_bound(zs.begin()+1, zs.end()-1, z) - zs.begin();
    return area[k-1] + (area[k] - area[k-1]) / (zs[k] - zs[k-1]) * (z - zs[k-1]);
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

    // Buoyant water: dp/dz = -rho_f0*g*(1 - alpha_f*(T - T0)).
    const double alpha_f = param.mat.fluid_alpha[param.mat.mattype_ref];
    double_vec zs, area;
    if (hydraulic_buoyancy(param) && alpha_f != 0)
        buoyancy_profile(param, var, alpha_f, zs, area);
    const double table = profile_integral(zs, area, param.ic.water_table_z);

    for (int i = 0; i < var.nnode; ++i) {
        ConstArrayAccessor coord = (*var.coord)[i];
        double z = coord[NDIMS - 1];

        const int side = pressure_boundary(i, var);
        if (side >= 0 && var.hbc_types[side] == 2) continue;
        // Pressure is compression-positive; z increases upward.
        ppressure[i] = -rho_fluid * gravity * (z - param.ic.water_table_z);
        if (!zs.empty())
            ppressure[i] -= rho_fluid * gravity * alpha_f * (table - profile_integral(zs, area, z));

        // Add excess pore pressure for non-boundary nodes
        if (side < 0) ppressure[i] += loading_efficiency * loading;
        dppressure[i] = 0.0;
    }
}

void prepare_pore_recovery(const Variables &var)
{
    if (!var.pore_dilation) return;
    int invalid = 0;
#ifndef ACC
    #pragma omp parallel for default(none) shared(var) reduction(|:invalid)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector reduction(|:invalid)
    for (int e = 0; e < var.nelem; ++e) {
        const double rate = var.mat->pore_healing_rate(e);
        const double change = (*var.pore_dilation)[e]*(-std::expm1(-var.dt*rate));
        const double loss = (var.mat->alpha_biot(e)-var.mat->phi(e))*std::expm1(change);
        (*var.pore_recovery)[e] = change;
        if (!(rate >= 0 && std::isfinite(rate) && change >= 0 && std::isfinite(change) &&
              loss >= 0 && std::isfinite(loss))) invalid = 1;
    }
    if (invalid) die(EXIT_RUNTIME_NAN,"Invalid pore recovery temperature or rate.");
}

void update_porosity(const Variables &var)
{
    if (!var.pore_strain) return;
    int invalid = 0;
#ifndef ACC
    #pragma omp parallel for default(none) shared(var) reduction(|:invalid)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector reduction(|:invalid)
    for (int e = 0; e < var.nelem; ++e) {
        ConstConnAccessor conn = (*var.connectivity)[e];
        double dp = 0;
        for (int i = 0; i < NODES_PER_ELEM; ++i)
            dp -= (*var.dppressure)[conn[i]] / NODES_PER_ELEM;
        const double eps = (*var.pore_strain)[e] + var.dt * trace((*var.strain_rate)[e])
                         + var.mat->beta_mineral(e) * dp
                         - var.mat->alpha_solid(e) * hydraulic_temperature_change(var, e)
                         - (var.pore_recovery ? (*var.pore_recovery)[e] : 0);
        const double phi = biot_porosity(var.mat->phi0(e), var.mat->alpha_biot(e), eps);
        (*var.etmp)[e] = eps;
        if (!biot_porosity_valid(*var.mat, e, phi) || !std::isfinite(eps))
            invalid = 1;
        if (var.pore_plastic_increment && !std::isfinite((*var.pore_plastic_increment)[e])) invalid = 1;
    }
    if (invalid)
        die(EXIT_RUNTIME_NAN, "Porosity leaves the Biot model range; reduce loading or the physical timestep.");
#ifndef ACC
    #pragma omp parallel for default(none) shared(var)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector
    for (int e = 0; e < var.nelem; ++e)
    {
        (*var.pore_strain)[e] = (*var.etmp)[e];
        if (var.pore_dilation) update_pore_dilation(var, e, (*var.pore_recovery)[e]);
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

    double diff_max_local = 0, source = 0, exchanged = 0;
    const bool has_retention = var.mat->has_retention;
    const bool strain_form = hydraulic_strain_form(param);
    const bool plastic_source = has_plastic_source(param);
    const bool strain_source = strain_form || plastic_source;
    const int mpi_ee = mpi_element_end(var), mpi_ne = mpi_node_end(var);

#ifndef ACC
    #pragma omp parallel for default(none) shared(var, ppressure, tmp_result, stress, old_mean_stress, param) \
        firstprivate(dt, has_retention, strain_form, plastic_source, strain_source, mpi_ee) reduction(max:diff_max_local) reduction(+:source, exchanged)
#endif
        #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector firstprivate(dt, has_retention, strain_form, plastic_source, strain_source, mpi_ee) reduction(max:diff_max_local) reduction(+:source, exchanged) async
    for (int e = 0; e < mpi_ee; e++) {
        ConstConnAccessor conn = (*var.connectivity)[e];
        ElemCacheAccessor tr = tmp_result[e];
        double current_mean_stress = trace(stress[e])/NDIMS;
        double mean_stress_change = current_mean_stress - old_mean_stress[e];

        double alpha_b = var.mat->alpha_biot(e);         // Biot coefficient
        double bulk_comp = hydraulic_mean_stress_compliance(*var.mat, e);


        // Darcy mobility k/mu in pressure form, well-defined when gravity is zero;
        // it is algebraically equivalent to the head form
        //   (k * gamma_w / mu) * grad(p / gamma_w + z).
        double mobility = hydraulic_mobility(*var.mat, e);
        // Bishop chi_e: the mean of the element's nodal Se, the same chi as in the
        // effective stress (rheology.cxx); exactly 1 when saturated.
        double chi = 1;
        if (has_retention) {
            // kr is the mean of the nodal values; it is exactly 1 when saturated.
            const int_vec &markers = (*var.elemmarkers)[e];
            int total = 0;
            mobility *= element_mean_kr(*var.mat, markers, conn, ppressure);
            element_mean_saturation(*var.mat, markers, conn, ppressure, &chi);
            for (int m = 0; m < var.mat->nmat; ++m) total += markers[m];
            tr[pore_slot] = (*var.volume)[e] * var.mat->phi(e) / (NODES_PER_ELEM * total);
        }
        const double kv_pressure = mobility * (*var.volume)[e];

        const double diff_e = hydraulic_diffusivity(param, *var.mat, e);
        diff_max_local = std::max(diff_max_local, diff_e);

        // volume term (poroelastic effect)
        double pe = 0;
        if (param.control.has_poroelastic_pressure_feedback && strain_source)
            // Match the volumetric rate used by the anti-locking constitutive update.
            pe = alpha_b * (*var.edvoldt)[e] * (*var.volume)[e] / NODES_PER_ELEM;
        else if (param.control.has_poroelastic_pressure_feedback)
            pe = alpha_b * mean_stress_change * bulk_comp * (*var.volume)[e] / NODES_PER_ELEM / var.dt;
        if (chi != 1) pe *= chi;
        // Bishop fixed-stress correction L*(chi*E - D), frozen for the physical step
        // (pressure_increments). An exact zero is skipped, so a saturated element
        // adds no term.
        if (plastic_source) {
            // The strain source already includes the plastic volume change;
            // only the previous raw pressure increment belongs to the buffer.
            pe -= hydraulic_fixed_stress_lag(*var.mat, e) * tr[pressure_d_slot] * (*var.volume)[e] /
                  (NODES_PER_ELEM * var.dt);
        } else if (has_retention && param.control.has_poroelastic_pressure_feedback && !strain_form) {
            const double correction = chi * tr[bishop_e_slot] - tr[pressure_d_slot];
            if (correction != 0)
                pe += hydraulic_fixed_stress_lag(*var.mat, e) * correction * (*var.volume)[e] /
                      NODES_PER_ELEM / var.dt;
        }
        if (var.old_temperature) {
            const double phi = var.mat->phi(e), solid = var.mat->alpha_solid(e);
            const double beta = pore_weighted(phi, alpha_b, var.mat->alpha_fluid(e), solid);
            const double thermal = strain_source ? beta :
                beta - alpha_b * bulk_comp * var.mat->bulkm(e) * solid;
            pe -= thermal * hydraulic_temperature_change(var, e) * (*var.volume)[e] /
                  (NODES_PER_ELEM * var.dt);
        }
        if (var.pore_recovery) {
            // Intrinsic pore closure expels water; it is not skeleton strain.
            const double loss = (alpha_b-var.mat->phi(e))*std::expm1((*var.pore_recovery)[e]);
            pe -= loss*(*var.volume)[e]/(NODES_PER_ELEM*var.dt);
        }
        source += pe * NODES_PER_ELEM * dt;
        if (var.matrix_pressure) {
            double pressure = 0;
            for (int j = 0; j < NODES_PER_ELEM; ++j) pressure += ppressure[conn[j]] / NODES_PER_ELEM;
            const double h = var.mat->exchange(e), capacity = var.mat->matrix_capacity(e);
            if (h > 0) {
                const double q = h * (pressure - (*var.matrix_pressure)[e]);
                (*var.matrix_pressure)[e] += dt*q/capacity;
                exchanged += dt*q*(*var.volume)[e];
                pe += q*(*var.volume)[e]/NODES_PER_ELEM;
            }
        }
        double grad[NDIMS][NODES_PER_ELEM];
        get_shape_gradients(var, e, grad);
        // One potential gradient and tensor product per element, shared by its nodes.
        double drive[NDIMS], flux[NDIMS] = {};
        hydraulic_drive(param, var, ppressure, e, grad, drive);
        for (int d = 0; d < NDIMS; ++d)
            for (int k = 0; k < NDIMS; ++k)
                flux[d] += kv_pressure * var.mat->perm_tensor[d][k] * drive[k];
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
    flux_facets(var,[&](int side, int e, int f, const int *, double area) {
        const double q = boundary_ratio * param.bc.hbc_flux[side] * area;
        if (e < mpi_ee) prescribed_flux += q * dt;
        for (int k = 0; k < NODES_PER_FACET; ++k)
            tmp_result[e][NODE_OF_FACET[f][k]] += q / NODES_PER_FACET;
    });
    double stored = exchanged, wells = 0, boundary = -prescribed_flux;
    int invalid = 0;

#ifndef ACC
    #pragma omp parallel for default(none) shared(param, var, tdot, ppressure, dppressure, tmp_result) \
        firstprivate(has_wells, has_points, well_ratio, pressure_ratio, dt, has_retention, mpi_ne) \
        reduction(+:stored, wells, boundary) reduction(|:invalid)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector reduction(+:stored, wells, boundary) reduction(|:invalid) firstprivate(has_points, mpi_ne)
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
        double next = p0 - increment, change = -mass * increment;
        if (has_retention) {
            // Below the patch's saturation threshold the content, not the
            // pressure, is advanced; saturated nodes keep the line above.
            const PorePatch pore = {*var.mat, *var.elemmarkers, tmp_result,
                                    var.support.patch(n), var.support.size(n)};
            if (std::min(p0, next) < patch_p_sat(pore)) {
                if (side < 0 && mass > 0) {
                    next = advance_content(pore, mass, p0, -(tdot[n] * dt), invalid);
                    increment = p0 - next;
                }
                change = content_change(pore, mass, p0, next);
            }
        }
        ppressure[n] = next;
        dppressure[n] += increment;
        stored += change;
        if (has_wells) wells += well_ratio * (*var.fluid_source)[n] * dt;
        if (has_points) wells -= (*var.fluid_point_source)[n]*dt;
        if (side >= 0) boundary += change + tdot[n]*dt;

    }

    #pragma acc wait
    if (mpi_any(var, invalid))
        die_together(EXIT_RUNTIME_NAN, "Unsaturated content update did not converge; check the pressure state.");
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
    // Gershgorin bound on the lumped operator, including rotated tensors.
    // Reuse assembly scratch; no pressure or material history is copied.
    double limit = std::numeric_limits<double>::max();
    const int mpi_ee = mpi_element_end(var);
#ifndef ACC
    #pragma omp parallel for default(none) shared(var, work) firstprivate(mpi_ee) reduction(min:limit)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector firstprivate(mpi_ee) reduction(min:limit)
    for (int e = 0; e < mpi_ee; ++e) {
        double g[NDIMS][NODES_PER_ELEM];
        get_shape_gradients(var, e, g);
        const double kv = hydraulic_mobility(*var.mat, e) * (*var.volume)[e];
        const double exchange = var.matrix_pressure ? var.mat->exchange(e) : 0;
        if (exchange > 0) limit = std::min(limit, var.mat->matrix_capacity(e)/exchange);
        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            double row = 0;
            for (int j = 0; j < NODES_PER_ELEM; ++j)
                row += std::abs(permeability_coupling(*var.mat, g, i, j));
            work[e][i] = kv * row + exchange*(*var.volume)[e]/NODES_PER_ELEM;
        }
    }
    mpi_sync(var, work, MpiEntity::Element, NODES_PER_ELEM);
    return nodal_step_limit(var, work, *var.hmass, limit);
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
    if (has_bishop_correction(param, var) || has_plastic_source(param))
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

enum FluidReference { fluid_phi, fluid_pressure, fluid_temperature, fluid_volume,
                      fluid_log_volume, fluid_pore_strain, fluid_recovery };

double effective_fluid_volume(const Param &param, const Variables &var, int e)
{
    return (*var.volume)[e]*(param.control.has_moving_mesh ? 1 : std::exp(trace((*var.strain)[e])));
}

#pragma acc routine seq
bool valid_fluid(const MatProps &mat, int e, const FluidState &state)
{
    return biot_porosity_valid(mat,e,state.porosity) &&
           state.density > 0 && std::isfinite(state.density) && state.pressure_slope >= 0 &&
           std::isfinite(state.pressure_slope) && std::isfinite(state.temperature_slope);
}

FluidState remap_fluid_state(const Variables &var, int e, int n, double p, double t)
{
    const auto state = nodal_fluid_state(var,e,n,p,t);
    if (!valid_fluid(*var.mat,e,state)) die(EXIT_RUNTIME_NAN,"Invalid fluid remap state.");
    return state;
}

double schedule_peak(const double_vec &time, const double_vec &ratio, double a, double b)
{
    double peak = std::max(std::fabs(schedule_mean(time,ratio,a,a)),std::fabs(schedule_mean(time,ratio,b,b)));
    for (std::size_t i = 0; i < time.size(); ++i)
        if (time[i]*YEAR2SEC >= a && time[i]*YEAR2SEC <= b) peak = std::max(peak,std::fabs(ratio[i]));
    return peak;
}

struct FluidSource {
    const Param &param;
    const Variables &var;
    double_vec &mass, &heat;
    bool bound;
};

void fluid_point_source(int point, int n, double volume_rate, void *context)
{
    auto &source = *static_cast<FluidSource *>(context);
    const auto &var = source.var;
    const int e = var.support.patch(n)[0];
    const double t = volume_rate > 0 && !source.param.injection.temperature.empty() ?
        source.param.injection.temperature[point] : (*var.temperature)[n];
    const double rate = volume_rate*fluid_density(*var.mat,e,(*var.ppressure)[n],t,
                                                 var.fluid_reference_temperature);
    source.mass[n] += source.bound ? std::fabs(rate) : rate;
    const double heat = var.mat->cp_fluid(e)*t*rate;
    source.heat[n] += source.bound ? std::fabs(heat) : heat;
}

} // namespace

namespace {
// Fluid mass of node n's patch at pressure p and temperature t, with its
// pressure and (optional) temperature slopes, or NaN for an invalid state.
// reference[e] is element e's FluidReference row; fraction scales its volume
// and recovery changes.
template <class Reference>
double patch_fluid_content(const Variables &var, const Reference &reference, int n, double p,
                           double t, double fraction, double &slope, double *thermal_slope)
{
    double value = 0;
    slope = 0;
    if (thermal_slope) *thermal_slope = 0;
    const int *patch = var.support.patch(n);
    for (int k = 0; k < var.support.size(n); ++k) {
        const int e = patch[k];
        const auto r = reference[e];
        const auto state = fluid_state(*var.mat,e,p,t,var.fluid_reference_temperature,r[fluid_phi],
            r[fluid_pressure],r[fluid_temperature],fraction*r[fluid_log_volume],fraction*r[fluid_recovery]);
        if (!valid_fluid(*var.mat,e,state)) return std::numeric_limits<double>::quiet_NaN();
        const double weight = r[fluid_volume]*std::exp(fraction*r[fluid_log_volume])/NODES_PER_ELEM;
        value += weight*state.density*state.porosity;
        slope += weight*state.pressure_slope;
        if (thermal_slope) *thermal_slope += weight*state.temperature_slope;
    }
    return value;
}

// FluidReference rows of the current state, without volume or recovery change.
struct CurrentReference
{
    const Param &param;
    const Variables &var;
    std::array<double,fluid_recovery+1> operator[](int e) const
    {
        double p, t; fluid_means(var,e,p,t);
        return {var.mat->phi(e),p,t,effective_fluid_volume(param,var,e),0,0,0};
    }
};

double current_fluid_content(const Param &param, const Variables &var, int n, double pressure,
                             double temperature, double &mp, double &mt)
{
    return patch_fluid_content(var,CurrentReference{param,var},n,pressure,temperature,0,mp,&mt);
}

double_vec lumped_volume(const Variables &var)
{
    double_vec volume(var.nnode,0);
    for (int e = 0; e < var.nelem; ++e)
        for (int j = 0; j < NODES_PER_ELEM; ++j)
            volume[(*var.connectivity)[e][j]] += (*var.volume)[e]/NODES_PER_ELEM;
    return volume;
}
} // namespace

double nodal_fluid_mass(const Param &param, const Variables &var, int n)
{
    double mp, mt;
    const double mass = current_fluid_content(param,var,n,(*var.ppressure)[n],(*var.temperature)[n],mp,mt);
    if (!(mass > 0 && std::isfinite(mass))) die(EXIT_RUNTIME_NAN,"Invalid saturated fluid content.");
    return mass;
}

void initialize_fluid_content(const Param &param, Variables &var)
{
    if (!var.fluid_split) return;
    #pragma acc wait
    for (int n = 0; n < var.nnode; ++n) {
        double mass = 0, capacity = 0;
        const int *patch = var.support.patch(n);
        for (int k = 0; k < var.support.size(n); ++k) {
            const int e = patch[k];
            double p, t; fluid_means(var,e,p,t);
            const auto state = nodal_fluid_state(var,e,n,p,t);
            if (!valid_fluid(*var.mat,e,state) || !((*var.temperature)[n] > 0))
                die(EXIT_RUNTIME_NAN,"Invalid initial finite-density fluid state.");
            const double solid = effective_fluid_volume(param,var,e)/NODES_PER_ELEM*
                                 (1-state.porosity)*var.mat->rho(e);
            if (!(solid > 0 && std::isfinite(solid) && var.mat->cp(e) >= 0 && std::isfinite(var.mat->cp(e))))
                die(EXIT_RUNTIME_NAN,"Invalid initial solid mass or heat capacity.");
            mass += solid;
            capacity += solid*var.mat->cp(e);
        }
        (*var.solid_mass)[n] = mass;
        (*var.solid_capacity)[n] = capacity;
        (*var.fluid_split)[n] = 0;
        (*var.pending_temperature)[n] = 0;
    }
}

double_vec hydraulic_reference(const Param &param)
{
    // These quantities define the saved mass/energy representation. Changing
    // a boundary, permeability or timestep must not redefine that inventory.
    double_vec law = {double(param.mat.nmat),double(param.control.has_moving_mesh),
        double(param.control.has_thermal_pressurization),param.mat.fluid_rho0[0],
        param.mat.fluid_bulk_modulus[0],param.mat.fluid_alpha[0],param.mat.fluid_heat_capacity[0],
        double(param.control.has_PT)};
    for (int m = 0; m < param.mat.nmat; ++m)
        law.insert(law.end(),{param.mat.porosity[m],param.mat.biot_coeff[m],
                             param.mat.bulk_modulus_s[m],param.mat.alpha[m]});
    return law;
}

#pragma acc routine seq
double HydraulicStep::content(const Param &param, const Variables &var, int n, double p,
                              double t, double fraction, double &slope, double *thermal_slope) const
{
    return patch_fluid_content(var,reference,n,p,t,fraction,slope,thermal_slope);
}

void HydraulicStep::begin(const Param &param, const Variables &var)
{
    if (!var.fluid_split) return;
    #pragma acc wait
    interval = var.dt;
    start_time = var.time-interval;
    pressure = *var.ppressure;
    temperature = *var.temperature;
    initial_split = *var.fluid_split;
    if (!param.control.has_PT && param.control.has_moving_mesh) {
        initial_velocity.resize(var.nnode,false);
        std::copy(var.vel->begin(),var.vel->end(),initial_velocity.begin());
        initial_inertia = *var.mass;
    }
    if (reference.size() != static_cast<std::size_t>(var.nelem)) {
        reference.resize(var.nelem);
        transport.resize(var.nelem);
        heat_transport.resize(var.nelem);
    }
    for (int e = 0; e < var.nelem; ++e) {
        double p, t; fluid_means(var,e,p,t);
        auto r = reference[e];
        r[fluid_phi] = var.mat->phi(e);
        r[fluid_pressure] = p;
        r[fluid_temperature] = t;
        r[fluid_volume] = effective_fluid_volume(param,var,e);
        r[fluid_log_volume] = 0;
        r[fluid_pore_strain] = (*var.pore_strain)[e];
        r[fluid_recovery] = var.pore_recovery ? (*var.pore_recovery)[e] : 0;
    }
    initial_mass.resize(var.nnode);
    initial_energy.resize(var.nnode);
    lag.assign(var.nnode,0);
    for (int n = 0; n < var.nnode; ++n) {
        double slope;
        initial_mass[n] = content(param,var,n,pressure[n],temperature[n],0,slope);
        initial_energy[n] = ((*var.solid_capacity)[n]+var.mat->cp_fluid(0)*initial_mass[n])*temperature[n];
        if (!(initial_mass[n] > 0 && std::isfinite(initial_mass[n]) && temperature[n] > 0))
            die(EXIT_RUNTIME_NAN,"Invalid finite-density physical-start state.");
        if (!param.control.has_PT) continue;
        const int *patch = var.support.patch(n);
        for (int k = 0; k < var.support.size(n); ++k) {
            const int e = patch[k];
            const double alpha = var.mat->alpha_biot(e);
            lag[n] += reference[e][fluid_volume]/NODES_PER_ELEM*
                fluid_density(*var.mat,e,pressure[n],temperature[n],var.fluid_reference_temperature)*
                alpha*alpha*hydraulic_mean_stress_compliance(*var.mat,e);
        }
    }
}

#pragma acc routine seq
double HydraulicStep::invert(const Param &param, const Variables &var, int n, double fraction,
                             double value, double guess, double t) const
{
    const int *patch = var.support.patch(n);
    const int first = patch[0];
    const double beta = var.mat->beta_fluid(first);
    // Keep exp itself normal as well as the resulting density. Some math
    // runtimes flush a subnormal exponential before multiplication by rho0.
    const double exponent_floor = std::log(std::numeric_limits<double>::min())+std::log(16.0)-
        std::min(0.0,std::log(var.mat->rho_fluid_reference(first)));
    double lower = (exponent_floor+
        var.mat->alpha_fluid(first)*(t-var.fluid_reference_temperature))/beta;
    double scale = 1;
    for (int k = 0; k < var.support.size(n); ++k) {
        const int e = patch[k];
        const auto r = reference[e];
        const double alpha = var.mat->alpha_biot(e), grain = 1/var.mat->beta_mineral(e);
        const double gap = alpha-r[fluid_phi];
        if (gap > 0)
            lower = std::max(lower,r[fluid_pressure]+grain*(std::log(gap/alpha)-
                fraction*r[fluid_log_volume]+fraction*r[fluid_recovery]+
                var.mat->alpha_solid(e)*(t-r[fluid_temperature])));
        scale = std::max(scale,grain);
    }
    const double roundoff = 64*std::numeric_limits<double>::epsilon();
    lower += roundoff*std::max(scale,std::fabs(lower));
    auto residual = [&](double p, double &slope) {
        const double m = content(param,var,n,p,t,fraction,slope);
        slope += lag[n];
        return m+lag[n]*(p-pressure[n])-value;
    };
    double slope, f = residual(lower,slope);
    const double tolerance = roundoff*std::max(std::fabs(value),initial_mass[n]);
    if (!std::isfinite(f) || f > tolerance) return std::numeric_limits<double>::quiet_NaN();
    double upper = std::max(pressure[n],lower+std::max(1e4,1e-3/beta));
    for (int i = 0; i < 100; ++i) {
        f = residual(upper,slope);
        if (!std::isfinite(f) || f >= 0) break;
        upper = lower+2*(upper-lower);
        if (i == 99) return std::numeric_limits<double>::quiet_NaN();
    }
    double p = std::max(lower,std::min(upper,guess));
    for (int i = 0; i < 100; ++i) {
        f = residual(p,slope);
        if (std::isfinite(f) && std::fabs(f) <= tolerance) return p;
        if (!std::isfinite(f) || f > 0) upper = p;
        else lower = p;
        const double next = p-f/slope;
        p = std::isfinite(next) && next > lower && next < upper ? next : lower+.5*(upper-lower);
    }
    return std::numeric_limits<double>::quiet_NaN();
}

#pragma acc routine seq
bool HydraulicStep::equilibrate(const Param &param, const Variables &var, int n, double fraction,
                                double mass_target, double energy_target, double &p, double &t) const
{
    const double cf = var.mat->cp_fluid(0), solid = (*var.solid_capacity)[n];
    const double tolerance = 128*std::numeric_limits<double>::epsilon()*
        std::max(std::fabs(energy_target),std::fabs(initial_energy[n]));
    if (!(energy_target > 0 && std::isfinite(energy_target))) return false;
    for (int i = 0; i < 80; ++i) {
        p = invert(param,var,n,fraction,mass_target,p,t);
        double mp, mt;
        const double m = content(param,var,n,p,t,fraction,mp,&mt);
        if (!(m > 0 && std::isfinite(m) && std::isfinite(p))) return false;
        const double capacity = solid+cf*m;
        const double residual = capacity*t-energy_target;
        const double tangent = capacity+cf*t*mt*lag[n]/(mp+lag[n]);
        if (!(t > 0 && tangent > 0 && std::isfinite(tangent))) return false;
        if (std::fabs(residual) <= tolerance) return true;
        double increment = -residual/tangent;
        bool accepted = false;
        for (int j = 0; j < 40; ++j) {
            const double next_t = t+increment;
            if (next_t > 0 && std::isfinite(next_t)) {
                const double next_p = invert(param,var,n,fraction,mass_target,p,next_t);
                double slope;
                const double next_mass = content(param,var,n,next_p,next_t,fraction,slope);
                if (next_mass > 0 && std::isfinite(next_mass) &&
                    std::fabs((solid+cf*next_mass)*next_t-energy_target) < std::fabs(residual)) {
                    t = next_t; p = next_p; accepted = true; break;
                }
            }
            increment *= .5;
        }
        if (!accepted) return false;
    }
    return false;
}

void HydraulicStep::advance(const Param &param, Variables &var)
{
    if (!var.fluid_split) return;
    #pragma acc wait
    var.fluid_pressure_consumed = false;
    if (pressure.size() != static_cast<std::size_t>(var.nnode) || reference.size() != static_cast<std::size_t>(var.nelem))
        die(EXIT_INTERNAL_ASSERT,"Finite-density physical step requires fixed topology until transport accepts.");
    double deformation = 0;
    for (int e = 0; e < var.nelem; ++e) {
        auto r = reference[e];
        r[fluid_log_volume] = std::log(effective_fluid_volume(param,var,e)/r[fluid_volume]);
        deformation = std::max(deformation,std::fabs(r[fluid_log_volume])+r[fluid_recovery]);
    }
    mass = target = initial_mass;
    energy = initial_energy;
    const bool thermal = has_temperature_transport(param);
    const double cf = var.mat->cp_fluid(0);
    pending_energy.assign(var.nnode,0);
    for (int n = 0; n < var.nnode; ++n)
        if (var.pending_heat) {
            const int *patch = var.support.patch(n);
            for (int k = 0; k < var.support.size(n); ++k)
                pending_energy[n] += (*var.pending_heat)[patch[k]]/NODES_PER_ELEM;
        }
    for (double_vec *work : {&rate,&row,&outgoing,&boundary_rate,&source_rate,&candidate,&candidate_mass,
                            &energy_rate,&heat_row,&heat_source,&source_heat,&boundary_heat,
                            &candidate_temperature,&candidate_energy})
        work->resize(var.nnode);
    double elapsed = 0, boundary = 0, extraction = 0;
    double supplied_heat = 0, well_heat = 0, carried_heat = 0, thermal_reaction = 0;
    std::uint64_t steps = 0;
    while (elapsed < interval) {
        const double fraction = elapsed/interval, time = start_time+elapsed;
        std::fill(transport.data(),transport.data()+transport.num_elements(),0.0);
        std::fill(heat_transport.data(),heat_transport.data()+heat_transport.num_elements(),0.0);
        for (int e = 0; e < var.nelem; ++e) {
            double g[NDIMS][NODES_PER_ELEM], q[NDIMS], rho[NODES_PER_ELEM], max_rho = 0, mean_rho = 0;
            get_shape_gradients(var,e,g);
            const auto conn = (*var.connectivity)[e];
            auto work = transport[e];
            auto heat = heat_transport[e];
            for (int i = 0; i < NODES_PER_ELEM; ++i) {
                rho[i] = nodal_fluid_density(var,e,conn[i]);
                max_rho = std::max(max_rho,rho[i]);
                mean_rho += rho[i]/NODES_PER_ELEM;
            }
            const double mobility = hydraulic_mobility(*var.mat,e);
            hydraulic_flux(param,var,e,g,mobility,mean_rho,q);
            double stiffness[NODES_PER_ELEM][NODES_PER_ELEM], all = 0;
            for (int i = 0; i < NODES_PER_ELEM; ++i)
                for (int j = 0; j < NODES_PER_ELEM; ++j) {
                    stiffness[i][j] = std::fabs(permeability_coupling(*var.mat,g,i,j));
                    all += stiffness[i][j];
                }
            const double volume = (*var.volume)[e], beta = var.mat->beta_fluid(e);
            for (int i = 0; i < NODES_PER_ELEM; ++i) {
                double local = 0;
                for (int j = 0; j < NODES_PER_ELEM; ++j) local += stiffness[i][j];
                work[NODES_PER_ELEM+i] += max_rho*mobility*volume*(local+all/NODES_PER_ELEM);
                for (int j = i+1; j < NODES_PER_ELEM; ++j) {
                    double flux = 0, gravity_flux = 0;
                    for (int d = 0; d < NDIMS; ++d) {
                        const double face = -volume/NODES_PER_ELEM*(g[d][i]-g[d][j]);
                        flux += face*q[d];
                        gravity_flux += face*(-mobility*var.mat->perm_tensor[d][NDIMS-1]*mean_rho*param.control.gravity);
                    }
                    flux *= flux >= 0 ? rho[i] : rho[j];
                    work[i] += flux; work[j] -= flux;
                    work[2*NODES_PER_ELEM+i] += std::max(0.0,flux);
                    work[2*NODES_PER_ELEM+j] += std::max(0.0,-flux);
                    const double density_row = beta*(std::fabs(flux)+max_rho*std::fabs(gravity_flux));
                    work[NODES_PER_ELEM+i] += density_row;
                    work[NODES_PER_ELEM+j] += density_row;
                    if (thermal) {
                        const double ti = (*var.temperature)[conn[i]], tj = (*var.temperature)[conn[j]];
                        const double advected = cf*flux*(flux >= 0 ? ti : tj);
                        double conductance = 0;
                        if (param.control.has_thermal_diffusion)
                            for (int d = 0; d < NDIMS; ++d)
                                conductance -= volume*var.mat->k(e)*g[d][i]*g[d][j];
                        const double heat_flux = advected+conductance*(ti-tj);
                        heat[i] += heat_flux; heat[j] -= heat_flux;
                        const double heat_bound = cf*std::fabs(flux)*(1+std::fabs(var.mat->alpha_fluid(e))*std::max(ti,tj))+
                                                  2*std::fabs(conductance);
                        heat[NODES_PER_ELEM+i] += heat_bound;
                        heat[NODES_PER_ELEM+j] += heat_bound;
                    }
                }
            }
        }
        const double end = start_time+interval;
        const bool wells = param.bc.num_wells > 0, points = param.injection.enabled;
        const bool scheduled = !param.bc.hbc_schedule_ratio.empty();
        const double well_peak = wells ? schedule_peak(param.bc.well_schedule_time_in_yr,param.bc.well_schedule_ratio,time,end) : 0;
        const double boundary_peak = scheduled ? schedule_peak(param.bc.hbc_schedule_time_in_yr,param.bc.hbc_schedule_ratio,time,end) : 1;
        std::fill(source_rate.begin(),source_rate.end(),0.0);
        std::fill(source_heat.begin(),source_heat.end(),0.0);
        if (points) {
            FluidSource source{param,var,source_rate,source_heat,true};
            deposit_point_sources(param.injection,var.fluid_points,var,time,end,fluid_point_source,&source,nullptr,true);
        }
        double limit = deformation > .05 ? interval*.05/deformation : interval;
        for (int n = 0; n < var.nnode; ++n) {
            double slope;
            double mt;
            content(param,var,n,(*var.ppressure)[n],(*var.temperature)[n],fraction,slope,&mt);
            rate[n] = row[n] = outgoing[n] = 0;
            energy_rate[n] = heat_row[n] = 0;
            const int *patch = var.support.patch(n), *local = var.support.local(n);
            for (int k = 0; k < var.support.size(n); ++k) {
                const auto work = transport[patch[k]];
                rate[n] -= work[local[k]];
                row[n] += work[NODES_PER_ELEM+local[k]];
                outgoing[n] += work[2*NODES_PER_ELEM+local[k]];
                energy_rate[n] -= heat_transport[patch[k]][local[k]];
                heat_row[n] += heat_transport[patch[k]][NODES_PER_ELEM+local[k]];
            }
            const int e = patch[0];
            const double density = nodal_fluid_density(var,e,n);
            const double source = density*(wells ? std::fabs((*var.fluid_source)[n])*well_peak : 0)+source_rate[n];
            outgoing[n] += source;
            row[n] += var.mat->beta_fluid(e)*source;
            if (!(std::isfinite(rate[n]) && row[n] >= 0 && std::isfinite(row[n]) &&
                  outgoing[n] >= 0 && std::isfinite(outgoing[n]) && slope >= 0 && std::isfinite(slope) && slope+lag[n] > 0))
                die(EXIT_RUNTIME_NAN,"Invalid finite-density transport coefficients.");
            if (row[n] > 0) limit = std::min(limit,.5*(slope+lag[n])/row[n]);
            if (outgoing[n] > 0) limit = std::min(limit,.5*mass[n]/outgoing[n]);
            if (thermal) {
                const double t = (*var.temperature)[n], capacity = (*var.solid_capacity)[n]+cf*mass[n];
                const double tangent = capacity+cf*t*mt*lag[n]/(slope+lag[n]);
                heat_row[n] += cf*source+std::fabs(var.mat->alpha_fluid(e))*source_heat[n];
                if (!(tangent > 0 && std::isfinite(tangent)))
                    die(EXIT_RUNTIME_NAN,"Finite-density sensible heat has nonpositive constrained capacity.");
                if (heat_row[n] > 0) limit = std::min(limit,.5*tangent/heat_row[n]);
            }
        }
        // Natural boundary fluxes use the same receiving/outgoing fluid density
        // as the mass reaction. Their peak bounds also cover schedule steps.
        flux_facets(var,[&](int side, int e, int, const int *nodes, double area) {
            const double volume_rate = std::fabs(param.bc.hbc_flux[side])*boundary_peak*area/NODES_PER_FACET;
            for (int j = 0; j < NODES_PER_FACET; ++j) {
                const int n = nodes[j];
                const double flux = volume_rate*nodal_fluid_density(var,e,n);
                if (flux > 0) {
                    double slope; content(param,var,n,(*var.ppressure)[n],(*var.temperature)[n],fraction,slope);
                    limit = std::min({limit,.5*mass[n]/flux,.5*(slope+lag[n])/(var.mat->beta_fluid(e)*flux)});
                }
            }
        });
        double dt = std::min(limit,interval-elapsed), next_fraction = 0;
        bool accepted = false;
        for (int attempt = 0; attempt < 50 && !accepted; ++attempt) {
            if (!(dt > 0 && std::isfinite(dt)) || elapsed+dt == elapsed)
                die(EXIT_RUNTIME_RESOURCE,"Finite-density substep exhausts physical-time resolution.");
            next_fraction = dt == interval-elapsed ? 1 : (elapsed+dt)/interval;
            const double stop = start_time+elapsed+dt;
            double well_ratio, boundary_ratio, pressure_ratio;
            schedule_ratios(param.bc,time,stop,well_ratio,boundary_ratio,pressure_ratio);
            std::fill(source_rate.begin(),source_rate.end(),0.0);
            std::fill(source_heat.begin(),source_heat.end(),0.0);
            if (points) {
                FluidSource source{param,var,source_rate,source_heat,false};
                deposit_point_sources(param.injection,var.fluid_points,var,time,stop,fluid_point_source,&source);
            }
            std::fill(heat_source.begin(),heat_source.end(),0.0);
            if (thermal) {
                if (param.heat_injection.enabled)
                    assemble_point_sources(param.heat_injection,var.heat_points,var,time,stop,heat_source);
                for (int n = 0; n < var.nnode; ++n) {
                    double power = 0, volume = 0;
                    const int *patch = var.support.patch(n);
                    for (int k = 0; k < var.support.size(n); ++k) {
                        const int e = patch[k];
                        power += (*var.volume)[e]*(*var.radiogenic_source)[e];
                        volume += (*var.volume)[e];
                    }
                    heat_source[n] += (*var.solid_mass)[n]*power/volume+pending_energy[n]/interval;
                }
            }
            std::fill(boundary_rate.begin(),boundary_rate.end(),0.0);
            std::fill(boundary_heat.begin(),boundary_heat.end(),0.0);
            flux_facets(var,[&](int side, int e, int, const int *nodes, double area) {
                const double q = boundary_ratio*param.bc.hbc_flux[side]*area/NODES_PER_FACET;
                for (int j = 0; j < NODES_PER_FACET; ++j) {
                    const int n = nodes[j];
                    const double t = (*var.temperature)[n];
                    const double flow = -q*nodal_fluid_density(var,e,n);
                    boundary_rate[n] += flow;
                    boundary_heat[n] += cf*t*flow;
                }
            });
            accepted = true;
            auto *const step = this;
#ifndef ACC
            #pragma omp parallel for default(none) shared(param,var) \
                firstprivate(step,cf,thermal,wells,well_ratio,dt,next_fraction,pressure_ratio) \
                reduction(&&:accepted) if(var.nnode >= 256)
#endif
            #pragma acc parallel loop copyin(step[0:1]) reduction(&&:accepted) if(des_acc_offload_enabled() && var.nnode >= 256)
            for (int n = 0; n < var.nnode; ++n) {
                const int e = var.support.patch(n)[0];
                const double t = (*var.temperature)[n];
                const double density = nodal_fluid_density(var,e,n);
                const double well = density*(wells ? -well_ratio*(*var.fluid_source)[n] : 0);
                step->source_rate[n] += well;
                step->source_heat[n] += cf*t*well;
                const double next_target = step->target[n]+dt*(step->rate[n]+step->source_rate[n]+step->boundary_rate[n]);
                const double predicted_energy = step->energy[n]+dt*(step->energy_rate[n]+step->source_heat[n]+step->boundary_heat[n]+step->heat_source[n]);
                const int side = pressure_boundary(n,var);
                double fixed_t;
                const bool fixed = thermal && prescribed_temperature(param,var,n,fixed_t);
                step->candidate_temperature[n] = fixed ? fixed_t : t;
                step->candidate[n] = (*var.ppressure)[n];
                const double augmented_target = next_target+next_fraction*step->initial_split[n];
                if (side >= 0) {
                    double fixed_p;
                    step->candidate[n] = prescribed_pressure(param,var,side,pressure_ratio,fixed_p) ? fixed_p : step->pressure[n];
                    if (thermal && !fixed) {
                        const double predicted_mass = step->mass[n]+dt*(step->rate[n]+step->source_rate[n]+step->boundary_rate[n]);
                        const double predicted_capacity = (*var.solid_capacity)[n]+cf*predicted_mass;
                        step->candidate_temperature[n] = predicted_capacity > 0 ? predicted_energy/predicted_capacity : -1;
                    }
                } else if (thermal && !fixed) {
                    if (!step->equilibrate(param,var,n,next_fraction,augmented_target,predicted_energy,
                                     step->candidate[n],step->candidate_temperature[n])) accepted = false;
                } else {
                    step->candidate[n] = step->invert(param,var,n,next_fraction,augmented_target,
                                          step->candidate[n],step->candidate_temperature[n]);
                }
                double slope;
                step->candidate_mass[n] = step->content(param,var,n,step->candidate[n],step->candidate_temperature[n],next_fraction,slope);
                step->candidate_energy[n] = ((*var.solid_capacity)[n]+cf*step->candidate_mass[n])*step->candidate_temperature[n];
                if (!(std::isfinite(step->candidate[n]) && step->candidate_mass[n] > 0 && std::isfinite(step->candidate_mass[n]) &&
                      step->candidate_temperature[n] > 0 && std::isfinite(step->candidate_temperature[n]) &&
                      std::isfinite(step->candidate_energy[n]))) accepted = false;
            }
            if (!accepted) dt *= .5;
        }
        if (!accepted) die(EXIT_RUNTIME_NAN,"Finite-density transport leaves its positive saturated state; reduce loading or the physical timestep.");
        for (int n = 0; n < var.nnode; ++n) {
            const double change = dt*(rate[n]+source_rate[n]+boundary_rate[n]);
            target[n] += change;
            boundary += dt*boundary_rate[n];
            extraction -= dt*source_rate[n];
            const double reaction = pressure_boundary(n,var) >= 0 ? candidate_mass[n]-(mass[n]+change) : 0;
            boundary += reaction;
            double accepted_energy = candidate_energy[n];
            if (thermal) {
                const double transported = dt*(energy_rate[n]+source_heat[n]+boundary_heat[n]+heat_source[n]);
                const double reaction_heat = cf*candidate_temperature[n]*reaction;
                supplied_heat += dt*heat_source[n];
                well_heat += dt*source_heat[n];
                carried_heat += dt*boundary_heat[n]+reaction_heat;
                double fixed_t;
                if (prescribed_temperature(param,var,n,fixed_t))
                    thermal_reaction += candidate_energy[n]-energy[n]-transported-reaction_heat;
                else
                    // Carry the integrated target, not each local Newton
                    // residual, through the remaining transport substeps.
                    accepted_energy = energy[n]+transported+reaction_heat;
            }
            mass[n] = candidate_mass[n];
            energy[n] = accepted_energy;
            (*var.ppressure)[n] = candidate[n];
            (*var.temperature)[n] = candidate_temperature[n];
        }
        for (int e = 0; e < var.nelem; ++e) {
            const auto r = reference[e];
            double p, t; fluid_means(var,e,p,t);
            (*var.pore_strain)[e] = r[fluid_pore_strain]+next_fraction*(r[fluid_log_volume]-r[fluid_recovery])+
                var.mat->beta_mineral(e)*(p-r[fluid_pressure])-var.mat->alpha_solid(e)*(t-r[fluid_temperature]);
        }
        elapsed += dt;
        if (++steps > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
            die(EXIT_RUNTIME_RESOURCE,"Finite-density transport needs too many explicit substeps.");
    }
    double stored = 0, buffer = 0, reset = 0, scale = 0;
    if (var.steps == 1 && steps > 1)
        std::printf("    Finite-density transport: %llu substeps in %.17g s\n",
                    static_cast<unsigned long long>(steps),interval);
    for (int n = 0; n < var.nnode; ++n) {
        const bool free = pressure_boundary(n,var) < 0;
        const double next_split = free ? lag[n]*((*var.ppressure)[n]-pressure[n]) : 0;
        if (!free) reset -= initial_split[n];
        (*var.fluid_split)[n] = next_split;
        (*var.dppressure)[n] = pressure[n]-(*var.ppressure)[n];
        (*var.pending_temperature)[n] = (*var.temperature)[n]-temperature[n];
        stored += nodal_fluid_mass(param,var,n)-initial_mass[n];
        buffer += next_split-initial_split[n];
        scale += initial_mass[n];
    }
    for (int e = 0; e < var.nelem; ++e)
        if (var.pore_dilation) update_pore_dilation(var,e,reference[e][fluid_recovery]);
    const double residual = stored+buffer+extraction-boundary-reset;
    const double tolerance = 512*std::numeric_limits<double>::epsilon()*scale+
        1e-10*(std::fabs(stored)+std::fabs(buffer)+std::fabs(extraction)+std::fabs(boundary)+std::fabs(reset));
    if (!std::isfinite(residual) || std::fabs(residual) > tolerance)
        die(EXIT_INTERNAL_ASSERT,"Finite-density transport fails its mass constraint.");
    *var.fluid_balance = {stored,buffer,extraction,boundary,reset,residual};
    if (thermal) {
        double change = 0, magnitude = 0;
        for (int n = 0; n < var.nnode; ++n) {
            const double actual = ((*var.solid_capacity)[n]+cf*nodal_fluid_mass(param,var,n))*(*var.temperature)[n];
            change += actual-initial_energy[n];
            magnitude += std::fabs(actual)+std::fabs(initial_energy[n]);
        }
        const double error = change-supplied_heat-well_heat-carried_heat-thermal_reaction;
        const double roundoff = 1024*std::numeric_limits<double>::epsilon()*magnitude+
            1e-10*(std::fabs(change)+std::fabs(supplied_heat)+std::fabs(well_heat)+std::fabs(carried_heat)+std::fabs(thermal_reaction));
        if (!std::isfinite(error) || std::fabs(error) > roundoff) {
            std::fprintf(stderr,"Sensible heat: change %.17g, supplied %.17g, source %.17g, carried %.17g, reaction %.17g, error %.17g, tolerance %.17g, substeps %llu\n",
                         change,supplied_heat,well_heat,carried_heat,thermal_reaction,error,roundoff,
                         static_cast<unsigned long long>(steps));
            die(EXIT_INTERNAL_ASSERT,"Finite-density transport fails its sensible-energy constraint.");
        }
        *var.energy_balance = {change,supplied_heat,well_heat,carried_heat,thermal_reaction,error};
        if (var.pending_heat) std::fill(var.pending_heat->begin(),var.pending_heat->end(),0.0);
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

// Mobile physical content, numerical fixed-stress buffer, immobile content.
void remap_content(const Param &param, const Variables &var, int e, int j, double value[3])
{
    const double alpha = var.mat->alpha_biot(e), phi = var.mat->phi(e);
    const bool feedback = param.control.has_poroelastic_pressure_feedback;
    double storage, lag;
    remap_coefficients(param,var,e,storage,lag);
    const int n = (*var.connectivity)[e][j];
    value[0] = storage*(*var.ppressure)[n];
    if (feedback) value[0] += alpha*trace((*var.strain)[e]);
    if (param.control.has_thermal_pressurization) {
        const double solid = var.mat->alpha_solid(e);
        double beta = pore_weighted(phi,alpha,var.mat->alpha_fluid(e),solid);
        if (!feedback && !hydraulic_strain_form(param))
            beta -= alpha*hydraulic_mean_stress_compliance(*var.mat,e)*var.mat->bulkm(e)*solid;
        value[0] -= beta*((*var.temperature)[n]-param.bc.surface_temperature);
    }
    value[1] = -lag*(*var.dppressure)[n];
    value[2] = var.matrix_pressure ? var.mat->matrix_capacity(e)*(*var.matrix_pressure)[e] : 0;
}

} // namespace

// Pressure transport commits once, on the accepted geometry. The explicit
// force uses that pressure in this interval; a moving-mesh predictor must not
// add a second momentum kick or use the next GVS interval/mass for this one.
void HydraulicStep::correct_mechanics(const Param &param, Variables &var) const
{
    #pragma acc wait
    var.dt = interval;
    if (param.control.has_moving_mesh) {
        std::copy(initial_velocity.begin(),initial_velocity.end(),var.vel->begin());
        *var.mass = initial_inertia;
    }
#ifndef ACC
    #pragma omp parallel for default(none) shared(var, param)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector
    for (int e = 0; e < var.nelem; ++e) {
        const auto conn = (*var.connectivity)[e];
        double increment = 0;
        for (int j = 0; j < NODES_PER_ELEM; ++j)
            increment += (*var.dppressure)[conn[j]]/NODES_PER_ELEM;
        increment *= var.mat->alpha_biot(e);
        for (int d = 0; d < NDIMS; ++d) (*var.stress)[e][d] += increment;
        if (NDIMS == 2 && var.mat->is_plane_strain) (*var.stressyy)[e] += increment;
        if (param.control.is_using_mixed_stress) (*var.dpressure)[e] += NDIMS*increment;
    }
    var.fluid_pressure_consumed = true;
}

ContentRemap::ContentRemap(const Param &param, const Variables &var) :
    linear(param.control.has_hydraulic_diffusion && !param.mat.porosity_model &&
           !var.mat->has_retention &&
           (!param.control.has_poroelastic_pressure_feedback ||
            (param.mat.rheol_type & MatProps::rh_elastic))),
    nonlinear(param.control.hydraulic_content_model == 1),
    active(linear || nonlinear || var.pending_heat)
{
    if (!active) return;
    #pragma acc wait
    if (var.pending_heat) {
        heat_density.resize(var.nelem);
        for (int e = 0; e < var.nelem; ++e) {
            heat_density[e] = (*var.pending_heat)[e]/(*var.volume)[e];
            heat_total += (*var.pending_heat)[e];
        }
    }
    if (nonlinear) {
        const auto volume = lumped_volume(var);
        inventory.resize(var.nelem);
        for (int e = 0; e < var.nelem; ++e) {
            double p, t; fluid_means(var,e,p,t);
            const double ratio = effective_fluid_volume(param,var,e)/(*var.volume)[e];
            for (int j = 0; j < NODES_PER_ELEM; ++j) {
                const int n = (*var.connectivity)[e][j];
                const auto state = remap_fluid_state(var,e,n,p,t);
                const double mass = ratio*state.density*state.porosity;
                const double solid = (*var.solid_mass)[n]/volume[n], capacity = (*var.solid_capacity)[n]/volume[n];
                const double values[6] = {mass,(*var.fluid_split)[n]/volume[n],solid,capacity,
                    (capacity+var.mat->cp_fluid(e)*mass)*(*var.temperature)[n],
                    (*var.temperature)[n]-(*var.pending_temperature)[n]};
                for (int k = 0; k < 6; ++k) {
                    inventory[e][k*NODES_PER_ELEM+j] = values[k];
                    physical_total[k] += (*var.volume)[e]*values[k]/NODES_PER_ELEM;
                }
            }
        }
    }
    if (!linear) return;
    old_content.resize(var.nelem);
    for (int e = 0; e < var.nelem; ++e)
        for (int j = 0; j < NODES_PER_ELEM; ++j) {
            double value[3];
            remap_content(param,var,e,j,value);
            for (int k = 0; k < 3; ++k) {
                old_content[e][k*NODES_PER_ELEM+j] = value[k];
                total[k] += (*var.volume)[e]*value[k]/NODES_PER_ELEM;
                magnitude[k] += (*var.volume)[e]*std::fabs(value[k])/NODES_PER_ELEM;
            }
        }
}

void ContentRemap::transfer(const Variables &var, const array_t &old_coordinate,
                         const conn_t &old_connectivity)
{
    if (!active) return;
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
    if (nonlinear) projected.resize(var.nnode,0.0);
    if (!heat_density.empty()) mapped_heat.assign(var.nelem,0.0);
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
            if (!heat_density.empty()) mapped_heat[e] += overlap.volume*heat_density[old];
            for (int i = 0; i < NODES_PER_ELEM; ++i)
                for (int j = 0; j < NODES_PER_ELEM; ++j) {
                    outside[e][i] -= overlap.product[i][j];
                    if (linear)
                        for (int k = 0; k < 3; ++k)
                            mapped[k] += overlap.product[i][j]*old_content[old][k*NODES_PER_ELEM+j];
                    if (nonlinear)
                        for (int k = 0; k < 6; ++k) {
                            const double value = overlap.product[i][j]*inventory[old][k*NODES_PER_ELEM+j];
                            projected[(*var.connectivity)[e][i]][k] += value;
                            physical_mapped[k] += value;
                        }
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
    if (!active) return;
    #pragma acc wait
    if (!heat_density.empty()) {
        *var.pending_heat = mapped_heat;
        double retained = 0;
        for (double energy : mapped_heat) retained += energy;
        std::printf("    Pending heat remap: retained/removed [J] = %.17g %.17g\n",retained,heat_total-retained);
    }
    if (nonlinear) { restore_inventory(param,var); return; }
    if (!linear) return;
    double current[3] = {}, current_magnitude[3] = {}, imported[3] = {}, capacity[3] = {};
    for (int e = 0; e < var.nelem; ++e) {
        const double volume = (*var.volume)[e];
        double storage, lag;
        remap_coefficients(param,var,e,storage,lag);
        for (int j = 0; j < NODES_PER_ELEM; ++j) {
            double value[3]; remap_content(param,var,e,j,value);
            for (int k = 0; k < 3; ++k) {
                current[k] += volume*value[k]/NODES_PER_ELEM;
                current_magnitude[k] += volume*std::fabs(value[k])/NODES_PER_ELEM;
                imported[k] += outside[e][j]*value[k];
            }
            if (pressure_boundary((*var.connectivity)[e][j],var) < 0) {
                capacity[0] += volume*storage/NODES_PER_ELEM;
                capacity[1] += volume*lag/NODES_PER_ELEM;
            }
        }
        if (var.matrix_pressure) capacity[2] += volume*var.mat->matrix_capacity(e);
    }
    double correction[3] = {}, reaction[3] = {};
    for (int k = 0; k < 3; ++k) {
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
        if (var.matrix_pressure && var.mat->matrix_capacity(e) > 0)
            (*var.matrix_pressure)[e] += correction[2];
        int free = 0;
        for (int j = 0; j < NODES_PER_ELEM; ++j)
            free += pressure_boundary((*var.connectivity)[e][j],var) < 0;
        const double mechanical = correction[0]+(hydraulic_strain_form(param) ? 0 : correction[1]);
        const double shift = -var.mat->alpha_biot(e)*mechanical*free/NODES_PER_ELEM;
        for (int d = 0; d < NDIMS; ++d) (*var.stress)[e][d] += shift;
        if (var.mat->is_plane_strain) (*var.stressyy)[e] += shift;
    }
    double remaining[3] = {};
    for (int e = 0; e < var.nelem; ++e)
        for (int j = 0; j < NODES_PER_ELEM; ++j) {
            double value[3]; remap_content(param,var,e,j,value);
            for (int k = 0; k < 3; ++k) remaining[k] += (*var.volume)[e]*value[k]/NODES_PER_ELEM;
        }
    for (int k = 0; k < 3; ++k)
        if (std::fabs(remaining[k]-mapped[k]-imported[k]-reaction[k]) >
            1e-10*std::max(magnitude[k]+current_magnitude[k],1e-30))
            die(EXIT_INTERNAL_ASSERT,"Fluid remapping fails its inventory constraint.");
    std::printf("    Biot content remap: removed/added/Dirichlet [m^%d] = %.17g %.17g %.17g; "
                "dp [Pa] = %.17g; buffer correction [m^%d] = %.17g\n",
                NDIMS,total[0]-mapped[0],imported[0],reaction[0],correction[0],NDIMS,
                capacity[1]*(-correction[1]));
}

void ContentRemap::restore_inventory(const Param &param, Variables &var) const
{
    // Only the global free-boundary inventories are constraints. Enforcing a
    // locally mapped mass at every node can create large interpolation spikes.
    Array2D<double,6> target(projected);
    const auto volume = lumped_volume(var);
    const double cf = var.mat->cp_fluid(0);
    const bool thermal = has_temperature_transport(param);
    double imported[6] = {};
    for (int e = 0; e < var.nelem; ++e) {
        double p, t; fluid_means(var,e,p,t);
        const double ratio = effective_fluid_volume(param,var,e)/(*var.volume)[e];
        for (int j = 0; j < NODES_PER_ELEM; ++j) {
            const int n = (*var.connectivity)[e][j];
            const auto state = remap_fluid_state(var,e,n,p,t);
            const double fluid = ratio*state.density*state.porosity;
            const double solid = ratio*(1-state.porosity)*var.mat->rho(e);
            const double capacity = solid*var.mat->cp(e);
            const double values[6] = {fluid,0,solid,capacity,(capacity+cf*fluid)*(*var.temperature)[n],
                                      (*var.temperature)[n]};
            for (int k = 0; k < 6; ++k) {
                const double value = outside[e][j]*values[k];
                target[n][k] += value;
                imported[k] += value;
            }
        }
    }
    uint_vec free_p(var.nnode), free_t(var.nnode);
    double mass_scale = 0, heat_scale = 0;
    int pressures = 0, temperatures = 0;
    for (int n = 0; n < var.nnode; ++n) {
        double fixed;
        free_p[n] = pressure_boundary(n,var) < 0;
        free_t[n] = thermal && !prescribed_temperature(param,var,n,fixed);
        pressures += free_p[n]; temperatures += free_t[n];
        if (!(target[n][0] > 0 && target[n][2] > 0 && target[n][3] >= 0 &&
              target[n][4] > 0 && volume[n] > 0))
            die(EXIT_RUNTIME_NAN,"Fluid remap has nonpositive physical inventory.");
        (*var.solid_mass)[n] = target[n][2];
        (*var.solid_capacity)[n] = target[n][3];
        mass_scale += target[n][0]; heat_scale += target[n][4];
    }
    const double mass_tolerance = 512*std::numeric_limits<double>::epsilon()*mass_scale;
    const double heat_tolerance = 512*std::numeric_limits<double>::epsilon()*heat_scale;
    auto evaluate = [&](double dp, double dt, double residual[2], double jacobian[2][2]) {
        residual[0] = residual[1] = 0;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j) jacobian[i][j] = 0;
        for (int n = 0; n < var.nnode; ++n) {
            const double p = (*var.ppressure)[n]+free_p[n]*dp, t = (*var.temperature)[n]+free_t[n]*dt;
            double mp, mt;
            const double mass = current_fluid_content(param,var,n,p,t,mp,mt);
            if (!(mass > 0 && std::isfinite(mass) && t > 0 && std::isfinite(t))) return false;
            if (free_p[n]) {
                residual[0] += mass-target[n][0];
                jacobian[0][0] += mp;
                if (free_t[n]) jacobian[0][1] += mt;
            }
            if (free_t[n]) {
                // A pressure reservoir carries cf*T*dM. Using its target mass
                // here incorporates that natural boundary enthalpy exactly.
                const double capacity = target[n][3]+cf*(free_p[n] ? mass : target[n][0]);
                residual[1] += capacity*t-target[n][4];
                jacobian[1][1] += capacity;
                if (free_p[n]) {
                    jacobian[1][0] += cf*t*mp;
                    jacobian[1][1] += cf*t*mt;
                }
            }
        }
        if (!pressures) jacobian[0][0] = 1;
        if (!temperatures) jacobian[1][1] = 1;
        return std::isfinite(residual[0]) && std::isfinite(residual[1]);
    };
    auto norm = [&](const double r[2]) {
        return std::max(std::fabs(r[0])/mass_tolerance,std::fabs(r[1])/heat_tolerance);
    };
    double dp = 0, dt = 0;
    bool converged = false;
    for (int iteration = 0; iteration < 100; ++iteration) {
        double residual[2], jacobian[2][2];
        if (!evaluate(dp,dt,residual,jacobian)) break;
        if (norm(residual) <= 1) { converged = true; break; }
        const double determinant = jacobian[0][0]*jacobian[1][1]-jacobian[0][1]*jacobian[1][0];
        if (!(determinant > 0 && std::isfinite(determinant))) break;
        double pstep = (jacobian[0][1]*residual[1]-jacobian[1][1]*residual[0])/determinant;
        double tstep = (jacobian[1][0]*residual[0]-jacobian[0][0]*residual[1])/determinant;
        bool accepted = false;
        for (int attempt = 0; attempt < 50; ++attempt) {
            double next[2], tangent[2][2];
            if (evaluate(dp+pstep,dt+tstep,next,tangent) && norm(next) < norm(residual)) {
                dp += pstep; dt += tstep; accepted = true; break;
            }
            pstep *= .5; tstep *= .5;
        }
        if (!accepted) break;
    }
    if (!converged) die(EXIT_RUNTIME_NAN,"Fluid remap cannot reconstruct positive mass and sensible energy.");
    double_vec mass(var.nnode);
    double boundary = 0, carried = 0, reaction = 0, numerical_reset = 0;
    double remaining_mass = 0, remaining_heat = 0;
    for (int n = 0; n < var.nnode; ++n) {
        double mp, mt;
        const double t = (*var.temperature)[n]+free_t[n]*dt;
        mass[n] = current_fluid_content(param,var,n,(*var.ppressure)[n]+free_p[n]*dp,t,mp,mt);
        const double water_reaction = free_p[n] ? 0 : mass[n]-target[n][0];
        const double heat_reaction = cf*t*water_reaction;
        const double energy = (target[n][3]+cf*mass[n])*t;
        boundary += water_reaction; carried += heat_reaction;
        if (!free_t[n]) reaction += energy-target[n][4]-heat_reaction;
        if (!free_p[n]) numerical_reset -= target[n][1];
        remaining_mass += mass[n]; remaining_heat += energy;
    }
    if (std::fabs(remaining_mass-physical_mapped[0]-imported[0]-boundary) > 2*mass_tolerance ||
        std::fabs(remaining_heat-physical_mapped[4]-imported[4]-carried-reaction) > 2*heat_tolerance)
        die(EXIT_INTERNAL_ASSERT,"Fluid remap fails its physical inventory constraints.");
    for (int n = 0; n < var.nnode; ++n) {
        (*var.ppressure)[n] += free_p[n]*dp;
        (*var.dppressure)[n] -= free_p[n]*dp;
        (*var.temperature)[n] += free_t[n]*dt;
        (*var.pending_temperature)[n] = (*var.temperature)[n]-target[n][5]/volume[n];
        (*var.fluid_split)[n] = free_p[n] ? target[n][1] : 0;
    }
    for (int e = 0; e < var.nelem; ++e) {
        double pressure_change = 0, temperature_change = 0;
        for (int j = 0; j < NODES_PER_ELEM; ++j) {
            const int n = (*var.connectivity)[e][j];
            pressure_change += free_p[n]*dp/NODES_PER_ELEM;
            temperature_change += free_t[n]*dt/NODES_PER_ELEM;
        }
        (*var.pore_strain)[e] += var.mat->beta_mineral(e)*pressure_change-
                                   var.mat->alpha_solid(e)*temperature_change;
    }
    double actual = 0;
    for (int n = 0; n < var.nnode; ++n) actual += nodal_fluid_mass(param,var,n);
    if (std::fabs(actual-remaining_mass) > 2*mass_tolerance)
        die(EXIT_INTERNAL_ASSERT,"Fluid remap reconstruction changes its anchored mass.");
    std::printf("    Fluid mass remap: removed/added/boundary [kg] = %.17g %.17g %.17g; numerical reset [kg] = %.17g\n",
                physical_total[0]-physical_mapped[0],imported[0],boundary,numerical_reset);
    std::printf("    Sensible heat remap: removed/added/carried/reaction [J] = %.17g %.17g %.17g %.17g; dp/dT = %.17g %.17g\n",
                physical_total[4]-physical_mapped[4],imported[4],carried,reaction,dp,dt);
    std::printf("    Solid remap: removed/added mass [kg] = %.17g %.17g; capacity [J/K] = %.17g %.17g\n",
                physical_total[2]-physical_mapped[2],imported[2],physical_total[3]-physical_mapped[3],imported[3]);
    std::printf("    Numerical buffer remap: removed/added/reset [kg] = %.17g %.17g %.17g\n",
                physical_total[1]-physical_mapped[1],imported[1],numerical_reset);
}
