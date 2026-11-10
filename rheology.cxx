#include <algorithm>
#include <cmath>
#include <iostream>

#include "3x3-C/dsyevc3.h"
#include "3x3-C/dsyevh3.h"

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "hydraulics.hpp"
#include "rheology.hpp"
#include "utils.hpp"

#ifdef THREED
// Half-width of the undecided band around the yield surface, relative to
// (|p0| + anphi*|p2| + |amc|). Covers dsyevc3's 6.6e-4 error relative to max|lambda|:
// anphi >= 1, so band >= 1e-2*max|lambda|; worst measured error/band 4e-2 over 2.6e7 tensors.
// Only the 3D yield test pre-filters, so 2D has no use for it.
const double YIELD_PREFILTER_MARGIN = 1e-2;
#endif

// Order p ascending, applying the same permutation to the columns of v when given.
// One ordering rule for both entry points; v == nullptr is the values-only case.
#pragma acc routine seq
static void sort_principal3(double p[3], double v[3][3])
{
    // Three compare-swaps: the (0,1),(1,2),(0,1) network sorts three elements.
    const int pairs[3][2] = {{0,1},{1,2},{0,1}};
    #pragma acc loop seq
    for (int k=0; k<3; ++k) {
        const int i = pairs[k][0], j = pairs[k][1];
        if (p[i] > p[j]) {
            const double tmp = p[i];
            p[i] = p[j];
            p[j] = tmp;
            if (v) {
                #pragma acc loop seq
                for (int r=0; r<3; ++r) {
                    const double b = v[r][i];
                    v[r][i] = v[r][j];
                    v[r][j] = b;
                }
            }
        }
    }
}

// Unflatten {XX, YY, ZZ, XY, XZ, YZ} into the upper triangle of a 3x3.
#pragma acc routine seq
template <typename T>
static void unflatten_stress3(T s, double a[3][3])
{
    a[0][0] = s[0];
    a[1][1] = s[1];
    a[2][2] = s[2];
    a[0][1] = s[3];
    a[0][2] = s[4];
    a[1][2] = s[5];
}

/* Eigenvalues ONLY, ascending, for the yield test in elasto_plastic(). dsyevc3 (Cardano)
 * is 4.8x cheaper (25.1 vs 121.1 ns/call) but inaccurate near degeneracy, so a decision
 * near the surface is undecided -- see YIELD_PREFILTER_MARGIN. */
#pragma acc routine seq
template <typename T>
static void principal_values3(T s, double p[3])
{
    double a[3][3];
    unflatten_stress3(s, a);
    dsyevc3(a, p);
    sort_principal3(p, nullptr);
}

/* Eigenvalues AND eigenvectors, p ascending, v holding eigenvectors as columns.
 * Vendored dsyevh3, reached only on the branch that needs vectors.
 */
#pragma acc routine seq
template <typename T>
static void principal_stresses3(T s, double p[3], double v[3][3])
{
    double a[3][3];
    unflatten_stress3(s, a);
    dsyevh3(a, v, p);
    sort_principal3(p, v);
}

#pragma acc routine seq
template <typename T>
static void principal_stresses2(T s, double p[2],
                                double& cos2t, double& sin2t)
{
    /* 's' is a flattened stress vector, with the components {XX, ZZ, XZ}.
     * Returns the eigenvalues 'p', and the direction cosine of the
     * eigenvectors in the X-Z plane.
     * The eigenvalues are ordered such that p[0] <= p[1].
     */

    // center and radius of Mohr circle
    double s0 = 0.5 * (s[0] + s[1]);
    double rad = second_invariant(s);

    // principal stresses in the X-Z plane
    p[0] = s0 - rad;
    p[1] = s0 + rad;

    {
        // direction cosine and sine of 2*theta
        const double eps = 1e-15;
        double a = 0.5 * (s[0] - s[1]);
        double b = - rad; // always negative
        if (b < -eps) {
            cos2t = a / b;
            sin2t = s[2] / b;
        }
        else {
            cos2t = 1;
            sin2t = 0;
        }
    }
}

#pragma acc routine seq
static void compute_slip_rate2(ConstTensorAccessor s, double& vx, double& vz, double& slip_rate)
{
    /* Computes the slip magnitude by projecting velocity (vx, vz) onto the maximum shear direction.
     *
     * Inputs:
     *  - s: Flattened stress tensor {σ_xx, σ_zz, σ_xz}
     *  - vx, vz: Velocity components
     * 
     * Output:
     *  - slip: The scalar magnitude of the slip vector (computed by projection)
     */

    // Compute center and radius of Mohr circle
    double s0 = 0.5 * (s[0] + s[1]);
    double rad = second_invariant(s);
    double cos2t, sin2t;

    {
        // Compute direction cosine and sine of 2*theta
        const double eps = 1e-15;
        double a = 0.5 * (s[0] - s[1]);
        double b = -rad; // Always negative
        if (b < -eps) {
            cos2t = a / b;
            sin2t = s[2] / b;
        }
        else {
            cos2t = 1;
            sin2t = 0;
        }
    }

    // Compute shear direction (rotated 45° from principal stress direction)
    double theta_shear = 0.5 * atan2(sin2t, cos2t);
    double shear_dir[2];
    shear_dir[0] = cos(theta_shear + M_PI_4); // Shear direction X-component
    shear_dir[1] = sin(theta_shear + M_PI_4); // Shear direction Z-component

    // slip magnitude
    slip_rate = std::abs(vx * shear_dir[0] + vz * shear_dir[1]);
}

#pragma acc routine seq
static void compute_slip_rate3(ConstTensorAccessor s, double vx, double vy, double vz, double& slip_rate)
{
    /* Computes the slip magnitude by projecting velocity (vx, vy, vz) onto the maximum shear direction in 3D.
     *
     * Inputs:
     *  - s: Flattened 3D stress tensor {σ_xx, σ_yy, σ_zz, σ_xy, σ_xz, σ_yz}
     *  - vx, vy, vz: Velocity components
     * 
     * Output:
     *  - slip_rate: Scalar magnitude of the slip vector (computed by projection)
     */

    double p[3], v[3][3];

    // Compute principal stresses and eigenvectors
    principal_stresses3(s, p, v);

    // Compute maximum shear stress plane
    double tau_max_1 = 0.5 * std::abs(p[2] - p[1]);  // Shear between σ3 and σ2
    double tau_max_2 = 0.5 * std::abs(p[2] - p[0]);  // Shear between σ3 and σ1
    double tau_max_3 = 0.5 * std::abs(p[1] - p[0]);  // Shear between σ2 and σ1

    int shear_plane_index = (tau_max_2 >= tau_max_1 && tau_max_2 >= tau_max_3) ? 1 :
                            (tau_max_3 >= tau_max_1 && tau_max_3 >= tau_max_2) ? 2 : 0;

    // Fault normal (perpendicular to max shear stress plane)
    double fault_normal[3] = {v[0][shear_plane_index], v[1][shear_plane_index], v[2][shear_plane_index]};

    // Two shear directions (perpendicular to fault normal)
    double shear_dir_1[3] = {v[0][(shear_plane_index + 1) % 3], 
                             v[1][(shear_plane_index + 1) % 3], 
                             v[2][(shear_plane_index + 1) % 3]};

    double shear_dir_2[3] = {v[0][(shear_plane_index + 2) % 3], 
                             v[1][(shear_plane_index + 2) % 3], 
                             v[2][(shear_plane_index + 2) % 3]};

    // Project velocity onto shear directions
    double slip_magnitude_1 = vx * shear_dir_1[0] + vy * shear_dir_1[1] + vz * shear_dir_1[2];
    double slip_magnitude_2 = vx * shear_dir_2[0] + vy * shear_dir_2[1] + vz * shear_dir_2[2];

    // Compute final slip magnitude
    slip_rate = std::sqrt(slip_magnitude_1 * slip_magnitude_1 + slip_magnitude_2 * slip_magnitude_2);
}

// Element width w for V = 2 w eps_II: the minimum altitude, matching compute_dt().
#pragma acc routine seq
static double elem_band_width(ConstArrayIndirectAccessor coord, double volume)
{
#ifdef THREED
    const int facet[4][3] = {{0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}};
    double max_area = 0.0;
    #pragma acc loop seq
    for (int i = 0; i < 4; ++i) {
        ConstArrayAccessor p = coord[facet[i][0]];
        ConstArrayAccessor q = coord[facet[i][1]];
        ConstArrayAccessor r = coord[facet[i][2]];
        const double u0 = q[0] - p[0];
        const double u1 = q[1] - p[1];
        const double u2 = q[2] - p[2];
        const double v0 = r[0] - p[0];
        const double v1 = r[1] - p[1];
        const double v2 = r[2] - p[2];
        const double c0 = u1 * v2 - u2 * v1;
        const double c1 = u2 * v0 - u0 * v2;
        const double c2 = u0 * v1 - u1 * v0;
        const double area = 0.5 * std::sqrt(c0*c0 + c1*c1 + c2*c2);
        max_area = std::max(max_area, area);
    }
    return (max_area > 0.0) ? 3.0 * volume / max_area : 0.0;
#else
    double max_length2 = 0.0;
    #pragma acc loop seq
    for (int i = 0; i < 3; ++i) {
        ConstArrayAccessor p = coord[i];
        ConstArrayAccessor q = coord[(i + 1) % 3];
        const double dx = p[0] - q[0];
        const double dz = p[1] - q[1];
        max_length2 = std::max(max_length2, dx*dx + dz*dz);
    }
    const double max_length = std::sqrt(max_length2);
    return (max_length > 0.0) ? 2.0 * volume / max_length : 0.0;
#endif
}

#pragma acc routine seq
static void select_slip_rate2(int projection_option, double continuum_rate,
                              ConstTensorAccessor stress,
                              double& vx, double& vz, double& slip_rate)
{
    if (projection_option == rsf_slip_rate_projection_total_strain_rate) {
        slip_rate = continuum_rate;
        return;
    }
    compute_slip_rate2(stress, vx, vz, slip_rate);
}

#pragma acc routine seq
static void select_slip_rate3(int projection_option, double continuum_rate,
                              ConstTensorAccessor stress,
                              double vx, double vy, double vz,
                              double& slip_rate)
{
    if (projection_option == rsf_slip_rate_projection_total_strain_rate) {
        slip_rate = continuum_rate;
        return;
    }
    compute_slip_rate3(stress, vx, vy, vz, slip_rate);
}

void refresh_rsf_friction(const Param& param, Variables& var,
                          double_vec& dyn_fric_coeff,
                          const double_vec& state_variable)
{
    if (!(param.mat.rheol_type & MatProps::rh_rsf)) return;

#ifndef ACC
    #pragma omp parallel for default(none) shared(param, var, dyn_fric_coeff, state_variable)
#endif
    for (int e = 0; e < var.nelem; ++e) {
        const int projection_option =
            param.control.rsf_slip_rate_projection_option;
        double vx = 0.0;
        double vy = 0.0;
#ifdef THREED
        double vz = 0.0;
#endif

        if (projection_option ==
            rsf_slip_rate_projection_maximum_shear) {
            ConstArrayIndirectAccessor v =
                var.vel->view_const((*var.connectivity)[e]);
            const double weight = 1.0 / NODES_PER_ELEM;
            for (int j = 0; j < NODES_PER_ELEM; ++j) {
                vx += v[j][0] * weight;
                vy += v[j][1] * weight;
#ifdef THREED
                vz += v[j][2] * weight;
#endif
            }
        }

        double continuum_rate = 0.0;
        if (projection_option ==
            rsf_slip_rate_projection_total_strain_rate) {
            const double width = elem_band_width(
                var.coord->view_const((*var.connectivity)[e]),
                (*var.volume)[e]);
            continuum_rate =
                2.0 * width * second_invariant((*var.strain_rate)[e]);
        }

        double slip_rate = 0.0;
#ifdef THREED
        select_slip_rate3(projection_option, continuum_rate,
                          (*var.stress)[e], vx, vy, vz, slip_rate);
#else
        select_slip_rate2(projection_option, continuum_rate,
                          (*var.stress)[e], vx, vy, slip_rate);
#endif
        var.mat->rsf_friction_from_state(e, (*var.plstrain)[e], slip_rate,
                                         state_variable[e], dyn_fric_coeff[e],
                                         param.mat.state_var_model);
    }
}

#pragma acc routine seq
template <typename T>
static void elastic(double bulkm, double shearm, const double* de, T s, double& syy)
{
    /* increment the stress s according to the incremental strain de */
    double lambda = bulkm - 2. /3 * shearm;
    double dev = trace(de);
    double volumetric_normal_increment = lambda * dev;

    for (int i=0; i<NDIMS; ++i)
        s[i] += 2 * shearm * de[i] + volumetric_normal_increment;
    for (int i=NDIMS; i<NSTR; ++i)
        s[i] += 2 * shearm * de[i];

    // Plane strain has zero out-of-plane strain.
    syy += volumetric_normal_increment;
}

#pragma acc routine seq
template <typename T>
static void elastic(double bulkm, double shearm, const double* de, T s)
{
    double unused_syy = 0;
    elastic(bulkm, shearm, de, s, unused_syy);
}

#pragma acc routine seq
template <typename T>
static void elastic_effective(double bulkm, double shearm, const double* de, T s,
                              double& dpp, double& syy)
{
    /* increment the stress s according to the incremental strain de */
    double lambda = bulkm - 2. /3 * shearm;
    double dev = trace(de);
    double volumetric_normal_increment = lambda * dev;

    for (int i=0; i<NDIMS; ++i)
        s[i] += 2 * shearm * de[i] + volumetric_normal_increment + dpp;

    for (int i=NDIMS; i<NSTR; ++i)
        s[i] += 2 * shearm * de[i];

    syy += volumetric_normal_increment + dpp;
}

#pragma acc routine seq
template <typename T>
static void elastic_effective(double bulkm, double shearm, const double* de, T s,
                              double& dpp)
{
    double unused_syy = 0;
    elastic_effective(bulkm, shearm, de, s, dpp, unused_syy);
}

#pragma acc routine seq
template <typename T>
static void maxwell(double bulkm, double shearm, double viscosity, double dt,
                    double dv, const double* de, T s, double* syy = nullptr, double dpp = 0)
{
    // non-dimensional parameter: dt/ relaxation time
    double tmp = 0.5 * dt * shearm / viscosity;
    double f1 = 1 - tmp;
    double f2 = 1 / (1  + tmp);

    const int normal_components = syy ? 3 : NDIMS;
    double dev = trace(de) / normal_components;
    double s0 = syy ? (trace(s) + *syy)/3 : trace(s)/NDIMS;

    if (syy)
        *syy = ((*syy - s0)*f1 - 2*shearm*dev)*f2 + s0 + bulkm*dv;

    // convert back to total stress
    for (int i=0; i<NDIMS; ++i)
        s[i] = ((s[i] - s0) * f1 + 2 * shearm * (de[i] - dev)) * f2 + s0 + bulkm * dv;
    for (int i=NDIMS; i<NSTR; ++i)
        s[i] = (s[i] * f1 + 2 * shearm * de[i]) * f2;
    if (dpp != 0) {
        for (int i=0; i<NDIMS; ++i) s[i] += dpp;
        if (syy) *syy += dpp;
    }
}


#pragma acc routine seq
static void viscous(double bulkm, double viscosity, double total_dv,
                    ConstTensorAccessor edot, TensorAccessor s, double* syy = nullptr, double pressure = 0)
{
    /* Viscous Model + incompressibility enforced by bulk modulus */

    double dev = trace(edot) / (syy ? 3 : NDIMS);

    for (int i=0; i<NDIMS; ++i)
        s[i] = 2 * viscosity * (edot[i] - dev) + bulkm * total_dv;
    for (int i=NDIMS; i<NSTR; ++i)
        s[i] = 2 * viscosity * edot[i];
    if (syy) *syy = -2*viscosity*dev + bulkm*total_dv;
    if (pressure != 0) {
        for (int i=0; i<NDIMS; ++i) s[i] -= pressure;
        if (syy) *syy -= pressure;
    }
}

// Recover only the irreversible increment with the inverse of the Hooke law
// actually used by this branch. dpp is a stress shift, not elastic strain.
// The start tensor is taken after reversible thermal expansion. Work uses
// effective stress, including the mean Bishop pressure over this interval.
#pragma acc routine seq
static double plastic_heat(double bulk, double shear, const double *start,
                           ConstTensorAccessor final, const double *de, double dpp,
                           bool has_syy, double syy_start, double syy_final,
                           double work_pressure)
{
    const double lambda = bulk - 2.0/3.0 * shear;
    const int trace_dim = has_syy ? 3 : NDIMS;
    const double ratio = lambda / (trace_dim * lambda + 2 * shear);
    const double t0 = trace(start) + (has_syy ? syy_start : 0);
    const double t1 = trace(final) - NDIMS*dpp + (has_syy ? syy_final-dpp : 0);
    double heat = 0;
    for (int d = 0; d < NSTR; ++d) {
        const bool normal = d < NDIMS;
        const double change = final[d] - start[d] - (normal ? dpp : 0);
        const double elastic = (change - (normal ? ratio*(t1-t0) : 0)) / (2*shear);
        const double work = 0.5*(start[d]+final[d]) + (normal ? work_pressure : 0);
        heat += (normal ? 1 : 2) * work * (de[d]-elastic);
    }
    if (has_syy) {
        const double elastic = (syy_final-dpp-syy_start-ratio*(t1-t0)) / (2*shear);
        heat -= (0.5*(syy_start+syy_final)+work_pressure) * elastic;
    }
    // Individual normal components are not dissipation: pressure contributions
    // cancel only in the sum. Never clamp each contribution separately.
    return heat < 0 ? 0 : heat;
}

// Maxwell uses a centered viscous strain increment. Memoryless flow uses
// its instantaneous deviatoric rate instead; neither dissipates bulk work.
#pragma acc routine seq
static double viscous_heat(double viscosity, double dt, ConstTensorAccessor rate,
                           const double *start, ConstTensorAccessor final,
                           bool maxwell_branch, bool has_syy,
                           double syy_start, double syy_final)
{
    const int trace_dim = has_syy ? 3 : NDIMS;
    const double mean = maxwell_branch ?
        (trace(start)+trace(final)+(has_syy ? syy_start+syy_final : 0))/(2*trace_dim) :
        trace(rate)/trace_dim;
    double squared = 0;
    for (int d = 0; d < NSTR; ++d) {
        const double value = (maxwell_branch ? 0.5*(start[d]+final[d]) : rate[d]) -
                             (d < NDIMS ? mean : 0);
        squared += (d < NDIMS ? 1 : 2) * value*value;
    }
    if (has_syy) {
        const double value = (maxwell_branch ? 0.5*(syy_start+syy_final) : 0) - mean;
        squared += value*value;
    }
    return dt * squared * (maxwell_branch ? 1/(2*viscosity) : 2*viscosity);
}

// Three solid principal stresses are used in both 3D and plane strain.
// A single-face return must not cross their ordering: try the two adjacent
// shear faces together, keeping the original trial axes and hardening tangent.
// Apex and combined tensile/shear returns retain their separate treatment.
#pragma acc routine seq
static bool shear_edge_return(const double trial[3], double p[3], double plastic[3],
                              double a1, double a2, double anphi, double anpsi,
                              double amc, double hardening, double ten_max)
{
    if (p[0] <= p[1] && p[1] <= p[2]) return false;
    const double scale = std::max(std::fabs(amc),
        std::max(std::fabs(trial[0]), std::max(std::fabs(trial[1]), std::fabs(trial[2]))));
    const double tolerance = 1e-12 * scale;
    for (int edge = 0; edge < 2; ++edge) {
        const int low[2] = {0, edge == 0 ? 1 : 0};
        const int high[2] = {2, edge == 0 ? 2 : 1};
        double flow[2][3], response[2][3];
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 3; ++k) {
                flow[j][k] = (k == high[j] ? anpsi : 0) - (k == low[j] ? 1 : 0);
                response[j][k] = a2 * (anpsi - 1) + (a1 - a2) * flow[j][k];
            }
        double matrix[2][2], rhs[2];
        for (int i = 0; i < 2; ++i) {
            rhs[i] = anphi * trial[high[i]] - trial[low[i]] - amc;
            for (int j = 0; j < 2; ++j)
                matrix[i][j] = anphi * response[j][high[i]] - response[j][low[i]] + hardening;
        }
        const double determinant = matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0];
        if (!(determinant > 0) || !std::isfinite(determinant)) continue;
        const double multiplier[2] = {
            (rhs[0]*matrix[1][1] - rhs[1]*matrix[0][1]) / determinant,
            (rhs[1]*matrix[0][0] - rhs[0]*matrix[1][0]) / determinant};
        if (!std::isfinite(multiplier[0]) || !std::isfinite(multiplier[1]) ||
            multiplier[0] < 0 || multiplier[1] < 0) continue;
        double candidate[3];
        for (int k = 0; k < 3; ++k)
            candidate[k] = trial[k] - multiplier[0]*response[0][k] - multiplier[1]*response[1][k];
        if (!std::isfinite(candidate[0]) || !std::isfinite(candidate[1]) || !std::isfinite(candidate[2]) ||
            candidate[0] > candidate[1] + tolerance || candidate[1] > candidate[2] + tolerance ||
            candidate[2] > ten_max + tolerance) continue;
        const double strength = amc + hardening * (multiplier[0] + multiplier[1]);
        if (candidate[0] - anphi*candidate[2] + strength < -tolerance) continue;
        // Enforce the repeated eigenvalue without letting solve roundoff pick
        // a preferred direction in its eigenspace.
        const int repeated = edge == 0 ? 0 : 1;
        candidate[repeated] = candidate[repeated+1] =
            0.5 * (candidate[repeated] + candidate[repeated+1]);
        for (int k = 0; k < 3; ++k) {
            p[k] = candidate[k];
            plastic[k] = multiplier[0]*flow[0][k] + multiplier[1]*flow[1][k];
        }
        return true;
    }
    return false;
}


// Return whether the shear correction satisfies ordering and the tensile bound.
#pragma acc routine seq
static bool shear_return(const double trial[3], double p[3], double &depls,
                         double a1, double a2, double anphi, double anpsi,
                         double amc, double hardening, double ten_max,
                         bool plane_strain, bool allow_edge = true)
{
    const double fs = trial[0] - trial[2]*anphi + amc;
    const double alam = fs / (a1 - a2*anpsi + a1*anphi*anpsi - a2*anphi + hardening);
    p[0] = trial[0] - alam * (a1 - a2*anpsi);
    p[1] = trial[1] - alam * (a2 - a2*anpsi);
    p[2] = trial[2] - alam * (a2 - a1*anpsi);
    double plastic[3] = {alam, 0, -alam*anpsi};
    const bool edge = allow_edge &&
        shear_edge_return(trial, p, plastic, a1, a2, anphi, anpsi, amc, hardening, ten_max);
    if (plane_strain) {
        depls = edge ? 0.5 * (std::max(plastic[0], std::max(plastic[1], plastic[2])) -
                             std::min(plastic[0], std::min(plastic[1], plastic[2]))) :
                      0.5 * std::fabs(alam + alam*anpsi);
    }
    else if (edge) {
        const double mean = (plastic[0] + plastic[1] + plastic[2]) / 3;
        double squared = mean*mean;
        for (int k = 0; k < 3; ++k) squared += (plastic[k] - mean)*(plastic[k] - mean);
        depls = std::sqrt(0.5 * squared);
    }
    else
        depls = std::fabs(alam) * std::sqrt((7 + 4*anpsi + 7*anpsi*anpsi) / 18);
    const double scale = std::max(std::fabs(amc),
        std::max(std::fabs(trial[0]), std::max(std::fabs(trial[1]), std::fabs(trial[2]))));
    const double tolerance = 1e-12 * scale;
    const double total_multiplier = (std::fabs(plastic[0]) + std::fabs(plastic[1]) +
                                     std::fabs(plastic[2])) / (1 + anpsi);
    return (edge || alam <= 0) && std::isfinite(depls) &&
           std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]) &&
           p[0] <= p[1] + tolerance && p[1] <= p[2] + tolerance && p[2] <= ten_max + tolerance &&
           p[0] - anphi*p[2] + amc + hardening*total_multiplier >= -tolerance;
}


#pragma acc routine seq
template <typename T>
static void rotate_stress2d(const double p[3], int n1, int n2, int n3,
                                double cos2t, double sin2t, T s, double &syy)
{
    const double dc2 = (p[n1] - p[n2]) * cos2t;
    const double dss = p[n1] + p[n2];
    s[0] = 0.5 * (dss + dc2);
    s[1] = 0.5 * (dss - dc2);
    s[2] = 0.5 * (p[n1] - p[n2]) * sin2t;
    syy = p[n3];
}


#pragma acc routine seq
template <typename T>
static void elasto_plastic(double bulkm, double shearm,
                           double amc, double anphi, double anpsi,
                           double hardn, double ten_max,
                           const double* de, double& depls, T s,
                           int &failure_mode,
                           bool has_hydraulic_diffusion,
                           double &dpp)
{
    /* Elasto-plasticity (Mohr-Coulomb criterion)
     *
     * failure_mode --
     *   0: no failure
     *   1: tensile failure
     *  10: shear failure
     */

    // elastic trial stress
    if (has_hydraulic_diffusion)
    {
        elastic_effective(bulkm, shearm, de, s, dpp);
    }
    else
    {
        elastic(bulkm, shearm, de, s);
    }
    depls = 0;
    failure_mode = 0;

    //
    // transform to principal stress coordinate system
    //
    // eigenvalues (principal stresses)
    double p[NDIMS];
#ifdef THREED
    // eigenvectors
    double v[3][3];

    // Eigenvalue-only PRE-FILTER: v is needed only to rotate a CORRECTED stress back, and
    // only 0.170% of elements yield, so the eigenvectors are almost never needed.
    // One-sided -- anything within the band falls through to the accurate solver.
    {
        double pf[3];
        principal_values3(s, pf);
        const double band = YIELD_PREFILTER_MARGIN
            * (std::fabs(pf[0]) + anphi * std::fabs(pf[2]) + std::fabs(amc));
        if (pf[0] - pf[2] * anphi + amc > band && pf[2] - ten_max < -band)
            return;   // no failure, and not close enough to it to need checking
    }

    principal_stresses3(s, p, v);
#else
    // In 2D, we only construct the eigenvectors from
    // cos(2*theta) and sin(2*theta) of Mohr circle
    double cos2t, sin2t;
    principal_stresses2(s, p, cos2t, sin2t);
#endif

    // composite (shear and tensile) yield criterion
    double fs = p[0] - p[NDIMS-1] * anphi + amc;
    double ft = p[NDIMS-1] - ten_max;

    if (fs > 0 && ft < 0) {
        // no failure
        return;
    }

    // yield, shear or tensile?
    double pa = std::sqrt(1 + anphi*anphi) + anphi;
    double ps = ten_max * anphi - amc;
    double h = p[NDIMS-1] - ten_max + pa * (p[0] - ps);
    double a1 = bulkm + 4. / 3 * shearm;
    double a2 = bulkm - 2. / 3 * shearm;

    double alam;
#ifdef THREED
    const double trial[3] = {p[0], p[1], p[2]};
    double candidate[3], increment;
    // A tensile elastic predictor does not imply tensile plastic flow when
    // the shear return itself finishes below the cap.
    if (ft > 0 && fs < 0 &&
        shear_return(trial, candidate, increment, a1, a2, anphi, anpsi, amc,
                     2*std::sqrt(anphi)*hardn, ten_max, false)) {
        for (int k = 0; k < 3; ++k) p[k] = candidate[k];
        depls = increment;
        failure_mode = 10;
    }
    else
#endif
    if (h < 0) {
        // shear failure
        failure_mode = 10;

#ifdef THREED
        shear_return(trial, p, depls, a1, a2, anphi, anpsi, amc,
                     2*std::sqrt(anphi)*hardn, ten_max, false);
#else
        alam = fs / (a1 - a2*anpsi + a1*anphi*anpsi - a2*anphi + 2*std::sqrt(anphi)*hardn);
        p[0] -= alam * (a1 - a2 * anpsi);
        p[NDIMS-1] -= alam * (a2 - a1 * anpsi);

        // 2nd invariant of plastic strain
        /* // plastic strain in the principal directions
         * double depls1 = alam;
         * double depls2 = -alam * anpsi;
         * double deplsm = (depls1 + depls2) / 2;
         * depls = std::sqrt(((depls1-deplsm)*(depls1-deplsm) +
         *                    (depls2-deplsm)*(depls2-deplsm) +
         *                    deplsm*deplsm) / 2);
         */
        // the equations above can be reduce to:
        depls = std::fabs(alam) * std::sqrt((3 + 2*anpsi + 3*anpsi*anpsi) / 8);
#endif
    }
    else {
        // tensile failure
        failure_mode = 1;

        alam = ft / a1;
        p[0] -= alam * a2;
#ifdef THREED
        p[1] -= alam * a2;
#endif
        p[NDIMS-1] -= alam * a1;

        // 2nd invariant of plastic strain
#ifdef THREED
        /* double depls1 = 0;
         * double depls3 = alam;
         * double deplsm = (depls1 + depls3) / 3;
         * depls = std::sqrt(((depls1-deplsm)*(depls1-deplsm) +
         *                    (-deplsm)*(-deplsm) +
         *                    (depls3-deplsm)*(depls3-deplsm) +
         *                    deplsm*deplsm) / 2);
         */
        depls = std::fabs(alam) * std::sqrt(7. / 18);
#else
        /* double depls1 = 0;
         * double depls3 = alam;
         * double deplsm = (depls1 + depls3) / 2;
         * depls = std::sqrt(((depls1-deplsm)*(depls1-deplsm) +
         *                    (depls2-deplsm)*(depls2-deplsm) +
         *                    deplsm*deplsm) / 2);
         */
        depls = std::fabs(alam) * std::sqrt(3. / 8);
#endif
    }

    // rotate the principal stresses back to global axes
    {
#ifdef THREED
        double ss[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        for(int m=0; m<3; m++) {
            for(int n=m; n<3; n++) {
                for(int k=0; k<3; k++) {
                    ss[m][n] += v[m][k] * v[n][k] * p[k];
                }
            }
        }
        s[0] = ss[0][0];
        s[1] = ss[1][1];
        s[2] = ss[2][2];
        s[3] = ss[0][1];
        s[4] = ss[0][2];
        s[5] = ss[1][2];
#else
        double dc2 = (p[0] - p[1]) * cos2t;
        double dss = p[0] + p[1];
        s[0] = 0.5 * (dss + dc2);
        s[1] = 0.5 * (dss - dc2);
        s[2] = 0.5 * (p[0] - p[1]) * sin2t;
#endif
    }
}

#pragma acc routine seq
template <typename T>
static void elasto_plastic2d(double bulkm, double shearm,
                             double amc, double anphi, double anpsi,
                             double hardn, double ten_max,
                             const double* de, double& depls,
                             T s, double &syy,
                             int &failure_mode,
                             bool has_hydraulic_diffusion,
                             double &dpp)
{
    /* Elasto-plasticity (Mohr-Coulomb criterion) */

    /* This function is derived from geoFLAC.
     * The original code in geoFLAC assumes 2D plane strain formulation,
     * i.e. there are 3 principal stresses (PSs) and only 2 principal strains
     * (Strain_yy, Strain_xy, and Strain_yz all must be 0).
     * Here, the code is adopted to pure 2D or 3D plane strain.
     *
     * failure_mode --
     *   0: no failure
     *   1: tensile failure, all PSs exceed tensile limit
     *   2: tensile failure, 2 PSs exceed tensile limit
     *   3: tensile failure, 1 PS exceeds tensile limit
     *  10: pure shear failure
     *  11, 12, 13: tensile + shear failure
     *  20, 21, 22, 23: shear + tensile failure
     */

    depls = 0;
    failure_mode = 0;

    // elastic trial stress
    double a1 = bulkm + 4. / 3 * shearm;
    double a2 = bulkm - 2. / 3 * shearm;
    double sxx = s[0] + de[1]*a2 + de[0]*a1;
    double szz = s[1] + de[0]*a2 + de[1]*a1;
    double sxz = s[2] + de[2]*2*shearm;
    syy += (de[0] + de[1]) * a2; // Stress YY component, plane strain

    // Apply the pore pressure effect if hydraulic diffusion is enabled
    if (has_hydraulic_diffusion)
    {
        sxx += dpp;
        syy += dpp;
        szz += dpp;
    }
    //
    // transform to principal stress coordinate system
    //
    // eigenvalues (principal stresses)
    double p[3];

    // In 2D, we only construct the eigenvectors from
    // cos(2*theta) and sin(2*theta) of Mohr circle
    double cos2t, sin2t;
    int n1, n2, n3;

    {
        // center and radius of Mohr circle
        double s0 = 0.5 * (sxx + szz);
        double rad = 0.5 * std::sqrt((sxx-szz)*(sxx-szz) + 4*sxz*sxz);

        // principal stresses in the X-Z plane
        double si = s0 - rad;
        double sii = s0 + rad;

        // direction cosine and sine of 2*theta
        const double eps = 1e-15;
        if (rad > eps) {
            cos2t = 0.5 * (szz - sxx) / rad;
            sin2t = -sxz / rad;
        }
        else {
            cos2t = 1;
            sin2t = 0;
        }

        // sort p.s.
#if 1
        //
        // 3d plane strain
        //
        if (syy > sii) {
            // syy is minor p.s.
            n1 = 0;
            n2 = 1;
            n3 = 2;
            p[0] = si;
            p[1] = sii;
            p[2] = syy;
        }
        else if (syy < si) {
            // syy is major p.s.
            n1 = 1;
            n2 = 2;
            n3 = 0;
            p[0] = syy;
            p[1] = si;
            p[2] = sii;
        }
        else {
            // syy is intermediate
            n1 = 0;
            n2 = 2;
            n3 = 1;
            p[0] = si;
            p[1] = syy;
            p[2] = sii;
        }
#else
        /* XXX: This case gives unreasonable result. Don't know why... */

        //
        // pure 2d case
        //
        n1 = 0;
        n2 = 2;
        p[0] = si;
        p[1] = syy;
        p[2] = sii;
#endif
    }

    // Test the uncut trial first: shear flow may return below the tensile cap.
    if (p[2] > ten_max && p[0] - p[2]*anphi + amc < 0) {
        double candidate[3], increment;
        if (shear_return(p, candidate, increment, a1, a2, anphi, anpsi,
                         amc, hardn, ten_max, true)) {
            rotate_stress2d(candidate, n1, n2, n3, cos2t, sin2t, s, syy);
            depls = increment;
            failure_mode = 10;
            return;
        }
    }

    // Possible tensile failure scenarios
    // 1. S1 (least tensional or greatest compressional principal stress) > ten_max:
    //    all three principal stresses must be greater than ten_max.
    //    Assign ten_max to all three and don't do anything further.
    if( p[0] >= ten_max ) {
        s[0] = s[1] = syy = ten_max;
        s[2] = 0.0;
        failure_mode = 1;
        return;
    }

    // 2. S2 (intermediate principal stress) > ten_max:
    //    S2 and S3 must be greater than ten_max.
    //    Assign ten_max to these two and continue to the shear failure block.
    if( p[1] >= ten_max ) {
        p[1] = p[2] = ten_max;
        failure_mode = 2;
    }

    // 3. S3 (most tensional or least compressional principal stress) > ten_max:
    //    Only this must be greater than ten_max.
    //    Assign ten_max to S3 and continue to the shear failure block.
    else if( p[2] >= ten_max ) {
        p[2] = ten_max;
        failure_mode = 3;
    }


    // shear yield criterion
    double fs = p[0] - p[2] * anphi + amc;
    if (fs >= 0.0) {
        // Tensile failure case S2 or S3 could have happened!!
        // XXX: Need to rationalize why exit without doing anything further.
        s[0] = sxx;
        s[1] = szz;
        s[2] = sxz;
        return;
    }

    failure_mode += 10;

    // shear failure
    const double trial[3] = {p[0], p[1], p[2]};
    shear_return(trial, p, depls, a1, a2, anphi, anpsi, amc, hardn, ten_max,
                 true, failure_mode == 10);

    //***********************************
    // The following seems redundant but... this is how it goes in geoFLAC.
    //
    // Possible tensile failure scenarios
    // 1. S1 (least tensional or greatest compressional principal stress) > ten_max:
    //    all three principal stresses must be greater than ten_max.
    //    Assign ten_max to all three and don't do anything further.
    if( p[0] >= ten_max ) {
        s[0] = s[1] = syy = ten_max;
        s[2] = 0.0;
        failure_mode += 20;
        return;
    }

    // 2. S2 (intermediate principal stress) > ten_max:
    //    S2 and S3 must be greater than ten_max.
    //    Assign ten_max to these two and continue to the shear failure block.
    if( p[1] >= ten_max ) {
        p[1] = p[2] = ten_max;
        failure_mode += 20;
    }

    // 3. S3 (most tensional or least compressional principal stress) > ten_max:
    //    Only this must be greater than ten_max.
    //    Assign ten_max to S3 and continue to the shear failure block.
    else if( p[2] >= ten_max ) {
        p[2] = ten_max;
        failure_mode += 20;
    }
    //***********************************


    // rotate the principal stresses back to global axes
    rotate_stress2d(p, n1, n2, n3, cos2t, sin2t, s, syy);
}

// Stored stress is total: sigma' = sigma + q I. Translating the yield
// intercept and tensile cap evaluates skeleton strength with the same
// return algorithm, without shifting and rounding the stored tensor.
#pragma acc routine seq
static void pore_strength(double q, double anphi, double &amc, double &ten_max)
{
    if (q == 0) return;
    amc -= (anphi - 1) * q;
    ten_max -= q;
}

void update_stress(const Param& param, Variables& var, tensor_t& stress,
                   double_vec& stressyy, double_vec& dpressure, double_vec& viscosity,
                   tensor_t& strain, double_vec& plstrain,
                   double_vec& delta_plstrain, tensor_t& strain_rate,
                   double_vec& ppressure, double_vec& dppressure, array_t& vel,
                   double_vec& dyn_fric_coeff, double_vec& state_variable,
                   bool trial_reference_geometry, bool initial_equilibrium)
{
#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif

    // Loop-invariant, so the per-element gathers they gate are decided once.
    // Instantaneous initialization equilibrates the existing stress without
    // consuming a pressure increment belonging to the next physical step.
    const bool has_hydraulic_diffusion =
        param.control.has_hydraulic_diffusion && !initial_equilibrium;
    const bool shear_heating = param.control.has_shear_heating && !initial_equilibrium;
    // Only option 0 projects the centroid velocity; option 1 derives its rate
    // from the strain-rate invariant.
    const bool needs_slip_rate = (param.mat.rheol_type == MatProps::rh_ep_rsf)
                              || (param.mat.rheol_type == MatProps::rh_evp_rsf);
    const bool needs_projected_velocity =
        needs_slip_rate &&
        param.control.rsf_slip_rate_projection_option ==
            rsf_slip_rate_projection_maximum_shear;

#ifndef ACC
    #pragma omp parallel for default(none) shared(param, var, ppressure, dppressure, \
        vel, stress, stressyy, dpressure, viscosity, strain, plstrain, delta_plstrain, \
        strain_rate, dyn_fric_coeff, state_variable) \
        firstprivate(has_hydraulic_diffusion, needs_slip_rate, trial_reference_geometry, initial_equilibrium, \
                     needs_projected_velocity, shear_heating)
#endif
    #pragma acc parallel loop gang vector firstprivate(shear_heating) async // TODO: ACC: CPU and GPU results are differet because of using 3x3 in elasto_plastic
    for (int e = 0; e < var.nelem; e++) {
        ConstConnAccessor conn = (*var.connectivity)[e];

        // Centroid interpolations, each gated on the branch that reads it: ungated they cost
        // 2 + NDIMS nodal gathers per element per step in every run. Absolute pressure
        // is needed by Bishop increments and the effective-stress strength criterion.
        double dpp = 0.0;
        double pore_strength_pressure = 0;
        const bool increment = has_hydraulic_diffusion && !var.fluid_pressure_consumed;
        const bool strength = param.control.has_hydraulic_diffusion &&
                              (var.mat->rheol_type & MatProps::rh_plastic);
        if (increment || strength) {
            const auto pressure = element_pressure<false>(var,e,ppressure,dppressure,increment,strength);
            // Biot-weighted pore-pressure INCREMENT, not a rate: update_pore_pressure()
            // folded dt in, and stores dppressure with the sign OPPOSITE to its change
            // to ppressure. elastic_effective() adds this to the diagonal components.
            // Initial equilibrium must test skeleton strength too; it only
            // suppresses consumption of pending pressure increments.
            const double alpha = var.mat->alpha_biot(e);
            if (increment) dpp = pressure.bishop * alpha;
            if (strength) pore_strength_pressure = pressure.strength * alpha;
        }

        // No vz in 2D: compute_slip_rate2 takes the vertical component as its second
        // argument, so vy carries it and a vz here would only ever be written.
        double vx = 0.0, vy = 0.0;
#ifdef THREED
        double vz = 0.0;
#endif
        if (needs_projected_velocity) {
            #pragma acc loop seq
            for (int j = 0; j < NODES_PER_ELEM; ++j) {
                vx += vel[conn[j]][0] / double(NODES_PER_ELEM);
                vy += vel[conn[j]][1] / double(NODES_PER_ELEM);
#ifdef THREED
                vz += vel[conn[j]][2] / double(NODES_PER_ELEM);
#endif
            }
        }

        // stress, strain and strain_rate of this element
        TensorAccessor s = stress[e];
        double& syy = stressyy[e];
        TensorAccessor es = strain[e];
        TensorAccessor edot = strain_rate[e];
        double old_s = trace(s);
        if (trial_reference_geometry && var.mat->is_plane_strain) old_s += syy;
        if (has_hydraulic_diffusion && var.old_temperature) {
            // alpha_solid is volumetric, consistent with the density law.
            const double thermal = var.mat->bulkm(e) * var.mat->alpha_solid(e) *
                                   hydraulic_temperature_change(var, e);
            for (int d = 0; d < NDIMS; ++d) s[d] -= thermal;
            if (var.mat->is_plane_strain) syy -= thermal;
        }

        double heat_start[NSTR];
        const double heat_syy_start = syy;
        const double pore_start_trace = var.pore_plastic_increment ?
            trace(s)+(var.mat->is_plane_strain ? syy : 0) : 0;
        if (shear_heating)
            for (int d = 0; d < NSTR; ++d) heat_start[d] = s[d];
        bool heat_maxwell = param.mat.rheol_type == MatProps::rh_maxwell;
        bool heat_has_syy = var.mat->is_plane_strain;

        // anti-mesh locking correction on strain rate
        {
            double div = trace(edot);
            #pragma acc loop seq
            for (int i=0; i<NDIMS; ++i) {
                edot[i] += ((*var.edvoldt)[e] - div) / NDIMS;  // XXX: should NDIMS -> 3 in plane strain?
            }
        }

        double continuum_rate = 0.0;
        if (needs_slip_rate &&
            param.control.rsf_slip_rate_projection_option ==
                rsf_slip_rate_projection_total_strain_rate) {
            const double width = elem_band_width(
                var.coord->view_const(conn), (*var.volume)[e]);
            continuum_rate = 2.0 * width * second_invariant(edot);
        }

        // update strain with strain rate
        #pragma acc loop seq
        for (int i=0; i<NSTR; ++i) {
            es[i] += edot[i] * var.dt;
        }

        // modified strain increment
        double de[NSTR];
        #pragma acc loop seq
        for (int i=0; i<NSTR; ++i) {
            de[i] = edot[i] * var.dt;
        }

        // delta_plstrain is the plastic-strain increment of THIS step (output as
        // "plastic strain-rate" after /dt, and read by the earthquake-state check). Only the
        // plastic branches write it, so clear it here: an evp/evp_rsf element that takes the
        // viscous branch, or an elastic/viscous/maxwell element, would otherwise keep the value
        // from the last time it yielded forever. That stale value is then smeared across the
        // mesh by the nearest-neighbour remap at each remesh, growing a spurious nonzero halo
        // that diffuses out of the yielding zones. Resetting to 0 makes it a true per-step rate.
        delta_plstrain[e] = 0.;
        if (var.pore_plastic_increment) (*var.pore_plastic_increment)[e] = 0;

        switch (param.mat.rheol_type) {
        case MatProps::rh_elastic:
            {
                double bulkm = var.mat->bulkm(e);
                double shearm = var.mat->shearm(e);
                if (has_hydraulic_diffusion)
                {
                    if (var.mat->is_plane_strain)
                        elastic_effective(bulkm, shearm, de, s, dpp, syy);
                    else
                        elastic_effective(bulkm, shearm, de, s, dpp);
                }
                else
                {
                    if (var.mat->is_plane_strain)
                        elastic(bulkm, shearm, de, s, syy);
                    else
                        elastic(bulkm, shearm, de, s);
                }
            }
            break;
        case MatProps::rh_viscous:
            {
                heat_has_syy = trial_reference_geometry && var.mat->is_plane_strain;
                double bulkm = var.mat->bulkm(e);
                viscosity[e] = var.mat->visc(e);
                double total_dv = trace(es);
                double pressure = 0;
                if (trial_reference_geometry && has_hydraulic_diffusion) {
                    for (int j=0; j<NODES_PER_ELEM; ++j)
                        pressure += ppressure[conn[j]] / NODES_PER_ELEM;
                    pressure *= var.mat->alpha_biot(e);
                }
                viscous(bulkm, viscosity[e], total_dv, edot, s,
                        trial_reference_geometry && var.mat->is_plane_strain ? &syy : nullptr, pressure);
            }
            break;
        case MatProps::rh_maxwell:
            {
                heat_has_syy = trial_reference_geometry && var.mat->is_plane_strain;
                double bulkm = var.mat->bulkm(e);
                double shearm = var.mat->shearm(e);
                viscosity[e] = var.mat->visc(e);
                double dv = trial_reference_geometry ? trace(edot)*var.dt
                    : (*var.volume)[e] / (*var.volume_old)[e] - 1;
                maxwell(bulkm, shearm, viscosity[e], var.dt, dv, de, s,
                        trial_reference_geometry && var.mat->is_plane_strain ? &syy : nullptr,
                        trial_reference_geometry ? dpp : 0);
            }
            break;
        case MatProps::rh_ep:
            {
                double depls = 0;
                double bulkm = var.mat->bulkm(e);
                double shearm = var.mat->shearm(e);
                double amc, anphi, anpsi, hardn, ten_max;
                var.mat->plastic_props(e, plstrain[e],
                                       amc, anphi, anpsi, hardn, ten_max);
                pore_strength(pore_strength_pressure, anphi, amc, ten_max);
                int failure_mode;
                if (var.mat->is_plane_strain) {
                    elasto_plastic2d(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                     de, depls, s, syy, failure_mode, 
                                     has_hydraulic_diffusion, dpp);
                }
                else {
                    elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                   de, depls, s, failure_mode, 
                                   has_hydraulic_diffusion, dpp);
                }
                plstrain[e] += depls;
                delta_plstrain[e] = depls;
            }
            break;
        case MatProps::rh_evp:
            {
                double depls = 0;
                double bulkm = var.mat->bulkm(e);
                double shearm = var.mat->shearm(e);
                viscosity[e] = var.mat->visc(e);
                double dv = trial_reference_geometry ? trace(edot)*var.dt
                    : (*var.volume)[e] / (*var.volume_old)[e] - 1;
                // stress due to maxwell rheology
                double sv[NSTR];
                #pragma acc loop seq
                for (int i=0; i<NSTR; ++i) sv[i] = s[i];
                double svyy = syy;
                const bool plane_trial = trial_reference_geometry && var.mat->is_plane_strain;
                maxwell(bulkm, shearm, viscosity[e], var.dt, dv, de, sv,
                        plane_trial ? &svyy : nullptr, trial_reference_geometry ? dpp : 0);
                double svII = plane_trial ? second_invariant2(sv, svyy) : second_invariant2(sv);

                double amc, anphi, anpsi, hardn, ten_max;
                var.mat->plastic_props(e, plstrain[e],
                                       amc, anphi, anpsi, hardn, ten_max);
                pore_strength(pore_strength_pressure, anphi, amc, ten_max);
                // stress due to elasto-plastic rheology
                double sp[NSTR], spyy;
                #pragma acc loop seq
                for (int i=0; i<NSTR; ++i) sp[i] = s[i];
                int failure_mode;
                if (var.mat->is_plane_strain) {
                    spyy = syy;
                    elasto_plastic2d(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                     de, depls, sp, spyy, failure_mode, 
                                     has_hydraulic_diffusion, dpp);
                }
                else {
                    elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                   de, depls, sp, failure_mode, 
                                   has_hydraulic_diffusion, dpp);
                }
                double spII = plane_trial ? second_invariant2(sp, spyy) : second_invariant2(sp);

                // use the smaller as the final stress
                if (svII < spII) {
                    heat_maxwell = true;
                    heat_has_syy = plane_trial;
                    #pragma acc loop seq
                    for (int i=0; i<NSTR; ++i) s[i] = sv[i];
                    if (plane_trial) syy = svyy;
                }
                else {
                    #pragma acc loop seq
                    for (int i=0; i<NSTR; ++i) s[i] = sp[i];
                    plstrain[e] += depls;
                    delta_plstrain[e] = depls;
                    if (var.mat->is_plane_strain) syy = spyy;
                }
            }
            break;
        case MatProps::rh_ep_rsf: // rate-and-state frition model
            {
                double depls = 0;
                double bulkm = var.mat->bulkm(e);
                double shearm = var.mat->shearm(e);

                // calculate the shear direction in maximum shear stress
                double slip_rate;
                
#ifdef THREED
                select_slip_rate3(
                    param.control.rsf_slip_rate_projection_option,
                    continuum_rate, s, vx, vy, vz, slip_rate);
#else
                select_slip_rate2(
                    param.control.rsf_slip_rate_projection_option,
                    continuum_rate, s, vx, vy, slip_rate);
#endif
        
                double amc, anphi, anpsi, hardn, ten_max;
                var.mat->plastic_props_rsf(e, plstrain[e],
                                       amc, anphi, anpsi, hardn, ten_max, slip_rate,
                                       dyn_fric_coeff[e], state_variable[e],
                                       var.dt, param.mat.state_var_model, initial_equilibrium);
                pore_strength(pore_strength_pressure, anphi, amc, ten_max);
                int failure_mode;
                if (var.mat->is_plane_strain) {
                    elasto_plastic2d(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                     de, depls, s, syy, failure_mode, 
                                     has_hydraulic_diffusion, dpp);
                }
                else {
                    elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                   de, depls, s, failure_mode, 
                                   has_hydraulic_diffusion, dpp);
                }
                plstrain[e] += depls;
                delta_plstrain[e] = depls;
            }
            break;
        case MatProps::rh_evp_rsf: // rate-and-state frition model
            {
                double depls = 0;
                double bulkm = var.mat->bulkm(e);
                double shearm = var.mat->shearm(e);
                viscosity[e] = var.mat->visc(e);
                double dv = trial_reference_geometry ? trace(edot)*var.dt
                    : (*var.volume)[e] / (*var.volume_old)[e] - 1;
                // stress due to maxwell rheology
                double sv[NSTR];
                #pragma acc loop seq
                for (int i=0; i<NSTR; ++i) sv[i] = s[i];
                double svyy = syy;
                const bool plane_trial = trial_reference_geometry && var.mat->is_plane_strain;
                maxwell(bulkm, shearm, viscosity[e], var.dt, dv, de, sv,
                        plane_trial ? &svyy : nullptr, trial_reference_geometry ? dpp : 0);
                double svII = plane_trial ? second_invariant2(sv, svyy) : second_invariant2(sv);

                 // calculate the shear direction in maximum shear stress
                double slip_rate;
                
#ifdef THREED
                select_slip_rate3(
                    param.control.rsf_slip_rate_projection_option,
                    continuum_rate, s, vx, vy, vz, slip_rate);
#else
                select_slip_rate2(
                    param.control.rsf_slip_rate_projection_option,
                    continuum_rate, s, vx, vy, slip_rate);
#endif

                double amc, anphi, anpsi, hardn, ten_max;
                var.mat->plastic_props_rsf(e, plstrain[e],
                                       amc, anphi, anpsi, hardn, ten_max, slip_rate,
                                       dyn_fric_coeff[e], state_variable[e],
                                       var.dt, param.mat.state_var_model, initial_equilibrium);
                pore_strength(pore_strength_pressure, anphi, amc, ten_max);
                // stress due to elasto-plastic rheology
                double sp[NSTR], spyy;
                #pragma acc loop seq
                for (int i=0; i<NSTR; ++i) sp[i] = s[i];
                int failure_mode;
                if (var.mat->is_plane_strain) {
                    spyy = syy;
                    elasto_plastic2d(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                     de, depls, sp, spyy, failure_mode, 
                                     has_hydraulic_diffusion, dpp);
                }
                else {
                    elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                   de, depls, sp, failure_mode, 
                                   has_hydraulic_diffusion, dpp);
                }
                double spII = plane_trial ? second_invariant2(sp, spyy) : second_invariant2(sp);

                // use the smaller as the final stress
                if (svII < spII) {
                    heat_maxwell = true;
                    heat_has_syy = plane_trial;
                    #pragma acc loop seq
                    for (int i=0; i<NSTR; ++i) s[i] = sv[i];
                    if (plane_trial) syy = svyy;
                }
                else {
                    #pragma acc loop seq
                    for (int i=0; i<NSTR; ++i) s[i] = sp[i];
                    plstrain[e] += depls;
                    delta_plstrain[e] = depls;
                    if (var.mat->is_plane_strain) syy = spyy;
                }
            }
            break;
        default:
//            std::cerr << "Error: unknown rheology type: " << rheol_type << "\n";
//            std::exit(1);
            break;
        }
        if (var.pore_plastic_increment && !initial_equilibrium && !heat_maxwell) {
            // The signed plastic volume is total strain minus the elastic
            // bulk response. Plane strain includes the hidden normal stress.
            const int dim = var.mat->is_plane_strain ? 3 : NDIMS;
            const double final_trace = trace(s)+(var.mat->is_plane_strain ? syy : 0);
            const double bulk = dim*(var.mat->bulkm(e)-2*var.mat->shearm(e)/3)+2*var.mat->shearm(e);
            (*var.pore_plastic_increment)[e] = trace(de)-(final_trace-pore_start_trace-dim*dpp)/bulk;
        }
        if (shear_heating) {
            double heat = 0;
            if (param.mat.rheol_type == MatProps::rh_viscous || heat_maxwell)
                heat = viscous_heat(viscosity[e], var.dt, edot, heat_start, s,
                                    heat_maxwell, heat_has_syy, heat_syy_start, syy);
            else if (var.mat->rheol_type & MatProps::rh_plastic)
                heat = plastic_heat(var.mat->bulkm(e), var.mat->shearm(e), heat_start, s,
                                    de, dpp, heat_has_syy, heat_syy_start, syy,
                                    pore_strength_pressure + 0.5*dpp);
            (*var.shear_heat)[e] = heat / var.dt;
        }
        if (param.control.is_using_mixed_stress) {
            const double new_trace = trial_reference_geometry && var.mat->is_plane_strain
                ? trace(s)+syy : trace(s);
            dpressure[e] = new_trace-old_s;
        }
        // std::cerr << "stress " << e << ": ";
        // print(std::cerr, s, NSTR);
        // std::cerr << '\n';
    }
#ifdef NPROF
    nvtxRangePop();
#endif
}
