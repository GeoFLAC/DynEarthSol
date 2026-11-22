#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>

#include "3x3-C/dsyevc3.h"
#include "3x3-C/dsyevh3.h"

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
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

// The stored plastic history is also the independent variable of hardening.
#pragma acc routine seq
static double plastic_increment(const double plastic[3], bool plane_strain)
{
    if (plane_strain)
        return 0.5 * (std::max(plastic[0], std::max(plastic[1], plastic[2])) -
                      std::min(plastic[0], std::min(plastic[1], plastic[2])));
    double mean = 0;
    for (int k = 0; k < NDIMS; ++k) mean += plastic[k];
    mean /= NDIMS;
    double squared = mean*mean;
    for (int k = 0; k < NDIMS; ++k) squared += (plastic[k]-mean)*(plastic[k]-mean);
    return std::sqrt(0.5*squared);
}

#pragma acc routine seq
static bool admissible_shear(const double p[3], double increment, double anphi,
                             double amc, double hardening, double ten_max, double tolerance)
{
    const double strength = amc + hardening*increment;
    return std::isfinite(increment) && std::isfinite(strength) && strength >= 0 &&
        std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]) &&
        p[0] <= p[1]+tolerance && p[1] <= p[2]+tolerance && p[2] <= ten_max+tolerance &&
        std::fabs(p[0]-anphi*p[2]+strength) <= tolerance;
}

// Solve both adjacent faces. Hardening is per unit stored plastic increment,
// rather than per unit sum of face multipliers.
#pragma acc routine seq
static bool shear_edge_return(const double trial[3], double p[3], double& depls,
                              double a1, double a2, double anphi, double anpsi,
                              double amc, double hardening, double ten_max,
                              bool plane_strain, double tolerance)
{
    for (int edge = 0; edge < 2; ++edge) {
        const int low[2] = {0, edge == 0 ? 1 : 0};
        const int high[2] = {2, edge == 0 ? 2 : 1};
        double flow[2][3], response[2][3];
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 3; ++k) {
                flow[j][k] = (k == high[j] ? anpsi : 0) - (k == low[j] ? 1 : 0);
                response[j][k] = a2*(anpsi-1) + (a1-a2)*flow[j][k];
            }
        double matrix[2][2], rhs[2];
        for (int i = 0; i < 2; ++i) {
            rhs[i] = anphi*trial[high[i]] - trial[low[i]] - amc;
            for (int j = 0; j < 2; ++j)
                matrix[i][j] = anphi*response[j][high[i]] - response[j][low[i]];
        }
        const double determinant = matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0];
        if (!(determinant > 0) || !std::isfinite(determinant)) continue;
        double multiplier[2];
        if (hardening == 0) {
            multiplier[0] = (rhs[0]*matrix[1][1] - rhs[1]*matrix[0][1])/determinant;
            multiplier[1] = (rhs[1]*matrix[0][0] - rhs[0]*matrix[1][0])/determinant;
        }
        else {
            // Isotropic elasticity gives a symmetric 2x2 matrix. In sum/difference
            // coordinates, d is fixed and E*s + H*depls(s,d) = R remains scalar.
            const double difference = (rhs[0]-rhs[1])/(matrix[0][0]-matrix[0][1]);
            const double elastic = 0.5*(matrix[0][0]+matrix[0][1]);
            const double r = 0.5*(rhs[0]+rhs[1])/elastic;
            const double h = hardening/elastic;
            double sum;
            if (plane_strain) {
                const double a = edge == 0 ? 0.5*anpsi+0.25 : 0.5+0.25*anpsi;
                const double b = edge == 0 ? 0.25 : 0.25*anpsi;
                if (!(1+h*a > 0)) continue;
                sum = (r-h*b*std::fabs(difference))/(1+h*a);
            }
            else {
                // depls^2 = A*s^2 + B*d^2. Solve its positive quadratic root;
                // reject loss of the local softening tangent instead of rotating
                // a crossed single-face correction back into the mesh.
                double along[3], across[3];
                for (int k = 0; k < 3; ++k) {
                    along[k] = 0.5*(flow[0][k]+flow[1][k]);
                    across[k] = 0.5*(flow[0][k]-flow[1][k]);
                }
                const double qa = plastic_increment(along,false);
                const double qb = plastic_increment(across,false);
                const double a = qa*qa, c = qb*qb*difference*difference;
                if (hardening < 0 && !(1-a*h*h > 0)) continue;
                const double squared = a*r*r+c;
                const double discriminant = squared-a*h*h*c;
                if (!(discriminant >= 0)) continue;
                const double denominator = std::sqrt(discriminant)+a*r*h;
                if (!(denominator > 0)) continue;
                sum = r-h*squared/denominator;
            }
            multiplier[0] = 0.5*(sum+difference);
            multiplier[1] = 0.5*(sum-difference);
        }
        if (!std::isfinite(multiplier[0]) || !std::isfinite(multiplier[1]) ||
            multiplier[0] < 0 || multiplier[1] < 0) continue;
        double candidate[3], plastic[3];
        for (int k = 0; k < 3; ++k) {
            candidate[k] = trial[k]-multiplier[0]*response[0][k]-multiplier[1]*response[1][k];
            plastic[k] = multiplier[0]*flow[0][k]+multiplier[1]*flow[1][k];
        }
        const double increment = plastic_increment(plastic,plane_strain);
        if (!admissible_shear(candidate,increment,anphi,amc,hardening,ten_max,tolerance)) continue;
        const int repeated = edge == 0 ? 0 : 1;
        candidate[repeated] = candidate[repeated+1] =
            0.5*(candidate[repeated]+candidate[repeated+1]);
        for (int k = 0; k < 3; ++k) p[k] = candidate[k];
        depls = increment;
        return true;
    }
    return false;
}

// Absolute tolerance of the yield, ordering and tensile-cap checks over the
// nd principal stresses in use (2 in pure 2D, else 3).
#pragma acc routine seq
static double return_tolerance(const double p[3], double amc, int nd = 3)
{
    double scale = std::fabs(amc);
    for (int k = 0; k < nd; ++k) scale = std::max(scale, std::fabs(p[k]));
    return 1e-12*scale;
}

// Rejected candidates never replace the caller's stress or plastic increment.
#pragma acc routine seq
static bool shear_return(const double trial[3], double p[3], double& depls,
                         double a1, double a2, double anphi, double anpsi,
                         double amc, double hardening, double ten_max,
                         bool plane_strain)
{
    const double tolerance = return_tolerance(trial, amc);
    const double fs = trial[0]-trial[2]*anphi+amc;
    const double unit[3] = {-1,0,anpsi};
    const double unit_increment = plastic_increment(unit,plane_strain);
    const double denominator = a1-a2*anpsi+a1*anphi*anpsi-a2*anphi+hardening*unit_increment;
    if (denominator > 0 && std::isfinite(denominator)) {
        const double alam = fs/denominator;
        const double candidate[3] = {trial[0]-alam*(a1-a2*anpsi),
                                    trial[1]-alam*(a2-a2*anpsi),
                                    trial[2]-alam*(a2-a1*anpsi)};
        const double increment = std::fabs(alam)*unit_increment;
        if (alam <= 0 && admissible_shear(candidate,increment,anphi,amc,hardening,ten_max,tolerance)) {
            for (int k = 0; k < 3; ++k) p[k] = candidate[k];
            depls = increment;
            return true;
        }
    }
    return shear_edge_return(trial,p,depls,a1,a2,anphi,anpsi,amc,hardening,ten_max,
                             plane_strain,tolerance);
}

// Pure 2D (two principal stresses): the single shear face. The strength
// changes by hardening*depls, with depls = |alam|*u the recorded measure, so
// the hardening term carries u. Accepted only below the tensile cap and in
// principal order; the caller otherwise uses consistent_return.
#pragma acc routine seq
static bool shear_return2d(const double trial[3], double p[3], double& depls,
                           double a1, double a2, double anphi, double anpsi,
                           double amc, double hardening, double ten_max)
{
    /* plastic strain in the principal directions: alam and -alam*anpsi, so
     * depls = sqrt(((d1-dm)^2 + (d2-dm)^2 + dm^2)/2) with dm = (d1+d2)/2
     * reduces to |alam|*u. */
    const double u = std::sqrt((3 + 2*anpsi + 3*anpsi*anpsi) / 8);
    const double denominator = a1 - a2*anpsi + a1*anphi*anpsi - a2*anphi + hardening*u;
    if (!(denominator > 0) || !std::isfinite(denominator)) return false;
    const double fs = trial[0] - trial[1]*anphi + amc;
    const double alam = fs / denominator;
    const double candidate[2] = {trial[0] - alam * (a1 - a2 * anpsi),
                                 trial[1] - alam * (a2 - a1 * anpsi)};
    const double increment = std::fabs(alam) * u;
    const double tolerance = return_tolerance(trial, amc, 2);
    if (!(alam <= 0) || !std::isfinite(candidate[0]) || !std::isfinite(candidate[1]) ||
        !(candidate[0] <= candidate[1] + tolerance) ||
        !(candidate[1] <= ten_max + tolerance) ||
        !(amc + hardening*increment >= 0))
        return false;
    p[0] = candidate[0];
    p[1] = candidate[1];
    depls = increment;
    return true;
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


// Why no admissible Mohr-Coulomb state was found; failure_mode = -reason.
enum ReturnRejection {
    reject_nonfinite = 1,    // trial stress is not finite
    reject_no_return = 2,    // no set of shear faces and tensile caps gives nonnegative
                             // multipliers and a stress inside every surface
    reject_no_root = 3,      // the law has more pieces than consistent_return searches
    reject_no_increment = 4, // no plastic increment returns the stress to the strength of the
                             // history it stores: the law jumps up there
    reject_budget = 5,       // the search used up return_budget before it finished
    reject_reason_count = 6
};

#pragma acc routine seq
static bool finite3(const double p[3])
{
    return std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]);
}

// Strength after a plastic increment delta: plastic_props at pls+delta.
// Rate-and-state strength (rsf) evaluates the same
// piecewise cohesion, friction, dilation and cap law with this step's slip
// rate, state and friction coefficient held fixed.
struct StrengthLaw {
    const MatProps* mat;
    int e;
    double pls;
    double amc, anphi, anpsi, ten_max;   // at delta = 0
    bool rsf;
    double slip_rate, state_variable, dyn_fric_coeff;   // rsf only
    int state_model;
    bool frozen_history;
    const MatProps::StrengthPiece* piece;   // the piece at pls, if already prepared
    // control.iterative_return_mapping_option 0: friction and dilation stay at
    // their values at pls, and only cohesion follows the increment.
    bool start_angles;
};

struct Strength {
    double amc, anphi, anpsi, ten_max;
};

// The strength law along this step's plastic increment delta, one linear
// piece of the weakening law at a time: each piece is prepared once and
// reused for every trial increment on it. The law declares through the piece
// which strength parameters change with delta; the returns below never ask
// which material or friction law produced it.
struct StrengthCurve {
    const StrengthLaw& law;
    MatProps::StrengthPiece piece;
    double from;   // increment at which the piece starts
    // Friction and dilation at delta = 0, for law.start_angles.
    double start_friction, start_dilation;
    MatProps::AngleFactors start_factors;
    // What the search of consistent_return has evaluated on this curve (spend).
    mutable int evaluations;
};

// For law.start_angles: the angles of every piece are those at delta = 0.
#pragma acc routine seq
static void hold_start_angles(StrengthCurve& curve)
{
    MatProps::StrengthPiece& piece = curve.piece;
    if (curve.from == 0) {
        curve.start_friction = piece.friction;
        curve.start_dilation = piece.dilation;
        curve.start_factors = piece.factors;
    }
    else {
        piece.friction = curve.start_friction;
        piece.dilation = curve.start_dilation;
        piece.factors = curve.start_factors;
    }
    piece.dfriction = piece.ddilation = 0;
    curve.law.mat->piece_strength(piece, 0, piece.amc, piece.anphi, piece.anpsi, piece.ten_max);
}

#pragma acc routine seq
static void prepare_piece(StrengthCurve& curve, double from)
{
    const StrengthLaw& law = curve.law;
    if (law.rsf)
        law.mat->rsf_strength_piece(law.e, law.pls+from, law.slip_rate, law.state_variable,
                                    law.dyn_fric_coeff, law.state_model, law.frozen_history,
                                    curve.piece);
    else
        law.mat->strength_piece(law.e, law.pls+from, curve.piece);
    curve.from = from;
    if (law.start_angles) hold_start_angles(curve);
}

// The curve at delta = 0, from the piece the caller prepared if there is one.
#pragma acc routine seq
static StrengthCurve start_curve(const StrengthLaw& law)
{
    StrengthCurve curve = {law, {}, 0};
    if (law.piece) {
        curve.piece = *law.piece;
        law.mat->complete_piece(law.e, curve.piece);
        if (law.start_angles) hold_start_angles(curve);
    }
    else prepare_piece(curve, 0);
    return curve;
}

// Strength after a plastic increment delta on the prepared piece; at delta =
// 0 it equals plastic_props at pls.
#pragma acc routine seq
static Strength strength_at(const StrengthCurve& curve, double delta)
{
    Strength r;
    curve.law.mat->piece_strength(curve.piece, delta - curve.from,
                                  r.amc, r.anphi, r.anpsi, r.ten_max);
    return r;
}

static constexpr double no_break = std::numeric_limits<double>::max();

// Plastic increment at which the prepared piece ends. Rounded up so that pls
// + increment reaches the breakpoint: next-pls alone can land one ulp short,
// and the next piece would then start on this one.
#pragma acc routine seq
static double next_break(const StrengthCurve& curve)
{
    const double next = curve.piece.end;
    if (next == no_break) return no_break;
    const double pls = curve.law.pls;
    double increment = next - pls;
    while (pls + increment < next)
        increment += std::numeric_limits<double>::epsilon()*next;
    return increment;
}

#pragma acc routine seq
static bool near_value(double x, double y)
{
    return std::fabs(x - y) <= 1e-12*std::max(std::fabs(x), std::fabs(y));
}

// Whether two strengths agree to the relative 1e-12 of the yield checks: the
// two pieces of a continuous law at their break differ only by rounding.
#pragma acc routine seq
static bool same_strength(const Strength& a, const Strength& b)
{
    return near_value(a.amc, b.amc) && near_value(a.anphi, b.anphi) &&
        near_value(a.anpsi, b.anpsi) && near_value(a.ten_max, b.ten_max);
}

// Whether the plastic strain that an increment stores, pls + recorded, reaches
// the end of the prepared piece where the law jumps (pls0 = pls1): its strength
// is then the next piece's, not this one's. Where the law is continuous there,
// the two differ only by rounding.
#pragma acc routine seq
static bool reaches_jump(const StrengthCurve& curve, double recorded)
{
    const double end = curve.piece.end;
    if (end == no_break || curve.law.pls + recorded < end) return false;
    const double at = next_break(curve);
    StrengthCurve next = curve;
    prepare_piece(next, at);
    return !same_strength(strength_at(curve, at), strength_at(next, at));
}

// The closed-form returns hold friction and dilation fixed, so they apply
// only where neither changes along the piece (always with law.start_angles).
// Elsewhere their result is the starting guess of consistent_return.
#pragma acc routine seq
static bool closed_form_applies(const StrengthCurve& curve)
{
    return !curve.piece.angles_change();
}

// The closed-form returns use the law's linear piece at pls: a fixed cohesion
// slope and fixed friction, dilation and tensile cap. Keep their result only
// while the increment stays on that piece and the cap of the stored history
// still bounds it, or still holds it when it was returned to the cap.
#pragma acc routine seq
static bool closed_form_consistent(const StrengthCurve& curve, const double p[3], double depls,
                                   int nd = 3)
{
    if (depls == 0) return true;
    if (!(depls <= next_break(curve)) || reaches_jump(curve, depls)) return false;
    const StrengthLaw& law = curve.law;
    const Strength end = strength_at(curve, depls);
    const double tolerance = return_tolerance(p, end.amc, nd);
    const double top = p[nd-1];
    const bool on_cap = top >= law.ten_max - tolerance;
    return end.anphi == law.anphi && end.anpsi == law.anpsi &&
        top <= end.ten_max + tolerance &&
        (!on_cap || top >= end.ten_max - tolerance);
}

// Surface s of the ordered region of nd principal stresses, grad.p - offset
// <= 0, with plastic flow direction flow. For nd = 3, 0..2 are Mohr-Coulomb
// faces with (low, high) = (0,2), (1,2), (0,1) and 3..5 the tensile caps on
// p2, p1, p0. For nd = 2 (pure 2D), 0 is the face (0,1) and 1..2 the caps on
// p1, p0.
#pragma acc routine seq
static int return_surface_count(int nd)
{
    return nd == 3 ? 6 : 3;
}

#pragma acc routine seq
static int return_face_count(int nd)
{
    return nd == 3 ? 3 : 1;
}

// The principal stresses (low, high) of face s.
#pragma acc routine seq
static void face_stresses(int s, int nd, int& low, int& high)
{
    low = s == 1 ? 1 : 0;
    high = nd == 2 ? 1 : (s == 2 ? 1 : 2);
}

#pragma acc routine seq
static double return_surface(int s, const Strength& r, int nd, double grad[3], double flow[3])
{
    const int faces = return_face_count(nd);
    for (int k = 0; k < 3; ++k) grad[k] = flow[k] = 0;
    if (s < faces) {
        int low, high;
        face_stresses(s, nd, low, high);
        grad[high] = r.anphi;
        grad[low] = -1;
        flow[high] = r.anpsi;
        flow[low] = -1;
        return r.amc;
    }
    const int cap = nd-1 - (s-faces);
    grad[cap] = flow[cap] = 1;
    return r.ten_max;
}

// How far the stress p lies outside surface s, grad.p - offset.
#pragma acc routine seq
static double surface_violation(int s, const Strength& r, int nd, const double p[3])
{
    double grad[3], flow[3];
    double f = -return_surface(s, r, nd, grad, flow);
    for (int k = 0; k < nd; ++k) f += grad[k]*p[k];
    return f;
}

// Gaussian elimination with partial pivoting; false if (near) singular.
#pragma acc routine seq
static bool solve_small(int n, double m[3][3], double b[3], double x[3])
{
    double scale = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) scale = std::max(scale, std::fabs(m[i][j]));
    if (!(scale > 0)) return false;
    for (int c = 0; c < n; ++c) {
        int pivot = c;
        for (int r = c+1; r < n; ++r)
            if (std::fabs(m[r][c]) > std::fabs(m[pivot][c])) pivot = r;
        if (!(std::fabs(m[pivot][c]) > 1e-12*scale)) return false;
        for (int k = 0; k < n; ++k) {
            const double t = m[c][k]; m[c][k] = m[pivot][k]; m[pivot][k] = t;
        }
        const double t = b[c]; b[c] = b[pivot]; b[pivot] = t;
        for (int r = c+1; r < n; ++r) {
            const double factor = m[r][c]/m[c][c];
            for (int k = c; k < n; ++k) m[r][k] -= factor*m[c][k];
            b[r] -= factor*b[c];
        }
    }
    for (int r = n-1; r >= 0; --r) {
        double v = b[r];
        for (int k = r+1; k < n; ++k) v -= m[r][k]*x[k];
        x[r] = v/m[r][r];
    }
    return true;
}

// Return for a fixed strength (Koiter's rule on the composite surface): the
// first active set, fewest surfaces first, whose multipliers are nonnegative
// and whose stress keeps the trial order and satisfies every surface. This
// covers faces, edges, the tensile cap, its edge and apex, and their corners
// (in pure 2D, nd = 2: the face, the cap, and the face-cap and cap-cap corners).
// nd is a template argument so that the loops over principal stresses and
// surfaces unroll; fixed_strength_return below dispatches on it.
#pragma acc routine seq
template <int ND>
static bool fixed_strength_return_nd(const double trial[3], double a1, double a2,
                                     const Strength& r, bool plane_strain,
                                     double p[3], double& depls, int& set)
{
    double grad[6][3], flow[6][3], response[6][3], offset[6], violation[6];
    const int surfaces = return_surface_count(ND);
    const double tolerance = return_tolerance(trial, r.amc, ND);
    bool admissible = true;
    for (int s = 0; s < surfaces; ++s) {
        offset[s] = return_surface(s, r, ND, grad[s], flow[s]);
        double volume = 0;
        for (int k = 0; k < ND; ++k) volume += flow[s][k];
        violation[s] = -offset[s];
        for (int k = 0; k < 3; ++k) response[s][k] = 0;
        for (int k = 0; k < ND; ++k) {
            response[s][k] = a2*volume + (a1-a2)*flow[s][k];
            violation[s] += grad[s][k]*trial[k];
        }
        admissible = admissible && violation[s] <= tolerance;
    }
    if (admissible) {
        for (int k = 0; k < 3; ++k) p[k] = trial[k];
        depls = 0;
        set = 0;
        return true;
    }
    for (int size = 1; size <= ND; ++size)
        for (int mask = 1; mask < (1 << surfaces); ++mask) {
            int active[3], n = 0;
            for (int s = 0; s < surfaces && n <= 3; ++s)
                if (mask & (1 << s)) {
                    if (n < 3) active[n] = s;
                    ++n;
                }
            if (n != size) continue;
            double m[3][3], b[3], x[3];
            for (int i = 0; i < n; ++i) {
                b[i] = violation[active[i]];
                for (int j = 0; j < n; ++j) {
                    m[i][j] = 0;
                    for (int k = 0; k < ND; ++k) m[i][j] += grad[active[i]][k]*response[active[j]][k];
                }
            }
            if (!solve_small(n, m, b, x)) continue;
            double largest = 0;
            for (int i = 0; i < n; ++i) largest = std::max(largest, std::fabs(x[i]));
            bool ok = std::isfinite(largest);
            for (int i = 0; i < n; ++i) ok = ok && x[i] >= -1e-12*largest;
            if (!ok) continue;
            double candidate[3], plastic[3];
            for (int k = 0; k < 3; ++k) {
                candidate[k] = trial[k];
                plastic[k] = 0;
                for (int i = 0; i < n; ++i) {
                    candidate[k] -= x[i]*response[active[i]][k];
                    plastic[k] += x[i]*flow[active[i]][k];
                }
            }
            for (int k = 0; k+1 < ND; ++k) ok = ok && candidate[k] <= candidate[k+1]+tolerance;
            for (int s = 0; s < surfaces; ++s) {
                double f = -offset[s];
                for (int k = 0; k < ND; ++k) f += grad[s][k]*candidate[k];
                ok = ok && f <= tolerance;
            }
            if (!ok) continue;
            for (int k = 0; k < 3; ++k) p[k] = candidate[k];
            depls = plastic_increment(plastic, plane_strain);
            set = mask;
            return true;
        }
    return false;
}

#pragma acc routine seq
static bool fixed_strength_return(const double trial[3], double a1, double a2,
                                  const Strength& r, bool plane_strain, int nd,
                                  double p[3], double& depls, int& set)
{
    return nd == 3 ? fixed_strength_return_nd<3>(trial, a1, a2, r, plane_strain, p, depls, set)
                   : fixed_strength_return_nd<2>(trial, a1, a2, r, plane_strain, p, depls, set);
}

// A fixed-strength return at plastic increment delta, its residual g = recorded - delta,
// and its active set: the surfaces it returns to, as a mask (0 if elastic).
struct ReturnPoint {
    double delta, g, recorded, p[3];
    int set;
};

// The fixed-strength returns and parts of strength ranges the search of one element's return
// may evaluate: a bound on its work, many times what returns take on the strength laws of the
// tests and examples. Past it the return is rejected (reject_budget), not taken as found.
static constexpr int return_budget = 4096;

#pragma acc routine seq
static bool spend(const StrengthCurve& curve)
{
    return ++curve.evaluations <= return_budget;
}

#pragma acc routine seq
static bool evaluate_return(const double trial[3], double a1, double a2,
                            const StrengthCurve& curve, bool plane_strain, int nd,
                            double delta, ReturnPoint& r)
{
    r.delta = delta;
    if (!spend(curve) ||
        !fixed_strength_return(trial, a1, a2, strength_at(curve, delta), plane_strain, nd,
                               r.p, r.recorded, r.set))
        return false;
    r.g = r.recorded - delta;
    return true;
}

#pragma acc routine seq
static int accept_return(const ReturnPoint& r, double p[3], double& depls)
{
    for (int k = 0; k < 3; ++k) p[k] = r.p[k];
    depls = r.recorded;
    return 0;
}

// Converged once g is below the relative 1e-12 of the yield checks
// (return_tolerance): the strength then misses the recorded increment by far
// less than those checks resolve.
#pragma acc routine seq
static bool root_found(const ReturnPoint& r)
{
    return std::fabs(r.g) <= 1e-12*std::max(r.delta, r.recorded);
}

static constexpr int keep_searching = -1;

// Go on to the next piece at the break hi, from the point end_point there. Where
// the law is continuous, the search goes on from end_point. Where it jumps, the
// break is evaluated again with the next piece's strength: a root there is the
// return, and otherwise the search goes on from there. If g is then below 0, no
// increment near the break records itself: less flow stores a history on this
// piece, more on the next.
#pragma acc routine seq
static int cross_break(const double trial[3], double a1, double a2, StrengthCurve& curve,
                       bool plane_strain, int nd, const ReturnPoint& end_point,
                       ReturnPoint& lo, double& hi, double p[3], double& depls)
{
    const double at = hi;
    const Strength end = strength_at(curve, at);
    prepare_piece(curve, at);
    hi = next_break(curve);
    lo = end_point;
    if (same_strength(end, strength_at(curve, at))) return keep_searching;
    if (!evaluate_return(trial, a1, a2, curve, plane_strain, nd, at, lo))
        return reject_no_return;
    if (root_found(lo)) return accept_return(lo, p, depls);
    return keep_searching;
}

// Accept a root on the prepared piece, unless the history it stores reaches a
// jump of the law; that one goes on at the break as above.
#pragma acc routine seq
static int finish_return(const double trial[3], double a1, double a2, StrengthCurve& curve,
                         bool plane_strain, int nd, const ReturnPoint& r,
                         ReturnPoint& lo, double& hi, double p[3], double& depls)
{
    if (!reaches_jump(curve, r.recorded)) return accept_return(r, p, depls);
    return cross_break(trial, a1, a2, curve, plane_strain, nd, r, lo, hi, p, depls);
}

// A collapsed bracket of g: within rounding of a root where g is continuous;
// a g still far from 0 means the strength jumps inside it.
#pragma acc routine seq
static bool near_root(const ReturnPoint& r)
{
    return std::fabs(r.g) <= 1e-9*std::max(r.delta, r.recorded);
}

// Between a and b on one piece, where g has the same sign at both and turns at most once,
// x is where g goes furthest the other way: its minimum if it is above 0 there, its maximum
// if below. Golden-section search finds it, and stops at the first point where g has the
// other sign. Returns false if a return fails.
#pragma acc routine seq
static bool find_turn(const double trial[3], double a1, double a2, const StrengthCurve& curve,
                      bool plane_strain, int nd, const ReturnPoint& a, const ReturnPoint& b,
                      ReturnPoint& x)
{
    const double s = a.g > 0 ? 1 : -1;  // the search lowers s*g
    const double r = 0.6180339887498949;  // (sqrt(5) - 1)/2
    double lo = a.delta, hi = b.delta;
    ReturnPoint u, v;
    if (!evaluate_return(trial, a1, a2, curve, plane_strain, nd, hi - r*(hi - lo), u) ||
        !evaluate_return(trial, a1, a2, curve, plane_strain, nd, lo + r*(hi - lo), v))
        return false;
    for (;;) {
        const bool left = s*u.g < s*v.g;
        x = left ? u : v;
        if (s*x.g < 0 || hi - lo <= 1e-9*hi) return true;
        if (left) {
            hi = v.delta;
            v = u;
            if (!evaluate_return(trial, a1, a2, curve, plane_strain, nd, hi - r*(hi - lo), u))
                return false;
        }
        else {
            lo = u.delta;
            u = v;
            if (!evaluate_return(trial, a1, a2, curve, plane_strain, nd, lo + r*(hi - lo), v))
                return false;
        }
    }
}

// Bounds of the largest violation of the trial's surfaces, v, and of the tolerance of the yield
// checks, tol, over the strengths of the prepared piece between increments x and y. Each
// surface is linear in anphi and in amc or ten_max (anpsi, of the flow only, plays no part),
// so it lies within its values at the corners of their ranges. Where the curvature of the
// strength is known, it also lies within curvature (y-x)^2/8 of the chord between its values
// at x and y, which follows cohesion and friction together: the bound narrows with the square
// of y-x even where their changes nearly cancel. Twice the lookup tables' error covers the
// difference between their strength and the smooth one.
#pragma acc routine seq
static void violation_range(const double trial[3], const StrengthCurve& curve, int nd,
                            double x, double y, double v[2], double tol[2])
{
    MatProps::StrengthRange range;
    curve.law.mat->piece_strength_range(curve.piece, x - curve.from, y - curve.from, range);
    const Strength ends[2] = {strength_at(curve, x), strength_at(curve, y)};
    const double chord = 0.125*(y - x)*(y - x);
    v[0] = v[1] = -std::numeric_limits<double>::max();
    for (int s = 0; s < return_surface_count(nd); ++s) {
        double lo = std::numeric_limits<double>::max(), hi = -lo;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j) {
                const Strength corner = {range.amc[j], range.anphi[i], range.anphi[i],
                                         range.ten_max[j]};
                const double f = surface_violation(s, corner, nd, trial);
                lo = std::min(lo, f);
                hi = std::max(hi, f);
            }
        double curvature = range.curvature_ten, error = range.error_ten;
        if (s < return_face_count(nd)) {
            int low, high;
            face_stresses(s, nd, low, high);
            const double scale = std::fabs(trial[high]);
            curvature = range.curvature_anphi < 0 || range.curvature_amc < 0 ? -1 :
                scale*range.curvature_anphi + range.curvature_amc;
            error = scale*range.error_anphi + range.error_amc;
        }
        if (curvature >= 0) {
            const double fx = surface_violation(s, ends[0], nd, trial);
            const double fy = surface_violation(s, ends[1], nd, trial);
            const double slack = curvature*chord + 2*error;
            lo = std::max(lo, std::min(fx, fy) - slack);
            hi = std::min(hi, std::max(fx, fy) + slack);
        }
        v[0] = std::max(v[0], lo);
        v[1] = std::max(v[1], hi);
    }
    tol[0] = return_tolerance(trial, range.amc[0], nd);
    tol[1] = return_tolerance(trial, range.amc[1], nd);
}

// Whether the trial yields at increment delta, as the first check of fixed_strength_return.
#pragma acc routine seq
static bool yields_at(const double trial[3], const StrengthCurve& curve, int nd, double delta)
{
    const Strength r = strength_at(curve, delta);
    const double tolerance = return_tolerance(trial, r.amc, nd);
    for (int s = 0; s < return_surface_count(nd); ++s)
        if (surface_violation(s, r, nd, trial) > tolerance) return true;
    return false;
}

// With friction and dilation fixed, the strength along a piece moves one way, and where the trial
// yields is one range at one end. Where they change, it can yield in several ranges: whether
// it starts or stops yielding between a and c, where it does the same, as the return y at the
// first such increment found. The search halves [a, c], left half first, down to a relative
// 1e-12 of c, except where violation_range shows that the trial yields throughout a part, or
// nowhere in it, and checks the middle of each part it halves as yields_at does. Each part
// counts against return_budget. Returns 1 if it finds a change, 0 once it has shown there is
// none, and -1 if a return fails or the budget runs out.
#pragma acc routine seq
static int yield_change(const double trial[3], double a1, double a2, const StrengthCurve& curve,
                        bool plane_strain, int nd, const ReturnPoint& a, const ReturnPoint& c,
                        ReturnPoint& y)
{
    if (!curve.piece.angles_change() || !(c.delta > a.delta)) return 0;
    const bool yields = a.set != 0;
    double size = c.delta - a.delta, i = 0;  // part i of [a, c] in parts of this size
    int level = 0;
    for (;;) {
        if (!spend(curve)) return -1;
        const double x = a.delta + i*size, z = std::min(c.delta, a.delta + (i + 1)*size);
        double v[2], tol[2];
        violation_range(trial, curve, nd, x, z, v, tol);
        if (!(yields ? v[0] > tol[1] : v[1] <= tol[0]) && size > 1e-12*c.delta) {
            const double m = 0.5*(x + z);
            if (yields_at(trial, curve, nd, m) != yields) {
                if (!evaluate_return(trial, a1, a2, curve, plane_strain, nd, m, y)) return -1;
                if ((y.set != 0) != yields) return 1;
            }
            size *= 0.5;
            i *= 2;
            ++level;
            continue;
        }
        // On to the next part to the right, up from every right half.
        for (; level > 0 && std::fmod(i, 2.) == 1; --level) {
            i = 0.5*(i - 1);
            size *= 2;
        }
        if (level == 0) return 0;
        i += 1;
    }
}

// The first fall of g to 0 between a and b on one piece, as the bracket [first, last] with
// g(first) > 0 >= g(last). Besides where friction or dilation change, g turns where the active
// set of the return changes, and between such changes it turns at most once: with friction
// and dilation fixed, the recorded increment is there the norm of an affine function of
// delta, and g is convex. [a, b] is split where the set changes, found by bisection, and where
// the trial starts or stops yielding (yield_change), and each part is searched as find_turn
// does. Returns 1 if g falls to 0, 0 if not, and -1 if a return fails.
#pragma acc routine seq
static int first_fall(const double trial[3], double a1, double a2, const StrengthCurve& curve,
                      bool plane_strain, int nd, ReturnPoint a, const ReturnPoint& b,
                      ReturnPoint& first, ReturnPoint& last)
{
    for (int part = 0; part < 16; ++part) {
        // [a, c] has a's active set, and d lies just past it; or c and d are b. Where the
        // trial starts or stops yielding within [a, c], d moves there and [a, d] is split again.
        // From an elastic a, the yield check alone tells the side of each step, and the return
        // is solved at d only.
        ReturnPoint c = b, d = b, x, y;
        bool elastic = false;  // the trial yields nowhere in [a, c]
        for (int split = 0; ; ++split) {
            if (a.set != d.set) {
                c = a;
                const double target = d.delta;
                for (int it = 0; it < 60 && d.delta - c.delta > 1e-12*d.delta; ++it) {
                    const double mid = 0.5*(c.delta + d.delta);
                    if (a.set == 0) {
                        if (!spend(curve)) return -1;
                        if (yields_at(trial, curve, nd, mid)) d.delta = mid;
                        else c.delta = mid;
                        continue;
                    }
                    ReturnPoint m;
                    if (!evaluate_return(trial, a1, a2, curve, plane_strain, nd, mid, m))
                        return -1;
                    if (m.set == a.set) c = m;
                    else d = m;
                }
                if (a.set == 0) {
                    c.g = c.recorded - c.delta;
                    if (d.delta != target &&
                        !evaluate_return(trial, a1, a2, curve, plane_strain, nd, d.delta, d))
                        return -1;
                }
            }
            if (split == 16) break;
            const int change = yield_change(trial, a1, a2, curve, plane_strain, nd, a, c, y);
            if (change < 0) return -1;
            if (!change) {
                elastic = a.set == 0;
                break;
            }
            d = y;
        }
        if (a.g > 0 && !(c.g > 0)) {
            first = a;
            last = c;
            return 1;
        }
        // Where the trial does not yield, nothing is recorded: g = -delta only falls.
        if (!elastic && (a.g > 0) == (c.g > 0)) {
            if (!find_turn(trial, a1, a2, curve, plane_strain, nd, a, c, x)) return -1;
            if (a.g > 0 && x.g < 0) {
                first = a;
                last = x;
                return 1;
            }
            if (!(a.g > 0) && x.g > 0) {
                first = x;
                last = c;
                return 1;
            }
        }
        if (c.delta == b.delta) return 0;
        if (c.g > 0 && !(d.g > 0)) {
            first = c;
            last = d;
            return 1;
        }
        a = d;
    }
    return 0;
}

// Illinois regula falsi on [a, b], g(a) > 0 >= g(b), and the return at its root, as
// finish_return gives it. Where the bracket collapses with g still far from 0, g jumps
// across 0 inside it, and the search goes on from below the jump (lo).
#pragma acc routine seq
static int solve_bracket(const double trial[3], double a1, double a2, StrengthCurve& curve,
                         bool plane_strain, int nd, ReturnPoint a, ReturnPoint b,
                         ReturnPoint& lo, double& hi, double p[3], double& depls)
{
    ReturnPoint x;
    int side = 0;
    if (b.g == 0) x = b;
    for (int it = 0; it < 200 && b.g != 0; ++it) {
        double xd = (a.delta*b.g - b.delta*a.g)/(b.g - a.g);
        if (!(xd > a.delta && xd < b.delta)) xd = 0.5*(a.delta + b.delta);
        if (!evaluate_return(trial, a1, a2, curve, plane_strain, nd, xd, x))
            return reject_no_return;
        if (root_found(x)) break;
        // The bracket collapsed: g jumps across 0 inside it if it is far from 0.
        if (b.delta - a.delta <= 1e-15*b.delta) break;
        if (x.g > 0) {
            const double gb = b.g;
            a = x;
            if (side == 1) b.g = 0.5*gb;
            side = 1;
        }
        else {
            const double ga = a.g;
            b = x;
            if (side == -1) a.g = 0.5*ga;
            side = -1;
        }
    }
    if (!near_root(x)) {
        lo = b;
        return keep_searching;
    }
    return finish_return(trial, a1, a2, curve, plane_strain, nd, x, lo, hi, p, depls);
}

// The return on the part of a piece from a to its break b: at the first fall of g to 0
// between them, or otherwise past the break, on the next piece.
#pragma acc routine seq
static int fall_or_cross(const double trial[3], double a1, double a2, StrengthCurve& curve,
                         bool plane_strain, int nd, ReturnPoint a, const ReturnPoint& b,
                         ReturnPoint& lo, double& hi, double p[3], double& depls)
{
    ReturnPoint first, last;
    const int fall = first_fall(trial, a1, a2, curve, plane_strain, nd, a, b, first, last);
    if (fall < 0) return reject_no_return;
    if (fall)
        return solve_bracket(trial, a1, a2, curve, plane_strain, nd, first, last, lo, hi, p,
                             depls);
    return cross_break(trial, a1, a2, curve, plane_strain, nd, b, lo, hi, p, depls);
}

// Where the active set changes between start and e, or the trial stops yielding between them,
// g may fall to 0 before the fall that ends at e: whether it does, with [first, last]
// bracketing that fall, or -1 if a return fails.
#pragma acc routine seq
static int earlier_fall(const double trial[3], double a1, double a2, const StrengthCurve& curve,
                        bool plane_strain, int nd, const ReturnPoint& start, const ReturnPoint& e,
                        ReturnPoint& first, ReturnPoint& last)
{
    if (start.set == e.set) {
        ReturnPoint y;
        const int change = yield_change(trial, a1, a2, curve, plane_strain, nd, start, e, y);
        if (change <= 0) return change;
    }
    const int fall = first_fall(trial, a1, a2, curve, plane_strain, nd, start, e, first, last);
    return fall > 0 && !(last.delta < e.delta) ? 0 : fall;
}

// Backward-Euler return for the whole strength law: the plastic increment
// delta whose fixed-strength return at strength(pls+delta) records delta.
// g(delta) = recorded - delta is continuous where the law is, and g(0) > 0
// when the trial yields. The law's linear pieces are searched in order for the
// first increment where g falls to 0; where it stays constant the root is
// exact. Where g jumps below 0 instead, at a rise of the law or inside a
// bracket, the search goes on to where it is back above 0; beyond the last
// break the strength stays and g only falls. Within a piece g can turn where
// friction or dilation change along it and where the active set of the return
// changes, from a face to an edge of the yield surface, say. Where friction or
// dilation change, the trial can also start and stop yielding more than once
// along a piece, and the range of the strength over parts of it finds where
// (yield_change). Of the sets of a yielding trial, the search assumes that each
// holds over one range of delta and that g turns at most once within it, as it
// does with friction and dilation fixed. Where the set changes within the part
// of a piece searched, or g has the same sign at both ends of it, first_fall
// splits it at the changes and looks within each range for a dip below 0 or a
// rise above 0. With law.start_angles
// only cohesion and the tensile cap change with delta, so the search runs only
// for returns the closed forms do not give: corners, the tensile cap and
// breakpoints.
//
// On each piece, secant steps come first: the strength changes little within
// one step, so g is nearly linear and they usually converge in a few
// evaluations. On the first piece they start from guess, the closed-form
// increment with friction and dilation held fixed, when it lies on that
// piece. A root the secant steps do not reach is bracketed and found by
// Illinois regula falsi (solve_bracket).
#pragma acc routine seq
static int search_return(const double trial[3], double a1, double a2,
                         StrengthCurve& curve, bool plane_strain, int nd,
                         double guess, double p[3], double& depls)
{
    if (!finite3(trial)) return reject_nonfinite;
    ReturnPoint lo;           // where the search stands: g > 0, or below 0 past a jump
    ReturnPoint start;        // where it entered the piece: g falls to 0 nowhere before
    if (!evaluate_return(trial, a1, a2, curve, plane_strain, nd, 0, start))
        return reject_no_return;
    double hi = next_break(curve), entered = curve.from;
    bool from_guess = guess > 0 && hi != no_break && guess < hi;
    bool have_lo = !from_guess;
    if (have_lo) lo = start;
    for (int piece = 0; piece < 128; ++piece) {
        if (curve.from != entered) {
            start = lo;
            entered = curve.from;
        }
        ReturnPoint x, b, first, last;
        if (have_lo && lo.g < 0 && !root_found(lo)) {
            // Below 0 past a jump: on from where g rises above 0 and falls back within the
            // piece, if it does, and otherwise across the break.
            if (hi == no_break) return reject_no_increment;
            if (!evaluate_return(trial, a1, a2, curve, plane_strain, nd, hi, b))
                return reject_no_return;
            const int status = fall_or_cross(trial, a1, a2, curve, plane_strain, nd, lo, b, lo,
                                             hi, p, depls);
            if (status != keep_searching) return status;
            continue;
        }
        if (have_lo && (!(lo.g > 0) || hi == no_break)) {
            // Root at lo, or a constant strength beyond lo: recorded(lo) is the root.
            const int status = finish_return(trial, a1, a2, curve, plane_strain, nd, lo, lo, hi,
                                             p, depls);
            if (status != keep_searching) return status;
            continue;
        }

        // Secant steps. The first point is the guess or lo; the next goes to
        // its recorded increment, where g = 0 if the strength did not change.
        ReturnPoint prev;
        bool have_prev = have_lo, bracketed = false, moved = false;
        if (have_lo) prev = lo;
        double next = from_guess ? guess : lo.recorded;
        from_guess = false;
        for (int step = 0; step < 5 && next > (have_lo ? lo.delta : 0) && next < hi; ++step) {
            if (!evaluate_return(trial, a1, a2, curve, plane_strain, nd, next, x))
                return reject_no_return;
            if (root_found(x)) {
                const int earlier = earlier_fall(trial, a1, a2, curve, plane_strain, nd, start, x,
                                                 first, last);
                if (earlier < 0) return reject_no_return;
                const int status = earlier
                    ? solve_bracket(trial, a1, a2, curve, plane_strain, nd, first, last, lo, hi,
                                    p, depls)
                    : finish_return(trial, a1, a2, curve, plane_strain, nd, x, lo, hi, p, depls);
                if (status != keep_searching) return status;
                have_lo = moved = true;
                break;
            }
            if (x.g < 0) {
                b = x;
                bracketed = true;
                break;
            }
            lo = x;
            have_lo = true;
            if (!have_prev) {
                prev = x;
                have_prev = true;
                next = x.recorded;
                continue;
            }
            const double slope = (x.g - prev.g)/(x.delta - prev.delta);
            prev = x;
            if (!(slope < 0)) break;
            next = x.delta - x.g/slope;
        }
        if (moved) continue;
        if (!have_lo) {
            // The guess overshot the root: bracket it from 0.
            lo = start;
            have_lo = true;
            if (!(lo.g > 0)) return accept_return(lo, p, depls);
        }
        if (!bracketed &&
            !evaluate_return(trial, a1, a2, curve, plane_strain, nd, hi, b))
            return reject_no_return;
        // The root in [lo, b] is the first on the piece unless g falls to 0 before it. Where
        // b.g > 0, g may dip below 0 within the piece, and otherwise the root lies past the
        // break.
        int status;
        if (b.g > 0)
            status = fall_or_cross(trial, a1, a2, curve, plane_strain, nd, start, b, lo, hi, p,
                                   depls);
        else {
            const int earlier = earlier_fall(trial, a1, a2, curve, plane_strain, nd, start, b,
                                             first, last);
            if (earlier < 0) return reject_no_return;
            status = solve_bracket(trial, a1, a2, curve, plane_strain, nd, earlier ? first : lo,
                                   earlier ? last : b, lo, hi, p, depls);
        }
        if (status != keep_searching) return status;
    }
    return reject_no_root;
}

// The return search_return finds; one whose search used up return_budget is rejected as
// such, not as one without an admissible return.
#pragma acc routine seq
static int consistent_return(const double trial[3], double a1, double a2,
                             StrengthCurve& curve, bool plane_strain, int nd,
                             double guess, double p[3], double& depls)
{
    const int status = search_return(trial, a1, a2, curve, plane_strain, nd, guess, p, depls);
    return status == reject_no_return && curve.evaluations > return_budget ? reject_budget
                                                                           : status;
}


#pragma acc routine seq
template <typename T>
static void elasto_plastic(double bulkm, double shearm,
                           double amc, double anphi, double anpsi,
                           double hardn, double ten_max,
                           const double* de, double& depls, T s,
                           int &failure_mode,
                           bool has_hydraulic_diffusion,
                           double &dpp, const StrengthLaw& law)
{
    /* Elasto-plasticity (Mohr-Coulomb criterion)
     *
     * failure_mode --
     *  <0: no admissible return, -ReturnRejection
     *   0: no failure
     *   1: tensile failure
     *  10: shear failure
     *  30: return for the whole strength law (see consistent_return)
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
    // eigenvalues (principal stresses); p[2] is unused in 2D
    double p[3] = {0, 0, 0};
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

    for (int i = 0; i < NSTR; ++i)
        if (!std::isfinite(s[i])) {
            failure_mode = -reject_nonfinite;
            return;
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
    const double shear_hardening = 2*std::sqrt(anphi)*hardn;

    bool shear_failure = h < 0;
    const double trial[3] = {p[0], p[1], p[2]};
    StrengthCurve curve = start_curve(law);
    bool closed_form = true;
    // Test the uncut trial first: shear flow may return below the tensile cap.
    if ((ft >= 0 && fs < 0) || shear_failure) {
        double candidate[3], increment;
#ifdef THREED
        const bool sheared = shear_return(trial, candidate, increment, a1, a2, anphi, anpsi,
                                          amc, shear_hardening, ten_max, false);
#else
        const bool sheared = shear_return2d(trial, candidate, increment, a1, a2, anphi, anpsi,
                                            amc, shear_hardening, ten_max);
#endif
        if (sheared) {
            for (int k = 0; k < NDIMS; ++k) p[k] = candidate[k];
            depls = increment;
            shear_failure = true;
        }
        else if (shear_failure)
            closed_form = false;
    }
    double alam;
    if (shear_failure) {
        // shear failure
        failure_mode = 10;
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

        // Biaxial tension (p[1] stays above the cap) and a return past the
        // shear strength need a corner return.
        const double tolerance = return_tolerance(p, amc);
        if (p[1] > p[2]+tolerance ||
            p[0]-anphi*p[2]+amc+shear_hardening*depls < -tolerance)
            closed_form = false;
#else
        /* double depls1 = 0;
         * double depls3 = alam;
         * double deplsm = (depls1 + depls3) / 2;
         * depls = std::sqrt(((depls1-deplsm)*(depls1-deplsm) +
         *                    (depls2-deplsm)*(depls2-deplsm) +
         *                    deplsm*deplsm) / 2);
         */
        depls = std::fabs(alam) * std::sqrt(3. / 8);

        // A return past the shear strength needs the face-cap corner.
        const double tolerance = return_tolerance(p, amc, 2);
        if (p[0] > p[1]+tolerance ||
            p[0]-anphi*p[1]+amc+shear_hardening*depls < -tolerance)
            closed_form = false;
#endif
    }

    if (!closed_form || !closed_form_applies(curve) ||
        !closed_form_consistent(curve, p, depls, NDIMS)) {
        const double guess = closed_form ? depls : 0;
        const int rejection = consistent_return(trial, a1, a2, curve, false, NDIMS, guess,
                                                p, depls);
        if (rejection) {
            failure_mode = -rejection;
            return;
        }
        failure_mode = 30;
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
                             double &dpp, const StrengthLaw& law)
{
    /* Elasto-plasticity (Mohr-Coulomb criterion) */

    /* This function is derived from geoFLAC.
     * The original code in geoFLAC assumes 2D plane strain formulation,
     * i.e. there are 3 principal stresses (PSs) and only 2 principal strains
     * (Strain_yy, Strain_xy, and Strain_yz all must be 0).
     * Here, the code is adopted to pure 2D or 3D plane strain.
     *
     * failure_mode --
     *  <0: no admissible return, -ReturnRejection
     *   0: no failure
     *  10: shear failure (closed-form face or edge return)
     *  30: return for the whole strength law from the uncut trial
     *      (see consistent_return): tension, tension-shear corners and
     *      shear returns the closed form does not cover
     */

    depls = 0;
    failure_mode = 0;

    // elastic trial stress
    double a1 = bulkm + 4. / 3 * shearm;
    double a2 = bulkm - 2. / 3 * shearm;
    const double shear_hardening = 2*std::sqrt(anphi)*hardn;
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
    if (!std::isfinite(sxx) || !std::isfinite(szz) || !std::isfinite(sxz) ||
        !std::isfinite(syy)) {
        failure_mode = -reject_nonfinite;
        return;
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

    const double trial[3] = {p[0], p[1], p[2]};
    const double fs = p[0] - p[2]*anphi + amc;
    const double ft = p[2] - ten_max;
    if (fs >= 0 && ft < 0) {
        // Elastic: keep the trial tensor instead of rebuilding it.
        s[0] = sxx;
        s[1] = szz;
        s[2] = sxz;
        return;
    }

    // A shear face or edge return that finishes below the tensile cap is the
    // closed form of the uncut-trial return. Tension always goes to
    // consistent_return, with associated flow on the cap and its plastic
    // strain in history, as in 3D.
    StrengthCurve curve = start_curve(law);
    bool closed_form = false;
    if (fs < 0) {
        failure_mode = 10;
        closed_form = shear_return(trial, p, depls, a1, a2, anphi, anpsi, amc,
                                   shear_hardening, ten_max, true);
    }
    if (!closed_form || !closed_form_applies(curve) || !closed_form_consistent(curve, p, depls)) {
        const double guess = closed_form ? depls : 0;
        const int rejection = consistent_return(trial, a1, a2, curve, true, 3, guess, p, depls);
        if (rejection) {
            failure_mode = -rejection;
            return;
        }
        failure_mode = 30;
    }

    // rotate the principal stresses back to global axes
    rotate_stress2d(p, n1, n2, n3, cos2t, sin2t, s, syy);
}

static const char* rejection_text(int reason)
{
    switch (reason) {
    case reject_nonfinite:
        return "non-finite trial stress";
    case reject_no_return:
        return "no set of shear faces and tensile caps gives an admissible return";
    case reject_no_root:
        return "the strength law has more linear pieces than the return searches";
    case reject_no_increment:
        return "no plastic increment returns the stress to the strength of the plastic strain "
               "it stores, where the strength law jumps";
    case reject_budget:
        return "the search for the return reached its limit of evaluations";
    default:
        return "unknown";
    }
}

static void report_rejected_returns(const Variables& var, int count,
                                    double first, int reasons)
{
    const long element = long(first) / reject_reason_count;
    const int reason = int(long(first) % reject_reason_count);
    std::cerr << "Error: no admissible Mohr-Coulomb return in " << count
              << " elements at step " << var.steps << ".\n"
              << "  First element " << element << ": " << rejection_text(reason) << ".\n";
    for (int r = 1; r < reject_reason_count; ++r)
        if (r != reason && (reasons & (1 << r)))
            std::cerr << "  Also: " << rejection_text(r) << ".\n";
    die(EXIT_RUNTIME_NONCONVERGENCE,
        reasons & (1 << reject_nonfinite) ? "Non-finite stress: check the time step." :
        reasons & (1 << reject_no_increment) ?
            "Spread the change of the strength law over a range of plastic strain." :
            "Please report this stress state.");
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
    // Only option 0 projects the centroid velocity; option 1 derives its rate
    // from the strain-rate invariant.
    const bool needs_slip_rate = (param.mat.rheol_type == MatProps::rh_ep_rsf)
                              || (param.mat.rheol_type == MatProps::rh_evp_rsf);
    const bool needs_projected_velocity =
        needs_slip_rate &&
        param.control.rsf_slip_rate_projection_option ==
            rsf_slip_rate_projection_maximum_shear;
    const bool start_angles = param.control.iterative_return_mapping_option == 0;

    // Count, lowest element*reject_reason_count+reason, and the set of reasons.
    int rejected_returns = 0;
    double first_rejection = std::numeric_limits<double>::max();
    int rejection_reasons = 0;
#ifndef ACC
    #pragma omp parallel for default(none) shared(param, var, ppressure, dppressure, \
        vel, stress, stressyy, dpressure, viscosity, strain, plstrain, delta_plstrain, \
        strain_rate, dyn_fric_coeff, state_variable) \
        firstprivate(has_hydraulic_diffusion, needs_slip_rate, trial_reference_geometry, initial_equilibrium, \
                     needs_projected_velocity, start_angles) reduction(+:rejected_returns) \
        reduction(min:first_rejection) reduction(|:rejection_reasons)
#endif
    #pragma acc parallel loop gang vector async reduction(+:rejected_returns) reduction(min:first_rejection) reduction(|:rejection_reasons) // TODO: ACC: CPU and GPU results are differet because of using 3x3 in elasto_plastic
    for (int e = 0; e < var.nelem; e++) {
        ConstConnAccessor conn = (*var.connectivity)[e];
        int rejection = 0;

        // Centroid interpolations, each gated on the branch that reads it: ungated they cost
        // 2 + NDIMS nodal gathers per element per step in every run. The `ppressure` LEVEL
        // is not interpolated at all; restore it here if a rheology comes to need it.
        double dpp = 0.0;
        if (has_hydraulic_diffusion) {
            #pragma acc loop seq
            for (int j = 0; j < NODES_PER_ELEM; ++j)
                dpp += dppressure[conn[j]] / double(NODES_PER_ELEM);
            // Biot-weighted pore-pressure INCREMENT, not a rate: update_pore_pressure()
            // folded dt in, and stores dppressure with the sign OPPOSITE to its change
            // to ppressure. elastic_effective() adds this to the diagonal components.
            dpp *= var.mat->alpha_biot(e);
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

        // // Calculate the center of the element
        // const int *conn = (*var.connectivity)[e];
        // double center_x = 0., center_z = 0., T = 0.;
        // for (int i=0; i<NODES_PER_ELEM; ++i){
        //     center_x += (*var.coord)[conn[i]][0];
        // center_z += (*var.coord)[conn[i]][NDIMS-1];
        //     T += (*var.temperature)[conn[i]];
        // }
        // T /= NODES_PER_ELEM;
        // center_x /= NODES_PER_ELEM;
        // center_z /= NODES_PER_ELEM;
        // // Find the most abundant marker mattype in this element
        // int_vec &a = (*var.elemmarkers)[e];
        // int material = std::distance(a.begin(), std::max_element(a.begin(), a.end()));
        
        // anti-mesh locking correction on strain rate
        if(1){
            double div = trace(edot);
            //double div2 = ((*var.volume)[e] / (*var.volume_old)[e] - 1) / var.dt;
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
                // The piece of the weakening law at pls gives the strength
                // here and is reused by the plastic returns.
                MatProps::StrengthPiece piece;
                var.mat->strength_here(e, plstrain[e], piece);
                const double amc = piece.amc, anphi = piece.anphi, anpsi = piece.anpsi;
                const double hardn = piece.dcohesion, ten_max = piece.ten_max;
                const StrengthLaw law = {var.mat, e, plstrain[e],
                                         amc, anphi, anpsi, ten_max, false, 0, 0, 0, 0, false,
                                         &piece, start_angles};
                int failure_mode;
                if (var.mat->is_plane_strain) {
                    elasto_plastic2d(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                     de, depls, s, syy, failure_mode, 
                                     has_hydraulic_diffusion, dpp, law);
                }
                else {
                    elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                   de, depls, s, failure_mode, 
                                   has_hydraulic_diffusion, dpp, law);
                }
                if (failure_mode < 0) rejection = -failure_mode;
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

                // The piece of the weakening law at pls gives the strength
                // here and is reused by the plastic returns.
                MatProps::StrengthPiece piece;
                var.mat->strength_here(e, plstrain[e], piece);
                const double amc = piece.amc, anphi = piece.anphi, anpsi = piece.anpsi;
                const double hardn = piece.dcohesion, ten_max = piece.ten_max;
                const StrengthLaw law = {var.mat, e, plstrain[e],
                                         amc, anphi, anpsi, ten_max, false, 0, 0, 0, 0, false,
                                         &piece, start_angles};
                // stress due to elasto-plastic rheology
                double sp[NSTR], spyy;
                #pragma acc loop seq
                for (int i=0; i<NSTR; ++i) sp[i] = s[i];
                int failure_mode;
                if (var.mat->is_plane_strain) {
                    spyy = syy;
                    elasto_plastic2d(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                     de, depls, sp, spyy, failure_mode, 
                                     has_hydraulic_diffusion, dpp, law);
                }
                else {
                    elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                   de, depls, sp, failure_mode, 
                                   has_hydraulic_diffusion, dpp, law);
                }
                double spII = plane_trial ? second_invariant2(sp, spyy) : second_invariant2(sp);

                // use the smaller as the final stress
                if (svII < spII) {
                    #pragma acc loop seq
                    for (int i=0; i<NSTR; ++i) s[i] = sv[i];
                    if (plane_trial) syy = svyy;
                }
                else {
                    #pragma acc loop seq
                    for (int i=0; i<NSTR; ++i) s[i] = sp[i];
                    // Only a stored plastic stress must be admissible.
                    if (failure_mode < 0) rejection = -failure_mode;
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
                const StrengthLaw law = {var.mat, e, plstrain[e],
                                         amc, anphi, anpsi, ten_max, true, slip_rate,
                                         state_variable[e], dyn_fric_coeff[e],
                                         param.mat.state_var_model, initial_equilibrium, nullptr,
                                         start_angles};
                int failure_mode;
                if (var.mat->is_plane_strain) {
                    elasto_plastic2d(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                     de, depls, s, syy, failure_mode, 
                                     has_hydraulic_diffusion, dpp, law);
                }
                else {
                    elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                   de, depls, s, failure_mode, 
                                   has_hydraulic_diffusion, dpp, law);
                }
                if (failure_mode < 0) rejection = -failure_mode;
                plstrain[e] += depls;
                delta_plstrain[e] = depls;
                // The iterative return used the friction of the stored
                // history; store that coefficient. The state was already
                // advanced this step, and frozen initial-equilibrium friction
                // stays as stored. Option 0 used the stored coefficient.
                if (depls > 0 && !initial_equilibrium && !start_angles)
                    var.mat->rsf_friction_from_state(e, plstrain[e], slip_rate,
                                                     state_variable[e], dyn_fric_coeff[e],
                                                     param.mat.state_var_model);
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
                const StrengthLaw law = {var.mat, e, plstrain[e],
                                         amc, anphi, anpsi, ten_max, true, slip_rate,
                                         state_variable[e], dyn_fric_coeff[e],
                                         param.mat.state_var_model, initial_equilibrium, nullptr,
                                         start_angles};
                // stress due to elasto-plastic rheology
                double sp[NSTR], spyy;
                #pragma acc loop seq
                for (int i=0; i<NSTR; ++i) sp[i] = s[i];
                int failure_mode;
                if (var.mat->is_plane_strain) {
                    spyy = syy;
                    elasto_plastic2d(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                     de, depls, sp, spyy, failure_mode, 
                                     has_hydraulic_diffusion, dpp, law);
                }
                else {
                    elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                   de, depls, sp, failure_mode, 
                                   has_hydraulic_diffusion, dpp, law);
                }
                double spII = plane_trial ? second_invariant2(sp, spyy) : second_invariant2(sp);

                // use the smaller as the final stress
                if (svII < spII) {
                    #pragma acc loop seq
                    for (int i=0; i<NSTR; ++i) s[i] = sv[i];
                    if (plane_trial) syy = svyy;
                }
                else {
                    #pragma acc loop seq
                    for (int i=0; i<NSTR; ++i) s[i] = sp[i];
                    // Only a stored plastic stress must be admissible.
                    if (failure_mode < 0) rejection = -failure_mode;
                    plstrain[e] += depls;
                    delta_plstrain[e] = depls;
                    // The iterative return used the friction of the stored
                    // history; store that coefficient. The state was already
                    // advanced this step, and frozen initial-equilibrium friction
                    // stays as stored. Option 0 used the stored coefficient.
                    if (depls > 0 && !initial_equilibrium && !start_angles)
                        var.mat->rsf_friction_from_state(e, plstrain[e], slip_rate,
                                                         state_variable[e], dyn_fric_coeff[e],
                                                         param.mat.state_var_model);
                    if (var.mat->is_plane_strain) syy = spyy;
                }
            }
            break;
        default:
//            std::cerr << "Error: unknown rheology type: " << rheol_type << "\n";
//            std::exit(1);
            break;
        }
        if (param.control.is_using_mixed_stress) {
            const double new_trace = trial_reference_geometry && var.mat->is_plane_strain
                ? trace(s)+syy : trace(s);
            dpressure[e] = new_trace-old_s;
        }
        // Updated without a branch or a select: NVHPC loses conditional updates of
        // these reductions on the GPU, and the run would go on without exit 54 or
        // report no element. The key is an exact integer in a double: NVHPC's GPU
        // min reduction of a long long is wrong.
        rejected_returns += rejection != 0;
        rejection_reasons |= (rejection != 0) << rejection;
        first_rejection = std::min(first_rejection, double(e)*int(reject_reason_count) + rejection
                                   + (rejection == 0)*std::numeric_limits<double>::max());
        // std::cerr << "stress " << e << ": ";
        // print(std::cerr, s, NSTR);
        // std::cerr << '\n';
    }
    #pragma acc wait
    if (rejected_returns != 0)
        report_rejected_returns(var, rejected_returns, first_rejection, rejection_reasons);
#ifdef NPROF
    nvtxRangePop();
#endif
}

void update_old_mean_stress(const Param& param, const Variables& var, tensor_t& stress,
                   double_vec& old_mean_stress)
{
#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif

#ifndef ACC
    #pragma omp parallel for default(none)                           \
        shared(param, var, stress, old_mean_stress)
#endif
    #pragma acc parallel loop gang vector async
    for (int e=0; e<var.nelem; ++e) {
        ConstTensorAccessor s = stress[e];
        old_mean_stress[e] =trace(s)/NDIMS;
    }
#ifdef NPROF
    nvtxRangePop();
#endif
}
