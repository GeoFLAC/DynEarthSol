#include <cmath>
#include <iostream>
#ifdef USE_NPROF
#include <nvToolsExt.h>
#endif

#include "3x3-C/dsyevh3.h"

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "rheology.hpp"
#include "utils.hpp"

#pragma acc routine seq
static void principal_stresses3(const double* s, double p[3], double v[3][3])
{
    /* s is a flattened stress vector, with the components {XX, YY, ZZ, XY, XZ, YZ}.
     * Returns the eigenvalues p and eignvectors v.
     * The eigenvalues are ordered such that p[0] <= p[1] <= p[2].
     *  - Maximum shear stress direction 'shear_dir'.
     */

    // unflatten s to a 3x3 tensor, only the upper part is needed.
    double a[3][3];
    a[0][0] = s[0];
    a[1][1] = s[1];
    a[2][2] = s[2];
    a[0][1] = s[3];
    a[0][2] = s[4];
    a[1][2] = s[5];

    dsyevh3(a, v, p);

    // reorder p and v
    if (p[0] > p[1]) {
        double tmp, b[3];
        tmp = p[0];
        p[0] = p[1];
        p[1] = tmp;
        for (int i=0; i<3; ++i)
            b[i] = v[i][0];
        for (int i=0; i<3; ++i)
            v[i][0] = v[i][1];
        for (int i=0; i<3; ++i)
            v[i][1] = b[i];
    }
    if (p[1] > p[2]) {
        double tmp, b[3];
        tmp = p[1];
        p[1] = p[2];
        p[2] = tmp;
        for (int i=0; i<3; ++i)
            b[i] = v[i][1];
        for (int i=0; i<3; ++i)
            v[i][1] = v[i][2];
        for (int i=0; i<3; ++i)
            v[i][2] = b[i];
    }
    if (p[0] > p[1]) {
        double tmp, b[3];
        tmp = p[0];
        p[0] = p[1];
        p[1] = tmp;
        for (int i=0; i<3; ++i)
            b[i] = v[i][0];
        for (int i=0; i<3; ++i)
            v[i][0] = v[i][1];
        for (int i=0; i<3; ++i)
            v[i][1] = b[i];
    }
}

#pragma acc routine seq
static void principal_stresses2(const double* s, double p[2],
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
static void compute_slip_rate2(const double* s, double& vx, double& vz, double& slip_rate)
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
static void compute_slip_rate3(const double* s, double vx, double vy, double vz, double& slip_rate)
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

#pragma acc routine seq
static void elastic(double bulkm, double shearm, const double* de, double* s)
{
    /* increment the stress s according to the incremental strain de */
    double lambda = bulkm - 2. /3 * shearm;
    double dev = trace(de);

    for (int i=0; i<NDIMS; ++i)
        s[i] += 2 * shearm * de[i] + lambda * dev;
    for (int i=NDIMS; i<NSTR; ++i)
        s[i] += 2 * shearm * de[i];
}

#pragma acc routine seq
static void elastic_effective(double bulkm, double shearm, const double* de, double* s,  double &dpp)
{
    /* increment the stress s according to the incremental strain de */
    double lambda = bulkm - 2. /3 * shearm;
    double dev = trace(de);

    for (int i=0; i<NDIMS; ++i)
        s[i] += 2 * shearm * de[i] + lambda * dev + dpp;

    for (int i=NDIMS; i<NSTR; ++i)
        s[i] += 2 * shearm * de[i];
}

#pragma acc routine seq
static void maxwell(double bulkm, double shearm, double viscosity, double dt,
                    double dv, const double* de, double* s)
{
    // non-dimensional parameter: dt/ relaxation time
    double tmp = 0.5 * dt * shearm / viscosity;
    double f1 = 1 - tmp;
    double f2 = 1 / (1  + tmp);

    double dev = trace(de) / NDIMS;
    double s0 = trace(s) / NDIMS;

    // convert back to total stress
    for (int i=0; i<NDIMS; ++i)
        s[i] = ((s[i] - s0) * f1 + 2 * shearm * (de[i] - dev)) * f2 + s0 + bulkm * dv;
    for (int i=NDIMS; i<NSTR; ++i)
        s[i] = (s[i] * f1 + 2 * shearm * de[i]) * f2;
}


#pragma acc routine seq
static void viscous(double bulkm, double viscosity, double total_dv,
                    const double* edot, double* s)
{
    /* Viscous Model + incompressibility enforced by bulk modulus */

    double dev = trace(edot) / NDIMS;

    for (int i=0; i<NDIMS; ++i)
        s[i] = 2 * viscosity * (edot[i] - dev) + bulkm * total_dv;
    for (int i=NDIMS; i<NSTR; ++i)
        s[i] = 2 * viscosity * edot[i];
}

#pragma acc routine seq
static void elasto_plastic(double bulkm, double shearm,
                           double amc, double anphi, double anpsi,
                           double hardn, double ten_max,
                           const double* de, double& depls, double* s,
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
    if (h < 0) {
        // shear failure
        failure_mode = 10;

        alam = fs / (a1 - a2*anpsi + a1*anphi*anpsi - a2*anphi + 2*std::sqrt(anphi)*hardn);
        p[0] -= alam * (a1 - a2 * anpsi);
#ifdef THREED
        p[1] -= alam * (a2 - a2 * anpsi);
#endif
        p[NDIMS-1] -= alam * (a2 - a1 * anpsi);

        // 2nd invariant of plastic strain
#ifdef THREED
        /* // plastic strain in the principal directions, depls2 is always 0
         * double depls1 = alam;
         * double depls3 = -alam * anpsi;
         * double deplsm = (depls1 + depls3) / 3;
         * depls = std::sqrt(((depls1-deplsm)*(depls1-deplsm) +
         *                    (-deplsm)*(-deplsm) +
         *                    (depls3-deplsm)*(depls3-deplsm) +
         *                    deplsm*deplsm) / 2);
         */
        // the equations above can be reduce to:
        depls = std::fabs(alam) * std::sqrt((7 + 4*anpsi + 7*anpsi*anpsi) / 18);
#else
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
static void elasto_plastic2d(double bulkm, double shearm,
                             double amc, double anphi, double anpsi,
                             double hardn, double ten_max,
                             const double* de, double& depls,
                             double* s, double &syy,
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
    const double alams = fs / (a1 - a2*anpsi + a1*anphi*anpsi - a2*anphi + hardn);
    p[0] -= alams * (a1 - a2 * anpsi);
    p[1] -= alams * (a2 - a2 * anpsi);
    p[2] -= alams * (a2 - a1 * anpsi);

    // 2nd invariant of plastic strain
    depls = 0.5 * std::fabs(alams + alams * anpsi);

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
    {
        double dc2 = (p[n1] - p[n2]) * cos2t;
        double dss = p[n1] + p[n2];
        s[0] = 0.5 * (dss + dc2);
        s[1] = 0.5 * (dss - dc2);
        s[2] = 0.5 * (p[n1] - p[n2]) * sin2t;
        syy = p[n3];
    }
}

void update_stress(const Param& param, Variables& var, tensor_t& stress,
                   double_vec& stressyy, double_vec& dpressure, double_vec& viscosity,
                   tensor_t& strain, double_vec& plstrain,
                   double_vec& delta_plstrain, tensor_t& strain_rate,
                   double_vec& ppressure, double_vec& dppressure, array_t& vel)

{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    // Interpolate pore pressure from nodes to element level
    double_vec ppressure_element(var.nelem, 0.0);
    double_vec dppressure_element(var.nelem, 0.0);
    // Define element velocity arrays
    double_vec velocity_x_element(var.nelem, 0.0);
    double_vec velocity_y_element(var.nelem, 0.0);
    double_vec velocity_z_element(var.nelem, 0.0); // Used only for 3D

    #pragma omp parallel for default(none) shared(var, ppressure, ppressure_element, dppressure, dppressure_element, \
        vel, velocity_x_element, velocity_y_element, velocity_z_element)
    #pragma acc parallel loop
    for (int e = 0; e < var.nelem; e++) {
        const int *conn = (*var.connectivity)[e];
        double pp_element = 0.0;
        double dpp_element = 0.0;

        const array_t& vel = *var.vel;
        double vx_element = 0.0, vy_element = 0.0, vz_element = 0.0;

        for (int j = 0; j < NODES_PER_ELEM; ++j) {
        #ifdef THREED
            pp_element += ppressure[conn[j]] / 4.0; // the centroid shape functions are 1/4 for each node in 3D
            dpp_element += dppressure[conn[j]] / 4.0; // the centroid shape functions are 1/4 for each node in 3D
            vx_element += vel[conn[j]][0] / 4.0; // the centroid shape functions are 1/4 for each node in 3D
            vy_element += vel[conn[j]][1] / 4.0; // the centroid shape functions are 1/4 for each node in 3D
            vz_element += vel[conn[j]][2] / 4.0; // the centroid shape functions are 1/4 for each node in 3D
        #else
            pp_element += ppressure[conn[j]] / 3.0; // the centroid shape functions are 1/3 for each node in 2D
            dpp_element += dppressure[conn[j]] / 3.0; // the centroid shape functions are 1/3 for each node in 2D
            vx_element += vel[conn[j]][0] / 3.0; // the centroid shape functions are 1/3 for each node in 2D
            vy_element += vel[conn[j]][1] / 3.0; // the centroid shape functions are 1/3 for each node in 2D

        #endif
        }
        ppressure_element[e] = pp_element;  // Store element-level pore pressure 
        dppressure_element[e] = dpp_element;  // Store element-level pore pressure rate
        velocity_x_element[e] = vx_element;  // Store element-level velocity in x direction
        velocity_y_element[e] = vy_element;  // Store element-level velocity in y direction
        #ifdef THREED
        velocity_z_element[e] = vz_element;  // Store element-level velocity in z direction
        #endif
    }

    #pragma omp parallel for default(none)                           \
        shared(param, var, stress, stressyy, dpressure, viscosity, strain, plstrain, delta_plstrain, \
               strain_rate, ppressure_element, dppressure_element, std::cerr, \
               velocity_x_element, velocity_y_element, velocity_z_element)
    #pragma acc parallel loop
    for (int e=0; e<var.nelem; ++e) {
        // stress, strain and strain_rate of this element
        double* s = stress[e];
        double& syy = stressyy[e];
        double* es = strain[e];
        double* edot = strain_rate[e];
    	double old_s = trace(s);

        // Calculate effective pore pressure using Biot's coefficient
        double alpha_b = var.mat->alpha_biot(e); // Biot coefficient
        double pp = ppressure_element[e];        // Use element-level interpolated pore pressure
        double dpp = dppressure_element[e];        // Use element-level interpolated pore pressure rate

        double vx = velocity_x_element[e];        // Use element-level interpolated velocity in x direction
        double vy = velocity_y_element[e];        // Use element-level interpolated velocity in y direction
        double vz = velocity_z_element[e];        // Use element-level interpolated velocity in z direction

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
        
        if (param.control.has_hydraulic_diffusion) {
            pp = alpha_b * pp; // Apply Biot coefficient to pore pressure
            dpp = alpha_b * dpp; // Apply Biot coefficient to pore pressure rate
        }

        // anti-mesh locking correction on strain rate
        if(1){
            double div = trace(edot);
            //double div2 = ((*var.volume)[e] / (*var.volume_old)[e] - 1) / var.dt;
            for (int i=0; i<NDIMS; ++i) {
                edot[i] += ((*var.edvoldt)[e] - div) / NDIMS;  // XXX: should NDIMS -> 3 in plane strain?
            }
        }

        // update strain with strain rate
        for (int i=0; i<NSTR; ++i) {
            es[i] += edot[i] * var.dt;
        }

        // modified strain increment
        double de[NSTR];
        for (int i=0; i<NSTR; ++i) {
            de[i] = edot[i] * var.dt;
        }

        switch (param.mat.rheol_type) {
        case MatProps::rh_elastic:
            {
                double bulkm = var.mat->bulkm(e);
                double shearm = var.mat->shearm(e);
                if (param.control.has_hydraulic_diffusion)
                {
                    elastic_effective(bulkm, shearm, de, s, dpp);
                }
                else
                {
                    elastic(bulkm, shearm, de, s);
                }
            }
            break;
        case MatProps::rh_viscous:
            {
                double bulkm = var.mat->bulkm(e);
                viscosity[e] = var.mat->visc(e);
                double total_dv = trace(es);
                viscous(bulkm, viscosity[e], total_dv, edot, s);
            }
            break;
        case MatProps::rh_maxwell:
            {
                double bulkm = var.mat->bulkm(e);
                double shearm = var.mat->shearm(e);
                viscosity[e] = var.mat->visc(e);
                double dv = (*var.volume)[e] / (*var.volume_old)[e] - 1;
                maxwell(bulkm, shearm, viscosity[e], var.dt, dv, de, s);
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
                int failure_mode;
                if (var.mat->is_plane_strain) {
                    elasto_plastic2d(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                     de, depls, s, syy, failure_mode, 
                                     param.control.has_hydraulic_diffusion, dpp);
                }
                else {
                    elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                   de, depls, s, failure_mode, 
                                   param.control.has_hydraulic_diffusion, dpp);
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
                double dv = (*var.volume)[e] / (*var.volume_old)[e] - 1;
                // stress due to maxwell rheology
                double sv[NSTR];
                for (int i=0; i<NSTR; ++i) sv[i] = s[i];
                maxwell(bulkm, shearm, viscosity[e], var.dt, dv, de, sv);
                double svII = second_invariant2(sv);

                double amc, anphi, anpsi, hardn, ten_max;
                var.mat->plastic_props(e, plstrain[e],
                                       amc, anphi, anpsi, hardn, ten_max);
                // stress due to elasto-plastic rheology
                double sp[NSTR], spyy;
                for (int i=0; i<NSTR; ++i) sp[i] = s[i];
                int failure_mode;
                if (var.mat->is_plane_strain) {
                    spyy = syy;
                    elasto_plastic2d(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                     de, depls, s, syy, failure_mode, 
                                     param.control.has_hydraulic_diffusion, dpp);
                }
                else {
                    elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                   de, depls, s, failure_mode, 
                                   param.control.has_hydraulic_diffusion, dpp);
                }
                double spII = second_invariant2(sp);

                // use the smaller as the final stress
                if (svII < spII)
                    for (int i=0; i<NSTR; ++i) s[i] = sv[i];
                else {
                    for (int i=0; i<NSTR; ++i) s[i] = sp[i];
                    plstrain[e] += depls;
                    delta_plstrain[e] = depls;
                    syy = spyy;
                }
            }
            case MatProps::rh_ep_rsf: // rate-and-state frition model
            {
                double depls = 0;
                double bulkm = var.mat->bulkm(e);
                double shearm = var.mat->shearm(e);

                // calculate the shear direction in maximum shear stress
                double slip_rate;
                
#ifdef THREED
                compute_slip_rate3(s, vx, vy, vz, slip_rate); // Calculate slip vector magnitude in 3D
#else
                compute_slip_rate2(s, vx, vy, slip_rate); // Calculate slip vector magnitude in 2D
#endif
        
                double amc, anphi, anpsi, hardn, ten_max;
                var.mat->plastic_props_rsf(e, plstrain[e],
                                       amc, anphi, anpsi, hardn, ten_max, slip_rate);
                int failure_mode;
                if (var.mat->is_plane_strain) {
                    elasto_plastic2d(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                     de, depls, s, syy, failure_mode, 
                                     param.control.has_hydraulic_diffusion, dpp);
                }
                else {
                    elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                   de, depls, s, failure_mode, 
                                   param.control.has_hydraulic_diffusion, dpp);
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
                double dv = (*var.volume)[e] / (*var.volume_old)[e] - 1;
                // stress due to maxwell rheology
                double sv[NSTR];
                for (int i=0; i<NSTR; ++i) sv[i] = s[i];
                maxwell(bulkm, shearm, viscosity[e], var.dt, dv, de, sv);
                double svII = second_invariant2(sv);

                 // calculate the shear direction in maximum shear stress
                double slip_rate;
                
#ifdef THREED
                compute_slip_rate3(s, vx, vy, vz, slip_rate); // Calculate slip vector magnitude in 3D
#else
                compute_slip_rate2(s, vx, vy, slip_rate); // Calculate slip vector magnitude in 2D
#endif

                double amc, anphi, anpsi, hardn, ten_max;
                var.mat->plastic_props_rsf(e, plstrain[e],
                                       amc, anphi, anpsi, hardn, ten_max, slip_rate);
                // stress due to elasto-plastic rheology
                double sp[NSTR], spyy;
                for (int i=0; i<NSTR; ++i) sp[i] = s[i];
                int failure_mode;
                if (var.mat->is_plane_strain) {
                    spyy = syy;
                    elasto_plastic2d(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                     de, depls, s, syy, failure_mode, 
                                     param.control.has_hydraulic_diffusion, dpp);
                }
                else {
                    elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                   de, depls, s, failure_mode, 
                                   param.control.has_hydraulic_diffusion, dpp);
                }
                double spII = second_invariant2(sp);

                // use the smaller as the final stress
                if (svII < spII)
                    for (int i=0; i<NSTR; ++i) s[i] = sv[i];
                else {
                    for (int i=0; i<NSTR; ++i) s[i] = sp[i];
                    plstrain[e] += depls;
                    delta_plstrain[e] = depls;
                    syy = spyy;
                }
            }
            break;
        default:
//            std::cerr << "Error: unknown rheology type: " << rheol_type << "\n";
//            std::exit(1);
            break;
        }
        if (param.control.is_using_mixed_stress)
            dpressure[e] = trace(s) - old_s;
        // std::cerr << "stress " << e << ": ";
        // print(std::cerr, s, NSTR);
        // std::cerr << '\n';
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void update_old_mean_stress(const Param& param, const Variables& var, tensor_t& stress,
                   double_vec& old_mean_stress)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    #pragma omp parallel for default(none)                           \
        shared(param, var, stress, old_mean_stress)
    #pragma acc parallel loop
    for (int e=0; e<var.nelem; ++e) {
        double* s = stress[e];
        old_mean_stress[e] =trace(s)/NDIMS;
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}