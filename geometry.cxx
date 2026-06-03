#include <cmath>
#include <limits>
#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "utils.hpp"
#include "geometry.hpp"
#include "bc.hpp"
#include "mesh.hpp"
#include "output.hpp"

/* Given two points, returns the distance^2 */
template <typename T>
double dist2(T a, T b)
{
    double sum = 0;
    for (int i=0; i<NDIMS; ++i) {
        double d = b[i] - a[i];
        sum += d * d;
    }
    return sum;
}

template
double dist2<ConstArrayAccessor>(ConstArrayAccessor a, ConstArrayAccessor b);
template
double dist2<double*>(double *a, double *b);


/* Given four 3D points, returns the (signed) volume of the enclosed
   tetrahedron */
#pragma acc routine seq
static double tetrahedron_volume(ConstArrayAccessor d0,
                                 ConstArrayAccessor d1,
                                 ConstArrayAccessor d2,
                                 ConstArrayAccessor d3)
{
    double x01 = d0[0] - d1[0];
    double x12 = d1[0] - d2[0];
    double x23 = d2[0] - d3[0];

    double y01 = d0[1] - d1[1];
    double y12 = d1[1] - d2[1];
    double y23 = d2[1] - d3[1];

    double z01 = d0[2] - d1[2];
    double z12 = d1[2] - d2[2];
    double z23 = d2[2] - d3[2];

    return (x01*(y23*z12 - y12*z23) +
            x12*(y01*z23 - y23*z01) +
            x23*(y12*z01 - y01*z12)) / 6;
}

#pragma acc routine seq
static double triangle_area2d(ConstArrayAccessor a,
                            ConstArrayAccessor b,
                            ConstArrayAccessor c)
{
    double ab0, ab1, ac0, ac1;

    // ab: vector from a to b
    ab0 = b[0] - a[0];
    ab1 = b[1] - a[1];
    // ac: vector from a to c
    ac0 = c[0] - a[0];
    ac1 = c[1] - a[1];

    return std::fabs(ab0*ac1 - ab1*ac0) / 2;
}

/* Given two points, returns the area of the enclosed triangle */
#pragma acc routine seq
static double triangle_area(ConstArrayAccessor a,
                            ConstArrayAccessor b,
                            ConstArrayAccessor c)
{
    double ab0, ab1, ac0, ac1;

    // ab: vector from a to b
    ab0 = b[0] - a[0];
    ab1 = b[1] - a[1];
    // ac: vector from a to c
    ac0 = c[0] - a[0];
    ac1 = c[1] - a[1];

#ifndef THREED
    // area = norm(cross product of ab and ac) / 2
    return std::fabs(ab0*ac1 - ab1*ac0) / 2;
#else
    double ab2, ac2;
    ab2 = b[2] - a[2];
    ac2 = c[2] - a[2];

    // vector components of ab x ac
    double d0, d1, d2;
    d0 = ab1*ac2 - ab2*ac1;
    d1 = ab2*ac0 - ab0*ac2;
    d2 = ab0*ac1 - ab1*ac0;

    // area = norm(cross product of ab and ac) / 2
    return std::sqrt(d0*d0 + d1*d1 + d2*d2) / 2;
#endif
}

double compute_area_facet(ConstArrayIndirectAccessor coord)
{
    ConstArrayAccessor a = coord[0];
    ConstArrayAccessor b = coord[1];

#ifndef THREED
    return std::fabs(a[0]-b[0]);
#else
    ConstArrayAccessor c = coord[2];

    return triangle_area2d(a, b, c);
#endif
}

double compute_volume(ConstArrayIndirectAccessor coord)
{
    ConstArrayAccessor a = coord[0];
    ConstArrayAccessor b = coord[1];
    ConstArrayAccessor c = coord[2];
#ifdef THREED
    ConstArrayAccessor d = coord[3];
    return tetrahedron_volume(a, b, c, d);
#else
    return triangle_area(a, b, c);
#endif
}

void compute_volume(const array_t &coord, const conn_t &connectivity,
                    double_vec &volume)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif

#ifndef ACC
    #pragma omp parallel for default(none)      \
        shared(coord, connectivity, volume)
#endif
    #pragma acc parallel loop gang vector async
    for (int e=0; e<volume.size(); ++e) {
        int n0 = connectivity[e][0];
        int n1 = connectivity[e][1];
        int n2 = connectivity[e][2];

        ConstArrayAccessor a = coord[n0];
        ConstArrayAccessor b = coord[n1];
        ConstArrayAccessor c = coord[n2];

#ifdef THREED
        int n3 = connectivity[e][3];
        ConstArrayAccessor d = coord[n3];
        volume[e] = tetrahedron_volume(a, b, c, d);
#else
        volume[e] = triangle_area(a, b, c);
#endif
    }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

void compute_volume(const Variables &var,
                    double_vec &volume)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif

#ifndef ACC
    #pragma omp parallel for default(none) shared(var, volume)
#endif
    #pragma acc parallel loop gang vector async
    for (int e=0; e<var.nelem; ++e) {
        int n0 = (*var.connectivity)[e][0];
        int n1 = (*var.connectivity)[e][1];
        int n2 = (*var.connectivity)[e][2];

        ConstArrayAccessor a = (*var.coord)[n0];
        ConstArrayAccessor b = (*var.coord)[n1];
        ConstArrayAccessor c = (*var.coord)[n2];

#ifdef THREED
        int n3 = (*var.connectivity)[e][3];
        ConstArrayAccessor d = (*var.coord)[n3];
        volume[e] = tetrahedron_volume(a, b, c, d);
#else
        volume[e] = triangle_area(a, b, c);
#endif
    }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

void compute_dvoldt(const Variables &var, double_vec &dvoldt, double_vec &etmp)
{
#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    /* dvoldt is the volumetric strain rate, weighted by the element volume,
     * lumped onto the nodes.
     */
//    std::fill_n(dvoldt.begin(), var.nnode, 0);

#ifndef ACC
    #pragma omp parallel for default(none)      \
        shared(var, etmp)
#endif
    #pragma acc parallel loop gang vector async
    for (int e=0;e<var.nelem;e++) {
        ConstTensorAccessor strain_rate = (*var.strain_rate)[e];
        // TODO: try another definition:
        // dj = (volume[e] - volume_old[e]) / volume_old[e] / dt
        double dj = trace(strain_rate);
        etmp[e] = dj * (*var.volume)[e];
    }

#ifndef ACC
    #pragma omp parallel for default(none)      \
        shared(var,dvoldt,etmp)
#endif
    #pragma acc parallel loop gang vector async
    for (int n=0;n<var.nnode;n++) {
        dvoldt[n] = 0.;
        for( auto e = (*var.support)[n].begin(); e < (*var.support)[n].end(); ++e)
	        dvoldt[n] += etmp[*e];
        dvoldt[n] /= (*var.volume_n)[n];
    }

    // std::cout << "dvoldt:\n";
    // print(std::cout, dvoldt);
    // std::cout << "\n";
#ifdef NPROF
    nvtxRangePop();
#endif
}


void compute_edvoldt(const Variables &var, double_vec &dvoldt,
                     double_vec &edvoldt)
{
#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    /* edvoldt is the averaged (i.e. smoothed) dvoldt on the element.
     * It is used in update_stress() to prevent mesh locking.
     */

#ifndef ACC
    #pragma omp parallel for default(none)      \
        shared(var, dvoldt, edvoldt)
#endif
    #pragma acc parallel loop gang vector async
    for (int e=0; e<var.nelem; ++e) {
        ConstConnAccessor conn = (*var.connectivity)[e];
        double dj = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            int n = conn[i];
            dj += dvoldt[n];
        }
        edvoldt[e] = dj / NODES_PER_ELEM;
    }
#ifdef NPROF
    nvtxRangePop();
#endif
    // std::cout << "edvoldt:\n";
    // print(std::cout, edvoldt);
    // std::cout << "\n";
}


void NMD_stress(const Variables &var, tensor_t& stress, double_vec &dp_nd, double_vec &etmp)
{
#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    // dp_nd is the pressure change, weighted by the element volume,
    // lumped onto the nodes.

#ifndef ACC
    #pragma omp parallel for default(none) shared(var,etmp)
#endif
    #pragma acc parallel loop gang vector async
    for (int e=0;e<var.nelem;e++) {
        etmp[e] = (*var.dpressure)[e] * (*var.volume)[e];
    }

#ifndef ACC
    #pragma omp parallel for default(none) shared(var,dp_nd,etmp)
#endif
    #pragma acc parallel loop gang vector async
    for (int n=0;n<var.nnode;n++) {
        dp_nd[n] = 0;
        for( auto e = (*var.support)[n].begin(); e < (*var.support)[n].end(); ++e)
            dp_nd[n] += etmp[*e];
        dp_nd[n] /= (*var.volume_n)[n];
    }

    // dp_el is the averaged (i.e. smoothed) dp_nd on the element.
#ifndef ACC
    #pragma omp parallel for default(none) shared(var, dp_nd, stress)
#endif
    #pragma acc parallel loop gang vector async
    for (int e=0; e<var.nelem; ++e) {
        ConstConnAccessor conn = (*var.connectivity)[e];
        double dp = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            int n = conn[i];
            dp += dp_nd[n];
        }
        double dp_el = dp / NODES_PER_ELEM;

    	TensorAccessor s = stress[e];

	    double dp_orig = (*var.dpressure)[e];
        double ddp = ( - dp_orig + dp_el ) / NDIMS;
	    for (int i=0; i<NDIMS; ++i)
            s[i] += ddp;
    }

#ifdef NPROF
    nvtxRangePop();
#endif
}

// ============================================================
// SPR (Superconvergent Patch Recovery) stress computation
// ============================================================

// Fit a linear polynomial sigma*(dx, dy) = a0 + a1*dx + a2*dy to the stress
// values at element centroids in a node's patch, using RELATIVE coordinates
// (centroid of the patch's element centroids as origin).  Evaluate at the
// node's relative coordinates.
//
// Using relative coordinates reduces the absolute magnitudes from ~50,000 m
// to ~300 m (one mesh element), dropping the 3×3 normal-matrix condition
// number from ~2.5e9 to ~2.3e4 — numerically safe for double precision.
//
// Returns true if the Cramer solve succeeded; false if degenerate (fallback).
#pragma acc routine seq
static bool spr_solve_centered(const double A[3][3], const double b[3],
                               const double dxi, const double dyi,
                               double &result)
{
    const double M00 = A[1][1]*A[2][2] - A[1][2]*A[1][2];
    const double M01 = A[0][1]*A[2][2] - A[1][2]*A[0][2];
    const double M02 = A[0][1]*A[1][2] - A[1][1]*A[0][2];
    const double det = A[0][0]*M00 - A[0][1]*M01 + A[0][2]*M02;

    // Scale-adaptive singularity threshold: A[1][1] ~ Σdx², A[2][2] ~ Σdz²
    // both O(h²) with patch-relative coords — well-conditioned.
    const double scale = A[0][0] * (A[1][1] > A[2][2] ? A[1][1] : A[2][2]);
    if (scale == 0.0 || fabs(det) < 1e-6 * scale)
        return false;

    const double inv_det = 1.0 / det;

    const double a0 = inv_det * (
          b[0]    * M00
        - A[0][1] * (b[1]*A[2][2] - A[1][2]*b[2])
        + A[0][2] * (b[1]*A[1][2] - A[1][1]*b[2]));

    const double a1 = inv_det * (
          A[0][0] * (b[1]*A[2][2] - A[1][2]*b[2])
        - b[0]    * M01
        + A[0][2] * (A[0][1]*b[2] - b[1]*A[0][2]));

    const double a2 = inv_det * (
          A[0][0] * (A[1][1]*b[2] - b[1]*A[1][2])
        - A[0][1] * (A[0][1]*b[2] - b[1]*A[0][2])
        + b[0]    * M02);

    result = a0 + a1*dxi + a2*dyi;
    return true;
}

#ifdef THREED
// Fit linear polynomial sigma*(dx,dy,dz) = a0 + a1*dx + a2*dy + a3*dz to element
// centroid values in the patch.  Uses Gaussian elimination with partial pivoting
// on the 4×4 normal system.  Returns true on success; false if the patch is
// degenerate (fewer than 4 independent directions → fallback).
#pragma acc routine seq
static bool spr_solve_centered_3d(const double A[4][4], const double b[4],
                                   const double dxi, const double dyi, const double dzi,
                                   double &result)
{
    // Per-column ∞-norms of the original A, used as per-step singularity thresholds.
    // Column 0 scale ~ n (count); columns 1–3 scale ~ n·h² (patch-relative coords).
    double col_scale[4] = {};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            const double v = fabs(A[i][j]);
            if (v > col_scale[j]) col_scale[j] = v;
        }
    if (col_scale[0] == 0.0) return false;

    // Augmented system [A | b]
    double M[4][5];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) M[i][j] = A[i][j];
        M[i][4] = b[i];
    }

    for (int col = 0; col < 4; ++col) {
        int pivot_row = col;
        double pval = fabs(M[col][col]);
        for (int row = col+1; row < 4; ++row) {
            const double v = fabs(M[row][col]);
            if (v > pval) { pval = v; pivot_row = row; }
        }
        // col_scale[col]==0 means that entire column of A was zero (degenerate patch,
        // e.g. all centroids at the same depth) → rank-deficient → fall back.
        if (col_scale[col] == 0.0 || pval < 1e-6 * col_scale[col]) return false;

        if (pivot_row != col)
            for (int j = 0; j <= 4; ++j) {
                double tmp = M[col][j]; M[col][j] = M[pivot_row][j]; M[pivot_row][j] = tmp;
            }

        const double inv_pivot = 1.0 / M[col][col];
        for (int row = col+1; row < 4; ++row) {
            const double factor = M[row][col] * inv_pivot;
            for (int j = col; j <= 4; ++j) M[row][j] -= factor * M[col][j];
        }
    }

    double x[4];
    for (int i = 3; i >= 0; --i) {
        double s = M[i][4];
        for (int j = i+1; j < 4; ++j) s -= M[i][j] * x[j];
        x[i] = s / M[i][i];
    }

    result = x[0] + x[1]*dxi + x[2]*dyi + x[3]*dzi;
    return true;
}
#endif // THREED

// Volume-weighted average of element stress values in the patch.
// Used as fallback when the SPR normal matrix is near-singular, and as a
// clamp range for the polynomial result to prevent out-of-range extrapolation.
#pragma acc routine seq
static double spr_volume_weighted_avg_strided(const int npatch, const int* patch, 
                                              const double_vec &elem_volume, 
                                              const double* ptr, 
                                              int stride) 
{
    double sum_w = 0.0, sum_ws = 0.0;
    for (int jp = 0; jp < npatch; ++jp) {
        int e = patch[jp];
        double w = elem_volume[e];
        sum_ws += w * ptr[e * stride]; 
        sum_w  += w;
    }
    return (sum_w > 0.0) ? sum_ws / sum_w : 0.0;
}

// # of max spr fields, 4 for 2D, 6 for 3D
#define MAX_SPR_FIELDS (NDIMS*2)

// Compute SPR nodal values
static void spr_fused_fields(const Variables &var, 
                             const double* const in_ptrs[], const int in_strides[],
                             double* const out_ptrs[], const int out_strides[], 
                             int num_fields)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif

#ifndef ACC
    #pragma omp parallel for default(none) shared(var, in_ptrs, in_strides, out_ptrs, \
                out_strides, num_fields)
#endif
    #pragma acc parallel loop gang vector async copyin(in_ptrs[0:num_fields], \
                in_strides[0:num_fields], out_ptrs[0:num_fields], out_strides[0:num_fields])
    for (int i = 0; i < var.nnode; ++i) {
        const int npatch = var.sup.size(i);

        double smin[MAX_SPR_FIELDS], smax[MAX_SPR_FIELDS];
        const int* patch = var.sup.patch(i);
        const int e0 = patch[0];
        
        #pragma acc loop seq
        for (int f = 0; f < num_fields; ++f) {
            smin[f] = smax[f] = in_ptrs[f][e0 * in_strides[f]];
        }

        double x0 = 0.0, y0 = 0.0;
#ifdef THREED
        double z0 = 0.0;
#endif

        // calculate patch centroid (x0, y0, z0) and field value extrema (smin, smax) for clamping
        #pragma acc loop seq
        for (int jp = 0; jp < npatch; ++jp) {
            const int e = patch[jp];
            ConstArrayIndirectAccessor ck = var.coord->view_const((*var.connectivity)[e]);
            double cx = 0.0, cy = 0.0;
#ifdef THREED
            double cz = 0.0;
#endif
            #pragma acc loop seq
            for (int k = 0; k < NODES_PER_ELEM; ++k) {
                cx += ck[k][0];
                cy += ck[k][1];
#ifdef THREED
                cz += ck[k][2];
#endif
            }
            cx /= NODES_PER_ELEM;
            cy /= NODES_PER_ELEM;
            x0 += cx;
            y0 += cy;
#ifdef THREED
            cz /= NODES_PER_ELEM;
            z0 += cz;
#endif

            #pragma acc loop seq
            for (int f = 0; f < num_fields; ++f) {
                double s = in_ptrs[f][e * in_strides[f]]; // use Stride
                if (s < smin[f]) smin[f] = s;
                if (s > smax[f]) smax[f] = s;
            }
        }
        x0 /= npatch; y0 /= npatch;
#ifdef THREED
        z0 /= npatch;
#endif
        double A[NODES_PER_ELEM][NODES_PER_ELEM] = {};
        double b_vec[MAX_SPR_FIELDS][NODES_PER_ELEM] = {};

        // build A and b_vec for all fields in the patch, 
        // with shared loops and Stride-aware access
        #pragma acc loop seq
        for (int jp = 0; jp < npatch; ++jp) {
            const int e = patch[jp];
            ConstArrayIndirectAccessor ck = var.coord->view_const((*var.connectivity)[e]);
            double cx = 0.0, cy = 0.0;
#ifdef THREED
            double cz = 0.0;
#endif
            #pragma acc loop seq
            for (int k = 0; k < NODES_PER_ELEM; ++k) {
                cx += ck[k][0];
                cy += ck[k][1];
#ifdef THREED
                cz += ck[k][2];
#endif
            }
            cx /= NODES_PER_ELEM;
            cy /= NODES_PER_ELEM;
#ifdef THREED
            cz /= NODES_PER_ELEM;
#endif

            const double dx = cx - x0;
            const double dy = cy - y0;
#ifdef THREED
            const double dz = cz - z0;
#endif

            A[0][0] += 1.0;
            A[0][1] += dx;    A[1][0] = A[0][1];
            A[0][2] += dy;    A[2][0] = A[0][2];
            A[1][1] += dx*dx;
            A[1][2] += dx*dy; A[2][1] = A[1][2];
            A[2][2] += dy*dy;
#ifdef THREED
            A[0][3] += dz;    A[3][0] = A[0][3];
            A[1][3] += dx*dz; A[3][1] = A[1][3];
            A[2][3] += dy*dz; A[3][2] = A[2][3];
            A[3][3] += dz*dz;
#endif

            #pragma acc loop seq
            for (int f = 0; f < num_fields; ++f) {
                double s = in_ptrs[f][e * in_strides[f]]; // use Stride
                b_vec[f][0] += s;
                b_vec[f][1] += dx * s;
                b_vec[f][2] += dy * s;
#ifdef THREED
                b_vec[f][3] += dz * s;
#endif
            }
        }

        const ConstArrayAccessor ci = (*var.coord)[i];
        const double dxi = ci[0] - x0;
        const double dyi = ci[1] - y0;
#ifdef THREED
        const double dzi = ci[2] - z0;
#endif

        // solve A for all fields in the patch, with fallback and clamping
        #pragma acc loop seq
        for (int f = 0; f < num_fields; ++f) {
            double temp_val;

#ifdef THREED
            bool solve_success = spr_solve_centered_3d(A, b_vec[f], dxi, dyi, dzi, temp_val);
#else
            bool solve_success = spr_solve_centered(A, b_vec[f], dxi, dyi, temp_val);
#endif

            temp_val = solve_success ?
                fmax(smin[f], fmin(smax[f], temp_val)) :
                // if the matrix is degenerate, the polynomial fit is unreliable. 
                // Fall back to volume-weighted average, which is stable but less accurate. 
                // The clamping range is still valid since it's based on the patch values.
                spr_volume_weighted_avg_strided(var.sup.size(i), var.sup.patch(i),
                                                *var.volume, in_ptrs[f], in_strides[f]);

            out_ptrs[f][i * out_strides[f]] = temp_val;
        }
    }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

void spr_elem_to_node(const Param &param, const Variables &var,
                      tensor_t *stress_n, double_vec *stressyy_n)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    // ----------------------------------------------------------------
    // Step A: SPR — recover smooth nodal stresses on the OLD mesh.
    // Must happen before prepare_interpolation (reads old var.stress,
    // *var.support, old_coord, old_connectivity — all still old here).
    //
    // Pressure-centering: subtract ref_pressure from each old element's
    // diagonal stress components before SPR, so the polynomial fits the
    // small deviatoric residual rather than the large lithostatic pressure.
    // Without centering, independent polynomial fits for σxx and σzz on a
    // one-sided surface-node patch differ by the amount of the actual
    // deviatoric stress, but the large-pressure background amplifies any
    // extrapolation error at boundary nodes. Cold surface ice (η≈visc_max)
    // has a Maxwell relaxation time of ~250 Myr, so even a tiny artifact
    // persists throughout the simulation.
    // ----------------------------------------------------------------

#ifndef ACC
    #pragma omp parallel for default(none) shared(param, var)
#endif
    #pragma acc parallel loop gang vector async    
    for (int e = 0; e < var.nelem; ++e) {
        ConstConnAccessor conn = (*var.connectivity)[e];
        double z = 0;
        for (int k = 0; k < NODES_PER_ELEM; ++k)
            z += (*var.coord)[conn[k]][NDIMS-1];
        z /= NODES_PER_ELEM;
        double p_ref_old = ref_pressure(param, z);
        TensorAccessor s = (*var.stress)[e];
        // Shift diagonal components: add p_ref so stress becomes small deviatoric
        for (int d = 0; d < NDIMS; ++d) s[d] += p_ref_old;

        if (param.mat.is_plane_strain)
            (*var.stressyy)[e] += p_ref_old;
    }

    // dynamic registration of all fields involved in SPR
    const double* in_ptrs[MAX_SPR_FIELDS];
    int in_strides[MAX_SPR_FIELDS];
    double* out_ptrs[MAX_SPR_FIELDS];
    int out_strides[MAX_SPR_FIELDS];
    int num_fields = 0;

    auto add_field_acc = [&](tensor_t::ConstComponentAccessor in_acc, tensor_t::ComponentAccessor out_acc) {
        in_ptrs[num_fields]     = in_acc.ptr_;
        in_strides[num_fields]  = in_acc.stride_;
        out_ptrs[num_fields]    = out_acc.ptr_;
        out_strides[num_fields] = out_acc.stride_;
        num_fields++;
    };

    auto add_field_ptr = [&](const double* in_p, double* out_p) {
        in_ptrs[num_fields]     = in_p;
        in_strides[num_fields]  = 1;
        out_ptrs[num_fields]    = out_p;
        out_strides[num_fields] = 1;
        num_fields++;
    };

    for (int d = 0; d < NSTR; ++d)
        add_field_acc(var.stress->component_const(d), stress_n->component(d));

    if (param.mat.is_plane_strain)
        add_field_ptr(var.stressyy->data(), stressyy_n->data());

    // Compute SPR-recovered nodal stresses from element-centroid stresses on the old mesh.
    // Uses Zienkiewicz-Zhu superconvergent patch recovery with linear polynomial basis [1, x, y].
    // call the fused SPR kernel for all registered fields
    spr_fused_fields(var, in_ptrs, in_strides, out_ptrs, out_strides, num_fields);

    #pragma acc wait

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

void spr_node_to_elem(const Param &param, const Variables &var, 
                      tensor_t *stress, double_vec *stressyy)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    // ----------------------------------------------------------------
    // Step C: Average SPR nodal stresses -> new element stresses.
    // ----------------------------------------------------------------
    for (int d = 0; d < NSTR; ++d)
        average_nodal_to_elem(var.stress_n->component_const(d), *var.connectivity,
                              var.nelem, stress->component(d));

    if (param.mat.is_plane_strain)
        average_nodal_to_elem(static_cast<const double*>(var.stressyy_n->data()), *var.connectivity,
                                var.nelem, stressyy->data());

    // Step C': Restore reference pressure at new element centroids.
    // The SPR operated on pressure-centered stress; add back p_ref at each
    // new element's depth to recover the correct total stress.
#ifndef ACC
    #pragma omp parallel for default(none) shared(param, var, stress, stressyy)
#endif
    #pragma acc parallel loop gang vector async
    for (int e = 0; e < var.nelem; ++e) {
        ConstConnAccessor conn = (*var.connectivity)[e];
        double z = 0;
        for (int k = 0; k < NODES_PER_ELEM; ++k)
            z += (*var.coord)[conn[k]][NDIMS-1];
        z /= NODES_PER_ELEM;
        double p_ref = ref_pressure(param, z);
        TensorAccessor s = (*stress)[e];
        for (int d = 0; d < NDIMS; ++d) s[d] -= p_ref;

        if (param.mat.is_plane_strain)
            (*stressyy)[e] -= p_ref;
    }

    #pragma acc wait

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

double compute_dt(const Param& param, Variables& var)

{
#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    // constant dt
    if (param.control.fixed_dt != 0) return param.control.fixed_dt;

    // dynamic dt
    double dt_maxwell = std::numeric_limits<double>::max();
    double dt_diffusion = std::numeric_limits<double>::max();
    double dt_hydro_diffusion = std::numeric_limits<double>::max();
    double minl = std::numeric_limits<double>::max();

    // Define element velocity arrays
    // double_vec velocity_x_element(var.nelem, 0.0);
    // double_vec velocity_y_element(var.nelem, 0.0);
    // double_vec velocity_z_element(var.nelem, 0.0); // Used only for 3D
    // double vx_element = 0.0, vy_element = 0.0, vz_element = 0.0;
    // Global max velocity for elements
    // double global_max_vem = std::max(std::max(var.max_vbc_val, 0.0), 1E-20);
    double global_max_vem = 0.0;  // 
    double global_dt_min = std::numeric_limits<double>::max(); // based on length and S wave velocity

#ifndef ACC
    #pragma omp parallel for reduction(min:minl, dt_maxwell, dt_diffusion, dt_hydro_diffusion, global_dt_min) \
        reduction(max: global_max_vem) \
        default(none) shared(param, var) //, velocity_x_element, velocity_y_element, velocity_z_element) \
        private(vx_element, vy_element, vz_element)
#endif
    #pragma acc parallel loop gang vector reduction(min:minl, dt_maxwell, dt_diffusion, dt_hydro_diffusion, global_dt_min) \
        reduction(max: global_max_vem) async
    for (int e=0; e<var.nelem; ++e) {

        double vx_element = 0.0, vy_element = 0.0, vz_element = 0.0;
        int n0 = (*var.connectivity)[e][0];
        int n1 = (*var.connectivity)[e][1];
        int n2 = (*var.connectivity)[e][2];

        // calculate maxium velocity in the element
        ConstArrayIndirectAccessor v = var.vel->view_const((*var.connectivity)[e]);
        double weight = 1.0 / NODES_PER_ELEM; 

        for (int j = 0; j < NODES_PER_ELEM; ++j) {
            vx_element += v[j][0] * weight;
            vy_element += v[j][1] * weight;
#ifdef THREED
            vz_element += v[j][2] * weight;
#endif
        }

//         velocity_x_element[e] = vx_element;
//         velocity_y_element[e] = vy_element;
// #ifdef THREED
//         velocity_z_element[e] = vz_element;
// #endif

        // min height of this element
#ifdef THREED
        double max_vem = std::sqrt(vx_element*vx_element + vy_element*vy_element + vz_element*vz_element);
#else
        double max_vem = std::sqrt(vx_element*vx_element + vy_element*vy_element);
#endif

        // Find global max velocity
        global_max_vem = std::max(global_max_vem, max_vem);

        // std::cout<< "Element: " << e << " max_vem: " << global_max_vem << std::scientific << std::setprecision(5) << std::endl;

        ConstArrayAccessor a = (*var.coord)[n0];
        ConstArrayAccessor b = (*var.coord)[n1];
        ConstArrayAccessor c = (*var.coord)[n2];

        // min height of this element
        double minh;
#ifdef THREED
        {
            int n3 = (*var.connectivity)[e][3];
            ConstArrayAccessor d = (*var.coord)[n3];

            // max facet area of this tet
            double maxa = std::max(std::max(triangle_area(a, b, c),
                                            triangle_area(a, b, d)),
                                   std::max(triangle_area(c, d, a),
                                            triangle_area(c, d, b)));
            minh = 3 * (*var.volume)[e] / maxa;
        }
#else
        {
            // max edge length of this triangle
            double maxl = std::sqrt(std::max(std::max(dist2(a, b),
                                                      dist2(b, c)),
                                             dist2(a, c)));
            minh = 2 * (*var.volume)[e] / maxl;
        }
#endif
        dt_maxwell = std::min(dt_maxwell,
                              0.5 * var.mat->visc_min / (1e-40 + var.mat->shearm(e)));
        if (param.control.has_thermal_diffusion)
            dt_diffusion = std::min(dt_diffusion,
                                    0.5 * minh * minh / var.mat->therm_diff_max);
        
        // Compute dt_hydro_diffusion (hydraulic)
        if (param.control.has_hydraulic_diffusion)
            if (var.mat->hydro_diff_max > 0) {
                dt_hydro_diffusion = std::min(dt_hydro_diffusion,
                                            0.5 * minh * minh / var.mat->hydro_diff_max);
            }
        minl = std::min(minl, minh);

        // Find global min delta t to meet CFL condition
        global_dt_min = std::min(global_dt_min, minl/std::sqrt(var.mat->shearm(e)/var.mat->rho(e)) /5.0);
    }

    #pragma acc wait

    double max_vbc_val;
    if (param.control.characteristic_speed == 0) {
        max_vbc_val = var.max_vbc_val;
        if (param.control.surface_process_option > 0)
            max_vbc_val = std::max(max_vbc_val, var.surfinfo.max_surf_vel*5e-1);
    }
    else
        max_vbc_val = param.control.characteristic_speed;

    double dt_advection = 0;
    double dt_elastic = 0;
    global_max_vem = std::max(global_max_vem, var.max_vbc_val);
    var.max_global_vel_mag = global_max_vem;
    var.global_dt_min = global_dt_min;

    // Calculate dt_advection and dt_elastic 
    if(param.control.use_global_velocity_scaling)
    {
        dt_advection = 0.5 * minl / global_max_vem;
        dt_elastic = (param.control.is_quasi_static) ?
        0.5 * minl / (global_max_vem * param.control.inertial_scaling) :
        0.5 * minl / std::sqrt(param.mat.bulk_modulus[0] / param.mat.rho0[0]);
        dt_elastic = std::max(dt_elastic, global_dt_min); // Ensure dt_elastic is not smaller than global_dt_min
    }
    else
    {
        dt_advection = 0.5 * minl / max_vbc_val;
        dt_elastic = (param.control.is_quasi_static) ?
        0.5 * minl / (max_vbc_val * param.control.inertial_scaling) :
        0.5 * minl / std::sqrt(param.mat.bulk_modulus[param.mat.mattype_ref] / param.mat.rho0[param.mat.mattype_ref]);
    }

    // Combine dt calculations and incorporate dt_hydro_diffusion
    double dt = std::min({dt_elastic, dt_maxwell, dt_advection, dt_diffusion, dt_hydro_diffusion}) * param.control.dt_fraction;
    // double dt = std::min({dt_elastic, dt_maxwell, dt_advection, dt_diffusion}) * param.control.dt_fraction;
    if (param.debug.dt) {
        std::cout << "step #" << var.steps << "  dt: " << dt_maxwell << " " << dt_diffusion << " " 
                  << dt_hydro_diffusion << " " << dt_advection << " " << dt_elastic << " sec\n";
    }
    if (dt <= 0) {
        std::cerr << "Error: dt <= 0!  " << dt_maxwell << " " << dt_diffusion
                  << " " << dt_hydro_diffusion << " " << dt_advection << " " << dt_elastic << "\n";
        var.output->write_exact_error(var);
        std::exit(11);
    }
    
#ifdef NPROF
    nvtxRangePop();
#endif
    return dt;
}

double compute_dt_PT(const Param& param, const Variables& var)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    // constant dt
    if (param.control.fixed_dt != 0) return param.control.fixed_dt;

    // dynamic dt
    double dt_maxwell = std::numeric_limits<double>::max();
    double dt_diffusion = std::numeric_limits<double>::max();
    double dt_hydro_diffusion = std::numeric_limits<double>::max();
    double minl = std::numeric_limits<double>::max();

    #pragma omp parallel for reduction(min:minl,dt_maxwell,dt_diffusion,dt_hydro_diffusion)    \
        default(none) shared(param,var)
    // #pragma acc parallel loop reduction(min:minl, dt_maxwell, dt_diffusion,dt_hydro_diffusion)
    for (int e=0; e<var.nelem; ++e) {
        int n0 = (*var.connectivity)[e][0];
        int n1 = (*var.connectivity)[e][1];
        int n2 = (*var.connectivity)[e][2];

        ConstArrayAccessor a = (*var.coord)[n0];
        ConstArrayAccessor b = (*var.coord)[n1];
        ConstArrayAccessor c = (*var.coord)[n2];

        // min height of this element
        double minh;
#ifdef THREED
        {
            int n3 = (*var.connectivity)[e][3];
            ConstArrayAccessor d = (*var.coord)[n3];

            // max facet area of this tet
            double maxa = std::max(std::max(triangle_area(a, b, c),
                                            triangle_area(a, b, d)),
                                   std::max(triangle_area(c, d, a),
                                            triangle_area(c, d, b)));
            minh = 3 * (*var.volume)[e] / maxa;
        }
#else
        {
            // max edge length of this triangle
            double maxl = std::sqrt(std::max(std::max(dist2(a, b),
                                                      dist2(b, c)),
                                             dist2(a, c)));
            minh = 2 * (*var.volume)[e] / maxl;
        }
#endif
        dt_maxwell = std::min(dt_maxwell,
                              0.5 * var.mat->visc_min / (1e-40 + var.mat->shearm(e)));
        // if (param.control.has_thermal_diffusion)
        //     dt_diffusion = std::min(dt_diffusion,
        //                             0.5 * minh * minh / var.mat->therm_diff_max);
        
        // // Compute dt_hydro_diffusion (hydraulic)
        // if (var.mat->hydro_diff_max > 0) {
        //     dt_hydro_diffusion = std::min(dt_hydro_diffusion,
        //                                   0.5 * minh * minh / var.mat->hydro_diff_max);
        // }
        minl = std::min(minl, minh);
    }


    // max_vbc_val is maximum boundary velocity
    double max_vbc_val;
    if (param.control.characteristic_speed == 0) {
        max_vbc_val = var.max_vbc_val; 

        if (param.control.surface_process_option > 0)
            max_vbc_val = std::max(max_vbc_val, var.surfinfo.max_surf_vel*5e-1);
    }
    else
        max_vbc_val = param.control.characteristic_speed;

    double dt_advection = 0.5 * minl / max_vbc_val;
    double dt_elastic = (param.control.is_quasi_static) ?
        0.5 * minl / (max_vbc_val * param.control.inertial_scaling) :
        0.5 * minl / std::sqrt(param.mat.bulk_modulus[param.mat.mattype_ref] / param.mat.rho0[param.mat.mattype_ref]);

    double dt = std::min({dt_elastic, dt_maxwell, dt_advection}) * param.control.dt_fraction;
    if (param.debug.dt) {
        std::cout << "step #" << var.steps << "  dt: " << dt_maxwell << " " << dt_advection << " " << dt_elastic << " sec\n";
    }
    if (dt <= 0) {
        std::cerr << "Error: dt <= 0!  " << dt_maxwell << " "  << dt_advection << " " << dt_elastic << "\n";
        var.output->write_exact_error(var);
        std::exit(11);
    }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
    return dt;
}

void compute_mass(const Param &param, const Variables &var,
                  double max_vbc_val, double_vec &volume_n,
                  double_vec &mass, double_vec &tmass, double_vec &hmass, double_vec &ymass, elem_cache &tmp_result)

{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    // volume_n is (node-averaged volume * NODES_PER_ELEM)
    // volume_n.assign(volume_n.size(), 0);
    // mass.assign(mass.size(), 0);
    // tmass.assign(tmass.size(), 0);

    const double pseudo_speed = max_vbc_val * param.control.inertial_scaling; // for non-ATP using max velocity on boundary
    const double pseudo_speed_ATP = var.max_global_vel_mag * param.control.inertial_scaling; // for ATP using global max velocity

    double diff_e;

    #pragma acc serial async
    {
        // Retrieve hydraulic properties for the element
        double perm_e = var.mat->perm(param.mat.mattype_ref);                // Intrinsic permeability 
        double mu_e = var.mat->mu_fluid(param.mat.mattype_ref);              // Fluid dynamic viscosity
        double alpha_b = var.mat->alpha_biot(param.mat.mattype_ref);         // Biot coefficient
        double rho_f = var.mat->rho_fluid(param.mat.mattype_ref);            // Fluid density
        double phi_e = var.mat->phi(param.mat.mattype_ref);        // Element porosity
        double comp_fluid = var.mat->beta_fluid(param.mat.mattype_ref);        // fluid comporessibility
        double bulkm = var.mat->bulkm(param.mat.mattype_ref);
        double shearm = var.mat->shearm(param.mat.mattype_ref);
        double matrix_comp = 1.0 / (bulkm +4.0*shearm/3.0);

        rho_f = 1000.0; 
        double gamma_w = rho_f * param.control.gravity; // specific weight
        
        // Hydraulic conductivity using permeability and viscosity
        double hydraulic_conductivity = perm_e * gamma_w / mu_e;
        
        // Compute element diffusivity and update max using reduction
        diff_e = hydraulic_conductivity / (phi_e * comp_fluid + alpha_b * matrix_comp) / gamma_w;
    }

    #pragma acc wait

    if (pseudo_speed < diff_e && param.control.has_hydraulic_diffusion)
    {
        std::cout << "pseudo speed is too slow, increase mass scaling" << std::endl;
        std::exit(11);
    }

#ifndef ACC
#ifdef GPP1X
    #pragma omp parallel default(none) shared(var, param, volume_n, mass, tmass, hmass, ymass, pseudo_speed, pseudo_speed_ATP, tmp_result)
#else
    #pragma omp parallel default(none) shared(var, param, volume_n, mass, tmass, hmass, ymass, tmp_result)
#endif
#endif
    {
#ifndef ACC
        #pragma omp for
#endif
        #pragma acc parallel loop gang vector async
        for (int e=0;e<var.nelem;e++) {
            ElemCacheAccessor tr = tmp_result[e];
            double rho;

            if(param.control.use_global_velocity_scaling)
            {
                double apprent_speed = std::min(pseudo_speed_ATP, std::sqrt(var.mat->shearm(e)/var.mat->rho(e))); // minimum speed

                rho = (param.control.is_quasi_static) ?
                (*var.mat).bulkm(e) / (apprent_speed * apprent_speed) :  // pseudo density for quasi-static sim
                (*var.mat).rho(e);  // true density for dynamic sim

                if (param.control.has_hydraulic_diffusion && (param.control.is_quasi_static == false)) {
                    // Modified density considering porosity for hydraulic diffusion
                        rho = (*var.mat).rho(e) * (1 - (*var.mat).phi(e)) + 1000.0 * (*var.mat).phi(e);
                    }
            }
            else
            {    
                rho = (param.control.is_quasi_static) ?
                (*var.mat).bulkm(e) / (pseudo_speed * pseudo_speed) :  // pseudo density for quasi-static sim
                (*var.mat).rho(e);  // true density for dynamic sim

                if (param.control.has_hydraulic_diffusion && (param.control.is_quasi_static == false)) {
                    // Modified density considering porosity for hydraulic diffusion
                        rho = (*var.mat).rho(e) * (1 - (*var.mat).phi(e)) + 1000.0 * (*var.mat).phi(e);
                    }

            }

            double bulk_comp = 1.0/(*var.mat).bulkm(e); // lambda + 2G/3
            if(NDIMS == 2) bulk_comp = 1.0/((*var.mat).bulkm(e) + (*var.mat).shearm(e)/3.0); // lambda + G 

            double hm_coeff = (*var.mat).alpha_biot(e) + (*var.mat).phi(e) - (*var.mat).alpha_biot(e) * (*var.mat).phi(e); 
            double m = rho * (*var.volume)[e] / NODES_PER_ELEM;
            double tm = (*var.mat).rho(e) * (*var.mat).cp(e) * (*var.volume)[e] / NODES_PER_ELEM;
            double hm = (hm_coeff * bulk_comp + (*var.mat).phi(e) * (*var.mat).beta_fluid(e)) * (*var.volume)[e] / NODES_PER_ELEM;
            double ym = 9 * (*var.mat).bulkm(e) * (*var.mat).shearm(e) / (3 * (*var.mat).bulkm(e) + (*var.mat).shearm(e)) / NODES_PER_ELEM; // Young's modulus

            tr[0] = (*var.volume)[e];
            tr[1] = m;
            if (param.control.has_thermal_diffusion)
                tr[2] = tm;
            tr[3] = hm; // check
            tr[4] = ym; 
        }

#ifndef ACC
        #pragma omp for
#endif
        #pragma acc parallel loop gang vector async
        for (int n=0;n<var.nnode;n++) {
            volume_n[n]=0;
            mass[n]=0;
            tmass[n]=0;
            hmass[n]=0;
            ymass[n]=0;
        
            for( auto e = (*var.support)[n].begin(); e < (*var.support)[n].end(); ++e) {
                ConstElemCacheAccessor tr = tmp_result[*e];
                volume_n[n] += tr[0];
                mass[n] += tr[1];
                if (param.control.has_thermal_diffusion)
                    tmass[n] += tr[2];
                hmass[n] += tr[3];
                ymass[n] += tr[4];
            }
        }
    }

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}


double elem_quality(const array_t &coord, const conn_t &connectivity,
                    const double_vec &volume, int e)
{
    /* This function returns the quality (0~1) of the element.
     * The quality of an equidistant (i.e. best quality) tetrahedron/triangle is 1.
     */
    double quality;
    double vol = volume[e];
    int n0 = connectivity[e][0];
    int n1 = connectivity[e][1];
    int n2 = connectivity[e][2];

    ConstArrayAccessor a = coord[n0];
    ConstArrayAccessor b = coord[n1];
    ConstArrayAccessor c = coord[n2];

#ifdef THREED
    {
        int n3 = connectivity[e][3];
        ConstArrayAccessor d = coord[n3];
        double normalization_factor = 216 * std::sqrt(3);

        double area_sum = (triangle_area(a, b, c) +
                           triangle_area(a, b, d) +
                           triangle_area(c, d, a) +
                           triangle_area(c, d, b));
        quality = normalization_factor * vol * vol / (area_sum * area_sum * area_sum);
    }
#else
    {
        double normalization_factor = 4 * std::sqrt(3);

        double dist2_sum = dist2(a, b) + dist2(b, c) + dist2(a, c);
        quality = normalization_factor * vol / dist2_sum;
    }
#endif

    return quality;
}


double worst_elem_quality(const array_t &coord, const conn_t &connectivity,
                          const double_vec &volume, int &worst_elem)
{
    double q = 1;
    worst_elem = 0;
    for (std::size_t e=0; e<volume.size(); e++) {
        double quality = elem_quality(coord, connectivity, volume, e);
        if (quality < q) {
            q = quality;
            worst_elem = e;
        }
    }
    return q;
}

