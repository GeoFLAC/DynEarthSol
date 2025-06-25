#include "iostream"
#include <numeric>
#include <stdexcept>
#ifdef USE_NPROF
#include <nvToolsExt.h> 
#endif

#include "constants.hpp"
#include "parameters.hpp"
#include "barycentric-fn.hpp"
#include "mesh.hpp"
#include "utils.hpp"
#include "knn.hpp"
#include "nn-interpolation.hpp"


namespace {

    void find_nearest_neighbor(const Variables &var, KNN &kdtree,
                               int_vec &idx, int_vec &is_changed,
                               int_vec &idx_changed, int_vec &changed)
    {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif

        double eps = 1e-15;

        array_t queries(var.nelem);

        elem_center(*var.coord, *var.connectivity, queries);

        neighbor_vec neighbors(var.nelem);

        kdtree.search(queries, neighbors, 1, 3.);

        #pragma acc parallel loop
        for (int e=0; e<var.nelem; e++) {
            idx[e] = int(neighbors[e].idx);
            is_changed[e] = (neighbors[e].dist2 < eps) ? 0 : 1;
        }

        int nchanged = 0;
        for (int e=0; e<var.nelem; e++) {
            if (is_changed[e] > 0) {
                changed.push_back(e);
                idx_changed[e] = nchanged;
                nchanged++;
            }
        }

#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }

    void find_acm_elem_ratios(const Variables &var,
                              const Barycentric_transformation &bary,
                              int_vec &is_changed,
                              KNN &kdtree,
                              int old_nelem,
                              int_vec &elems_vec,
                              double_vec &ratios_vec,
                              int_vec &changed)
    {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif
        const int neta0 = 10; // larger neta0, more accurate mapping
        const int neta1 = neta0 + 1; // different from neta0 to prevent the temporary point falling the edge of elements
        const int neta2 = neta0;
        const double spacing0 = 1.0 / neta0;
        const double spacing1 = 1.0 / neta1;
        const double spacing2 = 1.0 / neta2;
        const int max_el = std::min(32, old_nelem);
        const double eps = 0;

        int nchanged = changed.size();

        elems_vec.resize(nchanged*32,-1);
        ratios_vec.resize(nchanged*32);

        double_vec2D sample_eta;
        for (int i=0; i<neta0; i++)
            for (int j=0; j<neta1; j++) {
#ifdef THREED
                for (int k=0; k<neta2; k++) {
                    double eta[4] = {(i + 0.5) * spacing0,
                                    (j + 0.5) * spacing1,
                                    (k + 0.5) * spacing2,
                                    1 - (i + 0.5) * spacing0 - (j + 0.5) * spacing1 - (k + 0.5) * spacing2};
#else
                    double eta[3] = {(i + 0.5) * spacing0,
                                    (j + 0.5) * spacing1,
                                    1 - (i + 0.5) * spacing0 - (j + 0.5) * spacing1};
#endif
                    if (eta[NODES_PER_ELEM-1] < 0) continue;

                    sample_eta.push_back(double_vec(eta, eta + NODES_PER_ELEM));
#ifdef THREED
                }
#endif
            }
        int nsample = sample_eta.size();
        int nqueries = nchanged * nsample;

        // number of neighbors exceeding computational limit
        int queries_max = 1024 * 1024 * 16;

        neighbor_vec neighbors;

        int elems_per_block = queries_max / nsample;
        if (elems_per_block < 1) elems_per_block = 1;
        int nblocks = (nchanged + elems_per_block - 1) / elems_per_block;
        printf("  Using %d blocks, elements per block: %d, total queries: %d\n",
               nblocks, elems_per_block, nqueries);

        array_t queries(elems_per_block*nsample);

        for (int b=0; b<nblocks; b++) {
            int start = b * elems_per_block;
            int end = std::min((b + 1) * elems_per_block, nchanged);
            if (start >= end) continue;

            printf("    Block %3d: element from %7d to %7d", b, start, end);

            queries.resize((end-start) * nsample);

#ifndef ACC
            #pragma omp parallel for default(none) shared(var,start, end, \
                sample_eta, nsample, queries, changed)
#endif
            #pragma acc parallel loop async
            for (int i=start; i<end; i++) {
                int e = changed[i];
                int query_start = i - start;
                const int* conn = (*var.connectivity)[e];

                for (int j=0; j<nsample; j++) {
                    double *x = queries[query_start*nsample + j];
                    for (int d=0; d<NDIMS; d++) {
                        x[d] = 0;
                        for (int n=0; n<NODES_PER_ELEM; n++)
                            x[d] += (*var.coord)[ conn[n] ][d] * sample_eta[j][n];
                    }
                }
            }

            neighbors.resize((end-start) * nsample * max_el);

            #pragma acc wait

            kdtree.search(queries, neighbors, max_el, 3.);

#ifndef ACC
            #pragma omp parallel for default(none) shared(var, bary, is_changed, \
                elems_vec, ratios_vec, sample_eta, \
                nsample, nchanged, changed, queries, neighbors, start, end) firstprivate(max_el)
#endif
            #pragma acc parallel loop async
            for (int i=start; i<end; i++) {
                int e = changed[i];
                int elem_count_buf[32] = {0};
                int elem_keys[32] = {0};
                int elem_size = 0;
                /* Procedure:
                *   1. Create a bunch of temporary points, uniformly distributed in the element.
                *   2. Locate these points on old elements.
                *   3. The percentage of points in each old elements is used as (approximate)
                *      volume weighting for the mapping.
                */
                const int* conn = (*var.connectivity)[e];
                for (int j=0; j<nsample; j++) {

                    double x[NDIMS] = {0}; // coordinate of temporary point

                    // find the nearest point nn in old_center
                    int query_start = i - start;
                    for (int d=0; d<NDIMS; d++)
                        x[d] = queries[query_start * nsample + j][d];

                    neighbor *nn_idx = neighbors.data() + (query_start * nsample + j) * max_el;

                    // bool is_consist = true;
                    // for (int jj=0; jj<max_el; jj++) {
                    //     // compare the nn_idx[jj] with neighbors[jj].idx
                    //     if (nn_idx[jj] != nn_idx_ptr[jj].idx) {
                    //         is_consist = false;
                    //         break;
                    //     }
                    // }
                    // if (!is_consist) {
                    //     printf("Inconsistent neighbors for element %d, query %d:\n", e, count);
                    //     for (int jj=0; jj<max_el; jj++) {
                    //         printf("  nn_idx[%4d] = %6d, dist2 = %10.1f | ", jj, nn_idx[jj], out_dists_sqr[jj]);
                    //         printf("  neighbors[%4d] = %6d, dist2 = %10.1f", jj, nn_idx_ptr[jj].idx, nn_idx_ptr[jj].dist2);
                    //         if (nn_idx[jj] != nn_idx_ptr[jj].idx) {
                    //             printf(" *\n");
                    //         } else {
                    //             printf("  \n");
                    //         }
                    //     }
                    //     printf("\n");
                    // }

                    // std::cout << "  ";
                    // print(std::cout, eta, NODES_PER_ELEM);
                    // print(std::cout, x, NDIMS);

                    // find the old element that is enclosing x
                    double r[NDIMS];
                    int old_e;
                    for (int jj=0; jj<max_el; jj++) {
                        old_e = nn_idx[jj].idx;
                        bary.transform(x, old_e, r);
                        if (bary.is_inside(r)) {
                            bool found = false;
                            for (int ei = 0; ei < elem_size; ++ei) {
                                if (elem_keys[ei] == old_e) {
                                    elem_count_buf[ei]++;
                                    found = true;
                                    break;
                                }
                            }
                            if (!found && elem_size < 32) {
                                elem_keys[elem_size] = old_e;
                                elem_count_buf[elem_size] = 1;
                                elem_size++;
                            }
                            break;
                        }
                    }

                    /* not found, do nothing */
                    // std::cout << " not found\n";
                    continue;
                }

                // Count
                int total_count = 0;
                for (int k=0; k<=elem_size; ++k)
                    total_count += elem_count_buf[k];

                // std::cout << "  has " << total_count << " points\n";

                if (total_count == 0) {
                    // This happens when new material is added during remeshing,
                    // and this element is completely within the new material.
                    // Mark the element as unchanged instead to keep the result
                    // of nearest neighbor interpolation.
                    is_changed[e] = -1;
                    continue;
                }

                if (elem_size == 1) {
                    // This happens when the new is completely within the old element.
                    // Mark the element as unchanged instead to keep the result
                    // of nearest neighbor interpolation.
                    is_changed[e] = -1;
                    continue;
                }

                const double inv = 1.0 / total_count;
                for (int k=0; k<elem_size; ++k) {
                    elems_vec[i*32+k] = elem_keys[k];
                    ratios_vec[i*32+k] = elem_count_buf[k] * inv;
                }
            }
        } // end of for (int b=0; b<nblocks; b++)

#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }


    void prepare_interpolation(const Param& param, Variables &var,
                               const Barycentric_transformation &bary,
                               const array_t &old_coord,
                               const conn_t &old_connectivity,
                               int_vec &idx,
                               int_vec &is_changed,
                               int_vec &idx_changed,
                               int_vec &elems_vec,
                               double_vec &ratios_vec)
    {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif

        int_vec changed;
        int old_nelem = old_connectivity.size();

#ifdef USE_NPROF
        nvtxRangePushA("create kdtree for old elements");
#endif

        array_t points(old_nelem);
        elem_center(old_coord, old_connectivity, points);

#ifdef ACC
        array_t point_tmp(1);
        PointCloud cloud(point_tmp);
#else
        PointCloud cloud(points);
#endif

        NANOKDTree nano_kdtree(NDIMS, cloud);
        KNN kdtree(param, points, nano_kdtree);

#ifdef USE_NPROF
        nvtxRangePop();
#endif
        printf("    Finding nearest neighbor...\n");
        find_nearest_neighbor(var, kdtree, idx, is_changed, idx_changed, changed);

        printf("    Finding acm element ratios...\n");
        find_acm_elem_ratios(var, bary, is_changed, kdtree, old_nelem, elems_vec, ratios_vec, changed);

#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }


    void inject_field(const int_vec &idx,
                      const int_vec &is_changed,
                      const int_vec &idx_changed,
                      const int_vec &elems_vec,
                      const double_vec &ratios_vec,
                      const double_vec &source,
                      double_vec &target)
    {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif

        #pragma acc serial async
        int ntarget = target.size();

#ifndef ACC
        #pragma omp parallel default(none) shared(idx, source, \
            target, is_changed, idx_changed, ratios_vec, elems_vec, ntarget)
#endif
        {
#ifndef ACC
            #pragma omp for
#endif
            #pragma acc parallel loop async
            for (int i=0; i<ntarget; i++) {
                int n = idx[i];
                target[i] = source[n];
            }

#ifndef ACC
            #pragma omp for
#endif
            #pragma acc parallel loop async
            for (int i=0; i<ntarget; i++) {
                if (is_changed[i]>0) {
                    int n = idx_changed[i];

                    target[i] = 0;
                    for (int j=0; j<32; j++) {
                        if (elems_vec[n*32+j] < 0) break;
                        target[i] += ratios_vec[n*32+j] * source[ elems_vec[n*32+j] ];
                    }
                }
            }
        }
#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }


    void inject_field(const int_vec &idx,
                      const int_vec &is_changed,
                      const int_vec &idx_changed,
                      const int_vec &elems_vec,
                      const double_vec &ratios_vec,
                      const tensor_t &source,
                      tensor_t &target)
    {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif

        #pragma acc serial async
        int ntarget = target.size();

#ifndef ACC
        #pragma omp parallel default(none) shared(idx, source, \
            target, is_changed, idx_changed, elems_vec, ratios_vec, ntarget)
#endif
        {
#ifndef ACC
            #pragma omp for
#endif
            #pragma acc parallel loop async
            for (int i=0; i<ntarget; i++) {
                int n = idx[i];
                for (int d=0; d<NSTR; d++) {
                    target[i][d] = source[n][d];
                }
            }

#ifndef ACC
            #pragma omp for
#endif
            #pragma acc parallel loop async
            for (int i=0; i<ntarget; i++) {
                if (is_changed[i]>0) {
                    int n = idx_changed[i];

                for (int d=0; d<NSTR; d++) {
                    target[i][d] = 0;
                    for (int j=0; j<32; j++) {
                        if (elems_vec[n*32+j] < 0) break;
                        target[i][d] += ratios_vec[n*32+j] * source[ elems_vec[n*32+j] ][d];
                    }
                }
                }
            }
        }
#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }


    void nn_interpolate_elem_fields(Variables &var,
                                    const int_vec &idx,
                                    const int_vec &is_changed,
                                    const int_vec &idx_changed,
                                    const int_vec &elems_vec,
                                    const double_vec &ratios_vec)
    {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif
        const int e = var.nelem;

        double_vec *new_plstrain = new double_vec(e);
        inject_field(idx, is_changed, idx_changed, elems_vec, ratios_vec, *var.plstrain, *new_plstrain);

        double_vec *new_delta_pls = new double_vec(e);
        inject_field(idx, is_changed, idx_changed, elems_vec, ratios_vec, *var.delta_plstrain, *new_delta_pls);

        double_vec *new_eff_pls = new double_vec(e);
        inject_field(idx, is_changed, idx_changed, elems_vec, ratios_vec, *var.eff_pls, *new_eff_pls);

        tensor_t *new_strain = new tensor_t(e);
        inject_field(idx, is_changed, idx_changed, elems_vec, ratios_vec, *var.strain, *new_strain);

        tensor_t *new_stress = new tensor_t(e);
        inject_field(idx, is_changed, idx_changed, elems_vec, ratios_vec, *var.stress, *new_stress);

        double_vec *new_stressyy = new double_vec(e);
        inject_field(idx, is_changed, idx_changed, elems_vec, ratios_vec, *var.stressyy, *new_stressyy);

        double_vec *new_old_mean_stress = new double_vec(e);
        inject_field(idx, is_changed, idx_changed, elems_vec, ratios_vec, *var.old_mean_stress, *new_old_mean_stress);

        double_vec *new_radiogenic_source = new double_vec(e);
        inject_field(idx, is_changed, idx_changed, elems_vec, ratios_vec, *var.radiogenic_source, *new_radiogenic_source);

        double_vec *new_edvacc_surf = new double_vec(e);
        inject_field(idx, is_changed, idx_changed, elems_vec, ratios_vec, *var.surfinfo.edvacc_surf, *new_edvacc_surf);

        #pragma acc wait

        delete var.plstrain;
        var.plstrain = new_plstrain;

        delete var.delta_plstrain;
        var.delta_plstrain = new_delta_pls;

        delete var.strain;
        var.strain = new_strain;

        delete var.stress;
        var.stress = new_stress;

        delete var.stressyy;
        var.stressyy = new_stressyy;

        delete var.old_mean_stress;
        var.old_mean_stress = new_old_mean_stress;

        delete var.radiogenic_source;
        var.radiogenic_source = new_radiogenic_source;

        delete var.surfinfo.edvacc_surf;
        var.surfinfo.edvacc_surf = new_edvacc_surf;

        // b = new tensor_t(e);
        // inject_field(idx, is_changed, elems_vec, ratios_vec, *var.stress_old, *b);
        // delete var.stress_old;
        // var.stress_old = b;


#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }

} // anonymous namespace


void nearest_neighbor_interpolation(const Param& param, Variables &var,
                                    const Barycentric_transformation &bary,
                                    const array_t &old_coord,
                                    const conn_t &old_connectivity)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    int_vec idx(var.nelem); // nearest element
    int_vec is_changed(var.nelem); // is the element changed during remeshing?
    int_vec idx_changed(var.nelem); 

    int_vec elems_vec;
    double_vec ratios_vec;
    
    prepare_interpolation(param, var, bary, old_coord, old_connectivity, idx, is_changed, idx_changed, elems_vec, ratios_vec);

    std::cout << "    Interpolating fields...\n";
    nn_interpolate_elem_fields(var, idx, is_changed, idx_changed, elems_vec, ratios_vec);

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}
