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
#include "nn-interpolation.hpp"


namespace {

    void find_nearest_neighbor(const Variables &var, KDTree &kdtree,
                               int_vec &idx, int_vec &is_changed)
    {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif
        array_t *new_center = elem_center(*var.coord, *var.connectivity);

        const double eps = 1e-15;

        #pragma omp parallel for default(none) shared(var, kdtree, new_center, idx, is_changed)
        for(int e=0; e<var.nelem; e++) {

            double *q = (*new_center)[e];
            size_t nn_idx;
            double dd;

            KNNResultSet resultSet(1);
            resultSet.init(&nn_idx, &dd);  
            kdtree.findNeighbors(resultSet, q);

            idx[e] = static_cast<int>(nn_idx);
            is_changed[e] = (dd < eps) ? 0 : 1;
        }

        delete new_center;
#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }

    void find_acm_elem_ratios(const Variables &var,
                              const Barycentric_transformation &bary,
                              int_vec &is_changed,
                              KDTree &kdtree,
                              int old_nelem,
                              int_vec &elems_vec,
                              double_vec &ratios_vec,
                              int_vec &idx_changed)
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

        int nelem_changed = 0;

        for (int e=0; e<var.nelem; e++) {
            if (is_changed[e]) {
                idx_changed[e] = nelem_changed;
                nelem_changed++;
            }
        }

        elems_vec.resize(nelem_changed*32,-1);
        ratios_vec.resize(nelem_changed*32);

        #pragma omp parallel for default(none) shared(var, bary, is_changed, kdtree, elems_vec, ratios_vec, idx_changed) firstprivate(max_el)
        for(int e=0; e<var.nelem; e++) {
            if (is_changed[e]) {
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

                            double x[NDIMS] = {0}; // coordinate of temporary point
                            for (int d=0; d<NDIMS; d++)
                                for (int n=0; n<NODES_PER_ELEM; n++) {
                                    x[d] += (*var.coord)[ conn[n] ][d] * eta[n];
                                }

                            // find the nearest point nn in old_center
                            size_t_vec nn_idx(max_el);
                            double_vec out_dists_sqr(max_el);
                            KNNResultSet resultSet(max_el);
                            resultSet.init(nn_idx.data(), out_dists_sqr.data());

                            kdtree.findNeighbors(resultSet, x);

                            // std::cout << "  ";
                            // print(std::cout, eta, NODES_PER_ELEM);
                            // print(std::cout, x, NDIMS);

                            // find the old element that is enclosing x
                            double r[NDIMS];
                            int old_e;
                            for (int jj=0; jj<max_el; jj++) {
                                old_e = nn_idx[jj];
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
#ifdef THREED
                        }
#endif
                    }

                // Count
                int total_count = 0;
                for (int i=0; i<=elem_size; ++i)
                    total_count += elem_count_buf[i];

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
                for (int i=0; i<elem_size; ++i) {
                    elems_vec[idx_changed[e]*32+i] = elem_keys[i];
                    ratios_vec[idx_changed[e]*32+i] = elem_count_buf[i] * inv;
                }
            }
        }

#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }


    void prepare_interpolation(Variables &var,
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
        array_t *old_center = elem_center(old_coord, old_connectivity);
        int old_nelem = old_connectivity.size();

#ifdef USE_NPROF
        nvtxRangePushA("create kdtree for old elements");
#endif
        PointCloud cloud(*old_center);
        KDTree kdtree(NDIMS, cloud);
        kdtree.buildIndex();
#ifdef USE_NPROF
        nvtxRangePop();
#endif

        find_nearest_neighbor(var, kdtree, idx, is_changed);

        find_acm_elem_ratios(var, bary, is_changed, kdtree, old_nelem, elems_vec, ratios_vec, idx_changed);

        delete old_center;
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
        #pragma omp parallel default(none) shared(idx, source, target, is_changed, idx_changed, ratios_vec, elems_vec)
        {
            #pragma omp for
            for (std::size_t i=0; i<target.size(); i++) {
                int n = idx[i];
                target[i] = source[n];
            }

            #pragma omp for
            for (std::size_t i=0; i<target.size(); i++) {
                if (is_changed[i]>0) {
                    int n = idx_changed[i];

                    target[i] = 0;
                    for (std::size_t j=0; j<32; j++) {
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
        #pragma omp parallel default(none) shared(idx, source, target, is_changed, idx_changed, elems_vec, ratios_vec)
        {
            #pragma omp for
            for (std::size_t i=0; i<target.size(); i++) {
                int n = idx[i];
                for (int d=0; d<NSTR; d++) {
                    target[i][d] = source[n][d];
                }
            }

            #pragma omp for
            for (std::size_t i=0; i<target.size(); i++) {
                if (is_changed[i]>0) {
                    int n = idx_changed[i];

                    for (int d=0; d<NSTR; d++) {
                        target[i][d] = 0;
                        for (std::size_t j=0; j<32; j++) {
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
        const int n = var.nnode;
        const int e = var.nelem;

        double_vec *a;

        a = new double_vec(e);
        inject_field(idx, is_changed, idx_changed, elems_vec, ratios_vec, *var.plstrain, *a);
        delete var.plstrain;
        var.plstrain = a;

        a = new double_vec(e);
        inject_field(idx, is_changed, idx_changed, elems_vec, ratios_vec, *var.delta_plstrain, *a);
        delete var.delta_plstrain;
        var.delta_plstrain = a;

        tensor_t *b;
        b = new tensor_t(e);
        inject_field(idx, is_changed, idx_changed, elems_vec, ratios_vec, *var.strain, *b);
        delete var.strain;
        var.strain = b;

        b = new tensor_t(e);
        inject_field(idx, is_changed, idx_changed, elems_vec, ratios_vec, *var.stress, *b);
        delete var.stress;
        var.stress = b;

        a = new double_vec(e);
        inject_field(idx, is_changed, idx_changed, elems_vec, ratios_vec, *var.stressyy, *a);
        delete var.stressyy;
        var.stressyy = a;

        a = new double_vec(e);
        inject_field(idx, is_changed, idx_changed, elems_vec, ratios_vec, *var.old_mean_stress, *a);
        delete var.old_mean_stress;
        var.old_mean_stress = a;

        // b = new tensor_t(e);
        // inject_field(idx, is_changed, elems_vec, ratios_vec, *var.stress_old, *b);
        // delete var.stress_old;
        // var.stress_old = b;

        a = new double_vec(e);
        inject_field(idx, is_changed, idx_changed, elems_vec, ratios_vec, *var.radiogenic_source, *a);
        delete var.radiogenic_source;
        var.radiogenic_source = a;

        a = new double_vec(e);
        inject_field(idx, is_changed, idx_changed, elems_vec, ratios_vec, *var.surfinfo.edvacc_surf, *a);
        delete var.surfinfo.edvacc_surf;
        var.surfinfo.edvacc_surf = a;

        a = new double_vec(e);
        inject_field(idx, is_changed, elems_vec, ratios_vec, *var.surfinfo.edhacc_oc, *a);
        delete var.surfinfo.edhacc_oc;
        var.surfinfo.edhacc_oc = a;


#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }

} // anonymous namespace


void nearest_neighbor_interpolation(Variables &var,
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

    prepare_interpolation(var, bary, old_coord, old_connectivity, idx, is_changed, idx_changed, elems_vec, ratios_vec);

    nn_interpolate_elem_fields(var, idx, is_changed, idx_changed, elems_vec, ratios_vec);

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}
