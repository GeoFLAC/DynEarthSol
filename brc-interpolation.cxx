#include "algorithm"
#include "iostream"
#include "unordered_set"

#include "constants.hpp"
#include "parameters.hpp"

#include "barycentric-fn.hpp"
#include "utils.hpp"
#include "knn.hpp"
#include "brc-interpolation.hpp"

namespace { // anonymous namespace

typedef Array2D<double,NODES_PER_ELEM> brc_t;


void interpolate_field(const brc_t &brc, const int_vec &el, const conn_t &connectivity,
                       const double_vec &source, double_vec &target, int ntarget)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif

#ifndef ACC
    #pragma omp parallel for default(none)          \
        shared(brc, el, connectivity, source, target,ntarget)
#endif
    #pragma acc parallel loop gang vector async
    for (int i=0; i<ntarget; i++) {
        int e = el[i];
        ConstConnAccessor conn = connectivity[e];
        double result = 0;
        for (int j=0; j<NODES_PER_ELEM; j++) {
            result += source[conn[j]] * brc[i][j];
        }
        target[i] = result;
    }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}


void interpolate_field(const brc_t &brc, const int_vec &el, const conn_t &connectivity,
                       const array_t &source, array_t &target, int ntarget)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif

#ifndef ACC
    #pragma omp parallel for default(none)          \
        shared(brc, el, connectivity, source, target,ntarget)
#endif
    #pragma acc parallel loop gang vector async
    for (int i=0; i<ntarget; i++) {
        int e = el[i];
        ConstConnAccessor conn = connectivity[e];
        for (int d=0; d<NDIMS; d++) {
            double result = 0;
            for (int j=0; j<NODES_PER_ELEM; j++) {
                result += source[conn[j]][d] * brc[i][j];
            }
            target[i][d] = result;
        }
    }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}


void interpolate_field(const brc_t &brc, const int_vec &el, const conn_t &connectivity,
                       const tensor_t &source, tensor_t &target, int ntarget)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif

#ifndef ACC
    #pragma omp parallel for default(none) \
        shared(brc, el, connectivity, source, target, ntarget)
#endif
    #pragma acc parallel loop gang vector async
    for (int i = 0; i < ntarget; i++) {
        int e = el[i];
        ConstConnAccessor conn = connectivity[e];
        for (int d = 0; d < NSTR; d++) {
            double result = 0;
            for (int j = 0; j < NODES_PER_ELEM; j++)
                result += source[conn[j]][d] * brc[i][j];
            target[i][d] = result;
        }
    }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}


void prepare_interpolation(const Param& param, const Variables &var,
                           const Barycentric_transformation &bary,
                           const array_t &old_coord,
                           const conn_t &old_connectivity,
                           const Support &old_support,
                           brc_t &brc, int_vec &el)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    // for each new coord point, find the enclosing old element

#ifdef NPROF_DETAIL
    nvtxRangePush("create kdtree for coord");
#endif

#ifdef ACC
    array_t point_tmp(1);
    PointCloud cloud(point_tmp);
#else
    PointCloud cloud(old_coord);
#endif

    NANOKDTree nano_kdtree(NDIMS, cloud);
    KNN kdtree(param, old_coord, nano_kdtree, false);

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif

    int max_el = 1;
    int nodes_per_block_max = std::max(kdtree.max_batch_size(max_el), 1);

    int nodes_per_block = std::min(nodes_per_block_max, var.nnode);
    int nblocks = (var.nnode + nodes_per_block - 1) / nodes_per_block;

    {
        char knn_be[64]; kdtree.backend_str(knn_be, sizeof(knn_be));
        char mem_str[32] = "";
        size_t qm = KNN::query_mem_bytes(nodes_per_block, max_el);
        if (qm) snprintf(mem_str, sizeof(mem_str), ", %6.1f MB", qm / 1048576.0);

        printf("    Node:             knn: %10s pts (%s) for %11s q w/ %2d k%s",
                format_with_commas((unsigned long)kdtree.get_npoints()).c_str(), knn_be,
                format_with_commas((unsigned long)var.nnode).c_str(), max_el, mem_str
                ); fflush(stdout);

        if (nblocks > 1) {
            printf("/blk (%d x %sq)", nblocks,
                    format_with_commas((unsigned long)nodes_per_block).c_str());
            fflush(stdout);
        } else {
            printf("\n");
        }
    }

    array_t block_queries;

    // bary's coeff_ is filled by an async ACC kernel; it is read on the host below
    #pragma acc wait

    for (int b = 0; b < nblocks; b++) {
        int start = b * nodes_per_block;
        int end = std::min((b + 1) * nodes_per_block, var.nnode);
        if (start >= end) continue;

        if (nblocks > 1) {printf(" %d", b+1); fflush(stdout);}
        
        block_queries.resize(end - start);
        #pragma omp parallel for default(none) shared(var, block_queries, start, end)
        for (int i = start; i < end; i++) {
            ArrayAccessor q_dst = block_queries[i - start];
            ConstArrayAccessor q_src = (*var.coord)[i];
            for (int d = 0; d < NDIMS; d++) {
                q_dst[d] = q_src[d];
            }
        }

        neighbor* neighbors = kdtree.search(block_queries, (end - start), max_el);

        #pragma omp parallel for default(none) schedule(guided) \
            shared(var, bary, old_coord, old_connectivity, old_support, el, brc, neighbors, start, end)
        for (int i = start; i < end; i++) {
            
            int local_i = i - start;
            
            ConstArrayAccessor q = (*var.coord)[i];

            int nn = neighbors[local_i].idx;
            double dd = neighbors[local_i].dist2;

            // elements surrounding nn
            const int nn_nelem = old_support.size(nn);
            const int* nn_elem = old_support.patch(nn);

        // std::cout << i << " ";
        // print(std::cout, q, NDIMS);
        // std::cout << " " << nn << " " << dd[0] << '\n';

            double r[NDIMS];
            int e;

            // shortcut: q is exactly the same as nn
            if (dd == 0) {
                e = nn_elem[0];
                bary.transform(q, e, r);
                // r should be a permutation of [1, 0, 0]
                // normalize r to remove round-off error
                for (int d=0; d<NDIMS; d++)
                    r[d] = (r[d] > 0.9) ? 1 : 0;

                goto found;
            }

            // loop over (old) elements surrounding nn to find
            // the element that is enclosing q
            for (int j=0; j<nn_nelem; j++) {
                e = nn_elem[j];
                bary.transform(q, e, r);
                if (bary.is_inside(r)) {
                // std::cout << e << " ";
                // print(std::cout, r, NDIMS);
                // std::cout << '\n';
                    goto found;
                }
            }

            /* not_found */

            {
                /* Situation: q is in the upper element, but its nearest point is o!
                * we won't find the enclosing element with the method above
                *     x
                *    / \   <-- this is a large triangle
                *   / q                            \
                *  x---- x
                *   \-o-/   <-- this is a small triangle
                *
                * True Breadth-First Search (BFS) outward from nn's element support, 
                * expanding level by level.
                * The one-level expansion above fails when q lies in a large old element
                * whose nearest node nn is not one of its vertices (common after aggressive
                * 3D remeshing that places new nodes far from any existing node).
                * BFS visits neighbors level-by-level (capped by MAX_LAYERS below) and will
                * typically find the enclosing element even after aggressive remeshing. If
                * the capped search does not find q, it may be outside the old domain and
                * will fall through to the nearest-node fallback below.
                */
                std::unordered_set<int> visited(nn_elem, nn_elem + nn_nelem);
                int_vec frontier(nn_elem, nn_elem + nn_nelem);
                // Track the element where q is least outside (highest min bary coord).
                // Used post-BFS to distinguish boundary-face vs. clearly-outside vs. bug.
                int    best_e_bfs    = (nn_nelem == 0) ? 0 : nn_elem[0];
                double best_min_bary = -1e30;

                const int MAX_LAYERS = 3;
                for (int layer = 0; layer < MAX_LAYERS && !frontier.empty(); layer++) {
                    int_vec next_frontier;
                    for (int ee : frontier) {
                        ConstConnAccessor conn = old_connectivity[ee];
                        for (int m=0; m<NODES_PER_ELEM; m++) {
                            // np is a node close to q
                            int np = conn[m];
                            const int np_nelem = old_support.size(np);
                            const int* np_elem = old_support.patch(np);
                            for (int jj=0; jj<np_nelem; jj++) {
                                int candidate = np_elem[jj];
                                if (!visited.insert(candidate).second) continue;
                                bary.transform(q, candidate, r);
                                // min bary coord = how far inside (negative means outside)
                                double min_bary = r[0];
                                for (int d = 1; d < NDIMS; d++)
                                    if (r[d] < min_bary) min_bary = r[d];
                                double last = 1.0;
                                for (int d = 0; d < NDIMS; d++) last -= r[d];
                                if (last < min_bary) min_bary = last;
                                if (min_bary > best_min_bary) {
                                    best_min_bary = min_bary;
                                    best_e_bfs = candidate;
                                }
                                if (bary.is_inside(r)) {
                                    e = candidate;
                                    goto found;
                                }
                                next_frontier.push_back(candidate);
                            }
                        }
                    }
                    frontier = std::move(next_frontier);
                }
                // BFS exhausted without finding q.
                // best_min_bary > -1e-10: q is numerically on a face (FP tolerance
                //   missed it) — use best_e_bfs directly, no fallback needed.
                if (best_min_bary > -1e-10) {
                    e = best_e_bfs;
                    bary.transform(q, e, r);
                    goto found;
                }
                // For all other cases, distinguish by bcflag of the new node:
                //   bcflag == 0: interior node — might be inside old domain, BFS failure might be a bug.
                //   bcflag != 0: boundary node — may sit just outside the old mesh after
                //                remeshing; nearest-node fallback is correct.
                if ((*var.bcflag)[i] == 0) {
                    printf("Warning: prepare_interpolation: interior node %d (bcflag=0) "
                           "not found after capped BFS (MAX_LAYERS=%d; %zu/%d elements searched), "
                           "best_min_bary=%.3e. Mesh connectivity may be broken.\n",
                           i, MAX_LAYERS, visited.size(), (int)old_connectivity.size(), best_min_bary);
                }
            }
            {
            //std::cout << "New node is outside of the old domain. \n";

                // Situation: q must be outside the old domain
                // using nearest old_coord instead
                e = nn_elem[0];
                bary.transform(old_coord[nn], e, r);
            }
        found:
            el[i] = e;
            double sum = 0;
            for (int d=0; d<NDIMS; d++) {
                brc[i][d] = r[d];
                sum += r[d];
            }
            brc[i][NODES_PER_ELEM-1] = 1 - sum;
        }
    }

    if (nblocks > 1) printf("\n");

    // print(std::cout, *var.coord);
    // std::cout << '\n';
    // print(std::cout, el);
    // std::cout << '\n';
    // print(std::cout, bar);
    // std::cout << '\n';
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

} // anonymous namespace

void barycentric_node_interpolation(const Param& param, Variables &var,
                                    const Barycentric_transformation &bary,
                                    const array_t &old_coord,
                                    const conn_t &old_connectivity)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    const int n = var.nnode;

    int_vec el(n);
    brc_t brc(n);
    prepare_interpolation(param, var, bary, old_coord, old_connectivity, var.support, brc, el);

    double_vec *new_temperature = new double_vec(n);
    interpolate_field(brc, el, old_connectivity, *var.temperature, *new_temperature, n);

    double_vec *new_ppressure = new double_vec(n);
    interpolate_field(brc, el, old_connectivity, *var.ppressure, *new_ppressure, n);

    double_vec *new_dppressure = new double_vec(n);
    interpolate_field(brc, el, old_connectivity, *var.dppressure, *new_dppressure, n);

#ifdef USEMMG
    double_vec *new_init_elem_size_n = new double_vec(n);
    interpolate_field(brc, el, old_connectivity, *var.init_elem_size_n, *new_init_elem_size_n, n);
#endif

    #pragma acc wait

    array_t *new_vel = new array_t(n);
    interpolate_field(brc, el, old_connectivity, *var.vel, *new_vel, n);

    array_t *new_coord0 = new array_t(n);
    interpolate_field(brc, el, old_connectivity, *var.coord0, *new_coord0, n);

    // Null when remesh() switched the SPR chain off: nothing to carry across.
    tensor_t *new_stress_n = nullptr;
    if (var.stress_n) {
        new_stress_n = new tensor_t(n);
        interpolate_field(brc, el, old_connectivity, *var.stress_n, *new_stress_n, n);
    }

    double_vec *new_stressyy_n = nullptr;
    if (var.stressyy_n) {
        new_stressyy_n = new double_vec(n);
        interpolate_field(brc, el, old_connectivity, *var.stressyy_n, *new_stressyy_n, n);
    }

    delete var.temperature;
    var.temperature = new_temperature;

    delete var.ppressure;
    var.ppressure = new_ppressure;

    delete var.dppressure;
    var.dppressure = new_dppressure;

#ifdef USEMMG
    delete var.init_elem_size_n;
    var.init_elem_size_n = new_init_elem_size_n;
#endif

    #pragma acc wait

    delete var.vel;
    var.vel = new_vel;

    delete var.coord0;
    var.coord0 = new_coord0;

    delete var.stress_n;
    var.stress_n = new_stress_n;

    delete var.stressyy_n;
    var.stressyy_n = new_stressyy_n;

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}


void barycentric_node_interpolation_forT(const Param& param, const Variables &var,
                                         const Barycentric_transformation &bary,
                                         const array_t &input_coord,
                                         const conn_t &input_connectivity,
                                         const Support &input_support,
					 const double_vec &inputtemperature,
					 double_vec &outputtemperature)
{
    int_vec el(var.nnode);
    brc_t brc(var.nnode);
    prepare_interpolation(param, var, bary, input_coord, input_connectivity, input_support, brc, el);

    interpolate_field(brc, el, input_connectivity, inputtemperature, outputtemperature, var.nnode);

    #pragma acc wait
}
