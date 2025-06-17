#ifdef USE_NPROF
#include <nvToolsExt.h> 
#endif
#include "algorithm"
#include "iostream"

#include "constants.hpp"
#include "parameters.hpp"

#include "barycentric-fn.hpp"
#include "utils.hpp"
#include "brc-interpolation.hpp"

namespace { // anonymous namespace

typedef Array2D<double,NODES_PER_ELEM> brc_t;


void interpolate_field(const brc_t &brc, const int_vec &el, const conn_t &connectivity,
                       const double_vec &source, double_vec &target)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    #pragma acc serial async
    int ntarget = target.size();

#ifndef ACC
    #pragma omp parallel for default(none)          \
        shared(brc, el, connectivity, source, target,ntarget)
#endif
    #pragma acc parallel loop async
    for (int i=0; i<ntarget; i++) {
        int e = el[i];
        const int *conn = connectivity[e];
        double result = 0;
        for (int j=0; j<NODES_PER_ELEM; j++) {
            result += source[conn[j]] * brc[i][j];
        }
        target[i] = result;
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void interpolate_field(const brc_t &brc, const int_vec &el, const conn_t &connectivity,
                       const array_t &source, array_t &target)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    #pragma acc serial async
    int ntarget = target.size();

#ifndef ACC
    #pragma omp parallel for default(none)          \
        shared(brc, el, connectivity, source, target,ntarget)
#endif
    #pragma acc parallel loop async
    for (int i=0; i<ntarget; i++) {
        int e = el[i];
        const int *conn = connectivity[e];
        for (int d=0; d<NDIMS; d++) {
            double result = 0;
            for (int j=0; j<NODES_PER_ELEM; j++) {
                result += source[conn[j]][d] * brc[i][j];
            }
            target[i][d] = result;
        }
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void prepare_interpolation(const Variables &var,
                           const Barycentric_transformation &bary,
                           const array_t &old_coord,
                           const conn_t &old_connectivity,
                           const int_vec2D &old_support,
                           brc_t &brc, int_vec &el)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // for each new coord point, find the enclosing old element

#ifdef USE_NPROF
    nvtxRangePushA("create kdtree for coord");
#endif
    PointCloud cloud(old_coord);
    KDTree kdtree(NDIMS, cloud);
    kdtree.buildIndex();
#ifdef USE_NPROF
    nvtxRangePop();
#endif

    const int k = 1;

    #pragma omp parallel for default(none)          \
        shared(var, bary, old_coord, old_connectivity, old_support, kdtree, el, brc) \
        firstprivate(k)
    for (int i=0; i<var.nnode; i++) {
        double *q = (*var.coord)[i];
        size_t_vec nn_idx(k);
        double_vec dd(k);
        KNNResultSet resultSet(k);
        resultSet.init(nn_idx.data(), dd.data());

        // find the nearest point nn in old_coord
        kdtree.findNeighbors(resultSet, q);

        int nn = nn_idx[0];

        // elements surrounding nn
        const int_vec &nn_elem = old_support[nn];

        // std::cout << i << " ";
        // print(std::cout, q, NDIMS);
        // std::cout << " " << nn << " " << dd[0] << '\n';

        double r[NDIMS];
        int e;

        // shortcut: q is exactly the same as nn
        if (dd[0] == 0) {
            e = nn_elem[0];
            bary.transform(q, e, r);
            // r should be a permutation of [1, 0, 0]
            // normalize r to remove round-off error
            for (int d=0; d<NDIMS; d++) {
                if (r[d] > 0.9)
                    r[d] = 1;
                else
                    r[d] = 0;
            }
            goto found;
        }

        // loop over (old) elements surrounding nn to find
        // the element that is enclosing q
        for (std::size_t j=0; j<nn_elem.size(); j++) {
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
             */

            // this array contains the elements that have been searched so far
            int_vec searched(nn_elem);

            // search through elements that are neighbors of nn_elem
            for (std::size_t j=0; j<nn_elem.size(); j++) {
                int ee = nn_elem[j];
                const int *conn = old_connectivity[ee];
                for (int m=0; m<NODES_PER_ELEM; m++) {
                    // np is a node close to q
                    int np = conn[m];
                    const int_vec &np_elem = old_support[np];
                    for (std::size_t j=0; j<np_elem.size(); j++) {
                        e = np_elem[j];
                        auto it = std::find(searched.begin(), searched.end(), e);
                        if (it != searched.end()) {
                            // this element has been searched before
                            continue;
                        }
                        searched.push_back(e);
                        bary.transform(q, e, r);
                        // std::cout << e << " ";
                        // print(std::cout, r, NDIMS);
                        // std::cout << " ... \n";
                        if (bary.is_inside(r)) {
                            goto found;
                        }
                    }
                }
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

    // print(std::cout, *var.coord);
    // std::cout << '\n';
    // print(std::cout, el);
    // std::cout << '\n';
    // print(std::cout, bar);
    // std::cout << '\n';
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

} // anonymous namespace

void barycentric_node_interpolation(Variables &var,
                                    const Barycentric_transformation &bary,
                                    const array_t &old_coord,
                                    const conn_t &old_connectivity)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    int_vec el(var.nnode);
    brc_t brc(var.nnode);
    prepare_interpolation(var, bary, old_coord, old_connectivity, *var.support, brc, el);

    double_vec *new_temperature = new double_vec(var.nnode);
    interpolate_field(brc, el, old_connectivity, *var.temperature, *new_temperature);

    double_vec *new_ppressure = new double_vec(var.nnode);
    interpolate_field(brc, el, old_connectivity, *var.ppressure, *new_ppressure);

    double_vec *new_dppressure = new double_vec(var.nnode);
    interpolate_field(brc, el, old_connectivity, *var.dppressure, *new_dppressure);

    array_t *new_vel = new array_t(var.nnode);
    interpolate_field(brc, el, old_connectivity, *var.vel, *new_vel);

    array_t *new_coord0 = new array_t(var.nnode);
    interpolate_field(brc, el, old_connectivity, *var.coord0, *new_coord0);

    #pragma acc wait

    delete var.temperature;
    var.temperature = new_temperature;

    delete var.ppressure;
    var.ppressure = new_ppressure;

    delete var.dppressure;
    var.dppressure = new_dppressure;

    delete var.vel;
    var.vel = new_vel;

    delete var.coord0;
    var.coord0 = new_coord0;
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void barycentric_node_interpolation_forT(const Variables &var,
                                         const Barycentric_transformation &bary,
                                         const array_t &input_coord,
                                         const conn_t &input_connectivity,
                                         const int_vec2D &input_support,
					 const double_vec &inputtemperature,
					 double_vec &outputtemperature)
{
    int_vec el(var.nnode);
    brc_t brc(var.nnode);
    prepare_interpolation(var, bary, input_coord, input_connectivity, input_support, brc, el);

    interpolate_field(brc, el, input_connectivity, inputtemperature, outputtemperature);
}
