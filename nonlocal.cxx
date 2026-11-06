#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "constants.hpp"
#include "parameters.hpp"
#include "nonlocal.hpp"

namespace {

// Interaction cutoff in units of the length l: exp(-9) ~ 1e-4 of the central weight.
const double cutoff_in_length = 3.0;

// Element centroids and the CSR graph of normalized Gaussian weights within the cutoff.
void build_graph(const Param& param, Variables& var)
{
    NonlocalState& nl = var.nonlocal;
    const int nelem = var.nelem;
    const double l2 = param.mat.nonlocal_length * param.mat.nonlocal_length;
    const double rc = cutoff_in_length * param.mat.nonlocal_length;
    const double rc2 = rc * rc;

    #pragma acc wait   // coord and volume are written by device kernels

    array_t centroid(nelem);
    #pragma omp parallel for default(none) shared(var, centroid, nelem)
    for (int e = 0; e < nelem; ++e) {
        ConstConnAccessor conn = (*var.connectivity)[e];
        for (int d = 0; d < NDIMS; ++d) {
            double c = 0;
            for (int i = 0; i < NODES_PER_ELEM; ++i) c += (*var.coord)[conn[i]][d];
            centroid[e][d] = c / NODES_PER_ELEM;
        }
    }

    PointCloud cloud(centroid);
    NANOKDTree kdtree(NDIMS, cloud);
    kdtree.buildIndex();

    // Per-element neighbour lists first (counts are unknown), then packed into CSR.
    typedef nanoflann::ResultItem<NANOKDTree::IndexType, double> Item;
    std::vector<std::vector<Item>> found(nelem);
    #pragma omp parallel for default(none) shared(kdtree, centroid, found, nelem, rc2) schedule(dynamic, 64)
    for (int e = 0; e < nelem; ++e) {
        double x[NDIMS];   // centroid is SOA: copy the query point out
        for (int d = 0; d < NDIMS; ++d) x[d] = centroid[e][d];
        // L2_Simple_Adaptor: radius and returned distances are squared
        kdtree.radiusSearch(x, rc2, found[e]);
    }

    nl.idx.assign(nelem + 1, 0);
    for (int e = 0; e < nelem; ++e) nl.idx[e + 1] = nl.idx[e] + int(found[e].size());
    const int nnz = nl.idx[nelem];
    nl.arr.resize(nnz);
    nl.w.resize(nnz);

    #pragma omp parallel for default(none) shared(var, nl, found, nelem, l2)
    for (int e = 0; e < nelem; ++e) {
        const int start = nl.idx[e];
        double sum = 0;
        for (std::size_t k = 0; k < found[e].size(); ++k) {
            const int j = int(found[e][k].first);
            const double wv = std::exp(-found[e][k].second / l2) * (*var.volume)[j];
            nl.arr[start + k] = j;
            nl.w[start + k] = wv;
            sum += wv;
        }
        // the element itself is always found (r = 0), so sum > 0
        for (std::size_t k = 0; k < found[e].size(); ++k) nl.w[start + k] /= sum;
    }

    nl.built_nelem = nelem;
    nl.built_nremesh = var.nremesh;

    std::printf("  [nonlocal] graph: %d elements, %.1f neighbours per element "
                "(l = %g m, cutoff %g m)\n", nelem, double(nnz) / nelem,
                param.mat.nonlocal_length, rc);
}

} // anonymous namespace


void init_nonlocal_plstrain(const Param& param, Variables& var)
{
    if (param.mat.nonlocal_option == 0) return;
    #pragma acc wait
    *var.nonlocal.plstrain_nl = *var.plstrain;
}


void accumulate_nonlocal_plstrain(const Param& param, Variables& var)
{
    if (param.mat.nonlocal_option == 0) return;
    NonlocalState& nl = var.nonlocal;

    // The graph depends on the mesh only: build it at the first step and after each remesh.
    if (nl.built_nelem != var.nelem || nl.built_nremesh != var.nremesh)
        build_graph(param, var);

    #pragma acc wait   // delta_plstrain is written by update_stress on the device

    const int nelem = var.nelem;
    const double_vec& dpls = *var.delta_plstrain;
    double_vec& kbar = *nl.plstrain_nl;
    #pragma omp parallel for default(none) shared(nl, dpls, kbar, nelem)
    for (int e = 0; e < nelem; ++e) {
        double a = 0;
        for (int k = nl.idx[e]; k < nl.idx[e + 1]; ++k) a += nl.w[k] * dpls[nl.arr[k]];
        kbar[e] += a;
    }
}
