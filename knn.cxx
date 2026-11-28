#include <algorithm>
#include <iostream>
#include <vector>
#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "parameters.hpp"
#include "knn.hpp"
#include "utils.hpp"


void KNN::print_search(const char* stage, long nquery, int k_neig,
        [[maybe_unused]] int batch_queries) const
{
    printf("%sknn: %10s pts for %11s q w/ %2d k", stage,
            format_with_commas((unsigned long)numPoints).c_str(),
            format_with_commas((unsigned long)nquery).c_str(), k_neig);
#if defined(ACC) && !defined(USEMPI)
    // d_queries + d_results of one batch
    const size_t bytes = (size_t)batch_queries * (3 * sizeof(float) + (size_t)k_neig * sizeof(neighbor));
    printf(", %6.1f MB", bytes / 1048576.0);
#endif
    fflush(stdout);
}

KNN::KNN(const array_t& points_vec_, NANOKDTree& nano_kdtree_, [[maybe_unused]] int capacity) :
    points_vec(points_vec_), numPoints(points_vec_.size()),
    nano_kdtree(nano_kdtree_)
{
    h_results = nullptr;
    h_results_capacity = 0;

#if defined(ACC) && !defined(USEMPI)
    m_bvh_state = nullptr;

    std::vector<float3> h_points(numPoints);
    points_vec_.pack_to_xyz_float(h_points);
    if (capacity > numPoints) {
        m_bvh_state = knn_bvh_create_padded(
            reinterpret_cast<const float*>(h_points.data()), numPoints, capacity);
    } else {
        m_bvh_state = knn_bvh_create(
            reinterpret_cast<const float*>(h_points.data()), numPoints);
    }

#else
    // Mesh management in GPU MPI uses the same double-precision search as the
    // CPU solver. Float BVH coordinates can change marker and remap identities.
    #pragma acc wait
    nano_kdtree.buildIndex();
#endif
}

KNN::~KNN()
{
    if (h_results_capacity > 0) {
        delete[] h_results;
        h_results = nullptr;
    }
#if defined(ACC) && !defined(USEMPI)
    knn_bvh_destroy(static_cast<KNNBVHState*>(m_bvh_state));
#endif
}

int KNN::max_batch_size([[maybe_unused]] int k_neig) const
{
#if defined(ACC) && !defined(USEMPI)
    int batch_size = knn_bvh_max_batch(static_cast<KNNBVHState*>(m_bvh_state), k_neig);
    if (batch_size <= 0) {
        throw std::runtime_error("Insufficient GPU memory to process even a single query batch.");
    }
    return batch_size;
#else
    return 16 * 1024 * 1024;
#endif
}

neighbor* KNN::search(const array_t& queries, int nquery, int k_neig,
        bool is_sync_to_host, const float* d_guess_radii_sq)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif

    const size_t required_size = (size_t)nquery * k_neig;

#if defined(ACC) && !defined(USEMPI)
    if (is_sync_to_host && required_size > h_results_capacity) {
        if (h_results) delete[] h_results;
        h_results = new neighbor[required_size];
        h_results_capacity = required_size;
    }

    if (is_sync_to_host) {
        // CPU-built queries (brc-interpolation). Keep existing H2D-copy path.
        std::vector<float3> h_queries_v(nquery);
        queries.pack_to_xyz_float(h_queries_v, nquery);
#ifdef NPROF_DETAIL
        nvtxRangePop();
#endif
        return static_cast<neighbor*>(
            knn_bvh_search(
                static_cast<KNNBVHState*>(m_bvh_state),
                reinterpret_cast<const float*>(h_queries_v.data()),
                nquery, k_neig, h_results, is_sync_to_host,
                d_guess_radii_sq));
    }

    // GPU-direct path: queries are GPU-resident (managed memory).
    float3* d_qptr =
        knn_bvh_ensure_queries(static_cast<KNNBVHState*>(m_bvh_state), nquery);
    {
        #pragma acc parallel loop gang vector deviceptr(d_qptr)
        for (int i = 0; i < nquery; i++) {
            d_qptr[i].x = (float)queries[i][0];
            d_qptr[i].y = (float)queries[i][1];
#ifdef THREED
            d_qptr[i].z = (float)queries[i][2];
#else
            d_qptr[i].z = 0.0f;
#endif
        }
    }
    // Synchronous acc parallel above already flushes the default queue,
    // but explicit wait makes the OpenACC->CUDA handoff unambiguous.
    #pragma acc wait

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
    return static_cast<neighbor*>(
        knn_bvh_search_prepared(
            static_cast<KNNBVHState*>(m_bvh_state),
            nquery, k_neig,
            h_results,       // ignored when is_sync_to_host=false
            is_sync_to_host,
            d_guess_radii_sq));

#else
    #pragma acc wait

    if (required_size > h_results_capacity) {
        if (h_results) delete[] h_results;
        h_results = new neighbor[required_size];
        h_results_capacity = required_size;
    }

    // Pre-allocate per-thread scratch buffers to avoid heap contention in the loop.
#ifdef _OPENMP
    const int nthreads = omp_get_max_threads();
#else
    const int nthreads = 1;
#endif
    std::vector<size_t> nn_idx_buf(nthreads * k_neig);
    std::vector<double> dist_buf(nthreads * k_neig);

    #pragma omp parallel for default(none) \
        shared(queries, h_results, k_neig, nano_kdtree, nquery, nn_idx_buf, dist_buf)
    for (int i = 0; i < nquery; ++i) {
        neighbor *result = h_results + (size_t)i * k_neig;

#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        size_t* nn_idx        = nn_idx_buf.data() + tid * k_neig;
        double* out_dists_sqr = dist_buf.data()   + tid * k_neig;

        KNNResultSet resultSet(k_neig);
        resultSet.init(nn_idx, out_dists_sqr);

        double q[NDIMS];
        queries[i].copy_to(q);
        nano_kdtree.findNeighbors(resultSet, q);

        for (int j = 0; j < k_neig; ++j) {
            result[j].idx   = nn_idx[j];
            result[j].dist2 = out_dists_sqr[j];
        }
    }

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
    return h_results;

#endif
}
