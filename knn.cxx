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



KNN::KNN(const Param& param, const array_t& points_vec_, NANOKDTree& nano_kdtree_, int capacity) :
    resolution(param.mesh.resolution), points_vec(points_vec_),
    numPoints(points_vec_.size()),
    nano_kdtree(nano_kdtree_)
{
    h_results = nullptr;
    h_results_capacity = 0;

#ifdef ACC
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
    nano_kdtree.buildIndex();
#endif
}

KNN::~KNN()
{
    if (h_results_capacity > 0) {
        delete[] h_results;
        h_results = nullptr;
    }
#ifdef ACC
    knn_bvh_destroy(static_cast<KNNBVHState*>(m_bvh_state));
#endif
}

int KNN::max_batch_size(int k_neig) const
{
#ifdef ACC
    int batch_size = knn_bvh_max_batch(static_cast<KNNBVHState*>(m_bvh_state), k_neig);
    if (batch_size <= 0) {
        throw std::runtime_error("Insufficient GPU memory to process even a single query batch.");
    }
    return batch_size;
#else
    (void)k_neig;
    return 16 * 1024 * 1024;
#endif
}

neighbor* KNN::search(const array_t& queries, int nquery, int k_neig,
        bool is_sync_to_host, const float* d_guess_radii_sq)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    printf("      Running knn query on %d points ", numPoints);

    const size_t required_size = (size_t)nquery * k_neig;

#ifdef ACC
    printf("(lbvh GPU)\n");

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
    printf("(nano-kdtree)\n");

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
