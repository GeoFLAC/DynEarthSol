#include <iostream>
#ifdef USE_NPROF
#include <nvToolsExt.h> 
#endif

#include "parameters.hpp"
#include "knn.hpp"

#ifdef ACC

__device__ static double distance2_cuda(const double *a, const double *b) {
    double sum = 0.0;
    for (int i = 0; i < NDIMS; ++i) {
        double diff = a[i] - b[i];
        sum += diff * diff;
    }
    return sum;
}

__device__ static inline int clampi(int v, int a, int b) {
    return v < a ? a : (v > b ? b : v);
}

// CUDA kernel using spatial hash grid for KNN
__global__ static void knn_hashgrid_kernel(
    const double *points, int numPoints,
    const double *queries, int numQueries,
    neighbor *out_results, int k, int nheap, double radius2,
    const HashGrid grid)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numQueries) return;
    const double *query = queries + tid*NDIMS;
    // Find which cell query is in
    int c[NDIMS];
    for (int i = 0; i < NDIMS; ++i) {
        c[i] = int((query[i] - grid.params.origin[i]) / grid.params.cell_size);
        c[i] = clampi(c[i], 0, grid.params.grid_dim[i] - 1);
    }

    neighbor *results = out_results + tid * k;
    for (int i = 0; i < k; ++i) {
        results[i].idx = -1;
        results[i].dist2 = std::numeric_limits<double>::infinity(); // max distance
    }

    // Determine if on boundary
    bool is_on_boundary = (c[0] == 0 || c[1] == 0 ||
                           c[0] == grid.params.grid_dim[0] - 1 ||
                           c[1] == grid.params.grid_dim[1] - 1);

#ifdef THREED
    is_on_boundary = (is_on_boundary || c[2] == 0 || c[2] == grid.params.grid_dim[2] - 1);
#endif

    const int *DC = is_on_boundary ? grid.params.D5 : grid.params.D3;
    int range = is_on_boundary ? grid.params.ND5 : grid.params.ND3;

    int n_cand = 0;
    for (int di = 0; di < range; ++di) {
        int nc[NDIMS];
        for (int i = 0; i < NDIMS; ++i)
            nc[i] = clampi(c[i] + DC[di*NDIMS + i], 0, grid.params.grid_dim[i] - 1);
#ifdef THREED
        int cell_idx = (nc[2] * grid.params.grid_dim[1] + nc[1]) * grid.params.grid_dim[0] + nc[0];
#else
        int cell_idx = nc[1] * grid.params.grid_dim[0] + nc[0];
#endif
        int start = grid.cell_starts[cell_idx];
        int end = grid.cell_starts[cell_idx+1];
        for (int i = start; i < end; ++i) {
            int pi = grid.point_indices[i];
            double d2 = distance2_cuda(points+pi*NDIMS, query);

            if (d2 <= radius2) {
                if (n_cand >= k && results[k - 1].dist2 < d2)
                    continue; // No need to add if we already have k results with smaller distances
                bool is_duplicate = false;
                for (int j = 0; j < k; ++j) {
                    if (results[j].idx == pi) {
                        is_duplicate = true;
                        break;
                    }
                }
                if (is_duplicate) continue;

                for (int ii = 0; ii < k; ++ii) {
                    if (d2 < results[ii].dist2) {
                        // Shift the rest of the result down
                        for (int j = k - 1; j > ii; --j) {
                            results[j] = results[j - 1];
                        }
                        results[ii].idx = pi;
                        results[ii].dist2 = d2;
                        n_cand++;
                        break;
                    }
                }
            }
        }
    }
}

#endif

KNN::KNN(const Param& param, const array_t& points_vec, NANOKDTree& nano_kdtree_,
            double resoTimes_) : 
    resolution(param.mesh.resolution),
    points(points_vec.data()), numPoints(points_vec.size()),
    nano_kdtree(nano_kdtree_),
    resoTimes(resoTimes_)
{
#ifdef ACC
    build_hash_grid(resolution * 2.5);
#else
    nano_kdtree.buildIndex();
#endif
    // use managed memory
    // cudaMallocManaged(&d_grid, sizeof(HashGrid));
    // *d_grid = grid;

    // cudaMallocManaged(&d_points, numPoints * sizeof(double3));
    // cudaMemcpy(d_points, points, numPoints * sizeof(double3), cudaMemcpyHostToDevice);

    // cudaDeviceSynchronize();
}

KNN::~KNN()
{
#ifdef ACC
    delete [] grid.params.D3;
    delete [] grid.params.D5;
#endif
    // use managed memory
    // cudaFree(grid.cell_starts);
    // cudaFree(grid.point_indices);
    // cudaFree(d_grid);

    // cudaFree(d_points);
};

#ifdef ACC

// Host: build hash grid (for simplicity, on host then copy to device)
void KNN::build_hash_grid(double cell_size) {
    // Compute bounds
    double minp[NDIMS], maxp[NDIMS];
    for (int i = 0; i < NDIMS; ++i) {
        minp[i] = points[i];
        maxp[i] = points[i];
    }
    for (int i = 1; i < numPoints; ++i) {
        int idx = i*NDIMS;
        for (int j = 0; j < NDIMS; ++j) {
            minp[j] = std::min(minp[j], points[idx + j]);
            maxp[j] = std::max(maxp[j], points[idx + j]);
        }
    }
    // Small margin
    double eps = 1e-8;
    int dim[3];
    for (int i = 0; i < NDIMS; ++i) {
        minp[i] -= eps;
        maxp[i] += eps;
        dim[i] = int((maxp[i] - minp[i]) / cell_size) + 1;
        grid.params.origin[i] = minp[i];
        grid.params.grid_dim[i] = dim[i];
    }
    int num_cells = dim[0] * dim[1] * dim[2];
    grid.params.cell_size = cell_size;
    grid.num_cells = num_cells;

    // First, count points in each cell
    int_vec cell_counts(num_cells, 0);
    std::vector<unsigned int> cell_codes(numPoints);
    #pragma omp parallel for
    for (int i = 0; i < numPoints; ++i) {
        int idx = i * NDIMS;
        int c[NDIMS];
        for (int j = 0; j < NDIMS; ++j) {
            c[j] = int((points[idx + j] - minp[j]) / cell_size);
            if (c[j] < 0) c[j] = 0;
            if (c[j] >= dim[j]) c[j] = dim[j] - 1;
        }
#ifdef THREED
        unsigned int code = (c[2] * dim[1] + c[1]) * dim[0] + c[0];
#else
        unsigned int code = c[1] * dim[0] + c[0];
#endif
        cell_codes[i] = code;
        #pragma omp atomic update
        cell_counts[code]++;
    }
    // Prefix sum for cell_starts
    int_vec cell_starts(num_cells + 1, 0);
    for (int i = 0; i < num_cells; ++i) {
        cell_starts[i + 1] = cell_starts[i] + cell_counts[i];
    }
    // Fill point_indices (bucket sort)
    int_vec next_indices(num_cells, 0);
    int_vec point_indices(numPoints);
    for (int i = 0; i < num_cells; ++i) next_indices[i] = cell_starts[i];
    for (int i = 0; i < numPoints; ++i) {
        int c = cell_codes[i];
        point_indices[next_indices[c]++] = i;
    }
    // Allocate & copy to device
    cudaMallocManaged(&grid.cell_starts, sizeof(int) * (num_cells + 1));
    cudaMemcpy(grid.cell_starts, cell_starts.data(), sizeof(int) * (num_cells + 1), cudaMemcpyHostToDevice);
    cudaMallocManaged(&grid.point_indices, sizeof(int) * numPoints);
    cudaMemcpy(grid.point_indices, point_indices.data(), sizeof(int) * numPoints, cudaMemcpyHostToDevice);

    grid.params.D3 = new int[grid.params.ND3*NDIMS];
    grid.params.D5 = new int[grid.params.ND5*NDIMS];

    // Neighbor offsets for 3x3x3 and 5x5x5
    int idx = 0;
#ifdef THREED
    for (int dz = -1; dz <= 1; ++dz)
#endif
        for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx) {
                grid.params.D3[idx*NDIMS] = dx;
                grid.params.D3[idx*NDIMS+1] = dy;
#ifdef THREED
                grid.params.D3[idx*NDIMS+2] = dz;
#endif
                idx++;
            }

    // 5x5x5 neighbor offsets: from -2 to +2
    idx = 0;
#ifdef THREED
    for (int dz = -2; dz <= 2; ++dz)
#endif
        for (int dy = -2; dy <= 2; ++dy)
            for (int dx = -2; dx <= 2; ++dx) {
                grid.params.D5[idx*NDIMS] = dx;
                grid.params.D5[idx*NDIMS+1] = dy;
#ifdef THREED
                grid.params.D5[idx*NDIMS+2] = dz;
#endif
                idx++;
            }

}

void KNN::knnSearchCuda_hashgrid(const double *queries, int numQueries,
                      neighbor* results, int k, int nheap, 
                      double radius2, double cell_size) {
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    // use managed memory
    // points and queries already on device or managed
    // double3* d_points;
    // cudaMallocManaged(&d_points, numPoints * sizeof(double3));
    // cudaMemcpy(d_points, points, numPoints * sizeof(double3), cudaMemcpyHostToDevice);

    // double3* d_queries;
    // cudaMallocManaged(&d_queries, numQueries * sizeof(double3));
    // cudaMemcpy(d_queries, queries, numQueries * sizeof(double3), cudaMemcpyHostToDevice);
    // neighbor* d_results;
    // cudaMallocManaged(&results, numQueries * k * sizeof(neighbor));

    int threadsPerBlock = 256;
    int numBlocks = (numQueries + threadsPerBlock - 1) / threadsPerBlock;
    knn_hashgrid_kernel<<<numBlocks, threadsPerBlock>>>(
        points, numPoints, queries, numQueries, results, k, nheap, radius2, grid
    );

    // use managed memory
    // cudaCheckSync("knn_hashgrid_kernel");

    // cudaMemcpy(results, d_results, numQueries * k * sizeof(neighbor), cudaMemcpyDeviceToHost);
    // Cleanup
    // cudaFree(d_points);
    // cudaFree(d_queries);
    // cudaFree(d_results);

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

#endif // ACC

void KNN::search(const array_t& queries, neighbor_vec& neighbors, 
        int k, double resoTimes)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    printf("      Running knn query on %d points ", numPoints);

#ifdef ACC
    printf("(cuda spatial hash grid)\n");

    int heapSize = k * 100;
    double maxDist = resoTimes * resolution;
    double cell_size = maxDist;

    knnSearchCuda_hashgrid(queries.data(), queries.size(),
        neighbors.data(), k, heapSize, maxDist * maxDist, cell_size);

    // long max_size = 1024 * 1024 * 256;
    // int nqueries = queries.size();

    // int nblocks = (double)nqueries * k / (double)max_size;
    // if (nblocks < 1) nblocks = 1;
    // printf("        nqueries: %d, k: %d, npoints: %d, max_size: %d, nblocks: %d\n", nqueries, k, numPoints, max_size, nblocks);

    // int block_size = (nqueries + nblocks - 1) / nblocks;
    // for (int b=0; b<nblocks; b++) {
    //     int start = b * block_size;
    //     int end = std::min(start + block_size, nqueries);
    //     if (start >= end) continue;

    //     size_t free_mem, total_mem;
    //     cudaMemGetInfo(&free_mem, &total_mem);
    //     printf("          Block %3d: %10d to %10d\n", b, start, end);
    //     printf("            GPU memory: free = %zu MB, total = %zu MB\n", free_mem / (1024 * 1024), total_mem / (1024 * 1024));

    //     knnSearchCuda_hashgrid(queries.data() + start, end - start,
    //         neighbors.data() + start*k, k, heapSize, maxDist * maxDist, cell_size);
    // }

    cudaDeviceSynchronize();

#else
    printf("(nano-kdtree)\n");

    #pragma omp parallel for default(none) \
        shared(queries, neighbors, k, nano_kdtree)
    for (int i = 0; i < queries.size(); ++i) {
        neighbor *result = neighbors.data() + i * k;

        size_t_vec nn_idx(k);
        double_vec out_dists_sqr(k);
        KNNResultSet resultSet(k);
        resultSet.init(nn_idx.data(), out_dists_sqr.data());

        nano_kdtree.findNeighbors(resultSet, queries.data() + i * NDIMS);

        for (int j = 0; j < k; ++j) {
            result[j].idx = nn_idx[j];
            result[j].dist2 = out_dists_sqr[j];
        }
    }

#endif

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}
