#include <iostream>
#ifdef USE_NPROF
#include <nvToolsExt.h> 
#endif

#include "parameters.hpp"
#include "knn.cuh"

#ifdef ACC

__device__ static double distance2_cuda(const double3 &a, const double3 &b) {
    double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx*dx + dy*dy + dz*dz;
}

// Morton code or spatial hash for 3D grid cell coordinate
__host__ __device__ static inline unsigned int morton3D(int x, int y, int z) {
    // Interleave bits of x, y, z. Works for up to 10 bits per coordinate.
    unsigned int answer = 0;
    for (unsigned int i = 0; i < 10; ++i) {
        answer |= ((x >> i) & 1) << (3 * i + 0);
        answer |= ((y >> i) & 1) << (3 * i + 1);
        answer |= ((z >> i) & 1) << (3 * i + 2);
    }
    return answer;
}

__device__ static inline int clampi(int v, int a, int b) {
    return v < a ? a : (v > b ? b : v);
}

// Device: locate cell index for a point
__device__ static int get_cell_index(const double3& p, const HashGridParams& params) {
    int cx = int((p.x - params.origin.x) / params.cell_size);
    int cy = int((p.y - params.origin.y) / params.cell_size);
    int cz = int((p.z - params.origin.z) / params.cell_size);
    cx = clampi(cx, 0, params.grid_dim.x - 1);
    cy = clampi(cy, 0, params.grid_dim.y - 1);
    cz = clampi(cz, 0, params.grid_dim.z - 1);
    return (cz * params.grid_dim.y + cy) * params.grid_dim.x + cx;
}

// CUDA kernel using spatial hash grid for KNN
__global__ static void knn_hashgrid_kernel(
    const double3 *points, int numPoints,
    const double3 *queries, int numQueries,
    neighbor *out_results, int k, int nheap, double radius2,
    const HashGrid grid)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numQueries) return;
    double3 query = queries[tid];
    // Find which cell query is in
    int cx = int((query.x - grid.params.origin.x) / grid.params.cell_size);
    int cy = int((query.y - grid.params.origin.y) / grid.params.cell_size);
    int cz = int((query.z - grid.params.origin.z) / grid.params.cell_size);
    cx = clampi(cx, 0, grid.params.grid_dim.x - 1);
    cy = clampi(cy, 0, grid.params.grid_dim.y - 1);
    cz = clampi(cz, 0, grid.params.grid_dim.z - 1);

    neighbor *results = out_results + tid * k;
    for (int i = 0; i < k; ++i) {
        results[i].idx = -1;
        results[i].dist2 = std::numeric_limits<double>::infinity(); // max distance
    }

    // Neighbor offsets for 3x3x3 and 5x5x5
    const int DX3[27] = {
        0, -1, 0, 1, 0, 0, 0,
        1, 1, 1, 0, 0, 0,
        -1, -1, -1, 1, 1, 1,
        -1, -1, -1, 0, 0, 0, 0, 0
    };
    const int DY3[27] = {
        0, 0, -1, 0, 1, 0, 0,
        1, 0, -1, 1, 0, -1,
        1, 0, -1, 0, 1, -1,
        0, 1, -1, 1, 0, -1, 1, 0
    };
    const int DZ3[27] = {
        0, 0, 0, 0, 0, 1, -1,
        0, 1, -1, 1, -1, 0,
        0, 1, -1, 1, 0, -1,
        1, 0, -1, 1, 1, 1, -1, -1
    };
    // 5x5x5 neighbor offsets: from -2 to +2
    __shared__ int DX5[125], DY5[125], DZ5[125];
    // Only need to initialize once per block
    if (threadIdx.x == 0) {
        int idx = 0;
        for (int dz = -2; dz <= 2; ++dz)
            for (int dy = -2; dy <= 2; ++dy)
                for (int dx = -2; dx <= 2; ++dx) {
                    DX5[idx] = dx;
                    DY5[idx] = dy;
                    DZ5[idx] = dz;
                    idx++;
                }
    }
    __syncthreads();
    // Determine if on boundary
    bool is_on_boundary = (cx == 0 || cy == 0 || cz == 0 ||
                           cx == grid.params.grid_dim.x - 1 ||
                           cy == grid.params.grid_dim.y - 1 ||
                           cz == grid.params.grid_dim.z - 1);
    const int* DX = is_on_boundary ? DX5 : DX3;
    const int* DY = is_on_boundary ? DY5 : DY3;
    const int* DZ = is_on_boundary ? DZ5 : DZ3;
    int range = is_on_boundary ? 125 : 27;

    int n_cand = 0;
    for (int di = 0; di < range; ++di) {
        int ncx = clampi(cx + DX[di], 0, grid.params.grid_dim.x-1);
        int ncy = clampi(cy + DY[di], 0, grid.params.grid_dim.y-1);
        int ncz = clampi(cz + DZ[di], 0, grid.params.grid_dim.z-1);
        int cell_idx = (ncz * grid.params.grid_dim.y + ncy) * grid.params.grid_dim.x + ncx;
        int start = grid.cell_starts[cell_idx];
        int end = grid.cell_starts[cell_idx+1];
        for (int i = start; i < end; ++i) {
            int pi = grid.point_indices[i];
            double d2 = distance2_cuda(points[pi], query);

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

CudaKNN::CudaKNN(const Param& param, const double3_vec& points_vec,
            double resoTimes_) : 
    resolution(param.mesh.resolution),
    points(points_vec.data()), numPoints(points_vec.size()),
    resoTimes(resoTimes_)
{
    build_hash_grid(resolution * resoTimes);

    cudaMallocManaged(&d_grid, sizeof(HashGrid));
    *d_grid = grid;

    cudaMallocManaged(&d_points, numPoints * sizeof(double3));
    // cudaMemcpy(d_points, points, numPoints * sizeof(double3), cudaMemcpyHostToDevice);

    // cudaDeviceSynchronize();
}

CudaKNN::~CudaKNN()
{
    cudaFree(grid.cell_starts);
    cudaFree(grid.point_indices);
    cudaFree(d_grid);

    cudaFree(d_points);
};

// Host: build hash grid (for simplicity, on host then copy to device)
void CudaKNN::build_hash_grid(double cell_size) {
    // Compute bounds
    double3 minp = points[0], maxp = points[0];
    for (int i = 1; i < numPoints; ++i) {
        minp.x = std::min(minp.x, points[i].x);
        minp.y = std::min(minp.y, points[i].y);
        minp.z = std::min(minp.z, points[i].z);
        maxp.x = std::max(maxp.x, points[i].x);
        maxp.y = std::max(maxp.y, points[i].y);
        maxp.z = std::max(maxp.z, points[i].z);
    }
    // Small margin
    double eps = 1e-8;
    minp.x -= eps; minp.y -= eps; minp.z -= eps;
    maxp.x += eps; maxp.y += eps; maxp.z += eps;
    int3 dim;
    dim.x = int((maxp.x - minp.x) / cell_size) + 1;
    dim.y = int((maxp.y - minp.y) / cell_size) + 1;
    dim.z = int((maxp.z - minp.z) / cell_size) + 1;
    int num_cells = dim.x * dim.y * dim.z;
    grid.params.origin = minp;
    grid.params.cell_size = cell_size;
    grid.params.grid_dim = dim;
    grid.num_cells = num_cells;

    // First, count points in each cell
    std::vector<int> cell_counts(num_cells, 0);
    std::vector<unsigned int> cell_codes(numPoints);
    for (int i = 0; i < numPoints; ++i) {
        int cx = int((points[i].x - minp.x) / cell_size);
        int cy = int((points[i].y - minp.y) / cell_size);
        int cz = int((points[i].z - minp.z) / cell_size);
        if (cx < 0) cx = 0; if (cy < 0) cy = 0; if (cz < 0) cz = 0;
        if (cx >= dim.x) cx = dim.x - 1;
        if (cy >= dim.y) cy = dim.y - 1;
        if (cz >= dim.z) cz = dim.z - 1;
        unsigned int code = (cz * dim.y + cy) * dim.x + cx;
        cell_codes[i] = code;
        cell_counts[code]++;
    }
    // Prefix sum for cell_starts
    std::vector<int> cell_starts(num_cells + 1, 0);
    for (int i = 0; i < num_cells; ++i) {
        cell_starts[i + 1] = cell_starts[i] + cell_counts[i];
    }
    // Fill point_indices (bucket sort)
    std::vector<int> next_indices(num_cells, 0);
    std::vector<int> point_indices(numPoints);
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
}

void CudaKNN::knnSearchCuda_hashgrid(const double3* queries, int numQueries,
                      neighbor* results, int k, int nheap, 
                      double radius2, double cell_size) {
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    // points and queries already on device or managed
    // double3* d_points;
    // cudaMallocManaged(&d_points, numPoints * sizeof(double3));
    cudaMemcpy(d_points, points, numPoints * sizeof(double3), cudaMemcpyHostToDevice);

    double3* d_queries;
    cudaMallocManaged(&d_queries, numQueries * sizeof(double3));
    cudaMemcpy(d_queries, queries, numQueries * sizeof(double3), cudaMemcpyHostToDevice);
    neighbor* d_results;
    cudaMallocManaged(&d_results, numQueries * k * sizeof(neighbor));

    int threadsPerBlock = 256;
    int numBlocks = (numQueries + threadsPerBlock - 1) / threadsPerBlock;
    knn_hashgrid_kernel<<<numBlocks, threadsPerBlock>>>(
        d_points, numPoints, d_queries, numQueries, d_results, k, nheap, radius2, grid
    );

    cudaCheckSync("knn_hashgrid_kernel");

    cudaMemcpy(results, d_results, numQueries * k * sizeof(neighbor), cudaMemcpyDeviceToHost);
    // Cleanup
    // cudaFree(d_points);
    cudaFree(d_queries);
    cudaFree(d_results);

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void CudaKNN::search_grid(const double3_vec& queries, neighbor_vec& neighbors, 
        int k, double resoTimes)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    std::cout << "      Running knn query on " << numPoints
            << " points (spatial hash grid)" << std::endl;

    int heapSize = k * 100;
    double maxDist = resoTimes * resolution;
    double cell_size = maxDist;

    long max_size = 1024 * 1024 * 256;
    int nqueries = queries.size();

    int nblocks = (double)nqueries * k / (double)max_size;
    if (nblocks < 1) nblocks = 1;
    printf("        nqueries: %d, k: %d, npoints: %d, max_size: %d, nblocks: %d\n", nqueries, k, numPoints, max_size, nblocks);

    int block_size = (nqueries + nblocks - 1) / nblocks;
    for (int b=0; b<nblocks; b++) {
        int start = b * block_size;
        int end = std::min(start + block_size, nqueries);
        if (start >= end) continue;

        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);
        printf("          Block %3d: %10d to %10d\n", b, start, end);
        printf("            GPU memory: free = %zu MB, total = %zu MB\n", free_mem / (1024 * 1024), total_mem / (1024 * 1024));

        knnSearchCuda_hashgrid(queries.data() + start, end - start,
            neighbors.data() + start*k, k, heapSize, maxDist * maxDist, cell_size);
    }

    cudaDeviceSynchronize();
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}
#endif // ACC