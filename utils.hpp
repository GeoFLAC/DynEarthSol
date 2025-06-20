#ifndef DYNEARTHSOL3D_UTILS_HPP
#define DYNEARTHSOL3D_UTILS_HPP
#ifdef USE_NPROF
#include <nvToolsExt.h> 
#endif

#include "parameters.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <cfloat>
#include <math.h>
#include <iomanip>
#if defined(_WIN32)
#include <windows.h>
#else
#include <time.h>
#endif

static void print(std::ostream& os, const double& x)
{
  os << x;
}


static void print(std::ostream& os, const int& x)
{
  os << x;
}


static void print(std::ostream& os, const std::size_t& x)
{
  os << x;
}


template <typename T1, typename T2>
void print(std::ostream& os, const std::pair<T1,T2>& x)
{
    os << x.first << ':' << x.second;
}


template <typename T>
void print(std::ostream& os, const T& A, std::size_t size)
{
  os << "[";
  for (std::size_t i = 0; i != size; ++i) {
    print(os, A[i]);
    if (i+1 != size)
      os << ", ";
  }
  os << "]";
}


template <typename Array>
void print(std::ostream& os, const Array& A)
{
  typename Array::const_iterator i;
  os << "[";
  for (i = A.begin(); i != A.end(); ++i) {
    print(os, *i);
    os << ", ";
  }
  os << "]";
}

#ifdef USE_NPROF

#pragma acc routine seq
static inline double log_lookup(const double_vec& log_table, double x) {
    if (x <= 0.0)
        printf("Error: log_lookup called with x <= 0 (%g)\n", x);

    int exponent = 0;
    while (x < LOG_XMIN) {
        x *= 10.0;
        exponent -= 1;
    }
    while (x > LOG_XMAX) {
        x *= 0.1;
        exponent += 1;
    }

    int idx = static_cast<int>((x - LOG_XMIN) / LOG_XDELTA);
    double dx = x - (LOG_XMIN + idx * LOG_XDELTA);
    double slope = (log_table[idx + 1] - log_table[idx]) / LOG_XDELTA;
    return log_table[idx] + slope * dx + exponent * 2.302585092994046; // LN_10
}

#pragma acc routine seq
static inline double tan_lookup(const double_vec& tan_table, const double x) {
    if (x < TAN_XMIN || x > TAN_XMAX) {
        printf("Error: tan_lookup called with x out of range (%g)\n", x);
    }

    int idx = static_cast<int>((x - TAN_XMIN) / TAN_XDELTA);
    double dx = x - (TAN_XMIN + idx * TAN_XDELTA);
    double slope = (tan_table[idx + 1] - tan_table[idx]) / TAN_XDELTA;
    return tan_table[idx] + slope * dx;
}

#pragma acc routine seq
static inline double sin_lookup(const double_vec& sin_table, const double x) {
    if (x < SIN_XMIN || x > SIN_XMAX) {
        printf("Error: sin_lookup called with x out of range (%g)\n", x);
    }

    int idx = static_cast<int>((x - SIN_XMIN) / SIN_XDELTA);
    double dx = x - (SIN_XMIN + idx * SIN_XDELTA);
    double slope = (sin_table[idx + 1] - sin_table[idx]) / SIN_XDELTA;
    return sin_table[idx] + slope * dx;
}

#endif

#pragma acc routine seq
static inline double tan_safe(const double_vec& tan_table, const double x) {
#ifdef USE_NPROF
    return tan_lookup(tan_table, x);
#else
    return std::tan(x);
#endif
}

#pragma acc routine seq
static inline double log_safe(const double_vec& log_table, const double x) {
#ifdef USE_NPROF
    return log_lookup(log_table, x);
#else
    return std::log(x);
#endif
}

#pragma acc routine seq
static inline double sin_safe(const double_vec& sin_table, const double x) {
#ifdef USE_NPROF
    return sin_lookup(sin_table, x);
#else
    return std::sin(x);
#endif
}

#pragma acc routine seq
static double pow_safe(const double_vec& log_table, const double x, const double y) {
#ifdef USE_NPROF
    if (x < 0.) {
        printf("Error: pow_safe called with x < 0 (%g)\n", x);
    } else if (x == 0.) {
        return 0.0; // Avoid log(0) which is undefined
    } else
        return exp(y * log_lookup(log_table, x));
#else
    return std::pow(x, y);
#endif
}

#pragma acc routine seq
static double pow_1_5(const double x) {
    return x * sqrt(x);
}

#pragma acc routine seq
static double pow_2(const double x) {
    return x * x;
}

#pragma acc routine seq
static double trace(const double* s)
{
#ifdef THREED
    return s[0] + s[1] + s[2];
#else
    return s[0] + s[1];
#endif
}


static double second_invariant2(const double* t)
{
#ifdef THREED
    double a = (t[0] + t[1] + t[2]) / 3;
    return ( 0.5 * ((t[0]-a)*(t[0]-a) + (t[1]-a)*(t[1]-a) + (t[2]-a)*(t[2]-a))
             + t[3]*t[3] + t[4]*t[4] + t[5]*t[5] );
#else
    return 0.25*(t[0]-t[1])*(t[0]-t[1]) + t[2]*t[2];
#endif
}


static double second_invariant(const double* t)
{
    /* second invariant of the deviatoric part of tensor t
     * defined as: td = deviatoric(t); sqrt( td(i,j) * td(i,j) / 2)
     */
    return std::sqrt(second_invariant2(t));
}


static int findNearestNeighbourIndex( double x_new, const double_vec& x )
{
    /* find nearest neighbour index for interpolation
     * x vector only can be ascending
     */
    double dist = DBL_MAX;
    int idx = -1;
    for (size_t i = 0; i < x.size(); ++i ) {
        double newDist = x_new - x[i];
        if ( newDist >= 0 && newDist <= dist ) {
            dist = newDist;
            idx = i;
        }
    }

    return idx;
}


static double interp1(const double_vec& x, const double_vec& y, double x_new)
{
    int idx = findNearestNeighbourIndex( x_new, x);
    double slope = 0;

    if (idx < 0)
        idx = 0;
    else if ( idx < static_cast<int>(x.size()-1) )
        slope = (y[idx+1] - y[idx]) / (x[idx+1] - x[idx]);

    return slope * (x_new-x[idx]) + y[idx];
}

static int64_t get_nanoseconds() {
    #pragma acc wait

    #if defined(_WIN32)
    LARGE_INTEGER frequency, counter;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&counter);
    return (int64_t)((double)counter.QuadPart / frequency.QuadPart * 1e9);
    #else
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (int64_t)ts.tv_sec * 1e9 + ts.tv_nsec;
    #endif
}

static void print_time_ns(const int64_t duration) {
    int hours = duration / (int64_t)3600000000000;
    int minutes = (duration % (int64_t)3600000000000) / (int64_t)60000000000;
    double seconds = (duration % (int64_t)60000000000) / 1e9;
    std::cout << std::setw(3) << std::setfill('0') << hours << ":"
    << std::setw(2) << std::setfill('0') << minutes << ":"
    << std::setw(9) << std::fixed << std::setprecision(6) << std::setfill('0') << seconds;
}


#pragma acc routine seq
static int out_nan_error(const char* msg, const int idx0, const int idx1 = -1) {
    if (idx1 >= 0)
        printf("Error: %s[%d][%d] becomes NaN\n", msg, idx0, idx1);
    else
        printf("Error: %s[%d] becomes NaN\n", msg, idx0);
    return 1;
}

static void check_nan(const Variables& var) {
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    #pragma acc serial
    int is_nan = 0;

#ifndef ACC
    #pragma omp parallel default(none) shared(var,is_nan)
#endif
    {
#ifndef ACC
        #pragma omp for reduction(+:is_nan)
#endif
        #pragma acc parallel loop reduction(+:is_nan)
        for (int e=0; e<var.nelem;e++) {
            if (std::isnan((*var.volume)[e]))
                is_nan += out_nan_error("volume", e);
            
            if (std::isnan((*var.dpressure)[e]))
                is_nan += out_nan_error("dpressure", e);

            if (std::isnan((*var.viscosity)[e]))
               is_nan +=  out_nan_error("viscosity", e);
            
            for (int i=0; i<NODES_PER_ELEM;i++)
                if(std::isnan((*var.connectivity)[e][i]))
                    is_nan += out_nan_error("connectivity", e, i);

            for (int i=0; i<NSTR; i++)
                if (std::isnan((*var.stress)[e][i]))
                    is_nan += out_nan_error("stress", e, i);

            for (int i=0; i<NODES_PER_ELEM; i++)
                if(std::isnan((*var.shpdx)[e][i]))
                    is_nan += out_nan_error("shpdx", e, i);
            for (int i=0; i<NODES_PER_ELEM; i++)
                if(std::isnan((*var.shpdy)[e][i]))
                    is_nan += out_nan_error("shpdy", e, i);
            for (int i=0; i<NODES_PER_ELEM; i++)
                if(std::isnan((*var.shpdz)[e][i]))
                    is_nan += out_nan_error("shpdz", e, i);
        }

#ifndef ACC
        #pragma omp for reduction(+:is_nan)
#endif
        #pragma acc parallel loop reduction(+:is_nan)
        for (int n=0; n<var.nnode; n++) {
            if (std::isnan((*var.temperature)[n]))
                is_nan += out_nan_error("temperature", n);

            if (std::isnan((*var.tmass)[n]))
                is_nan += out_nan_error("tmass", n);

            for (int i=0; i<NDIMS; i++) {
                if (std::isnan((*var.force)[n][i]))
                    is_nan += out_nan_error("force", n, i);

                if (std::isnan((*var.vel)[n][i]))
                    is_nan += out_nan_error("vel", n, i);

                if (std::isnan((*var.coord)[n][i]))
                    is_nan += out_nan_error("coord", n, i);                
            }
        }
    }

    if (is_nan > 0) {
        std::cerr << "Error: " << is_nan << " NaN values found in the variables." << std::endl;
        std::exit(1);
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

#ifdef ACC

struct HashGridParams {
    double3 origin;
    double cell_size;
    int3 grid_dim;
};

struct HashGrid {
    // Flat arrays for CUDA-friendly access
    int* cell_starts; // size = num_cells+1
    int* point_indices; // size = npoints
    HashGridParams params;
    int num_cells;
};

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

static inline void cudaCheckSync(const char* msg) {
    cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Sync Error (%s): %s\n", msg, cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

// Host: build hash grid (for simplicity, on host then copy to device)
static void build_hash_grid(const double3* points, int npoints, double cell_size, HashGrid& grid) {
    // Compute bounds
    double3 minp = points[0], maxp = points[0];
    for (int i = 1; i < npoints; ++i) {
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
    std::vector<unsigned int> cell_codes(npoints);
    for (int i = 0; i < npoints; ++i) {
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
    std::vector<int> point_indices(npoints);
    for (int i = 0; i < num_cells; ++i) next_indices[i] = cell_starts[i];
    for (int i = 0; i < npoints; ++i) {
        int c = cell_codes[i];
        point_indices[next_indices[c]++] = i;
    }
    // Allocate & copy to device
    cudaMallocManaged(&grid.cell_starts, sizeof(int) * (num_cells + 1));
    cudaMemcpy(grid.cell_starts, cell_starts.data(), sizeof(int) * (num_cells + 1), cudaMemcpyHostToDevice);
    cudaMallocManaged(&grid.point_indices, sizeof(int) * npoints);
    cudaMemcpy(grid.point_indices, point_indices.data(), sizeof(int) * npoints, cudaMemcpyHostToDevice);
}

static void knnSearchCuda_hashgrid(const double3* points, int numPoints,
                      const double3* queries, int numQueries,
                      neighbor* results, int k, int nheap, 
                      double radius2, double cell_size) {
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // Build hash grid
    HashGrid grid;
    build_hash_grid(points, numPoints, cell_size, grid);
    // Copy grid struct to device (params is trivially copyable)
    HashGrid* d_grid;
    cudaMallocManaged(&d_grid, sizeof(HashGrid));
    *d_grid = grid;

    // points and queries already on device or managed
    double3* d_points;
    cudaMallocManaged(&d_points, numPoints * sizeof(double3));
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
    cudaFree(d_points);
    cudaFree(d_queries);
    cudaFree(d_results);
    cudaFree(grid.cell_starts);
    cudaFree(grid.point_indices);
    cudaFree(d_grid);
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

class CudaKNN
{
public:
    CudaKNN(const Param& param, const double3_vec& points) : 
        points(points), resolution(param.mesh.resolution) {};
    ~CudaKNN() {};

    void search_grid(const double3_vec& queries, neighbor_vec& neighbors, 
            int k, double resoTimes = 3)
    {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif
        std::cout << "      Running knn query on " << points.size()
                << " points (spatial hash grid)" << std::endl;

        int heapSize = k * 100;
        double maxDist = resoTimes * resolution;
        double cell_size = maxDist;

        long max_size = 1024 * 1024 * 512;
        int nqueries = queries.size();
        int npoints = points.size();

        int nblocks = (double)nqueries * k / (double)max_size + 1;
        if (nblocks < 1) nblocks = 1;
        printf("        nqueries: %d, k: %d, npoints: %d, max_size: %d, nblocks: %d\n", nqueries, k, npoints, max_size, nblocks);

        int block_size = (nqueries + nblocks - 1) / nblocks;
        for (int b=0; b<nblocks; b++) {
            int start = b * block_size;
            int end = std::min(start + block_size, nqueries);
            if (start >= end) continue;

            printf("          Block %3d: %10d to %10d\n", b, start, end);

            knnSearchCuda_hashgrid(points.data(), npoints, queries.data() + start, end - start,
                neighbors.data() + start*k, k, heapSize, maxDist * maxDist, cell_size);
        }

        cudaDeviceSynchronize();
#ifdef USE_NPROF
        nvtxRangePop();
#endif
    };
private:
    const double3_vec& points;
    const int resolution;
};
#endif // ACC

#endif // DYNEARTHSOL3D_UTILS_HPP