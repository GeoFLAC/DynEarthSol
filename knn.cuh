#ifndef DYNEARTHSOL3D_KNN_CU
#define DYNEARTHSOL3D_KNN_CU

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

class CudaKNN
{
public:
    CudaKNN(const Param& param, const double3_vec& points_vec_, 
            double resoTimes = 3);
    ~CudaKNN();

    void search_grid(const double3_vec& queries, neighbor_vec& neighbors, 
            int k, double resoTimes = 3);
private:
    const double3* points;
    int numPoints;

    const int resolution;
    const int resoTimes;
    HashGrid grid;
    HashGrid* d_grid;

    double3* d_points;

    void build_hash_grid(double cell_size);

    void knnSearchCuda_hashgrid(const double3* queries, int numQueries,
                        neighbor* results, int k, int nheap, 
                        double radius2, double cell_size);

    inline void cudaCheckSync(const char* msg) {
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            fprintf(stderr, "CUDA Sync Error (%s): %s\n", msg, cudaGetErrorString(err));
            exit(EXIT_FAILURE);
        }
    };

};

#endif

#endif // DYNEARTHSOL3D_KNN_CU