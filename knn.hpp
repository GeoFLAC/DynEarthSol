#ifndef DYNEARTHSOL3D_KNN_CU
#define DYNEARTHSOL3D_KNN_CU

struct HashGridParams {
    double origin[NDIMS];
    double cell_size;
    int grid_dim[NDIMS];
    int *D3;
    int *D5;
#ifdef THREED
    int ND5 = 5*5*5;
    int ND3 = 3*3*3;
#else
    int ND5 = 5*5;
    int ND3 = 3*3; 
#endif
};

struct HashGrid {
    // Flat arrays for CUDA-friendly access
    int* cell_starts; // size = num_cells+1
    int* point_indices; // size = npoints
    HashGridParams params;
    int num_cells;
};

class KNN
{
public:
    KNN(const Param& param, const array_t& points_vec_, NANOKDTree& nano_kdtree_,
            double resoTimes = 3);
    ~KNN();

    void search(const array_t& queries, neighbor_vec& neighbors, 
            int k, double resoTimes = 3);
private:
    const double* points;
    int numPoints;

    const int resolution;
    const int resoTimes;
    HashGrid grid;
    NANOKDTree& nano_kdtree;
    // HashGrid* d_grid;

    // double3* d_points;
#ifdef ACC
    void build_hash_grid(double cell_size);

    void knnSearchCuda_hashgrid(const double* queries, int numQueries,
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
#endif

};


#endif // DYNEARTHSOL3D_KNN_CU