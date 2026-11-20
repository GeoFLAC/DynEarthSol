#ifndef DYNEARTHSOL3D_KNN_CU
#define DYNEARTHSOL3D_KNN_CU


class KNN
{
public:
    KNN(const array_t& points_vec_, NANOKDTree& nano_kdtree_, int capacity = -1);
    ~KNN();

    // Search for k nearest neighbors.
    neighbor* search(const array_t& queries, const int nquery, const int k_neig,
            bool is_sync_to_host = true, const float* d_guess_radii_sq = nullptr);

    int max_batch_size(int k_neig) const;

    // Print the head of a remesh log line, without a newline:
    //   <stage>knn: <points> pts for <nquery> q w/ <k> k[, <MB of one batch on the GPU>]
    // The search structure is fixed at build time; the [runtime][device] line names it.
    void print_search(const char* stage, long nquery, int k_neig, int batch_queries) const;
private:
    const array_t& points_vec;
    int numPoints;
    size_t h_results_capacity;
    neighbor* h_results;

    NANOKDTree& nano_kdtree;

#ifdef ACC
    void* m_bvh_state;  // KNNBVHState*, managed by knn_bvh.cu
#endif
};


#endif // DYNEARTHSOL3D_KNN_CU
