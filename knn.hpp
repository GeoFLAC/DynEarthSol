#ifndef DYNEARTHSOL3D_KNN_CU
#define DYNEARTHSOL3D_KNN_CU


class KNN
{
public:
    KNN(const Param& param, const array_t& points_vec_, NANOKDTree& nano_kdtree_,
            int capacity = -1);
    ~KNN();

    // Search for k nearest neighbors.
    neighbor* search(const array_t& queries, const int nquery, const int k_neig,
            bool is_sync_to_host = true, const float* d_guess_radii_sq = nullptr);

    int max_batch_size(int k_neig) const;
private:
    const array_t& points_vec;
    int numPoints;
    size_t h_results_capacity;
    neighbor* h_results;

    const int resolution;
    NANOKDTree& nano_kdtree;

#ifdef ACC
    void* m_bvh_state;  // KNNBVHState*, managed by knn_bvh.cu
#endif
};


#endif // DYNEARTHSOL3D_KNN_CU
