#ifndef DYNEARTHSOL3D_KNN_CU
#define DYNEARTHSOL3D_KNN_CU


class KNN
{
public:
    KNN(const Param& param, const array_t& points_vec_, NANOKDTree& nano_kdtree_, bool is_msg_ = true,
            int capacity = -1);
    ~KNN();

    // Search for k nearest neighbors.
    neighbor* search(const array_t& queries, const int nquery, const int k_neig,
            bool is_sync_to_host = true, const float* d_guess_radii_sq = nullptr);

    int max_batch_size(int k_neig) const;

    const char* backend_name() const {
#ifdef ACC
        return "lbvh GPU";
#else
        return "nano-kdtree";
#endif
    }

    // Fill buf with backend name.  GPU: "lbvh GPU"  CPU: "nano-kdtree"
    void backend_str(char* buf, size_t n) const;

    // GPU memory consumed by d_queries + d_results for one batch.
    // Returns 0 for CPU builds (no GPU allocation).
    static size_t query_mem_bytes(int nquery, int k) {
#ifdef ACC
        return (size_t)nquery * (3 * sizeof(float) + (size_t)k * sizeof(neighbor));
#else
        (void)nquery; (void)k; return 0;
#endif
    }
    inline int get_npoints() const { return numPoints; }
private:
    const array_t& points_vec;
    int numPoints;
    size_t h_results_capacity;
    neighbor* h_results;

    const int resolution;
    const bool is_msg;
    NANOKDTree& nano_kdtree;

#ifdef ACC
    void* m_bvh_state;  // KNNBVHState*, managed by knn_bvh.cu
#endif
};


#endif // DYNEARTHSOL3D_KNN_CU
