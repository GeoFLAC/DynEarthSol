#ifndef DES_MPI_ACC_HPP
#define DES_MPI_ACC_HPP

// Only compact timestep commands offload in a GPU MPI build. Initialization,
// MMG, remapping and global output remain CPU mesh-management operations.
inline bool& mpi_gpu_offload_flag() {
    static bool enabled=false;
    return enabled;
}
inline bool des_acc_offload() {
#if defined(USEMPI) && defined(ACC)
    return mpi_gpu_offload_flag();
#else
    return true;
#endif
}

class MpiAccScope {
    bool previous_;
public:
    explicit MpiAccScope(bool enabled):previous_(mpi_gpu_offload_flag()) {
        mpi_gpu_offload_flag()=enabled;
    }
    ~MpiAccScope() {mpi_gpu_offload_flag()=previous_;}
    MpiAccScope(const MpiAccScope&)=delete;
    MpiAccScope& operator=(const MpiAccScope&)=delete;
};
#endif
