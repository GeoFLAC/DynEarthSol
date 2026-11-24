#ifndef DES_MPI_ACC_HPP
#define DES_MPI_ACC_HPP

// Only compact timestep commands offload in a GPU MPI build. Initialization,
// MMG, remapping and global output remain CPU mesh-management operations.
// des_acc_offload_enabled() is the if() clause of every OpenACC loop that can
// run on root's management mesh: false there in a GPU MPI build, else true.
// AccOffloadScope sets it for the duration of one MPI operation.
inline bool& des_acc_offload_flag() {
    static bool enabled=false;
    return enabled;
}
#if defined(USEMPI) && defined(ACC)
inline bool des_acc_offload_enabled() {
    return des_acc_offload_flag();
}
#else
// Every other build always offloads, decided at compile time.
constexpr bool des_acc_offload_enabled() {
    return true;
}
#endif

class AccOffloadScope {
    bool previous_;
public:
    explicit AccOffloadScope(bool enabled):previous_(des_acc_offload_flag()) {
        des_acc_offload_flag()=enabled;
    }
    ~AccOffloadScope() {des_acc_offload_flag()=previous_;}
    AccOffloadScope(const AccOffloadScope&)=delete;
    AccOffloadScope& operator=(const AccOffloadScope&)=delete;
};
#endif
