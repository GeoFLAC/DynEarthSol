#ifndef DES_MPI_SIMULATION_HPP
#define DES_MPI_SIMULATION_HPP
#include "parameters.hpp"
struct PTResult;
struct MpiNodeRows {const int* ids; int count;};
#ifndef USEMPI
enum class MpiEntity { Node, Element };
#endif
#ifdef USEMPI
#include <memory>
#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

enum class MpiEntity { Node, Element };
// Call before parsing parameters: managed allocations must use this rank's GPU.
void mpi_select_device();
// World rank set by the MPI launcher, readable before MPI_Init; 0 if unknown.
int mpi_launcher_rank();
enum class MpiOp {
    State, Temperature, PorePressure, StrainRate, Stress, Force, Velocity,
    Coordinate, Rotation, Volume, Dvoldt, Edvoldt, Nmd, Mass, TimeStep,
    OldMeanStress, VelocityBC, Residual, PhysicalPT, InitialPT, Surface, SurfaceIdle, SurfaceTopography, SwapVolume,
    Remesh, PhaseChange, Collect, Geometry, Boundary, Sample, Average, Metrics, Report, Stop, Count
};

// Root keeps the global management/output mesh. Each compute rank, including
// root, advances compact owned/ghost arrays and exchanges only neighbor halos.
class MpiSimulation {
public:
    explicit MpiSimulation(Param&);
    ~MpiSimulation();
    bool is_root() const;
    // Root has attached its management mesh to the compact backend.
    bool distributed() const;
    bool in_kernel() const;
    int worker_loop();
    void attach(Variables&);
    void synchronize(Variables&);
    bool dispatch(MpiOp, const Variables&, double argument = 0);
    double result() const;
    PTResult pt_result() const;
    int owned_nodes() const;
    int owned_elements() const;
    // The whole mesh's number of one of this rank's elements.
    int global_element(int element) const;
    int global_nodes() const;
    void exchange(MpiEntity, std::initializer_list<std::pair<double*,int>>);
    void begin_exchange(MpiEntity, std::initializer_list<std::pair<double*,int>>);
    void finish_exchange();
    MpiNodeRows gather_rows(int pass) const;
    void collect(const Variables&, bool averages = false);
    void geometry(const Variables&);
    void boundary(const Variables&);
    void sample(const Variables&, const int_vec& nodes, const int_vec& elements);
    bool metrics(const Variables&, double& speed, double& plastic_rate);
    void minimum(double* data, int count);
    double maximum(double value);
    double sum(double value);
    void sum(double* data, int count);
    void maximum(double* data, int count);
    int sum(int value);
    void maximum(int* data, int count);
    void report();
    // The manifest's [runtime.mpi] fields, gathered from every rank at start.
    const std::vector<std::pair<std::string, std::string>>& runtime_fields() const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    MpiSimulation(const MpiSimulation&) = delete;
    MpiSimulation& operator=(const MpiSimulation&) = delete;
};
#endif

// Whether this process prints a failure or warning that every rank reaches together
// inside an operation: rank zero alone, so it appears once, as in a serial run.
inline bool mpi_reports(const Variables& var) {
#ifdef USEMPI
    if (var.mpi && var.mpi->in_kernel()) return var.mpi->is_root();
#endif
    return true;
}
// A failure flag seen on some ranks, made common to all inside an operation, so that
// every rank takes the same branch (see die_together()).
inline int mpi_any(const Variables& var, int flag) {
#ifdef USEMPI
    if (var.mpi && var.mpi->in_kernel()) var.mpi->maximum(&flag, 1);
#endif
    return flag;
}

inline int mpi_gather_passes(const Variables& v) {
#ifdef USEMPI
    if(v.mpi && v.mpi->in_kernel()) return 2;
#endif
    return 1;
}
inline MpiNodeRows mpi_gather_nodes(const Variables& v,int pass) {
#ifdef USEMPI
    if(v.mpi && v.mpi->in_kernel()) return v.mpi->gather_rows(pass);
#endif
    return {nullptr,v.nnode};
}

inline int mpi_element_end(const Variables& var) {
#ifdef USEMPI
    if (var.mpi && var.mpi->in_kernel()) return var.mpi->owned_elements();
#endif
    return var.nelem;
}
inline int mpi_node_end(const Variables& var) {
#ifdef USEMPI
    if (var.mpi && var.mpi->in_kernel()) return var.mpi->owned_nodes();
#endif
    return var.nnode;
}
inline void mpi_sync(const Variables& var, double_vec& values, MpiEntity entity) {
#ifdef USEMPI
    if (var.mpi && var.mpi->in_kernel()) var.mpi->exchange(entity,{{values.data(),1}});
#endif
}
template<int N> inline void mpi_sync(const Variables& var, Array2D<double,N>& values,
                                     MpiEntity entity, int columns = N) {
#ifdef USEMPI
    if (var.mpi && var.mpi->in_kernel()) var.mpi->exchange(entity,{{values.data(),columns}});
#endif
}
template<int N> inline void mpi_begin_gather(const Variables& var,Array2D<double,N>& values,int columns) {
#ifdef USEMPI
    if(var.mpi && var.mpi->in_kernel())var.mpi->begin_exchange(MpiEntity::Element,{{values.data(),columns}});
#endif
}
// Completes an overlapped gather before its second pass; no-op otherwise.
inline void mpi_finish_gather(const Variables& var) {
#ifdef USEMPI
    if(var.mpi && var.mpi->in_kernel())var.mpi->finish_exchange();
#endif
}
#endif
