#ifndef DES_MPI_SIMULATION_HPP
#define DES_MPI_SIMULATION_HPP
#include "parameters.hpp"
#ifndef USEMPI
enum class MpiEntity { Node, Element };
#endif
#ifdef USEMPI
#include <memory>
#include <initializer_list>
#include <utility>

enum class MpiEntity { Node, Element };
enum class MpiOp {
    State, Temperature, PorePressure, StrainRate, Stress, Force, Velocity,
    Coordinate, Rotation, Volume, Dvoldt, Edvoldt, Nmd, Mass, TimeStep,
    OldMeanStress, Friction, VelocityBC, Residual, Surface, SwapVolume,
    Remesh, PhaseChange, Hydrous, Collect, Boundary, Sample, Average, Metrics, Report, Stop, Count
};

// Root keeps the global management/output mesh. Each compute rank, including
// root, advances compact owned/ghost arrays and exchanges only neighbor halos.
class MpiSimulation {
public:
    explicit MpiSimulation(Param&);
    ~MpiSimulation();
    bool is_root() const;
    bool distributed() const;
    bool in_kernel() const;
    int worker_loop();
    void attach(Variables&);
    void synchronize(Variables&);
    bool dispatch(MpiOp, const Variables&, double argument = 0);
    double result() const;
    int owned_nodes() const;
    int owned_elements() const;
    int global_nodes() const;
    void exchange(MpiEntity, std::initializer_list<std::pair<double*,int>>);
    void collect(const Variables&, bool averages = false);
    void boundary(const Variables&);
    void sample(const Variables&, const int_vec& nodes, const int_vec& elements);
    bool metrics(const Variables&, double& speed, double& plastic_rate);
    void minimum(double* data, int count);
    double maximum(double value);
    double sum(double value);
    int sum(int value);
    void maximum(int* data, int count);
    void report();
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    MpiSimulation(const MpiSimulation&) = delete;
    MpiSimulation& operator=(const MpiSimulation&) = delete;
};
#endif

inline int mpi_element_begin(const Variables&) { return 0; }
inline int mpi_node_begin(const Variables&) { return 0; }
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
#endif
