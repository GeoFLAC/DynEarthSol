#ifdef USEMPI
#include <mpi.h>
#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <new>
#include <numeric>
#include <set>
#include <sstream>
#include <utility>
#include "mpi_simulation.hpp"
#include "mpi_partition.hpp"
#include "mpi_local_fields.hpp"
#include "mpi_transport.hpp"
#include "mpi_halo_kernels.hpp"
#include "mpi_gpu_buffer.hpp"
#ifdef ACC
#include <openacc.h>
#endif
#include "bc.hpp"
#include "earthquake_state.hpp"
#include "fields.hpp"
#include "geometry.hpp"
#include "markerset.hpp"
#include "matprops.hpp"
#include "output.hpp"
#include "phasechanges.hpp"
#include "remeshing.hpp"
#include "rheology.hpp"
#include "pseudo_transient.hpp"

namespace {
#ifdef ACC
int selected_device = -1;
void prefetch_managed(void* pointer, std::size_t bytes) {
    if (!pointer || !bytes) return;
    cudaPointerAttributes info{};
    const cudaError_t status = cudaPointerGetAttributes(&info, pointer);
    if (status != cudaSuccess) {
        cudaGetLastError();
        return;
    }
    if (info.type == cudaMemoryTypeManaged)
        mpi_cuda_check(cudaMemPrefetchAsync(pointer, bytes, selected_device, 0));
}
#endif

// A list of field arrays packed into one message, column by column.
struct Fields {
    struct Field {
        double* data;
        int stride, columns;
    };
    std::array<Field,64> list;
    int nfields = 0, columns = 0;
    void add(double_vec* values) {
        if (values) add(values->data(), values->size(), 1);
    }
    template<int N> void add(Array2D<double,N>* values) {
        if (values) add(values->data(), values->size(), N);
    }
    void add(double* data, int rows, int cols) {
        if (nfields == int(list.size())) die(EXIT_INTERNAL_ASSERT, "too many MPI fields");
        list[nfields++] = {data, rows, cols};
        columns += cols;
    }
};

Fields node_fields(Variables& var) {
    Fields fields;
#define ADD_FIELD(x) fields.add(var.x);
    MPI_NODE_VECTORS(ADD_FIELD)
    MPI_NODE_ARRAYS(ADD_FIELD)
#undef ADD_FIELD
    return fields;
}

Fields element_fields(Variables& var) {
    Fields fields;
#define ADD_FIELD(x) fields.add(var.x);
    MPI_ELEMENT_VECTORS(ADD_FIELD)
    MPI_ELEMENT_ARRAYS(ADD_FIELD)
#undef ADD_FIELD
    return fields;
}

void pack(const Fields& fields, const int_vec& rows, double_vec& out) {
    out.resize(std::size_t(rows.size()) * fields.columns);
    std::size_t offset = 0;
    if (rows.empty()) return;
    for (int i=0; i<fields.nfields; ++i) {
        const auto& field = fields.list[i];
        mpi_pack_columns(field.data, field.stride, field.columns, rows.data(), rows.size(), out.data() + offset);
        offset += std::size_t(rows.size()) * field.columns;
    }
}

void unpack(const Fields& fields, const int* rows, int nrows, const double* in) {
    std::size_t offset = 0;
    if (!nrows) return;
    for (int i=0; i<fields.nfields; ++i) {
        const auto& field = fields.list[i];
        mpi_unpack_columns(field.data, field.stride, field.columns, rows, nrows, in + offset);
        offset += std::size_t(nrows) * field.columns;
    }
}

#ifdef ACC
void pack_device(const Fields& fields, const int_vec& rows, MpiGpuBuffer& out) {
    out.resize(std::size_t(rows.size()) * fields.columns);
    std::size_t offset = 0;
    if (rows.empty()) return;
    for (int i=0; i<fields.nfields; ++i) {
        const auto& field = fields.list[i];
        mpi_pack_columns<true>(field.data, field.stride, field.columns,
                               rows.data(), rows.size(), out.device() + offset);
        offset += std::size_t(rows.size()) * field.columns;
    }
    out.download();
}

void unpack_device(const Fields& fields, const int_vec& rows, MpiGpuBuffer& in) {
    if (rows.empty()) return;
    in.upload();
    std::size_t offset = 0;
    for (int i=0; i<fields.nfields; ++i) {
        const auto& field = fields.list[i];
        mpi_unpack_columns<true>(field.data, field.stride, field.columns,
                                 rows.data(), rows.size(), in.device() + offset);
        offset += std::size_t(rows.size()) * field.columns;
    }
}
#endif

int count(std::size_t entries) {
    if (entries > INT_MAX) die(EXIT_RUNTIME_RESOURCE, "MPI message exceeds INT_MAX entries");
    return entries;
}

void metadata(MpiState& state, Variables& var) {
#define ARCHIVE(x) state.value(var.x)
    ARCHIVE(time);
    ARCHIVE(last_remesh_time);
    ARCHIVE(dt);
    ARCHIVE(steps);
    ARCHIVE(nremesh);
    ARCHIVE(max_vbc_val);
    ARCHIVE(max_global_vel_mag);
    ARCHIVE(global_dt_min);
    ARCHIVE(compensation_pressure);
    ARCHIVE(bottom_temperature);
    ARCHIVE(nnode);
    ARCHIVE(nelem);
    ARCHIVE(vbc_val_z1_loading_period);
    ARCHIVE(PT_initial_residual_max);
    ARCHIVE(PT_skip_scale_update);
    ARCHIVE(PT_force_residual_valid);
    ARCHIVE(info_display_next_step);
    for (int i=0; i<nbdrytypes; ++i) {
        state.value(var.vbc_types[i]);
        state.value(var.vbc_values[i]);
        state.value(var.stress_bc_values[i]);
        state.value(var.bnodes[i]);
        state.value(var.bfacets[i]);
    }
    for (int i=0; i<nbdrytypes_hydro; ++i) {
        state.value(var.hbc_types[i]);
        state.value(var.stress_bc_types[i]);
    }
    for (int i=0; i<nbdrytypes*nbdrytypes; ++i) state.value(var.edge_slot[i]);
    ARCHIVE(bnormals);
    ARCHIVE(bcflag);
    ARCHIVE(edge_vec);
    ARCHIVE(vbc_vertical_div_x0);
    ARCHIVE(vbc_vertical_div_x1);
    ARCHIVE(vbc_vertical_ratio_x0);
    ARCHIVE(vbc_vertical_ratio_x1);
    state.value(var.surfinfo.max_surf_vel);
#undef ARCHIVE
}

void local_archive(MpiState& state, Variables& var) {
    metadata(state, var);
#define ARCHIVE(x) state.value(var.x);
    MPI_NODE_VECTORS(ARCHIVE)
    MPI_NODE_ARRAYS(ARCHIVE)
    MPI_ELEMENT_VECTORS(ARCHIVE)
    MPI_ELEMENT_ARRAYS(ARCHIVE)
    ARCHIVE(connectivity)
    ARCHIVE(elemmarkers)
    ARCHIVE(top_elems)
    ARCHIVE(ntop_elems)
    ARCHIVE(surfinfo.top_nodes)
    ARCHIVE(surfinfo.ntop)
#undef ARCHIVE
    state.value(var.support.arr_data);
    state.value(var.support.idx_data);
    state.value(var.support.lidx_data);
    var.support.rebind();
}

void release(Variables& var) {
    delete var.mat;
    var.mat = nullptr;
#define RELEASE(x)        \
    delete var.x;         \
    var.x = nullptr;
    MPI_NODE_VECTORS(RELEASE)
    MPI_NODE_ARRAYS(RELEASE)
    MPI_ELEMENT_VECTORS(RELEASE)
    MPI_ELEMENT_ARRAYS(RELEASE)
    RELEASE(connectivity)
    RELEASE(elemmarkers)
    RELEASE(bnormals)
    RELEASE(bcflag)
    RELEASE(top_elems)
    RELEASE(surfinfo.top_nodes)
    RELEASE(ntmp)
    RELEASE(etmp)
    RELEASE(etmp_int)
    RELEASE(tmp_result)
#undef RELEASE
    for (int i=0; i<nbdrytypes; ++i) {
        delete var.bnodes[i];
        delete var.bfacets[i];
        var.bnodes[i] = nullptr;
        var.bfacets[i] = nullptr;
    }
}

template<class T> void select(const std::vector<T>* src, std::vector<T>*& dst, const int_vec& ids) {
    if (!src) {
        delete dst;
        dst = nullptr;
        return;
    }
    if (!dst) dst = new std::vector<T>;
    dst->resize(ids.size());
    for (int i=0; i<int(ids.size()); ++i) (*dst)[i] = (*src)[ids[i]];
}

template<class T, int N> void select(const Array2D<T,N>* src, Array2D<T,N>*& dst, const int_vec& ids) {
    if (!src) {
        delete dst;
        dst = nullptr;
        return;
    }
    if (!dst) dst = new Array2D<T,N>;
    dst->resize(ids.size(), false);
    for (int d=0; d<N; ++d) {
        for (int i=0; i<int(ids.size()); ++i) (*dst)[i][d] = (*src)[ids[i]][d];
    }
}

template<class T, int N> void copy(const Array2D<T,N>* src, Array2D<T,N>*& dst) {
    if (!src) {
        delete dst;
        dst = nullptr;
        return;
    }
    if (!dst) dst = new Array2D<T,N>;
    dst->resize(src->size(), false);
    std::copy(src->data(), src->data() + src->num_elements(), dst->data());
}

// Build the compact state of one partition from the root's global mesh.
void make_local(const Variables& global, const MpiPartition& partition, Variables& var) {
#define COPY(x) var.x = global.x
    COPY(time);
    COPY(last_remesh_time);
    COPY(dt);
    COPY(steps);
    COPY(nremesh);
    COPY(max_vbc_val);
    COPY(max_global_vel_mag);
    COPY(global_dt_min);
    COPY(compensation_pressure);
    COPY(bottom_temperature);
    COPY(PT_initial_residual_max);
    COPY(PT_skip_scale_update);
    COPY(PT_force_residual_valid);
    COPY(info_display_next_step);
    COPY(vbc_val_z1_loading_period);
    COPY(edge_vec);
    COPY(vbc_vertical_div_x0);
    COPY(vbc_vertical_div_x1);
    COPY(vbc_vertical_ratio_x0);
    COPY(vbc_vertical_ratio_x1);
    COPY(surfinfo.max_surf_vel);
#undef COPY
    var.nnode = partition.nodes.size();
    var.nelem = partition.elements.size();
#define SELECT_NODES(x) select(global.x, var.x, partition.nodes);
    MPI_NODE_VECTORS(SELECT_NODES)
    MPI_NODE_ARRAYS(SELECT_NODES)
#undef SELECT_NODES
#define SELECT_ELEMENTS(x) select(global.x, var.x, partition.elements);
    MPI_ELEMENT_VECTORS(SELECT_ELEMENTS)
    MPI_ELEMENT_ARRAYS(SELECT_ELEMENTS)
    SELECT_ELEMENTS(connectivity)
    SELECT_ELEMENTS(elemmarkers)
#undef SELECT_ELEMENTS
    for (int e=0; e<var.nelem; ++e) {
        for (int j=0; j<NODES_PER_ELEM; ++j)
            (*var.connectivity)[e][j] = partition.node_local.at((*var.connectivity)[e][j]);
    }
    select(global.bcflag, var.bcflag, partition.nodes);
    copy(global.bnormals, var.bnormals);
    for (int i=0; i<nbdrytypes; ++i) {
        var.vbc_types[i] = global.vbc_types[i];
        var.vbc_values[i] = global.vbc_values[i];
        var.stress_bc_values[i] = global.stress_bc_values[i];
        if (!var.bnodes[i]) var.bnodes[i] = new int_vec;
        var.bnodes[i]->clear();
        for (int n : *global.bnodes[i]) {
            auto at = partition.node_local.find(n);
            if (at != partition.node_local.end()) var.bnodes[i]->push_back(at->second);
        }
        if (!var.bfacets[i]) var.bfacets[i] = new std::vector<std::pair<int,int>>;
        var.bfacets[i]->clear();
        for (auto facet : *global.bfacets[i]) {
            auto at = partition.element_local.find(facet.first);
            if (at != partition.element_local.end()) var.bfacets[i]->emplace_back(at->second, facet.second);
        }
    }
    for (int i=0; i<nbdrytypes_hydro; ++i) {
        var.hbc_types[i] = global.hbc_types[i];
        var.stress_bc_types[i] = global.stress_bc_types[i];
    }
    std::copy(global.edge_slot, global.edge_slot + nbdrytypes*nbdrytypes, var.edge_slot);
    // Copy each support in GLOBAL element order, including full stars of owned nodes.
    var.support.idx_data.assign(var.nnode+1, 0);
    var.support.arr_data.clear();
    var.support.lidx_data.clear();
    for (int n=0; n<var.nnode; ++n) {
        int global_node = partition.nodes[n];
        const int* patch = global.support.patch(global_node);
        const int* local_index = global.support.local(global_node);
        for (int i=0; i<global.support.size(global_node); ++i) {
            auto at = partition.element_local.find(patch[i]);
            if (at != partition.element_local.end()) {
                var.support.arr_data.push_back(at->second);
                var.support.lidx_data.push_back(local_index[i]);
            }
            else if (n < partition.owned_nodes) die(EXIT_INTERNAL_ASSERT, "owned node has an incomplete support");
        }
        var.support.idx_data[n+1] = var.support.arr_data.size();
    }
    var.support.rebind();
    if (!var.top_elems) var.top_elems = new int_vec;
    var.top_elems->clear();
    if (global.top_elems) {
        for (int id : *global.top_elems) {
            auto at = partition.element_local.find(id);
            if (at != partition.element_local.end()) var.top_elems->push_back(at->second);
        }
    }
    var.ntop_elems = var.top_elems->size();
    if (!var.surfinfo.top_nodes) var.surfinfo.top_nodes = new int_vec;
    var.surfinfo.top_nodes->clear();
    if (global.surfinfo.top_nodes) {
        for (int id : *global.surfinfo.top_nodes) {
            auto at = partition.node_local.find(id);
            if (at != partition.node_local.end()) var.surfinfo.top_nodes->push_back(at->second);
        }
    }
    var.surfinfo.ntop = var.surfinfo.top_nodes->size();
}

const char* op_name(MpiOp op) {
    switch (op) {
#define OP_NAME(x) case MpiOp::x: return #x
        OP_NAME(State);
        OP_NAME(Temperature);
        OP_NAME(PorePressure);
        OP_NAME(StrainRate);
        OP_NAME(Stress);
        OP_NAME(Force);
        OP_NAME(Velocity);
        OP_NAME(Coordinate);
        OP_NAME(Rotation);
        OP_NAME(Volume);
        OP_NAME(Dvoldt);
        OP_NAME(Edvoldt);
        OP_NAME(Nmd);
        OP_NAME(Mass);
        OP_NAME(TimeStep);
        OP_NAME(OldMeanStress);
        OP_NAME(VelocityBC);
        OP_NAME(Residual);
        OP_NAME(PhysicalPT);
        OP_NAME(InitialPT);
        OP_NAME(Surface);
        OP_NAME(SurfaceIdle);
        OP_NAME(SwapVolume);
        OP_NAME(Remesh);
        OP_NAME(PhaseChange);
        OP_NAME(Hydrous);
        OP_NAME(Collect);
        OP_NAME(Geometry);
        OP_NAME(Boundary);
        OP_NAME(Sample);
        OP_NAME(Average);
        OP_NAME(Metrics);
        OP_NAME(Report);
        OP_NAME(Stop);
        OP_NAME(Count);
#undef OP_NAME
    }
    return "unknown";
}

// Management operations run on root's global mesh, not in compact kernels.
bool management(MpiOp op) {
    return op == MpiOp::State || op == MpiOp::Remesh || op == MpiOp::Surface || op == MpiOp::PhaseChange ||
           op == MpiOp::Hydrous || op == MpiOp::Report || op == MpiOp::Collect || op == MpiOp::Geometry ||
           op == MpiOp::Boundary || op == MpiOp::Sample;
}

// Root broadcasts one command per operation; workers mirror its scalars.
struct Command {
    int operation, steps;
    double time, dt, max_vbc, max_vel, global_dt, surface_vel, argument;
    bool hydraulic, body_adjustment;
};
}

int mpi_launcher_rank() {
    for (const char* name : {"OMPI_COMM_WORLD_RANK", "PMIX_RANK", "PMI_RANK", "SLURM_PROCID"}) {
        const char* value = std::getenv(name);
        if (value) return std::atoi(value);
    }
    return 0;
}

void mpi_select_device() {
#ifdef ACC
    if (selected_device >= 0) return;
    const char* local = std::getenv("OMPI_COMM_WORLD_LOCAL_RANK");
    if (!local) local = std::getenv("MPI_LOCALRANKID");
    if (!local) local = std::getenv("PMI_LOCAL_RANK");
    const int index = local ? std::atoi(local) : 0;
    const int devices = acc_get_num_devices(acc_device_nvidia);
    if (index < 0 || index >= devices)
        die(EXIT_RUNTIME_RESOURCE, "GPU MPI requires one visible NVIDIA device per local rank");
    acc_set_device_num(index, acc_device_nvidia);
    acc_init(acc_device_nvidia);
    selected_device = index;
#endif
}

struct MpiSimulation::Impl {
    explicit Impl(Param& parameters) : param(parameters) {
        local.func_time = {};
        packing.func_time = {};
    }
    ~Impl() {
        release(local);
        release(packing);
    }
    Param& param;
    Variables* master = nullptr;
    Variables local{}, packing{};
    MpiPartition part;
    std::vector<MpiPartition> parts; // only rank zero retains global partition maps
    struct Buffers {
        double_vec send, recv;
#ifdef ACC
        MpiGpuBuffer gpu_send, gpu_recv;
#endif
    };
    std::vector<Buffers> peer;
    std::vector<MPI_Request> requests;
    std::vector<MpiPeerBuffer> transfers;
    int rank = 0, size = 1, global_n = 0, global_e = 0, generation = 0;
    bool active = false, inside = false;
    MpiOp current = MpiOp::State;
    double scalar = 0, measurements[2] = {};
    PTResult pt;
    std::uint64_t version = 0, root_version = 0, surface_version = ~std::uint64_t(0);
    std::uint64_t calls[int(MpiOp::Count)] = {}, payload_bytes = 0, state_bytes = 0, root_bytes = 0, all_payload_bytes = 0;
    double elapsed[int(MpiOp::Count)] = {}, communication_s = 0;
    double operation_comm[int(MpiOp::Count)] = {}, wait_s = 0, interior_window_s = 0;
    double_vec rank_profile;
    bool overlap_enabled = true, pending_active = false;
    int pending = 0;
    bool pending_nodes = false;
    Fields pending_fields;
    double posted_at = 0;
    std::uint64_t geometry_version = ~std::uint64_t(0), overlap_calls = 0;
    std::vector<unsigned char> bytes;
    double_vec send, recv, scatter;
    int_vec counts, offsets, queries_n, queries_e;
    array_t average_coord;
    tensor_t average_strain, average_stress;
    double_vec average_plastic;
    double average_time = 0;

    void command(const Command& received) {
        local.steps = received.steps;
        local.time = received.time;
        local.dt = received.dt;
        local.max_vbc_val = received.max_vbc;
        local.max_global_vel_mag = received.max_vel;
        local.global_dt_min = received.global_dt;
        local.surfinfo.max_surf_vel = received.surface_vel;
        param.control.has_hydraulic_diffusion = received.hydraulic;
        param.ic.has_body_force_adjustment = received.body_adjustment;
    }

    void trace() {
        if (rank != 0) return;
        const char* path = std::getenv("DES_MPI_PARTITION_TRACE");
        std::ofstream file;
        if (path) {
            std::ostringstream name;
            name << path << "." << std::setw(4) << std::setfill('0') << generation << ".csv";
            file.open(name.str());
            if (!file) die(EXIT_IO_OPEN, "cannot create MPI partition trace");
            file << "entity,id,owner,x,z";
#ifdef THREED
            file << ",y";
#endif
            for (int j=0; j<NODES_PER_ELEM; ++j) file << ",n" << j;
            file << "\n" << std::setprecision(17);
        }
        for (int r=0; r<size; ++r) {
            const auto& partition = parts[r];
            std::cout << "[mpi][partition] generation=" << generation << " step=" << master->steps << " rank=" << r
                      << " owned_nodes=" << partition.owned_nodes << " owned_elements=" << partition.owned_elements
                      << " owned_markers=" << partition.owned_markers
                      << " ghost_nodes=" << partition.nodes.size() - partition.owned_nodes
                      << " ghost_elements=" << partition.elements.size() - partition.owned_elements
                      << " neighbors=" << partition.halos.size()
                      << " interior_nodes=" << partition.interior_nodes.size()
                      << " interface_nodes=" << partition.interface_nodes.size() << '\n';
            if (!file) continue;
            for (int i=0; i<partition.owned_nodes; ++i) {
                int n = partition.nodes[i];
                file << "node," << n << ',' << r << ',' << (*master->coord)[n][0] << ',' << (*master->coord)[n][NDIMS-1];
#ifdef THREED
                file << ',' << (*master->coord)[n][1];
#endif
                for (int j=0; j<NODES_PER_ELEM; ++j) file << ',';
                file << '\n';
            }
            for (int i=0; i<partition.owned_elements; ++i) {
                int e = partition.elements[i];
                double center[NDIMS] = {};
                for (int j=0; j<NODES_PER_ELEM; ++j) {
                    for (int d=0; d<NDIMS; ++d)
                        center[d] += (*master->coord)[(*master->connectivity)[e][j]][d] / NODES_PER_ELEM;
                }
                file << "element," << e << ',' << r << ',' << center[0] << ',' << center[NDIMS-1];
#ifdef THREED
                file << ',' << center[1];
#endif
                for (int j=0; j<NODES_PER_ELEM; ++j) file << ',' << (*master->connectivity)[e][j];
                file << '\n';
            }
        }
    }

    void distribute(bool rebalance) {
        #pragma acc wait
        const double start = MPI_Wtime();
        if (rank == 0) {
            // Only root management changes marker composition. Refresh here,
            // before rebuilding the compact caches from the new marker counts.
            master->mat->refresh_elem_cache();
            // Marker operations stay on the management mesh. Compute ranks
            // need material counts, not copies of the individual records.
            if (rebalance && !master->markersets.empty())
                master->markersets[0]->check_marker_elem_consistency(*master);
            // Serial remeshing leaves this inactive hydraulic scratch array at
            // its old size. Compact selection/collection still visits it.
            if (!param.control.has_hydraulic_diffusion && master->old_mean_stress)
                master->old_mean_stress->resize(master->nelem, 0.0);
            global_n = master->nnode;
            global_e = master->nelem;
            if (rebalance || parts.empty()) {
                parts = mpi_partition_mesh(*master, size);
                ++generation;
                trace();
            }
        }
        int extent[] = {global_n, global_e};
        MPI_Bcast(extent, 2, MPI_INT, 0, MPI_COMM_WORLD);
        global_n = extent[0];
        global_e = extent[1];
        for (int r=0; r<size; ++r) {
            std::uint64_t nbytes = 0;
            if (rank == 0 && r) {
                // Root builds, serializes and sends the compact state of rank r.
                make_local(*master, parts[r], packing);
                MpiState state(bytes, true);
                parts[r].archive(state);
                local_archive(state, packing);
                nbytes = bytes.size();
                MPI_Send(&nbytes, 1, MPI_UINT64_T, r, 20, MPI_COMM_WORLD);
                MPI_Send(bytes.data(), count(nbytes), MPI_BYTE, r, 21, MPI_COMM_WORLD);
                state_bytes += nbytes;
            }
            if (rank == r) {
                if (r) {
                    MPI_Recv(&nbytes, 1, MPI_UINT64_T, 0, 20, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    bytes.resize(nbytes);
                    MPI_Recv(bytes.data(), count(nbytes), MPI_BYTE, 0, 21, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MpiState state(bytes, false);
                    part.archive(state);
                    local_archive(state, local);
                    state.finish();
                } else {
                    part = parts[0];
                    make_local(*master, part, local);
                }
                if (!local.ntmp) local.ntmp = new double_vec;
                if (!local.etmp) local.etmp = new double_vec;
                if (!local.etmp_int) local.etmp_int = new int_vec;
                if (!local.tmp_result) local.tmp_result = new elem_cache;
                local.ntmp->resize(local.nnode);
                local.etmp->resize(local.nelem);
                local.etmp_int->resize(local.nelem);
                local.tmp_result->resize(local.nelem, false);
                if (!local.mat) local.mat = new MatProps(param, local);
                else local.mat->refresh_elem_cache();
                local.mpi = master ? master->mpi : local.mpi;
            }
        }
        double hydro_diff_max = rank == 0 ? master->mat->hydro_diff_max : 0;
        MPI_Bcast(&hydro_diff_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        local.mat->hydro_diff_max = hydro_diff_max;
        peer.resize(part.halos.size());
        requests.resize(2*part.halos.size());
        transfers.resize(part.halos.size());
#ifdef ACC
        // Host MMG/management may resize every array. Rebind the compact state
        // first, then place managed payloads on this rank's GPU before timesteps.
        Fields local_node_fields = node_fields(local);
        Fields local_element_fields = element_fields(local);
        for (const Fields* fields : {&local_node_fields, &local_element_fields}) {
            for (int i=0; i<fields->nfields; ++i) {
                const auto& field = fields->list[i];
                prefetch_managed(field.data, std::size_t(field.stride) * field.columns * sizeof(double));
            }
        }
        prefetch_managed(local.connectivity->data(), local.connectivity->num_elements() * sizeof(int));
        prefetch_managed(local.bcflag->data(), local.bcflag->size() * sizeof(uint));
        for (const int_vec* ids : {&local.support.arr_data, &local.support.idx_data, &local.support.lidx_data,
                                   &part.interior_nodes, &part.interface_nodes, local.top_elems, local.surfinfo.top_nodes})
            prefetch_managed(const_cast<int*>(ids->data()), ids->size() * sizeof(int));
        mpi_cuda_check(cudaDeviceSynchronize());
#endif
        // Redistribution finishes a mutating command; root stays current when
        // execute() advances the command version on return.
        root_version = version + 1;
        surface_version = geometry_version = ~std::uint64_t(0);
        communication_s += MPI_Wtime() - start;
    }

    // Collect owned rows (all, surface only, or queried IDs) of each rank onto root.
    void gather_rows(bool nodes, const Fields& from, const Fields& to, bool surface=false, const int_vec* query=nullptr) {
        #pragma acc wait
        int_vec rows, global_rows;
        const int_vec* local_rows = nodes ? &part.node_rows : &part.element_rows;
        if (surface) local_rows = nodes ? &part.surface_owned_nodes : &part.surface_owned_elements;
        if (query) {
            const auto& local_ids = nodes ? part.node_local : part.element_local;
            const int owned = nodes ? part.owned_nodes : part.owned_elements;
            for (int id : *query) {
                auto at = local_ids.find(id);
                if (at != local_ids.end() && at->second < owned) rows.push_back(at->second);
            }
            local_rows = &rows;
        }
        pack(from, *local_rows, send);
        if (rank == 0) {
            counts.assign(size, 0);
            offsets.assign(size, 0);
            for (int r=0; r<size; ++r) {
                const auto& partition = parts[r];
                int nrows = nodes ? partition.owned_nodes : partition.owned_elements;
                if (surface && !query) {
                    const auto& selected = nodes ? partition.surface_owned_nodes : partition.surface_owned_elements;
                    const auto& ids = nodes ? partition.nodes : partition.elements;
                    nrows = selected.size();
                    for (int i : selected) global_rows.push_back(ids[i]);
                }
                if (query) {
                    nrows = 0;
                    const auto& local_ids = nodes ? partition.node_local : partition.element_local;
                    const int owned = nodes ? partition.owned_nodes : partition.owned_elements;
                    for (int id : *query) {
                        auto at = local_ids.find(id);
                        if (at != local_ids.end() && at->second < owned) {
                            ++nrows;
                            global_rows.push_back(id);
                        }
                    }
                }
                counts[r] = count(std::size_t(nrows) * from.columns);
                if (r) offsets[r] = offsets[r-1] + counts[r-1];
            }
            recv.resize(std::size_t(offsets.back()) + counts.back());
        }
        double start = MPI_Wtime();
        MPI_Gatherv(send.data(), count(send.size()), MPI_DOUBLE, recv.data(), counts.data(), offsets.data(),
                    MPI_DOUBLE, 0, MPI_COMM_WORLD);
        communication_s += MPI_Wtime() - start;
        std::size_t row_start = 0;
        if (rank == 0) {
            for (int r=0; r<size; ++r) {
                const auto& ids = nodes ? parts[r].nodes : parts[r].elements;
                const int nrows = counts[r] / from.columns;
                if (nrows) unpack(to, surface || query ? global_rows.data() + row_start : ids.data(), nrows,
                                  recv.data() + offsets[r]);
                row_start += nrows;
                if (r) root_bytes += std::uint64_t(counts[r]) * sizeof(double);
            }
        }
    }

    // Send root's surface rows back to every rank that holds them.
    void scatter_surface(bool nodes, const Fields& from, const Fields& to) {
        #pragma acc wait
        const auto& local_rows = nodes ? part.surface_nodes : part.surface_elements;
        if (rank == 0) {
            counts.assign(size, 0);
            offsets.assign(size, 0);
            scatter.clear();
            for (int r=0; r<size; ++r) {
                const auto& partition = parts[r];
                const auto& rows = nodes ? partition.surface_nodes : partition.surface_elements;
                const auto& ids = nodes ? partition.nodes : partition.elements;
                int_vec global_rows;
                for (int i : rows) global_rows.push_back(ids[i]);
                pack(from, global_rows, send);
                offsets[r] = scatter.size();
                counts[r] = count(send.size());
                scatter.insert(scatter.end(), send.begin(), send.end());
                if (r) root_bytes += send.size() * sizeof(double);
            }
        }
        recv.resize(std::size_t(local_rows.size()) * to.columns);
        double start = MPI_Wtime();
        MPI_Scatterv(scatter.data(), counts.data(), offsets.data(), MPI_DOUBLE, recv.data(), count(recv.size()),
                     MPI_DOUBLE, 0, MPI_COMM_WORLD);
        communication_s += MPI_Wtime() - start;
        unpack(to, local_rows.data(), local_rows.size(), recv.data());
    }

    void full() {
        if (root_version == version) return;
        Fields root_node_fields, root_element_fields;
        if (rank == 0) {
            root_node_fields = node_fields(*master);
            root_element_fields = element_fields(*master);
        }
        gather_rows(true, node_fields(local), root_node_fields);
        gather_rows(false, element_fields(local), root_element_fields);
        if (rank == 0) master->mat->hydro_diff_max = local.mat->hydro_diff_max;
        root_version = version;
    }

    void geometry() {
        if (geometry_version == version || root_version == version) return;
        Fields root_node_fields, root_element_fields, local_node_fields, local_element_fields;
        local_node_fields.add(local.coord);
        local_element_fields.add(local.volume);
        if (rank == 0) {
            root_node_fields.add(master->coord);
            root_element_fields.add(master->volume);
        }
        gather_rows(true, local_node_fields, root_node_fields);
        gather_rows(false, local_element_fields, root_element_fields);
        geometry_version = version;
    }

    Fields surface_node_fields(Variables& var) {
        Fields fields;
        fields.add(var.coord);
        fields.add(var.volume_n);
        return fields;
    }

    Fields surface_element_fields(Variables& var) {
        Fields fields;
        fields.add(var.volume);
        fields.add(var.plstrain);
        fields.add(var.stress);
        fields.add(var.strain);
        fields.add(var.strain_rate);
        return fields;
    }

    void boundary() {
        if (surface_version == version || root_version == version) return;
        Fields root_node_fields, root_element_fields;
        if (rank == 0) {
            root_node_fields = surface_node_fields(*master);
            root_element_fields = surface_element_fields(*master);
        }
        gather_rows(true, surface_node_fields(local), root_node_fields, true);
        gather_rows(false, surface_element_fields(local), root_element_fields, true);
        surface_version = version;
    }

    void query() {
        if (root_version == version) return;
        if (rank == 0) {
            MpiState state(bytes, true);
            state.value(queries_n);
            state.value(queries_e);
        }
        std::uint64_t nbytes = bytes.size();
        MPI_Bcast(&nbytes, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        if (rank) bytes.resize(nbytes);
        MPI_Bcast(bytes.data(), count(nbytes), MPI_BYTE, 0, MPI_COMM_WORLD);
        if (rank) {
            MpiState state(bytes, false);
            state.value(queries_n);
            state.value(queries_e);
            state.finish();
        }
        Fields root_node_fields, root_element_fields;
        if (rank == 0) {
            root_node_fields = node_fields(*master);
            root_element_fields = element_fields(*master);
        }
        gather_rows(true, node_fields(local), root_node_fields, false, &queries_n);
        gather_rows(false, element_fields(local), root_element_fields, false, &queries_e);
    }

    void collect_averages() {
        if (param.mesh.quality_check_step_interval <= 1) return;
        Fields root_node_fields, root_element_fields, local_node_fields, local_element_fields;
        local_node_fields.add(&average_coord);
        local_element_fields.add(&average_strain);
        local_element_fields.add(&average_stress);
        local_element_fields.add(&average_plastic);
        if (rank == 0) {
            Output& output = *master->output;
            output.time0 = average_time;
            output.coord0.resize(global_n, false);
            output.strain0.resize(global_e, false);
            output.stress_avg.resize(global_e, false);
            output.delta_plstrain_avg.resize(global_e);
            root_node_fields.add(&output.coord0);
            root_element_fields.add(&output.strain0);
            root_element_fields.add(&output.stress_avg);
            root_element_fields.add(&output.delta_plstrain_avg);
        }
        gather_rows(true, local_node_fields, root_node_fields);
        gather_rows(false, local_element_fields, root_element_fields);
    }

    void execute(MpiOp op, double arg) {
        current = op;
        inside = true;
        double started = MPI_Wtime();
        double comm_before = communication_s;
        Variables& var = local;
        MpiAccScope offload(!management(op));
        switch (op) {
        case MpiOp::State:
            distribute(true);
            break;
        case MpiOp::Temperature:
            update_temperature(param, var, *var.temperature, *var.tmp_result);
            exchange(MpiEntity::Node, {{var.temperature->data(), 1}});
            break;
        case MpiOp::PorePressure:
            update_pore_pressure(param, var, *var.ppressure, *var.dppressure, *var.ntmp, *var.tmp_result,
                                 *var.stress, *var.old_mean_stress);
            exchange(MpiEntity::Node, {{var.ppressure->data(), 1}, {var.dppressure->data(), 1}});
            break;
        case MpiOp::StrainRate:
            update_strain_rate(var, *var.strain_rate);
            break;
        case MpiOp::Dvoldt:
            compute_dvoldt(var, *var.ntmp, *var.etmp);
            exchange(MpiEntity::Node, {{var.ntmp->data(), 1}});
            break;
        case MpiOp::Edvoldt:
            compute_edvoldt(var, *var.ntmp, *var.edvoldt);
            break;
        case MpiOp::Stress:
            update_stress(param, var, *var.stress, *var.stressyy, *var.dpressure, *var.viscosity, *var.strain,
                          *var.plstrain, *var.delta_plstrain, *var.strain_rate, *var.ppressure, *var.dppressure,
                          *var.vel, *var.dyn_fric_coeff, *var.state_variable);
            break;
        case MpiOp::Nmd:
            NMD_stress(var, *var.stress, *var.ntmp, *var.etmp);
            break;
        case MpiOp::Force:
            update_force(param, var, *var.force, *var.force_residual, *var.tmp_result);
            break;
        case MpiOp::Velocity:
            update_velocity(var, *var.vel);
            exchange(MpiEntity::Node, {{var.vel->data(), NDIMS}});
            break;
        case MpiOp::VelocityBC:
            apply_vbcs(param, var, *var.vel);
            exchange(MpiEntity::Node, {{var.vel->data(), NDIMS}});
            break;
        case MpiOp::Coordinate:
            update_coordinate(var, *var.coord);
            exchange(MpiEntity::Node, {{var.coord->data(), NDIMS}});
            break;
        case MpiOp::Rotation:
            rotate_stress(var, *var.stress, *var.strain, arg);
            break;
        case MpiOp::Volume:
            compute_volume(var, *var.volume);
            exchange(MpiEntity::Element, {{var.volume->data(), 1}});
            break;
        case MpiOp::Mass:
            compute_mass(param, var, arg, *var.volume_n, *var.mass, *var.tmass, *var.hmass, *var.ymass,
                         *var.tmp_result);
            exchange(MpiEntity::Node, {{var.volume_n->data(), 1}, {var.mass->data(), 1}, {var.tmass->data(), 1},
                                       {var.hmass->data(), 1}, {var.ymass->data(), 1}});
            break;
        case MpiOp::OldMeanStress:
            update_old_mean_stress(param, var, *var.stress, *var.old_mean_stress);
            break;
        case MpiOp::TimeStep:
            scalar = compute_dt(param, var, arg != 0);
            if (rank == 0) {
                master->max_global_vel_mag = var.max_global_vel_mag;
                master->global_dt_min = var.global_dt_min;
            }
            break;
        case MpiOp::Residual:
            scalar = calculate_residual_force(var, *var.force_residual);
            if (rank == 0) master->l2_residual = scalar;
            break;
        case MpiOp::PhysicalPT:
        case MpiOp::InitialPT:
            pt = op == MpiOp::InitialPT ? run_initial_equilibrium_pt(param, var) : run_physical_step_pt(param, var);
            if (rank == 0) {
                master->l2_residual = var.l2_residual;
                master->PT_initial_residual_max = var.PT_initial_residual_max;
                master->PT_skip_scale_update = var.PT_skip_scale_update;
                master->PT_force_residual_valid = var.PT_force_residual_valid;
            }
            break;
        case MpiOp::SwapVolume:
            std::swap(var.volume, var.volume_old);
            if (rank == 0) std::swap(master->volume, master->volume_old);
            break;
        case MpiOp::SurfaceIdle:
            // Even without erosion, the original surface path corrects stress
            // for the volume change since the preceding coordinate update.
            correct_surface_element(var, *var.volume, *var.volume_n, *var.stress, *var.strain, *var.strain_rate,
                                    *var.plstrain);
            exchange(MpiEntity::Element, {{var.volume->data(), 1}, {var.plstrain->data(), 1},
                {var.stress->data(), NSTR}, {var.strain->data(), NSTR}, {var.strain_rate->data(), NSTR}});
            exchange(MpiEntity::Node, {{var.volume_n->data(), 1}});
            var.surfinfo.max_surf_vel = 0;
            break;
        case MpiOp::Surface: {
            bool whole = var.steps % param.mesh.quality_check_step_interval == 0 ||
                (param.control.surface_process_option != 0 && param.control.surface_process_option != 1);
            if (whole) full();
            else boundary();
            if (rank == 0)
                surface_processes(param, *master, *master->coord, *master->stress, *master->strain,
                                  *master->strain_rate, *master->plstrain, *master->volume, *master->volume_n,
                                  master->surfinfo, master->markersets, *master->elemmarkers,
                                  *master->markers_in_elem);
            if (whole) distribute(false);
            else {
                Fields root_node_fields, root_element_fields;
                if (rank == 0) {
                    root_node_fields = surface_node_fields(*master);
                    root_element_fields = surface_element_fields(*master);
                }
                scatter_surface(true, root_node_fields, surface_node_fields(var));
                scatter_surface(false, root_element_fields, surface_element_fields(var));
            }
            break;
        }
        case MpiOp::Remesh:
            full();
            if (rank == 0) remesh(param, *master, int(arg));
            distribute(true);
            break;
        case MpiOp::PhaseChange:
            full();
            if (rank == 0) phase_changes(param, *master);
            distribute(false);
            break;
        case MpiOp::Hydrous:
            full();
            if (rank == 0)
                advect_hydrous_markers(param, *master, arg, *master->markersets[master->hydrous_marker_index],
                                       *master->hydrous_elemmarkers);
            distribute(false);
            break;
        case MpiOp::Collect:
            full();
            if (arg != 0) collect_averages();
            break;
        case MpiOp::Geometry:
            geometry();
            break;
        case MpiOp::Boundary:
            boundary();
            break;
        case MpiOp::Sample:
            query();
            break;
        case MpiOp::Average:
            Output::average_fields(var, param.mesh.quality_check_step_interval,
                                   average_time, average_coord, average_strain, average_stress, average_plastic);
            break;
        case MpiOp::Metrics:
            compute_earthquake_metrics(var, measurements[0], measurements[1]);
            break;
        case MpiOp::Report: {
            MPI_Reduce(&payload_bytes, &all_payload_bytes, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
            double sample[3*int(MpiOp::Count)+3];
            for (int i=0; i<int(MpiOp::Count); ++i) {
                sample[3*i] = double(calls[i]);
                sample[3*i+1] = elapsed[i];
                sample[3*i+2] = operation_comm[i];
            }
            sample[3*int(MpiOp::Count)] = wait_s;
            sample[3*int(MpiOp::Count)+1] = interior_window_s;
            sample[3*int(MpiOp::Count)+2] = double(overlap_calls);
            if (rank == 0) rank_profile.resize(size*(3*int(MpiOp::Count)+3));
            MPI_Gather(sample, 3*int(MpiOp::Count)+3, MPI_DOUBLE, rank_profile.data(), 3*int(MpiOp::Count)+3,
                       MPI_DOUBLE, 0, MPI_COMM_WORLD);
            break;
        }
        default:
            die(EXIT_INTERNAL_ASSERT, "unknown MPI operation");
        }
        // Every command boundary permits host-side scalar/cache changes, output
        // or state redistribution. No outstanding GPU producer may cross it.
        #pragma acc wait
        // Scalar diagnostics do not invalidate current copies of field arrays.
        switch (op) {
        case MpiOp::Collect:
        case MpiOp::Geometry:
        case MpiOp::Boundary:
        case MpiOp::Sample:
        case MpiOp::Average:
        case MpiOp::Metrics:
        case MpiOp::TimeStep:
        case MpiOp::Residual:
        case MpiOp::Report:
            break;
        default:
            ++version;
        }
        ++calls[int(op)];
        elapsed[int(op)] += MPI_Wtime() - started;
        operation_comm[int(op)] += communication_s - comm_before;
        inside = false;
    }

    void begin_exchange(MpiEntity entity, std::initializer_list<std::pair<double*,int>> arrays) {
        if (pending_active) die(EXIT_INTERNAL_ASSERT, "overlapping halo buffers reused before completion");
        const bool nodes = entity == MpiEntity::Node;
        Fields fields;
        for (auto array : arrays) fields.add(array.first, nodes ? local.nnode : local.nelem, array.second);
        #pragma acc wait
        double start = MPI_Wtime();
        for (int i=0; i<int(part.halos.size()); ++i) {
            const auto& halo = part.halos[i];
            const auto& receiving = nodes ? halo.recv_nodes : halo.recv_elements;
#ifdef ACC
            const auto& sending = nodes ? halo.send_nodes : halo.send_elements;
            const int send_count = count(std::size_t(sending.size()) * fields.columns);
            const int recv_count = count(std::size_t(receiving.size()) * fields.columns);
            peer[i].gpu_recv.resize(recv_count);
            pack_device(fields, sending, peer[i].gpu_send);
            transfers[i] = {halo.rank, peer[i].gpu_send.host(), peer[i].gpu_recv.host(), send_count, recv_count};
            payload_bytes += std::uint64_t(send_count) * sizeof(double);
#else
            peer[i].recv.resize(std::size_t(receiving.size()) * fields.columns);
            pack(fields, nodes ? halo.send_nodes : halo.send_elements, peer[i].send);
            transfers[i] = {halo.rank, peer[i].send.data(), peer[i].recv.data(),
                            count(peer[i].send.size()), count(peer[i].recv.size())};
            payload_bytes += peer[i].send.size() * sizeof(double);
#endif
        }
        pending = mpi_post_buffers(transfers.data(), transfers.size(), nodes ? 30 : 31, requests.data());
        pending_nodes = nodes;
        pending_fields = fields;
        pending_active = true;
        posted_at = MPI_Wtime();
        communication_s += MPI_Wtime() - start;
    }

    void finish_exchange() {
        if (!pending_active) die(EXIT_INTERNAL_ASSERT, "halo completion without a posted exchange");
        const double started = MPI_Wtime();
        interior_window_s += started - posted_at;
        mpi_finish_buffers(pending, requests.data());
        wait_s += MPI_Wtime() - started;
        for (int i=0; i<int(part.halos.size()); ++i) {
#ifdef ACC
            unpack_device(pending_fields, pending_nodes ? part.halos[i].recv_nodes : part.halos[i].recv_elements,
                          peer[i].gpu_recv);
#else
            const auto& rows = pending_nodes ? part.halos[i].recv_nodes : part.halos[i].recv_elements;
            unpack(pending_fields, rows.data(), rows.size(), peer[i].recv.data());
#endif
        }
        communication_s += MPI_Wtime() - started;
        pending_active = false;
    }

    void exchange(MpiEntity entity, std::initializer_list<std::pair<double*,int>> arrays) {
        begin_exchange(entity, arrays);
        finish_exchange();
    }
};

MpiSimulation::MpiSimulation(Param& param) : impl_(new Impl(param)) {
    int provided;
    MPI_Init_thread(nullptr, nullptr, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED) die(EXIT_RUNTIME_RESOURCE, "MPI_THREAD_FUNNELED required");
    MPI_Comm_rank(MPI_COMM_WORLD, &impl_->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &impl_->size);
#ifdef ACC
    MPI_Comm local;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, impl_->rank, MPI_INFO_NULL, &local);
    int local_rank, local_size;
    MPI_Comm_rank(local, &local_rank);
    MPI_Comm_size(local, &local_size);
    if (selected_device != local_rank) die(EXIT_CONFIG_VALUE, "GPU device selection differs from MPI's local rank");
    int current = -1;
    mpi_cuda_check(cudaGetDevice(&current));
    if (current != selected_device) die(EXIT_INTERNAL_ASSERT, "CUDA/OpenACC device selection differs");
    cudaDeviceProp device{};
    mpi_cuda_check(cudaGetDeviceProperties(&device, current));
    std::vector<char> uuids(16*local_size);
    MPI_Allgather(device.uuid.bytes, 16, MPI_BYTE, uuids.data(), 16, MPI_BYTE, local);
    for (int i=0; i<local_size; ++i) {
        for (int j=i+1; j<local_size; ++j) {
            if (std::memcmp(uuids.data() + 16*i, uuids.data() + 16*j, 16) == 0)
                die(EXIT_CONFIG_VALUE, "GPU MPI ranks selected the same physical GPU");
        }
    }
    MPI_Comm_free(&local);
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int length;
    MPI_Get_processor_name(hostname, &length);
    const char* hex = "0123456789abcdef";
    char uuid[33] = {};
    for (int i=0; i<16; ++i) {
        unsigned char byte = static_cast<unsigned char>(device.uuid.bytes[i]);
        uuid[2*i] = hex[byte >> 4];
        uuid[2*i+1] = hex[byte & 15];
    }
    std::ostringstream identity;
    identity << "[mpi][device] rank=" << impl_->rank << " host=" << hostname
             << " local_rank=" << local_rank << " device=" << current << " uuid=" << uuid
             << " fields=managed halo_transport=host_staged\n";
    std::cout << identity.str() << std::flush;
#endif
    if (is_root()) {
        const char* option = std::getenv("DES_MPI_OVERLAP");
        if (option && std::string(option) != "0" && std::string(option) != "1")
            die(EXIT_CONFIG_VALUE, "DES_MPI_OVERLAP must be 0 or 1");
        impl_->overlap_enabled = !option || std::string(option) == "1";
    }
    MPI_Bcast(&impl_->overlap_enabled, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    impl_->local.mpi = this;
    if (is_root())
        std::cout << "[mpi] compact spatial partitions; ranks=" << impl_->size
                  << "; neighbor halos; centralized mesh management\n";
}

MpiSimulation::~MpiSimulation() {
    if (is_root()) {
        Command stop{};
        stop.operation = int(MpiOp::Stop);
        MPI_Bcast(&stop, sizeof(stop), MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    MPI_Finalize();
}

bool MpiSimulation::is_root() const {
    return impl_->rank == 0;
}

bool MpiSimulation::distributed() const {
    if (!impl_->active) return false;
#ifdef ACC
    return true; // A single GPU exercises the same compact backend without halos.
#else
    return impl_->size > 1;
#endif
}

bool MpiSimulation::in_kernel() const {
    return impl_->inside && !management(impl_->current);
}

int MpiSimulation::owned_nodes() const {
    return impl_->part.owned_nodes;
}

int MpiSimulation::owned_elements() const {
    return impl_->part.owned_elements;
}

int MpiSimulation::global_nodes() const {
    return impl_->global_n;
}

double MpiSimulation::result() const {
    return impl_->scalar;
}

PTResult MpiSimulation::pt_result() const {
    return impl_->pt;
}

void MpiSimulation::attach(Variables& var) {
    impl_->master = &var;
    impl_->active = true;
    synchronize(var);
}

void MpiSimulation::synchronize(Variables& var) {
    dispatch(MpiOp::State, var);
}

bool MpiSimulation::dispatch(MpiOp op, const Variables& var, double arg) {
    auto& impl = *impl_;
    if (!distributed() || impl.inside || impl.master != &var) return false;
    // Keep marker maintenance and requested surface reports on root. Other
    // erosion-free steps use the existing compact mechanical correction.
    if (op == MpiOp::Surface) {
        const bool surface_report = NDIMS == 2 && impl.param.control.is_reporting_terrigenous_info &&
            var.steps != 0 && var.steps % 10000 == 0;
        if (var.steps % impl.param.mesh.quality_check_step_interval != 0 && !surface_report) {
            bool idle = impl.param.control.surface_process_option == 0;
            if (idle) {
                std::fill(impl.master->surfinfo.dh->begin(), impl.master->surfinfo.dh->end(), 0.0);
                impl.master->surfinfo.max_surf_vel = 0;
                op = MpiOp::SurfaceIdle;
            }
        }
    }
    Command command{};
    command.operation = int(op);
    command.steps = var.steps;
    command.time = var.time;
    command.dt = var.dt;
    command.max_vbc = var.max_vbc_val;
    command.max_vel = var.max_global_vel_mag;
    command.global_dt = var.global_dt_min;
    command.surface_vel = var.surfinfo.max_surf_vel;
    command.argument = arg;
    command.hydraulic = impl.param.control.has_hydraulic_diffusion;
    command.body_adjustment = impl.param.ic.has_body_force_adjustment;
    MPI_Bcast(&command, sizeof(command), MPI_BYTE, 0, MPI_COMM_WORLD);
    impl.command(command);
    impl.execute(op, arg);
    return true;
}

int MpiSimulation::worker_loop() try {
    for (;;) {
        Command command{};
        MPI_Bcast(&command, sizeof(command), MPI_BYTE, 0, MPI_COMM_WORLD);
        if (MpiOp(command.operation) == MpiOp::Stop) return 0;
        impl_->command(command);
        impl_->execute(MpiOp(command.operation), command.argument);
    }
}
catch (const std::bad_alloc&) {
    die(EXIT_RUNTIME_RESOURCE, "MPI worker allocation failed");
}
catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    die(EXIT_INTERNAL_UNREACHABLE, "MPI worker exception");
}

void MpiSimulation::exchange(MpiEntity entity, std::initializer_list<std::pair<double*,int>> arrays) {
    impl_->exchange(entity, arrays);
}

void MpiSimulation::begin_exchange(MpiEntity entity, std::initializer_list<std::pair<double*,int>> arrays) {
    impl_->begin_exchange(entity, arrays);
    if (!impl_->overlap_enabled) impl_->finish_exchange();
    else ++impl_->overlap_calls;
}

void MpiSimulation::finish_exchange() {
    impl_->finish_exchange();
}

bool MpiSimulation::overlap() const {
    return impl_->overlap_enabled;
}

MpiNodeRows MpiSimulation::gather_rows(int pass) const {
    const auto& rows = pass == 0 ? impl_->part.interior_nodes : impl_->part.interface_nodes;
    return {rows.data(), int(rows.size())};
}

void MpiSimulation::collect(const Variables& var, bool averages) {
    dispatch(MpiOp::Collect, var, averages ? 1 : 0);
}

void MpiSimulation::geometry(const Variables& var) {
    dispatch(MpiOp::Geometry, var);
}

void MpiSimulation::boundary(const Variables& var) {
    dispatch(MpiOp::Boundary, var);
}

void MpiSimulation::sample(const Variables& var, const int_vec& nodes, const int_vec& elements) {
    impl_->queries_n = nodes;
    impl_->queries_e = elements;
    dispatch(MpiOp::Sample, var);
}

bool MpiSimulation::metrics(const Variables& var, double& speed, double& plastic_rate) {
    if (!dispatch(MpiOp::Metrics, var)) return false;
    speed = impl_->measurements[0];
    plastic_rate = impl_->measurements[1];
    return true;
}

void MpiSimulation::minimum(double* data, int count) {
    double start = MPI_Wtime();
    MPI_Allreduce(MPI_IN_PLACE, data, count, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    impl_->communication_s += MPI_Wtime() - start;
}

double MpiSimulation::maximum(double value) {
    maximum(&value, 1);
    return value;
}

double MpiSimulation::sum(double value) {
    sum(&value, 1);
    return value;
}

void MpiSimulation::sum(double* data, int count) {
    double start = MPI_Wtime();
    MPI_Allreduce(MPI_IN_PLACE, data, count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    impl_->communication_s += MPI_Wtime() - start;
}

void MpiSimulation::maximum(double* data, int count) {
    double start = MPI_Wtime();
    MPI_Allreduce(MPI_IN_PLACE, data, count, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    impl_->communication_s += MPI_Wtime() - start;
}

int MpiSimulation::sum(int value) {
    int total;
    double start = MPI_Wtime();
    MPI_Allreduce(&value, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    impl_->communication_s += MPI_Wtime() - start;
    return total;
}

void MpiSimulation::maximum(int* data, int count) {
    double start = MPI_Wtime();
    MPI_Allreduce(MPI_IN_PLACE, data, count, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    impl_->communication_s += MPI_Wtime() - start;
}

void MpiSimulation::report() {
    if (impl_->master) dispatch(MpiOp::Report, *impl_->master);
    const auto& impl = *impl_;
    std::cout << std::fixed << std::setprecision(9);
    for (int i=0; i<int(MpiOp::Count); ++i) {
        if (impl.calls[i])
            std::cout << "[mpi][kernel] name=" << op_name(MpiOp(i)) << " calls=" << impl.calls[i]
                      << " seconds=" << impl.elapsed[i] << '\n';
    }
    std::cout << "[mpi][communication] seconds=" << impl.communication_s << " halo_sent_bytes=" << impl.all_payload_bytes
              << " root_collection_bytes=" << impl.root_bytes << " state_bytes=" << impl.state_bytes << '\n';
    const int stride = 3*int(MpiOp::Count)+3;
    for (int r=0; r<impl.size && !impl.rank_profile.empty(); ++r) {
        const double* profile = impl.rank_profile.data() + r*stride;
        for (int i=0; i<int(MpiOp::Count); ++i) {
            if (profile[3*i])
                std::cout << "[mpi][rank-profile] rank=" << r << " operation=" << op_name(MpiOp(i))
                          << " calls=" << profile[3*i]
                          << " elapsed_s=" << profile[3*i+1] << " communication_s=" << profile[3*i+2]
                          << " outside_communication_s=" << std::max(0.0, profile[3*i+1] - profile[3*i+2]) << '\n';
        }
        std::cout << "[mpi][overlap] rank=" << r << " calls=" << profile[stride-1] << " wait_s=" << profile[stride-3]
                  << " work_window_s=" << profile[stride-2] << '\n';
    }
}
#endif
