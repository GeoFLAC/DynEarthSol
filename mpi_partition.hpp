#ifndef DES_MPI_PARTITION_HPP
#define DES_MPI_PARTITION_HPP
#ifdef USEMPI
#include <vector>
#include <unordered_map>
#include "parameters.hpp"
#include "mpi_state.hpp"

struct MpiHalo {
    int rank = 0;
    int_vec send_nodes, recv_nodes, send_elements, recv_elements;
    void archive(MpiState& a) {
        a.value(rank); a.value(send_nodes); a.value(recv_nodes);
        a.value(send_elements); a.value(recv_elements);
    }
};

// Owned rows precede ghost rows; each group follows global ID ordering.
// Node supports are rebuilt separately in the original global element order.
struct MpiPartition {
    int owned_nodes = 0, owned_elements = 0, owned_markers = 0;
    int_vec nodes, elements, node_rows, element_rows;
    int_vec surface_nodes, surface_elements; // all local rows touched by root surface work
    int_vec surface_owned_nodes, surface_owned_elements;
    std::vector<MpiHalo> halos;
    std::unordered_map<int,int> node_local, element_local;
    void index();
    void archive(MpiState&);
};

std::vector<MpiPartition> mpi_partition_mesh(const Variables&, int ranks);
#endif
#endif
