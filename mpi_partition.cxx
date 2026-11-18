#ifdef USEMPI
#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include "mpi_partition.hpp"
#include "utils.hpp"
#include <metis.h>

void MpiPartition::index() {
    node_local.clear();
    element_local.clear();
    for (int i=0; i<int(nodes.size()); ++i) node_local.emplace(nodes[i], i);
    for (int i=0; i<int(elements.size()); ++i) element_local.emplace(elements[i], i);
    node_rows.resize(owned_nodes);
    element_rows.resize(owned_elements);
    std::iota(node_rows.begin(), node_rows.end(), 0);
    std::iota(element_rows.begin(), element_rows.end(), 0);
}

void MpiPartition::archive(MpiState& state) {
    state.value(owned_nodes);
    state.value(owned_elements);
    state.value(owned_markers);
    state.value(nodes);
    state.value(elements);
    state.value(surface_nodes);
    state.value(surface_elements);
    state.value(interior_nodes);
    state.value(interface_nodes);
    state.value(surface_owned_nodes);
    state.value(surface_owned_elements);
    int nhalos = halos.size();
    state.value(nhalos);
    if (state.reading()) halos.resize(nhalos);
    for (auto& halo : halos) halo.archive(state);
    if (state.reading()) index();
}

namespace {
void metis_partition(const Variables& var, int ranks, int_vec& owner) {
    // A face-only dual graph omits dependencies in DES node-star gathers.
    // Connect elements sharing any vertex. The volume objective is a proxy;
    // actual node/element halo bytes still depend on the node-owner assignment.
    std::vector<idx_t> offsets(var.nelem+1, 0);
    std::vector<idx_t> adjacency;
    std::vector<idx_t> element_part(var.nelem);
    int_vec neighbors;
    for (int e=0; e<var.nelem; ++e) {
        neighbors.clear();
        for (int j=0; j<NODES_PER_ELEM; ++j) {
            const int n = (*var.connectivity)[e][j];
            const int* patch = var.support.patch(n);
            neighbors.insert(neighbors.end(), patch, patch + var.support.size(n));
        }
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
        for (int other : neighbors) {
            if (other != e) adjacency.push_back(other);
        }
        if (adjacency.size() > std::size_t(IDX_MAX))
            die(EXIT_RUNTIME_RESOURCE, "METIS adjacency exceeds idx_t; use a 64-bit METIS build");
        offsets[e+1] = adjacency.size();
    }
    idx_t vertices = var.nelem;
    idx_t constraints = 1;
    idx_t parts = ranks;
    idx_t objective = 0;
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0;
    options[METIS_OPTION_SEED] = 0;
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
    options[METIS_OPTION_UFACTOR] = 30; // allow 3% element-count imbalance
    const int status = METIS_PartGraphKway(&vertices, &constraints, offsets.data(), adjacency.data(),
                                           nullptr, nullptr, nullptr, &parts, nullptr, nullptr,
                                           options, &objective, element_part.data());
    if (status != METIS_OK) die(EXIT_MESH_QUALITY, "METIS_PartGraphKway failed");
    for (int e=0; e<var.nelem; ++e) {
        if (element_part[e] < 0 || element_part[e] >= ranks) die(EXIT_INTERNAL_ASSERT, "invalid METIS owner");
        owner[e] = element_part[e];
    }
    std::cout << "[mpi][metis] objective=volume value=" << objective << " imbalance_tolerance=0.03\n";
}

void bisect(int_vec& ids, int begin, int end, int rank_begin, int rank_end,
            const std::vector<std::array<double,NDIMS>>& centers, int_vec& owner) {
    if (rank_end - rank_begin == 1) {
        for (int i=begin; i<end; ++i) owner[ids[i]] = rank_begin;
        return;
    }
    std::array<double,NDIMS> lower{};
    std::array<double,NDIMS> upper{};
    if (begin < end) {
        lower = upper = centers[ids[begin]];
    }
    for (int i=begin; i<end; ++i) {
        for (int d=0; d<NDIMS; ++d) {
            lower[d] = std::min(lower[d], centers[ids[i]][d]);
            upper[d] = std::max(upper[d], centers[ids[i]][d]);
        }
    }
    int axis = 0;
    for (int d=1; d<NDIMS; ++d) {
        if (upper[d] - lower[d] > upper[axis] - lower[axis]) axis = d;
    }
    std::sort(ids.begin() + begin, ids.begin() + end, [&](int left, int right) {
        return centers[left][axis] < centers[right][axis] ||
              (centers[left][axis] == centers[right][axis] && left < right);
    });
    const int rank_middle = (rank_begin + rank_end) / 2;
    const int middle = begin + std::int64_t(end - begin) * (rank_middle - rank_begin) / (rank_end - rank_begin);
    bisect(ids, begin, middle, rank_begin, rank_middle, centers, owner);
    bisect(ids, middle, end, rank_middle, rank_end, centers, owner);
}
}

std::vector<MpiPartition> mpi_partition_mesh(const Variables& var, int nranks) {
    const auto started = std::chrono::steady_clock::now();
    int_vec elem_owner(var.nelem);
    int_vec node_owner(var.nnode, nranks);
    // METIS partitions every mesh. Its coarsening is inappropriate for tiny
    // meshes with empty partitions, which use recursive coordinate bisection.
    // One rank (a single GPU) owns every element; METIS 5.1 divides by zero
    // when asked for one part.
    const bool geometric = var.nelem <= nranks || nranks == 1;
    if (geometric) {
        int_vec ids(var.nelem);
        std::iota(ids.begin(), ids.end(), 0);
        std::vector<std::array<double,NDIMS>> centers(var.nelem);
        for (int e=0; e<var.nelem; ++e) {
            for (int d=0; d<NDIMS; ++d) {
                double sum = 0;
                for (int j=0; j<NODES_PER_ELEM; ++j) sum += (*var.coord)[(*var.connectivity)[e][j]][d];
                centers[e][d] = sum / NODES_PER_ELEM;
            }
        }
        bisect(ids, 0, var.nelem, 0, nranks, centers, elem_owner);
    }
    else metis_partition(var, nranks, elem_owner);
    int_vec node_count(nranks, 0);
    std::vector<std::pair<int,int_vec>> shared;
    for (int n=0; n<var.nnode; ++n) {
        std::set<int> candidates;
        const int* patch = var.support.patch(n);
        for (int i=0; i<var.support.size(n); ++i) candidates.insert(elem_owner[patch[i]]);
        if (candidates.empty()) candidates.insert(0);
        if (candidates.size() == 1) {
            node_owner[n] = *candidates.begin();
            ++node_count[node_owner[n]];
        }
        else shared.emplace_back(n, int_vec(candidates.begin(), candidates.end()));
    }
    for (const auto& shared_node : shared) {
        int least_loaded = shared_node.second.front();
        for (int candidate : shared_node.second) {
            if (node_count[candidate] < node_count[least_loaded]) least_loaded = candidate;
        }
        node_owner[shared_node.first] = least_loaded;
        ++node_count[least_loaded];
    }
    std::vector<char> surface_element(var.nelem, 0);
    std::vector<char> surface_node(var.nnode, 0);
    for (int n : *var.surfinfo.top_nodes) {
        const int* patch = var.support.patch(n);
        for (int j=0; j<var.support.size(n); ++j) surface_element[patch[j]] = 1;
    }
    for (int e=0; e<var.nelem; ++e) {
        if (surface_element[e])
            for (int j=0; j<NODES_PER_ELEM; ++j) surface_node[(*var.connectivity)[e][j]] = 1;
    }
    for (int boundary=0; boundary<nbdrytypes; ++boundary) {
        for (int n : *var.bnodes[boundary]) surface_node[n] = 1;
    }
    std::vector<MpiPartition> parts(nranks);
    for (int rank=0; rank<nranks; ++rank) {
        auto& partition = parts[rank];
        std::vector<char> needed_element(var.nelem, 0);
        std::vector<char> needed_node(var.nnode, 0);
        for (int e=0; e<var.nelem; ++e) {
            if (elem_owner[e] == rank) {
                partition.elements.push_back(e);
                needed_element[e] = 1;
                if (var.markers_in_elem) partition.owned_markers += (*var.markers_in_elem)[e].size();
            }
        }
        partition.owned_elements = partition.elements.size();
        for (int n=0; n<var.nnode; ++n) {
            if (node_owner[n] == rank) {
                partition.nodes.push_back(n);
                needed_node[n] = 1;
                const int* patch = var.support.patch(n);
                for (int j=0; j<var.support.size(n); ++j) needed_element[patch[j]] = 1;
            }
        }
        partition.owned_nodes = partition.nodes.size();
        for (int e=0; e<var.nelem; ++e) {
            if (needed_element[e]) {
                if (elem_owner[e] != rank) partition.elements.push_back(e);
                for (int j=0; j<NODES_PER_ELEM; ++j) needed_node[(*var.connectivity)[e][j]] = 1;
            }
        }
        for (int n=0; n<var.nnode; ++n) {
            if (needed_node[n] && node_owner[n] != rank) partition.nodes.push_back(n);
        }
        partition.index();
        for (int i=0; i<partition.owned_nodes; ++i) {
            const int n = partition.nodes[i];
            const int* patch = var.support.patch(n);
            bool interior = true;
            for (int j=0; j<var.support.size(n); ++j) {
                if (elem_owner[patch[j]] != rank) {
                    interior = false;
                    break;
                }
            }
            (interior ? partition.interior_nodes : partition.interface_nodes).push_back(i);
        }
        for (int i=0; i<int(partition.nodes.size()); ++i) {
            if (surface_node[partition.nodes[i]]) {
                partition.surface_nodes.push_back(i);
                if (i < partition.owned_nodes) partition.surface_owned_nodes.push_back(i);
            }
        }
        for (int i=0; i<int(partition.elements.size()); ++i) {
            if (surface_element[partition.elements[i]]) {
                partition.surface_elements.push_back(i);
                if (i < partition.owned_elements) partition.surface_owned_elements.push_back(i);
            }
        }
    }
    std::vector<std::map<int,MpiHalo>> links(nranks);
    for (int rank=0; rank<nranks; ++rank) {
        auto& partition = parts[rank];
        for (int i=partition.owned_nodes; i<int(partition.nodes.size()); ++i) {
            int owner = node_owner[partition.nodes[i]];
            links[rank][owner].rank = owner;
            links[owner][rank].rank = rank;
            links[rank][owner].recv_nodes.push_back(i);
            links[owner][rank].send_nodes.push_back(parts[owner].node_local.at(partition.nodes[i]));
        }
        for (int i=partition.owned_elements; i<int(partition.elements.size()); ++i) {
            int owner = elem_owner[partition.elements[i]];
            links[rank][owner].rank = owner;
            links[owner][rank].rank = rank;
            links[rank][owner].recv_elements.push_back(i);
            links[owner][rank].send_elements.push_back(parts[owner].element_local.at(partition.elements[i]));
        }
    }
    for (int rank=0; rank<nranks; ++rank) {
        for (auto& link : links[rank]) parts[rank].halos.push_back(std::move(link.second));
    }
    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    std::cout << "[mpi][partitioner] method=" << (geometric ? "rcb" : "metis")
              << " seconds=" << elapsed << '\n';
    return parts;
}
#endif
