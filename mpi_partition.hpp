#ifndef DES_MPI_PARTITION_HPP
#define DES_MPI_PARTITION_HPP
#ifdef USEMPI
#include <cstdint>
#include <cstring>
#include <limits>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>
#include "array2d.hpp"
#include "mpi_transport.hpp"
#include "parameters.hpp"
#include "utils.hpp"

// An in-memory state transfer between ranks of the same executable. No pointers
// or container internals are put into the message; marker records stay on root.
class MpiState {
public:
    MpiState(mpi_host_vector<unsigned char>& bytes, bool writing)
        : bytes_(bytes), writing_(writing), position_(0) { if (writing) bytes_.clear(); }
    bool reading() const { return !writing_; }
    void raw(void* data, std::size_t bytes) {
        if (!bytes) return;
        if (writing_) {
            const auto* first = static_cast<const unsigned char*>(data);
            bytes_.insert(bytes_.end(), first, first + bytes);
        } else {
            if (position_ > bytes_.size() || bytes > bytes_.size() - position_)
                die(EXIT_INTERNAL_ASSERT, "truncated MPI state");
            std::memcpy(data, bytes_.data() + position_, bytes);
            position_ += bytes;
        }
    }
    template<class T> typename std::enable_if<std::is_arithmetic<T>::value>::type value(T& item) {
        raw(&item, sizeof(T));
    }
    template<class T> void value(std::vector<T>& items) {
        std::uint64_t size = items.size(); value(size);
        if (reading()) items.resize(size);
        vector_data(items, typename std::is_arithmetic<T>::type());
    }
    template<class T, int N> void value(Array2D<T,N>& items) {
        int size = items.size(); value(size);
        if (size < 0) die(EXIT_INTERNAL_ASSERT, "negative MPI array size");
        if (reading()) items.resize(size, false);
        raw(items.data(), std::size_t(size) * N * sizeof(T));
    }
    template<class T> void value(T*& item) {
        bool present = item != nullptr; value(present);
        if (reading()) {
            if (!present) { delete item; item = nullptr; }
            else if (!item) item = new T;
        }
        if (present) value(*item);
    }
    template<class A, class B> void value(std::pair<A,B>& item) { value(item.first); value(item.second); }
    void finish() const {
        if (reading() && position_ != bytes_.size())
            die(EXIT_INTERNAL_ASSERT, "unused bytes in MPI state");
    }
private:
    template<class T> void vector_data(std::vector<T>& items, std::true_type) {
        raw(items.data(), items.size() * sizeof(T));
    }
    template<class T> void vector_data(std::vector<T>& items, std::false_type) {
        for (auto& item : items) value(item);
    }
    mpi_host_vector<unsigned char>& bytes_;
    bool writing_;
    std::size_t position_;
};

struct MpiHalo {
    int rank = 0;
    int_vec send_nodes, recv_nodes, send_elements, recv_elements;
    void archive(MpiState& a) {
        a.value(rank); a.value(send_nodes); a.value(recv_nodes);
        a.value(send_elements); a.value(recv_elements);
    }
};

// Owned rows precede ghost rows; each group follows global ID ordering.
// Node supports are rebuilt separately in global element order.
struct MpiPartition {
    int owned_nodes = 0, owned_elements = 0, owned_markers = 0;
    int_vec nodes, elements, node_rows, element_rows;
    int_vec interior_nodes, interface_nodes;
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
