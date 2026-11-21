#ifndef DES_MPI_STATE_HPP
#define DES_MPI_STATE_HPP
#ifdef USEMPI
#include <cstring>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>
#include "array2d.hpp"
#include "utils.hpp"

// An in-memory state transfer between ranks of the same executable. No pointers
// or container internals are put into the message; marker records stay on root.
class MpiState {
public:
    MpiState(std::vector<unsigned char>& bytes, bool writing)
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
    std::vector<unsigned char>& bytes_;
    bool writing_;
    std::size_t position_;
};
#endif
#endif
