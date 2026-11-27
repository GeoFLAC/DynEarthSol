#if !defined(DYNEARTHSOL3D_ARRAY2D_h)
#define DYNEARTHSOL3D_ARRAY2D_h

#include "mpi_acc.hpp"

#include <algorithm>
#include <cstddef>
#include <cstring>  // for memcpy
#include <type_traits>
#include <vector>

// n rows of N components. With SOA, component d of row i is at a_[d*n_ + i]; without, at
// a_[i*N + d]. The accessors hide the layout; data() does not, and the MPI halo packing reads
// it as SoA columns.
template <typename T, int N>
class Array2D {

#ifdef SOA
    static constexpr bool soa = true;
#else
    static constexpr bool soa = false;
#endif

    T* a_ = nullptr;
    int n_ = 0;

    // A view of row i: where its first component is, and how far apart its components are.
    T* row_start(std::size_t i) const { return soa ? a_ + i : a_ + i * N; }
    int row_stride() const { return soa ? n_ : 1; }

public:
    //
    // Accessors. A row accessor is a view of one row. The two are separate structs, not a
    // shared base: OpenACC device code built by nvc++ reads out of bounds through the
    // inherited members.
    //
    // Assigning to an Accessor copies the values into the row it views.
    struct Accessor {
        T* ptr_;      // the row's first component
        int stride_;  // distance between its components

        Accessor(T* p, int s) : ptr_(p), stride_(s) {}
        Accessor(const Accessor&) = default;  // copies the view; operator= copies the values

        T& operator[](int dim) const { return ptr_[dim * stride_]; }

        Accessor& operator=(const T& val) {
            for (int d = 0; d < N; ++d) (*this)[d] = val;
            return *this;
        }

        Accessor& operator=(const Accessor& other) {
            if (this != &other) copy_from(other);
            return *this;
        }

        void copy_to(T* dest) const {
            for (int d = 0; d < N; ++d) dest[d] = (*this)[d];
        }

        void copy_from(const T* src) {
            for (int d = 0; d < N; ++d) (*this)[d] = src[d];
        }

        template <typename OtherAccessor>
        void copy_from(const OtherAccessor& other) {
            for (int d = 0; d < N; ++d) (*this)[d] = other[d];
        }
    };

    // Assigning to a ConstAccessor points it at another row.
    struct ConstAccessor {
        const T* ptr_;
        int stride_;

        ConstAccessor(const T* p, int s) : ptr_(p), stride_(s) {}
        ConstAccessor(const Accessor& other) : ptr_(other.ptr_), stride_(other.stride_) {}

        const T& operator[](int dim) const { return ptr_[dim * stride_]; }

        void copy_to(T* dest) const {
            for (int d = 0; d < N; ++d) dest[d] = (*this)[d];
        }
    };

    // One component of every row.
    template <typename U>
    struct Component {
        U*  ptr_;    // component d of row 0
        int stride_; // distance between rows
        int n_;

        U& operator[](int i) const { return ptr_[i * stride_]; }
        int size() const { return n_; }
    };
    using ComponentAccessor = Component<T>;
    using ConstComponentAccessor = Component<const T>;

    // The rows named by an index list, such as the nodes of an element.
    struct ConstIndirectAccessor {
        const Array2D* source_array;
        const int* index_ptr;
        int index_stride;

        ConstAccessor operator[](int k) const {
            return (*source_array)[index_ptr[k * index_stride]];
        }
    };

    // indices is a pointer to row numbers, or an accessor of them (a row of conn_t).
    template <typename Indices>
    ConstIndirectAccessor view_const(Indices indices) const {
        if constexpr (std::is_pointer_v<Indices>)
            return ConstIndirectAccessor{ this, indices, 1 };
        else
            return ConstIndirectAccessor{ this, indices.ptr_, indices.stride_ };
    }

    //
    // constructors & destructor
    //
    Array2D() = default;

    explicit Array2D(int size) : a_(size > 0 ? new T[N * size] : nullptr), n_(size) {}

    Array2D(int size, const T& val) : Array2D(size) {
        if (a_) std::fill_n(a_, N * n_, val);
    }

    // From row-major (AoS) data; null aos_data gives zeros.
    Array2D(const T* aos_data, int size) : Array2D(size) {
        if (!a_) return;
        if (aos_data == nullptr) {
            std::fill_n(a_, N * n_, T(0));
        } else if constexpr (soa) {
#ifndef ACC
            #pragma omp parallel for default(none) shared(aos_data) \
                collapse(2) if(n_ > 10000)
#endif
            #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector collapse(2)
            for (int i = 0; i < n_; ++i) {
                for (int d = 0; d < N; ++d) {
                    a_[d * n_ + i] = aos_data[i * N + d];
                }
            }
        } else {
            std::memcpy(a_, aos_data, sizeof(T) * N * n_);
        }
    }

    Array2D(const Array2D& src) : Array2D(int(src.size())) {
        if (a_) std::memcpy(a_, src.data(), sizeof(T) * N * n_);
    }

    // A move takes the buffer and leaves the source empty.
    Array2D(Array2D&& src) noexcept : a_(src.a_), n_(src.n_) {
        src.a_ = nullptr;
        src.n_ = 0;
    }

    Array2D& operator=(Array2D&& rhs) noexcept {
        if (this != &rhs) {
            delete [] a_;
            a_ = rhs.a_;
            n_ = rhs.n_;
            rhs.a_ = nullptr;
            rhs.n_ = 0;
        }
        return *this;
    }

    // A copy is explicit: construct one. Assignment would hide an O(n) copy.
    Array2D& operator=(const Array2D&) = delete;

    ~Array2D() { delete [] a_; }

    //
    // methods
    //
    [[nodiscard]] T* data() noexcept { return a_; }
    [[nodiscard]] const T* data() const noexcept { return a_; }
    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] int num_elements() const noexcept { return N * n_; }

    void resize(int size, bool preserve_data = true) {
        if (size == n_) return;

        T* new_a = size > 0 ? new T[N * size] : nullptr;

        if (preserve_data && a_ != nullptr && size > 0) {
            const int copy_count = std::min(size, n_);
            if constexpr (soa) {
                // each component separately
                for (int d = 0; d < N; ++d)
                    std::memcpy(new_a + d * size, a_ + d * n_, sizeof(T) * copy_count);
            } else {
                std::memcpy(new_a, a_, sizeof(T) * N * copy_count);
            }
        }

        delete [] a_;
        a_ = new_a;
        n_ = size;
    }

    // Resize, and set every value to val. The new buffer comes first, so a failed allocation
    // leaves the array as it was.
    void resize(int size, const T& val) {
        if (size != n_) {
            T* new_a = size > 0 ? new T[N * size] : nullptr;
            delete [] a_;
            a_ = new_a;
            n_ = size;
        }
        if (a_) std::fill_n(a_, N * n_, val);
    }

    void fill_component(int d, const T& val = T(0)) {
        if (!a_ || n_ == 0 || d < 0 || d >= N) return;

        T* ptr = a_;
        // not const: g++ 8 refuses a const variable in shared(), later g++ require it there
        int stride = soa ? 1 : N;
        int offset = soa ? d * n_ : d;

#ifndef ACC
        #pragma omp parallel for default(none) shared(ptr, offset, stride, val) \
            if(n_ > 10000)
#endif
        #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector
        for (int i = 0; i < n_; ++i) {
            ptr[offset + i * stride] = val;
        }
    }

    void zero_component(int d) {
        fill_component(d, T(0));
    }

    //
    // I/O and transition helpers, in row-major (AoS) order
    //
    void load_from_buffer(const T* buffer, std::size_t count, bool should_strip = true) {
        if (a_ == nullptr) {
            a_ = new T[N * count];
            n_ = count;
        } else if (count > std::size_t(n_) || should_strip) {
            this->resize(count, false);
            n_ = count;
        }
#ifndef ACC
        #pragma omp parallel for default(none) shared(buffer, count) \
            collapse(2) if(count > 10000)
#endif
        #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector collapse(2)
        for (std::size_t i = 0; i < count; ++i) {
            for (int d = 0; d < N; ++d) {
                (*this)[i][d] = buffer[i * N + d];
            }
        }
    }

    void pack_to(std::vector<T>& buffer, std::size_t limit_size = 0) const {
        std::size_t count = (limit_size > 0 && limit_size <= std::size_t(n_)) ? limit_size : n_;
        std::size_t total_elements = count * N;

        buffer.resize(total_elements);

#ifndef ACC
        #pragma omp parallel for default(none) shared(buffer, count) \
            collapse(2) if(count > 10000)
#endif
        #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector collapse(2)
        for (std::size_t i = 0; i < count; ++i) {
            for (int d = 0; d < N; ++d) {
                buffer[i * N + d] = (*this)[i][d];
            }
        }
    }

#ifdef ACC
    void pack_to_xyz_float(std::vector<float3>& buffer, std::size_t limit_size = 0) const {
        std::size_t count = (limit_size > 0 && limit_size <= std::size_t(n_)) ? limit_size : n_;

        if (buffer.size() < count)
            buffer.resize(count);

        #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector
        for (std::size_t i = 0; i < count; ++i) {
            buffer[i].x = (float)(*this)[i][0];
            buffer[i].y = (float)(*this)[i][1];
#ifdef THREED
            buffer[i].z = (float)(*this)[i][2];
#else
            buffer[i].z = 0.0;
#endif
        }
    }
#endif

    //
    // index accessing
    //
    Accessor operator[](std::size_t i) { return Accessor{ row_start(i), row_stride() }; }
    ConstAccessor operator[](std::size_t i) const { return ConstAccessor{ row_start(i), row_stride() }; }

    ComponentAccessor component(int d) {
        return soa ? ComponentAccessor{ a_ + d * n_, 1, n_ } : ComponentAccessor{ a_ + d, N, n_ };
    }

    ConstComponentAccessor component_const(int d) const {
        return soa ? ConstComponentAccessor{ a_ + d * n_, 1, n_ } : ConstComponentAccessor{ a_ + d, N, n_ };
    }

    //
    // iterators, over the whole buffer in storage order
    //
    using iterator = T*;
    using const_iterator = const T*;
    iterator begin() { return a_; }
    const_iterator begin() const { return a_; }
    iterator end() { return a_ + N * n_; }
    const_iterator end() const { return a_ + N * n_; }
};

#endif
