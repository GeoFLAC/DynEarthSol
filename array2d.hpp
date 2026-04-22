#if !defined(DYNEARTHSOL3D_ARRAY2D_h)
#define DYNEARTHSOL3D_ARRAY2D_h

#include <vector>
#include <algorithm>
#include <cstring>  // for memcpy

template <typename T, int N>
class Array2D {

    T* a_;
    int n_;

public:
    // 
    // Accessor
    // 
    struct Accessor {
        T* ptr_;      // pointer to the start of the array (e.g., &a_[0])
        int stride_;  // stride = n_ (total elements) for SoA layout

        T& operator[](const int dim) const {
            return ptr_[dim * stride_];
        }

        Accessor& operator=(const T& val) {
            for (int d = 0; d < N; ++d) {
                (*this)[d] = val;
            }
            return *this;
        }

        Accessor& operator=(const Accessor& other) {
            if (this == &other) return *this;
            for (int d = 0; d < N; ++d) {
                (*this)[d] = other[d];
            }
            return *this;
        }

        void copy_to(T* dest) const {
            for (int d = 0; d < N; ++d) {
                dest[d] = ptr_[d * stride_];
            }
        }

        void copy_from(const T* src) {
            for (int d = 0; d < N; ++d) {
                ptr_[d * stride_] = src[d];
            }
        }

        template <typename OtherAccessor>
        void copy_from(const OtherAccessor& other) {
            for (int d = 0; d < N; ++d) {
                (*this)[d] = other[d];
            }
        }        
    };

    struct ConstAccessor {
        const T* ptr_;
        int stride_;

        ConstAccessor(const T* p, int s) : ptr_(p), stride_(s) {}

        ConstAccessor(const Accessor& other) 
            : ptr_(other.ptr_), stride_(other.stride_) 
        {}

        struct ConstStridedLocalIndirect {
            const T* ptr_;
            int stride_;
            const int* indices_;
            const T& operator[](int k) const {
                return ptr_[indices_[k] * stride_];
            }
        };
        
        ConstStridedLocalIndirect subset(const int* indices) const {
            return ConstStridedLocalIndirect{ptr_, stride_, indices};
        }

        const T& operator[](int dim) const {
            return ptr_[dim * stride_];
        }

        void copy_to(T* dest) const {
            for (int d = 0; d < N; ++d) {
                dest[d] = ptr_[d * stride_];
            }
        }
    };

    //
    // View
    //
    struct ConstIndirectAccessor {
        const Array2D<T, N>* source_array;
        const int* index_ptr;
        int index_stride;

        ConstAccessor operator[](int k) const {
            int real_idx = index_ptr[k * index_stride];
            return (*source_array)[real_idx];
        }
    };

    struct IndirectAccessor {
        Array2D<T, N>* source_array;
        const int* index_ptr;
        int index_stride;

        Accessor operator[](int k) const {
            int real_idx = index_ptr[k * index_stride];
            return (*source_array)[real_idx];
        }
    };

    IndirectAccessor view(int* indices) {
        return IndirectAccessor{ this, indices, 1 };
    }

    IndirectAccessor view(const int* indices) {
        return IndirectAccessor{ this, indices, 1 };
    }

    template <typename IntAccessor>
    IndirectAccessor view(IntAccessor indices_acc) {
        return IndirectAccessor{ this, indices_acc.ptr_, indices_acc.stride_ };
    }

    ConstIndirectAccessor view(const int* indices) const {
        return ConstIndirectAccessor{ this, indices, 1 };
    }
    
    template <typename IntAccessor>
    ConstIndirectAccessor view(IntAccessor indices_acc) const {
        return ConstIndirectAccessor{ this, indices_acc.ptr_, indices_acc.stride_ };
    }

    ConstIndirectAccessor view_const(int* indices) const {
        return ConstIndirectAccessor{ this, indices, 1 };
    }

    ConstIndirectAccessor view_const(const int* indices) const {
        return ConstIndirectAccessor{ this, indices, 1 };
    }

    template <typename IntAccessor>
    ConstIndirectAccessor view_const(IntAccessor indices_acc) const {
        return ConstIndirectAccessor{ this, indices_acc.ptr_, indices_acc.stride_ };
    }

    //
    // I/O and Transition Helpers
    //
    void load_from_buffer(const T* buffer, std::size_t count) {
        if (a_ == nullptr)
            a_ = new T[N * count];
        else
            this->resize(count, false);
        n_ = count;
#ifndef ACC
        #pragma omp parallel for collapse(2) if(count > 10000)
#endif
        #pragma acc parallel loop gang vector collapse(2)
        for (std::size_t i = 0; i < count; ++i) {
            for (int d = 0; d < N; ++d) {
                (*this)[i][d] = buffer[i * N + d]; 
            }
        }
    }

//     void copy_all_to(T* buffer) const {
//         if (!a_) return;
// #ifndef ACC
//         #pragma omp parallel for collapse(2) if(n_ > 10000)
// #endif
//         #pragma acc parallel loop gang vector collapse(2)
//         for(int i=0; i<n_; ++i) {
//             for(int d=0; d<N; ++d)
//                 buffer[i*N + d] = (*this)[i][d];
//         }
//     }

    void pack_to(std::vector<T>& buffer, std::size_t limit_size = 0) const {
        std::size_t count = (limit_size > 0 && limit_size <= n_) ? limit_size : n_;
        std::size_t total_elements = count * N;

        buffer.resize(total_elements);

#ifndef ACC
        #pragma omp parallel for collapse(2) if(n_ > 10000)
#endif
        #pragma acc parallel loop gang vector collapse(2)
        for (std::size_t i = 0; i < count; ++i) {
            for (int d = 0; d < N; ++d) {
                buffer[i * N + d] = (*this)[i][d];
            }
        }
    }

#ifdef ACC
    void pack_to_xyz_float(std::vector<float3>& buffer, std::size_t limit_size = 0) const {
        std::size_t count = (limit_size > 0 && limit_size <= n_) ? limit_size : n_;

        if (buffer.size() < count)
            buffer.resize(count);

        #pragma acc parallel loop gang vector
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
    // constructors & destructor
    //
    Array2D() : a_(nullptr), n_(0) {}
    
    explicit Array2D(int size) {
        n_ = size;
        if (n_ > 0) a_ = new T[N * n_];
        else a_ = nullptr;
    }

    Array2D(int size, const T& val) {
        n_ = size;
        if (n_ > 0) {
            a_ = new T[N * n_];
            std::fill_n(a_, N * n_, val);
        } else {
            a_ = nullptr;
        }
    }

    // AoS to SoA constructor
    Array2D(const T* aos_data, int size) {
        n_ = size;
        if (n_ > 0) {
            a_ = new T[N * n_];
            if (aos_data != nullptr) {
#ifdef SOA
#ifndef ACC
                #pragma omp parallel for collapse(2) if(n_ > 10000)
#endif
                #pragma acc parallel loop gang vector collapse(2)
                for (int i = 0; i < n_; ++i) {
                    for (int d = 0; d < N; ++d) {
                        // AoS(i, d) -> SoA(i, d)
                        a_[d * n_ + i] = aos_data[i * N + d];
                    }
                }
#else
                // If data is already in AoS layout, we can copy it directly
                std::memcpy(a_, aos_data, sizeof(T) * N * n_);
#endif
            } else {
                std::fill_n(a_, N * n_, T(0));
            }
        } else {
            a_ = nullptr;
        }
    }

    Array2D(const Array2D& src) {
        n_ = src.size();
        if (n_ > 0) {
            a_ = new T[N * n_];
            std::memcpy(a_, src.data(), sizeof(T) * N * n_);
        } else {
            a_ = nullptr;
        }
    }

    ~Array2D() { if (a_) delete [] a_; }

    //
    // methods
    //
    T* data() {return a_;}
    const T* data() const {return a_;}
    std::size_t size() const {return n_;}
    int num_elements() const {return N*n_;}

    void resize(int size, bool preserve_data = true) {
        if (size == n_) return;

        T* new_a = nullptr;
        if (size > 0) new_a = new T[N * size];

        if (preserve_data && a_ != nullptr && size > 0) {
            int copy_count = (size < n_) ? size : n_; 
#ifdef SOA
            // SoA Resize: move each dimension separately
            for (int d = 0; d < N; ++d) {
                T* src = a_ + d * n_;
                T* dst = new_a + d * size;
                std::memcpy(dst, src, sizeof(T) * copy_count);
            }
#else
            // AoS Resize: move all data at once
            std::memcpy(new_a, a_, sizeof(T) * N * copy_count);
#endif
        }

        if (a_) delete [] a_;
        a_ = new_a;
        n_ = size;
    }

    void resize(int size, const T& val) {
        if (size != n_) {
            T* new_a = nullptr;
            if (size > 0) new_a = new T[N * size];

            if (a_) delete [] a_;
            a_ = new_a;
            n_ = size;
        }
        if (a_) std::fill_n(a_, N * n_, val);
    }

    void steal_ref(Array2D& other) {
        if (a_) delete [] a_;
        a_ = other.a_;
        n_ = other.n_;
        other.a_ = nullptr;
        other.n_ = 0;
    }

    void reset(T* a, int n) {
        // Warning: this will take ownership of the pointer a, and delete it when destructed or resized
        if (a_) delete [] a_;
        a_ = a;
        n_ = n;
    }

    void nullify() {
        a_ = nullptr;
        n_ = 0;
    }

    //
    // index accessing
    //
    Accessor operator[](std::size_t i) {
        // pass n_ as stride
#ifdef SOA
        return Accessor{ a_ + i, n_ };
#else
        return Accessor{ a_ + i*N, 1 };
#endif
    }

    ConstAccessor operator[](std::size_t i) const {
#ifdef SOA
        return ConstAccessor{ a_ + i, n_ };
#else
        return ConstAccessor{ a_ + i*N, 1 };
#endif
    }

    Accessor at(std::size_t i) {
#ifdef SOA
        return Accessor{ a_ + i, n_ };
#else
        return Accessor{ a_ + i*N, 1 };
#endif
    }

    ConstAccessor at(std::size_t i) const {
#ifdef SOA
        return ConstAccessor{ a_ + i, n_ };
#else
        return ConstAccessor{ a_ + i*N, 1 };
#endif
    }

    //
    // iterators
    //
    typedef T* iterator;
    typedef const T* const_iterator;
    iterator begin() {return a_;}
    const_iterator begin() const {return a_;}
    iterator end() {return a_ + N*n_;}
    const_iterator end() const {return a_ + N*n_;}

    typedef T element;

private:
    // disable assignment operator
    Array2D<T,N>& operator=(const Array2D<T,N>& rhs);
};

#endif
