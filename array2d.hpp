#if !defined(DYNEARTHSOL3D_ARRAY2D_h)
#define DYNEARTHSOL3D_ARRAY2D_h

#include <algorithm>
#include <cstring>  // for memcpy


template <typename T, int N>
class Array2D {

    T* a_;
    int n_;

public:
    //
    // constructors & destructor
    //
    Array2D() {a_ = NULL; n_ = 0;}
    Array2D(T* a, int n) {a_ = a; n_ = n;}
    Array2D(const Array2D& src) {
        n_ = src.size();
        a_ = new T[src.num_elements()];
        std::memcpy(a_, src.data(), sizeof(T)*src.num_elements());
    }

    Array2D(int size, const T& val) {
        a_ = new T[N*size];
        n_ = size;
        std::fill_n(a_, N*n_, val);
    }

    explicit
    Array2D(int size) {
        a_ = new T[N*size];
        n_ = size;
    }

    ~Array2D() {delete [] a_;}

    //
    // methods
    //
    T* data() {return a_;}
    const T* data() const {return a_;}
    std::size_t size() const {return n_;}
    int num_elements() const {return N*n_;}

    void resize(int n) {
        if (n <= n_) {
            // shrink
            n_ = n;
        }
        else {
            // expand
            T* tmp = new T[N*n];
            std::memcpy(tmp, a_, sizeof(T)*N*n_);
            delete [] a_;
            a_ = tmp;
            n_ = n;
        }
    }

    // steal the pointer from other, leave a NULL to other
    void steal_ref(Array2D& other) {
        delete [] a_;
        a_ = other.a_;
        n_ = other.n_;
        other.a_ = NULL;
        other.n_ = 0;
    }

    void reset(T* a, int n) {
        delete [] a_;
        a_ = a;
        n_ = n;
    }

    void nullify() {
        a_ = NULL;
        n_ = 0;
    }

    //
    // index accessing
    //
    T* operator[](std::size_t i) {return a_ + N*i;}
    const T* operator[](std::size_t i) const {return a_ + N*i;}

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
