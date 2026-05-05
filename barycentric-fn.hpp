#ifndef DYNEARTHSOL3D_BARYCENTRIC_FN_HPP
#define DYNEARTHSOL3D_BARYCENTRIC_FN_HPP

#include "array2d.hpp"
#include "constants.hpp"
#include "parameters.hpp"

class Barycentric_transformation {

    /* Performing barycentric transformation to a point.
     *
     * The derivation of the formula can be found in
     * http://en.wikipedia.org/wiki/Barycentric_coordinate_system_(mathematics)
     */
    typedef Array2D<double,NODES_PER_ELEM*NDIMS> coeff_t;
    typedef coeff_t::Accessor CoeffAccessor;
    typedef coeff_t::ConstAccessor ConstCoeffAccessor;
    coeff_t coeff_;
    int nelem_;
    int elem_dim_;

public:

    Barycentric_transformation(const array_t &coord,
                               const conn_t &connectivity,
                               const double_vec &volume);
    Barycentric_transformation(const int_vec &elem,
                               const array_t &coord,
                               const conn_t &connectivity,
                               const double_vec &volume);
    Barycentric_transformation(const array_t &coord,
                               const conn_t&conn_surface,
                               const double_vec &volume,
                               const bool is_surface);
    Barycentric_transformation(ConstArrayIndirectAccessor coord,
                               const double volume);
    ~Barycentric_transformation();

    #pragma acc routine seq
    template <typename T>
    void transform(T point, int e, double *result) const;
    bool is_inside_elem(ConstArrayAccessor point, int elem) const;
    #pragma acc routine seq
    bool is_inside(const double *result) const;

private:

    inline int index1d(int node, int dim) const;
    inline int index2d(int node, int dim) const;
    inline int index3d(int node, int dim) const;

    void compute_coeff2d(ConstArrayAccessor a,
                         ConstArrayAccessor b,
                         ConstArrayAccessor c,
                         double area,
                         coeff_t::Accessor coeff_e);

#ifdef THREED
    void compute_coeff3d(ConstArrayAccessor a,
                         ConstArrayAccessor b,
                         ConstArrayAccessor c,
                         ConstArrayAccessor d,
                         double volume,
                         CoeffAccessor coeff_e);
#else
    void compute_coeff1d(ConstArrayAccessor a,
                         ConstArrayAccessor b,
                         double length,
                         CoeffAccessor coeff_e);
#endif

};

#endif
