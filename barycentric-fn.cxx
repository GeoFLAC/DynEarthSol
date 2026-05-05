#include "barycentric-fn.hpp"


Barycentric_transformation::Barycentric_transformation(const array_t &coord,
                                                       const conn_t &connectivity,
                                                       const double_vec &volume)
    : coeff_(connectivity.size()), nelem_(connectivity.size()), elem_dim_(NDIMS)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
#ifndef ACC
    #pragma omp parallel for default(none) \
        shared(coord, connectivity, volume)
#endif
    #pragma acc parallel loop gang vector async
    for (int e=0; e<nelem_; ++e) {
        int n0 = connectivity[e][0];
        int n1 = connectivity[e][1];
        int n2 = connectivity[e][2];

        ConstArrayAccessor a = coord[n0];
        ConstArrayAccessor b = coord[n1];
        ConstArrayAccessor c = coord[n2];

#ifdef THREED
        int n3 = connectivity[e][3];
        ConstArrayAccessor d = coord[n3];

        compute_coeff3d(a, b, c, d, volume[e], coeff_[e]);
#else
        compute_coeff2d(a, b, c, volume[e], coeff_[e]);
#endif
    }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

Barycentric_transformation::Barycentric_transformation(const int_vec &elem,
                                                       const array_t &coord,
                                                       const conn_t &connectivity,
                                                       const double_vec &volume)
    : coeff_(elem.size()), nelem_(elem.size()), elem_dim_(NDIMS)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
#ifndef ACC
    #pragma omp parallel for default(none) \
        shared(elem, coord, connectivity, volume)
#endif
    #pragma acc parallel loop gang vector async
    for (int i=0; i<nelem_; ++i) {
        int e = elem[i];
        int n0 = connectivity[e][0];
        int n1 = connectivity[e][1];
        int n2 = connectivity[e][2];

        ConstArrayAccessor a = coord[n0];
        ConstArrayAccessor b = coord[n1];
        ConstArrayAccessor c = coord[n2];

#ifdef THREED
        int n3 = connectivity[e][3];
        ConstArrayAccessor d = coord[n3];

        compute_coeff3d(a, b, c, d, volume[e], coeff_[i]);
#else
        compute_coeff2d(a, b, c, volume[e], coeff_[i]);
#endif
    }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}


Barycentric_transformation::Barycentric_transformation(const array_t &coord,
                                                       const conn_t &conn_surface,
                                                       const double_vec &area,
                                                       const bool is_surface)
    : coeff_(conn_surface.size()), nelem_(conn_surface.size()) ,elem_dim_(NDIMS-1)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
#ifndef ACC
    #pragma omp parallel for default(none) \
        shared(coord, conn_surface, area)
#endif
    #pragma acc parallel loop gang vector async
    for (int e=0; e<nelem_; ++e) {
        int n0 = conn_surface[e][0];
        int n1 = conn_surface[e][1];

        ConstArrayAccessor a = coord[n0];
        ConstArrayAccessor b = coord[n1];

#ifdef THREED
        int n2 = conn_surface[e][2];
        ConstArrayAccessor c = coord[n2];

        compute_coeff2d(a, b, c, area[e], coeff_[e]);
#else
        compute_coeff1d(a, b, area[e], coeff_[e]);
#endif
    }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

Barycentric_transformation::Barycentric_transformation(ConstArrayIndirectAccessor coord,
                                                       const double volume)
    : coeff_(1), nelem_(1), elem_dim_(NDIMS)
{
    ConstArrayAccessor a = coord[0];
    ConstArrayAccessor b = coord[1];
    ConstArrayAccessor c = coord[2];

#ifdef THREED
    ConstArrayAccessor d = coord[3];

    compute_coeff3d(a, b, c, d, volume, coeff_[0]);
#else
    compute_coeff2d(a, b, c, volume, coeff_[0]);
#endif
}


Barycentric_transformation::~Barycentric_transformation() {};

template <typename T>
void Barycentric_transformation::transform(T point, int e, double *result) const
{
    ConstCoeffAccessor cf = coeff_[e];
    if (elem_dim_ == 3) {
        for (int d=0; d<3; d++) {
            result[d] = cf[index3d(0,d)];
            for (int i=0; i<3; i++)
                result[d] += cf[index3d(i+1,d)]*point[i];
        }
    } else if (elem_dim_ == 2) {
        for (int d=0; d<2; d++) {
            result[d] = cf[index2d(0,d)];
            for (int i=0; i<2; i++)
                result[d] += cf[index2d(i+1,d)]*point[i];
        }
    } else if (elem_dim_ == 1) {
        for (int d=0; d<1; d++) {
            result[d] = cf[index1d(0,d)];
            for (int i=0; i<1; i++)
                result[d] += cf[index1d(i+1,d)]*point[i];
        }
    }
}

template
void Barycentric_transformation::transform<double*>(double* point, int e, double *result) const;
template
void Barycentric_transformation::transform<const double*>(const double* point, int e, double *result) const;
template
void Barycentric_transformation::transform<ArrayAccessor>(ArrayAccessor point, int e, double *result) const;
template
void Barycentric_transformation::transform<ConstArrayAccessor>(ConstArrayAccessor point, int e, double *result) const;


bool Barycentric_transformation::is_inside_elem(ConstArrayAccessor point, int elem) const
{
    double r[NDIMS];
    transform(point, elem, r);
    return is_inside(r);
}


bool Barycentric_transformation::is_inside(const double *r) const
{
    if (elem_dim_ == 3) {
        // 3D has larger round-off error in coeff_
        // => needs greater tolerance
        const double tolerance = 5e-11;

        if (r[0] >= -tolerance &&
            r[1] >= -tolerance &&
            r[2] >= -tolerance &&
            (r[0] + r[1] + r[2]) <= 1 + tolerance)
            return 1;
    } else if (elem_dim_ == 2) {
        const double tolerance = 1e-12;

        if (r[0] >= -tolerance &&
            r[1] >= -tolerance &&
            (r[0] + r[1]) <= 1 + tolerance)
            return 1;
    } else if (elem_dim_ == 1) {
        const double tolerance = 1e-12;

        if (r[0] >= -tolerance &&
            r[0] <= 1 + tolerance)
            return 1;
    }
    return 0;
}


inline int Barycentric_transformation::index1d(int node, int dim) const
{
    return node + dim;
}

inline int Barycentric_transformation::index2d(int node, int dim) const
{
    return node*2 + dim;
}
inline int Barycentric_transformation::index3d(int node, int dim) const
{
    return node*3 + dim;
}

void Barycentric_transformation::compute_coeff2d(ConstArrayAccessor a,
                                                 ConstArrayAccessor b,
                                                 ConstArrayAccessor c,
                                                 double area,
                                                 CoeffAccessor coeff_e)
{
    double det = 2 * area;

    coeff_e[index2d(0,0)] = (b[0]*c[1] - b[1]*c[0]) / det;
    coeff_e[index2d(0,1)] = (c[0]*a[1] - c[1]*a[0]) / det;
    coeff_e[index2d(1,0)] = (b[1] - c[1]) / det;
    coeff_e[index2d(1,1)] = (c[1] - a[1]) / det;
    coeff_e[index2d(2,0)] = (c[0] - b[0]) / det;
    coeff_e[index2d(2,1)] = (a[0] - c[0]) / det;
}

#ifdef THREED

void Barycentric_transformation::compute_coeff3d(ConstArrayAccessor a,
                                                 ConstArrayAccessor b,
                                                 ConstArrayAccessor c,
                                                 ConstArrayAccessor d,
                                                 double volume,
                                                 CoeffAccessor coeff_e)
{
    double det = 6 * volume;

    coeff_e[index3d(0,0)] = (b[0] * (c[1]*d[2] - d[1]*c[2]) +
                           c[0] * (d[1]*b[2] - b[1]*d[2]) +
                           d[0] * (b[1]*c[2] - c[1]*b[2])) / det;
    coeff_e[index3d(0,1)] = (a[0] * (d[1]*c[2] - c[1]*d[2]) +
                           c[0] * (a[1]*d[2] - d[1]*a[2]) +
                           d[0] * (c[1]*a[2] - a[1]*c[2])) / det;
    coeff_e[index3d(0,2)] = (a[0] * (b[1]*d[2] - d[1]*b[2]) +
                           b[0] * (d[1]*a[2] - a[1]*d[2]) +
                           d[0] * (a[1]*b[2] - b[1]*a[2])) / det;

    coeff_e[index3d(1,0)] = ((d[1] - b[1]) * (c[2] - b[2]) -
                           (c[1] - b[1]) * (d[2] - b[2])) / det;
    coeff_e[index3d(1,1)] = ((c[1] - a[1]) * (d[2] - c[2]) -
                           (d[1] - c[1]) * (c[2] - a[2])) / det;
    coeff_e[index3d(1,2)] = ((b[1] - d[1]) * (a[2] - d[2]) -
                           (a[1] - d[1]) * (b[2] - d[2])) / det;

    coeff_e[index3d(2,0)] = ((d[2] - b[2]) * (c[0] - b[0]) -
                           (c[2] - b[2]) * (d[0] - b[0])) / det;
    coeff_e[index3d(2,1)] = ((c[2] - a[2]) * (d[0] - c[0]) -
                           (d[2] - c[2]) * (c[0] - a[0])) / det;
    coeff_e[index3d(2,2)] = ((b[2] - d[2]) * (a[0] - d[0]) -
                           (a[2] - d[2]) * (b[0] - d[0])) / det;

    coeff_e[index3d(3,0)] = ((d[0] - b[0]) * (c[1] - b[1]) -
                           (c[0] - b[0]) * (d[1] - b[1])) / det;
    coeff_e[index3d(3,1)] = ((c[0] - a[0]) * (d[1] - c[1]) -
                           (d[0] - c[0]) * (c[1] - a[1])) / det;
    coeff_e[index3d(3,2)] = ((b[0] - d[0]) * (a[1] - d[1]) -
                           (a[0] - d[0]) * (b[1] - d[1])) / det;
}

#else

void Barycentric_transformation::compute_coeff1d(ConstArrayAccessor a,
                                                 ConstArrayAccessor b,
                                                 double length,
                                                 CoeffAccessor coeff_e)
{
    const double det = length;
    coeff_e[index1d(0,0)] = b[0] / det;
    coeff_e[index1d(0,1)] = -1.0 / det;
}

#endif


