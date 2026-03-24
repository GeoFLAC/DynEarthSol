// Metal compute kernels for DynEarthSol
//
// Requires Apple GPU Family 8 (M2/A16 or later) or Mac2 GPU family
// for double-precision floating-point arithmetic.
//
// 2D kernels: NDIMS=2, NODES_PER_ELEM=3, NSTR=3
// 3D kernels: NDIMS=3, NODES_PER_ELEM=4, NSTR=6

#include <metal_stdlib>
using namespace metal;


// ============================================================
// 2D kernels  (triangle elements)
// ============================================================

// Compute triangle area (volume) for each 2D element.
// coord:   flat array, nnode * 2 doubles (x, z interleaved)
// conn:    flat array, nelem * 3 ints
// volume:  output, nelem doubles
kernel void compute_volume_2d(
    device const double *coord  [[buffer(0)]],
    device const int   *conn   [[buffer(1)]],
    device double      *volume [[buffer(2)]],
    uint gid                   [[thread_position_in_grid]])
{
    int n0 = conn[gid * 3];
    int n1 = conn[gid * 3 + 1];
    int n2 = conn[gid * 3 + 2];

    double ax = coord[n0 * 2],     az = coord[n0 * 2 + 1];
    double bx = coord[n1 * 2],     bz = coord[n1 * 2 + 1];
    double cx = coord[n2 * 2],     cz = coord[n2 * 2 + 1];

    double ab0 = bx - ax,  ab1 = bz - az;
    double ac0 = cx - ax,  ac1 = cz - az;

    // abs() matches the CPU triangle_area() in geometry.cxx which also uses
    // std::fabs().  The mesh guarantees consistent winding, but abs() makes the
    // volume computation robust to occasional numerical sign flips.
    volume[gid] = abs(ab0 * ac1 - ab1 * ac0) / 2.0;
}


// Compute strain rate for each 2D triangle element.
// vel:          flat array, nnode * 2 doubles (vx, vz interleaved)
// conn:         flat array, nelem * 3 ints
// shpdx/shpdz:  flat arrays, nelem * 3 doubles
// strain_rate:  output, nelem * 3 doubles (sxx, szz, sxz)
kernel void update_strain_rate_2d(
    device const double *vel         [[buffer(0)]],
    device const int   *conn        [[buffer(1)]],
    device const double *shpdx      [[buffer(2)]],
    device const double *shpdz      [[buffer(3)]],
    device double      *strain_rate [[buffer(4)]],
    uint gid                        [[thread_position_in_grid]])
{
    const int NPE = 3;
    int eb = gid * NPE;

    double vx[3], vz[3];
    for (int i = 0; i < NPE; i++) {
        int n = conn[eb + i];
        vx[i] = vel[n * 2];
        vz[i] = vel[n * 2 + 1];
    }

    // XX component
    double sxx = 0;
    for (int i = 0; i < NPE; i++)
        sxx += vx[i] * shpdx[eb + i];

    // ZZ component
    double szz = 0;
    for (int i = 0; i < NPE; i++)
        szz += vz[i] * shpdz[eb + i];

    // XZ shear component
    double sxz = 0;
    for (int i = 0; i < NPE; i++)
        sxz += 0.5 * (vx[i] * shpdz[eb + i] + vz[i] * shpdx[eb + i]);

    int sb = gid * 3;
    strain_rate[sb]     = sxx;
    strain_rate[sb + 1] = szz;
    strain_rate[sb + 2] = sxz;
}


// Update velocity for 2D (per node):  vel[i][j] += dt * force[i][j] / mass[i]
// vel/force: flat arrays, nnode * 2 doubles
// mass:      flat array,  nnode doubles
// dt:        scalar passed via a single-element constant buffer
kernel void update_velocity_2d(
    device double      *vel   [[buffer(0)]],
    device const double *force [[buffer(1)]],
    device const double *mass  [[buffer(2)]],
    constant double    &dt    [[buffer(3)]],
    uint gid                  [[thread_position_in_grid]])
{
    double m  = mass[gid];
    int    nb = gid * 2;
    vel[nb]     += dt * force[nb]     / m;
    vel[nb + 1] += dt * force[nb + 1] / m;
}


// Update coordinates for 2D (per node):  coord[i][j] += vel[i][j] * dt
kernel void update_coordinate_2d(
    device double      *coord [[buffer(0)]],
    device const double *vel  [[buffer(1)]],
    constant double    &dt   [[buffer(2)]],
    uint gid                 [[thread_position_in_grid]])
{
    int nb = gid * 2;
    coord[nb]     += vel[nb]     * dt;
    coord[nb + 1] += vel[nb + 1] * dt;
}


// ============================================================
// 3D kernels  (tetrahedral elements)
// ============================================================

// Compute tetrahedral volume for each 3D element.
// coord:   flat array, nnode * 3 doubles (x, y, z interleaved)
// conn:    flat array, nelem * 4 ints
// volume:  output, nelem doubles
kernel void compute_volume_3d(
    device const double *coord  [[buffer(0)]],
    device const int   *conn   [[buffer(1)]],
    device double      *volume [[buffer(2)]],
    uint gid                   [[thread_position_in_grid]])
{
    int n0 = conn[gid * 4];
    int n1 = conn[gid * 4 + 1];
    int n2 = conn[gid * 4 + 2];
    int n3 = conn[gid * 4 + 3];

    double x0 = coord[n0 * 3],     y0 = coord[n0 * 3 + 1],  z0 = coord[n0 * 3 + 2];
    double x1 = coord[n1 * 3],     y1 = coord[n1 * 3 + 1],  z1 = coord[n1 * 3 + 2];
    double x2 = coord[n2 * 3],     y2 = coord[n2 * 3 + 1],  z2 = coord[n2 * 3 + 2];
    double x3 = coord[n3 * 3],     y3 = coord[n3 * 3 + 1],  z3 = coord[n3 * 3 + 2];

    double x01 = x0 - x1,  x12 = x1 - x2,  x23 = x2 - x3;
    double y01 = y0 - y1,  y12 = y1 - y2,  y23 = y2 - y3;
    double z01 = z0 - z1,  z12 = z1 - z2,  z23 = z2 - z3;

    volume[gid] = (x01 * (y23 * z12 - y12 * z23) +
                   x12 * (y01 * z23 - y23 * z01) +
                   x23 * (y12 * z01 - y01 * z12)) / 6.0;
}


// Compute strain rate for each 3D tetrahedral element.
// vel:                  flat array, nnode * 3 doubles
// conn:                 flat array, nelem * 4 ints
// shpdx/shpdy/shpdz:    flat arrays, nelem * 4 doubles
// strain_rate:          output, nelem * 6 doubles
//                       order: sxx, syy, szz, sxy, sxz, syz
kernel void update_strain_rate_3d(
    device const double *vel         [[buffer(0)]],
    device const int   *conn        [[buffer(1)]],
    device const double *shpdx      [[buffer(2)]],
    device const double *shpdy      [[buffer(3)]],
    device const double *shpdz      [[buffer(4)]],
    device double      *strain_rate [[buffer(5)]],
    uint gid                        [[thread_position_in_grid]])
{
    const int NPE = 4;
    int eb = gid * NPE;

    double vx[4], vy[4], vz[4];
    for (int i = 0; i < NPE; i++) {
        int n = conn[eb + i];
        vx[i] = vel[n * 3];
        vy[i] = vel[n * 3 + 1];
        vz[i] = vel[n * 3 + 2];
    }

    double sxx = 0, syy = 0, szz = 0, sxy = 0, sxz = 0, syz = 0;
    for (int i = 0; i < NPE; i++) {
        double dx = shpdx[eb + i];
        double dy = shpdy[eb + i];
        double dz = shpdz[eb + i];
        sxx += vx[i] * dx;
        syy += vy[i] * dy;
        szz += vz[i] * dz;
        sxy += 0.5 * (vx[i] * dy + vy[i] * dx);
        sxz += 0.5 * (vx[i] * dz + vz[i] * dx);
        syz += 0.5 * (vy[i] * dz + vz[i] * dy);
    }

    int sb = gid * 6;
    strain_rate[sb]     = sxx;
    strain_rate[sb + 1] = syy;
    strain_rate[sb + 2] = szz;
    strain_rate[sb + 3] = sxy;
    strain_rate[sb + 4] = sxz;
    strain_rate[sb + 5] = syz;
}


// Update velocity for 3D (per node).
// vel/force: flat arrays, nnode * 3 doubles
// mass:      flat array,  nnode doubles
kernel void update_velocity_3d(
    device double      *vel   [[buffer(0)]],
    device const double *force [[buffer(1)]],
    device const double *mass  [[buffer(2)]],
    constant double    &dt    [[buffer(3)]],
    uint gid                  [[thread_position_in_grid]])
{
    double m  = mass[gid];
    int    nb = gid * 3;
    vel[nb]     += dt * force[nb]     / m;
    vel[nb + 1] += dt * force[nb + 1] / m;
    vel[nb + 2] += dt * force[nb + 2] / m;
}


// Update coordinates for 3D (per node).
kernel void update_coordinate_3d(
    device double      *coord [[buffer(0)]],
    device const double *vel  [[buffer(1)]],
    constant double    &dt   [[buffer(2)]],
    uint gid                 [[thread_position_in_grid]])
{
    int nb = gid * 3;
    coord[nb]     += vel[nb]     * dt;
    coord[nb + 1] += vel[nb + 1] * dt;
    coord[nb + 2] += vel[nb + 2] * dt;
}
