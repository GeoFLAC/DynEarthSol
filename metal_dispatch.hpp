#ifndef DYNEARTHSOL3D_METAL_DISPATCH_HPP
#define DYNEARTHSOL3D_METAL_DISPATCH_HPP

// Metal GPU dispatch interface for DynEarthSol
//
// Utilises Apple Silicon's unified memory: CPU and GPU share the same physical
// DRAM, so Metal buffers wrapping existing C++ heap pointers have zero-copy cost.
//
// Double-precision arithmetic in Metal shaders requires:
//   Apple GPU Family 8 (M2 / A16 Bionic) or later, or Mac 2 GPU family
//   (discrete AMD Radeon on Intel Macs).
// On M1 / older hardware metal_init() returns false and all dispatch functions
// become no-ops; callers fall back automatically to the CPU / OpenMP path.
//
// Build with:   make metal=1   (macOS only)
// This sets -DMETAL and links -framework Metal -framework Foundation.

#ifdef METAL

#ifdef __cplusplus
extern "C" {
#endif

// ---------------------------------------------------------------------------
// Lifecycle
// ---------------------------------------------------------------------------

// Initialise the Metal GPU backend.  Must be called once before any dispatch
// functions.  Compiles the Metal shaders from embedded source (Metal caches
// the result on disk, so only the first call incurs compilation overhead).
// Returns true if the GPU supports double-precision arithmetic, false otherwise
// (M1 and older); in that case all dispatch functions below are safe no-ops.
bool metal_init(void);

// Release all Metal objects.  Safe to call even if metal_init returned false.
void metal_cleanup(void);

// True when the GPU is ready and double-precision shaders are loaded.
bool metal_gpu_available(void);


// ---------------------------------------------------------------------------
// Geometry kernels
// ---------------------------------------------------------------------------

// Compute element volumes.
//   ndims = 2  =>  triangle area    (coord: nnode*2, conn: nelem*3)
//   ndims = 3  =>  tetrahedron vol  (coord: nnode*3, conn: nelem*4)
void metal_compute_volume(
    int ndims, int nelem, int nnode,
    const double *coord, const int *conn,
    double *volume);


// ---------------------------------------------------------------------------
// Field kernels
// ---------------------------------------------------------------------------

// Compute per-element strain rate.
//   vel:          nnode * ndims doubles
//   conn:         nelem * (ndims+1) ints
//   shpdx/shpdz:  nelem * (ndims+1) doubles  (shape-function derivatives)
//   shpdy:        only used when ndims == 3
//   strain_rate:  nelem * nstr doubles  (nstr = 3 for 2-D, 6 for 3-D)
void metal_update_strain_rate(
    int ndims, int nelem, int nnode,
    const double *vel, const int *conn,
    const double *shpdx, const double *shpdy, const double *shpdz,
    double *strain_rate);

// Update nodal velocities:  vel[i][j] += dt * force[i][j] / mass[i]
//   vel, force:  nnode * ndims doubles
//   mass:        nnode doubles
void metal_update_velocity(
    int ndims, int nnode,
    double *vel, const double *force, const double *mass,
    double dt);

// Update nodal coordinates:  coord[i][j] += vel[i][j] * dt
//   coord, vel:  nnode * ndims doubles
void metal_update_coordinate(
    int ndims, int nnode,
    double *coord, const double *vel,
    double dt);

#ifdef __cplusplus
}
#endif

#endif /* METAL */

#endif /* DYNEARTHSOL3D_METAL_DISPATCH_HPP */
