#ifndef DYNEARTHSOL3D_GEOMETRY_HPP
#define DYNEARTHSOL3D_GEOMETRY_HPP

template <typename T>
double dist2(T a, T b);
#pragma acc routine seq
double compute_volume(ConstArrayIndirectAccessor coord);
#pragma acc routine seq
double compute_area_facet(ConstArrayIndirectAccessor coord);
void compute_volume(const array_t &coord, const conn_t &connectivity,
                    double_vec &volume);
void compute_volume(const Variables &var, double_vec &volume);

void compute_dvoldt(const Variables &var, double_vec &dvoldt,
                    double_vec &etmp);

void compute_edvoldt(const Variables &var, double_vec &dvoldt,
                     double_vec &edvoldt);

void NMD_stress(const Variables &var, tensor_t& stress, double_vec &dp_nd,
                double_vec &etmp);

// Surface-relative reference pressure for SPR remeshing: ref_pressure is evaluated
// at the effective depth below the CURRENT surface instead of the fixed datum z = 0,
//     zeff(x,z) = z - h_eff(x,d),   h_eff = sum_m w_m a_m e^{-k_m d},
// where {a_m} is the DCT-I cosine series of the surface elevation z_surf(x)
// (k_m = pi*m/L) and d is the depth below the local surface. e^{-k d} is the
// mean-stress kernel of a harmonic load on an elastic half-space: the k = 0 mode
// (uniform shift) never attenuates, and d <= 0 keeps the full h(x) so the
// free-surface pin at ref_pressure(0) = 0 stays exact. The 3-D build rasterizes the
// top boundary onto a uniform (x,y) grid and uses a separable 2-D DCT-I with
// kernel e^{-|k| d}, |k| = pi*sqrt((m/Lx)^2 + (n/Ly)^2).
// No persistent state: consumers build local instances from the current mesh, so
// remeshing and restart need no extra plumbing. When inactive (no top nodes),
// elev() returns 0.0 and zeff() returns z -- fixed-datum behavior, bit-exact.
struct SurfaceTopo {
    bool on = false;
    bool atten = false;           // attenuation table below is filled
    double_vec x, z;              // 2-D: top-boundary node coordinates, sorted by x
    // h_eff table: 2-D on (x_i, d_j) [i*ndg + j]; 3-D on (x_i, y_j, d_l)
    // [(i*mgy + j)*ndg + l]. x, y uniform; depth sqrt-spaced d = dmax*(l/(ndg-1))^2
    // (fine near the surface, where the high-k modes die).
    int mg = 0, ndg = 0;
    double x0g = 0.0, Lg = 0.0, dmax = 0.0;
    double_vec heff;
#ifdef THREED
    int mgy = 0;
    double y0g = 0.0, Lyg = 0.0;
    double_vec gh;                // 3-D: rasterized surface elevation grid [i*mgy + j]
#endif
    void build(const Param& param, const Variables& var);
    // p is the full coordinate: double[NDIMS] or an array_t row (accessor overloads).
    // The 2-D/3-D branch is compile-time, so call sites are dimension-agnostic.
    double elev(const double* p) const;             // z_surf at p; 0.0 when !on
    double zeff(const double* p) const;             // reference-depth coordinate
    double heff_at(const double* p, double d) const;  // attenuation-table lookup
    double elev(ConstArrayAccessor p) const;
    double zeff(ConstArrayAccessor p) const;
    double heff_at(ConstArrayAccessor p, double d) const;
#ifndef THREED
    // Mean z_surf over [xa, xb], exact for the piecewise-linear profile. build()
    // resamples with this instead of elev() when the DCT grid is coarser than the
    // surface nodes, so unresolved relief is averaged rather than aliased.
    double elev_avg(double xa, double xb) const;
#endif
};

// The remesh stress remap, in call order. center_/restore_ are the reference
// round trip and always run; the spr_* pair is the recovery and is skipped when
// the rheology has no viscous component (then the remap is the NN interpolation
// plus the round trip, which is what an infinite Maxwell time asks for).
//
// Two orderings the compiler cannot enforce:
//  1. compute_spr_blend_weight BEFORE center_stress_to_ref -- visc() reads the
//     stress trace, so centering changes the weight it would produce.
//  2. restore_stress_from_ref LAST -- it undoes the centering for whatever
//     produced the final centered stress (the blend, or the NN remap alone).
// Each takes the SurfaceTopo of its own mesh, built once by the caller: old mesh
// for the weight/centering pair, new mesh for the recovery/restore pair.
void compute_spr_blend_weight(const Param& param, const Variables& var);

void center_stress_to_ref(const Param& param, const Variables& var,
                          const SurfaceTopo& topo);

void spr_elem_to_node(const Param& param, const Variables& var,
                      tensor_t* stress_n, double_vec* stressyy_n);

void spr_node_to_elem(const Param& param, const Variables& var,
                      const SurfaceTopo& topo,
                      tensor_t* stress, double_vec* stressyy);

void restore_stress_from_ref(const Param& param, const Variables& var,
                             const SurfaceTopo& topo,
                             tensor_t* stress, double_vec* stressyy);

double compute_dt(const Param& param, Variables& var);

// double compute_dt_PT(const Param& param, const Variables& var);

void compute_mass(const Param &param, const Variables &var,
                  double max_vbc_val, double_vec &volume_n,
                  double_vec &mass, double_vec &tmass, double_vec &hmass, double_vec &ymass,

                  elem_cache &tmp_result);


double elem_quality(const array_t &coord, const conn_t &connectivity,
                    const double_vec &volume, int e);

double worst_elem_quality(const array_t &coord, const conn_t &connectivity,
                          const double_vec &volume, int &worst_elem);

#endif
