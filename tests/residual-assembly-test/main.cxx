// Focused tests of real DES routines: force assembly in update_force and its
// undamped output, PT velocity constraints and residual norm, material trials and
// Mohr-Coulomb returns, the 3x3 eigensolver, and hydraulic, thermal and
// point-source updates.
// Modes: -DBASE (5-arg update_force) / -DB0 (6-arg, with the force_undamped output).
// Mesh: two elements sharing node 0, identical constant stress -> cancellation at node 0.
// Force cases: cancellation, Neumann traction, damping exclusion, body force.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>

#include "parameters.hpp"
#include "fields.hpp"
#include "hydraulics.hpp"
#include "point_sources.hpp"
#include "bc.hpp"
#include "rheology.hpp"
#include "pseudo_transient.hpp"
#include <limits>
#include "matprops.hpp"
#include "geometry.hpp"
#include "3x3-C/dsyevh3.h"

#ifndef ATOL
#define ATOL 1e-9
#endif
#ifndef RTOL
#define RTOL 1e-14
#endif

static int nfail = 0;
static int npass = 0;

static bool close(double x, double y) {
    return std::fabs(x - y) <= ATOL + RTOL * std::fabs(x) + RTOL * std::fabs(y);
}

static void check(const char* what, double got, double expected) {
    if (close(got, expected)) { ++npass; }
    else { ++nfail; std::printf("FAIL %s: got %.17g expected %.17g\n", what, got, expected); }
}

struct Fixture {
    Variables var;
    Param param{};
    array_t coord, vel, force, force_residual;
#ifdef B0
    array_t force_undamped;
#endif
    conn_t connectivity;
    tensor_t stress;
    elem_cache tmp_result;
    double_vec volume;
    double_vec temperature, pressure, pressure_increment, lookup;
    double_vec thermal_mass, heat_source, old_temperature, hydraulic_mass;
    double_vec volume_rate, balance, matrix_pressure;
    uint_vec boundary_flags;
    tensor_t strain_rate, strain;
    int_vec2D markers;
    int_vec etmp_int;

    Fixture() {}
};

static void build_mesh(Fixture& m) {
    // Node layout: global 0 is shared by both elements (local 0 in each).
    // Elem 0 non-shared nodes: global 1..NDIMS (unit simplex, positive orientation).
    // Elem 1 non-shared nodes: global NDIMS+1 .. 2*NDIMS (reflected through origin).
    const int nnode = 1 + 2 * NDIMS;
    const int nelem = 2;
    m.var.nnode = nnode;
    m.var.nelem = nelem;
    m.var.time = 0;
    m.var.steps = 0;

    m.coord.resize(nnode);
    m.vel.resize(nnode);
    m.force.resize(nnode);
    m.force_residual.resize(nnode);
#ifdef B0
    m.force_undamped.resize(nnode);
#endif
    m.connectivity.resize(nelem);
    m.stress.resize(nelem);
    m.tmp_result.resize(nelem);
    m.volume.assign(nelem, 0.0);

    // Element 0: unit simplex in positive octant.
    for (int i = 0; i < NODES_PER_ELEM; ++i)
        for (int d = 0; d < NDIMS; ++d)
            m.coord[i][d] = (i == d + 1) ? 1.0 : 0.0;
    for (int i = 0; i < NODES_PER_ELEM; ++i)
        m.connectivity[0][i] = i;
    m.volume[0] = (NDIMS == 2) ? 0.5 : (1.0 / 6.0);

    // Element 1: reflected through origin. In 3D swap local 1 and 2 to keep
    // positive orientation (mirror flips handedness).
    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        int g = (i == 0) ? 0 : NDIMS + i;
        for (int d = 0; d < NDIMS; ++d)
            m.coord[g][d] = -(i == d + 1 ? 1.0 : 0.0);
    }
    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        if (NDIMS == 3 && i == 1) m.connectivity[1][i] = NDIMS + 2;
        else if (NDIMS == 3 && i == 2) m.connectivity[1][i] = NDIMS + 1;
        else m.connectivity[1][i] = (i == 0) ? 0 : NDIMS + i;
    }
    m.volume[1] = (NDIMS == 2) ? 0.5 : (1.0 / 6.0);

    for (int e = 0; e < nelem; ++e)
        for (int s = 0; s < NSTR; ++s)
            m.stress[e][s] = 0.0;
    for (int i = 0; i < nnode; ++i)
        for (int d = 0; d < NDIMS; ++d) {
            m.vel[i][d] = 0.0;
            m.force[i][d] = 0.0;
            m.force_residual[i][d] = 0.0;
#ifdef B0
            m.force_undamped[i][d] = 0.0;
#endif
        }

    // Full CSR support graph from connectivity.
    m.var.support.arr_data.clear();
    m.var.support.idx_data.assign(nnode + 1, 0);
    m.var.support.lidx_data.clear();
    for (int e = 0; e < nelem; ++e)
        for (int l = 0; l < NODES_PER_ELEM; ++l)
            m.var.support.idx_data[m.connectivity[e][l] + 1]++;
    for (int n = 1; n <= nnode; ++n)
        m.var.support.idx_data[n] += m.var.support.idx_data[n - 1];
    m.var.support.arr_data.assign(m.var.support.idx_data[nnode], 0);
    m.var.support.lidx_data.assign(m.var.support.idx_data[nnode], 0);
    std::vector<int> fill(m.var.support.idx_data.begin(), m.var.support.idx_data.end() - 1);
    for (int e = 0; e < nelem; ++e)
        for (int l = 0; l < NODES_PER_ELEM; ++l) {
            int n = m.connectivity[e][l];
            int pos = fill[n]++;
            m.var.support.arr_data[pos] = e;
            m.var.support.lidx_data[pos] = l;
        }
    m.var.support.rebind();

    m.var.coord = &m.coord;
    m.var.connectivity = &m.connectivity;
    m.var.stress = &m.stress;
    m.var.volume = &m.volume;
    m.var.vel = &m.vel;
    m.var.force = &m.force;
    m.var.force_residual = &m.force_residual;
    m.var.tmp_result = &m.tmp_result;
    m.var.strain_rate = &m.strain_rate;

    for (int i = 0; i < nbdrytypes; ++i) {
        m.var.vbc_types[i] = 0;
        m.var.vbc_values[i] = 0.0;
    }
    for (int i = 0; i < nbdrytypes_hydro; ++i) {
        m.var.stress_bc_types[i] = 0;
        m.var.stress_bc_values[i] = 0.0;
    }
    m.param.control.gravity = 0.0;
    m.param.control.damping_option = 0;
    m.param.control.damping_factor = 0.0;
    m.param.ic.has_body_force_adjustment = true;
}

// Uniform temperature, pore pressure and pending pressure increment on the
// fixture's nodes, one marker of material 0 per element, and var bound to them
// and to the lookup tables.
static void bind_fields(Fixture& m, double temperature, double pressure = 0,
                        double pressure_increment = 0) {
    m.temperature.assign(m.var.nnode, temperature);
    m.pressure.assign(m.var.nnode, pressure);
    m.pressure_increment.assign(m.var.nnode, pressure_increment);
    m.markers.assign(m.var.nelem, int_vec(1, 1));
    m.var.temperature = &m.temperature;
    m.var.ppressure = &m.pressure;
    m.var.dppressure = &m.pressure_increment;
    m.var.elemmarkers = &m.markers;
    m.var.log_table = m.var.tan_table = m.var.sin_table = &m.lookup;
}

// Build the material of var again from param.mat.
static void reset_material(Fixture& m) {
    delete m.var.mat;
    m.var.mat = new MatProps(m.param, m.var);
}

static void run_update_force(Fixture& m) {
#ifdef BASE
    update_force(m.param, m.var, m.force, m.force_residual, m.tmp_result);
#else
    update_force(m.param, m.var, m.force, m.force_residual, m.tmp_result, &m.force_undamped);
    // Supplying force_undamped must not change force or force_residual: reevaluate
    // the same assembly with the default null argument and compare exact bits.
    std::vector<double> force_before, residual_before;
    for (int n=0; n<m.var.nnode; ++n)
        for (int d=0; d<NDIMS; ++d) {
            force_before.push_back(m.force[n][d]);
            residual_before.push_back(m.force_residual[n][d]);
        }
    update_force(m.param, m.var, m.force, m.force_residual, m.tmp_result);
    for (int n=0; n<m.var.nnode; ++n)
        for (int d=0; d<NDIMS; ++d) {
            const int k = n*NDIMS+d;
            const double f = m.force[n][d], r = m.force_residual[n][d];
            if (std::memcmp(&f, &force_before[k], sizeof(double)) ||
                std::memcmp(&r, &residual_before[k], sizeof(double))) {
                ++nfail;
                std::printf("FAIL optional output changed assembly at %d,%d\n", n, d);
            } else ++npass;
        }
#endif
}

// Case 1: identical constant stress in both elements -> net internal force at shared node 0 is 0.
static void test_cancellation() {
    std::printf("[cancellation]\n");
    Fixture* m = new Fixture();
    build_mesh(*m);
    const double sx = 0.7, sy = -0.3, sz = 0.5, sxy = 0.2, sxz = -0.1, syz = 0.4;
    m->stress[0][0] = sx; m->stress[0][1] = sy; m->stress[0][2] = sz;
    if (NSTR > 3) { m->stress[0][3] = sxy; m->stress[0][4] = sxz; m->stress[0][5] = syz; }
    m->stress[1][0] = sx; m->stress[1][1] = sy; m->stress[1][2] = sz;
    if (NSTR > 3) { m->stress[1][3] = sxy; m->stress[1][4] = sxz; m->stress[1][5] = syz; }

    run_update_force(*m);

    // BASE records force_residual at the shared node, which each incident element
    // overwrites: it equals the last element's term, not the assembled sum, so
    // with identical stress in both elements it is nonzero while the true
    // assembled internal force is exactly 0.
#ifdef BASE
    {
        std::printf("  [measure] legacy force_residual@shared = (%.17g", m->force_residual[0][0]);
        for (int d = 1; d < NDIMS; ++d)
            std::printf(", %.17g", m->force_residual[0][d]);
        std::printf(")\n");
        bool residual_nonzero = false;
        for (int d = 0; d < NDIMS; ++d)
            if (!close(m->force_residual[0][d], 0.0)) residual_nonzero = true;
        if (!residual_nonzero) {
            ++nfail;
            std::printf("FAIL cancel legacy-residual-nonzero: expected nonzero legacy defect, got ~0\n");
        } else { ++npass; }
    }
#endif

    for (int d = 0; d < NDIMS; ++d) {
        char name[64];
        snprintf(name, sizeof name, "cancel f0[%d]", d);
        check(name, m->force[0][d], 0.0);
#ifdef B0
        snprintf(name, sizeof name, "cancel u0[%d]", d);
        check(name, m->force_undamped[0][d], 0.0);
#endif
    }
    delete m;
}

// Case 2: Neumann traction on facet away from shared node.
// Facet 0 of element 0 (nodes 1..NDIMS): positive outward normal components.
// type 1 -> x-direction only. Each facet node gets T*normal[d]/NODES_PER_FACET.
static void test_traction() {
    std::printf("[traction]\n");
    Fixture* m = new Fixture();
    build_mesh(*m);
    const double T = 0.9;
    m->var.bfacets[iboundx0] = new std::vector<std::pair<int,int>>;
    m->var.bfacets[iboundx0]->push_back(std::make_pair(0, 0));
    m->var.stress_bc_types[iboundx0] = 1;
    m->var.stress_bc_values[iboundx0] = T;
    m->param.ic.has_body_force_adjustment = false;

    run_update_force(*m);

    // Facet 0 of the unit simplex: area-weighted outward normal.
    // 2D edge (1,0)-(0,1): normal = (1,1), nx = 1. 3D face: normal = (0.5,0.5,0.5).
    const double nx = (NDIMS == 2) ? 1.0 : 0.5;
    const double add = T / NODES_PER_FACET;
    for (int l = 1; l < NODES_PER_ELEM; ++l) {
        int g = l; // facet 0 nodes are global 1..NDIMS
        char name[64];
        snprintf(name, sizeof name, "traction f%d[0]", g);
        check(name, m->force[g][0], add * nx);
        for (int d = 1; d < NDIMS; ++d) {
            snprintf(name, sizeof name, "traction f%d[%d]", g, d);
            check(name, m->force[g][d], 0.0);
        }
#ifdef B0
        snprintf(name, sizeof name, "traction u%d[0]", g);
        check(name, m->force_undamped[g][0], add * nx);
#endif
    }
    char name[64];
    snprintf(name, sizeof name, "traction f0[0]");
    check(name, m->force[0][0], 0.0);
#ifdef B0
    snprintf(name, sizeof name, "traction u0[0]");
    check(name, m->force_undamped[0][0], 0.0);
#endif
    m->var.stress_bc_types[iboundx0] = 4;
    m->var.stress_bc_values[iboundx0] = -T;
    run_update_force(*m);
    for (int d=0; d<NDIMS; ++d) {
        double resultant=0;
        for (int n=0; n<m->var.nnode; ++n) resultant+=m->force[n][d];
        check("oblique facet normal-stress resultant",resultant,-T*nx);
        for (int l=1; l<NODES_PER_ELEM; ++l)
            check("normal stress distributed over facet nodes",m->force[l][d],-add*nx);
    }
    delete m->var.bfacets[iboundx0];
    delete m;
}

static void test_platen() {
    std::printf("[Rigid platen momentum, resultant and virtual work]\n");
    auto* m=new Fixture();build_mesh(*m);
    const int nn=m->var.nnode,last=nn-1;
    auto& mass=m->hydraulic_mass;
    auto& flags=m->boundary_flags;
    mass.assign(nn,1); flags.assign(nn,0);
    mass[0]=2;mass[last]=5;
    flags[0]=flags[last]=BOUNDZ1;
    m->var.mass=&mass;m->var.bcflag=&flags;m->param.bc.vbc_z1=8;
    for (int n=0;n<nn;++n)
        for (int d=0;d<NDIMS;++d) {
            m->vel[n][d]=n+d+1;
            m->force[n][d]=2*n-d;
        }
    const array_t old_velocity=m->vel,old_force=m->force;
    const double common=(2*old_velocity[0][NDIMS-1]+5*old_velocity[last][NDIMS-1])/7;
    const double resultant=old_force[0][NDIMS-1]+old_force[last][NDIMS-1];
    apply_platen_velocity(m->param,m->var,m->vel);
    apply_platen_force(m->param,m->var,m->force);
    double power_velocity=0,power_force=0;
    for (int n=0;n<nn;++n)
        for (int d=0;d<NDIMS;++d) {
            const bool on_plate=(flags[n]&BOUNDZ1) && d==NDIMS-1;
            check("platen velocity and tangential freedom",m->vel[n][d],
                  on_plate ? common:old_velocity[n][d]);
            check("platen force and tangential freedom",m->force[n][d],
                  on_plate ? mass[n]*resultant/7:old_force[n][d]);
            power_velocity+=old_force[n][d]*m->vel[n][d];
            power_force+=m->force[n][d]*old_velocity[n][d];
        }
    check("platen projections preserve virtual work",power_velocity,power_force);
    const array_t projected=m->vel;
    apply_platen_velocity(m->param,m->var,m->vel);
    for (int n=0;n<nn;++n)
        for (int d=0;d<NDIMS;++d) {
            const double before=projected[n][d],after=m->vel[n][d];
            check("platen velocity projection is bitwise idempotent",
                  std::memcmp(&before,&after,sizeof(double))==0,1);
        }
    delete m;
}

// Case 3: damping_option=2 -> force damped, undamped output not.
static void test_damping() {
    std::printf("[damping]\n");
    Fixture* m = new Fixture();
    build_mesh(*m);
    const double sx = 0.6;
    m->stress[0][0] = sx;
    m->param.control.damping_option = 2;
    m->param.control.damping_factor = 0.25;
    for (int i = 0; i < m->var.nnode; ++i)
        for (int d = 0; d < NDIMS; ++d)
            m->vel[i][d] = 1.0;

    run_update_force(*m);

    // Node 1 (local 1 of elem 0): shpdx[1]=1, others 0 -> tr = sx*vol, force = -sx*vol.
    const double vol = (NDIMS == 2) ? 0.5 : (1.0 / 6.0);
    const double undamped = -sx * vol;
    const double damped = undamped * (1.0 - 0.25);
    char name[64];
    snprintf(name, sizeof name, "damping f1[0]");
    check(name, m->force[1][0], damped);
#ifdef B0
    snprintf(name, sizeof name, "damping u1[0]");
    check(name, m->force_undamped[1][0], undamped);
#endif
    update_force(m->param, m->var, m->force, m->force_residual,
                 m->tmp_result, nullptr, 0, false);
    check("PT damping bypass", m->force[1][0], undamped);
    apply_damping(m->param, m->var, m->force);
    #pragma acc wait
    check("accepted legacy damping", m->force[1][0], damped);
    delete m;
}

static void test_body_force() {
    std::printf("[body force]\n");
    Fixture* m = new Fixture();
    build_mesh(*m);
    bind_fields(*m, 273.0);
    m->strain_rate.resize(m->var.nelem);
    m->etmp_int.resize(m->var.nelem);
    m->var.strain_rate = &m->strain_rate;
    m->var.etmp_int = &m->etmp_int;
    // Suppress gravity-induced boundary tractions to isolate the body load.
    // Neumann tractions are exercised separately, using the real BC routine.
    for (int i=0; i<nbdrytypes; ++i) m->var.vbc_types[i] = 1;
    auto& p = m->param.mat;
    p.nmat = 1;
    p.rheol_type = MatProps::rh_elastic;
    p.rho0 = {2500.0}; p.alpha = {0.0}; p.porosity = {0.2};
    p.bulk_modulus = {1e9}; p.shear_modulus = {1e9};
    p.visc_exponent = {1.0}; p.visc_coefficient = {1.0};
    p.biot_coeff = {1.0}; p.heat_capacity = {1.0};
    p.therm_cond = {1.0}; p.fluid_bulk_modulus = {2e9}; p.bulk_modulus_s = {37e9};
    p.fluid_rho0 = {1200};
    p.hydraulic_perm = {1e-12}; p.fluid_visc = {1e-3};
    m->param.control.gravity = 9.8;
    m->param.control.has_hydraulic_diffusion = true;
    m->param.control.damping_option = 2;
    m->param.control.damping_factor = 0.25;
    reset_material(*m);
    run_update_force(*m);
    // Each simplex distributes its weight equally among its vertices. The
    // shared vertex receives two contributions, all other vertices one.
    const double weight = (2500.0*0.8 + 1200.0*0.2)*9.8*m->volume[0];
    for (int n=0; n<m->var.nnode; ++n)
        for (int d=0; d<NDIMS; ++d) {
            const double raw = d == NDIMS-1 ? -weight*(n == 0 ? 2 : 1)/NODES_PER_ELEM : 0;
            check("damped gravity", m->force[n][d], raw*0.75);
#ifdef B0
            check("undamped gravity", m->force_undamped[n][d], raw);
#endif
        }
    m->param.control.has_hydraulic_diffusion=true;
    m->param.control.is_quasi_static=false;
    m->param.control.inertial_scaling=1;
    m->var.max_vbc_val=m->var.max_global_vel_mag=0.01;
    check("saturated reference lithostat",ref_pressure(m->param,-2.5),2240*9.8*2.5);
    double_vec volume_n(m->var.nnode),mass(m->var.nnode),tmass(m->var.nnode),
               hmass(m->var.nnode),ymass(m->var.nnode);
    for(bool global : {false,true}) {
        m->param.control.use_global_velocity_scaling=global;
        compute_mass(m->param,m->var,0.01,volume_n,mass,tmass,hmass,ymass,m->tmp_result);
        #pragma acc wait
        for(int n=0;n<m->var.nnode;++n)
            check("dynamic saturated nodal mass",mass[n],2240*m->volume[0]*(n==0 ? 2 : 1)/NODES_PER_ELEM);
    }
    m->param.control.has_hydraulic_diffusion=false;
    run_update_force(*m);
    for(int n=0;n<m->var.nnode;++n) {
        const double dry_weight=2500*9.8*m->volume[0]*(n==0 ? 2 : 1)/NODES_PER_ELEM;
        check("dry body weight with nonzero porosity",m->force[n][NDIMS-1],-0.75*dry_weight);
    }
#ifdef B0
    // A horizontal unit foundation: uniform upward displacement must reduce
    // the vertical supporting force by rho*g*area*displacement. This checks
    // the physical spring stiffness independently of the PT iteration scheme.
    std::printf("[PT trial foundation load]\n");
    m->param.control.damping_option = 0;
    m->param.bc.has_winkler_foundation = true;
    m->param.mesh.zlength = 1;
    m->var.compensation_pressure = 40000;
    m->var.vbc_types[iboundz0] = 0;
    m->var.bfacets[iboundz0] = new std::vector<std::pair<int,int>>{{0,NDIMS}};
    m->var.bnodes[iboundz0] = new int_vec;
    for (int l=0;l<NODES_PER_FACET;++l)
        m->var.bnodes[iboundz0]->push_back(NODE_OF_FACET[NDIMS][l]);
    std::vector<double> coords;
    for (int n=0;n<m->var.nnode;++n) {
        for (int d=0;d<NDIMS;++d) coords.push_back(m->coord[n][d]);
        m->vel[n][NDIMS-1] = 0.25;
    }
    for (bool wet : {false,true}) {
        m->param.control.has_hydraulic_diffusion = wet;
        update_force(m->param,m->var,m->force,m->force_residual,m->tmp_result,&m->force_undamped);
        double baseline=0;
        for (int n=0;n<m->var.nnode;++n) baseline+=m->force_undamped[n][NDIMS-1];
        // Boundary assembly reuses tmp_result. The PT rounding floor needs the
        // node sums of the absolute element forces as assembled before that.
        auto* magnitude = new array_t(m->var.nnode);
        update_force(m->param,m->var,m->force,m->force_residual,m->tmp_result,
                     &m->force_undamped,0,false,magnitude);
        // Saturated elements weigh their solid and water; dry ones their solid only.
        const double rho = wet ? 2240 : 2500;
        for (int n=0;n<m->var.nnode;++n)
            for (int d=0;d<NDIMS;++d)
                check("node force magnitude survives boundary scratch",(*magnitude)[n][d],
                      d==NDIMS-1 ? rho*9.8*m->volume[0]*(n==0 ? 2 : 1)/NODES_PER_ELEM : 0);
        delete magnitude;
        const double area = NDIMS==2 ? 1 : 0.5;
        for (double dt : {0.5,1.0,0.5,0.0}) {
            update_force(m->param,m->var,m->force,m->force_residual,m->tmp_result,&m->force_undamped,dt);
            double total=0;
            for (int n=0;n<m->var.nnode;++n) total+=m->force_undamped[n][NDIMS-1];
            check("foundation incremental force",total-baseline,-rho*9.8*area*dt*0.25);
            int k=0;
            for (int n=0;n<m->var.nnode;++n) for (int d=0;d<NDIMS;++d)
                check("trial loads preserve coordinates",m->coord[n][d],coords[k++]);
        }
    }
    delete m->var.bfacets[iboundz0];
    delete m->var.bnodes[iboundz0];
#endif
    delete m->var.mat;
    delete m;
}


// Independent geometric oracles: identity, one normal, two normal constraints.
// Also check projector algebra and prescribed velocities using the real dispatch.
static void test_constraints() {
    std::printf("[PT velocity constraints]\n");
    Fixture* m = new Fixture();
    build_mesh(*m);
    auto* flags = new uint_vec(m->var.nnode, 0);
    auto* normals = new array_t(nbdrytypes, 0.0);
    auto* c = new VelocityConstraints(m->var.nnode);
    m->var.bcflag = flags;
    m->var.bnormals = normals;
    m->var.time = 0;
    m->var.vbc_val_z1_loading_period = 1;
    m->param.bc.vbc_period_x0_time_in_yr = {0,1};
    m->param.bc.vbc_period_x1_time_in_yr = {0,1};
    m->param.bc.vbc_period_x0_ratio = {1,1};
    m->param.bc.vbc_period_x1_ratio = {1,1};
    for (int k=0;k<4;++k) {
        m->var.vbc_vertical_div_x0[k] = m->var.vbc_vertical_div_x1[k] = k/3.;
        m->var.vbc_vertical_ratio_x0[k] = m->var.vbc_vertical_ratio_x1[k] = 1;
    }
    const double a = 0.6, b = 0.8;
    (*normals)[iboundn0][0] = a;
    (*normals)[iboundn0][NDIMS-1] = b;
    for (int scenario=0;scenario<5;++scenario) {
        std::printf("  scenario %d\n", scenario);
        m->param.bc.vbc_x0 = scenario == 2 ? 1 : 0;
        m->param.bc.vbc_val_x0 = 2;
        m->var.vbc_types[iboundx0] = m->param.bc.vbc_x0;
        m->var.vbc_types[iboundn0] = scenario == 3 ? 11 : (scenario == 4 ? 3 : 1);
        m->var.vbc_values[iboundn0] = 3;
        for (int n=0;n<m->var.nnode;++n)
            (*flags)[n] = scenario == 0 ? 0 : (1U << iboundn0) | (scenario == 2 ? BOUNDX0 : 0);
        build_velocity_constraints(m->param,m->var,*c);
        for (int n=0;n<m->var.nnode;++n) {
            const int rank = scenario == 0 ? NDIMS : (scenario == 2 ? NDIMS-2 : (scenario == 4 ? 0 : NDIMS-1));
            check("free rank",c->free_rank[n],rank);
            for (int d=0;d<NDIMS;++d) {
                double po = 0;
                for(int j=0;j<NDIMS;++j) {
                    double pp = 0;
                    for(int k=0;k<NDIMS;++k) pp += c->projector[n][d*NDIMS+k]*c->projector[n][k*NDIMS+j];
                    check("P squared",pp,c->projector[n][d*NDIMS+j]);
                    check("P symmetry",c->projector[n][d*NDIMS+j],c->projector[n][j*NDIMS+d]);
                    po += c->projector[n][d*NDIMS+j]*c->prescribed[n][j];
                }
                check("P offset",po,0);
                m->vel[n][d] = d+5;
            }
        }
        project_free_vectors(m->var,*c,m->vel,true);
        for(int n=0;n<m->var.nnode;++n) {
            if(scenario == 0) for(int d=0;d<NDIMS;++d) check("free identity",m->vel[n][d],d+5);
            else if(scenario == 3) check("horizontal normal",m->vel[n][0],3);
            else check("normal velocity",a*m->vel[n][0]+b*m->vel[n][NDIMS-1],3);
            if(scenario == 2) check("earlier normal retained",m->vel[n][0],2);
            for(int d=0;d<NDIMS;++d) m->force[n][d] = m->vel[n][d];
        }
        apply_vbcs(m->param,m->var,m->vel,false,true);
        for(int n=0;n<m->var.nnode;++n)
            for(int d=0;d<NDIMS;++d) check("affine fixed point",m->vel[n][d],m->force[n][d]);
    }
#ifdef THREED
    // Period-scaled normal motion must coexist with PT's homogeneous map;
    // lateral shear remains unscaled by the normal-velocity period table.
    m->var.time = 0.5*YEAR2SEC;
    m->param.bc.vbc_x0 = 0;
    m->param.bc.vbc_y0 = 0;
    m->param.bc.vbc_period_x0_ratio = {2,4};
    m->param.bc.vbc_period_y0_time_in_yr = {0,1};
    m->param.bc.vbc_period_y0_ratio = {4,6};
    m->param.bc.vbc_val_x0 = m->param.bc.vbc_val_y0 = 2;
    m->param.bc.vbc_val_x0_l = m->param.bc.vbc_val_y0_l = 7;
    for (int axis=0; axis<NDIMS-1; ++axis) {
        int& code = axis == 0 ? m->param.bc.vbc_x0 : m->param.bc.vbc_y0;
        for (int type : {1,3,4,5,6,7}) {
            code = type;
            for (int n=0;n<m->var.nnode;++n) (*flags)[n] = axis == 0 ? BOUNDX0 : BOUNDY0;
            for (bool homogeneous : {false,true}) {
                for (int n=0;n<m->var.nnode;++n)
                    for (int d=0;d<NDIMS;++d) m->vel[n][d] = 11;
                apply_vbcs(m->param,m->var,m->vel,homogeneous,true);
                #pragma acc wait
                const double prescribed = homogeneous ? 0 : 1;
                const double normal = type == 4 ? 11 : (type == 5 ? 0 : prescribed*2*(axis == 0 ? 3 : 5));
                const double shear = type == 1 ? 11 : (type == 4 || type == 5 ? prescribed*2 :
                                                      (type == 6 ? prescribed*7 : 0));
                for (int n=0;n<m->var.nnode;++n) {
                    check("period-scaled side normal",m->vel[n][axis],normal);
                    check("unscaled lateral shear",m->vel[n][1-axis],shear);
                    check("vertical side constraint",m->vel[n][NDIMS-1],type == 1 || type == 6 || type == 7 ? 11 : 0);
                }
            }
        }
        code = 0;
    }
#endif
    delete c; delete normals; delete flags; delete m;
}


// control.iterative_return_mapping_option of the material fixtures; main runs
// the Mohr-Coulomb tests with each option.
static int return_mapping_option = 0;

// One-material element fields for constitutive updates on the test mesh.
struct MaterialFixture {
    Fixture* m;
    tensor_t* strain;
    double_vec *yy, *dp, *visc, *pls, *dpls, *friction, *state, *div;
    MaterialFixture(int model, double alpha) {
        m = new Fixture();
        build_mesh(*m);
        const int ne=m->var.nelem;
        strain = new tensor_t(ne,0.);
        yy = new double_vec(ne,0.);
        dp = new double_vec(ne,0.);
        visc = new double_vec(ne,1e9);
        pls = new double_vec(ne,0.);
        dpls = new double_vec(ne,0.);
        friction = new double_vec(ne,0.6);
        state = new double_vec(ne,10.);
        div = new double_vec(ne,0.);
        m->strain_rate.resize(ne,0.);
        bind_fields(*m,300.,10.,-10.);
        m->var.strain=strain; m->var.stressyy=yy; m->var.dpressure=dp;
        m->var.viscosity=visc; m->var.plstrain=pls; m->var.delta_plstrain=dpls;
        m->var.dyn_fric_coeff=friction; m->var.state_variable=state;
        m->var.edvoldt=div; m->var.strain_rate=&m->strain_rate;
        m->var.dt=0.1; m->var.time=1;
        auto& p=m->param.mat;
        p.nmat=1; p.rheol_type=model; p.is_plane_strain=NDIMS==2;
        p.rho0={2500}; p.alpha={0}; p.porosity={0.2};
        p.bulk_modulus={1e9}; p.shear_modulus={1e9};
        p.visc_exponent={1}; p.visc_coefficient={1};
        p.visc_activation_energy={0}; p.visc_activation_volume={0};
        p.visc_min=p.visc_max=1e9;
        p.biot_coeff={alpha}; p.heat_capacity={1}; p.therm_cond={1}; p.fluid_bulk_modulus={2e9}; p.bulk_modulus_s={37e9};
        p.pls0={1}; p.pls1={2}; p.cohesion0=p.cohesion1={1e8};
        p.friction_angle0=p.friction_angle1={0};
        p.dilation_angle0=p.dilation_angle1={0}; p.tension_max=1e8;
        m->param.control.has_hydraulic_diffusion=true;
        m->param.control.is_using_mixed_stress=true;
        m->param.control.iterative_return_mapping_option=return_mapping_option;
        reset_material(*m);
    }
    ~MaterialFixture() {
        delete m->var.mat;
        delete strain; delete yy; delete dp; delete visc; delete pls; delete dpls;
        delete friction; delete state; delete div; delete m;
    }
    // Uniform principal trial stress with zero plastic history.
    void set_principal(const double trial[3]) {
        for(int e=0;e<m->var.nelem;++e) {
            for(int d=0;d<NDIMS;++d) m->stress[e][d]=trial[d];
            for(int d=NDIMS;d<NSTR;++d) m->stress[e][d]=0;
            (*yy)[e]=trial[2];
            (*pls)[e]=(*dpls)[e]=0;
        }
    }
    void set_material() { reset_material(*m); }
    void update() {
        update_stress(m->param,m->var,m->stress,*yy,*dp,*visc,*strain,*pls,*dpls,
                      m->strain_rate,m->pressure,m->pressure_increment,m->vel,
                      *friction,*state,true);
        #pragma acc wait
    }
};


static void test_material_trials() {
    std::printf("[PT material pressure and repeated trials]\n");
    for (int model : {MatProps::rh_elastic, MatProps::rh_maxwell, MatProps::rh_viscous,
                      MatProps::rh_ep, MatProps::rh_evp}) {
        for (double alpha : {0.,0.6,1.}) {
            MaterialFixture f(model, alpha);
            Fixture* m = f.m;
            const int ne=m->var.nelem, nn=m->var.nnode;
            tensor_t* strain=f.strain;
            double_vec *yy=f.yy, *dp=f.dp, *visc=f.visc, *pls=f.pls, *dpls=f.dpls,
                       *friction=f.friction, *state=f.state;
            auto* heat = new double_vec(ne,7.);
            auto* pending_heat = new double_vec(ne,17.);
            m->var.shear_heat = heat;
            m->var.pending_heat = pending_heat;
            m->var.pore_plastic_increment = new double_vec(ne,19.);
            m->var.pore_dilation = new double_vec(ne,23.);
            m->param.control.has_shear_heating = true;
            auto& p=m->param.mat;
            auto* baseline=new MechanicalState(ne);
            baseline->capture(m->var);
            if (model==MatProps::rh_elastic || model==MatProps::rh_ep) {
                update_stress(m->param,m->var,m->stress,*yy,*dp,*visc,*strain,*pls,*dpls,
                              m->strain_rate,m->pressure,m->pressure_increment,m->vel,
                              *friction,*state,true,true);
                #pragma acc wait
                for(int e=0;e<ne;++e) {
                    for(int d=0;d<NSTR;++d)
                        check("initial equilibrium ignores pending pressure",m->stress[e][d],0);
                    check("initial pressure yy unchanged",(*yy)[e],0);
                }
                for(int n=0;n<nn;++n) {
                    check("initial pressure level preserved",m->pressure[n],10);
                    check("initial pending pressure preserved",m->pressure_increment[n],-10);
                }
            }
            {
                MatProps *original_mat = m->var.mat;
                p.alpha = {1e-5};
                m->var.mat = new MatProps(m->param,m->var);
                auto* old_temperature = new double_vec(nn,299.);
                m->var.old_temperature = old_temperature;
                baseline->restore(m->var);
                update_stress(m->param,m->var,m->stress,*yy,*dp,*visc,*strain,*pls,*dpls,
                              m->strain_rate,m->pressure,m->pressure_increment,m->vel,*friction,*state,true);
                #pragma acc wait
                for(int e=0;e<ne;++e) {
                    check("reversible thermal and pressure shifts do not heat",(*heat)[e],0);
                    check("thermal trial retains pending energy",(*pending_heat)[e],17);
                }
                delete m->var.mat;
                m->var.mat = original_mat;
                m->var.old_temperature = nullptr;
                delete old_temperature;
                p.alpha = {0};
            }
            for (int trial=0;trial<5;++trial) {
                baseline->restore(m->var);
                update_stress(m->param,m->var,m->stress,*yy,*dp,*visc,*strain,*pls,*dpls,
                              m->strain_rate,m->pressure,m->pressure_increment,m->vel,*friction,*state,true);
                #pragma acc wait
                for(int e=0;e<ne;++e) {
                    for(int d=0;d<NSTR;++d) {
                        check("pressure-only total stress",m->stress[e][d],d<NDIMS?-10*alpha:0);
                        check("zero strain",(*strain)[e][d],0);
                    }
                    if(NDIMS==2) check("plane strain pressure yy",(*yy)[e],-10*alpha);
                    check("pressure trace",(*dp)[e],-30*alpha);
                    check("plastic history",(*pls)[e],0);
                    check("pressure-only increment has no shear heat",(*heat)[e],0);
                    check("constitutive trial retains pending energy",(*pending_heat)[e],17);
                    check("plastic increment",(*dpls)[e],0);
                    check("state unchanged",(*state)[e],10);
                    // Deliberately contaminate every snapshot scalar to check restore.
                    (*m->var.pore_plastic_increment)[e]=99;
                    (*yy)[e]=(*dp)[e]=(*visc)[e]=(*pls)[e]=(*dpls)[e]=(*friction)[e]=(*state)[e]=99;
                }
                baseline->restore(m->var);
                for(int e=0;e<ne;++e) {
                    check("restore stressyy exact",(*m->var.stressyy)[e]==baseline->stressyy[e],1);
                    check("restore dpressure exact",(*m->var.dpressure)[e]==baseline->dpressure[e],1);
                    check("restore viscosity exact",(*m->var.viscosity)[e]==baseline->viscosity[e],1);
                    check("restore plstrain exact",(*m->var.plstrain)[e]==baseline->plstrain[e],1);
                    check("restore delta_plstrain exact",(*m->var.delta_plstrain)[e]==baseline->delta_plstrain[e],1);
                    check("restore dyn_fric_coeff exact",(*m->var.dyn_fric_coeff)[e]==baseline->dyn_fric_coeff[e],1);
                    check("restore state_variable exact",(*m->var.state_variable)[e]==baseline->state_variable[e],1);
                    check("restore plastic volume exact",(*m->var.pore_plastic_increment)[e]==baseline->pore_plastic_increment[e],1);
                    check("trial retains recoverable dilation",(*m->var.pore_dilation)[e],23);
                }
            }
            // Nonzero shear increment: a discarded trial must not accumulate
            // total strain, while Maxwell/EVP use one physical relaxation interval.
            for (int trial=0;trial<5;++trial) {
                baseline->restore(m->var);
                for(int e=0;e<ne;++e) m->strain_rate[e][NDIMS]=0.02;
                update_stress(m->param,m->var,m->stress,*yy,*dp,*visc,*strain,*pls,*dpls,
                              m->strain_rate,m->pressure,m->pressure_increment,m->vel,*friction,*state,true);
                #pragma acc wait
                double shear=4e6;
                if(model==MatProps::rh_viscous) shear=4e7;
                if(model==MatProps::rh_maxwell || model==MatProps::rh_evp) shear/=1.05;
                for(int e=0;e<ne;++e) {
                    check("single shear strain increment",(*strain)[e][NDIMS],0.002);
                    check("constitutive shear",m->stress[e][NDIMS],shear);
                    check("friction restored",(*friction)[e],0.6);
                    check("history restored",(*state)[e],10);
                }
            }
            if ((model==MatProps::rh_ep || model==MatProps::rh_evp) && alpha==0) {
                double first_plastic=0, first_shear=0;
                for(int trial=0;trial<5;++trial) {
                    baseline->restore(m->var);
                    for(int e=0;e<ne;++e) m->strain_rate[e][NDIMS]=2;
                    update_stress(m->param,m->var,m->stress,*yy,*dp,*visc,*strain,*pls,*dpls,
                                  m->strain_rate,m->pressure,m->pressure_increment,m->vel,*friction,*state,true);
                    #pragma acc wait
                    if(trial==0) { first_plastic=(*pls)[0]; first_shear=m->stress[0][NDIMS]; }
                    check("plastic branch exercised",first_plastic>0,1);
                    for(int e=0;e<ne;++e) {
                        check("plastic history not accumulated",(*pls)[e],first_plastic);
                        check("one plastic increment",(*dpls)[e],first_plastic);
                        check("repeated plastic stress",m->stress[e][NDIMS],first_shear);
                        check("single yielded strain increment",(*strain)[e][NDIMS],0.2);
                    }
                }
            }
            if (model==MatProps::rh_elastic && alpha==0) {
                double first_strength=0;
                for(double slip : {1e-9,1.0,100.0}) {
                    double mu=0.6,theta=10;
                    const auto p=m->var.mat->plastic_props_rsf(0,0,slip,mu,theta,0.1,0,true);
                    const double anphi=p.strength.anphi;
                    if(slip==1e-9) first_strength=anphi;
                    check("initial RSF frozen friction",mu,0.6);
                    check("initial RSF frozen state",theta,10);
                    check("initial RSF independent of pseudo velocity",anphi,first_strength);
                }
            }

            if (model==MatProps::rh_elastic && alpha==0.6) {
                p.alpha={3e-5};
                f.set_material();
                auto& old_temperature=m->old_temperature;
                old_temperature.assign(nn,290);
                m->var.old_temperature=&old_temperature;
                m->pressure_increment.assign(nn,0);
                for(int trial=0;trial<3;++trial) {
                    baseline->restore(m->var);
                    m->strain_rate.resize(ne,0.0);
                    for(int e=0;e<ne;++e) for(int d=0;d<NSTR;++d) m->strain_rate[e][d]=0;
                    m->pressure_increment.assign(nn,0);
                    update_stress(m->param,m->var,m->stress,*yy,*dp,*visc,*strain,*pls,*dpls,
                                  m->strain_rate,m->pressure,m->pressure_increment,m->vel,*friction,*state,true);
                    #pragma acc wait
                    for(int e=0;e<ne;++e) {
                        for(int d=0;d<NDIMS;++d) check("thermal trial fixed-strain stress",m->stress[e][d],-3e5);
                        if(NDIMS==2) check("thermal plane-strain out-of-plane stress",(*yy)[e],-3e5);
                    }
                }
                m->var.old_temperature=nullptr; p.alpha={0};
                f.set_material();
                baseline->restore(m->var);
                m->pressure_increment.assign(nn,-10);
            }

            // A real failed solve must restore history after numerical updates,
            // including rejected yielded trials and the pending pressure input.
            if (alpha==0.6) {
                auto* ntmp=new double_vec(nn,0.);
                auto* etmp=new double_vec(ne,0.);
                auto* volume_n=new double_vec(nn,0.);
                auto* flags=new uint_vec(nn,0);
                auto* normals=new array_t(nbdrytypes,0.);
                m->var.ntmp=ntmp; m->var.etmp=etmp; m->var.volume_n=volume_n;
                m->var.bcflag=flags; m->var.bnormals=normals;
                m->var.tmp_result=&m->tmp_result;
                for(int n=0;n<nn;++n)
                    for(int k=0;k<m->var.support.size(n);++k)
                        (*volume_n)[n]+=m->volume[m->var.support.patch(n)[k]];
                m->var.vbc_val_z1_loading_period=1;
                m->param.bc.vbc_period_x0_time_in_yr={0,1};
                m->param.bc.vbc_period_x1_time_in_yr={0,1};
                m->param.bc.vbc_period_x0_ratio={1,1};
                m->param.bc.vbc_period_x1_ratio={1,1};
                for(int k=0;k<4;++k) {
                    m->var.vbc_vertical_div_x0[k]=m->var.vbc_vertical_div_x1[k]=k/3.;
                    m->var.vbc_vertical_ratio_x0[k]=m->var.vbc_vertical_ratio_x1[k]=1;
                }
                m->param.mesh.xlength=m->param.mesh.ylength=m->param.mesh.zlength=2;
                m->param.control.PT_CFL=0.25;
                m->param.control.PT_Re=14.90188239869415;
                m->param.control.PT_relative_tolerance=0;
                m->param.control.PT_absolute_tolerance=0;
                m->param.control.PT_stagnation_window=0;
                // Enabling diffusion must not turn mechanical trials into
                // physical surface steps, including rejected candidates.
                m->param.control.surface_process_option=1;
                auto* dh=new double_vec(nn,0.125);
                auto* dhacc=new double_vec(nn,0.375);
                m->var.surfinfo.ntop=nn;
                m->var.surfinfo.dh=dh;
                m->var.surfinfo.dhacc=dhacc;
                auto* after=new MechanicalState(ne);
                for(int option : {0,1}) for(bool running : {false,true}) for(int cap : {1,3}) {
                    m->param.control.PT_option=option;
                    m->param.control.PT_use_running_scale=running;
                    m->var.PT_initial_residual_max=0;
                    m->var.PT_skip_scale_update=true;
                    baseline->restore(m->var);
                    for(int n=0;n<nn;++n) for(int d=0;d<NDIMS;++d)
                        m->vel[n][d]=d==0 ? 2*m->coord[n][NDIMS-1] : 0;
                    m->param.control.PT_max_iter=cap;
                    m->var.PT_force_residual_valid=true;
                    const PTResult result=run_physical_step_pt(m->param,m->var);
                    #pragma acc wait
                    check("failed solve status",result.status==PTStatus::max_iterations,1);
                    check("failed solve invalidates diagnostic",!m->var.PT_force_residual_valid,1);
                    check("failed solve iterations",result.iterations,cap);
                    check("failed solve has imbalance",result.residual>0,1);
                    // Both criteria keep a rounding floor; the solve still fails above it.
                    check("rounding floor is finite and positive",std::isfinite(result.threshold) && result.threshold>0,1);
                    check("failed solve stays above its threshold",result.residual>result.threshold,1);
                    check("failed solve retains reference",m->var.PT_initial_residual_max,0);
                    check("failed solve retains remesh exclusion",m->var.PT_skip_scale_update,1);
                    after->capture(m->var);
                    for(int e=0;e<ne;++e) {
                        for(int d=0;d<NSTR;++d) {
                            check("rollback stress exact",after->stress[e][d]==baseline->stress[e][d],1);
                            check("rollback strain exact",after->strain[e][d]==baseline->strain[e][d],1);
                        }
                        check("rollback stressyy exact",after->stressyy[e]==baseline->stressyy[e],1);
                        check("rollback dpressure exact",after->dpressure[e]==baseline->dpressure[e],1);
                        check("rollback viscosity exact",after->viscosity[e]==baseline->viscosity[e],1);
                        check("rollback plstrain exact",after->plstrain[e]==baseline->plstrain[e],1);
                        check("rollback delta_plstrain exact",after->delta_plstrain[e]==baseline->delta_plstrain[e],1);
                        check("rollback dyn_fric_coeff exact",after->dyn_fric_coeff[e]==baseline->dyn_fric_coeff[e],1);
                        check("rollback state_variable exact",after->state_variable[e]==baseline->state_variable[e],1);
                        check("rollback heating power exact",after->shear_heat[e]==baseline->shear_heat[e],1);
                        check("failed PT does not deposit heat",(*pending_heat)[e],17);
                    }
                    for(int n=0;n<nn;++n) {
                        for(int d=0;d<NDIMS;++d)
                            check("rollback velocity exact",m->vel[n][d]==(d==0 ? 2*m->coord[n][NDIMS-1] : 0),1);
                        check("pending pressure retained",m->pressure_increment[n]==-10,1);
                        check("pressure retained",m->pressure[n]==10,1);
                        check("PT leaves surface increment untouched",(*dh)[n],0.125);
                        check("PT leaves accumulated surface change untouched",(*dhacc)[n],0.375);
                    }
                    check("physical time retained",m->var.time==1,1);
                    check("physical dt retained",m->var.dt==0.1,1);
                    check("stress allocation retained",m->var.stress==&m->stress,1);
                    check("strain allocation retained",m->var.strain==strain,1);
                    check("viscosity allocation retained",m->var.viscosity==visc,1);
                }
                delete dh; delete dhacc;
                delete after; delete ntmp; delete etmp; delete volume_n;
                delete flags; delete normals;
            }
            delete baseline; delete heat; delete pending_heat;
        }
    }
}

// The plastic-strain measure stored in history, computed independently of
// rheology.cxx: the invariant of the nd principal components (3D, and pure 2D
// with nd = 2, the two in-plane ones), or in plane strain half their range.
static double recorded_plastic_increment(const double plastic[3], int nd=3) {
    if(NDIMS==2 && nd==3)
        return 0.5*(std::max(plastic[0],std::max(plastic[1],plastic[2]))-
                    std::min(plastic[0],std::min(plastic[1],plastic[2])));
    double mean=0;
    for(int k=0;k<nd;++k) mean+=plastic[k];
    mean/=nd;
    double squared=mean*mean;
    for(int k=0;k<nd;++k) squared+=(plastic[k]-mean)*(plastic[k]-mean);
    return std::sqrt(squared/2);
}

// Plastic strain of a principal stress correction under the isotropic compliance of the
// moduli of material 0, K and G (lambda = K - 2G/3).
static void compliance_plastic(const Param& param, const double trial[3], const double result[3],
                               double plastic[3]) {
    const double K=param.mat.bulk_modulus[0], G=param.mat.shear_modulus[0];
    const double volume=(trial[0]+trial[1]+trial[2]-result[0]-result[1]-result[2])/(3*K);
    for(int k=0;k<3;++k) plastic[k]=(trial[k]-result[k]-(K-2*G/3)*volume)/(2*G);
}

// Pure 2D: the same with the in-plane moduli a1 = K + 4G/3 and a2 = K - 2G/3,
// inverted (determinant a1^2 - a2^2).
static void compliance_plastic_2d(const Param& param, const double trial[3], const double result[3],
                                  double plastic[3]) {
    const double K=param.mat.bulk_modulus[0], G=param.mat.shear_modulus[0];
    const double a1=K+4*G/3, a2=K-2*G/3, determinant=a1*a1-a2*a2;
    const double d0=trial[0]-result[0], d1=trial[1]-result[1];
    plastic[0]=(a1*d0-a2*d1)/determinant;
    plastic[1]=(a1*d1-a2*d0)/determinant;
    plastic[2]=0;
}

// The strength of element 0 at plastic strain pls, as update_stress takes it
// before a return.
static MatProps::StrengthPiece element_strength(const MaterialFixture& f, double pls) {
    return f.m->var.mat->strength_here(0,pls);
}

// The principal stresses of element 0 after a principal trial: in pure 2D the
// two of its in-plane tensor (and 0), otherwise its normal stresses.
static void element_principal(const MaterialFixture& f, double r[3]) {
    const auto& stress=f.m->stress;
    if(NDIMS==2 && !f.m->param.mat.is_plane_strain) {
        const double mean=0.5*(stress[0][0]+stress[0][1]);
        const double radius=std::hypot(0.5*(stress[0][0]-stress[0][1]),stress[0][2]);
        r[0]=mean-radius; r[1]=mean+radius; r[2]=0;
        return;
    }
    r[0]=stress[0][0]; r[1]=stress[0][1]; r[2]=NDIMS==2 ? (*f.yy)[0] : stress[0][2];
}

// Elastic-plastic material with c=100 kPa, phi=30 deg, K=65/3 GPa, G=10 GPa.
static void set_mc_material(MaterialFixture& f, double dilation) {
    auto& p=f.m->param.mat;
    p.bulk_modulus={65e9/3}; p.shear_modulus={1e10};
    p.cohesion0=p.cohesion1={1e5};
    p.friction_angle0=p.friction_angle1={30};
    p.dilation_angle0=p.dilation_angle1={dilation};
    p.tension_max=1e9;
    f.m->param.control.has_hydraulic_diffusion=false;
    for(int e=0;e<f.m->var.nelem;++e) for(int d=0;d<NSTR;++d) f.m->strain_rate[e][d]=0;
    f.div->assign(f.m->var.nelem,0.);
    f.set_material();
}

static void check_law_return(MaterialFixture& f, const char* name, const double trial[3],
                             double minimum_history=0, double initial_history=0,
                             bool start_angles=false);

static void test_mohr_coulomb_returns() {
    std::printf("[Mohr-Coulomb returns]\n");
    MaterialFixture f(MatProps::rh_ep, 0);
    Fixture* m=f.m;
    const int ne=m->var.nelem;
    auto& p=m->param.mat;
    auto& stress=m->stress;
    auto &yy=*f.yy, &pls=*f.pls, &dpls=*f.dpls;

    // Independent symmetric-edge states from the two active MC equations and
    // isotropic compliance.
    for(int dilation : {0,30}) {
        set_mc_material(f, dilation);
        for(int edge : {0,1}) for(int orientation=0;orientation<4;++orientation) {
            const double base = edge == 0 ? -1e6 : 1e5;
            const double unique = edge == 0 ? -1e5 : -1.4e5;
            const double lower = edge == 0 ?
                (dilation == 0 ? -949487.1659305394 : -1028287.187078898) :
                (dilation == 0 ? -102564.06460551017 : -145505.28461683675);
            const double upper = edge == 0 ?
                (dilation == 0 ? -201025.66813892123 : -227292.3418550404) :
                (dilation == 0 ? 81282.03230275508 : 66968.29229897956);
            const double target_base = edge == 0 ? lower : upper;
            const double target_unique = edge == 0 ? upper : lower;
            double normal[3]={0,0,0};
            if(orientation<3) normal[orientation]=1;
            else for(int d=0;d<NDIMS;++d) normal[d]=1/std::sqrt(double(NDIMS));
            for(int e=0;e<ne;++e) {
                for(int d=0;d<NDIMS;++d)
                    stress[e][d]=base+(unique-base)*normal[d]*normal[d];
                stress[e][NDIMS]=(unique-base)*normal[0]*normal[1];
#ifdef THREED
                stress[e][4]=(unique-base)*normal[0]*normal[2];
                stress[e][5]=(unique-base)*normal[1]*normal[2];
#endif
                yy[e]=base+(unique-base)*normal[2]*normal[2];
                pls[e]=dpls[e]=0;
            }
            f.update();
            const double edge_trial[3]={base,base,unique};
            const double edge_target[3]={target_base,target_base,target_unique};
            double plastic[3];
            compliance_plastic(m->param,edge_trial,edge_target,plastic);
            const double increment=recorded_plastic_increment(plastic);
            for(int e=0;e<ne;++e) {
                for(int d=0;d<NDIMS;++d)
                    check("MC edge normal stress",stress[e][d],
                        target_base+(target_unique-target_base)*normal[d]*normal[d]);
                check("MC edge shear",stress[e][NDIMS],(target_unique-target_base)*normal[0]*normal[1]);
#ifdef THREED
                check("MC edge rotated xz",stress[e][4],(target_unique-target_base)*normal[0]*normal[2]);
                check("MC edge rotated yz",stress[e][5],(target_unique-target_base)*normal[1]*normal[2]);
#else
                check("MC edge plane-strain yy",yy[e],target_base+(target_unique-target_base)*normal[2]*normal[2]);
#endif
                check("MC edge plastic history",pls[e],increment);
                check("MC edge plastic increment",dpls[e],increment);
            }
        }
    }

    // The elastic trial exceeds c/tan(phi), but its admissible associated
    // shear return lies below it; no tensile flow.
    {
        set_mc_material(f, 30);
        const double trial[3]={169e3,185e3,200e3};
        const double target[3]={165347.65655963385,170585.93935780309,170585.93935780309};
        double plastic[3];
        compliance_plastic(m->param,trial,target,plastic);
        const double increment=recorded_plastic_increment(plastic);
        f.set_principal(trial);
        f.update();
        for(int e=0;e<ne;++e) {
            for(int d=0;d<NDIMS;++d) check("shear return before the tensile cap",stress[e][d],target[d]);
            for(int d=NDIMS;d<NSTR;++d) check("uncut shear has no rotation",stress[e][d],0);
            check("uncut shear plastic history",pls[e],increment);
            check("uncut shear plastic increment",dpls[e],increment);
#ifndef THREED
            check("uncut shear yy",yy[e],target[2]);
#endif
        }
    }

    // Reconstruct trial stresses from prescribed plastic flow and final
    // stresses on yield, then rebuild strength from the stored history and
    // check a second, unloaded constitutive update.
    for(double slope : {-5e8,2e9}) for(int dilation : {0,30}) {
        p.pls0={0}; p.pls1={1e-4};
        set_mc_material(f, dilation);
        p.cohesion1={1e5+slope*1e-4};
        f.set_material();
        check("MC nonzero cohesion slope",element_strength(f,0).dcohesion,slope);
        const double psi=dilation==0 ? 1. : 3.;
        // face=0 is a single face; 1/2 repeat the low/high stress.
        for(int face=0;face<3;++face) for(double ratio : {1.,2.}) {
            if(face==0 && ratio!=1) continue;
            const double first=1e-6, second=ratio*first;
            double plastic[3]={-first,0,psi*first};
            if(face==1) { plastic[1]=-second; plastic[2]=psi*(first+second); }
            if(face==2) { plastic[0]=-(first+second); plastic[1]=psi*second; }
            const double increment=recorded_plastic_increment(plastic);
            const double strength=2*std::sqrt(3.)*(1e5+slope*increment);
            double target[3]={3*(-2e5)-strength,0,-2e5};
            target[1]=face==1 ? target[0] : face==2 ? target[2] : 0.5*(target[0]+target[2]);
            const double volume=plastic[0]+plastic[1]+plastic[2];
            double trial[3];
            for(int k=0;k<3;++k) trial[k]=target[k]+2e10*plastic[k]+15e9*volume;
            f.set_principal(trial);
            f.update();
            for(int e=0;e<ne;++e) {
                for(int d=0;d<NDIMS;++d)
                    check("MC tangent normal stress",stress[e][d],target[d]);
                for(int d=NDIMS;d<NSTR;++d) check("MC tangent has no rotation",stress[e][d],0);
#ifndef THREED
                check("MC tangent plane-strain yy",yy[e],target[2]);
#endif
                check("MC tangent plastic history",pls[e],increment);
                check("MC tangent plastic increment",dpls[e],increment);
                check("MC tangent stays in linear strength range",pls[e]<1e-4,1);
            }
            const MatProps::StrengthPiece stored=element_strength(f,pls[0]);
            const double high=NDIMS==2 ? yy[0] : stress[0][2];
            check("MC stored history reproduces yield",stress[0][0]-stored.strength.anphi*high+stored.strength.amc,0);
            f.update();
            for(int e=0;e<ne;++e) {
                for(int d=0;d<NDIMS;++d)
                    check("MC unloaded stress stays on yield",stress[e][d],target[d]);
#ifndef THREED
                check("MC unloaded plane-strain yy",yy[e],target[2]);
#endif
                check("MC unloaded plastic history",pls[e],increment);
                check("MC unloaded plastic increment",dpls[e],0);
            }
        }
    }

#ifndef THREED
    // Plane-strain trials at and on either side of the tensile cap 0. While the
    // uncut-trial edge return stays below the cap, zero-dilation shear flow
    // conserves trace: p1 = p2 = (trace + 2c*sqrt(3))/5, p0 = 3*p2 - 2c*sqrt(3).
    // Above it the cap is active; check_law_return checks those by properties.
    p.pls0={0}; p.pls1={1};
    set_mc_material(f, 0);
    p.tension_max=0;
    f.set_material();
    const double cases[][3]={{-1e6,2e6,2e6},{-1e6,-1e4,2e6},
        {-1e6,-1e4,0},{-1e6,-1e4,-1e-6},{-1e6,-1e4,1e-6}};
    for(const auto& trial : cases) {
        check_law_return(f,"plane-strain tension",trial);
        const double upper=(trial[0]+trial[1]+trial[2]+2e5*std::sqrt(3.))/5;
        if(upper>0) continue;
        for(int e=0;e<ne;++e) {
            check("MC uncut edge below the cap xx",stress[e][0],3*upper-2e5*std::sqrt(3.));
            check("MC uncut edge below the cap zz",stress[e][1],upper);
            check("MC uncut edge below the cap yy",yy[e],upper);
        }
    }
#endif
}

// Whether plastic is a nonnegative combination of the flows of the surfaces
// active at result: shear faces (low, high) with flow -e_low + Npsi*e_high,
// and tensile caps with flow e_k. In principal space three flows suffice.
// nd = 2 (pure 2D) has the one face (0,1) and caps on p0, p1.
static bool in_active_flow_cone(const double result[3], const double plastic[3],
                                double amc, double anphi, double anpsi, double cap,
                                double scale, int nd=3) {
    double flows[6][3];
    int n=0;
    const int faces[3][2]={{0,2},{1,2},{0,1}};
    const double tol=1e-9*scale;
    for(int i=nd==3 ? 0 : 2;i<3;++i) {
        const auto& face=faces[i];
        if(std::fabs(result[face[0]]-anphi*result[face[1]]+amc)<=tol) {
            double* m=flows[n++];
            m[0]=m[1]=m[2]=0; m[face[0]]=-1; m[face[1]]=anpsi;
        }
    }
    for(int k=0;k<nd;++k)
        if(result[k]>=cap-tol) {
            double* m=flows[n++];
            m[0]=m[1]=m[2]=0; m[k]=1;
        }
    const auto dot=[](const double* a, const double* b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; };
    const double size=std::sqrt(dot(plastic,plastic));
    if(size==0) return true;
    for(int mask=1;mask<(1<<n);++mask) {
        int idx[3], k=0;
        for(int i=0;i<n && k<=3;++i) if(mask>>i&1) { if(k<3) idx[k]=i; ++k; }
        if(k>3) continue;
        // Least squares on the chosen flows: normal equations, Gaussian elimination.
        double a[3][4];
        for(int i=0;i<k;++i) {
            for(int j=0;j<k;++j) a[i][j]=dot(flows[idx[i]],flows[idx[j]]);
            a[i][k]=dot(flows[idx[i]],plastic);
        }
        bool singular=false;
        for(int c=0;c<k && !singular;++c) {
            int pivot=c;
            for(int r=c+1;r<k;++r) if(std::fabs(a[r][c])>std::fabs(a[pivot][c])) pivot=r;
            if(std::fabs(a[pivot][c])<=1e-12) { singular=true; break; }
            for(int j=0;j<=k;++j) std::swap(a[c][j],a[pivot][j]);
            for(int r=0;r<k;++r) if(r!=c) {
                const double factor=a[r][c]/a[c][c];
                for(int j=c;j<=k;++j) a[r][j]-=factor*a[c][j];
            }
        }
        if(singular) continue;
        double x[3], residual[3]={plastic[0],plastic[1],plastic[2]}, largest=0;
        for(int i=0;i<k;++i) { x[i]=a[i][k]/a[i][i]; largest=std::max(largest,std::fabs(x[i])); }
        bool nonnegative=true;
        for(int i=0;i<k;++i) {
            nonnegative=nonnegative && x[i]>=-1e-9*largest;
            for(int d=0;d<3;++d) residual[d]-=x[i]*flows[idx[i]][d];
        }
        if(nonnegative && std::sqrt(dot(residual,residual))<=1e-6*size) return true;
    }
    return false;
}

// Properties of any accepted return from a principal trial with no plastic
// history, checked against the law at the stored history: inside the shear
// surface and tensile cap and on one of them, principal order kept, history
// equal to the measure of the plastic strain implied by isotropic compliance,
// and no change from a second, unloaded update. With start_angles (option 0
// with friction or dilation weakening and constant cohesion) the strength is
// the one at the initial history instead, and a second update may yield again.
static void check_law_return(MaterialFixture& f, const char* name, const double trial[3],
                             double minimum_history, double initial_history,
                             bool start_angles) {
    const auto flag=[name](const char* property, bool ok) {
        char what[256];
        std::snprintf(what,sizeof what,"%s: %s",name,property);
        check(what,ok,1);
    };
    // Pure 2D has two principal stresses, from the in-plane tensor.
    const int nd=NDIMS==2 && !f.m->param.mat.is_plane_strain ? 2 : 3;
    f.set_principal(trial);
    for(int e=0;e<f.m->var.nelem;++e) (*f.pls)[e]=initial_history;
    f.update();
    double result[3];
    element_principal(f,result);
    double plastic[3];
    if(nd==2) compliance_plastic_2d(f.m->param,trial,result,plastic);
    else compliance_plastic(f.m->param,trial,result,plastic);
    const double recorded=recorded_plastic_increment(plastic,nd);
    const double history=(*f.pls)[0];
    const MatProps::StrengthPiece piece=element_strength(f,start_angles ? initial_history : history);
    const double amc=piece.strength.amc, anphi=piece.strength.anphi, anpsi=piece.strength.anpsi,
                 cap=piece.strength.ten_max;
    double scale=std::fabs(amc);
    for(int k=0;k<nd;++k) scale=std::max(scale,std::fabs(result[k]));
    const double shear=result[0]-anphi*result[nd-1]+amc, tension=cap-result[nd-1];
    bool ordered=true;
    for(int k=0;k+1<nd;++k) ordered=ordered && result[k]<=result[k+1]+1e-9*scale;
    flag("principal order kept",ordered);
    flag("inside the shear surface",shear>=-1e-9*scale);
    flag("inside the tensile cap",tension>=-1e-9*scale);
    flag("on the composite surface",std::min(shear,tension)<=1e-9*scale);
    const double increment=(*f.dpls)[0];
    flag("history is the initial history plus the increment",history==initial_history+increment);
    flag("increment is the recorded plastic strain",std::fabs(increment-recorded)<=1e-6*recorded);
    flag("history reaches the expected range",history>=minimum_history);
    flag("plastic strain in the cone of active flows",
         in_active_flow_cone(result,plastic,amc,anphi,anpsi,cap,scale,nd));
    const double size=std::sqrt(plastic[0]*plastic[0]+plastic[1]*plastic[1]+plastic[2]*plastic[2]);
    flag("nonnegative plastic work",
         result[0]*plastic[0]+result[1]*plastic[1]+result[2]*plastic[2]>=-1e-9*scale*size);
    if(start_angles) return;
    f.update();
    double again[3];
    element_principal(f,again);
    for(int k=0;k<nd;++k) flag("unloaded stress unchanged",std::fabs(again[k]-result[k])<=1e-9*scale);
    flag("unloaded increment is zero",(*f.dpls)[0]<=1e-12*increment);
}

// Cohesion that jumps at one plastic strain (pls0 = pls1), for a trial whose return
// at 1e5 Pa records inc: a drop to 1e4 at 1e-5 (the return goes on past it), and a
// drop and a rise between 1e4 and 1e5 at inc itself, where the return lands on the
// jump. The strength there is the next piece's: past the drop, on the stronger side
// of the rise with the increment inc.
static void check_cohesion_jumps(MaterialFixture& f, const char* prefix, const double trial[3]) {
    auto& p=f.m->param.mat;
    char name[128];
    set_mc_material(f, 0);
    f.set_principal(trial);
    for(int e=0;e<f.m->var.nelem;++e) (*f.pls)[e]=0;
    f.update();
    const double inc=(*f.dpls)[0];
    p.pls0=p.pls1={1e-5};
    p.cohesion1={1e4};
    f.set_material();
    std::snprintf(name,sizeof name,"%scohesion drop at one plastic strain",prefix);
    check_law_return(f,name,trial,1e-5);
    p.pls0=p.pls1={inc};
    f.set_material();
    std::snprintf(name,sizeof name,"%scohesion drop where the return lands",prefix);
    check_law_return(f,name,trial,inc);
    p.cohesion0={1e4}; p.cohesion1={1e5};
    f.set_material();
    std::snprintf(name,sizeof name,"%scohesion rise where the return lands",prefix);
    check_law_return(f,name,trial,inc);
    std::snprintf(name,sizeof name,"%scohesion rise where the return lands: history",prefix);
    check(name,(*f.pls)[0],inc);
}

// Two materials of the Mohr-Coulomb test material, mixed 1:1 in every element.
static void mix_two_materials(MaterialFixture& f, bool pure_2d) {
    auto& p=f.m->param.mat;
    if(pure_2d) p.is_plane_strain=false;
    set_mc_material(f, 0);
    p.nmat=2;
    for(auto* v : {&p.rho0,&p.alpha,&p.porosity,&p.bulk_modulus,&p.shear_modulus,
                   &p.visc_exponent,&p.visc_coefficient,&p.visc_activation_energy,
                   &p.visc_activation_volume,&p.biot_coeff,&p.heat_capacity,&p.therm_cond,
                   &p.fluid_bulk_modulus,&p.friction_angle0,&p.friction_angle1,
                   &p.dilation_angle0,&p.dilation_angle1})
        v->resize(2,v->front());
    for(auto& marks : f.m->markers) marks={1,1};
}

// Two materials mixed 1:1 whose cohesions jump the other way: the mean rises from
// 4e5 to 6.5e5 Pa at 1e-5 and falls to 3.25e5 Pa at 1.4e-5, at once or over
// [1.1e-5, 1.4e-5], for a trial whose return records more than 1.4e-5 at the residual
// strength. No increment records itself at the rise; the return lies past the fall.
static void check_mixed_cohesions(const char* prefix, const double trial[3], bool pure_2d) {
    MaterialFixture f(MatProps::rh_ep, 0);
    auto& p=f.m->param.mat;
    mix_two_materials(f, pure_2d);
    char name[128];
    for(double fall : {1.4e-5,1.1e-5}) {
        p.pls0={1e-5,fall}; p.pls1={1e-5,1.4e-5};
        p.cohesion0={5e4,7.5e5}; p.cohesion1={5.5e5,1e5};
        f.set_material();
        std::snprintf(name,sizeof name,"%smixed cohesions that rise and then fall%s",prefix,
                      fall<1.4e-5 ? " over a range" : "");
        check_law_return(f,name,trial,1.4e-5);
    }
}

// g = recorded - delta can turn along a piece. With option 1, friction follows the increment.
// Two materials mixed 1:1: past a cohesion jump at 1.04e-6, one hardens from 3e5 to 2.1e6 Pa
// while its friction falls from 40 to 0 deg by 8.79e-6. For the trial rise, g is below 0 at
// both ends of that piece and above 0 between, where the return lies. Past a jump of one from
// 5e3 to 1e5 Pa at 1e-9, the other hardens from 1.59e5 to 1.84214e6 Pa while its friction
// falls from 58 to 26 deg by 1.2e-8: the trial rise yields only from 8.72e-9 to 1.10e-8, not
// at the ends of that piece nor where a golden-section search starts, and the return lies
// there. For ten times the trial rise, past a jump of one from 5e4 to 1e6 Pa, both at 44 deg,
// the other hardens from 1.62952e6 to 1.841807e7 Pa while its friction falls to 12 deg by
// 3e-8: the trial yields only from 2.38e-8 to 2.71e-8, where cohesion and friction nearly
// cancel in its yield function. Both soften friction from 40 to 15 deg and cohesion from 1.5e6
// to 5e5 Pa over 1e-4: for the trial dip, g is above 0 at both ends and below 0 between, and
// the return is its first root. With K = 48.18e9 Pa, a cap of 1.098e6 Pa and dilation, past a
// jump of one from 17e3 to 842e3 Pa at 6.694e-5, the other softens from 1.998e6 to 9.4e4 Pa
// while its friction rises from 22.75 to 34.33 deg and its dilation falls from 18.59 to
// 2.48 deg by 1.511e-4 (3D): for the trial apex, the return goes from two faces and the cap to
// two faces, all three faces and back; g is below 0 at both ends of that piece and above 0
// between, and the return lies at the apex of the faces. With either option, g turns where the
// return moves from a face to an edge (3D and plane strain): past a cohesion rise from 6.5e5
// to 7.5e5 Pa at 5.2e-6, the other material falls from 7.5e5 Pa at 4.94e10 Pa per unit plastic
// strain to 1.7e-5 or 2e-5. For the trial edge, g rises above 0 on the face, falls to 0 on the
// edge and rises again, above 0 at 2e-5.
static void check_turning_returns(const char* prefix, const double rise[3], const double dip[3],
                                  const double edge[3], bool pure_2d) {
    MaterialFixture f(MatProps::rh_ep, 0);
    auto& p=f.m->param.mat;
    mix_two_materials(f, pure_2d);
    // The return lies on the piece from start to end.
    const auto check_on_piece=[&](const char* what, const double trial[3], double start, double end) {
        char name[128];
        std::snprintf(name,sizeof name,"%s%s",prefix,what);
        check_law_return(f,name,trial,start);
        std::snprintf(name,sizeof name,"%s%s: on that piece",prefix,what);
        check(name,(*f.pls)[0]<end,1);
    };
    if(return_mapping_option==1) {
        p.pls0={1.04e-6,1.04e-6}; p.pls1={1.04e-6,8.79e-6};
        p.cohesion0={5e4,3e5}; p.cohesion1={3e5,2.1e6};
        p.friction_angle0={40,40}; p.friction_angle1={40,0};
        f.set_material();
        check_on_piece("return where g rises above 0 past a jump",rise,1.04e-6,8.79e-6);
        p.pls0={1e-9,1e-9}; p.pls1={1e-9,1.2e-8};
        p.cohesion0={5e3,1.59e5}; p.cohesion1={1e5,1.84214e6};
        p.friction_angle0={30,58}; p.friction_angle1={30,26};
        f.set_material();
        check_on_piece("return where the trial yields again within a piece",rise,1e-9,1.2e-8);
        const double strong[3]={10*rise[0],10*rise[1],10*rise[2]};
        p.pls1={1e-9,3e-8};
        p.cohesion0={5e4,1.62952e6}; p.cohesion1={1e6,1.841807e7};
        p.friction_angle0={44,44}; p.friction_angle1={44,12};
        f.set_material();
        check_on_piece("return where cohesion and friction nearly cancel",strong,1e-9,3e-8);
        p.pls0={0,0}; p.pls1={1e-4,1e-4};
        p.cohesion0={1.5e6,1.5e6}; p.cohesion1={5e5,5e5};
        p.friction_angle0={40,40}; p.friction_angle1={15,15};
        f.set_material();
        check_on_piece("first root where g dips below 0",dip,0,1e-4);
#ifdef THREED
        {
            const auto mixed=p;
            const double start=6.694e-5, end=1.511e-4, apex[3]={-1.373e6,3.531e6,3.801e6};
            p.bulk_modulus={48.18e9,48.18e9}; p.tension_max=1.098e6;
            p.pls0={start,start}; p.pls1={start,end};
            p.cohesion0={17e3,1998e3}; p.cohesion1={842e3,94e3};
            p.friction_angle0={22.75,22.75}; p.friction_angle1={22.75,34.33};
            p.dilation_angle0={2.48,18.59}; p.dilation_angle1={2.48,2.48};
            f.set_material();
            check_on_piece("return where the active set changes and comes back",apex,start,end);
            p=mixed;
        }
#endif
    }
    if(!edge) return;
    for(double end : {1.7e-5,2e-5}) {
        p.pls0={5.2e-6,5.2e-6}; p.pls1={5.2e-6,end};
        p.cohesion0={6.5e5,7.5e5}; p.cohesion1={7.5e5,7.5e5-4.94e10*(end-5.2e-6)};
        p.friction_angle0=p.friction_angle1={30,30};
        f.set_material();
        char what[96];
        std::snprintf(what,sizeof what,"first root where g turns at an edge, to %g",end);
        check_on_piece(what,edge,5.2e-6,end);
    }
}

#ifndef THREED
// Pure 2D (not plane strain): the single-face return softens by exactly the
// increment it stores, so the strength rebuilt from history keeps it on yield.
static void test_pure_2d_softening_return() {
    std::printf("[Pure 2D Mohr-Coulomb softening return]\n");
    MaterialFixture f(MatProps::rh_ep, 0);
    auto& p=f.m->param.mat;
    p.is_plane_strain=false;
    for(int dilation : {0,30}) {
        p.pls0={0}; p.pls1={1e-4};
        set_mc_material(f, dilation);
        p.cohesion1={1e5-5e8*1e-4};
        f.set_material();
        const double trial[3]={-1e6,-1e5,0};
        f.set_principal(trial);
        f.update();
        const double history=(*f.pls)[0];
        const MatProps::StrengthPiece piece=element_strength(f,history);
        const double amc=piece.strength.amc, anphi=piece.strength.anphi;
        double r[3];
        element_principal(f,r);
        check("pure 2D softening stays in the linear range",history>0 && history<1e-4,1);
        check("pure 2D history is the increment",history,(*f.dpls)[0]);
        check("pure 2D stored history reproduces yield",(r[0]-anphi*r[1]+amc)/amc,0);
    }
}
#endif

// Returns the closed form cannot give: a cohesion drop faster than elastic
// unloading (also from a history one ulp short of the breakpoint), friction
// and dilation softening within the step, and the tensile cap with its corners:
// in 3D the tensile edge, apex and tension-shear corners, in pure 2D (two
// principal stresses, not plane strain) tension, the face-cap corner and the
// cap-cap corner. A pure 2D trial keeps the lowest and highest principal stress
// of the 3D one.
static void test_mohr_coulomb_law_returns(bool pure_2d) {
    std::printf("[%sMohr-Coulomb returns over the strength law]\n",pure_2d ? "Pure 2D " : "");
    MaterialFixture f(MatProps::rh_ep, 0);
    auto& p=f.m->param.mat;
    if(pure_2d) p.is_plane_strain=false;
    using Trial=std::array<double,3>;
    const auto in_plane=[pure_2d](const Trial& t) { return pure_2d ? Trial{t[0],t[2],0} : t; };
    const char* prefix=pure_2d ? "pure 2D " : "";
    char name[128];
    const auto named=[&](const char* what) {
        std::snprintf(name,sizeof name,"%s%s",prefix,what);
        return name;
    };
    const Trial edge={-1e6,-1e6,-1e5}, face=in_plane({-1e6,-5e5,-1e5});
    // Cohesion drops 1e5 -> 9e4 over 1e-8: the return reaches the residual strength.
    for(int dilation : {0,30}) {
        p.pls0={0}; p.pls1={1e-8};
        set_mc_material(f, dilation);
        p.cohesion1={9e4};
        f.set_material();
        if(!pure_2d) check_law_return(f,"softening past the segment, edge",edge.data(),1e-8);
        check_law_return(f,named("softening past the segment, face"),face.data(),1e-8);
    }
    // Nonzero initial history whose pls + (pls1 - pls) rounds one ulp below
    // pls1: the search must still cross the breakpoint to the residual root.
    p.pls0={0}; p.pls1={0.0009287493630057339};
    set_mc_material(f, 0);
    p.cohesion1={9e4};
    f.set_material();
    const Trial deep=in_plane({-1e8,-1e8,-1e7});
    check_law_return(f,named("softening past the segment from nonzero history"),deep.data(),
                     p.pls1[0],0.0003829357836281204);
    const Trial far=in_plane({-6e6,-2e6,-1e6});
    check_cohesion_jumps(f,prefix,far.data());
    check_mixed_cohesions(prefix,far.data(),pure_2d);
    const Trial dip=in_plane({-7e6,-3e6,-1e5}), near_edge={-6e6,-5.8e6,-1e6};
    check_turning_returns(prefix,far.data(),dip.data(),pure_2d ? nullptr : near_edge.data(),
                          pure_2d);
    // Friction 30 -> 20 deg and dilation 10 -> 0 deg over 1e-5: the step
    // ends inside the segment with strength and flow of the stored history
    // (option 1) or of the start of the step (option 0).
    p.pls0={0}; p.pls1={1e-5};
    set_mc_material(f, 10);
    p.friction_angle1={20}; p.dilation_angle1={0};
    f.set_material();
    const bool start_angles=return_mapping_option==0;
    check_law_return(f,named("friction and dilation soften within the step"),face.data(),0,0,
                     start_angles);
    if(!pure_2d) {
        check_law_return(f,"friction and dilation soften within the step, edge",edge.data(),0,0,
                         start_angles);
        // Tension past the cap 0 with an endpoint friction angle 1e-8 deg away
        // from constant: the two returns differ by O(1e-8), not by a jump.
        for(int type=0;type<2;++type) {
            const double phi=type==0 ? 30 : 10, cohesion=type==0 ? 1e5 : 1e4;
            const double trial[3]={-1e6,8e5,9e5};
            double returned[2][3];
            for(int perturbed=0;perturbed<2;++perturbed) {
                p.pls0={0}; p.pls1={1e-3};
                set_mc_material(f, 0);
                p.cohesion0=p.cohesion1={cohesion};
                p.friction_angle0={phi}; p.friction_angle1={phi-(perturbed ? 1e-8 : 0)};
                p.tension_max=0;
                f.set_material();
                check_law_return(f,"tension with a constant or nearly constant friction angle",
                                 trial,0,0,perturbed && start_angles);
                returned[perturbed][0]=f.m->stress[0][0];
                returned[perturbed][1]=f.m->stress[0][1];
                returned[perturbed][2]=NDIMS==2 ? (*f.yy)[0] : f.m->stress[0][2];
            }
            for(int k=0;k<3;++k)
                check("tension return continuous in the friction angle",
                      std::fabs(returned[0][k]-returned[1][k])<=1e-6*cohesion,1);
        }
    }
    else {
        // Tensile cap at zero: tension, the face-cap corner and the cap-cap corner.
        for(int dilation : {0,30}) {
            p.pls0={1}; p.pls1={2};
            set_mc_material(f, dilation);
            p.tension_max=0;
            f.set_material();
            const double tension[][3]={{-1e5,3e5,0},{-1e6,3e5,0},{-3e5,3e5,0},{2e5,3e5,0}};
            for(const auto& trial : tension) check_law_return(f,"pure 2D tensile return",trial);
            check("pure 2D biaxial tension returns to the cap corner",
                  std::fabs(f.m->stress[0][0])<=1e-9*3e5 && std::fabs(f.m->stress[0][1])<=1e-9*3e5,1);
        }
    }
#ifdef THREED
    // Tensile cap at zero: biaxial and triaxial tension and tension-shear corners.
    p.pls0={1}; p.pls1={2};
    set_mc_material(f, 0);
    p.tension_max=0;
    f.set_material();
    const double tension[][3]={{0,2e5,3e5},{1e5,2.9e5,3e5},{-1e6,2e5,3e5},
                               {-1e6,-1e6,3e5},{2e5,2.5e5,3e5}};
    for(const auto& trial : tension) check_law_return(f,"tensile corner",trial);
    for(int d=0;d<3;++d) check("triaxial tension returns to the cap apex",f.m->stress[0][d],0);
#endif
    // The default cap c/tan(phi) meets the face (pure 2D) or passes through the
    // cone apex (3D) at c*sqrt(3) in every principal stress.
    if(pure_2d || NDIMS==3) {
        set_mc_material(f, 0);
        const Trial hydrostatic={3e5,3.2e5,pure_2d ? 0 : 3.5e5};
        check_law_return(f,pure_2d ? "pure 2D face-cap apex" : "cone apex",hydrostatic.data());
        for(int d=0;d<(pure_2d ? 2 : 3);++d)
            check(pure_2d ? "pure 2D apex stress" : "cone apex stress",f.m->stress[0][d],
                  1e5*std::sqrt(3.));
    }
}

// Rate-and-state returns evaluate the same piecewise cohesion law, with rate
// and state held at this step's values.
static void test_rsf_law_returns() {
    std::printf("[Rate-and-state Mohr-Coulomb returns over the strength law]\n");
    MaterialFixture f(MatProps::rh_ep_rsf, 0);
    auto& p=f.m->param.mat;
    const double edge[3]={-1e6,-1e6,-1e5};
    // With a = b = 0 and state model 0 the friction is the static 30 deg, so
    // cohesion 1e5 -> 9e4 over 1e-8 must reach the 90 kPa residual as in EP.
    p.pls0={0}; p.pls1={1e-8};
    set_mc_material(f, 0);
    p.cohesion1={9e4};
    p.direct_a={0}; p.evolution_b={0};
    p.characteristic_velocity={1e-6}; p.characteristic_distance={0.01};
    p.state_var_model=0;
    f.set_material();
    check_law_return(f,"rate-and-state softening past the segment",edge,1e-8);
    // Aging state law with a != b: a plastic update advances the state once,
    // exactly as an elastic update from the same state does.
    p.direct_a={0.01}; p.evolution_b={0.015}; p.state_var_model=1;
    f.set_material();
    const double elastic[3]={-1e6,-9e5,-8e5};
    double state[2];
    for(int plastic=0;plastic<2;++plastic) {
        f.set_principal(plastic ? edge : elastic);
        f.state->assign(f.m->var.nelem,10.);
        f.update();
        check(plastic ? "rate-and-state plastic update" : "rate-and-state elastic update",
              (*f.dpls)[0]>0,plastic);
        state[plastic]=(*f.state)[0];
    }
    check("rate-and-state state advanced once per step",state[1],state[0]);
    // Friction 30 -> 20 deg (within the step and past pls1): the returned
    // stress and the stored friction coefficient follow the final history
    // (option 1) or the start of the step (option 0). With a = b = 0 the
    // coefficient is tan(phi(pls)); with the aging law the stress must lie on
    // the face that the stored coefficient defines.
    const bool start_angles=return_mapping_option==0;
    const double face[3]={-1e6,-5e5,-1e5};
    for(double end : {1e-5,1e-8}) for(int aging=0;aging<2;++aging) {
        p.pls0={0}; p.pls1={end};
        set_mc_material(f, 0);
        p.friction_angle1={20};
        p.direct_a={aging ? 0.01 : 0.}; p.evolution_b={aging ? 0.015 : 0.};
        p.state_var_model=aging;
        f.set_material();
        if(!aging) check_law_return(f,"rate-and-state friction softening",face,0,0,start_angles);
        // One update: the next one would recompute the coefficient anyway.
        f.set_principal(face);
        f.friction->assign(f.m->var.nelem,0.6);
        f.state->assign(f.m->var.nelem,10.);
        f.update();
        const double history=(*f.pls)[0], mu=(*f.friction)[0];
        const MatProps::StrengthPiece piece=element_strength(f,start_angles ? 0 : history);
        const double amc=piece.strength.amc, anphi=piece.strength.anphi;
        const double cohesion=amc/(2*std::sqrt(anphi));
        const double static_mu=std::tan(std::asin((anphi-1)/(anphi+1)));
        check("rate-and-state friction softening yields",history>0,1);
        if(!aging) check(start_angles ? "stored friction coefficient of the start of the step"
                                      : "stored friction coefficient of the final history",
                         mu,static_mu);
        else check("stored friction coefficient differs from the static one",
                   std::fabs(mu-static_mu)>1e-3,1);
        const double sine=std::sin(std::atan(mu)), nphi=(1+sine)/(1-sine);
        auto& stress=f.m->stress;
        const double high=NDIMS==2 ? (*f.yy)[0] : stress[0][2];
        check("stress on the face of the stored friction coefficient",
              (stress[0][0]-nphi*high+2*cohesion*std::sqrt(nphi))/(2*cohesion),0);
    }
}

// Laws that jump up where the return would land, so that no increment records
// itself: the run stops with status 54 and that cause.
static void test_rejected_jump(bool cap) {
    // Option 1 for the cap: the fixture reads it when it is built.
    return_mapping_option=cap ? 1 : 0;
    MaterialFixture f(MatProps::rh_ep, 0);
    set_mc_material(f, 0);
    auto& p=f.m->param.mat;
    if(cap) {
        // No cohesion and friction 30 -> 0 deg over 1.5e-6, with the end-of-step
        // angles: the tensile cap opens from 0 to 1e5 Pa where the friction reaches 0.
        // Hydrostatic tension 2e5 Pa returns past that on the closed cap and short of
        // it on the open one. In 2D, pure 2D: plane strain measures plastic strain by
        // its largest shear, which this flow leaves at 0.
        p.is_plane_strain=false;
        p.pls0={0}; p.pls1={1.5e-6};
        p.cohesion0=p.cohesion1={0};
        p.friction_angle0={30}; p.friction_angle1={0};
        p.tension_max=1e5;
        f.set_material();
        const double trial[3]={2e5,2e5,2e5};
        f.set_principal(trial);
    }
    else {
        // Cohesion rises 1e4 -> 1e5 at 1e-5: the weaker strength returns past the rise
        // (1.26e-5) and the stronger short of it (7.4e-6).
        p.pls0=p.pls1={1e-5};
        p.cohesion0={1e4}; p.cohesion1={1e5};
        f.set_material();
        const double trial[3]={-1.2e6,-2e5,-1e5};
        f.set_principal(trial);
    }
    f.update();
}

// A non-finite trial has no admissible return. update_stress must stop the
// process with status 54; check.sh runs this case separately, in 2D both in
// plane strain and in pure 2D.
static void test_rejected_mc_return(bool pure_2d) {
    MaterialFixture f(MatProps::rh_ep, 0);
    if (pure_2d) f.m->param.mat.is_plane_strain=false;
    set_mc_material(f, 0);
    const double trial[3]={std::nan(""),-1e6,-1e5};
    f.set_principal(trial);
    f.update();
}

#ifdef THREED
// Two equal principal stresses split only by shears near underflow, the rounding a
// confined column leaves, send dsyevh3 to its QL fallback with a first row far below
// the diagonal. The return must match the same state without them.
static void test_tiny_shear_return() {
    std::printf("[Mohr-Coulomb return with near-underflow shear]\n");
    MaterialFixture f(MatProps::rh_ep, 0);
    set_mc_material(f, 0);
    const double trial[3]={-1e6,-1e6,-1e5};
    f.set_principal(trial);
    f.update();
    double reference[NSTR];
    for(int d=0;d<NSTR;++d) reference[d]=f.m->stress[0][d];
    f.set_principal(trial);
    for(int e=0;e<f.m->var.nelem;++e) {
        f.m->stress[e][4]=-4.9e-88;
        f.m->stress[e][5]=-6e-16;
    }
    f.update();
    for(int d=0;d<NSTR;++d)
        check("near-underflow shear leaves the return unchanged",f.m->stress[0][d],reference[d]);
}
#endif

// dsyevh3 on one symmetric matrix scaled from 1e-300 to 1e30: the eigenvalues scale
// with it and Q diag(w) Q^T rebuilds it. Below ~1e-7 the eigenvalues come from the QL
// fallback, whose first Householder reflection is built from the small first row.
static void test_eigen_scales() {
    std::printf("[3x3 eigensystem across input scales]\n");
    const double B[3][3]={{1,0.5,2},{0.5,3,-0.7},{2,-0.7,2}};
    auto solve=[&](double scale,double w[3],double& rebuilt) {
        double A[3][3],Q[3][3];
        for(int i=0;i<3;++i) for(int j=0;j<3;++j) A[i][j]=scale*B[i][j];
        // A NaN would drop out of std::max and break the sort below: check
        // finiteness first and stop there.
        bool solved=dsyevh3(A,Q,w)==0;
        for(int i=0;i<3;++i) {
            solved=solved && std::isfinite(w[i]);
            for(int j=0;j<3;++j) solved=solved && std::isfinite(Q[i][j]);
        }
        rebuilt=0;
        if(!solved) return false;
        for(int i=0;i<3;++i) for(int j=0;j<3;++j) {
            double a=0;
            for(int m=0;m<3;++m) a+=Q[i][m]*w[m]*Q[j][m];
            rebuilt=std::max(rebuilt,std::fabs(a-A[i][j])/(3*scale));
        }
        std::sort(w,w+3);
        return true;
    };
    double unit[3],rebuilt;
    check("unit-scale eigensystem is finite",solve(1,unit,rebuilt),1);
    for(int p=-300;p<=30;p+=10) {
        const double scale=std::pow(10.0,p);
        double w[3];
        const bool solved=solve(scale,w,rebuilt);
        double eigen=0;
        for(int i=0;i<3;++i) eigen=std::max(eigen,std::fabs(w[i]/scale-unit[i])/3);
        const bool ok=solved && eigen<1e-12 && rebuilt<1e-12;
        if(!ok) std::printf("  scale 1e%d: %s, eigenvalue error %.3g, rebuilt error %.3g\n",
                            p,solved?"finite":"failed or not finite",eigen,rebuilt);
        check("eigensystem scales with the matrix",ok,1);
    }
}

static void test_pt_norm() {
    std::printf("[PT projected residual reduction]\n");
    auto* m = new Fixture();
    build_mesh(*m);
    auto* constraints = new VelocityConstraints(m->var.nnode);
    for(int n=0;n<m->var.nnode;++n) constraints->free_rank[n]=NDIMS;
    for (int first : {0, m->var.nnode-2})
    for(double scale : {1.,1e200,1e-200}) {
        for(int n=0;n<m->var.nnode;++n)
            for(int d=0;d<NDIMS;++d) m->force[n][d]=0;
        m->force[first][0]=3*scale;
        m->force[first+1][NDIMS-1]=-4*scale;
        const double expected=scale*std::sqrt(25./(m->var.nnode*NDIMS));
        check("scaled residual ratio",pt_residual_rms(m->var,*constraints,m->force)/expected,1);
    }
    for(int n=0;n<m->var.nnode;++n) constraints->free_rank[n]=0;
    check("no free residual",pt_residual_rms(m->var,*constraints,m->force),0);
    m->force[0][0]=std::numeric_limits<double>::quiet_NaN();
    check("nonfinite force rejected",std::isnan(pt_residual_rms(m->var,*constraints,m->force)),1);
    m->force[0][0]=0;
    m->vel[0][0]=std::numeric_limits<double>::infinity();
    check("nonfinite velocity rejected",std::isnan(pt_residual_rms(m->var,*constraints,m->force)),1);
    delete constraints; delete m;
}

static void test_darcy_gravity_limit() {
    std::printf("[Darcy zero-gravity and hydrostatic invariants]\n");
    auto* m = new Fixture();
    build_mesh(*m);
    m->param.control.has_PT = true; // this fixture prescribes total-stress increments
    const int nn=m->var.nnode, ne=m->var.nelem;
    bind_fields(*m,273.0);
    auto* mass=new double_vec(nn,1.0);
    auto* rate=new double_vec(nn,0.0);
    auto* old_mean=new double_vec(ne,0.0);
    auto* flags=new uint_vec(nn,0);
    m->var.hmass=mass; m->var.bcflag=flags; m->var.dt=1;
    for(int i=0;i<nbdrytypes_hydro;++i) m->var.hbc_types[i]=0;
    auto& p=m->param.mat;
    p.nmat=1; p.rheol_type=MatProps::rh_elastic;
    p.rho0={2500}; p.alpha={0}; p.porosity={0.2};
    p.bulk_modulus={1e9}; p.shear_modulus={1e9};
    p.visc_exponent={1}; p.visc_coefficient={1};
    p.biot_coeff={1}; p.heat_capacity={1}; p.therm_cond={1};
    p.fluid_bulk_modulus={2e9}; p.hydraulic_perm={1e-12}; p.bulk_modulus_s={37e9};
    p.fluid_visc={1e-3}; p.fluid_rho0={1200}; p.fluid_alpha={0};
    reset_material(*m);
    double_vec reference(nn,0.0);
    double reference_diffusivity=0;
    for(int driven=0;driven<2;++driven) for(double gravity : {0.0,10.0}) {
        m->param.control.gravity=gravity;
        for(int n=0;n<nn;++n) {
            m->pressure[n]=2e5+driven*1e5*m->coord[n][0]-1200*gravity*m->coord[n][NDIMS-1];
            m->pressure_increment[n]=0;
        }
        update_pore_pressure(m->param,m->var,m->pressure,m->pressure_increment,
                             *rate,m->tmp_result,m->stress,*old_mean);
        #pragma acc wait
        double net=0,work=0;
        for(int n=0;n<nn;++n) {
            check("finite Darcy pressure",std::isfinite(m->pressure[n]),1);
            check("finite Darcy increment",std::isfinite(m->pressure_increment[n]),1);
            net+=(*mass)[n]*m->pressure_increment[n];
            work+=driven*1e5*m->coord[n][0]*m->pressure_increment[n];
            if(!driven) check("constant potential stationary",m->pressure_increment[n],0);
            if(gravity==0) reference[n]=m->pressure_increment[n];
            else check("hydrostatic background leaves flux unchanged",m->pressure_increment[n],reference[n]);
        }
        check("closed-system fluid conservation",net,0);
        if(driven) check("nonconstant pressure diffuses",work>0,1);
        const double diffusivity=m->var.mat->hydro_diff_max;
        check("finite positive hydraulic diffusivity",std::isfinite(diffusivity)&&diffusivity>0,1);
        if(gravity==0) reference_diffusivity=diffusivity;
        else check("gravity independent hydraulic diffusivity",diffusivity,reference_diffusivity);
    }
    // Prescribe an undrained elastic compression, independently enforcing
    // zeta = alpha*eps_v + p/M = 0 and sigma = C:eps - alpha*p*I.
    // The production pressure update must recover this pressure from stress.
    m->param.control.gravity=0;
    m->param.control.has_hydraulic_diffusion=true;
    for (double alpha : {0.2,0.6,1.0}) for (double grain_bulk : {2.5e9,1e30}) {
        p.biot_coeff[0]=alpha; p.bulk_modulus_s[0]=grain_bulk;
        reset_material(*m);
        const double inv_M=0.2/2e9+(alpha-0.2)/grain_bulk;
        const double eps_v=-1e-5;
        const double pressure=-alpha*eps_v/inv_M;
        const double mean_modulus=1e9+(NDIMS==2 ? 1e9/3 : 0);
        const double mean_stress=mean_modulus*eps_v-alpha*pressure;
        std::fill(mass->begin(),mass->end(),0);
        for(int e=0;e<ne;++e) {
            for(int i=0;i<NODES_PER_ELEM;++i)
                (*mass)[m->connectivity[e][i]]+=hydraulic_storage(*m->var.mat,e)*m->volume[e]/NODES_PER_ELEM;
            for(int s=0;s<NSTR;++s) m->stress[e][s]=(s<NDIMS ? mean_stress : 0);
        }
        for(double dt : {0.1,1.0}) {
            m->var.dt=dt;
            std::fill(m->pressure.begin(),m->pressure.end(),0);
            update_pore_pressure(m->param,m->var,m->pressure,m->pressure_increment,
                                 *rate,m->tmp_result,m->stress,*old_mean);
            for(int n=0;n<nn;++n) {
                check("undrained pressure from accepted stress",m->pressure[n]/pressure,1);
                check("undrained fluid content",(alpha*eps_v+inv_M*m->pressure[n])/std::abs(alpha*eps_v),0);
                check("pending increment is negative pressure change",m->pressure_increment[n]/pressure,-1);
            }
        }
        m->var.mat->hydro_diff_max=1e30;
        refresh_hydraulic_diffusivity(m->param,m->var);
        check("diffusivity seeded before pressure update",m->var.mat->hydro_diff_max/reference_diffusivity>0,1);
        check("refresh matches accepted transport bound",m->var.mat->hydro_diff_max,
              1e-12/1e-3/(inv_M+alpha*alpha/mean_modulus));
        m->param.control.gravity=10;
        m->param.ic.excess_pore_pressure=1e5;
        initial_hydrostatic_state(m->param,m->var,m->pressure,m->pressure_increment);
        const double H=1e9+4e9/3;
        const double initial_pressure=alpha*1e5/(H*inv_M+alpha*alpha);
        for(int n=0;n<nn;++n)
            check("uniaxial undrained initial loading",(m->pressure[n]+1200*10*m->coord[n][NDIMS-1])/initial_pressure,1);
        m->param.control.gravity=0;
    }
    {
        p.alpha={3e-5}; p.fluid_alpha={2e-4}; p.biot_coeff={0.6};
        reset_material(*m);
        auto& old_temperature=m->old_temperature;
        old_temperature.assign(nn,273);
        m->var.old_temperature=&old_temperature;
        m->temperature.assign(nn,283); m->pressure.assign(nn,0);
        std::fill(mass->begin(),mass->end(),0.0);
        std::fill(m->stress.data(),m->stress.data()+ne*NSTR,0.0);
        const double storage=hydraulic_storage(*m->var.mat,0);
        for(int e=0;e<ne;++e) for(int j=0;j<NODES_PER_ELEM;++j)
            (*mass)[m->connectivity[e][j]]+=storage*m->volume[e]/NODES_PER_ELEM;
        m->var.dt=0.1; m->param.control.gravity=0;
        update_pore_pressure(m->param,m->var,m->pressure,m->pressure_increment,
                             *rate,m->tmp_result,m->stress,*old_mean);
        const double beta=0.2*2e-4+(0.6-0.2)*3e-5;
        const double thermal=beta-0.6e9*3e-5*hydraulic_mean_stress_compliance(*m->var.mat,0);
        for(int n=0;n<nn;++n)
            check("sealed fixed-stress thermal pressure",m->pressure[n],thermal*10/storage);
        auto& balance=m->balance; balance.assign(5,0); m->var.fluid_balance=&balance;
        m->pressure.assign(nn,0);
        std::fill(flags->begin(),flags->end(),BOUNDZ1); m->var.hbc_types[iboundz1]=1;
        update_pore_pressure(m->param,m->var,m->pressure,m->pressure_increment,
                             *rate,m->tmp_result,m->stress,*old_mean);
        double volume=0;
        for(int e=0;e<ne;++e) volume+=m->volume[e];
        check("drained thermal pressure remains prescribed",m->pressure[0],0);
        check("drained thermal outflow",balance[3],-thermal*10*volume);
        m->var.fluid_balance=nullptr;
        std::fill(flags->begin(),flags->end(),0); m->var.hbc_types[iboundz1]=0;
        m->var.old_temperature=nullptr;
        p.alpha={0}; p.fluid_alpha={0};
    }
    std::fill(mass->begin(),mass->end(),1.0);
    for (double angle : {0.0, 0.37, 1.1}) {
        const double c=std::cos(angle), s=std::sin(angle);
        p.permeability_tensor={1+2*c*c, 1, 1+2*s*s, 0, 2*c*s, 0};
        reset_material(*m);
        double grad[NDIMS];
        for(int d=0;d<NDIMS;++d) grad[d]=1e5*(d+1);
        double dissipation=0, work=0, net=0;
        for(int n=0;n<nn;++n) {
            m->pressure[n]=0;
            for(int d=0;d<NDIMS;++d) m->pressure[n]+=grad[d]*m->coord[n][d];
        }
        for(int e=0;e<ne;++e)
            for(int i=0;i<NDIMS;++i) for(int j=0;j<NDIMS;++j)
                dissipation+=m->volume[e]*1e-9*grad[i]*m->var.mat->perm_tensor[i][j]*grad[j];
        reference=m->pressure;
        std::fill(m->stress.data(),m->stress.data()+ne*NSTR,0.0);
        m->param.control.gravity=0; m->var.dt=0.1;
        update_pore_pressure(m->param,m->var,m->pressure,m->pressure_increment,
                             *rate,m->tmp_result,m->stress,*old_mean);
        for(int n=0;n<nn;++n) {
            net+=m->pressure_increment[n];
            work+=reference[n]*m->pressure_increment[n];
        }
        check("rotated anisotropic Darcy dissipation",work/m->var.dt,dissipation);
        check("rotated tensor conserves fluid",net,0);
    }
    p.permeability_tensor.clear();
    reset_material(*m);
    {
        p.hydraulic_perm={0}; p.matrix_storage={1}; p.matrix_exchange={0.1};
        reset_material(*m);
        auto& matrix=m->matrix_pressure; matrix.assign(ne,0); m->var.matrix_pressure=&matrix;
        std::fill(mass->begin(),mass->end(),0.0);
        for(int e=0;e<ne;++e) for(int j=0;j<NODES_PER_ELEM;++j)
            (*mass)[m->connectivity[e][j]]+=0.5*m->volume[e]/NODES_PER_ELEM;
        double previous=0;
        for(double dt : {0.02,0.01}) {
            std::fill(matrix.begin(),matrix.end(),0);
            m->pressure.assign(nn,1); m->var.dt=dt;
            for(int k=0;k<int(10/dt);++k)
                update_pore_pressure(m->param,m->var,m->pressure,m->pressure_increment,
                                     *rate,m->tmp_result,m->stress,*old_mean);
            const double exact=1.0/3+2.0/3*std::exp(-3.0);
            const double error=std::abs(m->pressure[0]-exact);
            check("exchange first-order time convergence",previous==0 || error<0.51*previous,1);
            previous=error;
            for(int n=0;n<nn;++n) check("exchange homogeneous transport pressure",m->pressure[n],m->pressure[0]);
            for(int e=0;e<ne;++e) check("exchange conserves total reservoir fluid",0.5*m->pressure[0]+matrix[e],0.5);
        }
        m->pressure.assign(nn,1); std::fill(matrix.begin(),matrix.end(),0); m->var.dt=50;
        update_pore_pressure(m->param,m->var,m->pressure,m->pressure_increment,
                             *rate,m->tmp_result,m->stress,*old_mean);
        check("large fixed step uses stable explicit exchange substeps",m->pressure[0]>0 && m->pressure[0]<1,1);
        for(int n=0;n<nn;++n)
            check("pending pressure increment spans every substep",m->pressure_increment[n],1-m->pressure[n]);
        for(int e=0;e<ne;++e) check("subcycled exchange fluid balance",0.5*m->pressure[0]+matrix[e],0.5);
        m->var.matrix_pressure=nullptr;
        p.matrix_exchange.clear(); p.matrix_storage.clear(); p.hydraulic_perm={1e-12};
    }
    {
        const bool feedback=m->param.control.has_poroelastic_pressure_feedback;
        m->param.control.has_poroelastic_pressure_feedback=true;
        p.biot_coeff={.6}; p.matrix_storage={1e-9};
        auto& strain=m->strain; strain.resize(ne,0.0);
        auto& matrix=m->matrix_pressure; matrix.assign(ne,2e5);
        m->var.strain=&strain; m->var.matrix_pressure=&matrix;
        m->pressure.assign(nn,1e5); m->pressure_increment.assign(nn,-100);
        for(int e=0;e<ne;++e) strain[e][0]=.001;
        reset_material(*m);
        const double storage=hydraulic_inverse_biot_modulus(*m->var.mat,0);
        ContentRemap remap(m->param,m->var);
        remap.transfer(m->var,m->coord,m->connectivity);
        m->pressure.assign(nn,2e5); m->pressure_increment.assign(nn,-500);
        matrix.assign(ne,3e5);
        for(int e=0;e<ne;++e) {
            strain[e][0]=.002;
            for(int d=0;d<NSTR;++d) m->stress[e][d]=0;
        }
        remap.restore(m->param,m->var);
        for(int n=0;n<nn;++n) {
            check("remap physical content independent of strain interpolation",
                  (storage*(m->pressure[n]-1e5)+.6*.001)/(.6*.001),0);
            check("remap numerical buffer independent of physical water",m->pressure_increment[n],-100);
        }
        for(int e=0;e<ne;++e) {
            check("remap preserves immobile reservoir separately",matrix[e],2e5);
            const int n=m->connectivity[e][0];
            check("remap retains mechanical effective-stress reference",
                  (m->stress[e][0]+.6*(m->pressure[n]+m->pressure_increment[n]))/1e5,
                  .6*(2e5-500)/1e5);
        }
        ContentRemap fixed(m->param,m->var);
        fixed.transfer(m->var,m->coord,m->connectivity);
        const double fixed_pressure=m->pressure[0];
        for(int e=0;e<ne;++e) strain[e][0]+=.0001;
        std::fill(flags->begin(),flags->end(),BOUNDZ1); m->var.hbc_types[iboundz1]=1;
        fixed.restore(m->param,m->var);
        for(int n=0;n<nn;++n)
            check("all-Dirichlet remap preserves imposed pressure",m->pressure[n],fixed_pressure);
        std::fill(flags->begin(),flags->end(),0); m->var.hbc_types[iboundz1]=0;
        m->var.strain=nullptr; m->var.matrix_pressure=nullptr;
        p.matrix_storage.clear();
        m->param.control.has_poroelastic_pressure_feedback=feedback;
    }
    p.biot_coeff[0]=0.6; p.bulk_modulus_s[0]=2.5e9;
    p.permeability_model=1;
    m->var.pore_strain=new double_vec(ne,0);
    m->var.etmp=new double_vec(ne);
    m->strain_rate.resize(ne,0.0);
    m->var.strain_rate=&m->strain_rate;
    reset_material(*m);
    double_vec phi(ne,0.2);
    for(int step=0;step<10;++step) {
        m->var.dt=0.1;
        std::fill(m->pressure_increment.begin(),m->pressure_increment.end(),-1000*(step+1));
        for(int e=0;e<ne;++e) {
            m->strain_rate[e][0]=(e+1)*1e-3;
            const double deps=0.1*(e+1)*1e-3+1000*(step+1)/2.5e9;
            // Independent RK4 integration of dphi/deps = alpha - phi.
            for(int j=0;j<100;++j) {
                const double h=deps/100;
                const double a=0.6-phi[e], b=0.6-(phi[e]+h*a/2);
                const double c=0.6-(phi[e]+h*b/2), d=0.6-(phi[e]+h*c);
                phi[e]+=h*(a+2*b+2*c+d)/6;
            }
        }
        update_porosity(m->var);
        for(int e=0;e<ne;++e) {
            const double actual=m->var.mat->phi(e);
            check("Biot porosity against integrated pore-volume law",actual,phi[e]);
            const double resistance=m->var.mat->perm(e)/1e-12*(1-actual)*(1-actual)/(actual*actual*actual);
            check("Kozeny-Carman pore resistance",resistance,0.8*0.8/(0.2*0.2*0.2));
        }
    }
    p.permeability_model=2;
    for (double exponent : {0.5, 3.0, 5.0}) {
        p.permeability_exponent={exponent};
        reset_material(*m);
        double previous=1e-12;
        for (double eps : {0.0, 0.05, 0.1, 0.05, 0.0, -0.1}) {
            std::fill(m->var.pore_strain->begin(),m->var.pore_strain->end(),eps);
            const double phi=0.6-0.4*std::exp(-eps);
            const double k=m->var.mat->perm(0);
            if (eps == 0)
                check("power law returns to reference permeability",k/1e-12,1);
            else
                check("power law logarithmic sensitivity",std::log(k/1e-12)/std::log(phi/0.2),exponent);
            if (eps < 0) check("power law compaction reduces permeability",k < previous,1);
            previous=k;
        }
    }
    p.hydraulic_perm={0};
    p.permeability_exponent={1e300};
    reset_material(*m);
    std::fill(m->var.pore_strain->begin(),m->var.pore_strain->end(),0.1);
    check("impermeable power law avoids zero times overflow",m->var.mat->perm(0),0);
    delete m->var.pore_strain; delete m->var.etmp;
    delete m->var.mat; delete flags; delete old_mean; delete rate; delete mass; delete m;
}

static void test_thermal_subcycling() {
    std::printf("[Explicit heat eigenmode, conservation and capacity changes]\n");
    auto* m = new Fixture();
    build_mesh(*m);
    const int nn = m->var.nnode, ne = m->var.nelem;
    bind_fields(*m, 273);
    auto& source = m->heat_source;
    auto& capacity = m->thermal_mass;
    auto& old_temperature = m->old_temperature;
    auto& flags = m->boundary_flags;
    source.assign(ne, 0); capacity.resize(nn);
    old_temperature.assign(nn, 273); flags.assign(nn, 0);
    m->var.radiogenic_source = &source;
    m->var.tmass = &capacity;
    m->var.bcflag = &flags;
    m->var.old_temperature = &old_temperature;
    auto& p = m->param.mat;
    p.nmat=1; p.rheol_type=MatProps::rh_elastic;
    p.rho0={1}; p.alpha={0}; p.porosity={0.2};
    p.bulk_modulus={1e9}; p.shear_modulus={1e9};
    p.visc_exponent={1}; p.visc_coefficient={1};
    p.biot_coeff={1}; p.heat_capacity={1}; p.therm_cond={1};
    p.fluid_bulk_modulus={2e9}; p.hydraulic_perm={1e-12}; p.bulk_modulus_s={37e9};
    p.fluid_visc={1e-3}; p.fluid_rho0={1};
    m->param.control.has_thermal_diffusion = true;
    reset_material(*m);
    for (int n=0; n<nn; ++n)
        capacity[n] = m->volume[0] * (n == 0 ? 2 : 1) / NODES_PER_ELEM;
    // On two reflected unit simplices (-NDIMS,1,...,1) is an exact
    // lumped heat eigenvector with eigenvalue (NDIMS+1)^2.
    const double lambda = NODES_PER_ELEM * NODES_PER_ELEM;
    const double bound = 1.0 / (2 * NDIMS * NODES_PER_ELEM);
    m->var.dt = 3.5 * bound;
    for (int n=0; n<nn; ++n) m->temperature[n] += n == 0 ? -NDIMS : 1;
    update_temperature(m->param, m->var, m->temperature, m->tmp_result);
    #pragma acc wait
    check("heat row bound", m->var.mat->thermal_step_limit, bound);
    const double decay = std::pow(1 - lambda * m->var.dt / 4, 4);
    double heat = 0;
    for (int n=0; n<nn; ++n) {
        check("subcycled heat eigenmode", m->temperature[n], 273 + decay*(n == 0 ? -NDIMS : 1));
        check("thermal pressurization baseline untouched", old_temperature[n], 273);
        heat += capacity[n] * (m->temperature[n] - 273);
    }
    check("insulated heat conservation", heat, 0);
    double previous_error = 0;
    for (int steps : {10, 20, 40}) {
        m->var.dt = 0.01 / steps;
        for (int n=0; n<nn; ++n) m->temperature[n] = 273 + (n == 0 ? -NDIMS : 1);
        for (int i=0; i<steps; ++i)
            update_temperature(m->param, m->var, m->temperature, m->tmp_result);
        #pragma acc wait
        const double error = std::abs(m->temperature[0] - (273 - NDIMS*std::exp(-lambda*0.01)));
        if (previous_error > 0) check("heat first-order temporal convergence", previous_error/error > 1.9 && previous_error/error < 2.1, 1);
        previous_error = error;
    }
    // The same invalidation used after remesh/material refresh must rebuild
    // the capacity-dependent bound; it must not be serialized as history.
    double_vec volume_n(nn), mass(nn), hmass(nn), ymass(nn);
    compute_mass(m->param, m->var, 0.01, volume_n, mass, capacity, hmass, ymass, m->tmp_result);
    #pragma acc wait
    check("heat cache invalidated by mass assembly", m->var.mat->thermal_step_limit, 0);
    for (auto& value : capacity) value *= 2;
    flags[0] = BOUNDZ1;
    m->param.bc.surface_temperature = 280;
    m->var.dt = 3.5 * bound;
    update_temperature(m->param, m->var, m->temperature, m->tmp_result);
    #pragma acc wait
    check("heat bound follows changed capacity", m->var.mat->thermal_step_limit, 2*bound);
    check("subcycled prescribed surface temperature", m->temperature[0], 280);
    delete m->var.mat;
    delete m;
}

static void test_heat_remap() {
    std::printf("[Pending heat follows retained material]\n");
    for (int mode = 0; mode < 3; ++mode) {
        auto *m = new Fixture(); build_mesh(*m);
        m->param.control.has_hydraulic_diffusion = false;
        const array_t old_coordinate(m->coord);
        const conn_t old_connectivity(m->connectivity);
        double_vec energy(m->var.nelem,0);
        energy[0] = 2; energy[1] = 7;
        m->var.pending_heat = &energy;
        ContentRemap remap(m->param,m->var);
        m->var.nelem = 1;
        m->connectivity.resize(1);
        if (mode == 1)
            for (int n = 0; n < m->var.nnode; ++n)
                for (int d = 0; d < NDIMS; ++d) m->coord[n][d] *= .5;
        if (mode == 2)
            for (int n = 0; n < m->var.nnode; ++n)
                for (int d = 0; d < NDIMS; ++d) m->coord[n][d] += 10;
        remap.transfer(m->var,old_coordinate,old_connectivity);
        remap.restore(m->param,m->var);
        check("pending heat uses current cell count",energy.size(),1);
        const double retained = mode == 0 ? 2 : (mode == 1 ? 2*std::pow(.5,NDIMS) : 0);
        check("removed material does not return its heat",energy[0],retained);
        m->var.pending_heat = nullptr;
        delete m;
    }
}

static void test_fluid_state() {
    std::printf("[Finite-density Biot mass and thermodynamic tangents]\n");
    auto* m = new Fixture(); build_mesh(*m);
    bind_fields(*m,273.);
    auto& p=m->param.mat;
    p.nmat=1; p.rheol_type=MatProps::rh_elastic;
    p.rho0={2650}; p.alpha={3e-5}; p.porosity={.2}; p.biot_coeff={.6};
    p.bulk_modulus={2e7}; p.shear_modulus={1e7};
    p.visc_exponent={1}; p.visc_coefficient={1};
    p.heat_capacity={1000}; p.therm_cond={1}; p.fluid_bulk_modulus={1e8};
    p.bulk_modulus_s={1e8}; p.fluid_rho0={1000}; p.fluid_alpha={2e-4};
    p.hydraulic_perm={0}; p.fluid_visc={1e-3};
    reset_material(*m);
    const MatProps& mat=*m->var.mat;
    for (double eps : {-.01,0.,.01})
        for (double theta : {0.,10.,30.}) {
            // Equal constituent moduli give a closed sealed-cell pressure.
            const double pressure=1e8*(std::log((.2*std::exp(2e-4*theta)+
                .4*std::exp(3e-5*theta))/.6)-eps);
            const double t=273+theta;
            const auto state=fluid_state(mat,0,pressure,t,273.,.2,0.,273.,eps);
            const double j=std::exp(eps), mass=j*state.density*state.porosity;
            check("finite EOS sealed mass",mass/200,1);
            auto mass_at=[&](double pp,double tt,double ee) {
                const auto s=fluid_state(mat,0,pp,tt,273.,.2,0.,273.,ee);
                return std::exp(ee)*s.density*s.porosity;
            };
            const double dp=(mass_at(pressure+1000,t,eps)-mass_at(pressure-1000,t,eps))/2000;
            const double dt=(mass_at(pressure,t+.01,eps)-mass_at(pressure,t-.01,eps))/.02;
            const double de=(mass_at(pressure,t,eps+1e-5)-mass_at(pressure,t,eps-1e-5))/2e-5;
            check("mass pressure tangent",dp/(j*state.pressure_slope),1);
            check("mass temperature tangent",dt/(j*state.temperature_slope),1);
            check("Biot mass volume tangent",de/(j*state.density*.6),1);
            const auto anchored=fluid_state(mat,0,pressure,t,273.,state.porosity,pressure,t);
            check("constitutive anchor exact",anchored.porosity==state.porosity,1);
        }
    delete m->var.mat; delete m;
}

static void test_strain_pressure(bool scaled = false) {
    std::printf("[Strain-driven storage, heating and subcycling: %s mass]\n",
                scaled ? "scaled" : "physical");
    auto* m = new Fixture(); build_mesh(*m);
    const int nn=m->var.nnode, ne=m->var.nelem;
    bind_fields(*m,273.01,5e4,999);
    auto& p=m->param.mat;
    p.nmat=1; p.rheol_type=MatProps::rh_elastic;
    p.rho0={2500}; p.alpha={3e-5}; p.porosity={.2}; p.biot_coeff={.6};
    p.bulk_modulus={1e9}; p.shear_modulus={1e9};
    p.visc_exponent={1}; p.visc_coefficient={1};
    p.heat_capacity={1}; p.therm_cond={1}; p.fluid_bulk_modulus={2e9};
    p.bulk_modulus_s={37e9}; p.fluid_rho0={1000}; p.fluid_alpha={2e-4};
    p.hydraulic_perm={0}; p.fluid_visc={1e-3};
    m->param.control.has_hydraulic_diffusion=true;
    m->param.control.has_thermal_pressurization=true;
    m->param.control.is_quasi_static=scaled;
    m->param.control.inertial_scaling=1e5;
    p.visc_min=1e30;
    if (scaled) p.bulk_modulus={2e9};
    reset_material(*m);
    double_vec mass(nn),volume_n(nn),tmass(nn),ymass(nn),rate(nn),old_mean(ne,-2e5);
    auto& hmass=m->hydraulic_mass;
    auto& old_temperature=m->old_temperature;
    auto& vol_rate=m->volume_rate;
    auto& balance=m->balance;
    auto& flags=m->boundary_flags;
    hmass.resize(nn); old_temperature.assign(nn,273);
    vol_rate.resize(ne); balance.resize(5); flags.assign(nn,0);
    m->var.hmass=&hmass; m->var.bcflag=&flags; m->var.dt=2;
    m->var.edvoldt=&vol_rate; m->var.old_temperature=&old_temperature;
    m->var.fluid_balance=&balance;
    for(int e=0;e<ne;++e) {
        vol_rate[e]=(e+1)*1e-7;
        for(int d=0;d<NDIMS;++d) m->stress[e][d]=1e6;
    }
    compute_mass(m->param,m->var,.01,volume_n,mass,tmass,hmass,ymass,m->tmp_result);
    const double storage=.2/2e9+.4/37e9;
    const double beta=.2*2e-4+.4*3e-5;
    const double heat=beta*(m->temperature[0]-273);
    check("physical pressure capacity",hydraulic_storage(m->param,*m->var.mat,0)/storage,1);
    double expected_total=0;
    double_vec nodal_change(nn,0),nodal_volume(nn,0);
    for(int e=0;e<ne;++e) {
        const double change=(heat-.6*m->var.dt*vol_rate[e])*m->volume[e];
        expected_total+=change;
        for(int j=0;j<NODES_PER_ELEM;++j) {
            const int n=m->connectivity[e][j];
            nodal_change[n]+=change/NODES_PER_ELEM;
            nodal_volume[n]+=m->volume[e]/NODES_PER_ELEM;
        }
    }
    for(int substeps : {1,7}) {
        std::fill(m->pressure.begin(),m->pressure.end(),5e4);
        std::fill(m->pressure_increment.begin(),m->pressure_increment.end(),999);
        m->var.mat->hydro_step_limit=m->var.dt/substeps;
        update_pore_pressure(m->param,m->var,m->pressure,m->pressure_increment,
                             rate,m->tmp_result,m->stress,old_mean);
        #pragma acc wait
        for(int n=0;n<nn;++n) {
            const double dp=nodal_change[n]/(storage*nodal_volume[n]);
            check("physical pressure strain/heating increment",m->pressure[n],5e4+dp);
            check("raw increment reset and integrated once",m->pressure_increment[n],-dp);
        }
        check("physical fluid content",balance[0]/expected_total,1);
        check("physical fluid balance",balance[4]/expected_total,0);
    }
    delete m->var.mat; delete m;
}

static void test_fluid_heat_transport() {
    std::printf("[Darcy heat transport, mixed capacity and bounded substeps]\n");
    auto* m = new Fixture(); build_mesh(*m);
    const int nn=m->var.nnode, ne=m->var.nelem;
    bind_fields(*m,300);
    auto& capacity=m->thermal_mass;
    auto& sources=m->heat_source;
    auto& flags=m->boundary_flags;
    capacity.resize(nn); sources.assign(ne,0); flags.assign(nn,0);
    double_vec volume_n(nn),mass(nn),hmass(nn),ymass(nn);
    m->var.tmass=&capacity; m->var.bcflag=&flags; m->var.radiogenic_source=&sources;
    auto& p=m->param.mat;
    p.nmat=1; p.rheol_type=MatProps::rh_elastic;
    p.rho0={2500}; p.alpha={0}; p.porosity={.2}; p.biot_coeff={1};
    p.bulk_modulus={1e9}; p.shear_modulus={1e9};
    p.visc_exponent={1}; p.visc_coefficient={1};
    p.heat_capacity={1000}; p.therm_cond={0}; p.fluid_bulk_modulus={2e9};
    p.bulk_modulus_s={37e9}; p.fluid_rho0={1200}; p.fluid_heat_capacity={4000};
    p.fluid_alpha={0}; p.hydraulic_perm={1e-12}; p.fluid_visc={1e-3};
    p.permeability_tensor={2,1,1,0,.3,0};
    m->param.control.has_hydraulic_diffusion=true;
    m->param.control.has_fluid_thermal_advection=true;
    reset_material(*m);
    compute_mass(m->param,m->var,.01,volume_n,mass,capacity,hmass,ymass,m->tmp_result);
    #pragma acc wait
    const double cv=.8*2500*1000+.2*1200*4000;
    check("rock/water thermal capacity",m->var.mat->thermal_capacity(0)/cv,1);
    double total_volume=0;
    for(double v:m->volume) total_volume+=v;
    double total_capacity=0;
    for(double v:capacity) total_capacity+=v;
    check("assembled mixed thermal mass",total_capacity/(cv*total_volume),1);
    for(double gravity : {0.,10.}) {
        m->param.control.gravity=gravity;
        for(int n=0;n<nn;++n) {
            const auto x=m->coord[n];
            m->pressure[n]=1e5*x[0]+(2e4-1200*gravity)*x[NDIMS-1];
            m->temperature[n]=300+2*x[0]+3*x[NDIMS-1];
        }
        double expected_q[NDIMS]={};
        double gradient[NDIMS]={};gradient[0]=1e5;gradient[NDIMS-1]=2e4;
        for(int d=0;d<NDIMS;++d) for(int k=0;k<NDIMS;++k)
            expected_q[d]-=1e-9*m->var.mat->perm_tensor[d][k]*gradient[k];
        for(int e=0;e<ne;++e) {
            double q[NDIMS];hydraulic_element_flux(m->param,m->var,e,q);
            for(int d=0;d<NDIMS;++d)
                check("tensor pressure/gravity physical Darcy flux",q[d],expected_q[d]);
        }
        const double_vec before=m->temperature;
        m->var.dt=.25;
        update_temperature(m->param,m->var,m->temperature,m->tmp_result);
        #pragma acc wait
        double energy=0;
        for(int n=0;n<nn;++n) energy+=capacity[n]*(m->temperature[n]-before[n]);
        const double reference=-m->var.dt*1200*4000*(2*expected_q[0]+3*expected_q[NDIMS-1])*total_volume;
        check("integrated affine heat-advection energy",energy/reference,1);
        m->temperature.assign(nn,300);
        m->var.dt=10*m->var.mat->thermal_step_limit;
        update_temperature(m->param,m->var,m->temperature,m->tmp_result);
        #pragma acc wait
        for(double t:m->temperature) check("advection preserves constant temperature",t,300);
        for(int n=0;n<nn;++n) m->temperature[n]=n%2 ? 320:280;
        update_temperature(m->param,m->var,m->temperature,m->tmp_result);
        #pragma acc wait
        for(double t:m->temperature)
            check("advection substeps preserve temperature bounds",t>=280-1e-10&&t<=320+1e-10,1);
    }
    // An inactive mantle-temperature setting cannot drive uniform, static water.
    // Also exercise explicit thermal boundary values that differ from defaults.
    m->param.mesh.zlength=1;
    m->param.bc.surface_temperature=300;
    m->param.bc.mantle_temperature=400;
    m->param.control.gravity=10;
    p.fluid_alpha={.001};
    reset_material(*m);
    for (bool fixed_ends : {false,true}) {
        m->param.bc.tbc[iboundz0]=m->param.bc.tbc[iboundz1]=fixed_ends ? 1:0;
        m->param.bc.tbc_value[iboundz0]=360;
        m->param.bc.tbc_value[iboundz1]=320;
        const double slope=fixed_ends ? 40:0;
        const double top=fixed_ends ? 320:300;
        const double weight=1200*10;
        for (int n=0;n<nn;++n) {
            const double z=m->coord[n][NDIMS-1];
            m->temperature[n]=top-slope*z;
            m->pressure[n]=-weight*(1-.001*(top-300))*z-.5*weight*.001*slope*z*z;
        }
        for (int e=0;e<ne;++e) {
            double q[NDIMS];
            hydraulic_element_flux(m->param,m->var,e,q);
            for (int d=0;d<NDIMS;++d)
                check(fixed_ends ? "active thermal boundary conduction at rest":
                                  "inactive thermal boundary hydrostatic rest",q[d]*1e9,0);
        }
    }

    // At zero flow, enabling advection must retain the ordinary conduction
    // operator, including positive off-diagonals on an obtuse simplex. Hold
    // the mixed thermal mass fixed when comparing the two update paths.
    p.fluid_alpha={0}; p.therm_cond={10};
    m->param.control.gravity=0;
    m->param.control.has_thermal_diffusion=true;
    m->param.bc.tbc[iboundz0]=m->param.bc.tbc[iboundz1]=0;
    reset_material(*m);
    m->pressure.assign(nn,0);
    for (int n=0;n<nn;++n) {
        m->coord[n][0]-=1.5*m->coord[n][NDIMS-1];
        m->temperature[n]=300+2*m->coord[n][0]+3*m->coord[n][NDIMS-1];
    }
    compute_mass(m->param,m->var,.01,volume_n,mass,capacity,hmass,ymass,m->tmp_result);
    const double_vec initial=m->temperature;
    m->var.dt=1;
    m->param.control.has_fluid_thermal_advection=false;
    update_temperature(m->param,m->var,m->temperature,m->tmp_result);
    #pragma acc wait
    const double_vec conduction=m->temperature;
    m->temperature=initial;
    m->param.control.has_fluid_thermal_advection=true;
    update_temperature(m->param,m->var,m->temperature,m->tmp_result);
    #pragma acc wait
    for (int n=0;n<nn;++n)
        check("zero-flow heat operator on obtuse elements",m->temperature[n],conduction[n]);
    delete m->var.mat;delete m;
}

static void test_point_sources() {
    std::printf("[Point-source shape weights, signed schedules and intervals]\n");
    auto* m=new Fixture();build_mesh(*m);
    Injection source;
    source.enabled=true;source.num_points=3;source.rate_model=injection_rate_constant;
    const double c=1.0/NODES_PER_ELEM;
    source.points_x={c,-c,0};source.points_y={c,-c,0};source.points_z={c,-c,0};
    source.rate={5,-2,7};source.total_amount={3,-1.4,2.8};
    source.start_time={.2,.5,.3};source.end_time={.8,1.2,.7};
    double_vec rate(m->var.nnode,0),integrated(m->var.nnode,0);
    PointSourceOwners owners;
    place_point_sources(source,m->var,owners);
    assemble_point_sources(source,owners,m->var,0,1,rate);
    for(int n=0;n<m->var.nnode;++n) {
        const double expected=n==0 ? 2*c+2.8 : n<=NDIMS ? 3*c : -c;
        check("point barycentric rates and shared-node ownership",rate[n],expected);
    }
    const double_vec inflow={1,0,1}, temperature={300,0,400};
    assemble_point_sources(source,owners,m->var,0,1,rate,&inflow,true);
    for(int n=0;n<m->var.nnode;++n)
        check("active inlet rate bounds partial source windows",rate[n],
              n==0 ? 5*c+7 : n<=NDIMS ? 5*c : 0);
    assemble_point_sources(source,owners,m->var,0,1,rate,&temperature,true);
    for(int n=0;n<m->var.nnode;++n)
        check("temperature-weighted active inlet rate",rate[n],
              n==0 ? 1500*c+2800 : n<=NDIMS ? 1500*c : 0);
    source.temperature={300,273,400};
    double_vec heat(m->var.nnode,0),separate_rate(m->var.nnode,0),separate_heat(m->var.nnode,0);
    for (bool peak : {false,true}) {
        assemble_point_sources(source,owners,m->var,0,1,separate_rate,&inflow,peak);
        assemble_point_sources(source,owners,m->var,0,1,separate_heat,&temperature,peak);
        assemble_point_sources(source,owners,m->var,0,1,rate,&inflow,peak,&heat,&source.temperature);
        for(int n=0;n<m->var.nnode;++n) {
            check("shared point map preserves inlet volume",rate[n],separate_rate[n]);
            check("shared point map preserves weighted inlet",heat[n],separate_heat[n]);
        }
    }
    assemble_point_sources(source,owners,m->var,.81,1.4,rate,&inflow,true);
    for(double v:rate) check("inactive inlets do not constrain heat steps",v,0);
    const double breaks[]={0,.13,.47,.81,1.4};
    for(int i=0;i<4;++i) {
        assemble_point_sources(source,owners,m->var,breaks[i],breaks[i+1],rate);
        for(int n=0;n<m->var.nnode;++n) integrated[n]+=rate[n]*(breaks[i+1]-breaks[i]);
    }
    double total=0;for(double v:integrated) total+=v;
    check("subintervals conserve signed scheduled amount",total,4.4);
    source.rate_model=injection_rate_total_amount;
    assemble_point_sources(source,owners,m->var,0,1.4,rate);
    for(int n=0;n<m->var.nnode;++n)
        check("rate/total-amount equivalence",rate[n]*1.4,integrated[n]);
    assemble_point_sources(source,owners,m->var,0,1.4,rate,&inflow,true);
    for(int n=0;n<m->var.nnode;++n)
        check("total-amount inlet has the same active rate bound",rate[n],
              n==0 ? 5*c+7 : n<=NDIMS ? 5*c : 0);
    // The points are fixed in the material: they move with the mesh, away from
    // their configured coordinates.
    for(int n=0;n<m->var.nnode;++n) for(int d=0;d<NDIMS;++d) m->coord[n][d]+=10;
    assemble_point_sources(source,owners,m->var,0,1.4,rate);
    for(int n=0;n<m->var.nnode;++n)
        check("point sources move with the mesh",rate[n]*1.4,integrated[n]);
    source.enabled=false;
    place_point_sources(source,m->var,owners);
    assemble_point_sources(source,owners,m->var,0,1.4,rate);
    for(double v:rate) check("disabled point source clears derived rates",v,0);
    delete m;
}

int main(int argc, char** argv) {
    if(argc==2 && (std::strcmp(argv[1],"--reject-mc-return")==0 ||
                   std::strcmp(argv[1],"--reject-mc-return-pure-2d")==0)) {
        test_rejected_mc_return(NDIMS==2 && std::strcmp(argv[1],"--reject-mc-return-pure-2d")==0);
        return 1;  // check.sh requires status 54 from update_stress
    }
    if(argc==2 && (std::strcmp(argv[1],"--reject-mc-jump")==0 ||
                   std::strcmp(argv[1],"--reject-mc-cap-jump")==0)) {
        test_rejected_jump(std::strcmp(argv[1],"--reject-mc-cap-jump")==0);
        return 1;  // check.sh requires status 54 from update_stress
    }
    test_heat_remap();
    test_fluid_state();
    test_point_sources();
    test_fluid_heat_transport();
    test_strain_pressure();
    test_strain_pressure(true);
    test_thermal_subcycling();
    test_darcy_gravity_limit();
    test_pt_norm();
    test_material_trials();
    for(int option : {0,1}) {
        return_mapping_option=option;
        std::printf("[control.iterative_return_mapping_option = %d]\n",option);
        test_mohr_coulomb_returns();
        test_mohr_coulomb_law_returns(false);
        test_rsf_law_returns();
#ifndef THREED
        test_pure_2d_softening_return();
        test_mohr_coulomb_law_returns(true);
#endif
    }
    return_mapping_option=0;
#ifdef THREED
    test_tiny_shear_return();
#endif
    test_eigen_scales();
    test_constraints();
    test_cancellation();
    test_traction();
    test_platen();
    test_damping();
    test_body_force();
    std::printf("\n%d passed, %d failed\n", npass, nfail);
    return nfail == 0 ? 0 : 1;
}
