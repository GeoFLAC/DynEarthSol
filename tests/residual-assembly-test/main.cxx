// B0 focused test: undamped force exposure in update_force.
// Modes: -DBASE (5-arg API) / -DB0 (6-arg API with force_undamped output).
// Mesh: two elements sharing node 0, identical constant stress -> cancellation at node 0.
// Cases: cancellation, Neumann traction, damping exclusion, body force.

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

static void run_update_force(Fixture& m) {
#ifdef BASE
    update_force(m.param, m.var, m.force, m.force_residual, m.tmp_result);
#else
    update_force(m.param, m.var, m.force, m.force_residual, m.tmp_result, &m.force_undamped);
    // The optional output cannot change either existing output. Reevaluate the
    // identical assembly with the default null argument and compare exact bits.
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

    // Measurement only (legacy defect record): force_residual at the shared node
    // is overwritten per incident element, so it equals the last element's term,
    // not the assembled sum. With identical stress in both elements it is nonzero
    // while the true assembled internal force is exactly 0.
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
    m->temperature.assign(m->var.nnode, 273.0);
    m->pressure.assign(m->var.nnode, 0.0);
    m->pressure_increment.assign(m->var.nnode, 0.0);
    m->strain_rate.resize(m->var.nelem);
    m->markers.assign(m->var.nelem, int_vec(1, 1));
    m->etmp_int.resize(m->var.nelem);
    m->var.temperature = &m->temperature;
    m->var.ppressure = &m->pressure;
    m->var.dppressure = &m->pressure_increment;
    m->var.strain_rate = &m->strain_rate;
    m->var.elemmarkers = &m->markers;
    m->var.etmp_int = &m->etmp_int;
    m->var.log_table = m->var.tan_table = m->var.sin_table = &m->lookup;
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
    m->var.mat = new MatProps(m->param, m->var);
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
        const double rho = wet ? 2240 : 2500;
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


static void test_material_trials() {
    std::printf("[PT material pressure and repeated trials]\n");
    for (int model : {MatProps::rh_elastic, MatProps::rh_maxwell, MatProps::rh_viscous,
                      MatProps::rh_ep, MatProps::rh_evp}) {
        for (double alpha : {0.,0.6,1.}) {
            Fixture* m = new Fixture();
            build_mesh(*m);
            const int ne=m->var.nelem, nn=m->var.nnode;
            auto* strain = new tensor_t(ne,0.);
            auto* yy = new double_vec(ne,0.);
            auto* dp = new double_vec(ne,0.);
            auto* visc = new double_vec(ne,1e9);
            auto* pls = new double_vec(ne,0.);
            auto* dpls = new double_vec(ne,0.);
            auto* friction = new double_vec(ne,0.6);
            auto* state = new double_vec(ne,10.);
            auto* div = new double_vec(ne,0.);
            auto* heat = new double_vec(ne,7.);
            auto* pending_heat = new double_vec(ne,17.);
            m->var.shear_heat = heat;
            m->var.pending_heat = pending_heat;
            m->var.pore_plastic_increment = new double_vec(ne,19.);
            m->var.pore_dilation = new double_vec(ne,23.);
            m->param.control.has_shear_heating = true;
            m->strain_rate.resize(ne,0.);
            m->temperature.assign(nn,300.);
            m->pressure.assign(nn,10.);
            m->pressure_increment.assign(nn,-10.);
            m->markers.assign(ne,int_vec(1,1));
            m->var.strain=strain; m->var.stressyy=yy; m->var.dpressure=dp;
            m->var.viscosity=visc; m->var.plstrain=pls; m->var.delta_plstrain=dpls;
            m->var.dyn_fric_coeff=friction; m->var.state_variable=state;
            m->var.edvoldt=div; m->var.strain_rate=&m->strain_rate;
            m->var.temperature=&m->temperature; m->var.ppressure=&m->pressure;
            m->var.dppressure=&m->pressure_increment; m->var.elemmarkers=&m->markers;
            m->var.log_table=m->var.tan_table=m->var.sin_table=&m->lookup;
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
            m->var.mat=new MatProps(m->param,m->var);
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
                    double amc,anphi,anpsi,hardn,tenmax,mu=0.6,theta=10;
                    m->var.mat->plastic_props_rsf(0,0,amc,anphi,anpsi,hardn,tenmax,
                                                 slip,mu,theta,0.1,0,true);
                    if(slip==1e-9) first_strength=anphi;
                    check("initial RSF frozen friction",mu,0.6);
                    check("initial RSF frozen state",theta,10);
                    check("initial RSF independent of pseudo velocity",anphi,first_strength);
                }
            }

            if (model==MatProps::rh_elastic && alpha==0.6) {
                p.alpha={3e-5};
                delete m->var.mat; m->var.mat=new MatProps(m->param,m->var);
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
                delete m->var.mat; m->var.mat=new MatProps(m->param,m->var);
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
                    if (!running) check("failed solve threshold",result.threshold,0);
                    else check("running floor is finite and positive",std::isfinite(result.threshold) && result.threshold>0,1);
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
            if (model == MatProps::rh_ep && alpha == 0) {
                // Independent symmetric-edge states from the two active MC
                // equations and isotropic compliance (c=100 kPa, phi=30 deg).
                p.bulk_modulus={65e9/3}; p.shear_modulus={1e10};
                p.cohesion0=p.cohesion1={1e5};
                p.friction_angle0=p.friction_angle1={30};
                p.tension_max=1e9;
                m->param.control.has_hydraulic_diffusion=false;
                m->strain_rate.resize(ne,0.);
                for(int e=0;e<ne;++e) for(int d=0;d<NSTR;++d) m->strain_rate[e][d]=0;
                div->assign(ne,0.);
                for(int dilation : {0,30}) {
                    p.dilation_angle0=p.dilation_angle1={double(dilation)};
                    delete m->var.mat; m->var.mat=new MatProps(m->param,m->var);
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
                                m->stress[e][d]=base+(unique-base)*normal[d]*normal[d];
                            m->stress[e][NDIMS]=(unique-base)*normal[0]*normal[1];
#ifdef THREED
                            m->stress[e][4]=(unique-base)*normal[0]*normal[2];
                            m->stress[e][5]=(unique-base)*normal[1]*normal[2];
#endif
                            (*yy)[e]=base+(unique-base)*normal[2]*normal[2];
                            (*pls)[e]=(*dpls)[e]=0;
                        }
                        update_stress(m->param,m->var,m->stress,*yy,*dp,*visc,*strain,*pls,*dpls,
                                      m->strain_rate,m->pressure,m->pressure_increment,m->vel,
                                      *friction,*state,true);
                        #pragma acc wait
                        const double volume = (2*(base-target_base)+unique-target_unique)/(65e9);
                        const double plastic_base=(base-target_base-15e9*volume)/2e10;
                        const double plastic_unique=(unique-target_unique-15e9*volume)/2e10;
                        double increment;
#ifdef THREED
                        const double mean=(2*plastic_base+plastic_unique)/3;
                        increment=std::sqrt((2*(plastic_base-mean)*(plastic_base-mean)+
                            (plastic_unique-mean)*(plastic_unique-mean)+mean*mean)/2);
#else
                        increment=0.5*std::fabs(plastic_base-plastic_unique);
#endif
                        for(int e=0;e<ne;++e) {
                            for(int d=0;d<NDIMS;++d)
                                check("MC edge normal stress",m->stress[e][d],
                                    target_base+(target_unique-target_base)*normal[d]*normal[d]);
                            check("MC edge shear",m->stress[e][NDIMS],(target_unique-target_base)*normal[0]*normal[1]);
#ifdef THREED
                            check("MC edge rotated xz",m->stress[e][4],(target_unique-target_base)*normal[0]*normal[2]);
                            check("MC edge rotated yz",m->stress[e][5],(target_unique-target_base)*normal[1]*normal[2]);
#else
                            check("MC edge plane-strain yy",(*yy)[e],target_base+(target_unique-target_base)*normal[2]*normal[2]);
#endif
                            check("MC edge plastic history",(*pls)[e],increment);
                            check("MC edge plastic increment",(*dpls)[e],increment);
                        }
                    }
                }
            }
            if (model == MatProps::rh_ep && alpha == 0) {
                // The elastic trial exceeds c/tan(phi), but its admissible
                // associated shear return lies below it; no tensile flow.
                const double trial[3]={169e3,185e3,200e3};
                const double target[3]={165347.65655963385,170585.93935780309,170585.93935780309};
                const double volume=(trial[0]+trial[1]+trial[2]-target[0]-target[1]-target[2])/65e9;
                double plastic[3];
                for(int k=0;k<3;++k) plastic[k]=(trial[k]-target[k]-15e9*volume)/2e10;
                double increment;
#ifdef THREED
                const double mean=(plastic[0]+plastic[1]+plastic[2])/3;
                double squared=mean*mean;
                for(int k=0;k<3;++k) squared+=(plastic[k]-mean)*(plastic[k]-mean);
                increment=std::sqrt(squared/2);
#else
                increment=0.5*(std::max(plastic[0],std::max(plastic[1],plastic[2]))-
                               std::min(plastic[0],std::min(plastic[1],plastic[2])));
#endif
                for(int e=0;e<ne;++e) {
                    for(int d=0;d<NDIMS;++d) m->stress[e][d]=trial[d];
                    for(int d=NDIMS;d<NSTR;++d) m->stress[e][d]=0;
                    (*yy)[e]=trial[2];
                    (*pls)[e]=(*dpls)[e]=0;
                }
                update_stress(m->param,m->var,m->stress,*yy,*dp,*visc,*strain,*pls,*dpls,
                              m->strain_rate,m->pressure,m->pressure_increment,m->vel,*friction,*state,true);
                #pragma acc wait
                for(int e=0;e<ne;++e) {
                    for(int d=0;d<NDIMS;++d) check("shear return before tensile clipping",m->stress[e][d],target[d]);
                    for(int d=NDIMS;d<NSTR;++d) check("uncut shear has no rotation",m->stress[e][d],0);
                    check("uncut shear plastic history",(*pls)[e],increment);
                    check("uncut shear plastic increment",(*dpls)[e],increment);
#ifndef THREED
                    check("uncut shear yy",(*yy)[e],target[2]);
#endif
                }
            }
            delete baseline; delete m->var.mat;
            delete strain; delete yy; delete dp; delete visc; delete pls; delete dpls;
            delete friction; delete state; delete div;
            delete heat; delete pending_heat; delete m;
        }
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
    m->temperature.assign(nn,273.0);
    m->pressure.assign(nn,0.0);
    m->pressure_increment.assign(nn,0.0);
    m->markers.assign(ne,int_vec(1,1));
    m->var.temperature=&m->temperature;
    m->var.ppressure=&m->pressure;
    m->var.dppressure=&m->pressure_increment;
    m->var.elemmarkers=&m->markers;
    m->var.log_table=m->var.tan_table=m->var.sin_table=&m->lookup;
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
    m->var.mat=new MatProps(m->param,m->var);
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
        delete m->var.mat;
        m->var.mat=new MatProps(m->param,m->var);
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
        delete m->var.mat; m->var.mat=new MatProps(m->param,m->var);
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
        delete m->var.mat; m->var.mat=new MatProps(m->param,m->var);
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
    delete m->var.mat; m->var.mat=new MatProps(m->param,m->var);
    {
        p.hydraulic_perm={0}; p.matrix_storage={1}; p.matrix_exchange={0.1};
        delete m->var.mat; m->var.mat=new MatProps(m->param,m->var);
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
        delete m->var.mat; m->var.mat=new MatProps(m->param,m->var);
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
    delete m->var.mat;
    m->var.mat=new MatProps(m->param,m->var);
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
        delete m->var.mat;
        m->var.mat=new MatProps(m->param,m->var);
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
    delete m->var.mat;
    m->var.mat=new MatProps(m->param,m->var);
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
    m->temperature.assign(nn, 273);
    m->pressure.assign(nn, 0);
    m->markers.assign(ne, int_vec(1, 1));
    m->var.temperature = &m->temperature;
    m->var.ppressure = &m->pressure;
    m->var.elemmarkers = &m->markers;
    m->var.log_table = m->var.tan_table = m->var.sin_table = &m->lookup;
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
    m->var.mat = new MatProps(m->param, m->var);
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
    const int nn=m->var.nnode, ne=m->var.nelem;
    m->temperature.assign(nn,273.); m->markers.assign(ne,int_vec(1,1));
    m->var.temperature=&m->temperature; m->var.elemmarkers=&m->markers;
    m->var.log_table=m->var.tan_table=m->var.sin_table=&m->lookup;
    auto& p=m->param.mat;
    p.nmat=1; p.rheol_type=MatProps::rh_elastic;
    p.rho0={2650}; p.alpha={3e-5}; p.porosity={.2}; p.biot_coeff={.6};
    p.bulk_modulus={2e7}; p.shear_modulus={1e7};
    p.visc_exponent={1}; p.visc_coefficient={1};
    p.heat_capacity={1000}; p.therm_cond={1}; p.fluid_bulk_modulus={1e8};
    p.bulk_modulus_s={1e8}; p.fluid_rho0={1000}; p.fluid_alpha={2e-4};
    p.hydraulic_perm={0}; p.fluid_visc={1e-3};
    m->var.mat=new MatProps(m->param,m->var);
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
    m->temperature.assign(nn,273.01); m->pressure.assign(nn,5e4);
    m->pressure_increment.assign(nn,999); m->markers.assign(ne,int_vec(1,1));
    m->var.temperature=&m->temperature; m->var.ppressure=&m->pressure;
    m->var.dppressure=&m->pressure_increment; m->var.elemmarkers=&m->markers;
    m->var.log_table=m->var.tan_table=m->var.sin_table=&m->lookup;
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
    m->var.mat=new MatProps(m->param,m->var);
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
    m->temperature.assign(nn,300); m->pressure.assign(nn,0);
    m->markers.assign(ne,int_vec(1,1));
    m->var.temperature=&m->temperature; m->var.ppressure=&m->pressure;
    m->var.elemmarkers=&m->markers;
    m->var.log_table=m->var.tan_table=m->var.sin_table=&m->lookup;
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
    m->var.mat=new MatProps(m->param,m->var);
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
    delete m->var.mat;
    m->var.mat=new MatProps(m->param,m->var);
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
    delete m->var.mat;
    m->var.mat=new MatProps(m->param,m->var);
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
    assemble_point_sources(source,m->var,0,1,rate);
    for(int n=0;n<m->var.nnode;++n) {
        const double expected=n==0 ? 2*c+2.8 : n<=NDIMS ? 3*c : -c;
        check("point barycentric rates and shared-node ownership",rate[n],expected);
    }
    const double_vec inflow={1,0,1}, temperature={300,0,400};
    assemble_point_sources(source,m->var,0,1,rate,&inflow,true);
    for(int n=0;n<m->var.nnode;++n)
        check("active inlet rate bounds partial source windows",rate[n],
              n==0 ? 5*c+7 : n<=NDIMS ? 5*c : 0);
    assemble_point_sources(source,m->var,0,1,rate,&temperature,true);
    for(int n=0;n<m->var.nnode;++n)
        check("temperature-weighted active inlet rate",rate[n],
              n==0 ? 1500*c+2800 : n<=NDIMS ? 1500*c : 0);
    source.temperature={300,273,400};
    double_vec heat(m->var.nnode,0),separate_rate(m->var.nnode,0),separate_heat(m->var.nnode,0);
    for (bool peak : {false,true}) {
        assemble_point_sources(source,m->var,0,1,separate_rate,&inflow,peak);
        assemble_point_sources(source,m->var,0,1,separate_heat,&temperature,peak);
        assemble_point_sources(source,m->var,0,1,rate,&inflow,peak,&heat,&source.temperature);
        for(int n=0;n<m->var.nnode;++n) {
            check("shared point map preserves inlet volume",rate[n],separate_rate[n]);
            check("shared point map preserves weighted inlet",heat[n],separate_heat[n]);
        }
    }
    assemble_point_sources(source,m->var,.81,1.4,rate,&inflow,true);
    for(double v:rate) check("inactive inlets do not constrain heat steps",v,0);
    const double breaks[]={0,.13,.47,.81,1.4};
    for(int i=0;i<4;++i) {
        assemble_point_sources(source,m->var,breaks[i],breaks[i+1],rate);
        for(int n=0;n<m->var.nnode;++n) integrated[n]+=rate[n]*(breaks[i+1]-breaks[i]);
    }
    double total=0;for(double v:integrated) total+=v;
    check("subintervals conserve signed scheduled amount",total,4.4);
    source.rate_model=injection_rate_total_amount;
    assemble_point_sources(source,m->var,0,1.4,rate);
    for(int n=0;n<m->var.nnode;++n)
        check("rate/total-amount equivalence",rate[n]*1.4,integrated[n]);
    assemble_point_sources(source,m->var,0,1.4,rate,&inflow,true);
    for(int n=0;n<m->var.nnode;++n)
        check("total-amount inlet has the same active rate bound",rate[n],
              n==0 ? 5*c+7 : n<=NDIMS ? 5*c : 0);
    for(int n=0;n<m->var.nnode;++n) for(int d=0;d<NDIMS;++d) m->coord[n][d]+=10;
    for(int i=0;i<source.num_points;++i) {
        source.points_x[i]+=10;source.points_y[i]+=10;source.points_z[i]+=10;
    }
    assemble_point_sources(source,m->var,0,1.4,rate);
    for(int n=0;n<m->var.nnode;++n)
        check("mapping uses current translated mesh",rate[n]*1.4,integrated[n]);
    source.enabled=false;
    assemble_point_sources(source,m->var,0,1.4,rate);
    for(double v:rate) check("disabled point source clears derived rates",v,0);
    delete m;
}

int main() {
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
    test_constraints();
    test_cancellation();
    test_traction();
    test_platen();
    test_damping();
    test_body_force();
    std::printf("\n%d passed, %d failed\n", npass, nfail);
    return nfail == 0 ? 0 : 1;
}
