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
#include "bc.hpp"
#include "rheology.hpp"
#include "pseudo_transient.hpp"
#include <limits>
#include "matprops.hpp"

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
    tensor_t strain_rate;
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
    delete m->var.bfacets[iboundx0];
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
    p.therm_cond = {1.0}; p.fluid_bulk_modulus = {2e9};
    m->param.control.gravity = 9.8;
    m->param.control.damping_option = 2;
    m->param.control.damping_factor = 0.25;
    m->var.mat = new MatProps(m->param, m->var);
    run_update_force(*m);
    // Each simplex distributes its weight equally among its vertices. The
    // shared vertex receives two contributions, all other vertices one.
    const double weight = (2500.0*0.8 + 1000.0*0.2)*9.8*m->volume[0];
    for (int n=0; n<m->var.nnode; ++n)
        for (int d=0; d<NDIMS; ++d) {
            const double raw = d == NDIMS-1 ? -weight*(n == 0 ? 2 : 1)/NODES_PER_ELEM : 0;
            check("damped gravity", m->force[n][d], raw*0.75);
#ifdef B0
            check("undamped gravity", m->force_undamped[n][d], raw);
#endif
        }
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
            p.biot_coeff={alpha}; p.heat_capacity={1}; p.therm_cond={1}; p.fluid_bulk_modulus={2e9};
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
                    check("plastic increment",(*dpls)[e],0);
                    check("state unchanged",(*state)[e],10);
                    // Deliberately contaminate every snapshot scalar to check restore.
                    (*yy)[e]=(*dp)[e]=(*visc)[e]=(*pls)[e]=(*dpls)[e]=(*friction)[e]=(*state)[e]=99;
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
                auto* after=new MechanicalState(ne);
                for(int cap : {1,3}) {
                    baseline->restore(m->var);
                    for(int n=0;n<nn;++n) for(int d=0;d<NDIMS;++d)
                        m->vel[n][d]=d==0 ? 2*m->coord[n][NDIMS-1] : 0;
                    m->param.control.PT_max_iter=cap;
                    const PTResult result=run_physical_step_pt(m->param,m->var);
                    #pragma acc wait
                    check("failed solve status",result.status==PTStatus::max_iterations,1);
                    check("failed solve iterations",result.iterations,cap);
                    check("failed solve has imbalance",result.residual>0,1);
                    check("failed solve threshold",result.threshold,0);
                    after->capture(m->var);
                    for(int e=0;e<ne;++e) {
                        for(int d=0;d<NSTR;++d) {
                            check("rollback stress exact",after->stress[e][d]==baseline->stress[e][d],1);
                            check("rollback strain exact",after->strain[e][d]==baseline->strain[e][d],1);
                        }
                        for(int k=0;k<7;++k)
                            check("rollback history exact",after->scalars[e][k]==baseline->scalars[e][k],1);
                    }
                    for(int n=0;n<nn;++n) {
                        for(int d=0;d<NDIMS;++d)
                            check("rollback velocity exact",m->vel[n][d]==(d==0 ? 2*m->coord[n][NDIMS-1] : 0),1);
                        check("pending pressure retained",m->pressure_increment[n]==-10,1);
                        check("pressure retained",m->pressure[n]==10,1);
                    }
                    check("physical time retained",m->var.time==1,1);
                    check("physical dt retained",m->var.dt==0.1,1);
                    check("stress allocation retained",m->var.stress==&m->stress,1);
                    check("strain allocation retained",m->var.strain==strain,1);
                    check("viscosity allocation retained",m->var.viscosity==visc,1);
                }
                delete after; delete ntmp; delete etmp; delete volume_n;
                delete flags; delete normals;
            }
            delete baseline; delete m->var.mat;
            delete strain; delete yy; delete dp; delete visc; delete pls; delete dpls;
            delete friction; delete state; delete div; delete m;
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

int main() {
    test_pt_norm();
    test_material_trials();
    test_constraints();
    test_cancellation();
    test_traction();
    test_damping();
    test_body_force();
    std::printf("\n%d passed, %d failed\n", npass, nfail);
    return nfail == 0 ? 0 : 1;
}
