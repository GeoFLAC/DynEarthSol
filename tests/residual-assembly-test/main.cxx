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
        const double rho = wet ? 2200 : 2500;
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


// One-material element fields for constitutive updates on the test mesh.
struct MaterialFixture {
    Fixture* m;
    tensor_t* strain;
    double_vec *yy, *dp, *visc, *pls, *dpls, *friction, *state, *div;
    MaterialFixture(int model, double alpha) {
        m = new Fixture();
        build_mesh(*m);
        const int ne=m->var.nelem, nn=m->var.nnode;
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
    void set_material() {
        delete m->var.mat;
        m->var.mat=new MatProps(m->param,m->var);
    }
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
                baseline->restore(m->var);
                for(int e=0;e<ne;++e) {
                    check("restore stressyy exact",(*m->var.stressyy)[e]==baseline->stressyy[e],1);
                    check("restore dpressure exact",(*m->var.dpressure)[e]==baseline->dpressure[e],1);
                    check("restore viscosity exact",(*m->var.viscosity)[e]==baseline->viscosity[e],1);
                    check("restore plstrain exact",(*m->var.plstrain)[e]==baseline->plstrain[e],1);
                    check("restore delta_plstrain exact",(*m->var.delta_plstrain)[e]==baseline->delta_plstrain[e],1);
                    check("restore dyn_fric_coeff exact",(*m->var.dyn_fric_coeff)[e]==baseline->dyn_fric_coeff[e],1);
                    check("restore state_variable exact",(*m->var.state_variable)[e]==baseline->state_variable[e],1);
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
            delete baseline; delete m->var.mat;
            delete strain; delete yy; delete dp; delete visc; delete pls; delete dpls;
            delete friction; delete state; delete div; delete m;
        }
    }
}

// The plastic-strain measure stored in history (3D invariant, plane-strain
// half range), computed independently of rheology.cxx.
static double recorded_plastic_increment(const double plastic[3]) {
#ifdef THREED
    const double mean=(plastic[0]+plastic[1]+plastic[2])/3;
    double squared=mean*mean;
    for(int k=0;k<3;++k) squared+=(plastic[k]-mean)*(plastic[k]-mean);
    return std::sqrt(squared/2);
#else
    return 0.5*(std::max(plastic[0],std::max(plastic[1],plastic[2]))-
                std::min(plastic[0],std::min(plastic[1],plastic[2])));
#endif
}

// Pure 2D (two principal stresses, not plane strain): the same invariant of
// the two in-plane components.
static double recorded_plastic_increment_2d(const double plastic[3]) {
    const double mean=(plastic[0]+plastic[1])/2;
    double squared=mean*mean;
    for(int k=0;k<2;++k) squared+=(plastic[k]-mean)*(plastic[k]-mean);
    return std::sqrt(squared/2);
}

// Plastic strain of a principal stress correction under the compliance of
// set_mc_material: 3K = 65e9, lambda = K - 2G/3 = 15e9, 2G = 2e10.
static void compliance_plastic(const double trial[3], const double result[3], double plastic[3]) {
    const double volume=(trial[0]+trial[1]+trial[2]-result[0]-result[1]-result[2])/65e9;
    for(int k=0;k<3;++k) plastic[k]=(trial[k]-result[k]-15e9*volume)/2e10;
}

// Pure 2D: the in-plane moduli a1 = K + 4G/3 = 35e9 and a2 = K - 2G/3 = 15e9
// of set_mc_material, inverted (determinant a1^2 - a2^2 = 1e21).
static void compliance_plastic_2d(const double trial[3], const double result[3], double plastic[3]) {
    const double d0=trial[0]-result[0], d1=trial[1]-result[1];
    plastic[0]=(35e9*d0-15e9*d1)/1e21;
    plastic[1]=(35e9*d1-15e9*d0)/1e21;
    plastic[2]=0;
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
                             double minimum_history=0, double initial_history=0);

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
            compliance_plastic(edge_trial,edge_target,plastic);
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
        compliance_plastic(trial,target,plastic);
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
        double amc,anphi,anpsi,hardn,tenmax;
        m->var.mat->plastic_props(0,0,amc,anphi,anpsi,hardn,tenmax);
        check("MC nonzero cohesion slope",hardn,slope);
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
            double new_amc,new_phi,new_psi,new_slope,new_cap;
            m->var.mat->plastic_props(0,pls[0],new_amc,new_phi,new_psi,new_slope,new_cap);
            const double high=NDIMS==2 ? yy[0] : stress[0][2];
            check("MC stored history reproduces yield",stress[0][0]-new_phi*high+new_amc,0);
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

// Properties of any accepted return from a principal trial with no plastic
// history, checked against the law at the stored history: inside the shear
// surface and tensile cap and on one of them, principal order kept, history
// equal to the measure of the plastic strain implied by isotropic compliance,
// and no change from a second, unloaded update.
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

static void check_law_return(MaterialFixture& f, const char* name, const double trial[3],
                             double minimum_history, double initial_history) {
    const auto flag=[name](const char* property, bool ok) {
        char what[256];
        std::snprintf(what,sizeof what,"%s: %s",name,property);
        check(what,ok,1);
    };
    auto& stress=f.m->stress;
    // Pure 2D has two principal stresses, from the in-plane tensor.
    const int nd=NDIMS==2 && !f.m->param.mat.is_plane_strain ? 2 : 3;
    const auto principal=[&](double r[3]) {
        if(nd==2) {
            const double mean=0.5*(stress[0][0]+stress[0][1]);
            const double radius=std::hypot(0.5*(stress[0][0]-stress[0][1]),stress[0][2]);
            r[0]=mean-radius; r[1]=mean+radius; r[2]=0;
            return;
        }
        r[0]=stress[0][0]; r[1]=stress[0][1]; r[2]=NDIMS==2 ? (*f.yy)[0] : stress[0][2];
    };
    f.set_principal(trial);
    for(int e=0;e<f.m->var.nelem;++e) (*f.pls)[e]=initial_history;
    f.update();
    double result[3];
    principal(result);
    double plastic[3];
    if(nd==2) compliance_plastic_2d(trial,result,plastic);
    else compliance_plastic(trial,result,plastic);
    const double recorded=nd==2 ? recorded_plastic_increment_2d(plastic)
                                : recorded_plastic_increment(plastic);
    const double history=(*f.pls)[0];
    double amc,anphi,anpsi,hardn,cap;
    f.m->var.mat->plastic_props(0,history,amc,anphi,anpsi,hardn,cap);
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
    f.update();
    double again[3];
    principal(again);
    for(int k=0;k<nd;++k) flag("unloaded stress unchanged",std::fabs(again[k]-result[k])<=1e-9*scale);
    flag("unloaded increment is zero",(*f.dpls)[0]<=1e-12*increment);
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
        auto& stress=f.m->stress;
        const double history=(*f.pls)[0];
        double amc,anphi,anpsi,hardn,cap;
        f.m->var.mat->plastic_props(0,history,amc,anphi,anpsi,hardn,cap);
        const double low=0.5*(stress[0][0]+stress[0][1])-std::hypot(0.5*(stress[0][0]-stress[0][1]),stress[0][2]);
        const double high=stress[0][0]+stress[0][1]-low;
        check("pure 2D softening stays in the linear range",history>0 && history<1e-4,1);
        check("pure 2D history is the increment",history,(*f.dpls)[0]);
        check("pure 2D stored history reproduces yield",(low-anphi*high+amc)/amc,0);
    }
}

// Pure 2D returns the closed form cannot give, checked by the properties of
// check_law_return: softening past the segment (also from a history one ulp
// short of the breakpoint), friction and dilation softening within the step,
// and tension at the cap, its corner with the shear face and the cap-cap corner.
static void test_pure_2d_law_returns() {
    std::printf("[Pure 2D Mohr-Coulomb returns over the strength law]\n");
    MaterialFixture f(MatProps::rh_ep, 0);
    auto& p=f.m->param.mat;
    p.is_plane_strain=false;
    const double face[3]={-1e6,-1e5,0};
    for(int dilation : {0,30}) {
        p.pls0={0}; p.pls1={1e-8};
        set_mc_material(f, dilation);
        p.cohesion1={9e4};
        f.set_material();
        check_law_return(f,"pure 2D softening past the segment",face,1e-8);
    }
    p.pls0={0}; p.pls1={0.0009287493630057339};
    set_mc_material(f, 0);
    p.cohesion1={9e4};
    f.set_material();
    const double deep[3]={-1e8,-1e7,0};
    check_law_return(f,"pure 2D softening past the segment from nonzero history",deep,
                     p.pls1[0],0.0003829357836281204);
    p.pls0={0}; p.pls1={1e-5};
    set_mc_material(f, 10);
    p.friction_angle1={20}; p.dilation_angle1={0};
    f.set_material();
    check_law_return(f,"pure 2D friction and dilation soften within the step",face);
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
    // The default cap c/tan(phi) meets the face at p0 = p1 = c*sqrt(3).
    set_mc_material(f, 0);
    const double hydrostatic[3]={3e5,3.2e5,0};
    check_law_return(f,"pure 2D face-cap apex",hydrostatic);
    for(int d=0;d<2;++d) check("pure 2D apex stress",f.m->stress[0][d],1e5*std::sqrt(3.));
}
#endif

// Returns the closed form cannot give: a cohesion drop faster than elastic
// unloading, friction and dilation softening within the step, and (3D) the
// tensile edge, apex and tension-shear corners.
static void test_mohr_coulomb_law_returns() {
    std::printf("[Mohr-Coulomb returns over the strength law]\n");
    MaterialFixture f(MatProps::rh_ep, 0);
    auto& p=f.m->param.mat;
    const double edge[3]={-1e6,-1e6,-1e5}, face[3]={-1e6,-5e5,-1e5};
    // Cohesion drops 1e5 -> 9e4 over 1e-8: the return reaches the residual strength.
    for(int dilation : {0,30}) {
        p.pls0={0}; p.pls1={1e-8};
        set_mc_material(f, dilation);
        p.cohesion1={9e4};
        f.set_material();
        check_law_return(f,"softening past the segment, edge",edge,1e-8);
        check_law_return(f,"softening past the segment, face",face,1e-8);
    }
    // Nonzero initial history whose pls + (pls1 - pls) rounds one ulp below
    // pls1: the search must still cross the breakpoint to the residual root.
    p.pls0={0}; p.pls1={0.0009287493630057339};
    set_mc_material(f, 0);
    p.cohesion1={9e4};
    f.set_material();
    const double deep[3]={-1e8,-1e8,-1e7};
    check_law_return(f,"softening past the segment from nonzero history",deep,
                     p.pls1[0],0.0003829357836281204);
    // Friction 30 -> 20 deg and dilation 10 -> 0 deg over 1e-5: the step
    // ends inside the segment with strength and flow of the stored history.
    p.pls0={0}; p.pls1={1e-5};
    set_mc_material(f, 10);
    p.friction_angle1={20}; p.dilation_angle1={0};
    f.set_material();
    check_law_return(f,"friction and dilation soften within the step",face);
    check_law_return(f,"friction and dilation soften within the step, edge",edge);
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
            check_law_return(f,"tension with a constant or nearly constant friction angle",trial);
            returned[perturbed][0]=f.m->stress[0][0];
            returned[perturbed][1]=f.m->stress[0][1];
            returned[perturbed][2]=NDIMS==2 ? (*f.yy)[0] : f.m->stress[0][2];
        }
        for(int k=0;k<3;++k)
            check("tension return continuous in the friction angle",
                  std::fabs(returned[0][k]-returned[1][k])<=1e-6*cohesion,1);
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
    // The default cap c/tan(phi) passes through the cone apex.
    set_mc_material(f, 0);
    const double hydrostatic[3]={3e5,3.2e5,3.5e5};
    check_law_return(f,"cone apex",hydrostatic);
    for(int d=0;d<3;++d) check("cone apex stress",f.m->stress[0][d],1e5*std::sqrt(3.));
#endif
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
    // stress and the stored friction coefficient follow the final history.
    // With a = b = 0 the coefficient is tan(phi(pls)); with the aging law the
    // stress must lie on the face that the stored coefficient defines.
    const double face[3]={-1e6,-5e5,-1e5};
    for(double end : {1e-5,1e-8}) for(int aging=0;aging<2;++aging) {
        p.pls0={0}; p.pls1={end};
        set_mc_material(f, 0);
        p.friction_angle1={20};
        p.direct_a={aging ? 0.01 : 0.}; p.evolution_b={aging ? 0.015 : 0.};
        p.state_var_model=aging;
        f.set_material();
        if(!aging) check_law_return(f,"rate-and-state friction softening",face);
        // One update: the next one would recompute the coefficient anyway.
        f.set_principal(face);
        f.friction->assign(f.m->var.nelem,0.6);
        f.state->assign(f.m->var.nelem,10.);
        f.update();
        const double history=(*f.pls)[0], mu=(*f.friction)[0];
        double amc,anphi,anpsi,hardn,cap;
        f.m->var.mat->plastic_props(0,history,amc,anphi,anpsi,hardn,cap);
        const double cohesion=amc/(2*std::sqrt(anphi));
        const double static_mu=std::tan(std::asin((anphi-1)/(anphi+1)));
        check("rate-and-state friction softening yields",history>0,1);
        if(!aging) check("stored friction coefficient of the final history",mu,static_mu);
        else check("stored friction coefficient differs from the static one",
                   std::fabs(mu-static_mu)>1e-3,1);
        const double sine=std::sin(std::atan(mu)), nphi=(1+sine)/(1-sine);
        auto& stress=f.m->stress;
        const double high=NDIMS==2 ? (*f.yy)[0] : stress[0][2];
        check("stress on the face of the stored friction coefficient",
              (stress[0][0]-nphi*high+2*cohesion*std::sqrt(nphi))/(2*cohesion),0);
    }
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
    p.fluid_bulk_modulus={2e9}; p.hydraulic_perm={1e-12};
    p.fluid_visc={1e-3}; p.fluid_rho0={1000}; p.fluid_alpha={0};
    m->var.mat=new MatProps(m->param,m->var);
    double_vec reference(nn,0.0);
    double reference_diffusivity=0;
    for(int driven=0;driven<2;++driven) for(double gravity : {0.0,10.0}) {
        m->param.control.gravity=gravity;
        for(int n=0;n<nn;++n) {
            m->pressure[n]=2e5+driven*1e5*m->coord[n][0]-1000*gravity*m->coord[n][NDIMS-1];
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
    delete m->var.mat; delete flags; delete old_mean; delete rate; delete mass; delete m;
}

int main(int argc, char** argv) {
    if(argc==2 && (std::strcmp(argv[1],"--reject-mc-return")==0 ||
                   std::strcmp(argv[1],"--reject-mc-return-pure-2d")==0)) {
        test_rejected_mc_return(NDIMS==2 && std::strcmp(argv[1],"--reject-mc-return-pure-2d")==0);
        return 1;  // check.sh requires status 54 from update_stress
    }
    test_darcy_gravity_limit();
    test_pt_norm();
    test_material_trials();
    test_mohr_coulomb_returns();
    test_mohr_coulomb_law_returns();
    test_rsf_law_returns();
#ifndef THREED
    test_pure_2d_softening_return();
    test_pure_2d_law_returns();
#endif
    test_constraints();
    test_cancellation();
    test_traction();
    test_damping();
    test_body_force();
    std::printf("\n%d passed, %d failed\n", npass, nfail);
    return nfail == 0 ? 0 : 1;
}
