#include <iostream>
#include <limits>

#include "constants.hpp"
#include "parameters.hpp"
#include "bc.hpp"
#include "matprops.hpp"
#include "fields.hpp"
#include "mpi_simulation.hpp"
#include "hydraulics.hpp"
#include "utils.hpp"
#include "point_sources.hpp"

// Nodal point-source fields, sized to the current mesh.
static void allocate_source_fields(const Param &param, Variables &var)
{
    const int n = var.nnode;
    delete var.fluid_point_source;
    var.fluid_point_source = param.injection.enabled ? new double_vec(n,0) : nullptr;
    delete var.heat_source;
    var.heat_source = param.heat_injection.enabled ? new double_vec(n,0) : nullptr;
    delete var.inflow_rate;
    delete var.inflow_heat;
    var.inflow_rate = param.injection.temperature.empty() ? nullptr : new double_vec(n,0);
    var.inflow_heat = param.injection.temperature.empty() ? nullptr : new double_vec(n,0);
}

void allocate_variables(const Param &param, Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    var.volume = new double_vec(e);
    var.volume_n = new double_vec(n);

    var.mass = new double_vec(n);
    var.tmass = new double_vec(n);
    var.ymass = new double_vec(n);
    var.edvoldt = new double_vec(e);

    allocate_hydraulic_state(var);
    if (param.control.hydraulic_content_model == 1) {
        var.fluid_reference_temperature = param.bc.surface_temperature;
        var.fluid_split = new double_vec(n,0);
        var.solid_mass = new double_vec(n,0);
        var.solid_capacity = new double_vec(n,0);
        var.pending_temperature = new double_vec(n,0);
        if (has_temperature_transport(param)) var.energy_balance = new double_vec(6,0);
        var.fluid_balance->assign(6,0);
    }
    allocate_source_fields(param, var);
    if (std::any_of(param.mat.matrix_exchange.begin(), param.mat.matrix_exchange.end(),
                    [](double h) { return h > 0; }))
        var.matrix_pressure = new double_vec(e, 0);
    if (param.control.has_thermal_pressurization)
        var.old_temperature = new double_vec(var.nnode);
    if (param.control.has_shear_heating) {
        var.shear_heat = new double_vec(e, 0);
        var.pending_heat = new double_vec(e, 0);
    }
    if (param.mat.porosity_model)
        var.pore_strain = new double_vec(e, 0);
    if (std::any_of(param.mat.pore_healing_time.begin(), param.mat.pore_healing_time.end(),
                    [](double t) { return t > 0; })) {
        var.pore_dilation = new double_vec(e,0);
        var.pore_plastic_increment = new double_vec(e,0);
        var.pore_recovery = new double_vec(e,0);
    }

    {
        // these fields are reallocated during remeshing interpolation
        var.stress = new tensor_t(e, 0);
        var.stressyy = new double_vec(e, 0);
        var.volume_old = new double_vec(e); // for dv remeshing interpolation
        var.temperature = new double_vec(n);
        var.coord0 = new array_t(n);
        var.plstrain = new double_vec(e);
        var.delta_plstrain = new double_vec(e);
        var.plstrain_remesh = new double_vec(e);
        var.vel = new array_t(n, 0);
        var.strain = new tensor_t(e, 0);
        var.radiogenic_source = new double_vec(e, 0);
        var.dyn_fric_coeff = new double_vec(e);
        var.state_variable = new double_vec(e);

    }

    var.ntmp = new double_vec(n);
    var.dpressure = new double_vec(e, 0);

    var.viscosity = new double_vec(e,param.mat.visc_max);

    var.force = new array_t(n, 0);
    var.force_residual = new array_t(n, 0);
    if (param.control.has_PT) var.PT_force_residual_valid = false;

    var.strain_rate = new tensor_t(e, 0);

    var.mat = new MatProps(param, var);

    var.tmp_result = new elem_cache(e);
    var.etmp = new double_vec(e);
    var.etmp_int = new int_vec(e);
}


void reallocate_tmp(const Param& param, Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    delete var.ntmp;
    var.ntmp = new double_vec(n);

    delete var.tmp_result;
    var.tmp_result = new elem_cache(e);
    delete var.etmp;
    var.etmp = new double_vec(e);
    delete var.etmp_int;
    var.etmp_int = new int_vec(e);
}


void reallocate_variables(const Param& param, Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    if (var.fluid_split) {
        var.fluid_split->resize(n);
        var.solid_mass->resize(n);
        var.solid_capacity->resize(n);
        var.pending_temperature->resize(n);
    }

    delete var.volume;
    delete var.volume_n;
    var.volume = new double_vec(e);
    var.volume_n = new double_vec(n);

    delete var.mass;
    delete var.tmass;
    var.mass = new double_vec(n);
    var.tmass = new double_vec(n);

    reset_hydraulic_assembly(param, var);
    allocate_source_fields(param, var);
    delete var.ymass;
    var.ymass = new double_vec(n);
  
    delete var.edvoldt;
    var.edvoldt = new double_vec(e);
    if (var.pore_recovery) var.pore_recovery->assign(e, 0.0);


    delete var.dpressure;
    var.dpressure = new double_vec(e, 0);
    delete var.viscosity;
    var.viscosity = new double_vec(e,param.mat.visc_max);
    delete var.force;
    var.force = new array_t(n, 0);

    delete var.force_residual;
    var.force_residual = new array_t(n, 0);
    if (param.control.has_PT) var.PT_force_residual_valid = false;

    delete var.strain_rate;
    var.strain_rate = new tensor_t(e, 0);

    // var.stress and var.stressyy are NOT reallocated here: the NN interpolation
    // already swapped in remapped copies sized to the new mesh, which
    // spr_node_to_elem reads before overwriting with the SPR average.

    // Must run AFTER remap_markers(), which deletes and recreates *var.elemmarkers:
    // MatProps binds it by REFERENCE and its constructor seeds the property-mean cache.
    // Called earlier the cache is built from empty counts -- wrong means, not a crash.
    // It also binds strain_rate, recreated above. remesh() freed it (null) when the old
    // mesh was moved out.
    delete var.mat;
    var.mat = new MatProps(param, var);

    apply_temperature_bcs(param, var, *var.temperature);

}


namespace {
// Galerkin transport with symmetric graph diffusion for advection; rows sum
// to zero. Zero flow retains ordinary Galerkin conduction on obtuse elements.
#pragma acc routine seq
void thermal_transport_operator(const Param &param, const Variables &var, int e,
                                double a[NODES_PER_ELEM][NODES_PER_ELEM])
{
    double g[NDIMS][NODES_PER_ELEM];
    get_shape_gradients(var, e, g);
    double q[NDIMS];
    hydraulic_element_flux(param, var, e, q, g);
    const double kv = param.control.has_thermal_diffusion ? var.mat->k(e)*(*var.volume)[e] : 0;
    const double cv = var.mat->rho_fluid_reference(e)*var.mat->cp_fluid(e)*
        (*var.volume)[e]/NODES_PER_ELEM;
    double c[NODES_PER_ELEM];  // lumped advection; the same for every test node
    for (int j = 0; j < NODES_PER_ELEM; ++j) {
        c[j] = 0;
        for (int d = 0; d < NDIMS; ++d) c[j] += cv*q[d]*g[d][j];
    }
    for (int i = 0; i < NODES_PER_ELEM; ++i)
        for (int j = 0; j < NODES_PER_ELEM; ++j) {
            a[i][j] = c[j];
            for (int d = 0; d < NDIMS; ++d) a[i][j] += kv*g[d][i]*g[d][j];
        }
    for (int i = 0; i < NODES_PER_ELEM; ++i)
        for (int j = i+1; j < NODES_PER_ELEM; ++j) {
            // Clip positive couplings, but by no more than advection adds.
            const double diffusion = std::max(0.0, std::min(std::max(a[i][j], a[j][i]),
                                                            std::max(c[i], c[j])));
            a[i][j] -= diffusion; a[j][i] -= diffusion;
            a[i][i] += diffusion; a[j][j] += diffusion;
        }
}
} // namespace

#pragma acc routine seq
bool prescribed_temperature(const Param &param, const Variables &var, int n, double &value)
{
    const uint flags[6] = {BOUNDX0, BOUNDX1, BOUNDY0, BOUNDY1, BOUNDZ0, BOUNDZ1};
    for (int i = 0; i < 6; ++i)
        if (param.bc.tbc[i] == 1 && ((*var.bcflag)[n] & flags[i])) {
            value = thermal_boundary_value(param,i);
            return true;
        }
    return false;
}

void apply_temperature_bcs(const Param &param, const Variables &var, double_vec &temperature)
{
    if (!param.bc.custom_temperature_boundary) return;
    const uint flags[6] = {BOUNDX0, BOUNDX1, BOUNDY0, BOUNDY1, BOUNDZ0, BOUNDZ1};
    for (int n = 0; n < var.nnode; ++n) {
        bool fixed = false;
        double value = 0;
        for (int i = 0; i < 6; ++i) {
            if (param.bc.tbc[i] != 1 || !((*var.bcflag)[n] & flags[i])) continue;
            const double next = thermal_boundary_value(param, i);
            if (!std::isfinite(next) || (fixed && std::abs(next-value) > 1e-12*std::max(1.0, std::abs(value))))
                die(EXIT_CONFIG_VALUE, "Intersecting temperature boundaries require equal finite values.");
            fixed = true; value = next;
        }
        if (fixed) temperature[n] = value;
    }
}

// Volumetric heat capacity of the injected fluid, from the reference material.
static double inflow_capacity(const Param &param)
{
    const int m = param.mat.mattype_ref;
    return param.mat.fluid_rho0[m] * param.mat.fluid_heat_capacity[m];
}

// Fluid injected at a prescribed temperature over [start, stop]: the nodal
// injected volume rate Q_n and its temperature-weighted sum. Extraction leaves
// at the local temperature, which the advective transport already carries, so
// extracting points are left out.
static void assemble_inflow(const Param &param, const Variables &var, double start, double stop,
                            bool peak_rate = false)
{
    const Injection &source = param.injection;
    double_vec injecting(source.num_points);
    for (int p = 0; p < source.num_points; ++p) {
        const double strength = source.rate_model == injection_rate_total_amount ?
            source.total_amount[p] : source.rate[p];
        injecting[p] = strength > 0 ? 1.0 : 0.0;
    }
    assemble_point_sources(source, var.fluid_points, var, start, stop, *var.inflow_rate, &injecting, peak_rate,
                           var.inflow_heat, &source.temperature);
}

static double temperature_step_limit(const Param &param, const Variables &var, elem_cache &work)
{
    // Absolute element rows bound the assembled lumped diffusion operator.
    // Half the forward-Euler spectral bound leaves a stability margin.
    const int mpi_ee = mpi_element_end(var);
#ifndef ACC
    #pragma omp parallel for default(none) shared(param, var, work) firstprivate(mpi_ee)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector firstprivate(mpi_ee)
    for (int e = 0; e < mpi_ee; ++e) {
        if (param.control.has_fluid_thermal_advection) {
            double a[NODES_PER_ELEM][NODES_PER_ELEM];
            thermal_transport_operator(param, var, e, a);
            for (int i = 0; i < NODES_PER_ELEM; ++i) {
                double row = 0;
                for (int j = 0; j < NODES_PER_ELEM; ++j) row += std::abs(a[i][j]);
                work[e][i] = row;
            }
            continue;
        }
        double g[NDIMS][NODES_PER_ELEM];
        get_shape_gradients(var, e, g);
        const double kv = param.control.has_thermal_diffusion ? var.mat->k(e)*(*var.volume)[e] : 0;
        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            double row = 0;
            for (int j = 0; j < NODES_PER_ELEM; ++j) {
                double a = 0;
                for (int d = 0; d < NDIMS; ++d) a += g[d][i] * g[d][j];
                row += std::abs(a);
            }
            work[e][i] = kv * row;
        }
    }
    // Ghost elements take their rows from the ranks that own them.
    mpi_sync(var, work, MpiEntity::Element, NODES_PER_ELEM);
    // Injected fluid relaxes the node towards its temperature at rate rho_f c_f Q_n.
    const double capacity = var.inflow_rate ? inflow_capacity(param) : 0.0;
    return nodal_step_limit(var, work, *var.tmass, std::numeric_limits<double>::max(),
                            var.inflow_rate, capacity);
}

double nodal_step_limit(const Variables &var, const elem_cache &work, const double_vec &mass,
                        double limit, const double_vec *inflow, double capacity)
{
    const int mpi_ne = mpi_node_end(var);
#ifndef ACC
    #pragma omp parallel for default(none) shared(var, work, mass) firstprivate(inflow, capacity, mpi_ne) reduction(min:limit)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector firstprivate(mpi_ne) reduction(min:limit)
    for (int n = 0; n < mpi_ne; ++n) {
        double row = 0;
        const int *patch = var.support.patch(n), *local = var.support.local(n);
        for (int k = 0; k < var.support.size(n); ++k) row += work[patch[k]][local[k]];
        if (inflow) row += capacity * (*inflow)[n];
        if (row > 0 && mass[n] > 0)
            limit = std::min(limit, mass[n] / row);
    }
    #pragma acc wait
#ifdef USEMPI
    if (var.mpi && var.mpi->in_kernel()) var.mpi->minimum(&limit, 1);
#endif
    return limit;
}

static void temperature_step(const Param &param, const Variables &var,
                             double_vec &temperature, elem_cache &tmp_result,
                             double dt, double time)
{
    const int mpi_eb = 0, mpi_ee = mpi_element_end(var);
    const int mpi_nb = 0, mpi_ne = mpi_node_end(var);
    if (param.heat_injection.enabled) {
        assemble_point_sources(param.heat_injection,var.heat_points,var,time-dt,time,*var.heat_source);
        for (int n=0; n<var.nnode; ++n) {
            double prescribed = 0;
            if ((*var.heat_source)[n] != 0 &&
                (prescribed_temperature(param,var,n,prescribed) ||
                 !((*var.tmass)[n]>0 && std::isfinite((*var.tmass)[n]))))
                die(EXIT_RUNTIME_LOOKUP, "A heat point source reaches a fixed-temperature node or nonpositive heat capacity.");
        }
    }
    const double capacity = var.inflow_rate ? inflow_capacity(param) : 0.0;
    if (var.inflow_rate) {
        assemble_inflow(param, var, time-dt, time);
        for (int n=0; n<var.nnode; ++n) {
            double prescribed = 0;
            if ((*var.inflow_rate)[n] != 0 && prescribed_temperature(param,var,n,prescribed))
                die(EXIT_RUNTIME_LOOKUP, "Fluid injected at a temperature reaches a fixed-temperature node.");
        }
    }
#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif
#ifndef ACC
    #pragma omp parallel default(none) shared(param,var,temperature,tmp_result) firstprivate(mpi_eb,mpi_ee,mpi_nb,mpi_ne,dt,capacity)
#endif
    {
#ifndef ACC
        #pragma omp for
#endif
        #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector async firstprivate(mpi_eb,mpi_ee,mpi_nb,mpi_ne)
        for (int e=mpi_eb; e<mpi_ee;e++) {
            // diffusion matrix

            ConstConnAccessor conn = (*var.connectivity)[e];
            ElemCacheAccessor tr = tmp_result[e];
            double kv = param.control.has_thermal_diffusion ? var.mat->k(e)*(*var.volume)[e] : 0;
            double rh = (*var.radiogenic_source)[e] * (*var.volume)[e] * var.mat->rho(e) / NODES_PER_ELEM;
            if (param.control.has_fluid_thermal_advection) rh *= 1-var.mat->phi(e);
            // Integrate the previous accepted energy, not its old power times
            // this step's possibly different dt. Each substep consumes its share.
            if (var.pending_heat) rh += (*var.pending_heat)[e] / (NODES_PER_ELEM * var.dt);
            double transport[NODES_PER_ELEM][NODES_PER_ELEM];
            if (param.control.has_fluid_thermal_advection)
                thermal_transport_operator(param, var, e, transport);
#ifdef THREED
            double shpdx[NODES_PER_ELEM], shpdy[NODES_PER_ELEM], shpdz[NODES_PER_ELEM];
            if (!param.control.has_fluid_thermal_advection)
                get_local_shape_fn(var, e, shpdx, shpdy, shpdz);
#else
            double shpdx[NODES_PER_ELEM], shpdz[NODES_PER_ELEM];
            if (!param.control.has_fluid_thermal_advection)
                get_local_shape_fn(var, e, shpdx, shpdz);
#endif
            for (int i=0; i<NODES_PER_ELEM; ++i) {
                if (param.control.has_fluid_thermal_advection) {
                    double value = -rh;
                    for (int j = 0; j < NODES_PER_ELEM; ++j)
                        value += transport[i][j]*temperature[conn[j]];
                    tr[i] = value;
                    continue;
                }
                double diffusion = 0.;
                for (int j=0; j<NODES_PER_ELEM; ++j) {
#ifdef THREED
                    diffusion += (shpdx[i] * shpdx[j] + \
                                shpdy[i] * shpdy[j] + \
                                shpdz[i] * shpdz[j]) * temperature[conn[j]];
#else
                    diffusion += (shpdx[i] * shpdx[j] + \
                                shpdz[i] * shpdz[j]) * temperature[conn[j]];
#endif
                }
                tr[i] = diffusion * kv - rh;
            }
        }

#if !defined(ACC) && defined(USEMPI)
        #pragma omp master
#endif
        {
            mpi_begin_gather(var,tmp_result,NODES_PER_ELEM);
        }
#if !defined(ACC) && defined(USEMPI)
        #pragma omp barrier
#endif
        const int passes=mpi_gather_passes(var);
        for(int pass=0;pass<passes;++pass) {
            if(pass==1) {
#if !defined(ACC) && defined(USEMPI)
                #pragma omp master
#endif
                mpi_finish_gather(var);
#if !defined(ACC) && defined(USEMPI)
                #pragma omp barrier
#endif
            }
            const MpiNodeRows rows=mpi_gather_nodes(var,pass);
#ifndef ACC
            #pragma omp for
#endif
            #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector async firstprivate(rows, dt, capacity)
            for(int row=0;row<rows.count;++row) {
                const int n=rows.ids?rows.ids[row]:row;
                double fixed_temperature = 0;
                if (prescribed_temperature(param, var, n, fixed_temperature))
                    temperature[n] = fixed_temperature;
                else {
                    double tdot = 0;

                    // Sum the element rows through the node's support patch (element, local index).
                    const int npatch = var.support.size(n);
                    const int* patch = var.support.patch(n);
                    const int* lpatch = var.support.local(n);
                    for (int k=0; k<npatch; ++k)
                        tdot += tmp_result[patch[k]][lpatch[k]];
                    if (param.heat_injection.enabled) tdot -= (*var.heat_source)[n];
                    // Injected enthalpy relative to the node: rho_f c_f (sum Q T_in - T sum Q).
                    if (var.inflow_rate)
                        tdot -= capacity * ((*var.inflow_heat)[n] - temperature[n] * (*var.inflow_rate)[n]);
                    temperature[n] -= dt * tdot / (*var.tmass)[n];
                }
            }
        }
    }

#ifdef NPROF
    nvtxRangePop();
#endif
}

int explicit_substeps(double dt, double limit, const char *message)
{
    const double count = std::ceil(dt / limit);
    if (!std::isfinite(count) || count > std::numeric_limits<int>::max())
        die(EXIT_RUNTIME_RESOURCE, message);
    return std::max(1, static_cast<int>(count));
}

void update_temperature(const Param &param, const Variables &var,
                        double_vec &temperature, elem_cache &work)
{
#ifdef USEMPI
    if (var.mpi && var.mpi->dispatch(MpiOp::Temperature, var, 0.0)) return;
#endif
    // A short pulse can be faster than its whole-step mean. Bound every
    // thermal substep by the active rates, then integrate its actual overlap.
    if (var.inflow_rate) assemble_inflow(param, var, var.time - var.dt, var.time, true);
    if (var.mat->thermal_step_limit == 0 || param.control.has_fluid_thermal_advection)
        var.mat->thermal_step_limit = temperature_step_limit(param, var, work);
    const int steps = explicit_substeps(var.dt, var.mat->thermal_step_limit,
                                        "Thermal timestep needs too many explicit substeps.");
    for (int i = 0; i < steps; ++i) {
        // Each substep reads its neighbours' ghost temperatures; the dispatcher exchanges the last.
        if (i > 0) mpi_sync(var, temperature, MpiEntity::Node);
        temperature_step(param, var, temperature, work, var.dt / steps,
            i+1 == steps ? var.time : var.time-var.dt+(i+1)*(var.dt/steps));
    }
    if (var.pending_heat) {
        #pragma acc wait
        std::fill(var.pending_heat->begin(), var.pending_heat->end(), 0.0);
    }
}

void accumulate_shear_heat(const Variables &var)
{
    if (!var.pending_heat) return;
    #pragma acc wait
    int invalid = 0;
#ifndef ACC
    #pragma omp parallel for default(none) shared(var) reduction(|:invalid)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector reduction(|:invalid)
    for (int e = 0; e < var.nelem; ++e) {
        const double energy = (*var.shear_heat)[e] * var.dt * (*var.volume)[e];
        (*var.pending_heat)[e] += energy;
        if (!(energy >= 0 && std::isfinite((*var.pending_heat)[e]))) invalid = 1;
    }
    if (invalid) die(EXIT_RUNTIME_NAN, "Non-finite or negative accepted shear heat.");
}

void update_strain_rate(const Variables& var, tensor_t& strain_rate)
{
#ifdef USEMPI
    if (var.mpi && var.mpi->dispatch(MpiOp::StrainRate, var, 0.0)) return;
#endif
    const int mpi_eb = 0, mpi_ee = mpi_element_end(var);
    const int mpi_nb = 0, mpi_ne = mpi_node_end(var);

#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif

#ifndef ACC
    #pragma omp parallel for default(none) shared(var, strain_rate) firstprivate(mpi_eb,mpi_ee,mpi_nb,mpi_ne)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector async firstprivate(mpi_eb,mpi_ee,mpi_nb,mpi_ne)
    for (int e=mpi_eb; e<mpi_ee; ++e) {
#ifdef THREED
        double shpdx[NODES_PER_ELEM], shpdy[NODES_PER_ELEM], shpdz[NODES_PER_ELEM];
        get_local_shape_fn(var, e, shpdx, shpdy, shpdz);
#else
        double shpdx[NODES_PER_ELEM], shpdz[NODES_PER_ELEM];
        get_local_shape_fn(var, e, shpdx, shpdz);
#endif
        TensorAccessor s = strain_rate[e];
        
        ConstArrayIndirectAccessor v = var.vel->view_const((*var.connectivity)[e]);

        // XX component
        int n = 0;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][0] * shpdx[i];

#ifdef THREED
        // YY component
        n = 1;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][1] * shpdy[i];
#endif

        // ZZ component
#ifdef THREED
        n = 2;
#else
        n = 1;
#endif
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][NDIMS-1] * shpdz[i];

#ifdef THREED
        // XY component
        n = 3;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][0] * shpdy[i] + v[i][1] * shpdx[i]);
#endif

        // XZ component
#ifdef THREED
        n = 4;
#else
        n = 2;
#endif
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][0] * shpdz[i] + v[i][NDIMS-1] * shpdx[i]);

#ifdef THREED
        // YZ component
        n = 5;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][1] * shpdz[i] + v[i][2] * shpdy[i]);
#endif
    }
#ifdef NPROF
    nvtxRangePop();
#endif
}


void apply_damping(const Param& param, const Variables& var, array_t& force)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    // flatten 2d arrays to simplify indexing
    const double small_vel = 1e-13;

    switch (param.control.damping_option) {
    case 0:
        // no damping, stress field can become very noisy
        break;
    case 1:
        // damping when force and velocity are parallel
        // acclerating when force and velocity are anti-parallel
#ifndef ACC
#ifdef GPP1X
        #pragma omp parallel for default(none) \
            shared(var, param, force, small_vel)
#else
        #pragma omp parallel for default(none) \
            shared(var, param, force)  firstprivate(small_vel)
#endif
#endif
        #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector collapse(2) async
        for (int i=0; i<var.nnode; ++i) {
            for (int j=0;j<NDIMS;j++)
                if (std::fabs((*var.vel)[i][j]) > small_vel) {
                    force[i][j] -= param.control.damping_factor * std::copysign(force[i][j], (*var.vel)[i][j]);
                }
        }
        break;
    case 2:
        // damping prop. to force
#ifndef ACC
        #pragma omp parallel for default(none) \
            shared(var, param, force)
#endif
        #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector collapse(2) async
        for (int i=0; i<var.nnode; ++i) {
            for (int j=0;j<NDIMS;j++)
                force[i][j] -= param.control.damping_factor * force[i][j];
        }
        break;
    case 3:
        // damping when force and velocity are parallel
        // weakly acclerating when force and velocity are anti-parallel
#ifndef ACC
        #pragma omp parallel for default(none) shared(var, param, force)
#endif
        #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector collapse(2) async
        for (int i=0; i<var.nnode; ++i) {
            for (int j=0;j<NDIMS;j++) {
                if ((force[i][j]<0) == ((*var.vel)[i][j]<0)) {
                    // strong damping
                    force[i][j] -= param.control.damping_factor * force[i][j];
                }
                else {
                    // weak acceleration
                    force[i][j] += (1 - param.control.damping_factor) * force[i][j];
                }
            }
        }
        break;
    case 4:
        // rayleigh damping
#ifndef ACC
#ifdef GPP1X
        #pragma omp parallel for default(none) \
            shared(var, param, force, small_vel)
#else
        #pragma omp parallel for default(none) \
            shared(var, param, force)  firstprivate(small_vel)
#endif
#endif
        #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector async
        for (int i=0; i<var.nnode; ++i) {
            double mass = (*var.mass)[i];
            double young = (*var.ymass)[i];
            double critical_coeff = 2.0 * std::sqrt(mass * young);
            for (int j=0;j<NDIMS;j++)
                if (std::fabs((*var.vel)[i][j]) > small_vel) {
                    double f_C = param.control.damping_factor * std::copysign(force[i][j], (*var.vel)[i][j]);
                    double f_V = critical_coeff * (*var.vel)[i][j];
                    double f_damping = (std::fabs(f_C) < std::fabs(f_V)) ? f_V : f_C;
                    force[i][j] -= f_damping;
                }
        }
        break;
    default:
        if (mpi_reports(var))
            std::cerr << "Error: unknown damping_option: " << param.control.damping_option << '\n';
        die_together(EXIT_CONFIG_VALUE);
    }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

/*
#pragma acc routine seq
static double rho(const conn_t &var_connectivity, \
    const double_vec &var_temperature, const int_vec2D &var_elemmarkers, \
    const double_vec &rho0, const double_vec &alpha, int nmat,  int e)
{
    const double celsius0 = 273;

    // average temperature of this element
    double T = 0;
    const int *conn = var_connectivity[e];
    for (int i=0; i<NODES_PER_ELEM; ++i) {
        T += var_temperature[conn[i]];
    }
    T /= NODES_PER_ELEM;

    double TinCelsius = T - celsius0;
    double result = 0;
    int n = 0;
    for (int m=0; m<nmat; m++) {
        // TODO: compressibility
        result += rho0[m] * (1 - alpha[m] * TinCelsius) * var_elemmarkers[e][m];
        n += var_elemmarkers[e][m];
    }
    return result / n;
}
*/

void update_force(const Param& param, const Variables& var, array_t& force, array_t& force_residual,
                  elem_cache& tmp_result, array_t* force_undamped, double displacement_dt, bool damp_force,
                  array_t* node_force_magnitude)
{
#ifdef USEMPI
    if (var.mpi && var.mpi->dispatch(MpiOp::Force, var, 0.0)) return;
#endif
    const int mpi_eb = 0, mpi_ee = mpi_element_end(var);
    const int mpi_nb = 0, mpi_ne = mpi_node_end(var);

#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif

#ifndef ACC
    #pragma omp parallel default(none) shared(var,param,force,force_residual,tmp_result,node_force_magnitude) firstprivate(mpi_eb,mpi_ee,mpi_nb,mpi_ne)
#endif
    {
#ifndef ACC
        #pragma omp for
#endif
        #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector async firstprivate(mpi_eb,mpi_ee,mpi_nb,mpi_ne)
        for (int e=mpi_eb; e<mpi_ee;e++) {
#ifdef THREED
            double shpdx[NODES_PER_ELEM], shpdy[NODES_PER_ELEM], shpdz[NODES_PER_ELEM];
            get_local_shape_fn(var, e, shpdx, shpdy, shpdz);
#else
            double shpdx[NODES_PER_ELEM], shpdz[NODES_PER_ELEM];
            get_local_shape_fn(var, e, shpdx, shpdz);
#endif
            TensorAccessor s = (*var.stress)[e];
            double vol = (*var.volume)[e];
            ElemCacheAccessor tr = tmp_result[e];

            double buoy = 0;
            if (param.control.gravity != 0) {
                // Unsaturated pores carry only S of their water weight.
                double sat = 1;
                if (param.control.has_hydraulic_diffusion && var.mat->has_retention)
                    sat = element_mean_saturation(*var.mat, (*var.elemmarkers)[e],
                                                  (*var.connectivity)[e], *var.ppressure);
                const double density = param.control.has_hydraulic_diffusion
                    ? hydraulic_bulk_density(param,var,e,sat) : var.mat->rho(e);
                buoy = density * param.control.gravity / NODES_PER_ELEM;
            }


            for (int i=0; i<NODES_PER_ELEM; ++i) {
#ifdef THREED
                tr[i] = (s[0]*shpdx[i] + s[3]*shpdy[i] + s[4]*shpdz[i]) * vol;
                tr[i+NODES_PER_ELEM] = (s[3]*shpdx[i] + s[1]*shpdy[i] + s[5]*shpdz[i]) * vol;
                tr[i+NODES_PER_ELEM*2] = (s[4]*shpdx[i] + s[5]*shpdy[i] + s[2]*shpdz[i] + buoy) * vol;
#else
                tr[i] = (s[0]*shpdx[i] + s[2]*shpdz[i]) * vol;
                tr[i+NODES_PER_ELEM] = (s[2]*shpdx[i] + s[1]*shpdz[i] + buoy) * vol;
#endif
            }
        }

#if !defined(ACC) && defined(USEMPI)
        #pragma omp master
#endif
        {
            mpi_begin_gather(var,tmp_result,NODES_PER_ELEM*NDIMS);
        }
#if !defined(ACC) && defined(USEMPI)
        #pragma omp barrier
#endif
        const int passes=mpi_gather_passes(var);
        for(int pass=0;pass<passes;++pass) {
            if(pass==1) {
#if !defined(ACC) && defined(USEMPI)
                #pragma omp master
#endif
                mpi_finish_gather(var);
#if !defined(ACC) && defined(USEMPI)
                #pragma omp barrier
#endif
            }
            const MpiNodeRows rows=mpi_gather_nodes(var,pass);
#ifndef ACC
            #pragma omp for
#endif
            #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector async firstprivate(rows)
            for(int row=0;row<rows.count;++row) {
                const int n=rows.ids?rows.ids[row]:row;
                double magnitude[NDIMS] = {};
                ArrayAccessor f = force[n];
                f = 0;
                ArrayAccessor f_residual = force_residual[n];
                f_residual = 0;
                const int npatch = var.support.size(n);
                const int* patch = var.support.patch(n);
                const int* lpatch = var.support.local(n);
                for (int k=0;k<npatch;++k) {
                    ConstElemCacheAccessor tr = tmp_result[patch[k]];
                    const int i = lpatch[k];
                    for (int j=0;j<NDIMS;j++)
                    {
                        if (node_force_magnitude) magnitude[j] += std::abs(tr[i+NODES_PER_ELEM*j]);
                        f[j] -= tr[i+NODES_PER_ELEM*j];
                        f_residual[j] = tr[i+NODES_PER_ELEM*j];
                    }
                }
                if (node_force_magnitude)
                    for (int d=0; d<NDIMS; ++d) (*node_force_magnitude)[n][d] = magnitude[d];
            }
        }
    }

    apply_stress_bcs(param, var, force, displacement_dt);

    if (!param.ic.has_body_force_adjustment) apply_stress_bcs_neumann(param, var, force);
    apply_platen_force(param, var, force);
    apply_platen_force(param, var, force_residual);
    if (force_undamped) {
        // Same accelerator queue as assembly and damping; the wait below makes
        // both outputs visible to the caller. Constraint projection is separate.
#ifndef ACC
        #pragma omp parallel for default(none) shared(var, force, force_undamped) firstprivate(mpi_ne)
#endif
        #pragma acc parallel loop gang vector async firstprivate(mpi_ne)
        for (int n=0; n<mpi_ne; ++n)
            for (int d=0; d<NDIMS; ++d)
                (*force_undamped)[n][d] = force[n][d];
    }
    if (damp_force) {
        apply_damping(param, var, force);
        apply_platen_force(param, var, force);
    }

    #pragma acc wait

#ifdef NPROF
    nvtxRangePop();
#endif
}

double calculate_residual_force(const Variables& var, array_t& force_residual)
{
#ifdef USEMPI
    if (var.mpi && var.mpi->dispatch(MpiOp::Residual, var, 0.0)) return var.mpi->result();
#endif
    const int mpi_eb = 0, mpi_ee = mpi_element_end(var);
    const int mpi_nb = 0, mpi_ne = mpi_node_end(var);

#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    double l2 = 0.0;
    double num = var.nnode * NDIMS;
#ifdef USEMPI
    if (var.mpi && var.mpi->in_kernel()) num = var.mpi->global_nodes() * NDIMS;
#endif

#ifndef ACC
    #pragma omp parallel for default(none) shared(var, force_residual, num) reduction(+:l2) firstprivate(mpi_eb,mpi_ee,mpi_nb,mpi_ne)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector collapse(2) reduction(+:l2) async firstprivate(mpi_eb,mpi_ee,mpi_nb,mpi_ne)
    for (int i=mpi_nb; i<mpi_ne; ++i)
        for (int j = 0; j < NDIMS; ++j)
            l2 += std::pow(force_residual[i][j], 2) / num;

    #pragma acc wait

#ifdef NPROF
    nvtxRangePop();
#endif
#ifdef USEMPI
    if (var.mpi && var.mpi->in_kernel()) l2 = var.mpi->sum(l2);
#endif
    return std::sqrt(l2);
}


void update_velocity(const Variables& var, array_t& vel)
{
#ifdef USEMPI
    if (var.mpi && var.mpi->dispatch(MpiOp::Velocity, var, 0.0)) return;
#endif
    const int mpi_eb = 0, mpi_ee = mpi_element_end(var);
    const int mpi_nb = 0, mpi_ne = mpi_node_end(var);

#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif

#ifndef ACC
    #pragma omp parallel for default(none) shared(var, vel) firstprivate(mpi_eb,mpi_ee,mpi_nb,mpi_ne)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector collapse(2) async firstprivate(mpi_eb,mpi_ee,mpi_nb,mpi_ne)
    for (int i=mpi_nb; i<mpi_ne; ++i)
        for (int j=0;j<NDIMS;j++)
            vel[i][j] += var.dt * (*var.force)[i][j] / (*var.mass)[i];

#ifdef NPROF
    nvtxRangePop();
#endif
}

void update_coordinate(const Variables& var, array_t& coord)
{
#ifdef USEMPI
    if (var.mpi && var.mpi->dispatch(MpiOp::Coordinate, var, 0.0)) return;
#endif
    const int mpi_eb = 0, mpi_ee = mpi_element_end(var);
    const int mpi_nb = 0, mpi_ne = mpi_node_end(var);

#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    // double* x = var.coord->data();
    // const double* v = var.vel->data();

    // for gpu parallelization dt and bound need to be sent to 

#ifndef ACC
    #pragma omp parallel for default(none) shared(var) firstprivate(mpi_eb,mpi_ee,mpi_nb,mpi_ne)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector collapse(2) async firstprivate(mpi_eb,mpi_ee,mpi_nb,mpi_ne)
    for (int i=mpi_nb; i<mpi_ne; ++i) {
        for (int j=0 ; j<NDIMS; ++j){
            (*var.coord)[i][j] += (*var.vel)[i][j] * var.dt;
        }
        
    }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}


namespace {

#ifdef THREED
    #pragma acc routine seq
    void jaumann_rate_3d(TensorAccessor s, double dt, double w3, double w4, double w5)
    {
        double s_inc[NSTR];

        s_inc[0] =-2.0 * s[3] * w3 - 2.0 * s[4] * w4;
        s_inc[1] = 2.0 * s[3] * w3 - 2.0 * s[5] * w5;
        s_inc[2] = 2.0 * s[4] * w4 + 2.0 * s[5] * w5;
        s_inc[3] = s[0] * w3 - s[1] * w3 - s[4] * w5 - s[5] * w4;
        s_inc[4] = s[0] * w4 - s[2] * w4 + s[3] * w5 - s[5] * w3;
        s_inc[5] = s[1] * w5 - s[2] * w5 + s[3] * w4 + s[4] * w3;

        for(int i=0; i<NSTR; ++i)  {
            s[i] += dt * s_inc[i];
        }
    }

#else
    #pragma acc routine seq
    void jaumann_rate_2d(TensorAccessor s, double dt, double w2)
    {
        double s_inc[NSTR];

        s_inc[0] =-2.0 * s[2] * w2;
        s_inc[1] = 2.0 * s[2] * w2;
        s_inc[2] = s[0] * w2 - s[1] * w2;

        for(int i=0; i<NSTR; ++i)  {
            s[i] += dt * s_inc[i];
        }
    }

#endif

}


void rotate_stress(const Variables &var, tensor_t &stress, tensor_t &strain,
                   double dt)
{
#ifdef USEMPI
    if (var.mpi && var.mpi->dispatch(MpiOp::Rotation, var, dt)) return;
#endif
    const int mpi_eb = 0, mpi_ee = mpi_element_end(var);
    const int mpi_nb = 0, mpi_ne = mpi_node_end(var);

#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    // The spin rate tensor, W, and the Cauchy stress tensor, S, are
    //     [  0  w3  w4]     [s0 s3 s4]
    // W = [-w3   0  w5],  S=[s3 s1 s5].
    //     [-w4 -w5   0]     [s4 s5 s2]
    //
    // Stress (and strain) increment based on the Jaumann rate is
    // dt*(S*W-W*S). S*W-W*S is also symmetric.
    // So, following the indexing of stress tensor,
    // we get
    // sj[0] = dt * ( -2 * s3 * w3 - 2 * s4 * w4)
    // sj[1] = dt * (  2 * s3 * w3 - 2 * s5 * w5)
    // sj[2] = dt * (  2 * s4 * w4 + 2 * s5 * w5)
    // sj[3] = dt * ( s0 * w3 - s1 * w3 - s4 * w5 - s5 * w4)
    // sj[4] = dt * ( s0 * w4 - s2 * w4 + s3 * w5 - s5 * w3)
    // sj[5] = dt * ( s1 * w5 - s2 * w5 + s3 * w4 + s4 * w3)

#ifndef ACC
    #pragma omp parallel for default(none) shared(var, stress, strain) \
        firstprivate(dt) firstprivate(mpi_eb,mpi_ee,mpi_nb,mpi_ne)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector firstprivate(dt) async firstprivate(mpi_eb,mpi_ee,mpi_nb,mpi_ne)
    for (int e=mpi_eb; e<mpi_ee; ++e) {
        ConstConnAccessor conn = (*var.connectivity)[e];

#ifdef THREED

        double w3, w4, w5;
        {
            double shpdx[NODES_PER_ELEM], shpdy[NODES_PER_ELEM], shpdz[NODES_PER_ELEM];
            get_local_shape_fn(var, e, shpdx, shpdy, shpdz);

            ConstArrayIndirectAccessor v = var.vel->view_const(conn);

            w3 = 0;
            for (int i=0; i<NODES_PER_ELEM; ++i)
                w3 += 0.5 * (v[i][0] * shpdy[i] - v[i][1] * shpdx[i]);

            w4 = 0;
            for (int i=0; i<NODES_PER_ELEM; ++i)
                w4 += 0.5 * (v[i][0] * shpdz[i] - v[i][NDIMS-1] * shpdx[i]);

            w5 = 0;
            for (int i=0; i<NODES_PER_ELEM; ++i)
                w5 += 0.5 * (v[i][1] * shpdz[i] - v[i][NDIMS-1] * shpdy[i]);
        }

        jaumann_rate_3d(stress[e], dt, w3, w4, w5);
        jaumann_rate_3d(strain[e], dt, w3, w4, w5);

#else

        double w2;
        {
            double shpdx[NODES_PER_ELEM], shpdz[NODES_PER_ELEM];
            get_local_shape_fn(var, e, shpdx, shpdz);

            ConstArrayIndirectAccessor v = var.vel->view_const(conn);

            w2 = 0;
            for (int i=0; i<NODES_PER_ELEM; ++i)
                w2 += 0.5 * (v[i][NDIMS-1] * shpdx[i] - v[i][0] * shpdz[i]);
        }

        jaumann_rate_2d(stress[e], dt, w2);
        jaumann_rate_2d(strain[e], dt, w2);

#endif
    }
#ifdef NPROF
    nvtxRangePop();
#endif
}
