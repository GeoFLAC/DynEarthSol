#include <iostream>
#include <limits>
#include <cstdlib>
#include <fstream>


#include "constants.hpp"
#include "parameters.hpp"
#include "ats_output_scheduler.hpp"
#include "bc.hpp"
#include "binaryio.hpp"
#include "earthquake_state.hpp"
#include "fields.hpp"
#include "geometry.hpp"
#include "ic.hpp"
#include "input.hpp"
#include "matprops.hpp"
#include "markerset.hpp"
#include "mesh.hpp"
#include "monitor.hpp"
#include "output.hpp"
#include "phasechanges.hpp"
#include "remeshing.hpp"
#include "rheology.hpp"
#include "runtime_info.hpp"
#include "utils.hpp"
#include "gospl_driver/gospl-driver.hpp"

#ifdef WIN32
#ifdef _MSC_VER
#define snprintf _snprintf
#endif // _MSC_VER
namespace std { using ::snprintf; }
#endif // WIN32

void init_var(const Param& param, Variables& var)
{
    var.time = 0;
    var.last_remesh_time = 0;
    var.steps = 0;
    var.nremesh = 0;
    var.nremesh_prior = 0;
    var.noutput = 0;
    var.nphasechange = 0;
    var.func_time.output_time = 0;
    var.func_time.remesh_time = 0;
    var.func_time.phase_change_time = 0;
    var.func_time.start_time = get_nanoseconds();
    var.init_elem_size_n = new double_vec(0);
    var.refined_elem_size_n = new double_vec(0);
    var.reference_frame_time = 0.0;

    for (int i=0;i<nbdrytypes;++i)
        var.bfacets[i] = new int_pair_vec;
    for (int i=0;i<nbdrytypes;++i)
        var.bnodes[i] = new int_vec;
    var.bnormals = new array_t(nbdrytypes);

    if (param.control.characteristic_speed == 0)
        var.max_vbc_val = find_max_vbc(param.bc);
        // todo max_vbc_val change with boundary period.
    else
        var.max_vbc_val = param.control.characteristic_speed;

    // XXX: Hard coded boundary flag. If the order of ibound?? is changed
    //      in the future, the following lines have to be updated as well.
    var.vbc_types[0] = param.bc.vbc_x0;
    var.vbc_types[1] = param.bc.vbc_x1;
    var.vbc_types[2] = param.bc.vbc_y0;
    var.vbc_types[3] = param.bc.vbc_y1;
    var.vbc_types[4] = param.bc.vbc_z0;
    var.vbc_types[5] = param.bc.vbc_z1;
    var.vbc_types[6] = param.bc.vbc_n0;
    var.vbc_types[7] = param.bc.vbc_n1;
    var.vbc_types[8] = param.bc.vbc_n2;
    var.vbc_types[9] = param.bc.vbc_n3;

    var.vbc_values[0] = param.bc.vbc_val_x0;
    var.vbc_values[1] = param.bc.vbc_val_x1;
    var.vbc_values[2] = param.bc.vbc_val_y0;
    var.vbc_values[3] = param.bc.vbc_val_y1;
    var.vbc_values[4] = param.bc.vbc_val_z0;
    var.vbc_values[5] = param.bc.vbc_val_z1;
    var.vbc_values[6] = param.bc.vbc_val_n0;
    var.vbc_values[7] = param.bc.vbc_val_n1;
    var.vbc_values[8] = param.bc.vbc_val_n2;
    var.vbc_values[9] = param.bc.vbc_val_n3;

    var.vbc_vertical_div_x0[0] = 0.;
    var.vbc_vertical_div_x0[1] = param.bc.vbc_val_division_x0_min;
    var.vbc_vertical_div_x0[2] = param.bc.vbc_val_division_x0_max;
    var.vbc_vertical_div_x0[3] = 1.;
    var.vbc_vertical_div_x1[0] = 0.;
    var.vbc_vertical_div_x1[1] = param.bc.vbc_val_division_x1_min;
    var.vbc_vertical_div_x1[2] = param.bc.vbc_val_division_x1_max;
    var.vbc_vertical_div_x1[3] = 1.;

    var.vbc_vertical_ratio_x0[0] = param.bc.vbc_val_x0_ratio0;
    var.vbc_vertical_ratio_x0[1] = param.bc.vbc_val_x0_ratio1;
    var.vbc_vertical_ratio_x0[2] = param.bc.vbc_val_x0_ratio2;
    var.vbc_vertical_ratio_x0[3] = param.bc.vbc_val_x0_ratio3;
    var.vbc_vertical_ratio_x1[0] = param.bc.vbc_val_x1_ratio0;
    var.vbc_vertical_ratio_x1[1] = param.bc.vbc_val_x1_ratio1;
    var.vbc_vertical_ratio_x1[2] = param.bc.vbc_val_x1_ratio2;
    var.vbc_vertical_ratio_x1[3] = param.bc.vbc_val_x1_ratio3;

    var.hbc_types[0] = param.bc.hbc_x0;
    var.hbc_types[1] = param.bc.hbc_x1;
    var.hbc_types[2] = param.bc.hbc_y0;
    var.hbc_types[3] = param.bc.hbc_y1;
    var.hbc_types[4] = param.bc.hbc_z0;
    var.hbc_types[5] = param.bc.hbc_z1;

    var.stress_bc_types[0] = param.bc.stress_bc_x0;
    var.stress_bc_types[1] = param.bc.stress_bc_x1;
    var.stress_bc_types[2] = param.bc.stress_bc_y0;
    var.stress_bc_types[3] = param.bc.stress_bc_y1;
    var.stress_bc_types[4] = param.bc.stress_bc_z0;
    var.stress_bc_types[5] = param.bc.stress_bc_z1;

    var.stress_bc_values[0] = param.bc.stress_val_x0;
    var.stress_bc_values[1] = param.bc.stress_val_x1;
    var.stress_bc_values[2] = param.bc.stress_val_y0;
    var.stress_bc_values[3] = param.bc.stress_val_y1;
    var.stress_bc_values[4] = param.bc.stress_val_z0;
    var.stress_bc_values[5] = param.bc.stress_val_z1;

    var.vbc_val_z1_loading_period = param.bc.vbc_val_z1_loading_period;

    // for regular mesh
    var.nx = std::round(param.mesh.xlength/param.mesh.resolution) + 1;
    var.nz = std::round(param.mesh.zlength/param.mesh.resolution) + 1;
    var.ncell = (var.nx-1) * (var.nz-1);
    var.nnode = var.nx * var.nz;
    var.nelem = 2 * var.ncell;
    var.nseg = 2 * (var.nx + var.nz - 2);
#ifdef THREED
    var.ny = std::round(param.mesh.ylength/param.mesh.resolution) + 1;
    var.ncell *= (var.ny-1);
    var.nnode *= var.ny;
    var.nelem = 5 * var.ncell;
    var.nseg = 4 * ( (var.nx-1) * (var.ny-1) + \
                     (var.ny-1) * (var.nz-1) + \
                     (var.nz-1) * (var.nx-1) );
#endif

#ifdef NPROF
    var.log_table = new double_vec(LOG_TABLE_SIZE);
    for (int i=0; i<LOG_TABLE_SIZE; ++i)
        (*var.log_table)[i] = std::log(LOG_XDELTA * i + LOG_XMIN);

    var.tan_table = new double_vec(TAN_TABLE_SIZE);
    for (int i=0; i<TAN_TABLE_SIZE; ++i)
        (*var.tan_table)[i] = std::tan(TAN_XDELTA * i + TAN_XMIN);

    var.sin_table = new double_vec(SIN_TABLE_SIZE);
    for (int i=0; i<SIN_TABLE_SIZE; ++i)
        (*var.sin_table)[i] = std::sin(SIN_XDELTA * i + SIN_XMIN);
#endif
}


void init(const Param& param, Variables& var)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    std::cout << "Initializing mesh and field data...\n";

    create_new_mesh(param, var);
    create_boundary_flags(var);
    create_boundary_nodes(var);
    create_boundary_facets(var);
    create_support(var);
    // create_neighbor(var);
    create_elemmarkers(param, var);
    create_markers(param, var);

    allocate_variables(param, var);

    create_top_elems(var);
    create_surface_info(param,var,var.surfinfo);

    for(int i=0; i<var.nnode; i++)
        for(int d=0; d<NDIMS; d++)
            (*var.coord0)[i][d] = (*var.coord)[i][d];

    compute_volume(*var.coord, *var.connectivity, *var.volume);

    #pragma acc wait

    *var.volume_old = *var.volume;
    // Must precede apply_vbcs: it is the only writer of edge_vec/edge_slot, which
    // apply_vbcs reads for any node on two boundaries at once.
    create_boundary_normals(var, *var.bnormals, var.edge_vec, var.edge_slot);
    apply_vbcs(param, var, *var.vel); // Global-velocity scaling needs boundary conditions before compute_mass.
    // The RSF state variable is initialized below; skip its bound for this
    // bootstrap call. main() recomputes dt after state initialization.
    var.dt = compute_dt(param, var, false);
    // compute_mass runs at the END of init(): it reads rho(e), a function of var.temperature
    // and elemmarkers, which the initial_* calls below are the first to write.




    // temperature should be init'd before stress and strain
    initial_temperature(param, var, *var.temperature, *var.radiogenic_source, var.bottom_temperature, *var.markersets[0], *var.elemmarkers, *var.markers_in_elem);
    // initial_temperature() reassigns mantle to asthenosphere, moving elemmarkers.
    var.mat->refresh_elem_cache();
    initial_stress_state(param, var, *var.stress, *var.stressyy, *var.old_mean_stress, *var.strain, var.compensation_pressure);
    // initial_stress_state_1d_load(param, var, *var.stress, *var.stressyy, *var.old_mean_stress, *var.strain, var.compensation_pressure);
    if(param.control.has_hydraulic_diffusion)
        initial_hydrostatic_state(param, var, *var.ppressure, *var.dppressure);

    initial_weak_zone(param, var, *var.plstrain);
    *var.plstrain_remesh = *var.plstrain;   // a prescribed weak zone is not R2 growth
    if (param.mat.rheol_type & MatProps::rh_rsf) {
        initial_friction_coeff(param, var, *var.dyn_fric_coeff);
        initial_state_variable(param, var, *var.state_variable);
    }

    // init viscosity
    #pragma omp parallel for default(none) shared(var)
    for (int e=0; e<var.nelem; ++e)
        (*var.viscosity)[e] = var.mat->visc(e);

    // Must run after the initial_* block and refresh_elem_cache(): rho(e) reads var.temperature
    // and elemmarkers, which those calls are the first to write -- before them tmass was built at
    // T = 0 K, too large by rho0*(1 + 273*alpha)/rho(e) (5.4 % on hot silicate). On a static mesh
    // the time loop gates update_mesh and remesh on has_moving_mesh, so only its rate-and-state
    // branch and its rebuild after phase changes rebuild the masses; isostasy_adjustment() calls
    // update_mesh() ungated before the loop.
    compute_mass(param, var, var.max_vbc_val, *var.volume_n, *var.mass, *var.tmass, *var.hmass, *var.ymass, *var.tmp_result);

#ifdef USEMMG
    // After compute_mass: it fills var.volume_n, which this divides by.
    initialize_elem_size_n(var, *var.init_elem_size_n);
    // No refinement reached yet; the wait covers initialize_elem_size_n's async node kernel.
    #pragma acc wait
    *var.refined_elem_size_n = *var.init_elem_size_n;
#endif

    report_mesh_info(var, "initial");

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}


void restart(const Param& param, Variables& var)
{
    std::cout << "Initializing mesh and field data from checkpoints...\n";

    char filename_save[256];
#ifdef HDF5
    std::snprintf(filename_save, 255, "%s.save.%06d.vtkhdf",
                  param.sim.restarting_from_modelname.c_str(), param.sim.restarting_from_frame);
    HDF5Input bin_save(filename_save);
#else
    std::snprintf(filename_save, 255, "%s.save.%06d",
                  param.sim.restarting_from_modelname.c_str(), param.sim.restarting_from_frame);
    BinaryInput bin_save(filename_save);
#endif
    std::cout << "  Reading " << filename_save << "...\n";

    char filename[256];
    std::snprintf(filename, 255, "%s.info", param.sim.restarting_from_modelname.c_str());
    bool got_meta = false;
    std::FILE *f = std::fopen(filename, "r");
    if (f != NULL) {
        int frame, steps, nnode, nelem, nseg;
        while (1) {
            int n = std::fscanf(f, "%d %d %*f %*f %*f %d %d %d\n",
                                &frame, &steps, &nnode, &nelem, &nseg);
            if (n != 5) break;   // EOF or malformed line: stop scanning, fall back
            if (frame == param.sim.restarting_from_frame) {
                var.steps = steps;
                var.nnode = nnode;
                var.nelem = nelem;
                var.nseg = nseg;
                got_meta = true;
                break;
            }
        }
        std::fclose(f);
        if (!got_meta) {
            std::cerr << "Warning: frame " << param.sim.restarting_from_frame
                      << " not found in " << filename << "; using metadata embedded in "
                      << filename_save << ".\n";
        }
    } else {
        std::cerr << "Warning: cannot open info file " << filename
                << "; using metadata embedded in " << filename_save << ".\n";
    }

    if (!got_meta) {
        if (bin_save.has_array("steps") && bin_save.has_array("nseg")) {
            bin_save.read_scalar(var.steps, "steps");
            bin_save.read_scalar(var.nnode, "nnode");
            bin_save.read_scalar(var.nelem, "nelem");
            bin_save.read_scalar(var.nseg, "nseg");
        } else {
            std::cerr << "Error: cannot read frame metadata from " << filename
                      << " and " << filename_save << " has none embedded.\n";
            die(EXIT_IO_RESTART);
        }
    }

    char filename_chkpt[256];
#ifdef HDF5
    std::snprintf(filename_chkpt, 255, "%s.chkpt.%06d.vtkhdf",
                  param.sim.restarting_from_modelname.c_str(), param.sim.restarting_from_frame);
    HDF5Input bin_chkpt(filename_chkpt);
#else
    std::snprintf(filename_chkpt, 255, "%s.chkpt.%06d",
                  param.sim.restarting_from_modelname.c_str(), param.sim.restarting_from_frame);
    BinaryInput bin_chkpt(filename_chkpt);
#endif
    std::cout << "  Reading " << filename_chkpt << "...\n";

    //
    // Following the same procedure in init()
    //

    // Reading mesh, replacing create_new_mesh()
    {
        var.coord = new array_t(var.nnode);
        bin_save.read_array(*var.coord, "coordinate");
        var.connectivity = new conn_t(var.nelem);
        bin_save.read_array(*var.connectivity, "connectivity");

        var.segment = new segment_t(var.nseg);
        bin_chkpt.read_array(*var.segment, "segment");
        var.segflag = new segflag_t(var.nseg);
        bin_chkpt.read_array(*var.segflag, "segflag");
        // Note: regattr is not needed for restarting
        // var.regattr = new regattr_t(var.nelem);
        // bin_chkpt.read_array(*var.regattr, "regattr", var.nelem);
    }

    create_boundary_flags(var);
    create_boundary_nodes(var);
    create_boundary_facets(var);
    create_support(var);
    // create_neighbor(var);
    create_elemmarkers(param, var);

    // Replacing create_markers()
    var.markersets.push_back(new MarkerSet(param, var, bin_save, bin_chkpt, std::string("markerset")));
    if (param.control.has_hydration_processes) {
        var.hydrous_marker_index = var.markersets.size();
        var.markersets.push_back(new MarkerSet(param, var, bin_save, bin_chkpt, std::string("hydrous-markerset")));
    }

    allocate_variables(param, var);
#ifdef USEMMG
    var.init_elem_size_n->resize(var.nnode);
    var.refined_elem_size_n->resize(var.nnode);
#endif

    create_top_elems(var);
    create_surface_info(param,var,var.surfinfo);

    bin_save.read_array(*var.coord0, "coord0");

    // Misc. items
    {
        bin_chkpt.read_scalar(var.time, "time");
        bin_chkpt.read_scalar(var.info_display_next_step, "info_display_next_step");
        bin_chkpt.read_scalar(var.compensation_pressure, "compensation_pressure");
        bin_chkpt.read_scalar(var.bottom_temperature, "bottom_temperature");
        bin_chkpt.read_scalar(var.dt, "dt");
        bin_chkpt.read_scalar(var.max_global_vel_mag, "max_global_vel_mag");
        bin_chkpt.read_scalar(var.reference_frame_time, "reference_frame_time");
        bin_chkpt.read_scalar(var.last_remesh_time, "last_remesh_time");
        bin_chkpt.read_scalar(var.nremesh_prior, "nremesh");
    }

    if (var.steps % param.mesh.quality_check_step_interval == 0 &&
        var.steps >= var.info_display_next_step)
        var.info_display_next_step = var.steps + param.sim.info_display_step_interval;

    // Restore the cached marker-cell index (mat.phase_change_option 3, see
    // VirtualTopology in mesh.hpp) so the restart carries forward the EXACT
    // lagged binning the checkpoint captured, instead of a fresher rebuild.
    // Only when uses_serpentinization_grid(), the one case that builds it: a
    // checkpoint restarted with another option, or with serpentinization off,
    // leaves var.vtopo null, so the hooks and checkpoints do not carry a dead
    // index. has_array-gated: a checkpoint written by a run without the grid
    // has no vtopo arrays -> var.vtopo stays null and the first
    // phase_changes rebuilds it lazily.
    if (uses_serpentinization_grid(param) && bin_chkpt.has_array("vtopo_int")) {
        // read_array does not auto-resize; pre-size from the known layout.
        // vtopo_int is variable length (VirtualTopology::serialize()), so
        // vtopo_len records that length up front.
        int_vec vt_len(1);
        bin_chkpt.read_array(vt_len, "vtopo_len");
        if (vt_len[0] < 1)
            die(EXIT_IO_RESTART, "vtopo_len is not a length");
        int_vec vt_i(vt_len[0]);
        bin_chkpt.read_array(vt_i, "vtopo_int");
        // vtopo_dbl is variable length too (VirtualTopology::serialize()),
        // recorded in vtopo_dbl_len like vtopo_len.
        int_vec vt_dlen(1);
        bin_chkpt.read_array(vt_dlen, "vtopo_dbl_len");
        if (vt_dlen[0] < 1)
            die(EXIT_IO_RESTART, "vtopo_dbl_len is not a length");
        double_vec vt_d(vt_dlen[0]);
        bin_chkpt.read_array(vt_d, "vtopo_dbl");
        delete var.vtopo;
        var.vtopo = VirtualTopology::deserialize(vt_i, vt_d);
        // The input checks validated the cfg's cell size, not the checkpoint's:
        // a restart that changed mesh.resolution or mesh.vt_size drops the
        // restored grid, and the next phase_changes or remesh builds one at the
        // cfg's size. The run then no longer continues the checkpointed one exactly.
        // A restart that raises mat.serpentinization_thickness keeps the
        // checkpoint's grid bottom until the next rebuild, at most one quality
        // check later, so a look-down near the depth cap may stop at its last row.
        const double cell_size = param.mesh.resolution * param.mesh.vt_size;
        if (var.vtopo->cell_size() != cell_size) {
            std::cerr << "Warning: the checkpoint's serpentinization grid has " << var.vtopo->cell_size()
                      << " m cells, but mesh.resolution * mesh.vt_size is " << cell_size
                      << " m; the grid is rebuilt at the next phase change or remesh.\n";
            delete var.vtopo;
            var.vtopo = nullptr;
        }
    }

    // Initializing field variables
    {
        bin_save.read_array(*var.vel, "velocity");
        bin_save.read_array(*var.temperature, "temperature");
        bin_save.read_array(*var.strain, "strain");
        bin_save.read_array(*var.stress, "stress");
        bin_save.read_array(*var.plstrain, "plastic strain");
        // the R2 baseline, so the first remesh after a restart refines as in a continuous run
        bin_chkpt.read_array(*var.plstrain_remesh, "plstrain_remesh");
        bin_save.read_array(*var.radiogenic_source, "radiogenic source");
        bin_save.read_array(*var.ppressure, "pore pressure");
        // previous-step volume for volumetric strain rate.
        bin_chkpt.read_array(*var.volume_old, "volume_old");

        bin_chkpt.read_array(*var.surfinfo.edvacc_surf, "dv surface acc");
        // for surface marker correction after surface processes
        bin_chkpt.read_array(*var.surfinfo.dhacc, "dhacc");
#ifdef USEMMG
        bin_chkpt.read_array(*var.init_elem_size_n, "init_elem_size_n");
        bin_chkpt.read_array(*var.refined_elem_size_n, "refined_elem_size_n");
        {   // the last surface-process dh, stored per node (Output::write_checkpoint)
            double_vec dh_n(var.nnode);
            bin_chkpt.read_array(dh_n, "surface dh");
            for (int i = 0; i < var.surfinfo.ntop; ++i)
                (*var.surfinfo.dh)[i] = dh_n[(*var.surfinfo.top_nodes)[i]];
        }
#endif
        if (param.mat.is_plane_strain)
            bin_chkpt.read_array(*var.stressyy, "stressyy");
    }

    // Restored side-wall columns and wall shifts (remeshing_option 13, SideProfile); no record = no
    // column and no shift at the checkpoint. remesh()'s side-profile lifecycle handles a side that
    // changed at the restart. A restart without option 13 drops them: nothing would erase them,
    // and the marker replenishment would keep refilling from them.
    if (param.mesh.remeshing_option == 13 && bin_chkpt.has_array("side_prof_sizes")) {
        int_vec psz(3 * NSIDEWALL);
        bin_chkpt.read_array(psz, "side_prof_sizes");
        char aname[64];
        double_vec shifts(NSIDEWALL, 0.), rises(NSIDEWALL, 0.);   // wall-restore shifts, top rises
        bin_chkpt.read_array(shifts, "side_prof_wall_shift");
        bin_chkpt.read_array(rises, "side_prof_top_rise");
        for (int t = 0; t < NSIDEWALL; ++t) {
            SideProfile &prof = var.side_profile[t];
            prof.clear_column();
            prof.wall_shift = shifts[t];
            prof.top_rise = rises[t];
            const int nn = psz[t], ne = psz[NSIDEWALL + t], ns = psz[2 * NSIDEWALL + t];
            if (nn <= 0) continue;
            const char *sfx_t = SIDEWALL_NAME[t];
            auto rd = [&](double_vec &a, const char *name, int size) {
                a.resize(size);
                std::snprintf(aname, 64, "side_prof_%s.%s", name, sfx_t);
                bin_chkpt.read_array(a, aname);
            };
            rd(prof.node_reldepth, "node_reldepth", nn);
            rd(prof.node_coord0_z, "node_coord0_z", nn);
            rd(prof.node_temperature, "node_temperature", nn);
#ifdef USEMMG
            // written only under USEMMG; without MMG it must stay empty, which restore_side_fields
            // tests before it indexes var.init_elem_size_n
            rd(prof.node_init_elem_size, "node_init_elem_size", nn);
#endif
#ifdef THREED
            rd(prof.node_strike, "node_strike", nn);
            prof.facet_node.resize(ns * NODES_PER_FACET);
            std::snprintf(aname, 64, "side_prof_facet_node.%s", sfx_t);
            bin_chkpt.read_array(prof.facet_node, aname);
            prof.facet_elem.resize(ns);
            std::snprintf(aname, 64, "side_prof_facet_elem.%s", sfx_t);
            bin_chkpt.read_array(prof.facet_elem, aname);
#else
            prof.node_support_idx.resize(nn + 1);
            std::snprintf(aname, 64, "side_prof_node_support_idx.%s", sfx_t);
            bin_chkpt.read_array(prof.node_support_idx, aname);
            prof.node_support_arr.resize(ns);
            std::snprintf(aname, 64, "side_prof_node_support_arr.%s", sfx_t);
            bin_chkpt.read_array(prof.node_support_arr, aname);
            rd(prof.elem_reldepth, "elem_reldepth", ne);
#endif
            prof.elem_mattype.resize(ne);
            std::snprintf(aname, 64, "side_prof_elem_mattype.%s", sfx_t);
            bin_chkpt.read_array(prof.elem_mattype, aname);
            rd(prof.elem_radiogenic_source, "elem_radiogenic_source", ne);
            side_profile_index_elems(prof);
        }
        std::cout << "  Restored incoming side profiles from checkpoint:";
        for (int t = 0; t < NSIDEWALL; ++t)
            if (psz[t] > 0 || shifts[t] != 0.)
                std::cout << " " << SIDEWALL_NAME[t] << " (" << psz[t] << " nodes, "
                          << psz[NSIDEWALL + t] << " elements, wall shift " << shifts[t]
                          << " m, top rise " << rises[t] << " m)";
        std::cout << ".\n";
    }

    // the following fields are not required for restarting, yet
    {
        // for shear heating
        bin_save.read_array(*var.strain_rate, "strain-rate");
        // for tidal heating
        bin_save.read_array(*var.viscosity, "viscosity");
        bin_save.read_array(*var.force, "force");
        // delta_plstrain carries per-element state across steps:
        bin_save.read_array(*var.delta_plstrain, "plastic strain-rate");
    }

    compute_volume(*var.coord, *var.connectivity, *var.volume);
    // require max_global_vel_mag, var.volume, and var.temperature to be loaded before
    compute_mass(param, var, var.max_vbc_val, *var.volume_n, *var.mass, *var.tmass, *var.hmass, *var.ymass, *var.tmp_result);

    create_boundary_normals(var, *var.bnormals, var.edge_vec, var.edge_slot);

    apply_vbcs(param, var, *var.vel);

    if (param.ic.is_restarting_weakzone) {
        std::cout << "  Creating new weakzone...\n";
        const double_vec old_plstrain(*var.plstrain);
        initial_weak_zone(param, var, *var.plstrain);
        // Shift the R2 baseline by the re-seed only; growth since the last remesh still counts.
        for (int e = 0; e < var.nelem; ++e)
            (*var.plstrain_remesh)[e] += (*var.plstrain)[e] - old_plstrain[e];
    }

    // Restore the RSF state variable when available, then rebuild the derived
    // dynamic friction coefficient from the restored state and current fields.
    if (param.mat.rheol_type & MatProps::rh_rsf) {
        bool restored_state_variable = false;
        if (bin_chkpt.has_array("friction state variable")) {
            bin_chkpt.read_array(*var.state_variable, "friction state variable");
            restored_state_variable = true;
        } else if (bin_save.has_array("friction state variable")) {
            bin_save.read_array(*var.state_variable, "friction state variable");
            restored_state_variable = true;
        } else {
            initial_state_variable(param, var, *var.state_variable);
        }

        if (!restored_state_variable) {
            std::cout << "  RSF restart fallback applied for missing state variable dataset.\n";
        }
        if (param.control.rsf_slip_rate_projection_option ==
            rsf_slip_rate_projection_total_strain_rate) {
            update_strain_rate(var, *var.strain_rate);
            #pragma acc wait
        }
        refresh_rsf_friction(param, var, *var.dyn_fric_coeff, *var.state_variable);
    }

    // For some reason, the following is added by Denis
    // However, it is not clear why this is needed.
    if (param.control.use_global_velocity_scaling) {
        var.dt = compute_dt(param, var);
        compute_mass(param, var, var.max_vbc_val, *var.volume_n, *var.mass, *var.tmass, *var.hmass, *var.ymass, *var.tmp_result);
    }

    report_mesh_info(var, "restart");

}


void end(Variables& var) {
    delete var.output;

    for (size_t i=0; i<var.markersets.size(); i++)
        delete var.markersets[i];
}


void update_mesh(const Param& param, Variables& var)
{
#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif

    update_coordinate(var, *var.coord);

    if(!param.control.PT_jump)
    {
        // option 13: each profiled wall's top rise (SideProfile::top_rise)
        double top_before[NSIDEWALL] = {};
        const bool track_tops = param.mesh.remeshing_option == 13 &&
                                param.control.surface_process_option != 0;
        if (track_tops) {
            for (int t = 0; t < NSIDEWALL; ++t)
                if (!var.side_profile[t].empty())
                    top_before[t] = side_wall_top_z(var, t);
        }
        surface_processes(param, var, *var.coord, *var.stress, *var.strain, *var.strain_rate, \
                      *var.plstrain, *var.volume, *var.volume_n, \
                      var.surfinfo, var.markersets, *var.elemmarkers, *var.markers_in_elem);
        if (track_tops) {
            for (int t = 0; t < NSIDEWALL; ++t) {
                SideProfile &prof = var.side_profile[t];
                if (!prof.empty())
                    prof.top_rise += side_wall_top_z(var, t) - top_before[t];
            }
        }
    }

#ifdef NPROF_DETAIL
    nvtxRangePush("swap vectors");
#endif
    {
        double_vec *tmp = var.volume;
        var.volume = var.volume_old;
        var.volume_old = tmp;
    }
//    var.volume->swap(*var.volume_old);
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif

    compute_volume(var, *var.volume);

    // surface_processes() above can move elemmarkers via correct_surface_marker().
    var.mat->refresh_elem_cache();

    if (param.control.use_global_velocity_scaling) {
        var.dt = compute_dt(param, var);
    }

    compute_mass(param, var, var.max_vbc_val, *var.volume_n, *var.mass, *var.tmass, *var.hmass, *var.ymass, *var.tmp_result);


    #pragma acc wait

#ifdef NPROF
    nvtxRangePop();
#endif
}


void isostasy_adjustment(const Param &param, Variables &var)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    std::cout << "Adjusting isostasy for " << param.ic.isostasy_adjustment_time_in_yr << " yrs...\n";

    var.dt = compute_dt(param, var);
    int iso_steps = param.ic.isostasy_adjustment_time_in_yr*YEAR2SEC / var.dt;

    for (int n=0; n<iso_steps; n++) {
        update_strain_rate(var, *var.strain_rate);
        compute_dvoldt(var, *var.ntmp, *var.etmp);
        compute_edvoldt(var, *var.ntmp, *var.edvoldt);
        update_stress(param, var, *var.stress, *var.stressyy, *var.dpressure,
            *var.viscosity, *var.strain, *var.plstrain, *var.delta_plstrain,
            *var.strain_rate,
            *var.ppressure, *var.dppressure, *var.vel,
            *var.dyn_fric_coeff, *var.state_variable);

        update_force(param, var, *var.force, *var.force_residual, *var.tmp_result);
        update_velocity(var, *var.vel);

        #pragma acc wait

        // do not apply vbc to allow free boundary

        // displacment is vertical only
        #pragma omp parallel for default(none)          \
            shared(var, param)
        for (int i=0; i<var.nnode; ++i) {
            for (int j=0; j<NDIMS-1; ++j) {
                (*var.vel)[i][j] = 0;
            }
            if (param.bc.has_winkler_foundation == false &&
                (*var.bcflag)[i] & BOUNDZ0) {
                // holding bottom surface fixed
                (*var.vel)[i][NDIMS-1] = 0;
            }
        }

        update_mesh(param, var);

    }
    std::cout << "Adjusted isostasy for " << iso_steps << " steps.\n";
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

void initial_body_force_adjustment(const Param &param, Variables &var)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    
    double residual_old = std::numeric_limits<double>::max();
    double relative_change = 1.0;
    
    // pseudo transient (PT) loop
    var.l2_residual = calculate_residual_force(var, *var.force_residual);
    residual_old = var.l2_residual;
    if (param.control.has_PT)
    {   
        // var.dt = compute_dt_PT(param, var);
        param.control.PT_jump = true;
        for (int pt_step = 0; pt_step < param.control.PT_max_iter; ++pt_step) 
        {
            apply_vbcs(param, var, *var.vel);
            if (param.control.has_moving_mesh)
                update_mesh(param, var);
            update_strain_rate(var, *var.strain_rate);
            compute_dvoldt(var, *var.ntmp, *var.etmp);
            compute_edvoldt(var, *var.ntmp, *var.edvoldt);
            update_stress(param, var, *var.stress, *var.stressyy, *var.dpressure,
                *var.viscosity, *var.strain, *var.plstrain, *var.delta_plstrain,
                *var.strain_rate,
                *var.ppressure, *var.dppressure, *var.vel,
                *var.dyn_fric_coeff, *var.state_variable);
            update_force(param, var, *var.force, *var.force_residual, *var.tmp_result);
            // update_velocity_PT(var, *var.vel);
            update_velocity(var, *var.vel);
            var.l2_residual = calculate_residual_force(var, *var.force_residual);
            double relative_change = std::fabs((var.l2_residual - residual_old) / residual_old);
            if (relative_change < param.control.PT_relative_tolerance) {
            break;  // Exit the loop if relative change is small enough
            }
            residual_old = var.l2_residual;
        }
        param.control.PT_jump = false;
    }

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

int main(int argc, const char* argv[])
{
    // Only main() can tell its own default from the operator's (omp_wait_policy_src).
    bool wait_policy_from_des = false;
#if defined(__APPLE__) && defined(_OPENMP)
    // macOS/LLVM libomp sets blocktime=0 on Apple Silicon (hybrid CPU detection),
    // causing threads to immediately yield between parallel regions. This hurts
    // throughput in tight time loops with many short parallel regions.
    // OMP_WAIT_POLICY=active restores spin-wait behavior. No-op if already set.
    if (!env_sets_wait_policy()) {
        setenv("OMP_WAIT_POLICY", "active", 1);   // 1: overwrite an empty export
        wait_policy_from_des = true;
        std::cout << "[OpenMP] macOS: OMP_WAIT_POLICY=active set to avoid libomp "
                     "zero-blocktime regression on Apple Silicon.\n"
                     "         Override: export OMP_WAIT_POLICY=passive | KMP_BLOCKTIME=<ms>\n";
    }
#endif

    //
    // read command line
    //
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " config_file\n";
        std::cout << "       " << argv[0] << " -h or --help\n";
        return -1;
    }

    Param param;
    get_input_parameters(argv[1], param);

    const BuildInfo build = probe_build_info();
    const CpuInfo cpu = probe_cpu_info(wait_policy_from_des);
    // Selects the offload device, so it must precede any compute.
    const DeviceInfo dev = init_offload_device();
    const Manifest manifest = compose_manifest(param, build, cpu, dev);
    if (param.sim.has_runtime_info_display)
        report_build_and_runtime_info(manifest);
    // An -acc=gpu binary cannot compute on the host (its managed allocations need a
    // device): stop here, not in the first Array2D with a bare exit 1.
    if (dev.acc_build && !dev.using_gpu) {
        std::cerr << "Error: this is an OpenACC (GPU) build and no usable device was "
                     "found; check CUDA_VISIBLE_DEVICES and ACC_DEVICE_TYPE, or run a "
                     "CPU build\n";
        write_manifest(param, manifest, true);
        die(EXIT_RUNTIME_RESOURCE);
    }

    //
    // run simulation
    //
    static Variables var; // declared as static to silence valgrind's memory leak detection
    init_var(param, var);

    var.output = new Output(param, build, cpu, dev, manifest, var.func_time.start_time,
                  (param.sim.is_restarting) ? param.sim.restarting_from_frame : 0);

    if (! param.sim.is_restarting) {
        init(param, var);
        var.info_display_next_step = param.sim.info_display_step_interval;

        if (param.ic.isostasy_adjustment_time_in_yr > 0) {
            // output.write_exact(var);
            isostasy_adjustment(param, var);
        }

        var.dt = compute_dt(param, var);

        // Option 13: record the inflow walls' columns (detect_side_profile) after mesh and markers
        // exist and before the initial checkpoint. A restart reads the checkpointed ones instead.
        if (param.mesh.remeshing_option == 13) {
            const double t_now = var.time / YEAR2SEC;
            uint bits = 0;
            for (int t = 0; t < NSIDEWALL; ++t)
                if (side_wall_vbc_is_inflow(param.bc, t, t_now)) bits |= SIDEWALL_FLAG[t];
            if (bits)
                detect_side_profile(param, var, bits);
        }

        // Build the serpentinization grid (mat.phase_change_option 3 with
        // mat.has_serpentinization) up front rather than at the first
        // phase_changes step.  Built after isostasy adjustment so the
        // bounding box sees the settled mesh.
        build_serpentinization_grid(param, var);

        if (param.sim.has_initial_checkpoint)
            var.output->write_checkpoint(param, var);
    }
    else {
        // The serpentinization grid is not built here: restart() restores it
        // or leaves var.vtopo null (see the gate there).
        restart(param, var);
    }

    // var.dt_PT = var.dt;

#ifdef HAS_GOSPL_CPP_INTERFACE
    // Initialize GoSPL driver if surface process option is 11
    // This runs for both fresh starts and restarts
    if (param.control.surface_process_option == 11) {
        std::cout << "Initializing GoSPL driver for surface process option 11..." << std::endl;
        var.gospl_driver = new GoSPLDriver();
        
        // Step 1: Initialize Python/gospl_extensions (needed for generate_mesh)
        if (!var.gospl_driver->init_python()) {
            std::cerr << "Failed to initialize Python for GoSPL" << std::endl;
            delete var.gospl_driver;
            var.gospl_driver = nullptr;
            return -1;
        }
        
        // Step 2: Generate GoSPL mesh from DynEarthSol surface nodes (or reuse existing on restart)
        std::string mesh_file = "gospl_mesh.npz";
        bool reuse_mesh = false;
        
        if (param.sim.is_restarting) {
            // Only check for existing mesh when restarting
            std::ifstream mesh_check(mesh_file);
            reuse_mesh = mesh_check.good();
            mesh_check.close();
        }
        
        if (reuse_mesh) {
            std::cout << "Restarting: Found existing GoSPL mesh file: " << mesh_file << ". Reusing it." << std::endl;
        } else {
            std::cout << "Generating new GoSPL mesh..." << std::endl;
            std::vector<double> x_coords;
            std::vector<double> y_coords;
            
            for (int i = 0; i < var.nnode; ++i) {
                if ((*var.bcflag)[i] & BOUNDZ1) {  // Top surface nodes
                    x_coords.push_back((*var.coord)[i][0]);
#ifdef THREED
                    y_coords.push_back((*var.coord)[i][1]);
#else
                    y_coords.push_back(0.0);  // 2D case
#endif
                }
            }
            
            if (!x_coords.empty()) {
                var.gospl_driver->generate_mesh(x_coords, y_coords, mesh_file,
                                                param.control.gospl_mesh_resolution,
                                                param.control.gospl_mesh_perturbation,
                                                param.control.gospl_mesh_padding);
            }
        }
        
        // Step 3: Initialize GoSPL model (loads the mesh)
        if (var.gospl_driver->initialize(param.control.surface_process_gospl_config_file)) {
            std::cout << "GoSPL driver successfully initialized with config: " 
                      << param.control.surface_process_gospl_config_file << std::endl;
            // Suppress verbose output from GoSPL
            var.gospl_driver->set_verbose(false);
            // Set coupling frequency from config (default: 1 = every step)
            var.gospl_driver->coupling_by_time   = (param.control.gospl_coupling_mode == "time");
            var.gospl_driver->coupling_frequency = param.control.gospl_coupling_frequency;
            var.gospl_driver->coupling_interval_in_yr  = param.control.gospl_coupling_interval_in_yr;
            var.gospl_driver->velocity_coupling  = param.control.gospl_velocity_coupling;
            if (var.gospl_driver->coupling_by_time)
                std::cout << "GoSPL coupling mode: every " << var.gospl_driver->coupling_interval_in_yr << " yr"
                          << (var.gospl_driver->velocity_coupling ? " [velocity coupling enabled]" : "") << std::endl;
            else
                std::cout << "GoSPL coupling mode: every " << var.gospl_driver->coupling_frequency << " step(s)"
                          << (var.gospl_driver->velocity_coupling ? " [velocity coupling enabled]" : "") << std::endl;
        } else {
            std::cerr << "Failed to initialize GoSPL driver with config: " 
                      << param.control.surface_process_gospl_config_file << std::endl;
            delete var.gospl_driver;
            var.gospl_driver = nullptr;
            return -1;
        }
    }
#endif

    int64_t init_time = get_nanoseconds() - var.func_time.start_time;

    var.output->write_exact(var);
    monitor_initialize(param, var);

    // int rheol_type_old = param.mat.rheol_type;

    const double starting_time = var.reference_frame_time;
    var.reference_frame_time = starting_time + param.sim.output_time_interval_in_yr * YEAR2SEC;
    const double starting_step = var.steps;
    int next_regular_frame = 1;  // excluding frames due to output_during_remeshing

    EarthquakeState earthquake;
    init_earthquake_state(param, earthquake);

    double residual_old = std::numeric_limits<double>::max();
    double relative_change = 1.0;
    double dt_copy = 0.0;
    bool hydraulic_diffusion_switch = false;

    if(param.ic.has_body_force_adjustment)
    {
        if(param.control.has_hydraulic_diffusion) {param.control.has_hydraulic_diffusion = false; hydraulic_diffusion_switch = true;}
        // this is similar to isostasy_adjustment(param, var); so maybe should be merged to it later.
        // Only works with PT loop
        initial_body_force_adjustment(param, var); 
        if(hydraulic_diffusion_switch) {param.control.has_hydraulic_diffusion = true;}
        param.ic.has_body_force_adjustment = false;
    }

    std::cout << "Starting simulation...\n";
    std::cout << "  Showing model progress every "
              << param.sim.info_display_step_interval
              << " steps.\n";

    do {
#ifdef NPROF_DETAIL
        nvtxRangePush("dynearthsol");
        nvtxRangePush("kernel");
#endif
        var.steps ++;
        var.time += var.dt;
        // Pick up what the previous step's phase changes and remeshing moved.
        var.mat->refresh_elem_cache();
        // dt_copy = 0.0; dt_copy += var.dt;
        if (param.control.has_thermal_diffusion)
            update_temperature(param, var, *var.temperature, *var.tmp_result);

        if (param.control.has_hydraulic_diffusion)
            update_old_mean_stress(param, var, *var.stress, *var.old_mean_stress);

        update_strain_rate(var, *var.strain_rate);
        compute_dvoldt(var, *var.ntmp, *var.etmp);
        compute_edvoldt(var, *var.ntmp, *var.edvoldt);
        update_stress(param, var, *var.stress, *var.stressyy, *var.dpressure,
            *var.viscosity, *var.strain, *var.plstrain, *var.delta_plstrain,
            *var.strain_rate,
            *var.ppressure, *var.dppressure, *var.vel,
            *var.dyn_fric_coeff, *var.state_variable);

        // Nodal Mixed Discretization For Stress
        if (param.control.is_using_mixed_stress)
            NMD_stress(var, *var.stress, *var.ntmp, *var.etmp);
            
        update_force(param, var, *var.force, *var.force_residual, *var.tmp_result);
        update_velocity(var, *var.vel);

        // pseudo transient (PT) loop
        var.l2_residual = calculate_residual_force(var, *var.force_residual);
        residual_old = var.l2_residual;
        if (param.control.has_PT)
        {   
            // var.dt = compute_dt_PT(param, var);
            if (param.control.has_hydraulic_diffusion) {
                param.control.has_hydraulic_diffusion = false;
                hydraulic_diffusion_switch = true;
            }

            param.control.PT_jump = true;
            for (int pt_step = 0; pt_step < param.control.PT_max_iter; ++pt_step) 
            {
                apply_vbcs(param, var, *var.vel);
                if (param.control.has_moving_mesh)
                    update_mesh(param, var);
                update_strain_rate(var, *var.strain_rate);
                compute_dvoldt(var, *var.ntmp, *var.etmp);
                compute_edvoldt(var, *var.ntmp, *var.edvoldt);
                update_stress(param, var, *var.stress, *var.stressyy, *var.dpressure,
                    *var.viscosity, *var.strain, *var.plstrain, *var.delta_plstrain,
                    *var.strain_rate,
                    *var.ppressure, *var.dppressure, *var.vel,
                    *var.dyn_fric_coeff, *var.state_variable);
                update_force(param, var, *var.force, *var.force_residual, *var.tmp_result);
                // update_velocity_PT(var, *var.vel);
                update_velocity(var, *var.vel);
                var.l2_residual = calculate_residual_force(var, *var.force_residual);
                double relative_change = std::fabs((var.l2_residual - residual_old) / residual_old);
                if (relative_change < param.control.PT_relative_tolerance) {
                    // std::cout << "tolerance reached " << pt_step << std::endl;
                break;  // Exit the loop if relative change is small enough
                }
                residual_old = var.l2_residual;
                // var.dt = std::min({var.dt*1.01, dt_copy});

                if (pt_step % param.mesh.quality_check_step_interval == 0) {
                    if (param.control.has_moving_mesh)
                    {
                        int quality_is_bad, bad_quality_index;
                        double min_quality;
                        quality_is_bad = bad_mesh_quality(param, var, bad_quality_index, min_quality);
                        if (quality_is_bad) {

                            if (param.sim.has_output_during_remeshing) {
                                var.output->write_exact(var);
                            }

                            monitor_before_remesh(param, var);
                            remesh(param, var, quality_is_bad);
                            monitor_remesh_update(param, var);

                            if (param.sim.has_output_during_remeshing) {
                                var.output->write_exact(var);
                            }
                        }
                    }
                }
            }
            if(hydraulic_diffusion_switch) {param.control.has_hydraulic_diffusion = true;}
            // var.dt = dt_copy;
            param.control.PT_jump = false;

        }


        if(param.control.has_hydraulic_diffusion)
            update_pore_pressure(param, var, *var.ppressure, *var.dppressure, *var.ntmp, *var.tmp_result, *var.stress, *var.old_mean_stress);

        // The functions of a slow step are expensive and their inputs change
        // slowly, so they run every slow_updates_interval steps.
        const int slow_updates_interval = 10;
        const bool slow_step = var.steps % slow_updates_interval == 0;
        // phase_changes() converts nothing with one material or option 0.
        const bool phase_change_step = slow_step && param.mat.nmat > 1 &&
                                       param.mat.phase_change_option != 0;
        if (phase_change_step) {
            // Before update_mesh(): its compute_mass() then sees the converted
            // composition, for the next step and for a restart from this step,
            // which rebuilds the masses from it.
            phase_changes(param, var);
            var.mat->refresh_elem_cache();
        }

        // Objective rotation still belongs to the current physical step, even
        // if the update below selects var.dt for the next step.
        const double rotation_dt = var.dt;
        apply_vbcs(param, var, *var.vel);
        if (param.control.has_moving_mesh)
            update_mesh(param, var);
        else if (param.control.rsf_dtheta_max > 0.0) {
            // Select the next step from the accepted velocity. The GVS mass
            // must follow the same current rate used by the state bound.
            var.dt = compute_dt(param, var);
            compute_mass(param, var, var.max_vbc_val, *var.volume_n,
                         *var.mass, *var.tmass, *var.hmass, *var.ymass,
                         *var.tmp_result);
        }
        else if (phase_change_step) {
            // A fixed mesh has no update_mesh() to rebuild them.
            compute_mass(param, var, var.max_vbc_val, *var.volume_n,
                         *var.mass, *var.tmass, *var.hmass, *var.ymass,
                         *var.tmp_result);
        }

        // Rebuild the serpentinization grid at each quality check. Markers of
        // markersets[0] are added or removed only by the surface processes in
        // update_mesh() above, under this same step test, and by a remesh, which
        // builds the grid at its end. So the grid always indexes every marker,
        // and it does not go stale in a run that seldom remeshes. A checkpoint
        // holds the grid the run reads.
        if (var.steps % param.mesh.quality_check_step_interval == 0)
            build_serpentinization_grid(param, var);

        // elastic stress/strain are objective (frame-indifferent)
        if (var.mat->rheol_type & MatProps::rh_elastic)
            rotate_stress(var, *var.stress, *var.strain, rotation_dt);

        monitor_write_if_due(param, var);

        #pragma acc wait
        if (slow_step) {
            if (param.control.has_hydration_processes)
                advect_hydrous_markers(param, var, 10*var.dt,
                                       *var.markersets[var.hydrous_marker_index],
                                       *var.hydrous_elemmarkers);
            var.dt = compute_dt(param, var);
        }

        #pragma acc wait

#ifdef NPROF_DETAIL
        nvtxRangePop();
#endif

        if (param.sim.is_outputting_averaged_fields)
            var.output->average_fields(var);

        update_earthquake_tracking(param, var, earthquake);

        if (param.control.use_global_velocity_scaling) {
            handle_ats_output(param, var, *var.output, earthquake,
                              starting_time, starting_step, next_regular_frame);
        } else {
            if (( (param.sim.output_step_interval != std::numeric_limits<int>::max() &&
               (var.steps - starting_step) == next_regular_frame * param.sim.output_step_interval)
              ||
              (param.sim.output_time_interval_in_yr != std::numeric_limits<double>::max() &&
               (var.time - starting_time) > next_regular_frame * param.sim.output_time_interval_in_yr * YEAR2SEC)
             )
            // time or step output requirements are met
            &&
            ((! param.sim.is_outputting_averaged_fields) ||
                (param.sim.is_outputting_averaged_fields &&
                    (var.steps % param.mesh.quality_check_step_interval == 0)))
            // When is_outputting_averaged_fields is turned on, the output cannot be
            // done at arbitrary time steps.
            ) {
                if (next_regular_frame % param.sim.checkpoint_frame_interval == 0)
                    var.output->write_checkpoint(param, var);

                var.output->write(var);

                next_regular_frame ++;
                var.reference_frame_time = starting_time +
                    next_regular_frame * param.sim.output_time_interval_in_yr * YEAR2SEC;

            }
        }

        if (var.steps % param.mesh.quality_check_step_interval == 0) {
            double min_quality = 1.0;
            if (param.control.has_moving_mesh)
            {
                int quality_is_bad, bad_quality_index;
                quality_is_bad = bad_mesh_quality(param, var, bad_quality_index, min_quality);
                if (quality_is_bad) {

                    if (param.sim.has_output_during_remeshing) {
                        var.output->write_exact(var);
                    }

                    monitor_before_remesh(param, var);
                    remesh(param, var, quality_is_bad);
                    monitor_remesh_update(param, var);

                    if (param.sim.has_output_during_remeshing) {
                        var.output->write_exact(var);
                    }
                }
            }

            if (var.steps >= var.info_display_next_step) {
                int64_t now_ns = get_nanoseconds();
                std::cout << "              Step = " << var.steps
                    << ", time = " << std::scientific << std::setprecision(5)
                    << var.time / YEAR2SEC << " yr"
                    << ", vmax = " << var.max_global_vel_mag << " m/s";
                if (param.control.has_moving_mesh && min_quality < 1.0)
                    std::cout << ", min_q = " << std::fixed << std::setprecision(4) << min_quality;
                std::cout << ", wt = ";
                print_time_ns(now_ns - var.func_time.start_time);
                std::cout << "\n";
#ifdef HAS_GOSPL_CPP_INTERFACE
                if (var.gospl_driver != nullptr)
                    var.gospl_driver->print_coupling_summary();
#endif


                var.info_display_next_step = var.steps + param.sim.info_display_step_interval;
            }
        }
#ifdef NPROF_DETAIL
        nvtxRangePop();
#endif

    } while (var.steps < param.sim.max_steps &&
             (var.time <= param.sim.max_time_in_yr * YEAR2SEC ||
              (param.sim.is_outputting_averaged_fields &&
               var.steps % param.mesh.quality_check_step_interval != 0)));

#ifdef HAS_GOSPL_CPP_INTERFACE
    // Report the coupling events since the last progress line, if any.
    if (var.gospl_driver != nullptr && var.gospl_driver->summary_events > 0)
        var.gospl_driver->print_coupling_summary();
#endif

    monitor_finalize(var);

    // Before end() deletes the Output that keeps the run's peak.
    const double peak_rss_gib = var.output->update_peak_rss_gib();

    // at end of code, clean up lost memory reported by valgrind
    end(var);

    std::cout << "Ending simulation.\n";
    int64_t duration_ns = get_nanoseconds() - var.func_time.start_time;
    int64_t computation_time = duration_ns - var.func_time.remesh_time - var.func_time.output_time
        - var.func_time.phase_change_time - init_time;
    double computation_time_avg = computation_time * 1.e-9 / var.steps;

    std::cout << "Time summary...\n  Execute : ";
    print_time_ns(duration_ns);
    std::cout << "\n  Initiate: ";
    print_time_ns(init_time);
    std::cout << " (" <<  std::setw(5) << std::fixed << std::setprecision(2) << std::setfill(' ')
        << 100.*init_time/duration_ns << "%)\n";
    std::cout << "  Compute : ";
    print_time_ns(computation_time);
    std::cout << " (" <<  std::setw(5) << std::fixed << std::setprecision(2) << std::setfill(' ')
        << 100.*computation_time/duration_ns << "%)/ " << var.steps
        << " = " << std::setprecision(6) << computation_time_avg << " s/step\n";
    if (var.nremesh > 0) {
        double remesh_time_avg = var.func_time.remesh_time * 1.e-9 / var.nremesh;
        std::cout << "  Remesh  : ";
        print_time_ns(var.func_time.remesh_time);
        std::cout << " (" <<  std::setw(5) << std::fixed << std::setprecision(2) << std::setfill(' ')
            << 100.*var.func_time.remesh_time/duration_ns << "%)/ " << var.nremesh
            << " = " << std::setprecision(2) << remesh_time_avg << " s/remesh\n";
    }
    if (var.nphasechange > 0) {
        double phase_change_time_avg = var.func_time.phase_change_time * 1.e-9 / var.nphasechange;
        std::cout << "  PhaseChg: ";
        print_time_ns(var.func_time.phase_change_time);
        std::cout << " (" <<  std::setw(5) << std::fixed << std::setprecision(2) << std::setfill(' ')
            << 100.*var.func_time.phase_change_time/duration_ns << "%)/ " << var.nphasechange
            << " = " << std::setprecision(6) << phase_change_time_avg << " s/call\n";
    }
    if (var.noutput > 0) {
        double output_time_avg = var.func_time.output_time * 1.e-9 / var.noutput;
        std::cout << "  Output  : ";
        print_time_ns(var.func_time.output_time);
        std::cout << " (" <<  std::setw(5) <<  std::fixed << std::setprecision(2) << std::setfill(' ')
            << 100.*var.func_time.output_time/duration_ns << "%)/ " << var.noutput
            << " = " << std::setprecision(2) << output_time_avg << " s/output\n";
#ifdef HAS_GOSPL_CPP_INTERFACE
        // Clean up GoSPL driver if it was created
        if (var.gospl_driver != nullptr) {
            std::cout << "Cleaning up GoSPL driver..." << std::endl;
            delete var.gospl_driver;
            var.gospl_driver = nullptr;
        }
#endif
    }

    write_manifest_end(param, manifest, var, var.steps - static_cast<int>(starting_step),
                       duration_ns, init_time, computation_time, peak_rss_gib);
    return 0;
}
