#include <algorithm>  // For std::is_sorted
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <sstream>

#include "options.hpp"
#include "parameters.hpp"
#include "fields.hpp"
#include "matprops.hpp"
#include "phasechanges.hpp"
#include "utils.hpp"


static void declare_point_source(opt::Registry &cfg, const char *section,
                                 Injection &source)
{
    const auto key = [section](const char *name) {
        return std::string(section) + "." + name;
    };
    cfg.add_options()
        (key("enabled").c_str(), opt::value<bool>(&source.enabled)->default_value(false),
         "Enable scheduled point sources.")
        (key("num_points").c_str(), opt::value<int>(&source.num_points)->default_value(0),
         "Number of source points.")
        (key("points_x").c_str(), opt::list()->default_value("[]"),
         "Point coordinates; in 2D points_y is the legacy alias for points_z.")
        (key("points_y").c_str(), opt::list()->default_value("[]"),
         "Point coordinates; in 2D points_y is the legacy alias for points_z.")
        (key("points_z").c_str(), opt::list()->default_value("[]"),
         "Point coordinates; in 2D points_y is the legacy alias for points_z.")
        (key("points_unit").c_str(), opt::value<std::string>(&source.points_unit)->default_value("m"),
         "Coordinate unit: mm, cm, m or km.")
        (key("rate_model").c_str(), opt::value<std::string>()->default_value("constant_rate"),
         "Source schedule: constant_rate or total_amount.")
        (key("rate").c_str(), opt::list()->default_value("[0]"),
         "Positive rate supplies fluid volume or heat; negative rate extracts.")
        (key("total_amount").c_str(), opt::list()->default_value("[0]"),
         "Total fluid volume or heat over the scheduled interval.")
        (key("start_time_in_yr").c_str(), opt::list()->default_value("[0]"),
         "Source start time in years.")
        (key("end_time_in_yr").c_str(), opt::list()->default_value("[]"),
         "Finite end times in years; alternative to duration.")
        (key("duration_in_yr").c_str(), opt::list()->default_value("[]"),
         "Positive schedule durations in years.")
        ;
}


static void declare_parameters(opt::Registry &cfg,
                               Param &p)
{
    /* To have a new input parameter declared as such in parameters.hpp,
     *
     *     struct SectionType {
     *         type name;
     *     }
     *     struct Param {
     *         SectionType section;
     *     }
     *
     * add this line in this function:
     *
     *     ("section.name", opt::value<type>(&p.section.name), "help string")
     *
     */

    declare_point_source(cfg, "injection", p.injection);
    declare_point_source(cfg, "heat_injection", p.heat_injection);
    cfg.add_options()
        ("injection.temperature", opt::list()->default_value("[]"),
         "Temperature (K) of the injected fluid of each point; requires has_fluid_thermal_advection. "
         "An injecting node gains rho_f*c_f*Q*(T_in - T); extraction leaves at the local temperature. "
         "Unset: the fluid source carries no enthalpy and a heat source must be given separately.")
        ;
    cfg.add_options()
        ("sim.modelname", opt::value<std::string>(&p.sim.modelname),
         "Prefix for the output files")

        ("sim.max_steps", opt::value<int>(&p.sim.max_steps),
         "Max. number of time steps")
        ("sim.max_time_in_yr", opt::value<double>(&p.sim.max_time_in_yr),
         "Max. time (in years)")
        ("sim.output_step_interval", opt::value<int>(&p.sim.output_step_interval),
         "Output step interval")
        ("sim.output_time_interval_in_yr", opt::value<double>(&p.sim.output_time_interval_in_yr),
         "Output time interval (in years)")
        ("sim.info_display_step_interval", opt::value<int>(&p.sim.info_display_step_interval)->default_value(0),
         "Step interval for showing model status on screen.\n"
         " 0: use default (= 100 * mesh.quality_check_step_interval).\n"
         ">0: rounded up to a multiple of mesh.quality_check_step_interval.")
        ("sim.has_runtime_info_display", opt::value<bool>(&p.sim.has_runtime_info_display)->default_value(true),
         "Print the .manifest's sections on screen at run start, one line each,\n"
         "[build] then [runtime]. On by default: it is the start-of-run report every\n"
         "run printed before this option. The .manifest is written regardless.")

        ("sim.checkpoint_frame_interval", opt::value<int>(&p.sim.checkpoint_frame_interval)->default_value(10),
         "How frequent to write checkpoint file (used for restarting simulation)?")
        ("sim.restarting_from_modelname", opt::value<std::string>(&p.sim.restarting_from_modelname),
         "Prefix for the checkpoint and save files for restarting.")
        ("sim.restarting_from_frame", opt::value<int>(&p.sim.restarting_from_frame),
         "Which output frame to read for restarting?")
        ("sim.is_restarting", opt::value<bool>(&p.sim.is_restarting)->default_value(false),
         "Restarting from previous checkpoint file?")

        ("sim.has_initial_checkpoint", opt::value<bool>(&p.sim.has_initial_checkpoint)->default_value(true),
         "Output checkpoint file at 0th step?")
        ("sim.has_marker_output", opt::value<bool>(&p.sim.has_marker_output)->default_value(true),
         "Output marker coordinate and material?")
        ("sim.has_output_during_remeshing", opt::value<bool>(&p.sim.has_output_during_remeshing)->default_value(false),
         "Output immediately before and after remeshing?")
        ("sim.is_outputting_averaged_fields", opt::value<bool>(&p.sim.is_outputting_averaged_fields)->default_value(true),
         "Output time-averaged (smoothed) field variables or not. These fields are: velocity, strain rate, and stress.\n"
         "no: output instantaneous fields. The velocity and strain-rate might oscillate temporally.\n"
         "yes: output field variables averaged over mesh.quality_check_step_interval time steps.\n")
        // For rate-and-state + earthquake output control
        ("sim.earthquake_output_step_interval", opt::value<int>(&p.sim.earthquake_output_step_interval)->default_value(50),
         "Output step interval while earthquake mode is active.")
        ("sim.earthquake_start_factor", opt::value<double>(&p.sim.earthquake_start_factor)->default_value(5.0),
         "Earthquake start threshold factor: vmax > factor * max_vbc_val.")
        ("sim.earthquake_end_factor", opt::value<double>(&p.sim.earthquake_end_factor)->default_value(2.0),
         "Earthquake end threshold factor: vmax < factor * max_vbc_val.")
        ("sim.seismic_moment_calculate_output", opt::value<bool>(&p.sim.seismic_moment_calculate_output)->default_value(false),
         "If true, the seismic moment, sum of 2 G V times the plastic-strain increment over\n"
         "all elements, is accumulated during earthquake mode and written with Mw at its end.\n"
         "In 2D V is the element area and M0 is per unit length out of the plane, without Mw.")
        ("sim.hdf5_compression_level", opt::value<int>(&p.sim.hdf5_compression_level)->default_value(4),
            "HDF5 compression level (0-9). 0: no compression; 9: max. compression. \n"
            "A higher compression level could slow down the output speed.")
        ;

    cfg.add_options()
        ("monitor.enabled", opt::value<bool>(&p.monitor.enabled)->default_value(false),
         "Enable point-based generic monitoring output?")
        ("monitor.step_interval", opt::value<int>(&p.monitor.step_interval)->default_value(10),
         "Monitoring output step interval.")
        ("monitor.num_points", opt::value<int>(&p.monitor.num_points)->default_value(0),
         "Number of monitoring query points.")
        ("monitor.points_x", opt::list()->default_value("[]"),
         "Monitoring x coordinates array '[x0, x1, ...]'.")
        ("monitor.points_y", opt::list()->default_value("[]"),
         "Monitoring y coordinates array '[y0, y1, ...]' (3D only; 2D legacy alias for monitor.points_z).")
        ("monitor.points_z", opt::list()->default_value("[]"),
         "Monitoring z coordinates array '[z0, z1, ...]' (2D second coordinate or 3D third coordinate).")
        ("monitor.points_unit", opt::value<std::string>(&p.monitor.points_unit)->default_value("m"),
         "Unit of monitor.points_x/z in 2D, or monitor.points_x/y/z in 3D: mm, cm, m, km.")
        ("monitor.point_type", opt::value<std::string>()->default_value("particle"),
         "Each point follows the material through deformation, remeshing and restart. particle: values "
         "interpolated at the point; grid: values of the node and element nearest it.")
        ("monitor.remesh_output", opt::value<std::string>()->default_value("before"),
         "Rows at a remesh, marked in the last CSV column, row. before: the state before each remesh "
         "(row 1) unless that step already has a regular row (row 0); both: also the remapped state "
         "after it (row 2).")
        ("monitor.output_prefix", opt::value<std::string>(&p.monitor.output_prefix)->default_value("monitor"),
         "Prefix of monitoring output files.")
        ("monitor.write_header", opt::value<bool>(&p.monitor.write_header)->default_value(true),
         "Write CSV header row for monitoring outputs?")
        ("monitor.output_coord", opt::value<bool>(&p.monitor.output_coord)->default_value(true),
         "Output node coordinate columns.")
        ("monitor.output_velocity", opt::value<bool>(&p.monitor.output_velocity)->default_value(true),
         "Output node velocity columns.")
        ("monitor.output_force", opt::value<bool>(&p.monitor.output_force)->default_value(false),
         "Output node force columns.")
        ("monitor.output_temperature", opt::value<bool>(&p.monitor.output_temperature)->default_value(false),
         "Output node temperature.")
        ("monitor.output_pore_pressure", opt::value<bool>(&p.monitor.output_pore_pressure)->default_value(false),
         "Output node pore pressure.")
        ("monitor.output_bcflag", opt::value<bool>(&p.monitor.output_bcflag)->default_value(false),
         "Output node boundary-flag value.")
        ("monitor.output_stress", opt::value<bool>(&p.monitor.output_stress)->default_value(false),
         "Output element stress components.")
        ("monitor.output_strain", opt::value<bool>(&p.monitor.output_strain)->default_value(false),
         "Output element strain components.")
        ("monitor.output_strain_rate", opt::value<bool>(&p.monitor.output_strain_rate)->default_value(false),
         "Output element strain-rate components.")
        ("monitor.output_plastic_strain", opt::value<bool>(&p.monitor.output_plastic_strain)->default_value(false),
         "Output element plastic strain.")
        ("monitor.output_plastic_strain_rate", opt::value<bool>(&p.monitor.output_plastic_strain_rate)->default_value(false),
         "Output element plastic strain-rate.")
        ("monitor.output_radiogenic_source", opt::value<bool>(&p.monitor.output_radiogenic_source)->default_value(false),
         "Output element radiogenic source.")
        ("monitor.output_density", opt::value<bool>(&p.monitor.output_density)->default_value(false),
         "Output element density.")
        ("monitor.output_mesh_quality", opt::value<bool>(&p.monitor.output_mesh_quality)->default_value(false),
         "Output element mesh quality.")
        ("monitor.output_viscosity", opt::value<bool>(&p.monitor.output_viscosity)->default_value(false),
         "Output element viscosity.")
        ("monitor.output_material", opt::value<bool>(&p.monitor.output_material)->default_value(false),
         "Output element material index.")
        ("monitor.output_dynamic_friction", opt::value<bool>(&p.monitor.output_dynamic_friction)->default_value(false),
         "Output dynamic friction value.")
        ("monitor.output_state_variable", opt::value<bool>(&p.monitor.output_state_variable)->default_value(false),
         "Output state variable value.")
        ;

    cfg.add_options()
        ("mesh.meshing_option", opt::value<int>(&p.mesh.meshing_option)->default_value(1),
         "How to create the new mesh?\n"
         "1: rectangular box with roughly uniform resolution\n"
         "2: rectangular box with rectangular zone of refined resolution\n"
         "90: read the bounding polygon from a .poly file\n"
         "91: same as 90, but the max. element size is normalized with resolution.\n"
         "95: read the mesh from from a .exo file\n"
         )
        ("mesh.meshing_elem_shape", opt::value<int>(&p.mesh.meshing_elem_shape)->default_value(0),
         "Mesh layout. 0: unstructured triangles (2D) or tetrahedra (3D); 1: regular rectangular "
         "cells split into triangles/tetrahedra; 2: equilateral triangles (2D) or a regular "
         "tetrahedral layout (3D). All solver elements are simplices.")
        ("mesh.meshing_verbosity", opt::value<int>(&p.mesh.meshing_verbosity)->default_value(-1),
         "Screen output of meshing and remeshing. -1: no stage log (default). 0: the stage log "
         "([metric], [collapse], [freeze], [unfreeze], [require], [corner], [mmg], "
         "[restore], [GPU mem] lines: what a remesh re-sized, collapsed, froze, freed and restored) "
         "and the mesher's standard report. 1-3: Triangle/TetGen -V levels on top. Warnings print "
         "at every level.")
        ("mesh.meshing_sediment", opt::value<bool>(&p.mesh.meshing_sediment)->default_value(false),
         "Refine the elements with sediment marker during remeshing.")
        ("mesh.tetgen_optlevel", opt::value<int>(&p.mesh.tetgen_optlevel)->default_value(3),
         "Optimization level for tetgen. 0: no optimization; 1: multiple edge filps; 2: 1 & free vertex deletion; 3: 2 & new vertex insertion. High optimization level could slow down the speed of mesh generation. For 3D only.")

        ("mesh.xlength", opt::value<double>(&p.mesh.xlength)->required(),
         "Length of x (in meters)")
        ("mesh.ylength", opt::value<double>(&p.mesh.ylength)->required(),
         "Length of y (in meters), for 3D only")
        ("mesh.zlength", opt::value<double>(&p.mesh.zlength)->required(),
         "Length of z (in meters)")

        ("mesh.resolution", opt::value<double>(&p.mesh.resolution)->required(),
         "Spatial resolution (in meters)")
        ("mesh.smallest_size", opt::value<double>(&p.mesh.smallest_size)->default_value(0.01),
         "The size of smallest element relative to an element with typical resolution. It sets the "
         "tiny-element remesh trigger smallest_vol = smallest_size * sizefactor * resolution^NDIMS "
         "(sizefactor 0.433 in 2D, 0.118 in 3D: a unit-edge equilateral element's measure) and two "
         "lengths, hmin = smallest_size^(1/NDIMS) * resolution and min_dist = (smallest_size * "
         "sizefactor)^(1/NDIMS) * resolution. hmin is the Triangle/TetGen node-spacing floor (a "
         "remesh deletes interior nodes closer than hmin to a retained node), with MMG in 2D the "
         "boundary-segment length below which a remesh merges the segment's two nodes (never a "
         "domain corner), and the necking-stop distance of remeshing_option 1, 11 and 13 (see "
         "remeshing_option). With options 1, 11 and 13 a remesh removes the nodes off a restored "
         "boundary that lie past it or within min_dist of it, and option 2 the interior nodes within "
         "min_dist of the new bottom; with options 1 and 10-12 Triangle keeps a boundary node whose "
         "two segments are both longer than min_dist, and in 3D TetGen's boundary edge flips move a "
         "boundary by at most min_dist; at a receiving wall, option 13 gives the side profile to new "
         "material up to min_dist inland of the old wall. With MMG it also sets MMG's size floors "
         "(see mmg_remesh_defensive_tiny_ratio).")
        ("mesh.largest_size", opt::value<double>(&p.mesh.largest_size)->default_value(30),
         "The size of largest element relative to an element with typical resolution, which also "
         "caps the coarse zone of a meshing_option = 2 mesh. With MMG it also sets the maximum edge "
         "length, hmax = largest_size^(1/NDIMS) * resolution")
        ("mesh.sediment_size", opt::value<double>(&p.mesh.sediment_size)->default_value(1.),
         "The size threshold of the refined sediment elements")
        // for 2D only
        ("mesh.min_angle", opt::value<double>(&p.mesh.min_angle)->default_value(32.),
         "Min. angle of all triangles (in degrees), for 2D only")
        ("mesh.max_steiner_factor", opt::value<double>(&p.mesh.max_steiner_factor)->default_value(30.),
         "Cap on the number of Steiner points Triangle may insert while (re)meshing, "
         "expressed as a multiple of the number of input nodes (2D only). This bounds the "
         "pathological, non-terminating segment splitting that Triangle can enter when two "
         "boundary segments meet at a very small angle (e.g. a highly distorted mesh at "
         "extreme strain), which would otherwise spin forever and exhaust memory. The cap is "
         "far above any healthy remesh, so normal meshing is unaffected. Set <= 0 to disable "
         "the cap (Triangle default: unlimited).")
        // for 3D only
        ("mesh.min_tet_angle", opt::value<double>(&p.mesh.min_tet_angle)->default_value(22.),
         "Min. dihedral angle of all tetrahedra (in degrees), for 3D only")
        ("mesh.max_ratio", opt::value<double>(&p.mesh.max_ratio)->default_value(2.),
         "Max. radius / length ratio of all tetrahedra, for 3D only")

        // for meshing_option = 2 only
        /* read the value as string, then parse as two numbers later */
        ("mesh.refined_zonex", opt::list()->default_value("[0.4, 0.6]"),
         "Refining portion of xlength ([d0,d1]; 0<=d0<=d1<=1), for meshing_option=2 only")
        ("mesh.refined_zoney", opt::list()->default_value("[0.4, 0.6]"),
         "Refining portion of ylength ([d0,d1]; 0<=d0<=d1<=1), for meshing_option=2 only, for 3D only")
        ("mesh.refined_zonez", opt::list()->default_value("[0.8, 1]"),
         "Refined depth interval [d0,d1] as fractions of zlength measured downward from the top "
         "(0 <= d0 <= d1 <= 1); meshing_option=2 only.")

        // for meshing_option = 90 only
        ("mesh.poly_filename", opt::value<std::string>(&p.mesh.poly_filename)->default_value("mesh.poly"),
         "Input polygon filename for meshing_option=90 or 91. Triangle .poly format in 2D; "
         "TetGen .poly format in 3D. Region attributes select materials when "
         "ic.mattype_option=0. Point attributes, point boundary markers and holes are not "
         "supported. Also supply mesh.xlength, mesh.ylength and mesh.zlength.")

        // for meshing_option = 95 only
        ("mesh.exo_filename", opt::value<std::string>(&p.mesh.exo_filename)->default_value("mesh.exo"),
         "Filename of the input polygon, for meshing_option=95 only.\n"
         "Meshing software like Trelis is supposed to be used for creating a mesh in the ExodusII format.\n"
         "Limitation: Only box-like domains with six sides are supported.\n"
         "Users still need to provide these parameters: mesh.xlength, mesh.ylength, and mesh.zlength.")

        ("mesh.quality_check_step_interval", opt::value<int>(&p.mesh.quality_check_step_interval)->default_value(100),
         "How often to check mesh quality?")
        ("mesh.min_quality", opt::value<double>(&p.mesh.min_quality)->default_value(0.4),
         "Min. mesh quality before remeshing (between 0 and 1)")
        ("mesh.max_boundary_distortion", opt::value<double>(&p.mesh.max_boundary_distortion)->default_value(0.25),
         "Max. distance of boundary distortion before remeshing (in unit of mesh.resolution).")

        ("mesh.remeshing_option", opt::value<int>(&p.mesh.remeshing_option)->default_value(0),
         "How to deal with the boundaries during remeshing?\n"
         " 0: no modification on any boundary.\n"
         " 1: move all bottom nodes to initial depth and delete the side-wall nodes left just "
         "above or below it; other boundaries are intact.\n"
         " 2: create a new bottom boundary at the initial depth, other boundaries are intact "
         "(2D only, not with MMG remeshing).\n"
         "10: no modification on any boundary, except small boundary segments might get merged "
         "(in 3D with TetGen, see smallest_size).\n"
         "11: move all bottom nodes to initial depth and delete the side-wall nodes left just "
         "above or below it; other boundaries are intact, except small boundary segments might "
         "get merged.\n"
         "12: move all x0 (left) wall nodes to the x of the wall's top corner and delete the "
         "nodes left outside, other boundaries are intact. For a wall pushed over a fixed bottom "
         "with its push tapered toward it (vbc_x0 = 3 with bottom_shear_zone_thickness), whose "
         "foot lags its top. With MMG remeshing, 2D only.\n"
         "13: move all bottom and side-wall nodes (left/right; also front/back in 3D) to initial "
         "dimensions, and give material entering through a side wall the structure recorded at that "
         "wall (its side profile). A wall receives while its nodes drift inward between remeshes (by "
         "any amount when its vbc at the remesh is inflow, else by more than min(1e-3, "
         "max_boundary_distortion) * resolution). It keeps the column (2D) or wall face (3D) "
         "recorded when it started receiving (at run start when its vbc is inflow, else at the first "
         "remesh at which it receives), and drops it when it neither receives nor has an inflow vbc. "
         "Incoming material gets the recorded temperature, layering, radiogenic source, coord0 "
         "reference and zero strain at its depth below the wall's top point less the rise surface "
         "processes gave it, and in 3D at its along-strike position (a point off the recorded face "
         "takes the values of the face's nearest point). A deposit on the wall top gets the profile "
         "top's values. Markers a remesh creates in an element at the wall take the layering at the "
         "element's centroid and time 0; in an element centred in the deposit, mattype_sediment and the "
         "time of a sediment marker of the element, else the remesh time. Temperatures are not "
         "re-conducted under a risen surface, so they can differ from the domain beside the wall by "
         "up to the near-surface temperature gradient times the rise. Requires MMG when "
         "meshing_elem_shape = 0: the restored wall accumulates sub-hmin elements that only MMG "
         "collapses. With MMG on an unstructured mesh, a Triangle/TetGen initial mesh first gets its "
         "restored-wall elements refined to ~resolution, whatever use_mmg_init says (in 3D under "
         "mattype_option 0 this can move a regional-attribute boundary inside the wall band), and "
         "every remesh caps the size target at every restored-wall node at resolution (see "
         "mmg_aniso_wall_ratio), so the walls can carry more elements than the cfg's sizes give.\n"
         "14: restore selected reference-boundary components (2D unstructured meshes only; "
         "mesh.remeshing_reset_boundaries).\n"
         "On an unstructured mesh (meshing_elem_shape = 0), options 1, 11 and 13 stop the run at a "
         "remesh (exit 42) when a top-surface node is within smallest_size^(1/NDIMS) * resolution "
         "of the restored bottom, or below it: the domain has necked through.\n")

        ("mesh.num_remeshing_reset_boundaries", opt::value<int>(&p.mesh.num_remeshing_reset_boundaries)->default_value(0),
         "Number of selected boundary-marker bits for remeshing_option=14.")
        ("mesh.remeshing_reset_boundaries", opt::list()->default_value("[]"),
         "Boundary-marker bits to restore, e.g. [1,16]. Single faces restore normal components; "
         "two selected faces restore their corner fully.")

        ("mesh.mmg_remesh_active_plstrain", opt::value<double>(&p.mesh.mmg_remesh_active_plstrain)->default_value(0.01),
         "Conservative MMG remeshing (always on with MMG): a remesh keeps quiet elements as they are "
         "(frozen) and rebuilds the rest. This key is the plastic-strain increment since the last "
         "remesh above which an element is refined: MMG may split it (not along an edge it shares "
         "with a frozen, MMG-required element) but keeps its existing nodes, so the band's plastic "
         "strain is "
         "carried verbatim. Distorted or tiny elements are repaired with node motion (see "
         "mmg_remesh_defensive_quality_ratio and mmg_remesh_defensive_tiny_ratio), and so is the "
         "mesh that restoring a boundary squashed or collapsed or that merging a short boundary "
         "segment (see smallest_size) reshaped. Elements far "
         "from their size target are freed (see mmg_remesh_size_recovery_ratio), as is the band beside "
         "a receiving or outward-drifting restored wall (remeshing_option 13). Boundary nodes stay "
         "pinned except the free-surface nodes a surface process moves by more than 1% of resolution "
         "per remesh interval (2D only, surface_process_option 1 or 102), the off-bottom nodes of "
         "unrestored side walls and slant boundaries, the nodes of repaired elements and the band "
         "nodes off an outward-drifting restored wall. Keep it > "
         "0: at 0, incidental yielding makes quiet far-field elements splittable at every remesh. "
         "Default 0.01, heuristic: above the increment of incidental yielding, below a localizing "
         "band's increment over a remesh interval.")
        ("mesh.mmg_remesh_defensive_tiny_ratio", opt::value<double>(&p.mesh.mmg_remesh_defensive_tiny_ratio)->default_value(1.5),
         "Gap, in element measure, between the tiny-element remesh trigger smallest_vol (see "
         "smallest_size) and the smallest elements MMG builds (usemmg=1, meshing_elem_shape = 0; "
         "remeshing and MMG-built initial meshes). MMG's minimum edge length is hmin * "
         "ratio^(1/NDIMS) (hmin: see smallest_size), so MMG aims to build every element at "
         "smallest_vol * ratio or more, which must shrink by this factor before it re-triggers, and "
         "a remesh rebuilds any element below "
         "smallest_vol * ratio. The ratio also sets MMG's metric floor one step higher, hmin * "
         "ratio^(2/NDIMS) (at most MMG's hmax), and the post-remesh quality gate at smallest_vol * "
         "sqrt(ratio). Must be >= 1 (at 1 a just-remeshed element re-triggers under any "
         "compression); a ratio at or above largest_size / smallest_size is capped just below it, "
         "with a warning. Default 1.5 (a 2/3 measure gap), heuristic.")
        ("mesh.mmg_remesh_defensive_quality_ratio", opt::value<double>(&p.mesh.mmg_remesh_defensive_quality_ratio)->default_value(1.5),
         "MMG remeshing repairs (unfreezes, with its nodes movable except where disjoint material "
         "bodies meet) every element of quality below "
         "min_quality * this ratio, so a marginal element is rebuilt before it crosses the "
         "min_quality remesh trigger instead of triggering another remesh a few steps later. The "
         "trigger itself stays at min_quality. Must be >= 1 (1: repair only below the trigger). "
         "Default 1.5, heuristic.")
        ("mesh.mmg_remesh_size_recovery_ratio", opt::value<double>(&p.mesh.mmg_remesh_size_recovery_ratio)->default_value(1.4),
         "Size recovery in MMG remeshing. A quiet element (e.g. one that entered at ~resolution "
         "through a restored inflow wall and advected inland) is freed to coarsen when it is finer "
         "by more than this factor than the gradation envelope, the smallest size MMG's gradation "
         "(hgrad) allows there given every intentionally fine size target; an element at least 2x "
         "its own plastic-strain size target is freed to split. Neither frees an element refined for "
         "plastic-strain growth (mmg_remesh_active_plstrain), one on the top surface, one at a "
         "material interface (by its own and its neighbours' markers) or one within 4 rings of a "
         "restored wall that receives no material, and neither unpins a restored-wall or bottom "
         "node; the remesh log's '[unfreeze] size recovery' line counts the top-surface, interface "
         "and wall-band skips. On a moving "
         "mesh, a Triangle/TetGen initial mesh also gets the 2x split in its weak zone before any "
         "field exists, except on a top-surface node or a region boundary ('[weakzone]' log lines). "
         "0 disables all three; otherwise must be > 1 (a margin against re-freeing the legal size "
         "transition). Default 1.4 ~ sqrt(2) in edge length (one 2D bisection level), heuristic.")

        ("mesh.remesh_deborah_min", opt::value<double>(&p.mesh.remesh_deborah_min)->default_value(1e0),
         "During remeshing the element stress is a Deborah-number-weighted blend of "
         "the NN-remapped stress and the SPR recovery, De = Maxwell time / time since "
         "the last remesh. The weight is a smoothstep -- the cubic 3t^2 - 2t^3, flat "
         "at both ends -- in log10(De): high-De (cold, effectively "
         "elastic) elements keep the NN stress, preserving elastic stress memory; "
         "low-De (weak, viscous) elements take the SPR average, which suppresses "
         "element-scale noise. This sets the De at or below which the remeshed "
         "stress is purely the SPR recovery.")

        ("mesh.remesh_deborah_max", opt::value<double>(&p.mesh.remesh_deborah_max)->default_value(1e2),
         "De at or above which the remeshed stress is purely the NN-remapped stress. "
         "The default puts the two limits two decades apart, which keeps them clear of "
         "the De = 1 crossover while remesh intervals themselves vary by a factor of a "
         "few between events; widen it to blend more of the domain, narrow it to make "
         "the switch sharper. The window width is otherwise untuned.")

        ("mesh.is_discarding_internal_segments", opt::value<bool>(&p.mesh.is_discarding_internal_segments)->default_value(true),
         "Discard the initial mesh's internal segments (refined-zone box, .poly internal facets)? "
         "'no' keeps them as Triangle/TetGen constraints across remeshes, but a remesh that deletes "
         "one of their nodes (tiny-element repair, node cleanup near a restored boundary) misroutes "
         "them, which can make the mesher fail or hang; in 3D with remeshing_option 1 or 10-12 the "
         "first remesh that deletes a boundary node drops them all. MMG remeshing (usemmg=1, "
         "meshing_elem_shape = 0) requires 'yes': MMG keeps no internal segment.")

        // for mesh optimization with MMG2D/3D
        ("mesh.mmg_debug", opt::value<int>(&p.mesh.mmg_debug)->default_value(0),
         "Run MMG remesher in debug mode? No:0; Yes:1\n")
        ("mesh.mmg_verbose", opt::value<int>(&p.mesh.mmg_verbose)->default_value(0),
         "Verbosity level of MMG remesher. For debugging, set a value greater than 4.\n")
        ("mesh.mmg_hausd_factor", opt::value<double>(&p.mesh.mmg_hausd_factor)->default_value(0.01),
         "Factor multiplied to param.mesh.resolution to set the Hausdorff distance between original and remeshed surfaces. "
         "At 0.01 MMG refines a sharply bent free surface down to its size floor and builds elements "
         "below the tiny-element trigger there, so the trigger fires at nearly every quality check: "
         "use 0.02-0.05 along a 3D fault scarp and 0.1 on the top/wall ridges next to restored side "
         "walls. MMG3D does not bound its top-surface edge swaps between pinned nodes by this "
         "distance.\n")
        ("mesh.mmg_init_coarsening_factor", opt::value<double>(&p.mesh.mmg_init_coarsening_factor)->default_value(10.0),
         "Coarsening factor for the two-stage TetGen+MMG init mesh. TetGen creates a mesh "
         "with element volumes scaled by factor^NDIMS; MMG then refines to target resolution. "
         "Higher values = smaller coarse mesh = faster TetGen but more MMG work. Default: 10.")
        ("mesh.use_mmg_init", opt::value<bool>(&p.mesh.use_mmg_init)->default_value(false),
         "Use two-stage TetGen+MMG initialization? Requires compilation with USEMMG. "
         "Set to yes for metric-aware refinement; no (default) uses plain triangle/tetgen init "
         "(MMG still refines its restored walls under remeshing_option 13 and, on a moving mesh, "
         "its weak zone: see mmg_remesh_size_recovery_ratio). Under "
         "mattype_option 0, which reads the regional attributes as mattypes, a boundary between "
         "them stays on the .poly line in 2D and can move by up to about half a coarse element in "
         "3D (see mmg_init_coarsening_factor).")
        ("mesh.mmg_metric_refine_coeff", opt::value<double>(&p.mesh.mmg_metric_refine_coeff)->default_value(5.0),
         "Plastic-strain sensitivity of the MMG remeshing metric: the target element volume is "
         "scaled by 1/(1 + coeff*plastic_strain), so the target edge length by that ratio^(1/NDIMS) "
         "(whether a node keeps its finest target: see mmg_metric_refine_recover). 0 disables "
         "plastic refinement from "
         "a fresh start (the metric floor, the surface-motion inflation of shrinking top elements "
         "and the option-13 restored-wall cap still apply). MMG starts from this target except at a "
         "node on a required edge (an edge of a frozen element or, in 2D, an option-13 "
         "material-interface edge) that belongs to no refined element (see "
         "mmg_remesh_active_plstrain): there it starts from the mean length of the required edges "
         "through the node, MMG's own rule. Larger values mesh a plastically strained zone finer, "
         "down to the metric floor (see mmg_remesh_defensive_tiny_ratio): a shear band can localize "
         "narrower, at more elements and, where an elastic, advection or diffusion limit sets dt, a "
         "shorter time step. Default 5.")
        ("mesh.mmg_metric_refine_recover",
         opt::value<bool>(&p.mesh.mmg_metric_refine_recover)->default_value(false),
         "no (default): each node keeps the finest plastic-strain size target (see "
         "mmg_metric_refine_coeff) it has had at a remesh, so a refined zone stays refined where the "
         "target coarsens, chiefly where plastic strain drops (the remap of rebuilt elements, the 2D "
         "top-element plastic-strain decay, a growing top element); the mesh can then hold more "
         "elements. Size recovery (mmg_remesh_size_recovery_ratio) measures against the kept size, and "
         "material entering through a restored side wall starts unrefined. yes: the target follows "
         "the current plastic strain. The record is checkpointed: with no, a restart that lowers "
         "mmg_metric_refine_coeff (even to 0) keeps the sizes recorded before it, and a restart that "
         "switches from yes to no holds from the last remesh's targets.")
        ("mesh.mmg_aniso_wall_ratio", opt::value<double>(&p.mesh.mmg_aniso_wall_ratio)->default_value(0.0),
         "Anisotropic MMG metric at inflow restored side walls (remeshing_option 13). 0 (default): "
         "off, isotropic. >= 1: aspect ratio; the inflow-wall elements are coarsened by this factor "
         "normal to the wall and kept fine along it, so the incoming column is carried by thin, wide "
         "elements. MMG's hmax (see largest_size) caps the aspect (~5.5 in 2D, "
         "~3.1 in 3D at the default largest_size). Other values below 1 are rejected. Requires "
         "USEMMG.")

        // for virtual topology
        ("mesh.vt_size", opt::value<double>(&p.mesh.vt_size)->default_value(0.5),
         "Cell size of the virtual-topology marker grid (in unit of mesh.resolution), used by the "
         "serpentinization grid (mat.phase_change_option 3 with mat.has_serpentinization). 0.5 = half an "
         "element edge: the one-element-thick slab crust "
         "gets its own wet-crust-dominant cell row instead of averaging with the mantle above, while a "
         "cell nominally (elements of edge mesh.resolution) holds ~2.3 markers in 2D (~4.2 in 3D) at the "
         "default 4 per element; a refined mesh puts more in each cell. Heuristic, untuned. With serpentinization "
         "on it must be > 0, and mesh.resolution * mesh.vt_size must not exceed mat.serpentinization_thickness. "
         "The grid is rebuilt at each mesh.quality_check_step_interval and after each remesh.")
        ;

    cfg.add_options()
        ("markers.init_marker_option", opt::value<int>(&p.markers.init_marker_option)->default_value(1),
         "How to generate markers?\n"
         "1: randomly distributed markers.\n"
         "2: regularly spaced markers.\n")
        ("markers.markers_per_element", opt::value<int>(&p.markers.markers_per_element)->default_value(4),
         "Number of markers per element. Used when init_marker_option=1.")
        ("markers.init_marker_spacing", opt::value<double>(&p.markers.init_marker_spacing)->default_value(0.3),
         "Spacing of markers (in unit of mesh.resolution). Used when init_marker_option=2. An "
         "element the grid misses gets min_num_markers_in_element (at least 1) markers of its "
         "initial material, with a warning.")
        ("markers.min_num_markers_in_element", opt::value<int>(&p.markers.min_num_markers_in_element)->default_value(3),
         "When the number of markers in an element is less than this number, a new marker will be replenished in the element.")
        ("markers.replenishment_option", opt::value<int>(&p.markers.replenishment_option)->default_value(2),
         "How to determine the mattype of replenished markers?\n"
         "0: always set to 0 (fastest option).\n"
         "1: by the probability of marker mattype of the element or surrounding elements.\n"
         "2: same as the mattype of the nearest marker (slowest option).\n"
         "12: if the element's surviving markers are all one mattype, use that mattype directly "
         "(deposit-time interpolated for mattype_sediment, as in option 2); otherwise fall back "
         "to option 2 (nearest marker). Exact in single-material elements.\n"
         "At a remesh, independently of this option, an element whose part inside the old mesh lies "
         "in one old element whose markers were all one mattype gets that mattype, and with "
         "remeshing_option 13 an element at a receiving wall gets the wall's recorded layering (see "
         "remeshing_option).")
        ("markers.random_seed", opt::value<uint>(&p.markers.random_seed)->default_value(1),
         "Random seed of marker position. If 0, the current time is used as the seed.")
        ;

    cfg.add_options()
        ("control.gravity", opt::value<double>(&p.control.gravity)->default_value(10),
         "Magnitude of the gravity (in m/s^2)")

        ("control.characteristic_speed",
         opt::value<double>(&p.control.characteristic_speed)->default_value(0),
         "Characteristic tectonic speed (in m/s). "
         "It is used to estimate the size of stable time step. "
         "If it is 0, its value is inferred from the imposed boundary velocity. "
         "Set it to other value if the imposed boundary velocity is 0 everywhere.")

        ("control.is_quasi_static", opt::value<bool>(&p.control.is_quasi_static)->default_value(true),
         "Is the simulation quasi-static or dynamic? If quasi-static, inertial scaling and strong damping is applied.\n")
        ("control.dt_fraction", opt::value<double>(&p.control.dt_fraction)->default_value(1.0),
         "Take dt as a fraction of max. stable time step size (0-1).\n")
        ("control.fixed_dt", opt::value<double>(&p.control.fixed_dt)->default_value(0),
         "Fixed dt size (in seconds). If 0, dt sized will be determined dynamically.\n")
        ("control.inertial_scaling", opt::value<double>(&p.control.inertial_scaling)->default_value(1e5),
         "Scaling factor for inertial (a large number)")

        ("control.damping_option", opt::value<int>(&p.control.damping_option)->default_value(1),
         "How to damp the elastic wave?\n"
         "0: no damping.\n"
         "1: damping/acceleration depends on force and direction of velocity.\n"
         "2: damping depends on force.\n"
         "3: damping/(weaker) acceleration depends on force and direction of velocity.\n"
         "4: Rayleigh damping.\n")
        ("control.damping_factor", opt::value<double>(&p.control.damping_factor)->default_value(0.8),
         "A factor for force damping (0-1)")

        ("control.ref_pressure_option", opt::value<int>(&p.control.ref_pressure_option)->default_value(0),
         "How to define reference pressure?\n"
         "0: using density of the mat.mattype_ref-th material to compute lithostatic pressure.\n"
         "1: computing reference pressure from the PREM model.\n"
         "2: computing reference pressure from the PREM model, modified for continent.\n"
         "Any other value is rejected at startup.\n")
        ("control.is_using_mixed_stress", opt::value<bool>(&p.control.is_using_mixed_stress)->default_value(true),
         "If use Nodal Mixed Discretization For Stress")

        ("control.surface_process_option", opt::value<int>(&p.control.surface_process_option)->default_value(0),
         "What kind of surface processes? 0: no surface processes. "
         "1: using simple diffusion to modify surface topography. "
         "2: using simple deposition. "
         "11: using GoSPL landscape evolution model. "
         "101: custom function. "
         "102: using simple deposition and diffusion.")
        ("control.surface_diffusivity", opt::value<double>(&p.control.surface_diffusivity)->default_value(1e-6),
         "Standard diffusition coefficient of surface topography (m^2/s)")
        ("control.surf_diff_ratio_terrig", opt::value<double>(&p.control.surf_diff_ratio_terrig)->default_value(1.),
         "Ratio between terrigenous and standard Diffusition coefficients.")
        ("control.surf_diff_ratio_marine", opt::value<double>(&p.control.surf_diff_ratio_marine)->default_value(1.),
         "Ratio between marine and standard Diffusition coefficients.")
        ("control.surf_depo_universal", opt::value<double>(&p.control.surf_depo_universal)->default_value(0.),
         "Sedimentation rate of universal suspended deposite (in m/s).")
        ("control.surf_base_level",opt::value<double>(&p.control.surf_base_level)->default_value(0.e0),
         "Base level of surface processes (m); also the sea level used by bc.has_water_loading.")
        ("control.surface_process_gospl_config_file", opt::value<std::string>(&p.control.surface_process_gospl_config_file)->default_value(""),
         "Configuration file path for GoSPL surface processes (used with surface_process_option = 11).")
        ("control.gospl_coupling_frequency", opt::value<int>(&p.control.gospl_coupling_frequency)->default_value(1),
         "Run GoSPL coupling every N DES steps (used when gospl_coupling_mode = steps).")
        ("control.gospl_coupling_interval_in_yr", opt::value<double>(&p.control.gospl_coupling_interval_in_yr)->default_value(1000.0),
         "Run GoSPL coupling every T years (used when gospl_coupling_mode = time). E.g. 1000 = couple once per 1000 yr.")
        ("control.gospl_coupling_mode", opt::value<std::string>(&p.control.gospl_coupling_mode)->default_value("steps"),
         "Coupling frequency mode: 'steps' (use gospl_coupling_frequency) or 'time' (use gospl_coupling_interval_in_yr).")
        ("control.gospl_velocity_coupling", opt::value<bool>(&p.control.gospl_velocity_coupling)->default_value(true),
         "Send all three DES surface velocity components (vx, vy, vz) to GoSPL each coupling step. "
         "vz drives uplift/subsidence; vx/vy drive semi-Lagrangian horizontal advection of GoSPL topography.")
        ("control.gospl_mesh_resolution", opt::value<double>(&p.control.gospl_mesh_resolution)->default_value(-1.0),
         "GoSPL mesh node spacing in meters (default: -1 = auto-size based on DES surface nodes). Set to desired resolution for control over GoSPL grid density.")
        ("control.gospl_mesh_perturbation", opt::value<double>(&p.control.gospl_mesh_perturbation)->default_value(0.3),
         "Fraction of grid spacing to randomly perturb node positions (0-1). 0 = regular grid, 0.3 = moderate perturbation to avoid grid artifacts.")
        ("control.gospl_mesh_padding", opt::value<double>(&p.control.gospl_mesh_padding)->default_value(0.1),
         "Fractional extension of the GoSPL mesh beyond the DES domain on each side (default: 0.1 = 10%). Keeps DES nodes in the interior of GoSPL's mesh, away from boundary artifacts.")
        ("control.terrig_sediment_diffusivity",opt::value<double>(&p.control.terrig_sediment_diffusivity)->default_value(3.17e-6),
         "Submarine diffusion coefficient (=100 m^2/yr, Kaufman et al., 1991) (in m^2/s).")
        ("control.terrig_depth_coefficient",opt::value<double>(&p.control.terrig_depth_coefficient)->default_value(5e-4),
         "Submarine diffusion coefficient decay (Kaufman et al., 1991) (in m^-1).")
        ("control.is_reporting_terrigenous_info",opt::value<bool>(&p.control.is_reporting_terrigenous_info)->default_value(false),
         "If report surface terrigenous information.")
        ("control.hemipelagic_sedimentation_rate",opt::value<double>(&p.control.hemipelagic_sedimentation_rate)->default_value(3.17e-11),
         "Hemipelagic sedimentation (Emmerich et al., 2009) (1 km/Myr ~= 3.17e-11 m/s)")
        ("control.hemipelagic_width",opt::value<double>(&p.control.hemipelagic_width)->default_value(90.),
         "Width of the productivity function (m)")
        ("control.hemipelagic_max_depth",opt::value<double>(&p.control.hemipelagic_max_depth)->default_value(40.),
         "Maximum accumulation bathymetry (m)")
        ("control.pelagic_sedimentation_rate",opt::value<double>(&p.control.pelagic_sedimentation_rate)->default_value(1.585e-12),
         "Pelagic sedimentation rate (Emmerich et al., 2009) (0.05 km/Myr ~= 1.585e-12 m/s)")
        ("control.pelagic_increasing_width",opt::value<double>(&p.control.pelagic_increasing_width)->default_value(30.),
         "Width of the suface increasing zone (m)")

        ("control.terrig_sediment_volume",opt::value<double>(&p.control.terrig_sediment_volume)->default_value(1.),
         "The volume of the sediment bring from continent source (in m^3/s). "
         "5e7 ton/yr ~= 1 m^3/s")
        ("control.terrig_sediment_area",opt::value<double>(&p.control.terrig_sediment_area)->default_value(3.17e-7),
         "The are of the sediment bring from continent source (=10 m^2/yr) (in m^2/s). "
         "Assuming the width of the channel is 50 km, 5e7 ton/yr ~= 2.e-5 m^2/s")

        ("control.has_thermal_pressurization", opt::value<bool>(&p.control.has_thermal_pressurization)->default_value(false),
         "Couple thermal expansion to solid stress and fluid content; requires hydraulic and thermal diffusion.")
        ("control.has_shear_heating", opt::value<bool>(&p.control.has_shear_heating)->default_value(false),
         "Convert accepted viscous/plastic work to heat. Pending energy is deposited at the next "
         "physical thermal update, independently of changes in timestep; initial PT equilibration "
         "does not generate heat. All computed nonelastic work is converted (unit heat fraction).")
        ("control.has_fluid_thermal_advection", opt::value<bool>(&p.control.has_fluid_thermal_advection)->default_value(false),
         "Transport heat with Darcy flow using combined rock/water heat capacity; requires hydraulics.")
        ("control.has_poroelastic_pressure_feedback", opt::value<bool>(&p.control.has_poroelastic_pressure_feedback)->default_value(true),
         "Include mechanical volumetric/mean-stress changes in the hydraulic source; transport remains active when false.")
        ("control.has_thermal_diffusion", opt::value<bool>(&p.control.has_thermal_diffusion)->default_value(true),
         "Does the model have thermal diffusion? If not, temperature is advected, but not diffused.\n")
        ("control.has_hydraulic_diffusion", opt::value<bool>(&p.control.has_hydraulic_diffusion)->default_value(false),
         "Does the model have hydraulic diffusion? If not, pore pressure is advected, but not diffused.\n") 
        ("control.hydraulic_content_model", opt::value<int>(&p.control.hydraulic_content_model)->default_value(0),
         "Fluid content: 0 linear Biot/retention, 1 finite-density saturated mass balance with evolving Biot porosity.")

        ("control.has_hydration_processes", opt::value<bool>(&p.control.has_hydration_processes)->default_value(false),
         "Does the model have hydration processes? It is required to model some types of phase changes.")
        ("control.hydration_migration_speed", opt::value<double>(&p.control.hydration_migration_speed)->default_value(3e-9),
         "The upward migration speed of hydrous fluid (in m/s).\n")

        ("control.has_PT", opt::value<bool>(&p.control.has_PT)->default_value(false),
         "Does the model have Pseudo-transient (PT) loop?\n")
        ("control.PT_use_running_scale", opt::value<bool>(&p.control.PT_use_running_scale)->default_value(false),
         "Use running physical initial-residual scale with a roundoff floor (epsilon=1e-6, C=1000)")
        ("control.PT_max_iter", opt::value<int>(&p.control.PT_max_iter)->default_value(5000),
         "Maximum iteration for PT loop")
        ("control.PT_relative_tolerance",opt::value<double>(&p.control.PT_relative_tolerance)->default_value(1e-6),
         "Relative tolerance on projected physical-force RMS, scaled by the first warm trial")
        ("control.PT_absolute_tolerance", opt::value<double>(&p.control.PT_absolute_tolerance)->default_value(0.0),
         "Absolute projected-force RMS tolerance (N/m in 2D, N in 3D)")
        ("control.PT_option", opt::value<int>(&p.control.PT_option)->default_value(0),
         "PT iteration scheme. Both accept the same physical residual.\n"
         "0: accelerated pseudo-transient relaxation (Rass et al. 2022), uses PT_CFL and PT_Re.\n"
         "1: adaptive dynamic relaxation: Gershgorin fictitious mass and damping from a\n"
         "   Rayleigh-quotient estimate of the lowest mode (Underwood 1983; Papadrakakis 1981).\n")
        ("control.PT_CFL", opt::value<double>(&p.control.PT_CFL)->default_value(0.25),
         "Pseudo-wave CFL factor for local dual-time relaxation (0 < CFL <= 0.5)")
        ("control.PT_Re", opt::value<double>(&p.control.PT_Re)->default_value(14.90188239869415),
         "Positive pseudo Reynolds number for dual-time relaxation")
        ("control.PT_stagnation_window", opt::value<int>(&p.control.PT_stagnation_window)->default_value(500),
         "Warn once per solve after this many iterations without best-residual improvement; 0 disables warning")

         ("control.has_moving_mesh", opt::value<bool>(&p.control.has_moving_mesh)->default_value(true),
         "Does the model update mesh coordinates (Lagrangian)?\n")
         ("control.use_global_velocity_scaling",
          opt::value<bool>(&p.control.use_global_velocity_scaling)->default_value(false),
          "Use the global maximum model velocity to scale both dt and pseudo-density/mass scaling.\n")
        ("control.mass_scaling_reference_speed",
         opt::value<std::string>()->default_value("shear"),
         "Elastic-speed ceiling used by global velocity scaling. It sets the floor of\n"
         "the fictitious density rho_fict = K / v_elastic^2.\n"
         "shear: sqrt(G/rho), giving rho_fict >= rho K/G (historical default).\n"
         "bulk : sqrt(K/rho), giving rho_fict >= rho and restoring physical density\n"
         "       when the pseudo-wave speed reaches the bulk-wave ceiling.\n")
        ("control.rsf_slip_rate_projection_option",
         opt::value<int>(&p.control.rsf_slip_rate_projection_option)
             ->default_value(rsf_slip_rate_projection_maximum_shear),
         "Velocity-dimensional rate supplied to the RSF update.\n"
         "0: project the element velocity onto the maximum-shear direction inferred from stress\n"
         "   (default; the historical behavior).\n"
         "1: use V = 2 w eps_II from the total deviatoric strain rate, where w is the\n"
         "   element's minimum altitude.\n")
        ("control.rsf_dtheta_max",
         opt::value<double>(&p.control.rsf_dtheta_max)->default_value(0.0),
         "For adaptive stepping with the aging law and total-strain rate option 1, apply\n"
         "dt <= f D_c / V and dt <= f theta over all elements. The two bounds limit\n"
         "slip within one characteristic distance and fractional healing per step.\n"
         "f must be in [0, 2); 0 (default) disables this state-update limit.\n")
        ;

    cfg.add_options()
        ("bc.tbc_x0", opt::value<int>(&p.bc.tbc[0])->default_value(0),
         "Thermal boundary: 0 insulated, 1 fixed temperature.")
        ("bc.tbc_val_x0", opt::value<double>(&p.bc.tbc_value[0])->default_value(std::numeric_limits<double>::quiet_NaN()),
         "Fixed boundary temperature in K; defaults to surface_temperature (mantle_temperature at z0). "
         "It is applied at initialization, so a value that differs from the initial temperature is part "
         "of the initial state, not a thermal load; to heat a boundary at t > 0, restart with the new value.")
        ("bc.tbc_x1", opt::value<int>(&p.bc.tbc[1])->default_value(0),
         "Thermal boundary: 0 insulated, 1 fixed temperature.")
        ("bc.tbc_val_x1", opt::value<double>(&p.bc.tbc_value[1])->default_value(std::numeric_limits<double>::quiet_NaN()),
         "Fixed boundary temperature in K; defaults to surface_temperature (mantle_temperature at z0). "
         "It is applied at initialization, so a value that differs from the initial temperature is part "
         "of the initial state, not a thermal load; to heat a boundary at t > 0, restart with the new value.")
        ("bc.tbc_y0", opt::value<int>(&p.bc.tbc[2])->default_value(0),
         "Thermal boundary: 0 insulated, 1 fixed temperature.")
        ("bc.tbc_val_y0", opt::value<double>(&p.bc.tbc_value[2])->default_value(std::numeric_limits<double>::quiet_NaN()),
         "Fixed boundary temperature in K; defaults to surface_temperature (mantle_temperature at z0). "
         "It is applied at initialization, so a value that differs from the initial temperature is part "
         "of the initial state, not a thermal load; to heat a boundary at t > 0, restart with the new value.")
        ("bc.tbc_y1", opt::value<int>(&p.bc.tbc[3])->default_value(0),
         "Thermal boundary: 0 insulated, 1 fixed temperature.")
        ("bc.tbc_val_y1", opt::value<double>(&p.bc.tbc_value[3])->default_value(std::numeric_limits<double>::quiet_NaN()),
         "Fixed boundary temperature in K; defaults to surface_temperature (mantle_temperature at z0). "
         "It is applied at initialization, so a value that differs from the initial temperature is part "
         "of the initial state, not a thermal load; to heat a boundary at t > 0, restart with the new value.")
        ("bc.tbc_z0", opt::value<int>(&p.bc.tbc[4])->default_value(0),
         "Thermal boundary: 0 insulated, 1 fixed temperature.")
        ("bc.tbc_val_z0", opt::value<double>(&p.bc.tbc_value[4])->default_value(std::numeric_limits<double>::quiet_NaN()),
         "Fixed boundary temperature in K; defaults to surface_temperature (mantle_temperature at z0). "
         "It is applied at initialization, so a value that differs from the initial temperature is part "
         "of the initial state, not a thermal load; to heat a boundary at t > 0, restart with the new value.")
        ("bc.tbc_z1", opt::value<int>(&p.bc.tbc[5])->default_value(1),
         "Thermal boundary: 0 insulated, 1 fixed temperature.")
        ("bc.tbc_val_z1", opt::value<double>(&p.bc.tbc_value[5])->default_value(std::numeric_limits<double>::quiet_NaN()),
         "Fixed boundary temperature in K; defaults to surface_temperature (mantle_temperature at z0). "
         "It is applied at initialization, so a value that differs from the initial temperature is part "
         "of the initial state, not a thermal load; to heat a boundary at t > 0, restart with the new value.")
        ("bc.hbc_val_x0", opt::value<double>(&p.bc.hbc_legacy_value[0])->default_value(std::numeric_limits<double>::quiet_NaN()),
         "Research-input compatibility: finite value with hbc=1 maps to prescribed pressure (hbc=2).")
        ("bc.hbc_val_x1", opt::value<double>(&p.bc.hbc_legacy_value[1])->default_value(std::numeric_limits<double>::quiet_NaN()),
         "Research-input compatibility: finite value with hbc=1 maps to prescribed pressure (hbc=2).")
        ("bc.hbc_val_y0", opt::value<double>(&p.bc.hbc_legacy_value[2])->default_value(std::numeric_limits<double>::quiet_NaN()),
         "Research-input compatibility: finite value with hbc=1 maps to prescribed pressure (hbc=2).")
        ("bc.hbc_val_y1", opt::value<double>(&p.bc.hbc_legacy_value[3])->default_value(std::numeric_limits<double>::quiet_NaN()),
         "Research-input compatibility: finite value with hbc=1 maps to prescribed pressure (hbc=2).")
        ("bc.hbc_val_z0", opt::value<double>(&p.bc.hbc_legacy_value[4])->default_value(std::numeric_limits<double>::quiet_NaN()),
         "Research-input compatibility: finite value with hbc=1 maps to prescribed pressure (hbc=2).")
        ("bc.hbc_val_z1", opt::value<double>(&p.bc.hbc_legacy_value[5])->default_value(std::numeric_limits<double>::quiet_NaN()),
         "Research-input compatibility: finite value with hbc=1 maps to prescribed pressure (hbc=2).")
        ("bc.surface_temperature", opt::value<double>(&p.bc.surface_temperature)->default_value(273),
         "Surface temperature (in Kelvin)")
        ("bc.mantle_temperature", opt::value<double>(&p.bc.mantle_temperature)->default_value(1600),
         "Mantle temperature (in Kelvin)")
        ("bc.has_winkler_foundation", opt::value<bool>(&p.bc.has_winkler_foundation)->default_value(true),
         "Using Winkler foundation for the bottom boundary?")
        ("bc.winkler_delta_rho", opt::value<double>(&p.bc.winkler_delta_rho)->default_value(0),
         "Excess density of the bottom Winkler foundation (in kg/m^3)")

        ("bc.has_elastic_foundation", opt::value<bool>(&p.bc.has_elastic_foundation)->default_value(false),
         "Using elastic foundation for the bottom boundary?")
        ("bc.elastic_foundation_constant", opt::value<double>(&p.bc.elastic_foundation_constant)->default_value(1e11),
         "Elastic constant for elastic foundation.")

        ("bc.has_water_loading", opt::value<bool>(&p.bc.has_water_loading)->default_value(true),
         "Applying water loading for top boundary that is below sea level?")

        ("bc.sea_water_density", opt::value<double>(&p.bc.sea_water_density)->default_value(1030),
         "Density of the loading water (in kg/m^3). Used for the top-boundary "
         "water load and for the free-surface stress pin that must match it.")

         // pore pressure boundary condition
        ("bc.hbc_x0", opt::value<int>(&p.bc.hbc_x0)->default_value(0),
         "Hydraulic boundary type. 0: no normal flow; 1: hold the initial nodal pressure; 2: "
         "prescribed bc.hbc_pressure entry (Pa), multiplied by the boundary schedule; 3: "
         "prescribed bc.hbc_flux entry (m/s, positive outward), multiplied by the boundary "
         "schedule. Lists use x0 x1 y0 y1 z0 z1 order; y sides apply only in 3D.")
        ("bc.hbc_x1", opt::value<int>(&p.bc.hbc_x1)->default_value(0),
         "Hydraulic boundary type. 0: no normal flow; 1: hold the initial nodal pressure; 2: "
         "prescribed bc.hbc_pressure entry (Pa), multiplied by the boundary schedule; 3: "
         "prescribed bc.hbc_flux entry (m/s, positive outward), multiplied by the boundary "
         "schedule. Lists use x0 x1 y0 y1 z0 z1 order; y sides apply only in 3D.")
        ("bc.hbc_y0", opt::value<int>(&p.bc.hbc_y0)->default_value(0),
         "Hydraulic boundary type. 0: no normal flow; 1: hold the initial nodal pressure; 2: "
         "prescribed bc.hbc_pressure entry (Pa), multiplied by the boundary schedule; 3: "
         "prescribed bc.hbc_flux entry (m/s, positive outward), multiplied by the boundary "
         "schedule. Lists use x0 x1 y0 y1 z0 z1 order; y sides apply only in 3D.")
        ("bc.hbc_y1", opt::value<int>(&p.bc.hbc_y1)->default_value(0),
         "Hydraulic boundary type. 0: no normal flow; 1: hold the initial nodal pressure; 2: "
         "prescribed bc.hbc_pressure entry (Pa), multiplied by the boundary schedule; 3: "
         "prescribed bc.hbc_flux entry (m/s, positive outward), multiplied by the boundary "
         "schedule. Lists use x0 x1 y0 y1 z0 z1 order; y sides apply only in 3D.")
        ("bc.hbc_z0", opt::value<int>(&p.bc.hbc_z0)->default_value(0),
         "Hydraulic boundary type. 0: no normal flow; 1: hold the initial nodal pressure; 2: "
         "prescribed bc.hbc_pressure entry (Pa), multiplied by the boundary schedule; 3: "
         "prescribed bc.hbc_flux entry (m/s, positive outward), multiplied by the boundary "
         "schedule. Lists use x0 x1 y0 y1 z0 z1 order; y sides apply only in 3D.")
        ("bc.hbc_z1", opt::value<int>(&p.bc.hbc_z1)->default_value(0),
         "Hydraulic boundary type. 0: no normal flow; 1: hold the initial nodal pressure; 2: "
         "prescribed bc.hbc_pressure entry (Pa), multiplied by the boundary schedule; 3: "
         "prescribed bc.hbc_flux entry (m/s, positive outward), multiplied by the boundary "
         "schedule. Lists use x0 x1 y0 y1 z0 z1 order; y sides apply only in 3D.")

        ("bc.hbc_pressure", opt::list()->default_value("[0,0,0,0,0,0]"),
         "Pressure in Pa for hbc=2, ordered x0 x1 y0 y1 z0 z1.")
        ("bc.hbc_flux", opt::list()->default_value("[0,0,0,0,0,0]"),
         "Outward Darcy flux in m/s for hbc=3, ordered x0 x1 y0 y1 z0 z1.")
        ("bc.num_hbc_schedule", opt::value<int>(&p.bc.num_hbc_schedule)->default_value(1),
         "Number of shared boundary schedule points.")
        ("bc.hbc_schedule_time_in_yr", opt::list()->default_value("[0]"),
         "Boundary schedule times in years, nondecreasing; repeated times are steps.")
        ("bc.hbc_schedule_ratio", opt::list()->default_value("[1]"),
         "Ratios multiplying prescribed pressure and flux.")

         // rate-controlled wells
        ("bc.num_wells", opt::value<int>(&p.bc.num_wells)->default_value(0),
         "Number of vertical rate-controlled wells. Requires control.has_hydraulic_diffusion.")
        ("bc.well_x", opt::list()->default_value("[]"),
         "x coordinate of each well (in m).")
        ("bc.well_y", opt::list()->default_value("[]"),
         "y coordinate of each well (in m, 3D only).")
        ("bc.well_z_top", opt::list()->default_value("[]"),
         "Top of each well screen (in m).")
        ("bc.well_z_bottom", opt::list()->default_value("[]"),
         "Bottom of each well screen (in m), below its top.")
        ("bc.well_rate", opt::list()->default_value("[]"),
         "Fluid volume extracted by each well at schedule ratio 1, spread uniformly over the "
         "screen: m^3/s in 3D, m^2/s per unit out-of-plane length in 2D. Negative injects.")
        ("bc.num_well_schedule", opt::value<int>(&p.bc.num_well_schedule)->default_value(1),
         "Number of points in the well rate schedule.")
        ("bc.well_schedule_time_in_yr", opt::list()->default_value("[0]"),
         "Non-decreasing well schedule times (in yr). A repeated time is a step.")
        ("bc.well_schedule_ratio", opt::list()->default_value("[1]"),
         "Rate multiplier at each schedule time: linear in between, held outside. Each step "
         "applies its exact time average, so the volume pumped is independent of dt.")

         // General stress (Neumann) boundary conditions 
        ("bc.stress_bc_x0", opt::value<int>(&p.bc.stress_bc_x0)->default_value(0),
         "Stress boundary type. 0: no added traction; 1: x normal-stress component; 2: y "
         "normal-stress component (3D only); 3: z normal-stress component; 4: facet-normal "
         "stress on any orientation. Use the corresponding bc.stress_val_* in Pa; tension is "
         "positive and compression negative.")
        ("bc.stress_bc_x1", opt::value<int>(&p.bc.stress_bc_x1)->default_value(0),
         "Stress boundary type. 0: no added traction; 1: x normal-stress component; 2: y "
         "normal-stress component (3D only); 3: z normal-stress component; 4: facet-normal "
         "stress on any orientation. Use the corresponding bc.stress_val_* in Pa; tension is "
         "positive and compression negative.")
        ("bc.stress_bc_y0", opt::value<int>(&p.bc.stress_bc_y0)->default_value(0),
         "Stress boundary type. 0: no added traction; 1: x normal-stress component; 2: y "
         "normal-stress component (3D only); 3: z normal-stress component; 4: facet-normal "
         "stress on any orientation. Use the corresponding bc.stress_val_* in Pa; tension is "
         "positive and compression negative.")
        ("bc.stress_bc_y1", opt::value<int>(&p.bc.stress_bc_y1)->default_value(0),
         "Stress boundary type. 0: no added traction; 1: x normal-stress component; 2: y "
         "normal-stress component (3D only); 3: z normal-stress component; 4: facet-normal "
         "stress on any orientation. Use the corresponding bc.stress_val_* in Pa; tension is "
         "positive and compression negative.")
        ("bc.stress_bc_z0", opt::value<int>(&p.bc.stress_bc_z0)->default_value(0),
         "Stress boundary type. 0: no added traction; 1: x normal-stress component; 2: y "
         "normal-stress component (3D only); 3: z normal-stress component; 4: facet-normal "
         "stress on any orientation. Use the corresponding bc.stress_val_* in Pa; tension is "
         "positive and compression negative.")
        ("bc.stress_bc_z1", opt::value<int>(&p.bc.stress_bc_z1)->default_value(0),
         "Stress boundary type. 0: no added traction; 1: x normal-stress component; 2: y "
         "normal-stress component (3D only); 3: z normal-stress component; 4: facet-normal "
         "stress on any orientation. Use the corresponding bc.stress_val_* in Pa; tension is "
         "positive and compression negative.")

        ("bc.stress_val_x0", opt::value<double>(&p.bc.stress_val_x0)->default_value(0.0),
         "Prescribed stress value (Pa): positive in tension, negative in compression; the "
         "corresponding bc.stress_bc_* selects the component or facet-normal direction.")
        ("bc.stress_val_x1", opt::value<double>(&p.bc.stress_val_x1)->default_value(0.0),
         "Prescribed stress value (Pa): positive in tension, negative in compression; the "
         "corresponding bc.stress_bc_* selects the component or facet-normal direction.")
        ("bc.stress_val_y0", opt::value<double>(&p.bc.stress_val_y0)->default_value(0.0),
         "Prescribed stress value (Pa): positive in tension, negative in compression; the "
         "corresponding bc.stress_bc_* selects the component or facet-normal direction.")
        ("bc.stress_val_y1", opt::value<double>(&p.bc.stress_val_y1)->default_value(0.0),
         "Prescribed stress value (Pa): positive in tension, negative in compression; the "
         "corresponding bc.stress_bc_* selects the component or facet-normal direction.")
        ("bc.stress_val_z0", opt::value<double>(&p.bc.stress_val_z0)->default_value(0.0),
         "Prescribed stress value (Pa): positive in tension, negative in compression; the "
         "corresponding bc.stress_bc_* selects the component or facet-normal direction.")
        ("bc.stress_val_z1", opt::value<double>(&p.bc.stress_val_z1)->default_value(0.0),
         "Prescribed stress value (Pa): positive in tension, negative in compression; the "
         "corresponding bc.stress_bc_* selects the component or facet-normal direction.")

        ("bc.vbc_x0", opt::value<int>(&p.bc.vbc_x0)->default_value(1),
         "Type of velocity boundary condition for the left/western side. "
         "Odd number indicates the normal component of the velocity is fixed. "
         "Possible type is \n"
         "0: all velocity components free;\n"
         "1: normal component fixed, shear components free;\n"
         "2: normal component free, shear components fixed at 0;\n"
         "3: normal component fixed, shear components fixed at 0;\n"
         "4: normal component free, shear component (not z) fixed, z component fixed at 0, only in 3D;\n"
         "5: normal component fixed at 0, shear component (not z) fixed, z component fixed at 0, only in 3D;\n"
         "6: normal component fixed (bc.vbc_val_*), shear component (not z) fixed (bc.vbc_val_*_l), z component free;\n"
         "7: normal component fixed, shear component (not z) fixed at 0, z component free, only in 3D;\n"
         "11: horizontal normal component fixed, z component and horizontal shear component free, only in 3D;\n"
         "13: horizontal normal component fixed, z component and horizontal shear component fixed at 0, only in 3D;\n")
        ("bc.vbc_x1", opt::value<int>(&p.bc.vbc_x1)->default_value(1),
         "Velocity boundary type; see bc.vbc_x0 for the type definitions and dimensional "
         "restrictions. The prescribed value is bc.vbc_val_x1.")
        ("bc.vbc_val_x0", opt::value<double>(&p.bc.vbc_val_x0)->default_value(-1e-9),
         "Value of boundary condition for left/western side (if velocity, unit is m/s; if stress, unit is Pa)")
        ("bc.vbc_val_x0_l", opt::value<double>(&p.bc.vbc_val_x0_l)->default_value(0),
         "Lateral shear velocity for left/western side (used with bc type 6, unit is m/s)")
        ("bc.vbc_val_division_x0_min", opt::value<double>(&p.bc.vbc_val_division_x0_min)->default_value(1.),
         "Ratio of first rate change division from the top to the right boundary width (for option = 8)")
        ("bc.vbc_val_division_x0_max", opt::value<double>(&p.bc.vbc_val_division_x0_max)->default_value(1.),
         "Ratio of second rate change division from the top to the right boundary width (for option = 8)")
        ("bc.vbc_val_x0_ratio0", opt::value<double>(&p.bc.vbc_val_x0_ratio0)->default_value(1.),
         "Ratio of the velocity at the surface of the left boundary width (for option = 1 and 3)")
        ("bc.vbc_val_x0_ratio1", opt::value<double>(&p.bc.vbc_val_x0_ratio1)->default_value(1.),
         "Ratio of the velocity at the vbc_val_division_x0_min of the left boundary width (for option = 1 and 3)")
        ("bc.vbc_val_x0_ratio2", opt::value<double>(&p.bc.vbc_val_x0_ratio2)->default_value(1.),
         "Ratio of the velocity at the vbc_val_division_x1_max of the left boundary width (for option = 1 and 3)")
        ("bc.vbc_val_x0_ratio3", opt::value<double>(&p.bc.vbc_val_x0_ratio3)->default_value(1.),
         "Ratio of the velocity at the bottom of the left boundary width (for option = 1 and 3)")

        ("bc.num_vbc_period_x0", opt::value<int>(&p.bc.num_vbc_period_x0)->default_value(1),
         "Number of points in the x0 velocity period table")
        ("bc.num_vbc_period_x1", opt::value<int>(&p.bc.num_vbc_period_x1)->default_value(1),
         "Number of points in the x1 velocity period table")
        ("bc.vbc_period_x0_time_in_yr", opt::list()->default_value("[1.e6]"),
         "Times of the x0 velocity period table '[t0, t1, ...]' (num_vbc_period_x0 of them, in yr, "
         "non-decreasing): the ratio is linear between them and constant outside; repeat a time "
         "for a step")
        ("bc.vbc_period_x1_time_in_yr", opt::list()->default_value("[1.e6]"),
         "Times of the x1 velocity period table (in yr), as for x0")
        ("bc.vbc_period_x0_ratio", opt::list()->default_value("[1.]"),
         "Ratio of vbc_val_x0 at each time of the x0 period table '[r0, r1, r2, ...]'. It scales "
         "the normal velocity only: tangential values (types 4 and 5, vbc_val_x0_l) stay constant")
        ("bc.vbc_period_x1_ratio", opt::list()->default_value("[1.]"),
         "Ratio of vbc_val_x1 at each time of the x1 period table, as for x0")

        ("bc.vbc_val_x1", opt::value<double>(&p.bc.vbc_val_x1)->default_value(1e-9),
         "Value of boundary condition for the right/eastern side (if velocity, unit is m/s; if stress, unit is Pa)")
        ("bc.vbc_val_x1_l", opt::value<double>(&p.bc.vbc_val_x1_l)->default_value(0),
         "Lateral shear velocity for right/eastern side (used with bc type 6, unit is m/s)")
        ("bc.vbc_val_division_x1_min", opt::value<double>(&p.bc.vbc_val_division_x1_min)->default_value(1.),
         "Ratio of first rate change division from the top to the right boundary width (for option = 8)")
        ("bc.vbc_val_division_x1_max", opt::value<double>(&p.bc.vbc_val_division_x1_max)->default_value(1.),
         "Ratio of second rate change division from the top to the right boundary width (for option = 8)")
        ("bc.vbc_val_x1_ratio0", opt::value<double>(&p.bc.vbc_val_x1_ratio0)->default_value(1.),
         "Ratio of the velocity at the surface of the right boundary width (for option = 1 and 3)")
        ("bc.vbc_val_x1_ratio1", opt::value<double>(&p.bc.vbc_val_x1_ratio1)->default_value(1.),
         "Ratio of the velocity at the vbc_val_division_x0_min of the right boundary width (for option = 1 and 3)")
        ("bc.vbc_val_x1_ratio2", opt::value<double>(&p.bc.vbc_val_x1_ratio2)->default_value(1.),
         "Ratio of the velocity at the vbc_val_division_x1_max of the right boundary width (for option = 1 and 3)")
        ("bc.vbc_val_x1_ratio3", opt::value<double>(&p.bc.vbc_val_x1_ratio3)->default_value(1.),
         "Ratio of the velocity at the bottom of the right boundary width (for option = 1 and 3)")

        ("bc.vbc_y0", opt::value<int>(&p.bc.vbc_y0)->default_value(0),
         "Velocity boundary type; see bc.vbc_x0 for the type definitions and dimensional "
         "restrictions. The prescribed value is bc.vbc_val_y0.")
        ("bc.vbc_y1", opt::value<int>(&p.bc.vbc_y1)->default_value(0),
         "Velocity boundary type; see bc.vbc_x0 for the type definitions and dimensional "
         "restrictions. The prescribed value is bc.vbc_val_y1.")
        ("bc.vbc_val_y0", opt::value<double>(&p.bc.vbc_val_y0)->default_value(0),
         "Value of boundary condition for the southern side (if velocity, unit is m/s; if stress, unit is Pa)")
        ("bc.vbc_val_y0_l", opt::value<double>(&p.bc.vbc_val_y0_l)->default_value(0),
         "Lateral shear velocity for southern side (used with bc type 6, unit is m/s)")
        ("bc.vbc_val_y1", opt::value<double>(&p.bc.vbc_val_y1)->default_value(0),
         "Value of boundary condition for the northern side (if velocity, unit is m/s; if stress, unit is Pa)")
        ("bc.vbc_val_y1_l", opt::value<double>(&p.bc.vbc_val_y1_l)->default_value(0),
         "Lateral shear velocity for northern side (used with bc type 6, unit is m/s)")
        ("bc.num_vbc_period_y0", opt::value<int>(&p.bc.num_vbc_period_y0)->default_value(1),
         "Number of points in the y0 velocity period table (3D only)")
        ("bc.num_vbc_period_y1", opt::value<int>(&p.bc.num_vbc_period_y1)->default_value(1),
         "Number of points in the y1 velocity period table (3D only)")
        ("bc.vbc_period_y0_time_in_yr", opt::list()->default_value("[1.e6]"),
         "Times of the y0 velocity period table (in yr, 3D only), as for x0")
        ("bc.vbc_period_y1_time_in_yr", opt::list()->default_value("[1.e6]"),
         "Times of the y1 velocity period table (in yr, 3D only), as for x0")
        ("bc.vbc_period_y0_ratio", opt::list()->default_value("[1.]"),
         "Ratio of vbc_val_y0 at each time of the y0 period table (3D only), as for x0")
        ("bc.vbc_period_y1_ratio", opt::list()->default_value("[1.]"),
         "Ratio of vbc_val_y1 at each time of the y1 period table (3D only), as for x0")

        ("bc.vbc_z0", opt::value<int>(&p.bc.vbc_z0)->default_value(0),
         "Bottom velocity boundary type. 0: all components free; 1: prescribed vertical velocity "
         "bc.vbc_val_z0, horizontal free; 2: vertical free, horizontal zero; 3: prescribed "
         "vertical velocity, horizontal zero; 4: prescribed x velocity and zero vertical "
         "velocity (2D only).")
        ("bc.vbc_z1", opt::value<int>(&p.bc.vbc_z1)->default_value(0),
         "Top velocity boundary type. 0: all components free; 1: prescribed vertical velocity "
         "bc.vbc_val_z1, horizontal free; 2: vertical free, horizontal zero; 3: prescribed "
         "vertical velocity, horizontal zero; 4: prescribed x velocity and zero vertical "
         "velocity (2D only); 8: rigid frictionless platen with one shared free vertical "
         "velocity; load it with bc.stress_bc_z1.")
        ("bc.vbc_val_z0", opt::value<double>(&p.bc.vbc_val_z0)->default_value(0),
         "Value of boundary condition for the bottom side (if velocity, unit is m/s; if stress, unit is Pa)")
        ("bc.vbc_val_z1", opt::value<double>(&p.bc.vbc_val_z1)->default_value(0),
         "Value of boundary condition for the top side (if velocity, unit is m/s; if stress, unit is Pa)")
        ("bc.vbc_val_z1_loading_period", opt::value<double>(&p.bc.vbc_val_z1_loading_period)->default_value(std::numeric_limits<double>::max()),
         "Loading period for the velocity on the top side")

        ("bc.vbc_n0", opt::value<int>(&p.bc.vbc_n0)->default_value(1),
         "Type of boundary condition for slant boundary #0 (only type 1, 3, 11, 13 are supported).")
        ("bc.vbc_val_n0", opt::value<double>(&p.bc.vbc_val_n0)->default_value(0),
         "Value of boundary condition for slant boundary #0 (if velocity, unit is m/s, "
         "outward normal direction is positive; if stress, unit is Pa)")
        ("bc.vbc_n1", opt::value<int>(&p.bc.vbc_n1)->default_value(1),
         "Type of boundary condition for slant boundary #1 (only type 1, 3, 11, 13 are supported).")
        ("bc.vbc_val_n1", opt::value<double>(&p.bc.vbc_val_n1)->default_value(0),
         "Value of boundary condition for slant boundary #1 (if velocity, unit is m/s, "
         "outward normal direction is positive; if stress, unit is Pa)")
        ("bc.vbc_n2", opt::value<int>(&p.bc.vbc_n2)->default_value(1),
         "Type of boundary condition for slant boundary #2 (only type 1, 3, 11, 13 are supported).")
        ("bc.vbc_val_n2", opt::value<double>(&p.bc.vbc_val_n2)->default_value(0),
         "Value of boundary condition for slant boundary #2 (if velocity, unit is m/s, "
         "outward normal direction is positive; if stress, unit is Pa)")
        ("bc.vbc_n3", opt::value<int>(&p.bc.vbc_n3)->default_value(1),
         "Type of boundary condition for slant boundary #3 (only type 1, 3, 11, 13 are supported).")
        ("bc.vbc_val_n3", opt::value<double>(&p.bc.vbc_val_n3)->default_value(0),
         "Value of boundary condition for slant boundary #3 (if velocity, unit is m/s, "
         "outward normal direction is positive; if stress, unit is Pa)")
        
        ("bc.bottom_shear_zone_thickness", opt::value<double>(&p.bc.bottom_shear_zone_thickness)->default_value(0.),
         "Thickness of the bottom shear zone for bc_x0 = 3 and ndims = 2, in m: the x0 wall velocity "
         "tapers linearly to 0 from this height above the bottom down to the bottom. 0 disables it.")
        ;

    cfg.add_options()
        ("ic.mattype_option", opt::value<int>(&p.ic.mattype_option)->default_value(0),
         "How to set the initial material type of markers?\n"
         "0: marker's mattype is determined by regional attribute (currently can be modified only in the code or in a poly file).\n"
         "1: marker's mattype is layered.\n"
         "101: custom mattype.")
        ("ic.num_mattype_layers", opt::value<int>(&p.ic.num_mattype_layers)->default_value(2),
         "Number of material layers")
        ("ic.layer_mattypes", opt::list()->default_value("[0,1]"),
         "Material type of each material layer '[d0, d1, d2, ...]', d0<=d1<=d2...")
        ("ic.mattype_layer_depths", opt::list()->default_value("[0.5]"),
         "Depths of the interfaces of each material layer '[d0, d1, d2, ...]', d0<=d1<=d2..., (in unit of zlength)")

        ("ic.weakzone_option", opt::value<int>(&p.ic.weakzone_option)->default_value(1),
         "How to set the initial weak zone?\n"
         "0: no weak zone.\n"
         "1: planar weak zone with specified azimuth, inclination, halfwidth, min/max depth range, and center location.\n"
         "2: ellipsoidal weak zone with specified center location and semi-axes.\n"
         "3: Gaussian distribution point weak zone with specified center location and standard deviation.\n"
         "4: multiple planar segments; each segment uses the same parameters as option 1 but given as '[v0, v1, ...]' arrays.\n")
        ("ic.is_restarting_weakzone", opt::value<bool>(&p.ic.is_restarting_weakzone)->default_value(false),
         "Create a new weakzone when restarting?")
        ("ic.weakzone_plstrain", opt::value<double>(&p.ic.weakzone_plstrain)->default_value(0.1),
         "Initial plastic strain in the weak zone")
        ("ic.weakzone_azimuth", opt::value<double>(&p.ic.weakzone_azimuth)->default_value(0),
         "Azimuth angle (relative to +y axis) (in degree)")
        ("ic.weakzone_inclination", opt::value<double>(&p.ic.weakzone_inclination)->default_value(90),
         "Inclination angle (relative to horizontal plane) (in degree)")
        ("ic.weakzone_halfwidth", opt::value<double>(&p.ic.weakzone_halfwidth)->default_value(1.5),
         "Half-width of the weak zone (in unit of mesh.resolution)")
        ("ic.weakzone_y_min", opt::value<double>(&p.ic.weakzone_y_min)->default_value(0),
         "Y-direction range (between 0 and 1, in unit of mesh.ylength)")
        ("ic.weakzone_y_max", opt::value<double>(&p.ic.weakzone_y_max)->default_value(1),
         "Y-direction range (between 0 and 1, in unit of mesh.ylength)")
        ("ic.weakzone_depth_min", opt::value<double>(&p.ic.weakzone_depth_min)->default_value(0),
         "Depth upper (shallower) range (between 0 and 1, in unit of mesh.zlength)")
        ("ic.weakzone_depth_max", opt::value<double>(&p.ic.weakzone_depth_max)->default_value(1),
         "Depth lower (deeper) range (between 0 and 1, in unit of mesh.zlength)")
        ("ic.weakzone_xcenter", opt::value<double>(&p.ic.weakzone_xcenter)->default_value(0.5),
         "Location of weak zone center in x direction (between 0 and 1, in unit of mesh.xlength)")
        ("ic.weakzone_ycenter", opt::value<double>(&p.ic.weakzone_ycenter)->default_value(0.5),
         "Location of weak zone center in y direction (between 0 and 1, in unit of mesh.ylength)")
        ("ic.weakzone_zcenter", opt::value<double>(&p.ic.weakzone_zcenter)->default_value(0.5),
         "Location of weak zone center in z direction (between 0 and 1, in unit of mesh.zlength)")
        ("ic.weakzone_xsemi_axis", opt::value<double>(&p.ic.weakzone_xsemi_axis)->default_value(1e3),
         "Length of weak zone semi-axis in x direction (in meters)")
        ("ic.weakzone_ysemi_axis", opt::value<double>(&p.ic.weakzone_ysemi_axis)->default_value(1e3),
         "Length of weak zone semi-axis in y direction (in meters)")
        ("ic.weakzone_zsemi_axis", opt::value<double>(&p.ic.weakzone_zsemi_axis)->default_value(1e3),
         "Length of weak zone semi-axis in z direction (in meters)\n")
        ("ic.weakzone_standard_deviation", opt::value<double>(&p.ic.weakzone_standard_deviation)->default_value(1e3),
         "Standard deviation of Gaussian weak zone (in meters).\n"
         "Used as point-Gaussian std dev (option 3) and as along-strike sigma (option 4).\n")
        ("ic.weakzone_gaussian_amplitude", opt::value<double>(&p.ic.weakzone_gaussian_amplitude)->default_value(0),
         "Amplitude of the Gaussian x-shift for weakzone_option=4 (in meters).\n"
         "The fault x-position is offset by A*exp(-(y-y0)^2/(2*sigma^2)) where\n"
         "A = weakzone_gaussian_amplitude and sigma = weakzone_standard_deviation.\n")

        // multi-segment planar weak zone (weakzone_option == 5)
        ("ic.weakzone_num_segments", opt::value<int>(&p.ic.weakzone_num_segments)->default_value(1),
         "Number of planar segments for multi-segment weak zone (weakzone_option=5)")
        ("ic.weakzone_segments_xcenter", opt::list()->default_value("[0.5]"),
         "X centers of each segment (between 0 and 1, in unit of mesh.xlength), '[v0, v1, ...]'")
        ("ic.weakzone_segments_ycenter", opt::list()->default_value("[0.5]"),
         "Y centers of each segment (between 0 and 1, in unit of mesh.ylength), '[v0, v1, ...]' (3D only)")
        ("ic.weakzone_segments_zcenter", opt::list()->default_value("[0.5]"),
         "Z centers of each segment (between 0 and 1, in unit of mesh.zlength), '[v0, v1, ...]'")
        ("ic.weakzone_segments_azimuth", opt::list()->default_value("[0]"),
         "Azimuth angles for each segment (in degrees), '[v0, v1, ...]'")
        ("ic.weakzone_segments_inclination", opt::list()->default_value("[90]"),
         "Inclination angles for each segment (in degrees), '[v0, v1, ...]'")
        ("ic.weakzone_segments_halfwidth", opt::list()->default_value("[1.5]"),
         "Half-widths for each segment (true perpendicular distance, in unit of mesh.resolution), '[v0, v1, ...]'")
        ("ic.weakzone_segments_x_min", opt::list()->default_value("[0]"),
         "X lower bound for each segment (between 0 and 1, in unit of mesh.xlength), '[v0, v1, ...]'")
        ("ic.weakzone_segments_x_max", opt::list()->default_value("[1]"),
         "X upper bound for each segment (between 0 and 1, in unit of mesh.xlength), '[v0, v1, ...]'")
        ("ic.weakzone_segments_depth_min", opt::list()->default_value("[0]"),
         "Depth upper (shallower) bound for each segment (between 0 and 1, in unit of mesh.zlength), '[v0, v1, ...]'")
        ("ic.weakzone_segments_depth_max", opt::list()->default_value("[1]"),
         "Depth lower (deeper) bound for each segment (between 0 and 1, in unit of mesh.zlength), '[v0, v1, ...]'")
        ("ic.weakzone_segments_y_min", opt::list()->default_value("[0]"),
         "Y lower bound for each segment (between 0 and 1, in unit of mesh.ylength), '[v0, v1, ...]' (3D only)")
        ("ic.weakzone_segments_y_max", opt::list()->default_value("[1]"),
         "Y upper bound for each segment (between 0 and 1, in unit of mesh.ylength), '[v0, v1, ...]' (3D only)\n")

        ("ic.temperature_option", opt::value<int>(&p.ic.temperature_option)->default_value(0),
         "How to set the initial temperature?\n"
         "0: uniform half-space cooling.\n"
         "1: continental thermal gradient.\n"
         "2: equilibrium temperature for crustal radiogenic heating.\n"
         "3: continental thermal gradient with adiabatic mantle.\n"
         "4: uniform ic.temp_background plus a Gaussian anomaly.\n"
         "5: linear conduction from bc.surface_temperature at the top to bc.mantle_temperature\n"
         "   at the bottom (mesh.zlength), plus the same Gaussian anomaly.\n"
         "90: temperature read from an external grid\n")

        // for temperature_option = 0
	    ("ic.oceanic_plate_age_in_yr", opt::value<double>(&p.ic.oceanic_plate_age_in_yr)->default_value(60e6),
         "Age of the oceanic plate (in years), used for the temperature profile on the plate.\n")
        // for temperature_option = 2
        ("ic.num_radiogenic_heat_layer", opt::value<int>(&p.ic.nhlayer)->default_value(1),
         "Number of radiogenic heat layers")
        ("ic.radiogenic_heat_boundry", opt::list()->default_value("[-1,-1]"),
        "The layer boundry of radiogenic heating from the top to the bottom (depth -> positive, in meters).\n"
        "Using -1 to indicate the top and bottom of the model.\n")
        ("ic.radiogenic_heat_mat_in_layer",  opt::list()->default_value("[0]"),
        "The material type of radiogenic heating in each layer (one value: the same in every layer).\n")
        // for temperature_option = 1 and 2
        ("ic.radiogenic_folding_depth", opt::value<double>(&p.ic.radiogenic_folding_depth)->default_value(10.e3),
         "Radiogenic folding depth (in meters), used for the temperature profile on the plate.\n")
        // for temperature_option = 1
        ("ic.radiogenic_heating_of_crust", opt::value<double>(&p.ic.radiogenic_heating_of_crust)->default_value(1.e-9),
         "Radiogenic heating of crust (in W/kg), used for the temperature profile on the plate.\n")
        ("ic.continental_plate_age_in_yr", opt::value<double>(&p.ic.continental_plate_age_in_yr)->default_value(100.e6),
         "Age of the continental plate (in years), used for the temperature profile on the plate.\n")
        ("ic.radiogenic_crustal_thickness", opt::value<double>(&p.ic.radiogenic_crustal_thickness)->default_value(30.e3),
         "Radiogenic crustal thickness (in meters), used for the temperature profile on the plate.\n")
        ("ic.lithospheric_thickness", opt::value<double>(&p.ic.lithospheric_thickness)->default_value(100.e3),
         "Lithospheric thickness (in meters), used for the temperature profile on the plate.\n")
         // for temperature_option = 3
        ("ic.radiogenic_heat_dome_center_x", opt::value<double>(&p.ic.radiogenic_heat_dome_center_x)->default_value(0.5),
         "x of center of geothermal dome, used for the temperature profile on the plate.\n")
        ("ic.radiogenic_heat_dome_center_y", opt::value<double>(&p.ic.radiogenic_heat_dome_center_y)->default_value(0.5),
         "y of center of geothermal dome, used for the temperature profile on the plate.\n")
        ("ic.surface_heat_flux", opt::value<double>(&p.ic.surface_heat_flux)->default_value(70e-3),
         "surface heat flux, used for the temperature profile on the plate.\n")
        ("ic.radiogenic_heat_dome_amplitude", opt::value<double>(&p.ic.radiogenic_heat_dome_amplitude)->default_value(0),
         "Dome amplitude, used for the temperature profile on the plate.\n")
        ("ic.radiogenic_heat_dome_width", opt::value<double>(&p.ic.radiogenic_heat_dome_width)->default_value(0),
         "Dome half width, used for the temperature profile on the plate.\n")
        ("ic.radiogenic_heat_dome_width_y", opt::value<double>(&p.ic.radiogenic_heat_dome_width_y)->default_value(0),
         "Dome half width in y-axis, used for the temperature profile on the plate (3D).\n"
         " 0: Isotropic (uses radiogenic_heat_dome_width for both axes).\n"
         ">0: Specific width for y-axis (elliptical dome).\n"
         "<0: Infinite width; thermal profile is uniform along y-axis (ridge/hinge).\n")
        ("ic.radiogenic_heat_dome_azimuth", opt::value<double>(&p.ic.radiogenic_heat_dome_azimuth)->default_value(0),
         "Dome azimuth (rotation angle), used for the temperature profile on the plate (in degree) (3D).\n")

        ("ic.temp_background", opt::value<double>(&p.ic.temp_background)->default_value(273),
         "Uniform background temperature in K for temperature_option=4.")
        ("ic.temp_anomaly_amplitude", opt::value<double>(&p.ic.temp_anomaly_amplitude)->default_value(0),
         "Gaussian temperature amplitude in K.")
        ("ic.temp_anomaly_xcenter", opt::value<double>(&p.ic.temp_anomaly_xcenter)->default_value(0),
         "Gaussian thermal anomaly x coordinate in m.")
        ("ic.temp_anomaly_ycenter", opt::value<double>(&p.ic.temp_anomaly_ycenter)->default_value(0),
         "Gaussian thermal anomaly y coordinate in m.")
        ("ic.temp_anomaly_zcenter", opt::value<double>(&p.ic.temp_anomaly_zcenter)->default_value(0),
         "Gaussian thermal anomaly z coordinate in m.")
        ("ic.temp_anomaly_xsemi_axis", opt::value<double>(&p.ic.temp_anomaly_xsemi_axis)->default_value(1),
         "Positive Gaussian x scale in m.")
        ("ic.temp_anomaly_ysemi_axis", opt::value<double>(&p.ic.temp_anomaly_ysemi_axis)->default_value(1),
         "Positive Gaussian y scale in m.")
        ("ic.temp_anomaly_zsemi_axis", opt::value<double>(&p.ic.temp_anomaly_zsemi_axis)->default_value(1),
         "Positive Gaussian z scale in m.")
        // for temperature_option = 90
        ("ic.Temp_filename", opt::value<std::string>(&p.ic.Temp_filename)->default_value("Thermal.dat"),
         "Filename of the input thermal field, for temperature_option=90.\n"
         "The file format is: <x-coord y-coord (if 3D) z-coord temperature>.\n"
	    "Can be generated by COMSOL.")
	    ("ic.Nodes_filename", opt::value<std::string>(&p.ic.Nodes_filename)->default_value("Coord.dat"),
         "Filename of the mesh coordinates, for temperature_option=90.\n"
         "The nodes may be ordered differently with the thermal file. The file format is: "
         "<x-coord y-coord (if 3D) z-coord>.\n"
         "Can be generated by COMSOL.")
	    ("ic.Connectivity_filename", opt::value<std::string>(&p.ic.Connectivity_filename)->default_value("Connectivity.dat"),
         "Filename of the mesh connectivity, for temperature_option = 90 only.\n"
         "The file format is: <#ofnode0 #ofnode1 #ofnode2 #ofnode3 (if 3D)>.\n"
         "Can be generated by COMSOL.\n")

        ("ic.isostasy_adjustment_time_in_yr", opt::value<double>(&p.ic.isostasy_adjustment_time_in_yr)->default_value(0),
         "Time for spinning up isostasy adjustment.\n")
        
        ("ic.excess_pore_pressure", opt::value<double>(&p.ic.excess_pore_pressure)->default_value(0.0),
         "Initial excess_pore_pressure except for boundary.\n")
        ("ic.pore_pressure_option", opt::value<int>(&p.ic.pore_pressure_option)->default_value(0),
         "Initial pressure: 0 hydrostatic/excess, 1 linear x profile between prescribed pressure faces, 2 uniform.")
        ("ic.pore_pressure_value", opt::value<double>(&p.ic.pore_pressure_value)->default_value(0),
         "Uniform initial pressure in Pa for pore_pressure_option=2.")
        ("ic.water_table_z", opt::value<double>(&p.ic.water_table_z)->default_value(0.0),
         "Height of the initial hydrostatic water table (m): p = -rho_f*g*(z - water_table_z), "
         "negative (suction) above it. The default 0 is the original surface.\n")
         ("ic.has_body_force_adjustment", opt::value<bool>(&p.ic.has_body_force_adjustment)->default_value(false),
         "Conducting PT loop to get initial stress field from inital guess")
        ("ic.initial_stress_option", opt::value<int>(&p.ic.initial_stress_option)->default_value(0),
         "How to initialize stress?\n"
         "0: use the legacy gravity-dependent initialization.\n"
         "1: prescribe a homogeneous absolute Cauchy stress tensor; requires gravity=0.\n")
#ifdef THREED
        ("ic.initial_stress", opt::list()->default_value("[0,0,0,0,0,0]"),
         "Homogeneous absolute Cauchy stress [sxx,syy,szz,sxy,sxz,syz] in Pa; "
         "compression is negative.")
#else
        ("ic.initial_stress", opt::list()->default_value("[0,0,0]"),
         "Homogeneous absolute Cauchy stress [sxx,szz,sxz] in Pa; compression is "
         "negative. Plane strain initializes syy=(sxx+szz)/2.")
#endif

        ;

    cfg.add_options()
        ("mat.rheology_type", opt::value<std::string>()->required(),
         "Type of rheology, either 'elastic', 'viscous' (experimental), 'maxwell', "
         "'elasto-plastic', 'elasto-visco-plastic', "
         "'elasto-plastic-rate-state-friction', or 'elasto-visco-plastic-rate-state-friction' "
         "(short '-rsf' aliases are also accepted)")
        ("mat.is_plane_strain", opt::value<bool>(&p.mat.is_plane_strain)->default_value(false),
         "Use plane-strain constitutive updates in 2D (elastic and elasto-plastic paths).\n")

        ("mat.phase_change_option", opt::value<int>(&p.mat.phase_change_option)->default_value(0),
         "What kind of phase changes?\n"
         "0: no phase changes.\n"
         "1: simple rules of subduction-related phase changes. See simple_subduction() in phasechanges.cxx for more details.\n"
         "3: subduction phase-change rule engine (see subduction_phase_changes() in phasechanges.cxx). A marker "
         "converts to a metamorphic product when it crosses a pressure-temperature reaction boundary (pressure: "
         "see mat.phase_change_pressure_option; depth: below the local "
         "surface; chlorite-out and the first two sediment steps read depth only); each transition "
         "activates only when its target mat.mattype_* is set to a distinct material and its switch is on "
         "(mat.has_eclogitization and mat.has_serpentinization default on, mat.has_hydrous_mantle_dehydration "
         "and mat.has_sediment_metamorphism off). Transitions: oceanic crust->"
         "eclogite (Hacker 2003), mantle<->serpentinite (Ulmer & Trommsdorff 1995; serpentinization also needs "
         "wet crust -- oceanic crust, sediment or sedimentary rock -- within mat.serpentinization_thickness below "
         "the marker, and no deeper than mat.serpentinization_max_depth), "
         "hydrous mantle->mantle (Grove 2009), and "
         "sediment->sedimentary rock->metasediment->schist (compaction uses geoflac's unsourced thresholds, "
         "the next step Montanez & Crossey 2017, the last Wei 2004).\n"
         "101: custom phase changes (custom_phase_change() in phasechanges.cxx).")
        ("mat.phase_change_pressure_option", opt::value<int>(&p.mat.phase_change_pressure_option)->default_value(0),
         "The pressure the mat.phase_change_option 3 rules read (options 1 and 101 always use ref_pressure()).\n"
         "0: the mean stress of the marker's element, compression positive, in 2D the in-plane mean. Tectonic "
         "over- and underpressure move the reaction boundaries.\n"
         "1: the reference pressure, ref_pressure() (control.ref_pressure_option), at the marker's depth below "
         "the local surface. The phase diagrams are drawn for a lithostat, and geoflac reads one.")
        ("mat.num_materials", opt::value<int>(&p.mat.nmat)->default_value(1),
         "Number of material types")
        ("mat.mattype_ref", opt::value<int>(&p.mat.mattype_ref)->default_value(0),
         "Index of reference material. For compute_dt(), compute_mass() and ref_pressure()")
        ("mat.mattype_mantle", opt::value<int>(&p.mat.mattype_mantle)->default_value(0),
         "Index of mantle material. For continental thermal gradient. With mat.phase_change_option 3: what "
         "serpentinizes above the slab and what serpentinite and hydrous mantle dehydrate to; no rule may "
         "start from it, and it must be set in the cfg when option 3 reads it.")
        ("mat.mattype_depleted_mantle", opt::value<int>(&p.mat.mattype_depleted_mantle)->default_value(0),
         "Index of depleted mantle material. Unused: no phase change reads it.")
        ("mat.mattype_partial_melting_mantle", opt::value<int>(&p.mat.mattype_partial_melting_mantle)->default_value(0),
         "Index of partial melting mantle material. Unused: no phase change reads it.")
        ("mat.mattype_continental_crust", opt::value<int>(&p.mat.mattype_continental_crust)->default_value(-1),
         "Index of continental crust. Read only by the continental geotherm (ic.temperature_option 1), "
         "which needs it in [0, mat.num_materials). Default -1 (unset).")
        ("mat.mattype_oceanic_crust", opt::value<int>(&p.mat.mattype_oceanic_crust)->default_value(-1),
         "Index of oceanic crust. Read only by mat.phase_change_option 3: the source of the eclogite "
         "rule (with mat.has_eclogitization) and wet crust for the serpentinization grid (with "
         "mat.has_serpentinization); either needs it in [0, mat.num_materials). Default -1 (unset).")
        ("mat.mattype_serpentinite", opt::value<int>(&p.mat.mattype_serpentinite)->default_value(-1),
         "Index of serpentinite material, for mat.phase_change_option 3 with mat.has_serpentinization: what "
         "mantle above the slab turns into, and the source of serpentinite->mantle dehydration. Default -1 "
         "(unset); with serpentinization on it must be in [0, mat.num_materials) and distinct from "
         "mattype_mantle and the wet-crust materials (mattype_oceanic_crust, mattype_sediment, and "
         "mattype_sediment_rock when set).")
        ("mat.mattype_eclogite", opt::value<int>(&p.mat.mattype_eclogite)->default_value(-1),
         "Index of eclogite material. For the mat.phase_change_option 3 rule engine: oceanic-crust markers past the "
         "basalt->eclogite P-T boundary (Hacker 2003) convert to this. Default -1 (unset); must be set while "
         "mat.has_eclogitization is on.")
        ("mat.mattype_hydrous_mantle", opt::value<int>(&p.mat.mattype_hydrous_mantle)->default_value(-1),
         "Index of hydrous (chlorite-bearing) mantle material. For mat.phase_change_option 3: dehydrates to "
         "mattype_mantle above the chlorite-out boundary (Grove et al. 2009). Default -1 (unset); must be set while "
         "mat.has_hydrous_mantle_dehydration is on.")
        ("mat.mattype_sediment_rock", opt::value<int>(&p.mat.mattype_sediment_rock)->default_value(-1),
         "Index of consolidated sedimentary-rock material. For mat.phase_change_option 3: mattype_sediment compacts to "
         "this, which then metamorphoses to mattype_metasediment; when set, it is also wet crust for the "
         "serpentinization grid. Default -1 (unset); must be set while mat.has_sediment_metamorphism is on.")
        ("mat.mattype_metasediment", opt::value<int>(&p.mat.mattype_metasediment)->default_value(-1),
         "Index of metamorphic-sediment material. For mat.phase_change_option 3: mattype_sediment_rock -> this -> "
         "mattype_schist. Default -1 (unset); must be set while mat.has_sediment_metamorphism is on.")
        ("mat.mattype_schist", opt::value<int>(&p.mat.mattype_schist)->default_value(-1),
         "Index of schist material. For mat.phase_change_option 3: mattype_metasediment dehydrates to this "
         "(Wei et al. 2004). Default -1 (unset); must be set while mat.has_sediment_metamorphism is on.")
        ("mat.mattype_sediment", opt::value<int>(&p.mat.mattype_sediment)->default_value(0),
         "Index of sediment material. Should be the same as mattype['sediment'] in 2vtk.py for sediment time convertion. "
         "With mat.phase_change_option 3: wet crust for the serpentinization grid (with mat.has_serpentinization) "
         "and the source of the sediment chain (with mat.has_sediment_metamorphism); with either it must be "
         "set in the cfg.")
        ("mat.mattype_mor_extrusion", opt::value<int>(&p.mat.mattype_mor_extrusion)->default_value(0),
         "Index of mid-ocean-ridge extrusion material. Unused: no phase change reads it.")
        ("mat.mattype_asthenosphere", opt::value<int>(&p.mat.mattype_asthenosphere)->default_value(0),
         "Index of asthenosphere.")
        ("mat.max_viscosity", opt::value<double>(&p.mat.visc_max)->default_value(1e24),
         "Max. value of viscosity (in Pa.s)")
        ("mat.min_viscosity", opt::value<double>(&p.mat.visc_min)->default_value(1e18),
         "Min. value of viscosity (in Pa.s)")
        ("mat.max_tension", opt::value<double>(&p.mat.tension_max)->default_value(1e9),
         "Max. value of tensile stress (in Pa)")
        ("mat.max_thermal_diffusivity", opt::value<double>(&p.mat.therm_diff_max)->default_value(5e-6),
         "Max. value of thermal diffusivity (in m^2/s)")
        ("mat.convert_rate_oceanic_crust", opt::value<double>(&p.mat.convert_rate_oceanic_crust)->default_value(1.e-13),
         "Converting rate from melting to oceanic crust. Unused: no phase change reads it.")
        ("mat.serpentinization_thickness", opt::value<double>(&p.mat.serpentinization_thickness)->default_value(10e3),
         "Thickness (in m) of the mantle layer above wet crust (oceanic crust, sediment or sedimentary rock) that turns into serpentinite "
         "(mat.phase_change_option is 3 with mat.has_serpentinization). Counted in whole background-grid "
         "cells below the marker's own "
         "cell (cell size mesh.resolution * mesh.vt_size), down to a wet-crust-dominant cell "
         "(mat.wet_crust_fraction_threshold), so the effective "
         "reach differs from it by up to two cells, and the thickness must be positive and at least one cell. "
         "Default 10e3: heuristic, untuned.")
        ("mat.serpentinization_wet_crust_margin", opt::value<double>(&p.mat.serpentinization_wet_crust_margin)->default_value(-1.0),
         "Horizontal margin (in m) for the serpentinization wet-crust-below test. When >= 0, a wet-crust-dominant "
         "grid cell triggers conversion only if the marker's horizontal position lies within the cell's "
         "wet-crust extent (recorded per cell at grid build) widened by this margin -- suppressing the "
         "column-quantization false positive at a slab bend, where a cell is wet-crust-dominant from wet crust "
         "up to a whole cell size away horizontally. < 0 (default) disables: whole-column behavior.")
        ("mat.serpentinization_max_depth", opt::value<double>(&p.mat.serpentinization_max_depth)->default_value(150e3),
         "mat.phase_change_option 3 with mat.has_serpentinization: no mantle deeper than this (in m below the "
         "local free surface) serpentinizes. It also ends the serpentinization grid, at this depth plus "
         "mat.serpentinization_thickness and one cell. Must be positive. Default 150e3, heuristic, untuned: slab "
         "dehydration, the water source, is essentially complete by ~150 km, shallower than the antigorite "
         "stability field's cold-slab extrapolation (> 200 km at 500 C).")
        ("mat.wet_crust_fraction_threshold", opt::value<double>(&p.mat.wet_crust_fraction_threshold)->default_value(0.5),
         "mat.phase_change_option 3 with mat.has_serpentinization: a serpentinization grid cell is wet-crust-dominant, "
         "so it can serpentinize the mantle above it, when more than this fraction of its markers is wet crust. "
         "In [0, 1). A cell holds markers.markers_per_element times (cell / element measure) markers, nominally "
         "~2.3 in 2D and ~4.2 in 3D at mesh.vt_size 0.5, so at 0.5 a two-marker cell needs both. Default 0.5, "
         "untuned: geoflac's 0.8 (of an element, with no source) leaves out the cell that straddles the slab top "
         "and cuts the serpentinite band of examples/tectonics/subduction.cfg by more than half.")
        ("mat.has_serpentinization", opt::value<bool>(&p.mat.has_serpentinization)->default_value(true),
         "mat.phase_change_option 3: enable mantle<->serpentinite (mantle above the slab serpentinizes inside the "
         "antigorite stability field; serpentinite above 550 C outside it dehydrates back to mantle). Default on. "
         "Needs mattype_serpentinite and mattype_oceanic_crust, and mattype_sediment set in the cfg. "
         "Off: no serpentinization grid is built, and the checks on its inputs are skipped: mattype_mantle, "
         "mattype_serpentinite and the wet-crust materials distinct, mattype_serpentinite in range and set, "
         "mattype_sediment set, serpentinization_thickness positive and at least one grid cell, mesh.vt_size positive. "
         "mattype_mantle and mattype_sediment stay range-checked with the other base material "
         "indices, and the rule "
         "checks still run (see the other has_* switches); every key an active rule reads must be set in "
         "the cfg.")
        ("mat.has_eclogitization", opt::value<bool>(&p.mat.has_eclogitization)->default_value(true),
         "mat.phase_change_option 3: enable oceanic crust->eclogite (Hacker 2003). Default on. Needs "
         "mattype_oceanic_crust and mattype_eclogite set, distinct, and not mattype_mantle, and "
         "mattype_sediment different from mattype_oceanic_crust (deposits would turn into eclogite). "
         "Like every has_* switch, no two active rules may start from one material.")
        ("mat.has_hydrous_mantle_dehydration", opt::value<bool>(&p.mat.has_hydrous_mantle_dehydration)->default_value(false),
         "mat.phase_change_option 3: enable hydrous (chlorite) mantle->mantle dehydration (Grove 2009). "
         "Default off. Needs mattype_hydrous_mantle.")
        ("mat.has_sediment_metamorphism", opt::value<bool>(&p.mat.has_sediment_metamorphism)->default_value(false),
         "mat.phase_change_option 3: enable the sediment chain sediment->sedimentary rock->metasediment->schist "
         "(compaction uses geoflac's unsourced thresholds, the next step Montanez & Crossey 2017, the last "
         "Wei 2004). Default off. "
         "Needs mattype_sediment_rock / mattype_metasediment / mattype_schist.")

        // these parameters need to parsed later
        ("mat.rho0", opt::list()->default_value("[3210]"),
         "Density of the materials at 0 Pa and 273 K '[d0, d1, d2, ...]' (in kg/m^3)")
        ("mat.alpha", opt::list()->default_value("[3e-5]"),
         "Volumetic thermal expansion of the materials '[d0, d1, d2, ...]' (in 1/Kelvin)")

        ("mat.bulk_modulus", opt::list()->default_value("[128.2e9]"),
         "Bulk modulus of the materials '[d0, d1, d2, ...]' (in Pa)")
        ("mat.shear_modulus", opt::list()->default_value("[80.5e9]"),
         "Shear modulus of the materials '[d0, d1, d2, ...]' (in Pa)")

        ("mat.visc_exponent", opt::list()->default_value("[3.05]"),
         "Exponents of non-linear viscosity of the materials'[d0, d1, d2, ...]'")
        ("mat.visc_coefficient", opt::list()->default_value("[1.25e-1]"),
         "Pre-exponent coefficient of non-linear viscosity of the materials '[d0, d1, d2, ...]'")
        ("mat.visc_activation_energy", opt::list()->default_value("[3.76e5]"),
         "Activation energy of non-linear viscosity of the materials '[d0, d1, d2, ...]' (in J/mol)")
        ("mat.visc_activation_volume", opt::list()->default_value("[0.e-6]"),
         "Activation volume of non-linear viscosity of the materials '[d0, d1, d2, ...]' (in m^3/mol)")

        ("mat.heat_capacity", opt::list()->default_value("[1e3]"),
         "Heat capacity (isobaric) of the materials '[d0, d1, d2, ...]' (in J/kg/Kelvin)")
        ("mat.therm_cond", opt::list()->default_value("[3]"),
         "Thermal conductivity of the materials '[d0, d1, d2, ...]' (in W/m/Kelvin)")
        ("mat.radiogenic_heat_prod", opt::list()->default_value("[0]"),
         "Radiogenic heat production rate of the materials '[d0, d1, d2, ...]' (in W/kg)"
         "Only for initial temperature condition")

        ("mat.pls0", opt::list()->default_value("[0]"),
         "Plastic strain of the materials where weakening starts '[d0, d1, d2, ...]' (no unit)")
        ("mat.pls1", opt::list()->default_value("[0.1]"),
         "Plastic strain of the materials where weakening saturates '[d0, d1, d2, ...]' (no unit)")
        ("mat.cohesion0", opt::list()->default_value("[4e7]"),
         "Cohesion of the materials when weakening starts '[d0, d1, d2, ...]' (in Pa)")
        ("mat.cohesion1", opt::list()->default_value("[4e6]"),
         "Cohesion of the materials when weakening saturates '[d0, d1, d2, ...]' (in Pa)")
        ("mat.friction_angle0", opt::list()->default_value("[30]"),
         "Friction angle of the materials when weakening starts '[d0, d1, d2, ...]' (in degree)")
        ("mat.friction_angle1", opt::list()->default_value("[5]"),
         "Friction angle of the materials when weakening saturates '[d0, d1, d2, ...]' (in degree)")
        ("mat.dilation_angle0", opt::list()->default_value("[0]"),
         "Dilation angle of the materials when weakening starts '[d0, d1, d2, ...]' (in degree)")
        ("mat.dilation_angle1", opt::list()->default_value("[0]"),
         "Dilation angle of the materials when weakening saturates '[d0, d1, d2, ...]' (in degree)")

        ("mat.porosity_model", opt::value<int>(&p.mat.porosity_model)->default_value(0),
         "Porosity: 0 fixed, 1 Eulerian Biot volume/grain-pressure evolution.")
        ("mat.permeability_model", opt::value<int>(&p.mat.permeability_model)->default_value(0),
         "Permeability: 0 fixed, 1 Kozeny-Carman, 2 power law from current porosity.")
        ("mat.permeability_exponent", opt::list()->default_value("[3]"),
         "Positive dimensionless exponent n in k=k0*(phi/phi0)^n for permeability_model=2. "
         "Material exponents are averaged by marker count.")
        ("mat.pore_healing_time", opt::list()->default_value("[0]"),
         "Reference recovery time of plastic pore dilation (s); 0 disables recovery. Requires porosity_model=1.")
        ("mat.pore_healing_temperature", opt::list()->default_value("[293]"),
         "Reference absolute temperature for pore recovery (K).")
        ("mat.pore_healing_energy", opt::list()->default_value("[0]"),
         "Nonnegative activation energy for pore recovery (J/mol); tau=tau0*(T/T0)*exp(E/R*(1/T-1/T0)).")
        ("mat.porosity", opt::list()->default_value("[0.0]"),
         "Porosity of the materials '[d0, d1, d2, ...]' (no unit)")
        ("mat.hydraulic_perm", opt::list()->default_value("[1e-14]"),
         "Intrinsic permeability of the materials '[d0, d1, d2, ...]' (m^2)")
        ("mat.matrix_storage", opt::list()->default_value("[0]"),
         "Immobile matrix reservoir storage per bulk volume in 1/Pa.")
        ("mat.matrix_exchange", opt::list()->default_value("[0]"),
         "Explicit exchange coefficient between transport and matrix pressure in 1/(Pa s).")
        ("mat.permeability_tensor", opt::list()->default_value("[1,1,1,0,0,0]"),
         "Dimensionless SPD permeability tensor [xx,yy,zz,xy,xz,yz], multiplying hydraulic_perm.")
        ("mat.fluid_rho0", opt::list()->default_value("[1000]"),
         "Reference fluid density for Darcy gravity and saturated mixture weight '[d0, d1, d2, ...]' (kg/m^3)")
        ("mat.fluid_heat_capacity", opt::list()->default_value("[4180]"),
         "Fluid heat capacity '[d0, d1, ...]' in J/kg/K, for groundwater heat transport.")
        ("mat.derive_biot_coeff_from_bulk_moduli", opt::value<bool>(&p.mat.derive_biot_coeff_from_bulk_moduli)->default_value(false),
         "Use biot_coeff=1-bulk_modulus/bulk_modulus_s with finite positive ordered moduli.")
        ("mat.fluid_alpha", opt::list()->default_value("[2.07e-4]"),
         "Thermal expansivity of the fluids at 0 Pa and 293 K '[d0, d1, d2, ...]' (in K^-1)")
        ("mat.fluid_bulk_modulus", opt::list()->default_value("[2.17e9]"),
         "Bulk modulus of the fluids at 0 Pa and 293 K '[d0, d1, d2, ...]' (in Pa)")
        ("mat.fluid_visc", opt::list()->default_value("[1.002e-3]"),
         "Dynamic viscosity of the fluids at 0 Pa and 293 K '[d0, d1, d2, ...]' (in Pa s).")
        ("mat.biot_coeff", opt::list()->default_value("[1.0]"),
         "Biot-Willis coefficient (HM coupling coeff) '[d0, d1, d2, ...]' (no unit)")
        ("mat.bulk_modulus_s", opt::list()->default_value("[37e9]"),
         "Bulk modulus of the solid grains or minerals '[d0, d1, d2, ...]' (Pa)")
        ("mat.retention_model", opt::list()->default_value("[0]"),
         "Retention law of the materials '[m0, m1, ...]': 0 saturated, "
         "1 van Genuchten-Mualem, 2 Brooks-Corey-Burdine, 3 Gardner exponential.")
        ("mat.retention_alpha", opt::list()->default_value("[1e-4]"),
         "Retention alpha of the materials '[d0, d1, ...]' (1/Pa); "
         "1/alpha is the Brooks-Corey air-entry pressure.")
        ("mat.retention_n", opt::list()->default_value("[2]"),
         "Retention exponent '[d0, d1, ...]': van Genuchten n > 1, "
         "or the Brooks-Corey pore-size index > 0. Unused by Gardner.")
        ("mat.residual_saturation", opt::list()->default_value("[0]"),
         "Residual saturation of the materials '[d0, d1, ...]', 0 <= Sr < 1.")
        // for rate-and-state friction
        ("mat.direct_a", opt::list()->default_value("[0.020]"),
         "rate-and-state friction parameter a '[d0, d1, d2, ...]' (-)")
        ("mat.evolution_b", opt::list()->default_value("[0.025]"),
         "rate-and-state friction parameter b '[d0, d1, d2, ...]' (-)")
        ("mat.characteristic_velocity", opt::list()->default_value("[1e-6]"),
         "rate-and-state friction parameter V0 '[d0, d1, d2, ...]' (-)")
        ("mat.characteristic_distance", opt::list()->default_value("[4e-3]"),
         "rate-and-state friction parameter Dc '[d0, d1, d2, ...]' (-)")
        ("mat.state_var_model", opt::value<int>(&p.mat.state_var_model)->default_value(0),
         "State variable evolution model for rate-and-state friction.\n"
         "0: steady-state friction (no state evolution),\n"
         "1: aging law,\n"
         "2: slip law.")
        ;
    /*
    Example of bulk modulus of rock-forming minerals
    Quartz (SiO₂) 37 GPa
    Feldspar 60-75 GPa
    Olivine ((Mg, Fe)₂SiO₄) 125-130 GPa
    Pyroxene (e.g., Enstatite MgSiO₃, Diopside CaMgSi₂O₆) 100-110 GPa
    Garnet (e.g., Almandine Fe₃Al₂(SiO₄)₃, Pyrope Mg₃Al₂(SiO₄)₃) 170-180 GPa
    Muscovite (a type of mica) 55 GPa
    Calcite (CaCO₃) 70-75 GPa
    Dolomite (CaMg(CO₃)₂) 90 GPa
    Amphibole (e.g., Hornblende) 80-90 GPa
    */
    /* These parameters will enable additional debugging output. DO NOT document these parameters
       in defaults.cfg. */
    cfg.add_options()
        ("debug.dt", opt::value<bool>(&p.debug.dt)->default_value(false),
         "Print all dt criteria")
        ;
}


static void read_parameters_from_file(const char* filename, opt::Registry &cfg)
{
    try {
        cfg.read_file(filename);
    }
    catch (const opt::ConfigError& e) {
        std::cerr << "Error reading config_file '" << filename << "'\n";
        std::cerr << e.what() << "\n";
        die(e.bad_value ? EXIT_CONFIG_VALUE : EXIT_CONFIG);
    }
    catch (const std::exception& e) {
        // Anything the parser did not foresee is still a problem with the file, not
        // a crash of the run.
        std::cerr << "Error reading config_file '" << filename << "'\n";
        std::cerr << e.what() << "\n";
        die(EXIT_CONFIG);
    }
}


template<class T>
static int read_numbers(const std::string &input, std::vector<T> &vec, int len)
{
    /* Read 'len' numbers from input.
     * The format of input must be '[n0, n1, n2]' or '[n0, n1, n2,]' (with a trailing ,),
     * for len=3.
     */

    std::istringstream stream(input);
    vec.resize(len);

    char sentinel;

    stream >> sentinel;
    if (sentinel != '[') return 1;

    for (int i=0; i<len; ++i) {
        stream >> vec[i];

        if (i == len-1) break;

        // consume ','
        char sep;
        stream >> sep;
        if (sep != ',') return 1;
    }

    stream >> sentinel;
    if (sentinel == ',') stream >> sentinel;
    if (sentinel != ']') return 1;

    if (! stream.good()) return 1;

    // success
    return 0;
}


// Read exactly len numbers from key name; len_name says what sets len.
// With one_for_all, a single number is also accepted and repeated len times.
template<class T>
static void get_numbers(const opt::Registry &vm, const char *name,
                        std::vector<T> &values, int len, const char *len_name,
                        bool one_for_all=false)
{
    if ( ! vm.count(name) ) {
        std::cerr << "Error: " << name << " is not provided.\n";
        die(EXIT_CONFIG_VALUE);
    }
    if (len < 0) {
        std::cerr << "Error: " << name << " cannot list " << len << " numbers ("
                  << len_name << ").\n";
        die(EXIT_CONFIG_VALUE);
    }

    std::string str = vm.str(name);
    int err = read_numbers(str, values, len);
    if (err && one_for_all) {
        err = read_numbers(str, values, 1);
        if (!err) {
            const T v = values[0];
            values.assign(len, v);
        }
    }

    if (err) {
        std::cerr << "Error: " << name << " must list " << len
                  << (len == 1 ? " number (" : " numbers (") << len_name << ")"
                  << (one_for_all && len != 1 ? " or 1 for all" : "")
                  << " as '[d0, d1, ...]', got '" << str << "'.\n";
        die(EXIT_CONFIG_VALUE);
    }
}


static void parse_point_source(const opt::Registry &vm, const Param &p,
                               const char *section, Injection &source, bool fluid)
{
    const auto key = [section](const char *name) { return std::string(section)+"."+name; };
    const auto fail = [section](const char *why) {
        const std::string text = std::string(section)+": "+why;
        die(EXIT_CONFIG_VALUE,text.c_str());
    };
    if (source.num_points < 0 || (source.enabled && source.num_points == 0))
        fail("enabled sources require a positive num_points.");
    if (source.enabled && fluid && !p.control.has_hydraulic_diffusion)
        fail("fluid sources require has_hydraulic_diffusion=true.");
    if (source.num_points == 0) return;
    const int n = source.num_points;
    const std::string group = key("num_points");
    const auto list = [&](const char *name, double_vec &v, bool broadcast=false) {
        const std::string name_key = key(name);
        get_numbers(vm,name_key.c_str(),v,n,group.c_str(),broadcast);
    };
    list("points_x",source.points_x);
#ifdef THREED
    list("points_y",source.points_y);
    list("points_z",source.points_z);
#else
    const bool y = !vm.defaulted(key("points_y")), z = !vm.defaulted(key("points_z"));
    if (y == z) fail("supply points_z in 2D; points_y is a mutually exclusive legacy alias.");
    list(z ? "points_z" : "points_y",source.points_z);
#endif
    double scale = 0;
    if (source.points_unit == "m") scale=1;
    else if (source.points_unit == "mm") scale=1e-3;
    else if (source.points_unit == "cm") scale=1e-2;
    else if (source.points_unit == "km") scale=1e3;
    else fail("points_unit must be mm, cm, m or km.");
    const bool box = p.mesh.meshing_option == 1 || p.mesh.meshing_option == 2;
    const auto coordinate = [&](double &v,double lo,double hi) {
        v *= scale;
        if (!std::isfinite(v)) fail("point coordinates must remain finite after unit conversion.");
        if (box) {
            const double tol = 5e-11*std::abs(hi-lo);
            if (v < lo-tol || v > hi+tol) fail("source point lies outside the box mesh.");
            v = std::max(lo,std::min(v,hi));
        }
    };
    for (int i=0;i<n;++i) {
        coordinate(source.points_x[i],0,p.mesh.xlength);
#ifdef THREED
        coordinate(source.points_y[i],0,p.mesh.ylength);
#endif
        coordinate(source.points_z[i],-p.mesh.zlength,0);
    }
    const std::string model = vm.str(key("rate_model"));
    if (model == "constant_rate") source.rate_model=injection_rate_constant;
    else if (model == "total_amount") source.rate_model=injection_rate_total_amount;
    else fail("rate_model must be constant_rate or total_amount.");
    list("rate",source.rate,true);list("total_amount",source.total_amount,true);
    list("start_time_in_yr",source.start_time,true);
    const bool end = !vm.defaulted(key("end_time_in_yr"));
    const bool duration = !vm.defaulted(key("duration_in_yr"));
    if (end && duration) fail("end_time_in_yr and duration_in_yr are mutually exclusive.");
    if (source.rate_model == injection_rate_total_amount && !end && !duration)
        fail("total_amount requires an end time or duration.");
    if (end) list("end_time_in_yr",source.end_time,true);
    else if (duration) {
        double_vec span;list("duration_in_yr",span,true);
        source.end_time.resize(n);
        for (int i=0;i<n;++i) {
            if (!(span[i]>0 && std::isfinite(span[i]))) fail("durations must be finite and positive.");
            source.end_time[i]=source.start_time[i]+span[i];
        }
    } else source.end_time.assign(n,std::numeric_limits<double>::infinity());
    if (fluid && !vm.defaulted(key("temperature"))) {
        if (!p.control.has_fluid_thermal_advection)
            fail("temperature requires control.has_fluid_thermal_advection=true.");
        list("temperature",source.temperature,true);
        for (int i=0;i<n;++i)
            if (!(source.temperature[i] > 0 && std::isfinite(source.temperature[i])))
                fail("injection temperatures must be finite and positive (K).");
    }
    for (int i=0;i<n;++i) {
        if (!std::isfinite(source.rate[i]) || !std::isfinite(source.total_amount[i]) ||
            !std::isfinite(source.start_time[i]) ||
            ((end || duration) && !std::isfinite(source.end_time[i])))
            fail("source strengths and supplied times must be finite.");
        source.start_time[i] *= YEAR2SEC;
        if (end || duration) source.end_time[i] *= YEAR2SEC;
        if (!std::isfinite(source.start_time[i]) ||
            ((end || duration) && !std::isfinite(source.end_time[i])) ||
            !(source.end_time[i]>source.start_time[i])) fail("invalid or overflowing source interval.");
    }
}

static void validate_parameters(const opt::Registry &vm, Param &p)
{
    parse_point_source(vm,p,"injection",p.injection,true);
    parse_point_source(vm,p,"heat_injection",p.heat_injection,false);
    std::cout << "Checking consistency of input parameters...\n";

    if (p.mesh.num_remeshing_reset_boundaries < 0)
        die(EXIT_CONFIG_VALUE, "mesh.num_remeshing_reset_boundaries must be nonnegative.");
    get_numbers(vm, "mesh.remeshing_reset_boundaries", p.mesh.remeshing_reset_boundaries,
                p.mesh.num_remeshing_reset_boundaries, "mesh.num_remeshing_reset_boundaries");
    p.mesh.remeshing_reset_boundary_mask = 0;
    for (int marker : p.mesh.remeshing_reset_boundaries) {
        if (marker <= 0 || marker > static_cast<int>(BOUNDN3) || (marker & (marker-1)))
            die(EXIT_CONFIG_VALUE, "mesh.remeshing_reset_boundaries requires single boundary-marker bits.");
        p.mesh.remeshing_reset_boundary_mask |= static_cast<uint>(marker);
    }
    if (p.mesh.remeshing_option == 14 && !p.mesh.remeshing_reset_boundary_mask)
        die(EXIT_CONFIG_VALUE, "mesh.remeshing_option=14 requires selected boundary-marker bits.");
#ifdef THREED
    if (p.mesh.remeshing_option == 14)
        die(EXIT_UNSUPPORTED_DIM, "mesh.remeshing_option=14 is available only in 2D.");
#endif
    if (p.mesh.remeshing_option == 14 && p.mesh.meshing_elem_shape != 0)
        die(EXIT_CONFIG_VALUE, "mesh.remeshing_option=14 requires an unstructured mesh (meshing_elem_shape=0).");

    //
    // Unspecified thermal side values use bc.surface_temperature and bc.mantle_temperature.
    {
        const char *side[6] = {"x0", "x1", "y0", "y1", "z0", "z1"};
        for (int i = 0; i < 6; ++i) {
            const std::string type = std::string("bc.tbc_") + side[i];
            const std::string value = std::string("bc.tbc_val_") + side[i];
            const bool supplied = !vm.defaulted(value);
            p.bc.custom_temperature_boundary |= supplied || !vm.defaulted(type);
            if (p.bc.tbc[i] < 0 || p.bc.tbc[i] > 1 ||
                (supplied && !std::isfinite(p.bc.tbc_value[i])))
                die(EXIT_CONFIG_VALUE, "Thermal boundaries require tbc=0/1 and finite supplied temperatures.");
        }
    }
    // stopping condition and output interval are based on either model time or step
    //
    if ( ! (vm.count("sim.max_steps") || vm.count("sim.max_time_in_yr")) ) {
        die(EXIT_CONFIG, "Must provide either sim.max_steps or sim.max_time_in_yr");
    }
    if ( ! vm.count("sim.max_steps") )
        p.sim.max_steps = std::numeric_limits<int>::max();
    if ( ! vm.count("sim.max_time_in_yr") )
        p.sim.max_time_in_yr = std::numeric_limits<double>::max();

#ifdef USEMMG
    // An aspect ratio: in (0,1) it would refine the wall-normal axis; a negative one means nothing.
    if ( p.mesh.mmg_aniso_wall_ratio != 0.0 && p.mesh.mmg_aniso_wall_ratio < 1.0 ) {
        std::cerr << "mesh.mmg_aniso_wall_ratio must be 0 (off) or >= 1 (aspect ratio); "
                     "a smaller ratio would refine the wall or mean nothing.\n";
        die(EXIT_CONFIG_VALUE);
    }
#endif

    if ( ! (vm.count("sim.output_step_interval") || vm.count("sim.output_time_interval_in_yr")) ) {
        die(EXIT_CONFIG, "Must provide either sim.output_step_interval or sim.output_time_interval_in_yr");
    }
    if ( ! vm.count("sim.output_step_interval") )
        p.sim.output_step_interval = std::numeric_limits<int>::max();
    if ( ! vm.count("sim.output_time_interval_in_yr") )
        p.sim.output_time_interval_in_yr = std::numeric_limits<double>::max();

    //
    // These parameters are required when restarting
    //
    if (p.sim.is_restarting) {
        if ( ! vm.count("sim.restarting_from_modelname") ) {
            die(EXIT_CONFIG, "Must provide sim.restarting_from_modelname when restarting.");
        }
        if ( ! vm.count("sim.restarting_from_frame") ) {
            die(EXIT_CONFIG, "Must provide sim.restarting_from_frame when restarting.");
        }
    }

    // Both are modulo divisors, the first one just below, so a zero is a SIGFPE.
    if (p.mesh.quality_check_step_interval < 1) {
        die(EXIT_CONFIG_VALUE, "mesh.quality_check_step_interval must be positive.");
    }
    if (p.sim.checkpoint_frame_interval < 1) {
        die(EXIT_CONFIG_VALUE, "sim.checkpoint_frame_interval must be positive.");
    }

    if (p.sim.is_outputting_averaged_fields == true)
        if (vm.count("sim.output_step_interval") &&
            p.sim.output_step_interval%p.mesh.quality_check_step_interval !=0) {
            die(EXIT_CONFIG_VALUE, "sim.output_step_interval must be a multiple of mesh.quality_check_step_interval!.");
    }

    // The status line prints only on quality-check steps, so round the interval up to a
    // multiple of quality_check_step_interval.
    if (p.sim.info_display_step_interval <= 0)
        p.sim.info_display_step_interval = p.mesh.quality_check_step_interval * 100;
    if (p.sim.info_display_step_interval % p.mesh.quality_check_step_interval != 0) {
        int q = p.mesh.quality_check_step_interval;
        p.sim.info_display_step_interval = ((p.sim.info_display_step_interval + q - 1) / q) * q;
        std::cout << "  info_display_step_interval rounded up to "
                  << p.sim.info_display_step_interval
                  << " (must be a multiple of quality_check_step_interval="
                  << q << ")\n";
    }
    if (p.sim.earthquake_output_step_interval < 1) {
        die(EXIT_CONFIG_VALUE, "sim.earthquake_output_step_interval must be >= 1.");
    }
    if (p.sim.earthquake_start_factor <= 0 || p.sim.earthquake_end_factor <= 0) {
        die(EXIT_CONFIG_VALUE, "sim.earthquake_start_factor and sim.earthquake_end_factor must be > 0.");
    }
    if (p.sim.earthquake_start_factor <= p.sim.earthquake_end_factor) {
        std::cerr << "Warning: sim.earthquake_start_factor <= sim.earthquake_end_factor; "
                  << "hysteresis may be weak.\n";
    }

    // these parameters are required in mesh.meshing_elem_shape >= 1
#ifdef THREED
    if (p.mesh.meshing_elem_shape == 2) {
        die(EXIT_UNSUPPORTED_DIM, "mesh.meshing_elem_shape == 2 is not available in 3D.");
    }
#endif
    if (p.mesh.meshing_elem_shape >= 1) {
        if ( p.mesh.meshing_option != 1) {
            die(EXIT_UNSUPPORTED_DIM, "mesh.meshing_elem_shape >= 1 is only for mesh.meshing_option == 1.");
        }
    }

    //
    // mesh.refined_zone{x,y,z}, read only for mesh.meshing_option == 2
    //
    if (p.mesh.meshing_option == 2) {
        /* get 2 numbers from the string */
        double_vec tmp;
        int err;
        std::string str;
        str = vm.str("mesh.refined_zonex");
        err = read_numbers(str, tmp, 2);
        if (err || tmp[0] < 0 || tmp[1] > 1 || tmp[0] > tmp[1]) {
            std::cerr << "Error: incorrect value for mesh.refined_zonex,\n"
                      << "       must be in this format '[d0, d1]', 0 <= d0 <= d1 <= 1.\n";
            die(EXIT_CONFIG_VALUE);
        }
        p.mesh.refined_zonex.first = tmp[0];
        p.mesh.refined_zonex.second = tmp[1];
#ifdef THREED
        str = vm.str("mesh.refined_zoney");
        err = read_numbers(str, tmp, 2);
        if (err || tmp[0] < 0 || tmp[1] > 1 || tmp[0] > tmp[1]) {
            std::cerr << "Error: incorrect value for mesh.refined_zoney,\n"
                      << "       must be in this format '[d0, d1]', 0 <= d0 <= d1 <= 1.\n";
            die(EXIT_CONFIG_VALUE);
        }
        p.mesh.refined_zoney.first = tmp[0];
        p.mesh.refined_zoney.second = tmp[1];
#endif
        str = vm.str("mesh.refined_zonez");
        err = read_numbers(str, tmp, 2);
        if (err || tmp[0] < 0 || tmp[1] > 1 || tmp[0] > tmp[1]) {
            std::cerr << "Error: incorrect value for mesh.refined_zonez,\n"
                      << "       must be in this format '[d0, d1]', 0 <= d0 <= d1 <= 1.\n";
            die(EXIT_CONFIG_VALUE);
        }
        p.mesh.refined_zonez.first = tmp[0];
        p.mesh.refined_zonez.second = tmp[1];
    }

    if (p.mesh.smallest_size > p.mesh.largest_size) {
        die(EXIT_CONFIG_VALUE, "mesh.smallest_size is greater than mesh.largest_size.");
    }

#ifdef USEMMG
    // Only MMG reads the mmg_* keys checked in this block.
    if (p.mesh.mmg_metric_refine_coeff < 0) {
        // 1 + coeff*plstrain divides the target volume: it must stay positive.
        std::cerr << "Error: mesh.mmg_metric_refine_coeff must be non-negative "
                     "(0 disables plastic-strain refinement).\n";
        die(EXIT_CONFIG_VALUE);
    }

    // ratio < 1 puts MMG's floor below the trigger: a just-remeshed element re-triggers at once.
    if (p.mesh.mmg_remesh_defensive_tiny_ratio < 1.0)
        die(EXIT_CONFIG_VALUE, "mesh.mmg_remesh_defensive_tiny_ratio must be >= 1 "
            "(1 = MMG builds at the trigger; larger = a wider gap above it).");

    // MMG needs hmin < hmax: cap the ratio just under largest_size/smallest_size. 1e-9 clears the
    // round-off of mmg_adapt's pow() products (~1e-16) and changes no size that matters. Equal
    // sizes cap it at 1 - 1e-9 < 1, so the >= 1 check must stay above this cap.
    if (p.mesh.meshing_elem_shape == 0 &&
        p.mesh.smallest_size * p.mesh.mmg_remesh_defensive_tiny_ratio >= p.mesh.largest_size) {
        const double cap = (1 - 1e-9) * p.mesh.largest_size / p.mesh.smallest_size;
        std::cerr << "Warning: mesh.mmg_remesh_defensive_tiny_ratio = "
                  << p.mesh.mmg_remesh_defensive_tiny_ratio << " puts MMG's size floor "
                  << "(smallest_size x ratio) at or above largest_size; capped just below "
                  << "largest_size / smallest_size = " << p.mesh.largest_size / p.mesh.smallest_size
                  << ". MMG builds closer to the tiny-element trigger: expect frequent "
                  << "tiny-element remeshes unless smallest_size < largest_size / ratio.\n";
        p.mesh.mmg_remesh_defensive_tiny_ratio = cap;
    }

    if (p.mesh.mmg_remesh_defensive_quality_ratio < 1.0) {
        // ratio < 1 can leave a below-trigger element frozen, re-triggering at every remesh.
        std::cerr << "Error: mesh.mmg_remesh_defensive_quality_ratio must be >= 1 (1 = repair only below the trigger).\n";
        die(EXIT_CONFIG_VALUE);
    }

    // Quality <= 1, so near a product of 1 R1 frees almost every element at every remesh.
    // Threshold 0.9, heuristic; written 0.875 because 0.6 * 1.5 rounds below 0.9.
    if (p.mesh.meshing_elem_shape == 0 &&
        p.mesh.min_quality * p.mesh.mmg_remesh_defensive_quality_ratio >= 0.875) {
        std::cerr << "Warning: mesh.min_quality * mesh.mmg_remesh_defensive_quality_ratio "
                     ">= 0.875: at every remesh the MMG freeze unfreezes every element of quality "
                     "below that product (all of them at >= 1) for MMG to repair, and an element "
                     "MMG rebuilds is re-interpolated; the remesh log's 'R1 repair "
                     "(low-quality/tiny)' count gives their number, tiny elements included.\n";
    }

    if (p.mesh.mmg_remesh_size_recovery_ratio != 0.0 && p.mesh.mmg_remesh_size_recovery_ratio <= 1.0) {
        // ratio <= 1 frees elements on the gradation envelope itself, at every remesh.
        std::cerr << "Error: mesh.mmg_remesh_size_recovery_ratio must be 0 (off) or > 1.\n";
        die(EXIT_CONFIG_VALUE);
    }
#endif

    if (p.mesh.remesh_deborah_min <= 0 ||
        p.mesh.remesh_deborah_min >= p.mesh.remesh_deborah_max) {
        die(EXIT_CONFIG_VALUE, "mesh.remesh_deborah_min must be positive and less than mesh.remesh_deborah_max.");
    }

    //
    // monitor
    //
    {
        if (p.monitor.step_interval < 1) {
            die(EXIT_CONFIG_VALUE, "monitor.step_interval must be >= 1.");
        }
        if (p.monitor.num_points < 0) {
            die(EXIT_CONFIG_VALUE, "monitor.num_points must be >= 0.");
        }

        const char *npts = "monitor.num_points";
        get_numbers(vm, "monitor.points_x", p.monitor.points_x, p.monitor.num_points, npts);
#ifdef THREED
        get_numbers(vm, "monitor.points_y", p.monitor.points_y, p.monitor.num_points, npts);
        get_numbers(vm, "monitor.points_z", p.monitor.points_z, p.monitor.num_points, npts);
#else
        // 2D legacy alias: points_y means z when points_z is not given
        const bool z_from_y = vm.str("monitor.points_z") == "[]" &&
                              vm.str("monitor.points_y") != "[]";
        get_numbers(vm, z_from_y ? "monitor.points_y" : "monitor.points_z", p.monitor.points_z,
                    p.monitor.num_points, npts);
#endif

        if (p.monitor.enabled && p.monitor.num_points <= 0) {
            die(EXIT_CONFIG_VALUE, "monitor.enabled=true requires monitor.num_points > 0.");
        }

        if (p.monitor.points_unit == "mm")
            p.monitor.points_scale_to_m = 1e-3;
        else if (p.monitor.points_unit == "cm")
            p.monitor.points_scale_to_m = 1e-2;
        else if (p.monitor.points_unit == "m")
            p.monitor.points_scale_to_m = 1.0;
        else if (p.monitor.points_unit == "km")
            p.monitor.points_scale_to_m = 1e3;
        else {
            die(EXIT_CONFIG_VALUE, "monitor.points_unit must be one of mm, cm, m, km.");
        }

        for (int i = 0; i < p.monitor.num_points; ++i) {
            p.monitor.points_x[i] *= p.monitor.points_scale_to_m;
#ifdef THREED
            p.monitor.points_y[i] *= p.monitor.points_scale_to_m;
            p.monitor.points_z[i] *= p.monitor.points_scale_to_m;
#else
            p.monitor.points_z[i] *= p.monitor.points_scale_to_m;
#endif
        }

        const std::string point_type = vm.str("monitor.point_type");
        if (point_type == "particle")
            p.monitor.point_type = monitor_point_particle;
        else if (point_type == "grid")
            p.monitor.point_type = monitor_point_grid;
        else
            die(EXIT_CONFIG_VALUE, "monitor.point_type must be particle or grid.");

        const std::string remesh_output = vm.str("monitor.remesh_output");
        if (remesh_output == "before")
            p.monitor.remesh_output = monitor_remesh_before;
        else if (remesh_output == "both")
            p.monitor.remesh_output = monitor_remesh_both;
        else
            die(EXIT_CONFIG_VALUE, "monitor.remesh_output must be before or both.");

        const bool any_monitor_output =
            p.monitor.output_coord ||
            p.monitor.output_velocity ||
            p.monitor.output_force ||
            p.monitor.output_temperature ||
            p.monitor.output_pore_pressure ||
            p.monitor.output_bcflag ||
            p.monitor.output_stress ||
            p.monitor.output_strain ||
            p.monitor.output_strain_rate ||
            p.monitor.output_plastic_strain ||
            p.monitor.output_plastic_strain_rate ||
            p.monitor.output_radiogenic_source ||
            p.monitor.output_density ||
            p.monitor.output_mesh_quality ||
            p.monitor.output_viscosity ||
            p.monitor.output_material ||
            p.monitor.output_dynamic_friction ||
            p.monitor.output_state_variable;

        if (p.monitor.enabled && !any_monitor_output) {
            die(EXIT_CONFIG_VALUE, "monitor.enabled=true requires at least one monitor.output_* = true.");
        }
    }

#ifdef THREED
    if (p.mesh.remeshing_option == 2) {
        die(EXIT_UNSUPPORTED_DIM, "mesh.remeshing_option=2 is not available in 3D.");
    }
#endif
#if defined(USEMMG) && defined(THREED)
    // optimize_mesh (3D MMG) has no flatten_x0_corner case.
    if (p.mesh.remeshing_option == 12 && p.mesh.meshing_elem_shape == 0)
        die(EXIT_UNSUPPORTED_DIM, "mesh.remeshing_option = 12 is not available with MMG remeshing "
            "in 3D (usemmg=1, mesh.meshing_elem_shape = 0); use a build without MMG.");
#endif
#ifdef USEMMG
    // Only new_mesh (Triangle) consumes option 2's new_bottom segment rewrite and deletion list.
    if (p.mesh.remeshing_option == 2 && p.mesh.meshing_elem_shape == 0)
        die(EXIT_CONFIG_VALUE, "mesh.remeshing_option = 2 is not available with MMG remeshing "
            "(usemmg=1, mesh.meshing_elem_shape = 0); use 1 or 11, or a build without MMG.");
#endif
    // MMG remeshing needs boundary-only segments (see drop_unflagged_facets).
    if (!p.mesh.is_discarding_internal_segments && p.mesh.meshing_elem_shape == 0) {
#ifdef USEMMG
        die(EXIT_CONFIG_VALUE, "mesh.is_discarding_internal_segments = no is not available with "
            "MMG remeshing (usemmg=1, mesh.meshing_elem_shape = 0): MMG keeps no internal "
            "segment. Use yes, or a build without MMG.");
#else
        std::cerr << "Warning: mesh.is_discarding_internal_segments = no: a remesh that deletes a "
                     "node of an internal segment (tiny-element repair, node cleanup near a "
                     "restored boundary) misroutes it, which can make the mesher fail or hang; "
                     "in 3D with remeshing_option 1 or 10-12, the first remesh that deletes a "
                     "boundary node drops them all.\n";
#endif
    }

#ifndef USEMMG
    if (p.mesh.use_mmg_init) {
        p.mesh.use_mmg_init = false;
        std::cerr << "Warning: MMG not compiled in (no USEMMG), mesh.use_mmg_init turned off.\n";
    }
    // Without USEMMG nothing reads the mmg_* keys: name the ones the cfg set, not the defaulted ones.
    {
        const char *mmg_only[] = {
            "mesh.mmg_remesh_active_plstrain", "mesh.mmg_remesh_defensive_quality_ratio",
            "mesh.mmg_remesh_size_recovery_ratio", "mesh.mmg_debug", "mesh.mmg_verbose",
            "mesh.mmg_hausd_factor", "mesh.mmg_init_coarsening_factor",
            "mesh.mmg_metric_refine_coeff", "mesh.mmg_metric_refine_recover",
            "mesh.mmg_aniso_wall_ratio",
            "mesh.mmg_remesh_defensive_tiny_ratio",
        };
        for (const char *name : mmg_only)
            if (vm.count(name) && ! vm.defaulted(name))
                std::cerr << "Warning: MMG not compiled in (no USEMMG), " << name
                          << " has no effect.\n";
    }
    // Option 13 re-flattens the side walls at every remesh, which breeds sub-hmin wall elements. Only
    // MMG collapses them: decimate_below_hmin keeps boundary nodes, and the tiny-element repair
    // skips corners.
    if (p.mesh.remeshing_option == 13 && p.mesh.meshing_elem_shape == 0)
        die(EXIT_UNSUPPORTED_LIB, "mesh.remeshing_option = 13 requires MMG (build with usemmg=1). "
            "Without it the restored side wall accumulates sub-hmin elements that the "
            "Triangle/TetGen path cannot collapse. Uniform meshes (mesh.meshing_elem_shape > 0) "
            "are unaffected.");
#endif

    //
    // bc
    //
    {
        if ( p.bc.has_winkler_foundation && p.control.gravity == 0 ) {
            p.bc.has_winkler_foundation = 0;
            std::cerr << "Warning: no gravity, Winkler foundation is turned off.\n";
        }
        if ( p.bc.has_winkler_foundation && p.bc.vbc_z0 != 0 ) {
            p.bc.vbc_z0 = 0;
            std::cerr << "Warning: Winkler foundation is turned on, setting bc.vbc_z0 to 0.\n";
        }
        if ( p.bc.has_water_loading && p.control.gravity == 0 ) {
            p.bc.has_water_loading = 0;
            std::cerr << "Warning: no gravity, water loading is turned off.\n";
        }
        if ( p.bc.has_water_loading && p.bc.vbc_z1 != 0 ) {
            p.bc.vbc_z1 = 0;
            std::cerr << "Warning: water loading is turned on, setting bc.vbc_z1 to 0.\n";
        }

#ifdef THREED
        if ( p.bc.vbc_z0 > 3) {
            die(EXIT_CONFIG_VALUE, "bc.vbc_z0 is not 0, 1, 2, or 3.");
        }
        if ( p.bc.vbc_z1 > 3 && p.bc.vbc_z1 != 8) {
            die(EXIT_CONFIG_VALUE, "bc.vbc_z1 is not 0, 1, 2, 3, or 8.");
        }
#else
        if ( p.bc.vbc_z0 > 4) {
            die(EXIT_CONFIG_VALUE, "bc.vbc_z0 is not 0, 1, 2, 3, or 4.");
        }
        if ( p.bc.vbc_z1 > 4 && p.bc.vbc_z1 != 8) {
            die(EXIT_CONFIG_VALUE, "bc.vbc_z1 is not 0, 1, 2, 3, 4, or 8.");
        }
#endif
        if ( p.bc.vbc_n0 != 1 && p.bc.vbc_n0 != 3 && p.bc.vbc_n0 != 11 && p.bc.vbc_n0 != 13 ) {
            die(EXIT_CONFIG_VALUE, "bc.vbc_n0 is not 1, 3, 11, or 13.");
        }
        if ( p.bc.vbc_n1 != 1 && p.bc.vbc_n1 != 3 && p.bc.vbc_n1 != 11 && p.bc.vbc_n1 != 13 ) {
            die(EXIT_CONFIG_VALUE, "bc.vbc_n1 is not 1, 3, 11, or 13.");
        }
        if ( p.bc.vbc_n2 != 1 && p.bc.vbc_n2 != 3 && p.bc.vbc_n2 != 11 && p.bc.vbc_n2 != 13 ) {
            die(EXIT_CONFIG_VALUE, "bc.vbc_n2 is not 1, 3, 11, or 13.");
        }
        if ( p.bc.vbc_n3 != 1 && p.bc.vbc_n3 != 3 && p.bc.vbc_n3 != 11 && p.bc.vbc_n3 != 13 ) {
            die(EXIT_CONFIG_VALUE, "bc.vbc_n3 is not 1, 3, 11, or 13.");
        }
    }

    //
    // control
    //
    {
        if ( p.control.dt_fraction < 0 || p.control.dt_fraction > 1 ) {
            die(EXIT_CONFIG_VALUE, "control.dt_fraction must be between 0 and 1.");
        }
        if ( p.control.damping_factor < 0 || p.control.damping_factor > 1 ) {
            die(EXIT_CONFIG_VALUE, "control.damping_factor must be between 0 and 1.");
        }
        // ref_pressure() is 'acc routine seq' and cannot reject an option itself; it
        // returns 0. Keep this range in step with its branches and the help above.
        if (p.control.ref_pressure_option < 0 || p.control.ref_pressure_option > 2) {
            std::cerr << "Error: control.ref_pressure_option must be 0, 1, or 2 (got "
                      << p.control.ref_pressure_option << ").\n";
            die(EXIT_CONFIG_VALUE);
        }
        if (p.control.rsf_slip_rate_projection_option !=
                rsf_slip_rate_projection_maximum_shear &&
            p.control.rsf_slip_rate_projection_option !=
                rsf_slip_rate_projection_total_strain_rate) {
            die(EXIT_CONFIG_VALUE,
                "control.rsf_slip_rate_projection_option must be 0 or 1.");
        }
        if (!std::isfinite(p.control.rsf_dtheta_max) ||
            p.control.rsf_dtheta_max < 0 ||
            p.control.rsf_dtheta_max >= 2) {
            die(EXIT_CONFIG_VALUE,
                "control.rsf_dtheta_max must be finite and in [0, 2).");
        }

    }

    //
    // ic
    //
    {
        if ( p.ic.weakzone_option == 5 ) {
            int n = p.ic.weakzone_num_segments;
            const char *nseg = "ic.weakzone_num_segments";
            get_numbers(vm, "ic.weakzone_segments_xcenter",     p.ic.weakzone_segments_xcenter,     n, nseg, true);
            get_numbers(vm, "ic.weakzone_segments_zcenter",     p.ic.weakzone_segments_zcenter,     n, nseg, true);
            get_numbers(vm, "ic.weakzone_segments_azimuth",     p.ic.weakzone_segments_azimuth,     n, nseg, true);
            get_numbers(vm, "ic.weakzone_segments_inclination", p.ic.weakzone_segments_inclination, n, nseg, true);
            get_numbers(vm, "ic.weakzone_segments_halfwidth",   p.ic.weakzone_segments_halfwidth,   n, nseg, true);
            get_numbers(vm, "ic.weakzone_segments_x_min",       p.ic.weakzone_segments_x_min,       n, nseg, true);
            get_numbers(vm, "ic.weakzone_segments_x_max",       p.ic.weakzone_segments_x_max,       n, nseg, true);
            get_numbers(vm, "ic.weakzone_segments_depth_min",   p.ic.weakzone_segments_depth_min,   n, nseg, true);
            get_numbers(vm, "ic.weakzone_segments_depth_max",   p.ic.weakzone_segments_depth_max,   n, nseg, true);
#ifdef THREED
            get_numbers(vm, "ic.weakzone_segments_ycenter",     p.ic.weakzone_segments_ycenter,     n, nseg, true);
            get_numbers(vm, "ic.weakzone_segments_y_min",       p.ic.weakzone_segments_y_min,       n, nseg, true);
            get_numbers(vm, "ic.weakzone_segments_y_max",       p.ic.weakzone_segments_y_max,       n, nseg, true);
#endif
        }

        if ( p.ic.mattype_option == 1) {
            if (p.ic.num_mattype_layers < 1) {
                die(EXIT_CONFIG_VALUE, "ic.num_mattype_layers must be >= 1.");
            }
            get_numbers(vm, "ic.layer_mattypes", p.ic.layer_mattypes, p.ic.num_mattype_layers, "ic.num_mattype_layers");
            get_numbers(vm, "ic.mattype_layer_depths", p.ic.mattype_layer_depths, p.ic.num_mattype_layers-1, "ic.num_mattype_layers - 1");
            // mattype_layer_depths must be already sorted
            if (! std::is_sorted(p.ic.mattype_layer_depths.begin(), p.ic.mattype_layer_depths.end())) {
                std::cerr << "Error: the content of ic.mattype_layer_depths is not ordered from"
                    " small to big values.\n";
                die(EXIT_CONFIG_VALUE);
            }
        }

        if (!std::isfinite(p.ic.water_table_z))
            die(EXIT_CONFIG_VALUE, "ic.water_table_z must be finite.");
        if (p.ic.temperature_option == 4 || p.ic.temperature_option == 5) {
            if (p.ic.temperature_option == 5 &&
                !(std::isfinite(p.bc.surface_temperature) && std::isfinite(p.bc.mantle_temperature)))
                die(EXIT_CONFIG_VALUE, "Linear conduction requires finite top and bottom temperatures.");
            const double base = p.ic.temperature_option == 4 ? p.ic.temp_background :
                std::min(p.bc.surface_temperature, p.bc.mantle_temperature);
            const double values[] = {base,p.ic.temp_anomaly_amplitude,
                p.ic.temp_anomaly_xcenter,p.ic.temp_anomaly_ycenter,p.ic.temp_anomaly_zcenter,
                p.ic.temp_anomaly_xsemi_axis,p.ic.temp_anomaly_ysemi_axis,p.ic.temp_anomaly_zsemi_axis};
            for (double v : values)
                if (!std::isfinite(v)) die(EXIT_CONFIG_VALUE, "Gaussian temperature parameters must be finite.");
            if (!(base+std::min(0.0,p.ic.temp_anomaly_amplitude) > 0 &&
                  p.ic.temp_anomaly_xsemi_axis > 0 && p.ic.temp_anomaly_ysemi_axis > 0 &&
                  p.ic.temp_anomaly_zsemi_axis > 0))
                die(EXIT_CONFIG_VALUE, "Gaussian temperature requires positive Kelvin values and scales.");
        }
        if (p.ic.initial_stress_option != 0 &&
            p.ic.initial_stress_option != 1) {
            die(EXIT_CONFIG_VALUE,
                "ic.initial_stress_option must be 0 or 1.");
        }
        get_numbers(vm, "ic.initial_stress", p.ic.initial_stress, NSTR, "one per stress component");
        if (p.ic.initial_stress_option == 1) {
            if (p.control.gravity != 0.0) {
                die(EXIT_CONFIG_VALUE,
                    "ic.initial_stress_option=1 prescribes an absolute stress "
                    "tensor and requires control.gravity=0.");
            }
            for (int i = 0; i < NSTR; ++i) {
                if (!std::isfinite(p.ic.initial_stress[i])) {
                    die(EXIT_CONFIG_VALUE,
                        "ic.initial_stress must contain only finite values.");
                }
            }
        }

        if (p.ic.temperature_option == 3) {
            if (p.ic.radiogenic_heat_dome_width == 0) {
                die(EXIT_CONFIG_VALUE, "ic.radiogenic_heat_dome_width must be greater than 0 for ic.temperature_option=3.");
            }
        }
    }

    //
    // marker
    //
    {

    }

    //
    // material properties
    //
    {
        const std::string str =
            vm.str("control.mass_scaling_reference_speed");
        if (str == std::string("shear"))
            p.control.mass_scaling_reference_speed =
                mass_scaling_speed_shear;
        else if (str == std::string("bulk"))
            p.control.mass_scaling_reference_speed =
                mass_scaling_speed_bulk;
        else {
            std::cerr
                << "Error: control.mass_scaling_reference_speed must be "
                   "'shear' or 'bulk', not '" << str << "'\n";
            die(EXIT_CONFIG_VALUE);
        }
    }

    {
        std::string str = vm.str("mat.rheology_type");
        if (str == std::string("elastic"))
            p.mat.rheol_type = MatProps::rh_elastic;
        else if (str == std::string("viscous"))
            p.mat.rheol_type = MatProps::rh_viscous;
        else if (str == std::string("maxwell"))
            p.mat.rheol_type = MatProps::rh_maxwell;
        else if (str == std::string("elasto-plastic"))
            p.mat.rheol_type = MatProps::rh_ep;
        else if (str == std::string("elasto-visco-plastic"))
            p.mat.rheol_type = MatProps::rh_evp;
        else if (str == std::string("elasto-plastic-rate-state-friction") ||
                 str == std::string("elasto-plastic-rsf"))
            p.mat.rheol_type = MatProps::rh_ep_rsf;
        else if (str == std::string("elasto-visco-plastic-rate-state-friction") ||
                 str == std::string("elasto-visco-plastic-rsf"))
            p.mat.rheol_type = MatProps::rh_evp_rsf;
        else {
            std::cerr << "Error: unknown rheology: '" << str << "'\n";
            die(EXIT_CONFIG_VALUE);
        }

        if ((p.mat.rheol_type & MatProps::rh_rsf) && !p.control.use_global_velocity_scaling) {
            p.control.use_global_velocity_scaling = true;
            std::cerr << "Warning: RSF rheology requires control.use_global_velocity_scaling=true. Forcing it on.\n";
        }
        if (p.control.mass_scaling_reference_speed ==
                mass_scaling_speed_bulk &&
            !p.control.use_global_velocity_scaling) {
            die(EXIT_CONFIG_VALUE,
                "control.mass_scaling_reference_speed=bulk requires "
                "control.use_global_velocity_scaling=true.");
        }

#ifdef THREED
        if ( p.mat.is_plane_strain ) {
            p.mat.is_plane_strain = false;
            std::cerr << "Warning: mat.is_plane_strain is not avaiable in 3D.\n";
        }
#endif

        // phase_changes() implements only these; any other value would run to its first
        // phase-change step and die there.
        switch (p.mat.phase_change_option) {
        case 0: case 1: case 3: case 101:
            break;
        default:
            die(EXIT_CONFIG_VALUE, "mat.phase_change_option must be 0, 1, 3 or 101 (see --help).");
        }
        if (p.mat.phase_change_pressure_option != 0 && p.mat.phase_change_pressure_option != 1) {
            die(EXIT_CONFIG_VALUE, "mat.phase_change_pressure_option must be 0 or 1 (see --help).");
        }
        if (p.mat.phase_change_option != 0 && p.mat.nmat == 1) {
            die(EXIT_CONFIG_VALUE, "mat.phase_change_option is chosen, but mat.num_materials is 1.");
        }
        if (p.mat.phase_change_option == 1 && p.mat.nmat < 8) {
            die(EXIT_CONFIG_VALUE, "mat.phase_change_option is 1, but mat.num_materials is less than 8.");
        }
        // A base key left out of the cfg keeps its default 0, which names a real
        // material, so a rule that reads it would convert material 0. Every key
        // option 3 reads must be set. Checked before the alias checks below, which
        // would blame the symptom; a key whose default -1 is unset is left to the
        // range checks and check_phase_transition_rules(), which name it.
        if (p.mat.phase_change_option == 3)
            for (const char* key : phase_change_keys(p)) {
                const std::string full = std::string("mat.") + key;
                const int v = vm.as<int>(full);
                if (vm.defaulted(full) && v >= 0) {
                    const std::string msg = std::string("mat.phase_change_option 3 reads mat.") + key
                        + ", but the cfg does not set it, and its default " + std::to_string(v)
                        + " would name material " + std::to_string(v) + ".";
                    die(EXIT_CONFIG_VALUE, msg.c_str());
                }
            }
        // The serpentinization grid's inputs, checked only where it is built.
        if (uses_serpentinization_grid(p)) {
            const std::pair<const char*, int> serp_types[] = {
                {"mattype_mantle", p.mat.mattype_mantle},
                {"mattype_serpentinite", p.mat.mattype_serpentinite},
                {"mattype_oceanic_crust", p.mat.mattype_oceanic_crust},
                {"mattype_sediment", p.mat.mattype_sediment},
            };
            for (const auto& t : serp_types)
                if (t.second < 0 || t.second >= p.mat.nmat) {
                    const std::string msg = std::string("mat.phase_change_option 3 with "
                        "mat.has_serpentinization needs mat.") + t.first + " in [0, mat.num_materials).";
                    die(EXIT_CONFIG_VALUE, msg.c_str());
                }
            // The grid counts these as wet crust (counts_as_wet_crust() in phasechanges.cxx);
            // the mantle or serpentinite among them would trigger itself.
            const int wet_crust_types[] = { p.mat.mattype_oceanic_crust, p.mat.mattype_sediment,
                                            p.mat.mattype_sediment_rock };
            bool aliased = p.mat.mattype_mantle == p.mat.mattype_serpentinite;
            for (int t : wet_crust_types)
                aliased = aliased || t == p.mat.mattype_mantle || t == p.mat.mattype_serpentinite;
            if (aliased) {
                die(EXIT_CONFIG_VALUE, "mat.phase_change_option is 3, but mattype_mantle, mattype_serpentinite "
                    "and the wet-crust materials (mattype_oceanic_crust, mattype_sediment, and "
                    "mattype_sediment_rock when set) are not distinct.");
            }
            if (p.mat.serpentinization_thickness <= 0) {
                die(EXIT_CONFIG_VALUE, "mat.serpentinization_thickness must be positive.");
            }
            if (p.mat.serpentinization_max_depth <= 0) {
                die(EXIT_CONFIG_VALUE, "mat.serpentinization_max_depth must be positive.");
            }
            // The dominance test is a strict '>': at 1 no cell could be dominant.
            if (p.mat.wet_crust_fraction_threshold < 0 || p.mat.wet_crust_fraction_threshold >= 1) {
                die(EXIT_CONFIG_VALUE, "mat.wet_crust_fraction_threshold must be in [0, 1).");
            }
            if (p.mesh.vt_size <= 0) {
                die(EXIT_CONFIG_VALUE, "mesh.vt_size must be positive.");
            }
            // The wet-crust-below search steps down in whole grid cells; a thickness
            // under one cell examines no cell at all and serpentinization never fires.
            if (p.mat.serpentinization_thickness < p.mesh.resolution * p.mesh.vt_size) {
                die(EXIT_CONFIG_VALUE, "mat.serpentinization_thickness is smaller than one virtual-topology "
                    "cell (mesh.resolution * mesh.vt_size); the wet-crust-below search would examine no cell. "
                    "Raise the thickness or lower mesh.vt_size.");
            }
        }
        if (p.mat.phase_change_option == 3) {
            // Metamorphic products for the phase-change rule engine
            // (subduction_phase_changes): -1 (the default) or an in-range
            // material index. check_phase_transition_rules() then requires
            // the ones an active switch reads to be set.
            const int metamorphic_types[] = {
                p.mat.mattype_eclogite, p.mat.mattype_hydrous_mantle,
                p.mat.mattype_sediment_rock, p.mat.mattype_metasediment,
                p.mat.mattype_schist };
            for (int t : metamorphic_types)
                if (t != -1 && (t < 0 || t >= p.mat.nmat)) {
                    die(EXIT_CONFIG_VALUE, "mat.phase_change_option is 3, but a metamorphic "
                        "mat.mattype_* index is neither -1 (unset) nor in [0, mat.num_materials).");
                }
            check_phase_transition_rules(p);
        }

        if (p.mat.nmat < 1) {
            die(EXIT_CONFIG_VALUE, "mat.num_materials must be greater than 0.");
        }

        // Each base mat.mattype_* below indexes the per-material arrays (rho0,
        // therm_cond, elemmarkers[e], ...) unchecked, so each must name an existing
        // material. The rule engine's indices (-1 = unset) are checked with
        // mat.phase_change_option 3 above, mattype_continental_crust with the
        // continental geotherm below.
        const std::pair<const char*, int> mattypes[] = {
            {"mattype_ref", p.mat.mattype_ref},
            {"mattype_mantle", p.mat.mattype_mantle},
            {"mattype_depleted_mantle", p.mat.mattype_depleted_mantle},
            {"mattype_partial_melting_mantle", p.mat.mattype_partial_melting_mantle},
            {"mattype_sediment", p.mat.mattype_sediment},
            {"mattype_mor_extrusion", p.mat.mattype_mor_extrusion},
            {"mattype_asthenosphere", p.mat.mattype_asthenosphere},
        };
        for (const auto& m : mattypes) {
            if (m.second < 0 || m.second >= p.mat.nmat) {
                std::cerr << "Error: mat." << m.first << " (" << m.second
                          << ") must be within [0, mat.num_materials-1] = [0, "
                          << p.mat.nmat - 1 << "].\n";
                die(EXIT_CONFIG_VALUE);
            }
        }
        // The continental geotherm indexes rho0 and therm_cond with it.
        if (p.ic.temperature_option == 1 &&
                (p.mat.mattype_continental_crust < 0 || p.mat.mattype_continental_crust >= p.mat.nmat)) {
            const std::string msg = "ic.temperature_option 1 needs mat.mattype_continental_crust ("
                + std::to_string(p.mat.mattype_continental_crust) + ") in [0, mat.num_materials).";
            die(EXIT_CONFIG_VALUE, msg.c_str());
        }

        if (p.mat.nmat == 1 && p.control.ref_pressure_option != 0) {
            p.control.ref_pressure_option = 0;
            std::cerr << "Warning: mat.num_materials is 1, using simplest control.ref_pressure_option.\n";
        }
        if (p.mat.nmat == 1 && p.markers.replenishment_option != 1) {
            p.markers.replenishment_option = 1;
            std::cerr << "Warning: mat.num_materials is 1, using simplest markers.replenishment_option.\n";
        }

        // Each period table lists exactly num_vbc_period_* (>= 1) points in non-decreasing time:
        // interp1 reads t[0] and ratio[i+1] for each time but the last; a repeated time is a step.
        if (std::min({p.bc.num_vbc_period_x0, p.bc.num_vbc_period_x1,
                      p.bc.num_vbc_period_y0, p.bc.num_vbc_period_y1}) < 1) {
            die(EXIT_CONFIG_VALUE, "bc.num_vbc_period_* must be positive.");
        }
        get_numbers(vm, "bc.vbc_period_x0_time_in_yr", p.bc.vbc_period_x0_time_in_yr, p.bc.num_vbc_period_x0, "bc.num_vbc_period_x0");
        get_numbers(vm, "bc.vbc_period_x1_time_in_yr", p.bc.vbc_period_x1_time_in_yr, p.bc.num_vbc_period_x1, "bc.num_vbc_period_x1");
        get_numbers(vm, "bc.vbc_period_y0_time_in_yr", p.bc.vbc_period_y0_time_in_yr, p.bc.num_vbc_period_y0, "bc.num_vbc_period_y0");
        get_numbers(vm, "bc.vbc_period_y1_time_in_yr", p.bc.vbc_period_y1_time_in_yr, p.bc.num_vbc_period_y1, "bc.num_vbc_period_y1");
        get_numbers(vm, "bc.vbc_period_x0_ratio", p.bc.vbc_period_x0_ratio, p.bc.num_vbc_period_x0, "bc.num_vbc_period_x0");
        get_numbers(vm, "bc.vbc_period_x1_ratio", p.bc.vbc_period_x1_ratio, p.bc.num_vbc_period_x1, "bc.num_vbc_period_x1");
        get_numbers(vm, "bc.vbc_period_y0_ratio", p.bc.vbc_period_y0_ratio, p.bc.num_vbc_period_y0, "bc.num_vbc_period_y0");
        get_numbers(vm, "bc.vbc_period_y1_ratio", p.bc.vbc_period_y1_ratio, p.bc.num_vbc_period_y1, "bc.num_vbc_period_y1");
        const std::pair<const char*, const double_vec*> period_times[] = {
            {"vbc_period_x0_time_in_yr", &p.bc.vbc_period_x0_time_in_yr},
            {"vbc_period_x1_time_in_yr", &p.bc.vbc_period_x1_time_in_yr},
            {"vbc_period_y0_time_in_yr", &p.bc.vbc_period_y0_time_in_yr},
            {"vbc_period_y1_time_in_yr", &p.bc.vbc_period_y1_time_in_yr},
        };
        for (const auto& t : period_times) {
            if (! std::is_sorted(t.second->begin(), t.second->end())) {
                std::cerr << "Error: bc." << t.first << " must be non-decreasing.\n";
                die(EXIT_CONFIG_VALUE);
            }
        }

        if (p.bc.num_wells < 0) {
            die(EXIT_CONFIG_VALUE, "bc.num_wells must be >= 0.");
        }
        get_numbers(vm, "bc.hbc_pressure", p.bc.hbc_pressure, 6, "six boundaries");
        {
            const char *side[6] = {"x0","x1","y0","y1","z0","z1"};
            int *type[6] = {&p.bc.hbc_x0,&p.bc.hbc_x1,&p.bc.hbc_y0,&p.bc.hbc_y1,&p.bc.hbc_z0,&p.bc.hbc_z1};
            for (int i = 0; i < 6; ++i) {
                if (vm.defaulted(std::string("bc.hbc_val_")+side[i])) continue;
                const double value = p.bc.hbc_legacy_value[i];
                if (!std::isfinite(value) || !vm.defaulted("bc.hbc_pressure"))
                    die(EXIT_CONFIG_VALUE, "Use either hbc_val_* compatibility values or hbc_pressure, with finite values.");
                p.bc.hbc_pressure[i] = value;
                if (*type[i] == 1) *type[i] = 2;
            }
            if (p.ic.pore_pressure_option < 0 || p.ic.pore_pressure_option > 2 ||
                !std::isfinite(p.ic.pore_pressure_value))
                die(EXIT_CONFIG_VALUE, "Invalid initial pressure option/value.");
            if (p.ic.pore_pressure_option == 1 &&
                !(p.bc.hbc_x0 == 2 && p.bc.hbc_x1 == 2 && p.mesh.xlength > 0))
                die(EXIT_CONFIG_VALUE, "Linear initial pressure needs prescribed pressures on x0 and x1.");
        }
        get_numbers(vm, "bc.hbc_flux", p.bc.hbc_flux, 6, "six boundaries");
        if (p.bc.num_hbc_schedule < 1)
            die(EXIT_CONFIG_VALUE, "bc.num_hbc_schedule must be positive.");
        get_numbers(vm, "bc.hbc_schedule_time_in_yr", p.bc.hbc_schedule_time_in_yr,
                    p.bc.num_hbc_schedule, "bc.num_hbc_schedule");
        get_numbers(vm, "bc.hbc_schedule_ratio", p.bc.hbc_schedule_ratio,
                    p.bc.num_hbc_schedule, "bc.num_hbc_schedule");
        if (!std::is_sorted(p.bc.hbc_schedule_time_in_yr.begin(), p.bc.hbc_schedule_time_in_yr.end()))
            die(EXIT_CONFIG_VALUE, "Boundary schedule times must be nondecreasing.");
        const int hbc[] = {p.bc.hbc_x0,p.bc.hbc_x1,p.bc.hbc_y0,p.bc.hbc_y1,p.bc.hbc_z0,p.bc.hbc_z1};
        for (int i = 0; i < 6; ++i) {
            if (hbc[i] < 0 || hbc[i] > 3)
                die(EXIT_CONFIG_VALUE, "Hydraulic boundary types must be 0, 1, 2 or 3.");
            for (int j = 0; j < i; ++j)
                if ((p.mesh.meshing_option == 1 || p.mesh.meshing_option == 2) &&
                    hbc[i] == 2 && hbc[j] == 2 && i/2 != j/2 &&
                    p.bc.hbc_pressure[i] != p.bc.hbc_pressure[j])
                    die(EXIT_CONFIG_VALUE, "Intersecting prescribed-pressure faces must agree at their shared nodes.");
        }
        for (const double_vec *v : {&p.bc.hbc_pressure, &p.bc.hbc_flux,
                                   &p.bc.hbc_schedule_time_in_yr, &p.bc.hbc_schedule_ratio})
            for (double x : *v)
                if (!std::isfinite(x)) die(EXIT_CONFIG_VALUE, "Hydraulic boundary data must be finite.");
        if (p.bc.num_wells > 0) {
            if (! p.control.has_hydraulic_diffusion) {
                die(EXIT_CONFIG_VALUE, "bc.num_wells > 0 requires control.has_hydraulic_diffusion.");
            }
            const char *nw = "bc.num_wells";
            get_numbers(vm, "bc.well_x", p.bc.well_x, p.bc.num_wells, nw);
#ifdef THREED
            get_numbers(vm, "bc.well_y", p.bc.well_y, p.bc.num_wells, nw);
#endif
            get_numbers(vm, "bc.well_z_top", p.bc.well_z_top, p.bc.num_wells, nw);
            get_numbers(vm, "bc.well_z_bottom", p.bc.well_z_bottom, p.bc.num_wells, nw);
            get_numbers(vm, "bc.well_rate", p.bc.well_rate, p.bc.num_wells, nw);
            for (int w = 0; w < p.bc.num_wells; ++w) {
                if (! (p.bc.well_z_bottom[w] < p.bc.well_z_top[w]) ||
                    ! std::isfinite(p.bc.well_z_top[w] - p.bc.well_z_bottom[w]) ||
                    ! std::isfinite(p.bc.well_rate[w])) {
                    std::cerr << "Error: well " << w << " needs well_z_bottom < well_z_top and a finite well_rate.\n";
                    die(EXIT_CONFIG_VALUE);
                }
            }
            if (p.bc.num_well_schedule < 1) {
                die(EXIT_CONFIG_VALUE, "bc.num_well_schedule must be positive.");
            }
            const char *ns = "bc.num_well_schedule";
            get_numbers(vm, "bc.well_schedule_time_in_yr", p.bc.well_schedule_time_in_yr, p.bc.num_well_schedule, ns);
            get_numbers(vm, "bc.well_schedule_ratio", p.bc.well_schedule_ratio, p.bc.num_well_schedule, ns);
            for (const double_vec *v : {&p.bc.well_x, &p.bc.well_y, &p.bc.well_z_top,
                                       &p.bc.well_z_bottom, &p.bc.well_rate,
                                       &p.bc.well_schedule_time_in_yr, &p.bc.well_schedule_ratio})
                for (double x : *v)
                    if (!std::isfinite(x)) die(EXIT_CONFIG_VALUE, "Well coordinates, rates and schedules must be finite.");
            if (! std::is_sorted(p.bc.well_schedule_time_in_yr.begin(), p.bc.well_schedule_time_in_yr.end())) {
                die(EXIT_CONFIG_VALUE, "bc.well_schedule_time_in_yr must be non-decreasing.");
            }
        }

        if (p.ic.nhlayer < 0) {
            die(EXIT_CONFIG_VALUE, "ic.num_radiogenic_heat_layer must be >= 0.");
        }
        get_numbers(vm, "ic.radiogenic_heat_boundry", p.ic.radiogenic_heat_boundry, p.ic.nhlayer+1, "ic.num_radiogenic_heat_layer + 1");
        if (p.ic.radiogenic_heat_boundry[0] == -1)
            p.ic.radiogenic_heat_boundry[0] = 0;
        if (p.ic.radiogenic_heat_boundry[p.ic.nhlayer] == -1)
            p.ic.radiogenic_heat_boundry[p.ic.nhlayer] = p.mesh.zlength;

        get_numbers(vm, "ic.radiogenic_heat_mat_in_layer", p.ic.radiogenic_heat_mat_in_layer, p.ic.nhlayer, "ic.num_radiogenic_heat_layer", true);

        const char *nmat_key = "mat.num_materials";
        get_numbers(vm, "mat.rho0", p.mat.rho0, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.alpha", p.mat.alpha, p.mat.nmat, nmat_key, true);

        get_numbers(vm, "mat.bulk_modulus", p.mat.bulk_modulus, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.shear_modulus", p.mat.shear_modulus, p.mat.nmat, nmat_key, true);

        get_numbers(vm, "mat.visc_exponent", p.mat.visc_exponent, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.visc_coefficient", p.mat.visc_coefficient, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.visc_activation_energy", p.mat.visc_activation_energy, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.visc_activation_volume", p.mat.visc_activation_volume, p.mat.nmat, nmat_key, true);

        get_numbers(vm, "mat.heat_capacity", p.mat.heat_capacity, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.therm_cond", p.mat.therm_cond, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.radiogenic_heat_prod", p.mat.radiogenic_heat_prod, p.mat.nmat, nmat_key, true);

        get_numbers(vm, "mat.pls0", p.mat.pls0, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.pls1", p.mat.pls1, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.cohesion0", p.mat.cohesion0, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.cohesion1", p.mat.cohesion1, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.friction_angle0", p.mat.friction_angle0, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.friction_angle1", p.mat.friction_angle1, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.dilation_angle0", p.mat.dilation_angle0, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.dilation_angle1", p.mat.dilation_angle1, p.mat.nmat, nmat_key, true);

        // Hydraulic parameters
        get_numbers(vm, "mat.porosity", p.mat.porosity, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.hydraulic_perm", p.mat.hydraulic_perm, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.permeability_exponent", p.mat.permeability_exponent, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.pore_healing_time", p.mat.pore_healing_time, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.pore_healing_temperature", p.mat.pore_healing_temperature, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.pore_healing_energy", p.mat.pore_healing_energy, p.mat.nmat, nmat_key, true);
        for (int m = 0; m < p.mat.nmat; ++m) {
            const double tau = p.mat.pore_healing_time[m], t0 = p.mat.pore_healing_temperature[m];
            const double energy = p.mat.pore_healing_energy[m];
            if (!(tau >= 0 && std::isfinite(tau) && t0 > 0 && std::isfinite(t0) &&
                  energy >= 0 && std::isfinite(energy)))
                die(EXIT_CONFIG_VALUE, "Pore recovery requires finite nonnegative time/activation energy and positive reference temperature.");
            if (tau > 0 && (p.mat.porosity_model != 1 ||
                            !(p.mat.rheol_type & MatProps::rh_plastic) ||
                            !p.control.has_poroelastic_pressure_feedback))
                die(EXIT_CONFIG_VALUE, "Pore recovery requires porosity_model=1, plastic rheology and coupled pressure feedback.");
        }
        for (double exponent : p.mat.permeability_exponent)
            if (!(exponent > 0 && std::isfinite(exponent)))
                die(EXIT_CONFIG_VALUE, "mat.permeability_exponent must be positive and finite.");
        get_numbers(vm, "mat.matrix_storage", p.mat.matrix_storage, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.matrix_exchange", p.mat.matrix_exchange, p.mat.nmat, nmat_key, true);
        for (int m = 0; m < p.mat.nmat; ++m) {
            const double c = p.mat.matrix_storage[m], h = p.mat.matrix_exchange[m];
            if (!(std::isfinite(c) && c >= 0 && std::isfinite(h) && h >= 0) ||
                (h > 0 && !(c > 0 && p.control.has_hydraulic_diffusion)))
                die(EXIT_CONFIG_VALUE, "Matrix exchange requires positive storage and active hydraulic transport.");
        }
        get_numbers(vm, "mat.permeability_tensor", p.mat.permeability_tensor, 6, "six tensor components");
        {
            const double_vec &k = p.mat.permeability_tensor;
            const double det = k[0]*k[1]*k[2] + 2*k[3]*k[4]*k[5]
                             - k[0]*k[5]*k[5] - k[1]*k[4]*k[4] - k[2]*k[3]*k[3];
            if (!std::all_of(k.begin(), k.end(), [](double x) { return std::isfinite(x); }) ||
                !(k[0] > 0 && k[0]*k[1] > k[3]*k[3] && det > 0))
                die(EXIT_CONFIG_VALUE, "mat.permeability_tensor must be symmetric positive definite.");
        }
        get_numbers(vm, "mat.fluid_rho0", p.mat.fluid_rho0, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.fluid_heat_capacity", p.mat.fluid_heat_capacity, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.fluid_alpha", p.mat.fluid_alpha, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.fluid_bulk_modulus", p.mat.fluid_bulk_modulus, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.fluid_visc", p.mat.fluid_visc, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.biot_coeff", p.mat.biot_coeff, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.bulk_modulus_s", p.mat.bulk_modulus_s, p.mat.nmat, nmat_key, true);
        if (p.mat.derive_biot_coeff_from_bulk_moduli)
            for (int m = 0; m < p.mat.nmat; ++m) {
                const double kd = p.mat.bulk_modulus[m], ks = p.mat.bulk_modulus_s[m];
                if (!(kd > 0 && std::isfinite(kd) && ks > 0 && std::isfinite(ks) && kd <= ks))
                    die(EXIT_CONFIG_VALUE, "Derived Biot coefficient requires 0 < bulk_modulus <= bulk_modulus_s, both finite.");
                p.mat.biot_coeff[m] = 1-kd/ks;
            }
        for (int m = 0; m < p.mat.nmat; ++m)
            if (p.control.has_thermal_pressurization &&
                !(std::isfinite(p.mat.alpha[m]) && std::isfinite(p.mat.fluid_alpha[m])))
                die(EXIT_CONFIG_VALUE, "Thermal expansivities must be finite.");
        if (p.control.has_thermal_pressurization && !(p.mat.rheol_type & MatProps::rh_elastic))
            die(EXIT_CONFIG_VALUE, "Thermal pressurization requires a rheology with elastic bulk response.");
        if (p.control.has_fluid_thermal_advection) {
            if (!p.control.has_hydraulic_diffusion)
                die(EXIT_CONFIG_VALUE, "Fluid thermal advection requires has_hydraulic_diffusion=true.");
            for (int m = 0; m < p.mat.nmat; ++m) {
                const double phi = p.mat.porosity[m], cs = p.mat.heat_capacity[m];
                const double cf = p.mat.fluid_heat_capacity[m];
                const double capacity = (1-phi)*p.mat.rho0[m]*cs + phi*p.mat.fluid_rho0[m]*cf;
                if (!(cs >= 0 && std::isfinite(cs) && cf > 0 && std::isfinite(cf) &&
                      capacity > 0 && std::isfinite(capacity) &&
                      p.mat.therm_cond[m] >= 0 && std::isfinite(p.mat.therm_cond[m])))
                    die(EXIT_CONFIG_VALUE, "Fluid heat transport requires nonnegative solid capacity/conductivity, positive fluid capacity and bulk storage.");
            }
        }
        // Buoyancy follows the material: rho_f0 is the density at bc.surface_temperature.
        if (p.control.has_hydraulic_diffusion &&
            (p.control.has_thermal_diffusion || p.control.has_fluid_thermal_advection))
            for (int m = 0; m < p.mat.nmat; ++m)
                if (!std::isfinite(p.mat.fluid_alpha[m]))
                    die(EXIT_CONFIG_VALUE, "mat.fluid_alpha must be finite with hydraulic and thermal transport.");
        if (p.heat_injection.enabled && !p.control.has_fluid_thermal_advection)
            for (int m=0; m<p.mat.nmat; ++m)
                if (!(p.mat.rho0[m]>0 && std::isfinite(p.mat.rho0[m]) &&
                      p.mat.heat_capacity[m]>0 && std::isfinite(p.mat.heat_capacity[m]) &&
                      p.mat.therm_cond[m]>=0 && std::isfinite(p.mat.therm_cond[m])))
                    die(EXIT_CONFIG_VALUE, "Heat sources require positive finite solid density/capacity and nonnegative conductivity.");
        if (p.control.has_thermal_pressurization &&
            !(p.control.has_hydraulic_diffusion && p.control.has_thermal_diffusion))
            die(EXIT_CONFIG_VALUE, "Thermal pressurization requires hydraulic and thermal diffusion.");
        if (p.mat.porosity_model < 0 || p.mat.porosity_model > 1 ||
            p.mat.permeability_model < 0 || p.mat.permeability_model > 2)
            die(EXIT_CONFIG_VALUE, "mat.porosity_model must be 0 or 1; mat.permeability_model must be 0, 1 or 2.");
        if (p.control.hydraulic_content_model < 0 || p.control.hydraulic_content_model > 1)
            die(EXIT_CONFIG_VALUE,"control.hydraulic_content_model must be 0 or 1.");
        if (p.control.hydraulic_content_model == 1) {
            if (!(p.control.has_hydraulic_diffusion && p.control.has_poroelastic_pressure_feedback &&
                  p.control.is_quasi_static && p.mat.porosity_model == 1 &&
                  (p.mat.rheol_type & MatProps::rh_elastic)))
                die(EXIT_CONFIG_VALUE,"Finite-density content requires coupled quasistatic hydraulics, porosity_model=1 and elastic bulk response.");
            if (has_temperature_transport(p) &&
                !(p.control.has_thermal_pressurization && p.control.has_fluid_thermal_advection))
                die(EXIT_CONFIG_VALUE,"Finite-density heat transport requires thermal pressurization and fluid thermal advection.");
            if (p.control.surface_process_option != 0 ||
                (p.mat.nmat > 1 && p.mat.phase_change_option != 0) ||
                std::any_of(p.mat.matrix_exchange.begin(),p.mat.matrix_exchange.end(),[](double h) { return h > 0; }))
                die(EXIT_CONFIG_VALUE,"Finite-density content does not yet support surface processes, phase conversion or matrix exchange.");
            for (int m = 0; m < p.mat.nmat; ++m) {
                if (!(p.mat.rho0[m] > 0 && std::isfinite(p.mat.rho0[m]) &&
                      p.mat.alpha[m] >= 0 && std::isfinite(p.mat.alpha[m]) &&
                      p.mat.fluid_alpha[m] >= 0 && std::isfinite(p.mat.fluid_alpha[m]) &&
                      p.mat.heat_capacity[m] >= 0 && std::isfinite(p.mat.heat_capacity[m]) &&
                      p.mat.fluid_heat_capacity[m] > 0 && std::isfinite(p.mat.fluid_heat_capacity[m])))
                    die(EXIT_CONFIG_VALUE,"Finite-density content requires positive solid density/fluid heat capacity and nonnegative solid heat capacity.");
                if (p.mat.fluid_rho0[m] != p.mat.fluid_rho0[0] ||
                    p.mat.fluid_bulk_modulus[m] != p.mat.fluid_bulk_modulus[0] ||
                    p.mat.fluid_alpha[m] != p.mat.fluid_alpha[0] ||
                    p.mat.fluid_heat_capacity[m] != p.mat.fluid_heat_capacity[0] ||
                    p.mat.fluid_visc[m] != p.mat.fluid_visc[0])
                    die(EXIT_CONFIG_VALUE,"Finite-density content requires one common fluid across solid materials.");
            }
        }
        if ((p.mat.porosity_model || p.mat.permeability_model) && !p.control.has_hydraulic_diffusion)
            die(EXIT_CONFIG_VALUE, "Evolving hydraulic properties require has_hydraulic_diffusion=true.");
        if (p.control.has_hydraulic_diffusion) {
            const auto positive = [](double x) { return x > 0 && std::isfinite(x); };
            for (int m = 0; m < p.mat.nmat; ++m) {
                const double phi = p.mat.porosity[m];
                const double alpha = p.mat.biot_coeff[m];
                if ((p.mat.porosity_model || p.mat.permeability_model) && !(phi > 0 && phi < 1))
                    die(EXIT_CONFIG_VALUE, "Evolving hydraulic properties require 0 < mat.porosity < 1.");
                if (!(phi >= 0 && phi <= alpha && alpha <= 1) ||
                    !positive(p.mat.bulk_modulus[m]) ||
                    !(p.mat.shear_modulus[m] >= 0 && std::isfinite(p.mat.shear_modulus[m])) ||
                    !positive(p.mat.bulk_modulus_s[m]) ||
                    !positive(p.mat.fluid_bulk_modulus[m]) ||
                    !positive(p.mat.fluid_visc[m]) ||
                    !positive(p.mat.fluid_rho0[m]) ||
                    !(p.mat.hydraulic_perm[m] >= 0 && std::isfinite(p.mat.hydraulic_perm[m])) ||
                    alpha == 0) {
                    std::cerr << "Invalid hydraulic properties for material " << m << ".\n";
                    die(EXIT_CONFIG_VALUE,
                        "Hydraulics requires 0 <= porosity <= biot_coeff <= 1, nonzero "
                        "biot_coeff, positive finite bulk moduli, fluid_visc and fluid_rho0, "
                        "and nonnegative finite shear_modulus and hydraulic_perm.");
                }
            }
        }
        get_numbers(vm, "mat.retention_model", p.mat.retention_model, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.retention_alpha", p.mat.retention_alpha, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.retention_n", p.mat.retention_n, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.residual_saturation", p.mat.residual_saturation, p.mat.nmat, nmat_key, true);
        {
            bool active = false;
            for (int m = 0; m < p.mat.nmat; ++m) {
                const int model = p.mat.retention_model[m];
                if (model < 0 || model > 3)
                    die(EXIT_CONFIG_VALUE, "mat.retention_model must be 0, 1, 2 or 3.");
                if (model == 0) continue;
                active = true;
                const double a = p.mat.retention_alpha[m], n = p.mat.retention_n[m];
                const double sr = p.mat.residual_saturation[m];
                if (!(a > 0 && std::isfinite(a)) || !(sr >= 0 && sr < 1) ||
                    (model == 1 && !(n > 1 && std::isfinite(n))) ||
                    (model == 2 && !(n > 0 && std::isfinite(n))) ||
                    !(p.mat.porosity[m] > 0)) {
                    std::cerr << "Invalid retention properties for material " << m << ".\n";
                    die(EXIT_CONFIG_VALUE,
                        "Retention laws require finite alpha > 0, 0 <= residual_saturation < 1, "
                        "porosity > 0, and retention_n > 1 (van Genuchten) or > 0 (Brooks-Corey).");
                }
            }
            if (active && !p.control.has_hydraulic_diffusion)
                die(EXIT_CONFIG_VALUE, "mat.retention_model > 0 requires has_hydraulic_diffusion=true.");
            // Each of these changes the pore volume or exchanges fluid outside
            // the unsaturated content balance, which is not accounted yet.
            const bool exchange = std::any_of(p.mat.matrix_exchange.begin(), p.mat.matrix_exchange.end(),
                                              [](double h) { return h > 0; });
            if (active && (p.mat.porosity_model || p.mat.permeability_model ||
                           p.control.has_thermal_pressurization || p.control.has_fluid_thermal_advection || exchange))
                die(EXIT_CONFIG_VALUE, "mat.retention_model > 0 does not yet support evolving "
                    "porosity/permeability, thermal coupling or matrix exchange.");
        }
        // Rate-and-state friction parameters
        get_numbers(vm, "mat.direct_a", p.mat.direct_a, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.evolution_b", p.mat.evolution_b, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.characteristic_velocity", p.mat.characteristic_velocity, p.mat.nmat, nmat_key, true);
        get_numbers(vm, "mat.characteristic_distance", p.mat.characteristic_distance, p.mat.nmat, nmat_key, true);
        if (p.mat.state_var_model < 0 || p.mat.state_var_model > 2) {
            die(EXIT_CONFIG_VALUE, "mat.state_var_model must be 0, 1, or 2.");
        }
        if (p.control.has_PT) {
            if (p.control.PT_use_running_scale &&
                (!vm.defaulted("control.PT_relative_tolerance") ||
                 !vm.defaulted("control.PT_absolute_tolerance")))
                std::cerr << "Warning: PT_use_running_scale ignores explicitly configured "
                          << "PT_relative_tolerance/PT_absolute_tolerance for physical steps; "
                          << "they still apply to initial equilibrium. Physical steps use "
                          << "epsilon=1e-6 and C=1000.\n";
            if (p.control.PT_max_iter < 1 || p.control.PT_stagnation_window < 0 ||
                p.control.PT_option < 0 || p.control.PT_option > 1 ||
                !(p.control.PT_CFL > 0 && p.control.PT_CFL <= 0.5) ||
                !(p.control.PT_Re > 0 && std::isfinite(p.control.PT_Re)) ||
                !(p.control.PT_relative_tolerance >= 0 && std::isfinite(p.control.PT_relative_tolerance)) ||
                !(p.control.PT_absolute_tolerance >= 0 && std::isfinite(p.control.PT_absolute_tolerance)) ||
                (p.control.PT_relative_tolerance == 0 && p.control.PT_absolute_tolerance == 0))
                die(EXIT_CONFIG_VALUE, "Invalid PT iteration, tolerance or relaxation parameters.");
        }
        if (p.control.rsf_dtheta_max > 0) {
            if (p.control.fixed_dt != 0.0) {
                die(EXIT_CONFIG_VALUE,
                    "control.rsf_dtheta_max requires control.fixed_dt=0.");
            }
            if (p.control.has_PT) {
                die(EXIT_CONFIG_VALUE,
                    "control.rsf_dtheta_max is not supported with "
                    "control.has_PT=true.");
            }
            if (p.ic.has_body_force_adjustment) {
                die(EXIT_CONFIG_VALUE,
                    "control.rsf_dtheta_max is not supported with "
                    "ic.has_body_force_adjustment=true.");
            }
            if (p.ic.isostasy_adjustment_time_in_yr > 0) {
                die(EXIT_CONFIG_VALUE,
                    "control.rsf_dtheta_max is not supported during "
                    "isostasy adjustment.");
            }
            if (!(p.mat.rheol_type & MatProps::rh_rsf)) {
                die(EXIT_CONFIG_VALUE,
                    "control.rsf_dtheta_max requires an RSF rheology.");
            }
            if (p.mat.state_var_model != 1) {
                die(EXIT_CONFIG_VALUE,
                    "control.rsf_dtheta_max requires the aging law "
                    "(mat.state_var_model=1).");
            }
            if (p.control.rsf_slip_rate_projection_option !=
                rsf_slip_rate_projection_total_strain_rate) {
                die(EXIT_CONFIG_VALUE,
                    "control.rsf_dtheta_max requires "
                    "control.rsf_slip_rate_projection_option=1.");
            }
        }
        if (p.mat.rheol_type & MatProps::rh_rsf) {
            for (int m = 0; m < p.mat.nmat; ++m) {
                if (!std::isfinite(p.mat.characteristic_velocity[m]) ||
                    p.mat.characteristic_velocity[m] <= 0.0) {
                    std::cerr << "Error: mat.characteristic_velocity must be > 0 for RSF materials "
                              << "(material index " << m << ").\n";
                    die(EXIT_CONFIG_VALUE);
                }
                if (!std::isfinite(p.mat.characteristic_distance[m]) ||
                    p.mat.characteristic_distance[m] <= 0.0) {
                    std::cerr << "Error: mat.characteristic_distance must be > 0 for RSF materials "
                              << "(material index " << m << ").\n";
                    die(EXIT_CONFIG_VALUE);
                }
            }
        }
    }

}


void get_input_parameters(const char* filename, Param& p)
{
    opt::Registry cfg;

    declare_parameters(cfg, p);
    // print help message
    if (std::strncmp(filename, "-h", 3) == 0 ||
        std::strncmp(filename, "--help", 7) == 0) {
        std::cout << "Config file options:\n" << cfg.help();
        const char* exe =
#ifdef THREED
            "dynearthsol3d"
#else
            "dynearthsol2d"
#endif
#ifdef ACC
            ".gpu"
#endif
            ;
        // The patterns are split across literals: a whole copy in .rodata would surface as a
        // stray match in the very extraction they describe.
        std::cout << "\nBuild snapshot:\n"
                     "  The build identity is embedded in the executable. A run\n"
                     "  prints it at start ([build][...] lines) unless\n"
                     "  sim.has_runtime_info_display = no; to read it from a binary\n"
                     "  without running it:\n"
                     "    strings " << exe << " | grep '^build"
                  << "\\.snapshot\\.'\n"
                     "  make snapshot_diff=1 also embeds the uncommitted code changes;\n"
                     "  extract them with:\n"
                     "    strings " << exe << " | sed -n '/^build\\.code-changes\\.beg"
                  << "in :$/,/^build\\.code-changes\\.e" << "nd   :$/p'\n";
        // Not die(): that prints an "[DES exit N] <category>" banner meant for
        // failures, and --help is a successful run with nothing to diagnose.
        std::exit(EXIT_OK);
    }
    read_parameters_from_file(filename, cfg);
    validate_parameters(cfg, p);
}
