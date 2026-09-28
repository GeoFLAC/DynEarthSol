#include <cmath>
#include <limits>

#include "pseudo_transient.hpp"
#include "bc.hpp"
#include "fields.hpp"
#include "geometry.hpp"
#include "monitor.hpp"
#include "output.hpp"
#include "remeshing.hpp"
#include "rheology.hpp"

void run_physical_step_pt(Param &param, Variables &var, bool &hydraulic_diffusion_switch)
{
    double residual_old = std::numeric_limits<double>::max();

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
}
