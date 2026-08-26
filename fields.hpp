#ifndef DYNEARTHSOL3D_FIELDS_HPP
#define DYNEARTHSOL3D_FIELDS_HPP

#include "parameters.hpp"

#ifdef ACC
#pragma acc routine seq
#endif
inline bool is_fixed_pore_pressure_node(int node, const Variables &var)
{
    const uint flag = (*var.bcflag)[node];
    return ((flag & BOUNDX0) && var.hbc_types[0] == 1) ||
           ((flag & BOUNDX1) && var.hbc_types[1] == 1) ||
           ((flag & BOUNDY0) && var.hbc_types[2] == 1) ||
           ((flag & BOUNDY1) && var.hbc_types[3] == 1) ||
           ((flag & BOUNDZ0) && var.hbc_types[4] == 1) ||
           ((flag & BOUNDZ1) && var.hbc_types[5] == 1);
}

void enforce_pore_pressure_bcs(const Param &param, const Variables &var,
                               double_vec &ppressure,
                               double_vec &dppressure);
void allocate_variables(const Param &param, Variables& var);
void reallocate_tmp(const Param &param, Variables& var);
void reallocate_variables(const Param &param, Variables& var);
void update_temperature(const Param &param, const Variables &var,
    double_vec &temperature, elem_cache& tmp_result);
void update_pore_pressure(const Param &param, const Variables &var,
    double_vec &ppressure, double_vec &dppressure, double_vec &tdot, elem_cache& tmp_result, tensor_t &stress, double_vec& old_mean_stress);
void update_strain_rate(const Variables& var, tensor_t& strain_rate);
void update_force(const Param& param, const Variables& var, array_t& force, array_t& force_residual, 
    elem_cache& tmp_result);
double calculate_residual_force(const Variables& var, array_t& vel);
void update_velocity(const Variables& var, array_t& vel);
// void update_velocity_PT(const Variables& var, array_t& vel);
void update_coordinate(const Variables& var, array_t& coord);
void rotate_stress(const Variables &var, tensor_t &stress, tensor_t &strain);

#endif
