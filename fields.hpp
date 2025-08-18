#ifndef DYNEARTHSOL3D_FIELDS_HPP
#define DYNEARTHSOL3D_FIELDS_HPP

void allocate_variables(const Param &param, Variables& var);
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
void update_velocity_PT(const Variables& var, array_t& vel);
void update_coordinate(const Variables& var, array_t& coord);
void rotate_stress(const Variables &var, tensor_t &stress, tensor_t &strain);

#endif
