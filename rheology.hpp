#ifndef DYNEARTHSOL3D_RHEOLOGY_HPP
#define DYNEARTHSOL3D_RHEOLOGY_HPP

void update_stress(const Param& param , Variables& var, tensor_t& stress, double_vec& stressyy,
                   double_vec& dpressure, double_vec& viscosity, tensor_t& strain, double_vec& plstrain,
                   double_vec& delta_plstrain, tensor_t& strain_rate,
                   double_vec& ppressure, double_vec& dppressure, array_t& vel);

void update_old_mean_stress(const Param& param ,const Variables& var, tensor_t& stress, double_vec& old_mean_stress);

// #ifdef RS
// void friction_variables(double &T, double &direct_a, double &evolution_b, double &characteristic_velocity, double &static_friction_coefficient);
// void update_state1(const Variables &var, double_vec &state1, double_vec &slip_velocity, int a);
// #endif

// void update_stress_old(const Param& param ,const Variables& var, tensor_t& stress, tensor_t& stress_old);

// void add_stress_old(const Param& param ,const Variables& var, tensor_t& stress, tensor_t& stress_old);

#endif
