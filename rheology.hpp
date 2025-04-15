#ifndef DYNEARTHSOL3D_RHEOLOGY_HPP
#define DYNEARTHSOL3D_RHEOLOGY_HPP

void update_stress(const Param& param , Variables& var, tensor_t& stress, double_vec& stressyy,
                   double_vec& dpressure, double_vec& viscosity, tensor_t& strain, double_vec& plstrain,
                   double_vec& delta_plstrain, tensor_t& strain_rate,
                   double_vec& ppressure, double_vec& dppressure, array_t& vel);

void update_old_mean_stress(const Param& param ,const Variables& var, tensor_t& stress, double_vec& old_mean_stress);
#endif
