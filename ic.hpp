#ifndef DYNEARTHSOL3D_IC_HPP
#define DYNEARTHSOL3D_IC_HPP

double get_prem_pressure(double depth);
double get_prem_pressure_modified(double depth);
void initial_stress_state(const Param &param, const Variables &var,
                          tensor_t &stress, double_vec &stressyy, double_vec &old_mean_stress, tensor_t &strain,
                          double &compensation_pressure);
void initial_stress_state_1d_load(const Param &param, const Variables &var,
                          tensor_t &stress, double_vec &stressyy, double_vec &old_mean_stress, tensor_t &strain,
                          double &compensation_pressure);
void initial_weak_zone(const Param &param, const Variables &var,
                       double_vec &plstrain);
void initial_temperature(const Param &param, const Variables &var,
                         double_vec &temperature, double_vec &radiogenic_source, double &bottom_temperature);
void initial_hydrostatic_state(const Param &param, const Variables &var,
                          double_vec &ppressure, double_vec &dppressure);

#endif
