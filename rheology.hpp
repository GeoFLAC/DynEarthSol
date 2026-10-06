#ifndef DYNEARTHSOL3D_RHEOLOGY_HPP
#define DYNEARTHSOL3D_RHEOLOGY_HPP

void update_stress(const Param& param , Variables& var, tensor_t& stress, double_vec& stressyy,
                   double_vec& dpressure, double_vec& viscosity, tensor_t& strain, double_vec& plstrain,
                   double_vec& delta_plstrain, tensor_t& strain_rate,
                   double_vec& ppressure, double_vec& dppressure, array_t& vel,
                   double_vec& dyn_fric_coeff, double_vec& state_variable,
                   bool trial_reference_geometry = false, bool initial_equilibrium = false);
void refresh_rsf_friction(const Param& param, Variables& var,
                          double_vec& dyn_fric_coeff,
                          const double_vec& state_variable);

void update_old_mean_stress(const Param& param ,const Variables& var, tensor_t& stress, double_vec& old_mean_stress);
// One physical-start baseline; owned by the mechanical solve, never remapped
// while trials are live. Restore in place to preserve MatProps aliases.
struct MechanicalState {
    int nnode = 0;
    double physical_dt = 0, physical_time = 0;
    const conn_t* connectivity = nullptr;
    tensor_t stress, strain;
    Array2D<double, 7> scalars;
    explicit MechanicalState(int nelem) : stress(nelem), strain(nelem), scalars(nelem) {}
    void capture(const Variables& var);
    void restore(Variables& var) const;
};
#endif
