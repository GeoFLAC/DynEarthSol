#ifndef DYNEARTHSOL3D_MATPROPS_HPP
#define DYNEARTHSOL3D_MATPROPS_HPP

#include "parameters.hpp"


#pragma acc routine seq
double ref_pressure(const Param& param, double z);


class MatProps
{
public:
    MatProps(const Param& param, const Variables& var);
    ~MatProps();

    const int rheol_type;
    const int nmat;

    #pragma acc routine seq
    double bulkm(int e) const;
    #pragma acc routine seq
    double shearm(int e) const;
    #pragma acc routine seq
    double visc(int e) const;

    #pragma acc routine seq
    double rho(int e) const;
    #pragma acc routine seq
    double cp(int e) const;
    #pragma acc routine seq
    double k(int e) const;

    // hydraulic parameter function
    #pragma acc routine seq
    double phi(int e) const;
    #pragma acc routine seq
    double perm(int e) const;
    #pragma acc routine seq
    double alpha_fluid(int e) const;
    #pragma acc routine seq
    double beta_fluid(int e) const;
    #pragma acc routine seq
    double rho_fluid(int e) const;
    #pragma acc routine seq
    double mu_fluid(int e) const;
    #pragma acc routine seq
    double alpha_biot(int e) const;
    #pragma acc routine seq
    double beta_mineral(int e) const;

    // rate-and-state friction parameters
    #pragma acc routine seq
    double d_a(int e) const;
    #pragma acc routine seq
    double e_b(int e) const;
    #pragma acc routine seq
    double c_v(int e) const;

    #pragma acc routine seq
    void plastic_props(int e, double pls,
                       double& amc, double& anphi, double& anpsi,
                       double& hardn, double& ten_max) const;
    #pragma acc routine seq
    void plastic_props_rsf(int e, double pls,
                       double& amc, double& anphi, double& anpsi,
                       double& hardn, double& ten_max, double& slip_rate) const;

    const bool is_plane_strain;
    const double visc_min;
    const double visc_max;
    const double tension_max;
    const double therm_diff_max;
    double hydro_diff_max;

    const static int rh_elastic = 1 << 0; // Decimal value 1
    const static int rh_viscous = 1 << 1; // Decimal value 2
    const static int rh_plastic = 1 << 2; // Decimal value 4
    const static int rh_plastic2d = rh_plastic | 1 << 3; // Decimal value 12
    const static int rh_rsf = 1 << 4; // Decimal value 16
    const static int rh_maxwell = rh_elastic | rh_viscous; // Decimal value 3
    const static int rh_ep = rh_elastic | rh_plastic; // Decimal value 5
    const static int rh_evp = rh_elastic | rh_viscous | rh_plastic; // Decimal value 7
    const static int rh_ep_rsf = rh_elastic | rh_plastic | rh_rsf;  // Decimal value 21
    const static int rh_evp_rsf = rh_elastic | rh_viscous | rh_plastic | rh_rsf; // Decimal value 23

private:

    // alias to field variables in var
    // ie. var.mat.temperature == var.temperature
    const array_t &coord;
    const conn_t &connectivity;
    const double_vec &temperature;
    const tensor_t &stress;
    const tensor_t &strain_rate;
    const int_vec2D &elemmarkers;
    const double_vec &log_table;
    const double_vec &tan_table;
    const double_vec &sin_table;

    double_vec rho0, alpha;
    double_vec bulk_modulus, shear_modulus;
    double_vec visc_exponent, visc_coefficient, visc_activation_energy;
    double_vec visc_activation_volume;
    double_vec heat_capacity, therm_cond;
    double_vec pls0, pls1;
    double_vec cohesion0, cohesion1;
    double_vec friction_angle0, friction_angle1;
    double_vec dilation_angle0, dilation_angle1;

    // hydraulic process
    const double_vec &ppressure;
    const double_vec &dppressure;
    double_vec porosity, hydraulic_perm, fluid_rho0;
    double_vec fluid_alpha, fluid_bulk_modulus, fluid_visc;
    double_vec biot_coeff, bulk_modulus_s;

    // rate-and-state friction
    double_vec direct_a, evolution_b, characteristic_velocity;
    // double_vec static_friction_coefficient;

    #pragma acc routine seq
    void plastic_weakening(int e, double pls,
                           double &cohesion, double &friction_angle,
                           double &dilation_angle, double &hardening) const;

    #pragma acc routine seq
    void plastic_weakening_rsf(int e, double pls,
                           double &cohesion, double &friction_angle,
                           double &dilation_angle, double &hardening, double &slip_rate) const;
};


#endif