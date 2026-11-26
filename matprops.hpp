#ifndef DYNEARTHSOL3D_MATPROPS_HPP
#define DYNEARTHSOL3D_MATPROPS_HPP

#include <cmath>
#include "parameters.hpp"
#include "unsaturated.hpp"

// Eulerian Biot porosity: eps is the accumulated volume/grain-pressure strain.
#pragma acc routine seq
inline double biot_porosity(double phi0, double alpha, double eps)
{
    return alpha == phi0 ? phi0 : phi0 - (alpha - phi0) * std::expm1(-eps);
}

#pragma acc routine seq
inline double kozeny_carman(double k0, double phi0, double phi)
{
    const double pore = phi / phi0;
    const double solid = (1 - phi0) / (1 - phi);
    return k0 * pore * pore * pore * solid * solid;
}

#pragma acc routine seq
double ref_pressure(const Param& param, double z);


class MatProps
{
public:
    MatProps(const Param& param, const Variables& var);
    ~MatProps();

    const int rheol_type;
    const int nmat;
    const int permeability_model;
    const bool fluid_thermal_advection;
    double perm_tensor[NDIMS][NDIMS]; // constant shape; material magnitude is perm(e)
    double perm_bound; // upper bound on the largest tensor eigenvalue
    double hydro_step_limit = 0; // rebuilt after storage/mesh changes, not checkpointed
    double thermal_step_limit = 0; // rebuilt after capacity/mesh changes, not checkpointed

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
    double cp_fluid(int e) const;
    #pragma acc routine seq
    double thermal_capacity(int e) const;
    #pragma acc routine seq
    double k(int e) const;

    // hydraulic parameter function
    #pragma acc routine seq
    double phi(int e) const;
    #pragma acc routine seq
    double phi0(int e) const;
    #pragma acc routine seq
    double perm(int e) const;
    #pragma acc routine seq
    double alpha_fluid(int e) const;
    #pragma acc routine seq
    double alpha_solid(int e) const;
    #pragma acc routine seq
    double matrix_capacity(int e) const;
    #pragma acc routine seq
    double exchange(int e) const;
    #pragma acc routine seq
    double pore_healing_rate(int e) const;
    #pragma acc routine seq
    double beta_fluid(int e) const;
    // Constant material reference density for Darcy gravity and mixture weight.
    #pragma acc routine seq
    double rho_fluid_reference(int e) const;
    #pragma acc routine seq
    double mu_fluid(int e) const;
    #pragma acc routine seq
    double alpha_biot(int e) const;
    #pragma acc routine seq
    double beta_mineral(int e) const;
    // Retention law of material m; has_retention when any law is unsaturated.
    #pragma acc routine seq
    const Retention &retention(int m) const { return retention_laws[m]; }
    const bool has_retention;

    // rate-and-state friction parameters
    #pragma acc routine seq
    double d_a(int e) const;
    #pragma acc routine seq
    double e_b(int e) const;
    #pragma acc routine seq
    double c_v(int e) const;
    #pragma acc routine seq
    double d_c(int e) const;

    #pragma acc routine seq
    void plastic_props(int e, double pls,
                       double& amc, double& anphi, double& anpsi,
                       double& hardn, double& ten_max) const;
    // Factors of the strength that depend only on the friction and dilation
    // angles, so that a strength with fixed angles needs no trigonometry.
    struct AngleFactors {
        double anphi, anpsi, sqrt_anphi, tan_phi;
        bool zero_friction;
    };
    // One linear piece of an element's weakening law, starting at plastic
    // strain pls: the cohesion and the friction and dilation angles at pls,
    // their change per unit plastic strain, and the plastic strain where the
    // piece ends (the next pls0 or pls1 of the element's materials, or the
    // largest double). A rate-and-state piece also keeps this step's rate and
    // state terms, which stay fixed along it.
    struct StrengthPiece {
        double pls, end;
        double cohesion, friction, dilation;
        double dcohesion, dfriction, ddilation;
        int rsf;              // 0: none; 1: rate and state; 2: frozen friction
        double rsf_terms[2];  // added to tan(friction), in this order
        AngleFactors factors;               // at pls
        double amc, anphi, anpsi, ten_max;  // strength at pls
        bool has_slopes;                    // dfriction, ddilation and end are set

        // Whether the friction or dilation angle changes along the piece.
        #pragma acc routine seq
        bool angles_change() const
        {
            return (rsf != 2 && dfriction != 0) || ddilation != 0;
        }
    };
    #pragma acc routine seq
    void strength_piece(int e, double pls, StrengthPiece& piece) const;
    // The strength at pls only, for every element; complete_piece adds the
    // angle slopes and the end once the element yields.
    #pragma acc routine seq
    void strength_here(int e, double pls, StrengthPiece& piece) const;
    #pragma acc routine seq
    void complete_piece(int e, StrengthPiece& piece) const;
    #pragma acc routine seq
    void plastic_props_rsf(int e, double pls,
                       double& amc, double& anphi, double& anpsi,
                       double& hardn, double& ten_max, double& slip_rate,
                       double& dyn_fric_coeff, double& state_variable,
                       double dt, int state_model, bool frozen_history = false) const;
    // The rate-and-state piece at pls with this step's rate and state held
    // fixed: no state update and no change to the stored friction coefficient.
    #pragma acc routine seq
    void rsf_strength_piece(int e, double pls, double slip_rate, double state_variable,
                            double dyn_fric_coeff, int state_model, bool frozen_history,
                            StrengthPiece& piece) const;
    // Strength at plastic strain piece.pls + offset, within the piece.
    #pragma acc routine seq
    void piece_strength(const StrengthPiece& piece, double offset,
                        double& amc, double& anphi, double& anpsi, double& ten_max) const;
    // The strength of a piece over the offsets from x to y, for a cohesion that is not
    // negative: the lowest and highest amc, anphi and ten_max, and bounds of the size of
    // their second derivatives per unit plastic strain (for ten_max, of cohesion/tan_phi
    // before the cap at tension_max; -1 where not known) and of the error of the lookup
    // tables in them.
    struct StrengthRange {
        double amc[2], anphi[2], ten_max[2];
        double curvature_amc, curvature_anphi, curvature_ten;
        double error_amc, error_anphi, error_ten;
    };
    #pragma acc routine seq
    void piece_strength_range(const StrengthPiece& piece, double x, double y,
                              StrengthRange& range) const;
    #pragma acc routine seq
    void rsf_friction_from_state(int e, double pls, double slip_rate,
                                 double state_variable, double& dyn_fric_coeff,
                                 int state_model) const;

    // Rebuild the per-element property means if elemmarkers changed since the last call.
    // Must be called after ANYTHING that mutates elemmarkers, before the next read.
    void refresh_elem_cache();

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

    double_vec rho0, alpha;
    double_vec bulk_modulus, shear_modulus;
    double_vec visc_exponent, visc_coefficient, visc_activation_energy;
    double_vec visc_activation_volume;

    // Material-only terms of the Chen & Morgan creep law, built by the constructor.
    // Pure function of visc_exponent / visc_coefficient, never mutated.
    double_vec visc_pow_edot, visc_coef_term, visc_nR;

    // The hot means, rebuilt by refresh_elem_cache(): pure functions of elemmarkers.
    // Not exhaustive -- rho/visc CANNOT be cached (T, stress, edot). Stride 8 keeps one
    // element's slots in a 64 B line; seven separate arrays measured 1.6x slower.
    static const int MP_STRIDE = 8;
    enum { MPC_BULKM = 0, MPC_SHEARM, MPC_PHI, MPC_ALPHA_BIOT,
           MPC_CP, MPC_K, MPC_BETA_FLUID };
    double_vec cache_props;

    // Flat copy of elemmarkers, kept ONLY to detect that it moved. A version counter
    // bumped by each mutation site would go stale the first time one is added.
    int_vec cache_markers;
    double_vec heat_capacity, therm_cond;
    double_vec fluid_heat_capacity;
    double_vec pls0, pls1;
    double_vec cohesion0, cohesion1;
    double_vec friction_angle0, friction_angle1;
    double_vec dilation_angle0, dilation_angle1;
    // Change of cohesion and friction and dilation angles per unit plastic
    // strain on each material's weakening segment, built by the constructor.
    double_vec cohesion_slope, friction_slope, dilation_slope;

    const double_vec &log_table;
    const double_vec &tan_table;
    const double_vec &sin_table;

    double_vec matrix_storage, matrix_exchange;
    double_vec porosity, hydraulic_perm, fluid_rho0;
    double_vec permeability_exponent;
    double_vec pore_healing_time, pore_healing_temperature, pore_healing_energy;
    const double_vec *pore_strain;
    double_vec fluid_alpha, fluid_bulk_modulus, fluid_visc;
    double_vec biot_coeff, bulk_modulus_s;
    std::vector<Retention> retention_laws;

    // rate-and-state friction
    double_vec direct_a, evolution_b, characteristic_velocity, characteristic_distance;
    // double_vec static_friction_coefficient;

    // Marker means of cohesion, friction and dilation angles and the cohesion
    // slope at pls over n markers.
    #pragma acc routine seq
    void weakening_means(int e, double pls,
                         double &cohesion, double &friction_angle,
                         double &dilation_angle, double &hardening, int &n) const;

    #pragma acc routine seq
    void plastic_weakening(int e, double pls,
                           double &cohesion, double &friction_angle,
                           double &dilation_angle, double &hardening) const;

    #pragma acc routine seq
    void plastic_weakening_rsf(int e, double pls,
                           double &cohesion, double &dynamic_friction_angle,
                           double &dilation_angle, double &hardening,
                           double &slip_rate, double& dyn_fric_coeff,
                           double& state_variable, int state_model, double dt, bool frozen_history = false) const;
    #pragma acc routine seq
    void update_state_variable(int e, double slip_rate, double& state_variable,
                               double dt, int state_model) const;
    // The rate and state terms added to tan(static friction).
    #pragma acc routine seq
    void rsf_terms(int e, double slip_rate, double state_variable, int state_model,
                   double terms[2]) const;
    // rsf selects the library trigonometry of rate-and-state instead of the
    // tables of Mohr-Coulomb.
    #pragma acc routine seq
    void angle_factors(double phi, double psi, bool rsf, AngleFactors& factors) const;
    // The cohesion and angle factors at plastic strain piece.pls + offset.
    #pragma acc routine seq
    void piece_point(const StrengthPiece& piece, double offset,
                     double& cohesion, AngleFactors& factors) const;
    #pragma acc routine seq
    void strength_from(double cohesion, const AngleFactors& factors,
                       double& amc, double& anphi, double& anpsi, double& ten_max) const;
    #pragma acc routine seq
    void finish_piece(StrengthPiece& piece) const;
    #pragma acc routine seq
    void ep_derived(double cohesion, double phi, double psi,
                    double& amc, double& anphi, double& anpsi, double& ten_max) const;
    #pragma acc routine seq
    void rsf_derived(double cohesion, double phi, double psi,
                     double& amc, double& anphi, double& anpsi, double& ten_max) const;
};


#endif
