#ifndef DYNEARTHSOL3D_FIELDS_HPP
#define DYNEARTHSOL3D_FIELDS_HPP

// Element shape gradients shared by field updates and the current-rate timestep bound.
#ifdef THREED
#pragma acc routine seq
inline void get_local_shape_fn(const Variables &var, int e, double shpdx[4], double shpdy[4], double shpdz[4])
{
    ConstArrayIndirectAccessor d = var.coord->view_const((*var.connectivity)[e]);

    double iv = 1.0 / (6.0 * (*var.volume)[e]);

    double x01 = d[0][0] - d[1][0]; double x02 = d[0][0] - d[2][0]; double x03 = d[0][0] - d[3][0];
    double x12 = d[1][0] - d[2][0]; double x13 = d[1][0] - d[3][0]; double x23 = d[2][0] - d[3][0];
    double y01 = d[0][1] - d[1][1]; double y02 = d[0][1] - d[2][1]; double y03 = d[0][1] - d[3][1];
    double y12 = d[1][1] - d[2][1]; double y13 = d[1][1] - d[3][1]; double y23 = d[2][1] - d[3][1];
    double z01 = d[0][2] - d[1][2]; double z02 = d[0][2] - d[2][2]; double z03 = d[0][2] - d[3][2];
    double z12 = d[1][2] - d[2][2]; double z13 = d[1][2] - d[3][2]; double z23 = d[2][2] - d[3][2];

    shpdx[0] = iv * (y13*z12 - y12*z13);
    shpdx[1] = iv * (y02*z23 - y23*z02);
    shpdx[2] = iv * (y13*z03 - y03*z13);
    shpdx[3] = iv * (y01*z02 - y02*z01);

    shpdy[0] = iv * (z13*x12 - z12*x13);
    shpdy[1] = iv * (z02*x23 - z23*x02);
    shpdy[2] = iv * (z13*x03 - z03*x13);
    shpdy[3] = iv * (z01*x02 - z02*x01);

    shpdz[0] = iv * (x13*y12 - x12*y13);
    shpdz[1] = iv * (x02*y23 - x23*y02);
    shpdz[2] = iv * (x13*y03 - x03*y13);
    shpdz[3] = iv * (x01*y02 - x02*y01);
}
#else
#pragma acc routine seq
inline void get_local_shape_fn(const Variables &var, int e, double shpdx[3], double shpdz[3])
{
    ConstArrayIndirectAccessor d = var.coord->view_const((*var.connectivity)[e]);

    double iv = 1.0 / (2.0 * (*var.volume)[e]);

    shpdx[0] = iv * (d[1][1] - d[2][1]);
    shpdx[1] = iv * (d[2][1] - d[0][1]);
    shpdx[2] = iv * (d[0][1] - d[1][1]);

    shpdz[0] = iv * (d[2][0] - d[1][0]);
    shpdz[1] = iv * (d[0][0] - d[2][0]);
    shpdz[2] = iv * (d[1][0] - d[0][0]);
}
#endif

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
void rotate_stress(const Variables &var, tensor_t &stress, tensor_t &strain,
                   double dt);

#endif
