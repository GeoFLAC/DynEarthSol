#ifndef DYNEARTHSOL3D_FIELDS_HPP
#define DYNEARTHSOL3D_FIELDS_HPP

// Element shape gradients shared by field updates and the current-rate timestep bound.
#ifdef THREED
#pragma acc routine seq
inline void get_local_shape_fn(const Variables &var, int e,
                               double shpdx[NODES_PER_ELEM], double shpdy[NODES_PER_ELEM],
                               double shpdz[NODES_PER_ELEM])
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
inline void get_local_shape_fn(const Variables &var, int e,
                               double shpdx[NODES_PER_ELEM], double shpdz[NODES_PER_ELEM])
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

// The same gradients as rows g[d][j] = dN_j/dx_d, for loops over the dimensions.
#pragma acc routine seq
inline void get_shape_gradients(const Variables &var, int e, double g[NDIMS][NODES_PER_ELEM])
{
#ifdef THREED
    get_local_shape_fn(var, e, g[0], g[1], g[2]);
#else
    get_local_shape_fn(var, e, g[0], g[1]);
#endif
}

// Number of explicit substeps of at most limit that cover dt; dies with
// message when it does not fit in an int.
int explicit_substeps(double dt, double limit, const char *message);
// Smallest mass/row over the owned nodes of all ranks, starting from limit.
// A node's row sums its element rows in work.
double nodal_step_limit(const Variables &var, const elem_cache &work, const double_vec &mass,
                        double limit);

void allocate_variables(const Param &param, Variables& var);
void reallocate_tmp(const Param &param, Variables& var);
void reallocate_variables(const Param &param, Variables& var);
void update_temperature(const Param &param, const Variables &var,
    double_vec &temperature, elem_cache& tmp_result);
void update_strain_rate(const Variables& var, tensor_t& strain_rate);
// force_undamped optionally receives a caller-owned nnode-sized copy of the force with
// the physical loads but without artificial damping. It must not alias force,
// force_residual or assembly inputs. damp_force=false skips artificial damping.
// node_force_magnitude optionally receives, at each node, the sum of the
// absolute element force contributions, before boundary assembly reuses tmp_result.
void update_force(const Param& param, const Variables& var, array_t& force, array_t& force_residual,
    elem_cache& tmp_result, array_t* force_undamped = nullptr, double displacement_dt = 0, bool damp_force = true,
    array_t* node_force_magnitude = nullptr);
// Artificial damping only; OpenACC work is queued asynchronously.
void apply_damping(const Param& param, const Variables& var, array_t& force);
double calculate_residual_force(const Variables& var, array_t& force_residual);
void update_velocity(const Variables& var, array_t& vel);
void update_coordinate(const Variables& var, array_t& coord);
void rotate_stress(const Variables &var, tensor_t &stress, tensor_t &strain,
                   double dt);

#endif
