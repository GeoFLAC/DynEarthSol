#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "monitor.hpp"
#include "geometry.hpp"
#include "matprops.hpp"

namespace {

struct MonitorPointState {
    int id = -1;
    double query_coord_initial[NDIMS] = {0.0};
    double query_coord_rebind[NDIMS] = {0.0};
    int node_id = -1;
    int elem_id = -1;
    double node_dist = std::numeric_limits<double>::infinity();
    double elem_dist = std::numeric_limits<double>::infinity();
    std::FILE* fp = nullptr;
};

struct MonitorManager {
    bool initialized = false;
    int last_written_step = std::numeric_limits<int>::min();
    std::vector<MonitorPointState> points;
};

MonitorManager g_monitor;

inline const char* axis_name(const int d)
{
#ifdef THREED
    static const char* kAxis[3] = {"x", "y", "z"};
#else
    static const char* kAxis[2] = {"x", "z"};
#endif
    return kAxis[d];
}

inline void csv_write_sep(std::FILE* fp, bool& first)
{
    if (!first) std::fputc(',', fp);
    first = false;
}

inline void csv_write_name(std::FILE* fp, const char* name, bool& first)
{
    csv_write_sep(fp, first);
    std::fputs(name, fp);
}

inline void csv_write_int(std::FILE* fp, int v, bool& first)
{
    csv_write_sep(fp, first);
    std::fprintf(fp, "%d", v);
}

inline void csv_write_double(std::FILE* fp, double v, bool& first)
{
    csv_write_sep(fp, first);
    if (std::isnan(v)) std::fputs("nan", fp);
    else if (std::isinf(v)) std::fputs((v > 0) ? "inf" : "-inf", fp);
    else std::fprintf(fp, "%.17e", v);
}

void write_component_header(std::FILE* fp, const char* prefix, int ncomp, bool& first)
{
    for (int i = 0; i < ncomp; ++i) {
        char name[64];
        std::snprintf(name, sizeof(name), "%s_%d", prefix, i);
        csv_write_name(fp, name, first);
    }
}

void write_csv_header(std::FILE* fp, const Param& param)
{
    bool first = true;
    csv_write_name(fp, "step", first);
    csv_write_name(fp, "time_s", first);
    for (int d = 0; d < NDIMS; ++d) {
        char name[32];
        std::snprintf(name, sizeof(name), "query_%s", axis_name(d));
        csv_write_name(fp, name, first);
    }
    csv_write_name(fp, "matched_node", first);
    csv_write_name(fp, "matched_elem", first);

    if (param.monitor.output_coord) {
        for (int d = 0; d < NDIMS; ++d) {
            char name[32];
            std::snprintf(name, sizeof(name), "coord_%s", axis_name(d));
            csv_write_name(fp, name, first);
        }
    }
    if (param.monitor.output_velocity) {
        for (int d = 0; d < NDIMS; ++d) {
            char name[32];
            std::snprintf(name, sizeof(name), "velocity_%s", axis_name(d));
            csv_write_name(fp, name, first);
        }
    }
    if (param.monitor.output_force) {
        for (int d = 0; d < NDIMS; ++d) {
            char name[32];
            std::snprintf(name, sizeof(name), "force_%s", axis_name(d));
            csv_write_name(fp, name, first);
        }
    }
    if (param.monitor.output_temperature) csv_write_name(fp, "temperature", first);
    if (param.monitor.output_pore_pressure) csv_write_name(fp, "pore_pressure", first);
    if (param.monitor.output_bcflag) csv_write_name(fp, "bcflag", first);

    if (param.monitor.output_stress) write_component_header(fp, "stress", NSTR, first);
    if (param.monitor.output_strain) write_component_header(fp, "strain", NSTR, first);
    if (param.monitor.output_strain_rate) write_component_header(fp, "strain_rate", NSTR, first);
    if (param.monitor.output_plastic_strain) csv_write_name(fp, "plastic_strain", first);
    if (param.monitor.output_plastic_strain_rate) csv_write_name(fp, "plastic_strain_rate", first);
    if (param.monitor.output_radiogenic_source) csv_write_name(fp, "radiogenic_source", first);
    if (param.monitor.output_density) csv_write_name(fp, "density", first);
    if (param.monitor.output_mesh_quality) csv_write_name(fp, "mesh_quality", first);
    if (param.monitor.output_viscosity) csv_write_name(fp, "viscosity", first);
    if (param.monitor.output_material) csv_write_name(fp, "material", first);
    if (param.monitor.output_dynamic_friction) csv_write_name(fp, "dynamic_friction", first);
    if (param.monitor.output_state_variable) csv_write_name(fp, "state_variable", first);

    std::fputc('\n', fp);
}

bool compute_elem_centroid(const Variables& var, int e, double c[NDIMS])
{
    if (e < 0 || e >= var.nelem) return false;
    ConstConnAccessor conn = (*var.connectivity)[e];
    for (int d = 0; d < NDIMS; ++d) c[d] = 0.0;
    for (int k = 0; k < NODES_PER_ELEM; ++k) {
        ConstArrayAccessor x = (*var.coord)[conn[k]];
        for (int d = 0; d < NDIMS; ++d) c[d] += x[d];
    }
    const double inv_n = 1.0 / static_cast<double>(NODES_PER_ELEM);
    for (int d = 0; d < NDIMS; ++d) c[d] *= inv_n;
    return true;
}

template <typename T>
double distance2_nd(T a, const double* b)
{
    double d2 = 0.0;
    for (int d = 0; d < NDIMS; ++d) {
        const double dv = a[d] - b[d];
        d2 += dv * dv;
    }
    return d2;
}

int find_nearest_node_id(const Variables& var, const double q[NDIMS], double& min_d2_out)
{
    min_d2_out = std::numeric_limits<double>::infinity();
    int best = -1;
    for (int n = 0; n < var.nnode; ++n) {
        const double d2 = distance2_nd((*var.coord)[n], q);
        if (d2 < min_d2_out) {
            min_d2_out = d2;
            best = n;
        }
    }
    return best;
}

int find_nearest_elem_id(const Variables& var, const double q[NDIMS], double& min_d2_out)
{
    min_d2_out = std::numeric_limits<double>::infinity();
    int best = -1;
    double c[NDIMS];
    for (int e = 0; e < var.nelem; ++e) {
        compute_elem_centroid(var, e, c);
        const double d2 = distance2_nd(c, q);
        if (d2 < min_d2_out) {
            min_d2_out = d2;
            best = e;
        }
    }
    return best;
}

void monitor_select_rebind_coord(const Param& param, const MonitorPointState& p, double out_q[NDIMS])
{
    if (param.monitor.remesh_rebind_mode == monitor_rebind_initial_coord) {
        for (int d = 0; d < NDIMS; ++d) out_q[d] = p.query_coord_initial[d];
        return;
    }
    for (int d = 0; d < NDIMS; ++d) out_q[d] = p.query_coord_rebind[d];
}

void monitor_match_all_points(const Param& param, const Variables& var, MonitorManager& manager)
{
    for (std::size_t i = 0; i < manager.points.size(); ++i) {
        MonitorPointState& p = manager.points[i];
        double q[NDIMS];
        monitor_select_rebind_coord(param, p, q);
        double node_d2 = 0.0;
        double elem_d2 = 0.0;
        p.node_id = find_nearest_node_id(var, q, node_d2);
        p.elem_id = find_nearest_elem_id(var, q, elem_d2);
        p.node_dist = std::sqrt(std::max(node_d2, 0.0));
        p.elem_dist = std::sqrt(std::max(elem_d2, 0.0));
    }
}

void monitor_open_files(const Param& param, MonitorManager& manager)
{
    for (std::size_t i = 0; i < manager.points.size(); ++i) {
        MonitorPointState& p = manager.points[i];
        char fname[512];
        std::snprintf(fname, sizeof(fname), "%s_point_%d.csv",
                      param.monitor.output_prefix.c_str(), p.id);
        p.fp = std::fopen(fname, "w");
        if (!p.fp) {
            std::cerr << "Error: cannot open monitor file '" << fname << "'\n";
            continue;
        }
        if (param.monitor.write_header) {
            write_csv_header(p.fp, param);
            std::fflush(p.fp);
        }
    }
}

int get_material_index(const Variables& var, int e)
{
    if (e < 0 || e >= var.nelem || var.elemmarkers == nullptr) return -1;
    const int_vec& a = (*var.elemmarkers)[e];
    if (a.empty()) return -1;
    return static_cast<int>(std::distance(a.begin(), std::max_element(a.begin(), a.end())));
}

double node_or_nan(const Variables& var, int node, int dim, const array_t* field)
{
    if (field == nullptr || node < 0 || node >= var.nnode) return std::numeric_limits<double>::quiet_NaN();
    return (*field)[node][dim];
}

double node_scalar_or_nan(const Variables& var, int node, const double_vec* field)
{
    if (field == nullptr || node < 0 || node >= var.nnode) return std::numeric_limits<double>::quiet_NaN();
    return (*field)[node];
}

double elem_scalar_or_nan(const Variables& var, int elem, const double_vec* field)
{
    if (field == nullptr || elem < 0 || elem >= var.nelem) return std::numeric_limits<double>::quiet_NaN();
    return (*field)[elem];
}

double elem_tensor_or_nan(const Variables& var, int elem, int comp, const tensor_t* field)
{
    if (field == nullptr || elem < 0 || elem >= var.nelem) return std::numeric_limits<double>::quiet_NaN();
    return (*field)[elem][comp];
}

void monitor_write_point_row(const Param& param,
                             const Variables& var,
                             const MonitorPointState& p,
                             std::FILE* fp)
{
    if (fp == nullptr) return;

    bool first = true;
    csv_write_int(fp, var.steps, first);
    csv_write_double(fp, var.time, first);

    double q[NDIMS];
    monitor_select_rebind_coord(param, p, q);
    for (int d = 0; d < NDIMS; ++d) csv_write_double(fp, q[d], first);

    csv_write_int(fp, p.node_id, first);
    csv_write_int(fp, p.elem_id, first);

    if (param.monitor.output_coord) {
        for (int d = 0; d < NDIMS; ++d) {
            csv_write_double(fp, node_or_nan(var, p.node_id, d, var.coord), first);
        }
    }
    if (param.monitor.output_velocity) {
        for (int d = 0; d < NDIMS; ++d) {
            csv_write_double(fp, node_or_nan(var, p.node_id, d, var.vel), first);
        }
    }
    if (param.monitor.output_force) {
        for (int d = 0; d < NDIMS; ++d) {
            csv_write_double(fp, node_or_nan(var, p.node_id, d, var.force), first);
        }
    }
    if (param.monitor.output_temperature) {
        csv_write_double(fp, node_scalar_or_nan(var, p.node_id, var.temperature), first);
    }
    if (param.monitor.output_pore_pressure) {
        csv_write_double(fp, node_scalar_or_nan(var, p.node_id, var.ppressure), first);
    }
    if (param.monitor.output_bcflag) {
        int bc = -1;
        if (var.bcflag != nullptr && p.node_id >= 0 && p.node_id < var.nnode) {
            bc = static_cast<int>((*var.bcflag)[p.node_id]);
        }
        csv_write_int(fp, bc, first);
    }

    if (param.monitor.output_stress) {
        for (int c = 0; c < NSTR; ++c) {
            csv_write_double(fp, elem_tensor_or_nan(var, p.elem_id, c, var.stress), first);
        }
    }
    if (param.monitor.output_strain) {
        for (int c = 0; c < NSTR; ++c) {
            csv_write_double(fp, elem_tensor_or_nan(var, p.elem_id, c, var.strain), first);
        }
    }
    if (param.monitor.output_strain_rate) {
        for (int c = 0; c < NSTR; ++c) {
            csv_write_double(fp, elem_tensor_or_nan(var, p.elem_id, c, var.strain_rate), first);
        }
    }
    if (param.monitor.output_plastic_strain) {
        csv_write_double(fp, elem_scalar_or_nan(var, p.elem_id, var.plstrain), first);
    }
    if (param.monitor.output_plastic_strain_rate) {
        csv_write_double(fp, elem_scalar_or_nan(var, p.elem_id, var.delta_plstrain), first);
    }
    if (param.monitor.output_radiogenic_source) {
        csv_write_double(fp, elem_scalar_or_nan(var, p.elem_id, var.radiogenic_source), first);
    }
    if (param.monitor.output_density) {
        double val = std::numeric_limits<double>::quiet_NaN();
        if (p.elem_id >= 0 && p.elem_id < var.nelem && var.mat != nullptr) {
            val = var.mat->rho(p.elem_id);
        }
        csv_write_double(fp, val, first);
    }
    if (param.monitor.output_mesh_quality) {
        double val = std::numeric_limits<double>::quiet_NaN();
        if (p.elem_id >= 0 && p.elem_id < var.nelem &&
            var.coord != nullptr && var.connectivity != nullptr && var.volume != nullptr) {
            val = elem_quality(*var.coord, *var.connectivity, *var.volume, p.elem_id);
        }
        csv_write_double(fp, val, first);
    }
    if (param.monitor.output_viscosity) {
        double val = std::numeric_limits<double>::quiet_NaN();
        if (p.elem_id >= 0 && p.elem_id < var.nelem && var.mat != nullptr) {
            val = var.mat->visc(p.elem_id);
        }
        csv_write_double(fp, val, first);
    }
    if (param.monitor.output_material) {
        csv_write_int(fp, get_material_index(var, p.elem_id), first);
    }
    if (param.monitor.output_dynamic_friction) {
        csv_write_double(fp, elem_scalar_or_nan(var, p.elem_id, var.dyn_fric_coeff), first);
    }
    if (param.monitor.output_state_variable) {
        csv_write_double(fp, elem_scalar_or_nan(var, p.elem_id, var.state_variable), first);
    }

    std::fputc('\n', fp);
    std::fflush(fp);
}

void monitor_write_all_points(const Param& param, const Variables& var, MonitorManager& manager)
{
    for (std::size_t i = 0; i < manager.points.size(); ++i) {
        monitor_write_point_row(param, var, manager.points[i], manager.points[i].fp);
    }
}

void monitor_capture_rebind_coords_before_remesh(const Variables& var, MonitorManager& manager)
{
    for (std::size_t i = 0; i < manager.points.size(); ++i) {
        MonitorPointState& p = manager.points[i];
        if (p.node_id >= 0 && p.node_id < var.nnode) {
            for (int d = 0; d < NDIMS; ++d) p.query_coord_rebind[d] = (*var.coord)[p.node_id][d];
            continue;
        }

        double c[NDIMS];
        if (compute_elem_centroid(var, p.elem_id, c)) {
            for (int d = 0; d < NDIMS; ++d) p.query_coord_rebind[d] = c[d];
        }
    }
}

} // namespace

void monitor_initialize(const Param& param, Variables& var)
{
    monitor_finalize(var);
    if (!param.monitor.enabled) return;

    g_monitor.points.resize(param.monitor.num_points);
    for (int i = 0; i < param.monitor.num_points; ++i) {
        MonitorPointState& p = g_monitor.points[i];
        p.id = i;
        p.query_coord_initial[0] = param.monitor.points_x[i];
        p.query_coord_initial[1] = param.monitor.points_z[i];
#ifdef THREED
        p.query_coord_initial[2] = param.monitor.points_z[i];
#endif
        for (int d = 0; d < NDIMS; ++d) {
            p.query_coord_rebind[d] = p.query_coord_initial[d];
        }
    }

    monitor_match_all_points(param, var, g_monitor);
    monitor_open_files(param, g_monitor);
    g_monitor.initialized = true;

    // Always dump the initial state at t=0 regardless of step_interval.
    if (var.steps == 0 && std::fabs(var.time) <= 1e-30) {
        monitor_write_all_points(param, var, g_monitor);
        g_monitor.last_written_step = var.steps;
    } else {
        g_monitor.last_written_step = std::numeric_limits<int>::min();
    }
}

void monitor_write_if_due(const Param& param, Variables& var)
{
    if (!param.monitor.enabled || !g_monitor.initialized) return;
    if (param.monitor.step_interval <= 0) return;
    if (var.steps % param.monitor.step_interval != 0) return;
    if (g_monitor.last_written_step == var.steps) return;

    monitor_write_all_points(param, var, g_monitor);
    g_monitor.last_written_step = var.steps;
}

void monitor_before_remesh(const Param& param, Variables& var)
{
    if (!param.monitor.enabled || !g_monitor.initialized) return;
    if (param.monitor.remesh_rebind_mode != monitor_rebind_pre_remesh_coord) return;
    monitor_capture_rebind_coords_before_remesh(var, g_monitor);
}

void monitor_remesh_update(const Param& param, Variables& var)
{
    if (!param.monitor.enabled || !g_monitor.initialized) return;
    if (param.monitor.remesh_rebind_mode == monitor_rebind_initial_coord) {
        for (std::size_t i = 0; i < g_monitor.points.size(); ++i) {
            for (int d = 0; d < NDIMS; ++d) {
                g_monitor.points[i].query_coord_rebind[d] = g_monitor.points[i].query_coord_initial[d];
            }
        }
    }
    monitor_match_all_points(param, var, g_monitor);
}

void monitor_finalize(Variables& /*var*/)
{
    for (std::size_t i = 0; i < g_monitor.points.size(); ++i) {
        if (g_monitor.points[i].fp != nullptr) {
            std::fclose(g_monitor.points[i].fp);
            g_monitor.points[i].fp = nullptr;
        }
    }
    g_monitor.points.clear();
    g_monitor.last_written_step = std::numeric_limits<int>::min();
    g_monitor.initialized = false;
}
