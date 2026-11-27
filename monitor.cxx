#include <algorithm>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

#include <sys/types.h>
#include <unistd.h>

#include "monitor.hpp"
#include "barycentric-fn.hpp"
#include "mpi_simulation.hpp"
#include "geometry.hpp"
#include "matprops.hpp"
#include "utils.hpp"

namespace {

// A point is located by its containing element and barycentric weights, as a marker is: it moves
// with the mesh (follows the material) and is relocated by its position after a remesh.
struct MonitorPointState {
    int id = -1;
    double position[NDIMS] = {0.0};          // material position before a remesh or from a restart
    int elem = -1;                           // containing element, -1 outside the mesh
    double weight[NODES_PER_ELEM] = {0.0};
    int node_id = -1;                        // node nearest the point
    int elem_id = -1;                        // element sampled: containing (particle) or nearest centroid (grid)
    bool warned_outside = false;
    std::FILE* fp = nullptr;
};

// The last CSV column, row: what a row records.
enum MonitorRowKind {
    monitor_row_regular = 0,        // the step_interval series
    monitor_row_before_remesh = 1,  // the state just before a remesh
    monitor_row_after_remesh = 2    // the remapped state just after it
};

struct MonitorManager {
    bool initialized = false;
    int last_written_step = std::numeric_limits<int>::min();
    std::vector<MonitorPointState> points;
};

MonitorManager g_monitor;

struct MonitorCsvRestartPlan {
    bool exists = false;
    bool has_checkpoint_row = false;
    off_t keep_bytes = 0;
    std::vector<std::tuple<int, int, double> > row_keys;   // step, row kind, time
};

void monitor_csv_restart_error(const std::string& filename,
                               const std::string& detail)
{
    std::cerr << "Error: monitor restart file '" << filename << "' "
              << detail << ".\n";
    die(EXIT_IO_RESTART);
}

void monitor_csv_io_error(const std::string& filename, const char* operation,
                          ExitCode code)
{
    std::cerr << "Error: cannot " << operation << " monitor file '"
              << filename << "'.\n";
    die(code);
}

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

inline void csv_append_name(std::string& header, const char* name, bool& first)
{
    if (!first) header.push_back(',');
    first = false;
    header += name;
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

void append_component_header(std::string& header, const char* prefix, int ncomp,
                             bool& first)
{
    for (int i = 0; i < ncomp; ++i) {
        char name[64];
        std::snprintf(name, sizeof(name), "%s_%d", prefix, i);
        csv_append_name(header, name, first);
    }
}

void append_axis_header(std::string& header, const char* prefix, bool& first)
{
    for (int d = 0; d < NDIMS; ++d) {
        char name[64];
        std::snprintf(name, sizeof(name), "%s_%s", prefix, axis_name(d));
        csv_append_name(header, name, first);
    }
}

int get_material_index(const Variables& var, int e)
{
    if (e < 0 || e >= var.nelem || var.elemmarkers == nullptr) return -1;
    const int_vec& a = (*var.elemmarkers)[e];
    if (a.empty()) return -1;
    return dominant_material(a);
}

// A nodal field at the point: interpolated (particle) or at the nearest node (grid).
template <typename ValueAtNode>
double sample_nodal(const Param& param, const Variables& var, const MonitorPointState& p,
                    ValueAtNode value)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (param.monitor.point_type == monitor_point_grid)
        return p.node_id < 0 ? nan : value(p.node_id);
    if (p.elem < 0) return nan;
    ConstConnAccessor conn = (*var.connectivity)[p.elem];
    double v = 0.0;
    for (int k = 0; k < NODES_PER_ELEM; ++k) v += p.weight[k] * value(conn[k]);
    return v;
}

double nodal(const Param& param, const Variables& var, const MonitorPointState& p, int dim, const array_t* field)
{
    if (field == nullptr) return std::numeric_limits<double>::quiet_NaN();
    return sample_nodal(param, var, p, [field, dim](int n) { return (*field)[n][dim]; });
}

double nodal_scalar(const Param& param, const Variables& var, const MonitorPointState& p, const double_vec* field)
{
    if (field == nullptr) return std::numeric_limits<double>::quiet_NaN();
    return sample_nodal(param, var, p, [field](int n) { return (*field)[n]; });
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

double elem_prop_or_nan(const Variables& var, int elem, double (MatProps::*prop)(int) const)
{
    if (var.mat == nullptr || elem < 0 || elem >= var.nelem) return std::numeric_limits<double>::quiet_NaN();
    return (var.mat->*prop)(elem);
}

double elem_quality_or_nan(const Variables& var, int elem)
{
    if (var.coord == nullptr || var.connectivity == nullptr || var.volume == nullptr ||
        elem < 0 || elem >= var.nelem) return std::numeric_limits<double>::quiet_NaN();
    return elem_quality(*var.coord, *var.connectivity, *var.volume, elem);
}

// A field a row can record, set by its monitor.output_* flag. A vector field has NDIMS columns
// named by axis, a tensor field NSTR columns named by component. An integer column has `integer`
// in place of `real`.
enum MonitorColumnShape { monitor_scalar, monitor_vector, monitor_tensor };

struct MonitorColumn {
    bool Monitor::*flag;
    const char* name;
    MonitorColumnShape shape;
    double (*real)(const Param&, const Variables&, const MonitorPointState&, int comp);
    int (*integer)(const Variables&, const MonitorPointState&);
};

// The optional columns, in the order the header and every row hold them.
const MonitorColumn kMonitorColumns[] = {
    {&Monitor::output_coord, "coord", monitor_vector,
     [](auto& param, auto& var, auto& p, int d) { return nodal(param, var, p, d, var.coord); }},
    {&Monitor::output_velocity, "velocity", monitor_vector,
     [](auto& param, auto& var, auto& p, int d) { return nodal(param, var, p, d, var.vel); }},
    {&Monitor::output_force, "force", monitor_vector,
     [](auto& param, auto& var, auto& p, int d) { return nodal(param, var, p, d, var.force); }},
    {&Monitor::output_temperature, "temperature", monitor_scalar,
     [](auto& param, auto& var, auto& p, int) { return nodal_scalar(param, var, p, var.temperature); }},
    {&Monitor::output_pore_pressure, "pore_pressure", monitor_scalar,
     [](auto& param, auto& var, auto& p, int) { return nodal_scalar(param, var, p, var.ppressure); }},
    {&Monitor::output_bcflag, "bcflag", monitor_scalar, nullptr,
     [](auto& var, auto& p) {
         return var.bcflag != nullptr && p.node_id >= 0 && p.node_id < var.nnode
             ? static_cast<int>((*var.bcflag)[p.node_id]) : -1;
     }},
    {&Monitor::output_stress, "stress", monitor_tensor,
     [](auto&, auto& var, auto& p, int c) { return elem_tensor_or_nan(var, p.elem_id, c, var.stress); }},
    {&Monitor::output_strain, "strain", monitor_tensor,
     [](auto&, auto& var, auto& p, int c) { return elem_tensor_or_nan(var, p.elem_id, c, var.strain); }},
    {&Monitor::output_strain_rate, "strain_rate", monitor_tensor,
     [](auto&, auto& var, auto& p, int c) { return elem_tensor_or_nan(var, p.elem_id, c, var.strain_rate); }},
    {&Monitor::output_plastic_strain, "plastic_strain", monitor_scalar,
     [](auto&, auto& var, auto& p, int) { return elem_scalar_or_nan(var, p.elem_id, var.plstrain); }},
    {&Monitor::output_plastic_strain_rate, "plastic_strain_rate", monitor_scalar,
     [](auto&, auto& var, auto& p, int) { return elem_scalar_or_nan(var, p.elem_id, var.delta_plstrain); }},
    {&Monitor::output_radiogenic_source, "radiogenic_source", monitor_scalar,
     [](auto&, auto& var, auto& p, int) { return elem_scalar_or_nan(var, p.elem_id, var.radiogenic_source); }},
    {&Monitor::output_density, "density", monitor_scalar,
     [](auto&, auto& var, auto& p, int) { return elem_prop_or_nan(var, p.elem_id, &MatProps::rho); }},
    {&Monitor::output_mesh_quality, "mesh_quality", monitor_scalar,
     [](auto&, auto& var, auto& p, int) { return elem_quality_or_nan(var, p.elem_id); }},
    {&Monitor::output_viscosity, "viscosity", monitor_scalar,
     [](auto&, auto& var, auto& p, int) { return elem_prop_or_nan(var, p.elem_id, &MatProps::visc); }},
    {&Monitor::output_material, "material", monitor_scalar, nullptr,
     [](auto& var, auto& p) { return get_material_index(var, p.elem_id); }},
    {&Monitor::output_dynamic_friction, "dynamic_friction", monitor_scalar,
     [](auto&, auto& var, auto& p, int) { return elem_scalar_or_nan(var, p.elem_id, var.dyn_fric_coeff); }},
    {&Monitor::output_state_variable, "state_variable", monitor_scalar,
     [](auto&, auto& var, auto& p, int) { return elem_scalar_or_nan(var, p.elem_id, var.state_variable); }},
};

int column_count(const MonitorColumn& column)
{
    return column.shape == monitor_vector ? NDIMS : column.shape == monitor_tensor ? NSTR : 1;
}

std::string make_csv_header(const Param& param)
{
    std::string header;
    bool first = true;
    csv_append_name(header, "step", first);
    csv_append_name(header, "time_s", first);
    append_axis_header(header, "query", first);
    csv_append_name(header, "matched_node", first);
    csv_append_name(header, "matched_elem", first);
    for (const MonitorColumn& column : kMonitorColumns) {
        if (!(param.monitor.*column.flag)) continue;
        if (column.shape == monitor_vector) append_axis_header(header, column.name, first);
        else if (column.shape == monitor_tensor) append_component_header(header, column.name, NSTR, first);
        else csv_append_name(header, column.name, first);
    }
    csv_append_name(header, "row", first);

    return header;
}

bool write_csv_header(std::FILE* fp, const Param& param)
{
    const std::string header = make_csv_header(param);
    return std::fwrite(header.data(), 1, header.size(), fp) == header.size() &&
           std::fputc('\n', fp) != EOF;
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

int find_nearest_node_id(const Variables& var, const double q[NDIMS])
{
    double min_d2 = std::numeric_limits<double>::infinity();
    int best = -1;
    for (int n = 0; n < var.nnode; ++n) {
        const double d2 = distance2_nd((*var.coord)[n], q);
        if (d2 < min_d2) {
            min_d2 = d2;
            best = n;
        }
    }
    return best;
}

// The element whose centroid is nearest q.
int find_nearest_elem_id(const Variables& var, const double q[NDIMS])
{
    double min_d2 = std::numeric_limits<double>::infinity();
    int best = -1;
    double c[NDIMS];
    for (int e = 0; e < var.nelem; ++e) {
        ConstArrayIndirectAccessor coord = var.coord->view_const((*var.connectivity)[e]);
        for (int d = 0; d < NDIMS; ++d) c[d] = 0.0;
        for (int k = 0; k < NODES_PER_ELEM; ++k) {
            for (int d = 0; d < NDIMS; ++d) c[d] += coord[k][d];
        }
        const double inv_n = 1.0 / static_cast<double>(NODES_PER_ELEM);
        for (int d = 0; d < NDIMS; ++d) c[d] *= inv_n;
        const double d2 = distance2_nd(c, q);
        if (d2 < min_d2) {
            min_d2 = d2;
            best = e;
        }
    }
    return best;
}

// Containing element and barycentric weights of q; elem = -1 if q is outside the mesh.
void locate(const Variables& var, const Barycentric_transformation& bary,
            const double q[NDIMS], MonitorPointState& p)
{
    p.elem = -1;
    for (int e = 0; e < var.nelem; ++e) {
        double eta[NDIMS];
        bary.transform(q, e, eta);
        if (!bary.is_inside(eta)) continue;
        double last = 1.0;
        for (int d = 0; d < NDIMS; ++d) {
            p.weight[d] = eta[d];
            last -= eta[d];
        }
        p.weight[NDIMS] = last;
        p.elem = e;
        return;
    }
}

// Current position of the material point, carried by the mesh; its last position once outside it.
void current_position(const Variables& var, const MonitorPointState& p, double x[NDIMS])
{
    if (p.elem < 0) {
        for (int d = 0; d < NDIMS; ++d) x[d] = p.position[d];
        return;
    }
    barycentric_point(var, p.elem, p.weight, x);
}

// Before a row is written: the node, and for grid the element, nearest the point.
void bind(const Param& param, const Variables& var, MonitorPointState& p)
{
    double x[NDIMS];
    current_position(var, p, x);
    p.node_id = find_nearest_node_id(var, x);
    p.elem_id = param.monitor.point_type == monitor_point_particle ? p.elem : find_nearest_elem_id(var, x);
    if (p.elem < 0 && !p.warned_outside) {
        std::cerr << "Warning: monitor point " << p.id << " has left the mesh; "
                  << (param.monitor.point_type == monitor_point_particle ? "its values are NaN.\n"
                                                                         : "it reports the nearest node.\n");
        p.warned_outside = true;
    }
}

std::string monitor_filename(const Param& param, int point_id)
{
    char suffix[64];
    std::snprintf(suffix, sizeof(suffix), "_point_%d.csv", point_id);
    return param.monitor.output_prefix + suffix;
}

bool read_monitor_csv_line(std::FILE* fp, const std::string& filename,
                           std::string& line,
                           off_t& end_offset)
{
    line.clear();
    int ch = 0;
    while ((ch = std::fgetc(fp)) != EOF) {
        line.push_back(static_cast<char>(ch));
        if (ch == '\n') break;
    }
    if (ch == EOF) {
        if (std::ferror(fp))
            monitor_csv_io_error(filename, "read", EXIT_IO_RW);
        if (line.empty()) return false;
        monitor_csv_restart_error(filename, "has an unterminated final row");
    }

    end_offset = ftello(fp);
    if (end_offset < 0)
        monitor_csv_io_error(filename, "locate rows in", EXIT_IO_RW);

    line.erase(line.size() - 1); // terminal '\n'
    if (!line.empty() && line[line.size() - 1] == '\r')
        line.erase(line.size() - 1);
    if (line.find('\r') != std::string::npos ||
        line.find('\n') != std::string::npos ||
        line.find('\0') != std::string::npos) {
        monitor_csv_restart_error(filename, "contains an invalid row delimiter");
    }
    return true;
}

std::size_t monitor_csv_column_count(const std::string& line)
{
    return 1 + static_cast<std::size_t>(
        std::count(line.begin(), line.end(), ','));
}

void parse_monitor_csv_key(const std::string& filename,
                           const std::string& line,
                           int& step, double& time)
{
    const char* begin = line.c_str();
    char* end = nullptr;
    errno = 0;
    const long parsed_step = std::strtol(begin, &end, 10);
    if (end == begin || *end != ',' || errno == ERANGE ||
        parsed_step < 0 || parsed_step > INT_MAX) {
        monitor_csv_restart_error(filename, "has a malformed step column");
    }

    begin = end + 1;
    const double parsed_time = std::strtod(begin, &end);
    if (end == begin || (*end != ',' && *end != '\0') ||
        !std::isfinite(parsed_time)) {
        monitor_csv_restart_error(filename, "has a malformed time_s column");
    }

    step = static_cast<int>(parsed_step);
    time = parsed_time;
}

int parse_monitor_csv_kind(const std::string& filename, const std::string& line)
{
    const std::size_t comma = line.rfind(',');
    const char* begin = line.c_str() + (comma == std::string::npos ? 0 : comma + 1);
    char* end = nullptr;
    const long kind = std::strtol(begin, &end, 10);
    if (end == begin || *end != '\0' || kind < monitor_row_regular || kind > monitor_row_after_remesh)
        monitor_csv_restart_error(filename, "has a malformed row column");
    return static_cast<int>(kind);
}

MonitorCsvRestartPlan inspect_monitor_restart_file(
    const Param& param, const Variables& var, const std::string& filename)
{
    MonitorCsvRestartPlan plan;
    errno = 0;
    std::FILE* fp = std::fopen(filename.c_str(), "rb");
    if (fp == nullptr) {
        if (errno == ENOENT) return plan;
        monitor_csv_io_error(filename, "open", EXIT_IO_OPEN);
    }
    plan.exists = true;

    const std::string expected_header = make_csv_header(param);
    const std::size_t expected_columns = monitor_csv_column_count(expected_header);
    off_t line_end = 0;
    std::string line;

    if (param.monitor.write_header) {
        if (!read_monitor_csv_line(fp, filename, line, line_end)) {
            monitor_csv_restart_error(filename, "is missing its CSV header");
        }
        if (line != expected_header) {
            monitor_csv_restart_error(
                filename, "has a header that does not match the current monitor schema");
        }
        plan.keep_bytes = line_end;
    }

    bool have_previous = false;
    int previous_step = -1;
    int previous_kind = -1;
    double previous_time = 0.0;

    // Keep the rows up to the checkpoint's regular row. The rest, including the
    // checkpoint step's remesh rows that the restart repeats, is rewritten unread.
    while (read_monitor_csv_line(fp, filename, line, line_end)) {
        int step = -1;
        double time = 0.0;
        parse_monitor_csv_key(filename, line, step, time);
        if (step > var.steps) break;
        if (monitor_csv_column_count(line) != expected_columns) {
            monitor_csv_restart_error(
                filename, "has a data row with the wrong number of columns");
        }
        const int kind = parse_monitor_csv_kind(filename,line);
        if (step == var.steps && kind != monitor_row_regular) break;

        if (kind == monitor_row_regular && step != 0 && step % param.monitor.step_interval != 0) {
            monitor_csv_restart_error(
                filename, "does not match the current monitor step interval");
        }
        if (have_previous) {
            if (step < previous_step || time < previous_time) {
                monitor_csv_restart_error(
                    filename, "has non-monotonic step/time rows");
            }
            // Rows of one step follow the order they are written in: regular, before, after.
            if (step == previous_step && (time != previous_time || kind <= previous_kind)) {
                monitor_csv_restart_error(
                    filename, "has duplicate or out-of-order rows for one step");
            }
        }
        if (step == var.steps) {
            if (time != var.time) {
                monitor_csv_restart_error(
                    filename, "has a checkpoint step with a different time");
            }
            plan.has_checkpoint_row = true;
        }
        else if (time > var.time) {
            monitor_csv_restart_error(
                filename, "has a pre-checkpoint step after the checkpoint time");
        }

        plan.keep_bytes = line_end;
        plan.row_keys.push_back(std::make_tuple(step, kind, time));
        have_previous = true;
        previous_step = step;
        previous_kind = kind;
        previous_time = time;
    }

    if (var.steps % param.monitor.step_interval == 0 &&
        !plan.has_checkpoint_row) {
        monitor_csv_restart_error(
            filename, "is missing the row corresponding to the restart checkpoint");
    }

    if (std::fclose(fp) != 0)
        monitor_csv_io_error(filename, "close", EXIT_IO_RW);
    return plan;
}

void open_monitor_files_fresh(const Param& param, MonitorManager& manager)
{
    for (std::size_t i = 0; i < manager.points.size(); ++i) {
        MonitorPointState& p = manager.points[i];
        const std::string filename = monitor_filename(param, p.id);
        p.fp = std::fopen(filename.c_str(), "w");
        if (!p.fp)
            monitor_csv_io_error(filename, "open", EXIT_IO_OPEN);
        if (param.monitor.write_header) {
            if (!write_csv_header(p.fp, param))
                monitor_csv_io_error(filename, "write", EXIT_IO_RW);
            if (std::fflush(p.fp) != 0)
                monitor_csv_io_error(filename, "flush", EXIT_IO_RW);
        }
    }
}

void open_monitor_files_restart(const Param& param, const Variables& var,
                                MonitorManager& manager)
{
    std::vector<MonitorCsvRestartPlan> plans;
    plans.reserve(manager.points.size());
    for (std::size_t i = 0; i < manager.points.size(); ++i) {
        const std::string filename = monitor_filename(param, manager.points[i].id);
        plans.push_back(inspect_monitor_restart_file(param, var, filename));
        if (i == 0) continue;

        const bool existence_matches = plans[i].exists == plans[0].exists;
        const bool rows_match = !plans[i].exists ||
            plans[i].row_keys == plans[0].row_keys;
        if (!existence_matches || !rows_match) {
            monitor_csv_restart_error(
                filename, "does not have the same history as the other monitor points");
        }
    }

    // Open every point before truncating any of them. A permission failure on
    // one point therefore leaves every existing history byte-for-byte intact.
    for (std::size_t i = 0; i < manager.points.size(); ++i) {
        MonitorPointState& p = manager.points[i];
        const std::string filename = monitor_filename(param, p.id);
        p.fp = std::fopen(filename.c_str(), plans[i].exists ? "r+b" : "w+b");
        if (p.fp == nullptr) {
            monitor_csv_io_error(filename, "open for restart", EXIT_IO_OPEN);
        }
    }

    for (std::size_t i = 0; i < manager.points.size(); ++i) {
        MonitorPointState& p = manager.points[i];
        const std::string filename = monitor_filename(param, p.id);
        if (plans[i].exists) {
            if (fseeko(p.fp, 0, SEEK_END) != 0)
                monitor_csv_io_error(filename, "seek in", EXIT_IO_RW);
            const off_t file_size = ftello(p.fp);
            if (file_size < 0)
                monitor_csv_io_error(filename, "measure", EXIT_IO_RW);
            if (plans[i].keep_bytes < file_size &&
                ftruncate(fileno(p.fp), plans[i].keep_bytes) != 0) {
                monitor_csv_io_error(filename, "trim stale rows from", EXIT_IO_RW);
            }
        }
        else if (param.monitor.write_header && !write_csv_header(p.fp, param)) {
            monitor_csv_io_error(filename, "write", EXIT_IO_RW);
        }

        std::clearerr(p.fp);
        if (fseeko(p.fp, 0, SEEK_END) != 0)
            monitor_csv_io_error(filename, "seek in", EXIT_IO_RW);
        if (std::fflush(p.fp) != 0)
            monitor_csv_io_error(filename, "flush", EXIT_IO_RW);
    }

    if (!plans.empty() && plans[0].exists &&
        plans[0].has_checkpoint_row) {
        manager.last_written_step = var.steps;
    }
}

void monitor_write_point_row(const Param& param,
                             const Variables& var,
                             const MonitorPointState& p, int kind)
{
    std::FILE* fp = p.fp;
    if (fp == nullptr) return;

    bool first = true;
    csv_write_int(fp, var.steps, first);
    csv_write_double(fp, var.time, first);

    double q[NDIMS];
    current_position(var, p, q);
    for (int d = 0; d < NDIMS; ++d) csv_write_double(fp, q[d], first);

    csv_write_int(fp, p.node_id, first);
    csv_write_int(fp, p.elem_id, first);

    for (const MonitorColumn& column : kMonitorColumns) {
        if (!(param.monitor.*column.flag)) continue;
        for (int c = 0; c < column_count(column); ++c) {
            if (column.integer) csv_write_int(fp, column.integer(var, p), first);
            else csv_write_double(fp, column.real(param, var, p, c), first);
        }
    }
    csv_write_int(fp,kind,first);

    std::fputc('\n', fp);
    std::fflush(fp);
}

void monitor_write_all_points(const Param& param, const Variables& var, MonitorManager& manager,
                              int kind = monitor_row_regular)
{
#ifdef USEMPI
    // The ranks hold the current mesh. Points are located on current
    // coordinates; then only the nodes and elements the rows read are gathered.
    if (var.mpi) var.mpi->geometry(var);
#endif
    for (std::size_t i = 0; i < manager.points.size(); ++i)
        bind(param, var, manager.points[i]);
#ifdef USEMPI
    if (var.mpi) {
        int_vec nodes;
        int_vec elements;
        for (const auto& p : manager.points) {
            if (p.node_id >= 0) nodes.push_back(p.node_id);
            for (int elem : {p.elem, p.elem_id}) {
                if (elem < 0) continue;
                elements.push_back(elem);
                for (int j = 0; j < NODES_PER_ELEM; ++j)
                    nodes.push_back((*var.connectivity)[elem][j]);
            }
        }
        std::sort(nodes.begin(), nodes.end());
        nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
        std::sort(elements.begin(), elements.end());
        elements.erase(std::unique(elements.begin(), elements.end()), elements.end());
        var.mpi->sample(var, nodes, elements);
    }
#endif
    for (std::size_t i = 0; i < manager.points.size(); ++i)
        monitor_write_point_row(param, var, manager.points[i], kind);
}

} // namespace

void monitor_initialize(const Param& param, Variables& var, const MonitorCheckpoint& saved)
{
    monitor_finalize();
    if (!param.monitor.enabled) return;

    const double_vec& positions = saved.positions;
    const int_vec& elements = saved.elements;
    const double_vec& weights = saved.weights;
    if (std::any_of(positions.begin(), positions.end(), [](double v) { return !std::isfinite(v); }))
        die(EXIT_IO_RESTART,"Invalid monitor checkpoint positions.");
    if (std::any_of(elements.begin(), elements.end(), [&var](int e) { return e < -1 || e >= var.nelem; }))
        die(EXIT_IO_RESTART,"Invalid monitor checkpoint element.");
    if (std::any_of(weights.begin(), weights.end(), [](double w) { return !std::isfinite(w); }))
        die(EXIT_IO_RESTART,"Invalid monitor checkpoint weight.");
    const bool resume = !positions.empty();
    g_monitor.points.resize(param.monitor.num_points);
    for (int i = 0; i < param.monitor.num_points; ++i) {
        MonitorPointState& p = g_monitor.points[i];
        p.id = i;
        // A restart continues the material point where the checkpoint left it.
        if (resume) {
            for (int d = 0; d < NDIMS; ++d)
                p.position[d] = positions[i*NDIMS + d];
            p.elem = elements[i];
            for (int k = 0; k < NODES_PER_ELEM; ++k)
                p.weight[k] = weights[i*NODES_PER_ELEM+k];
        }
        else {
            p.position[0] = param.monitor.points_x[i];
#ifdef THREED
            p.position[1] = param.monitor.points_y[i];
#endif
            p.position[NDIMS-1] = param.monitor.points_z[i];
        }
    }
    if (!resume) {
        Barycentric_transformation bary(*var.coord,*var.connectivity,*var.volume);
        #pragma acc wait
        for (MonitorPointState& p : g_monitor.points) locate(var,bary,p.position,p);
    }

    if (param.sim.is_restarting)
        open_monitor_files_restart(param, var, g_monitor);
    else
        open_monitor_files_fresh(param, g_monitor);
    g_monitor.initialized = true;

    // A fresh history always starts at t=0 regardless of step_interval.  A
    // frame-0 restart keeps the validated row that is already on disk.
    if (var.steps == 0 && std::fabs(var.time) <= 1e-30 &&
        g_monitor.last_written_step != var.steps) {
        monitor_write_all_points(param, var, g_monitor);
        g_monitor.last_written_step = var.steps;
    } else if (g_monitor.last_written_step != var.steps) {
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
#ifdef USEMPI
    // The rows and positions read the whole mesh, which only the root holds.
    if (var.mpi) var.mpi->collect(var);
#endif
    // A regular row of this step already records the state before the remesh.
    if (g_monitor.last_written_step != var.steps)
        monitor_write_all_points(param, var, g_monitor, monitor_row_before_remesh);
    for (MonitorPointState& p : g_monitor.points)
        current_position(var, p, p.position);
}

void monitor_remesh_update(const Param& param, Variables& var)
{
    if (!param.monitor.enabled || !g_monitor.initialized) return;
    Barycentric_transformation bary(*var.coord,*var.connectivity,*var.volume);
    #pragma acc wait
    for (MonitorPointState& p : g_monitor.points)
        locate(var,bary,p.position,p);
    if (param.monitor.remesh_output == monitor_remesh_both)
        monitor_write_all_points(param, var, g_monitor, monitor_row_after_remesh);
}

MonitorCheckpoint monitor_checkpoint(const Param& param, const Variables& var)
{
    MonitorCheckpoint saved;
    if (!param.monitor.enabled || !g_monitor.initialized) return saved;
    const std::size_t count = g_monitor.points.size();
    saved.positions.resize(count*NDIMS);
    saved.elements.resize(count);
    saved.weights.assign(count*NODES_PER_ELEM,0);
    for (std::size_t i = 0; i < count; ++i) {
        const MonitorPointState& p = g_monitor.points[i];
        current_position(var,p,&saved.positions[i*NDIMS]);
        saved.elements[i] = p.elem;
        if (p.elem >= 0)
            for (int k = 0; k < NODES_PER_ELEM; ++k) saved.weights[i*NODES_PER_ELEM+k] = p.weight[k];
    }
    return saved;
}

void monitor_finalize()
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
