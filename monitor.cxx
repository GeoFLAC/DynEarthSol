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
#include <utility>
#include <vector>

#ifdef _WIN32
#include <io.h>
#else
#include <sys/types.h>
#include <unistd.h>
#endif

#include "monitor.hpp"
#include "geometry.hpp"
#include "matprops.hpp"
#include "utils.hpp"

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

#ifdef _WIN32
typedef __int64 monitor_file_offset_t;

monitor_file_offset_t monitor_file_tell(std::FILE* fp)
{
    return _ftelli64(fp);
}

int monitor_file_seek(std::FILE* fp, monitor_file_offset_t offset, int origin)
{
    return _fseeki64(fp, offset, origin);
}

int monitor_file_truncate(std::FILE* fp, monitor_file_offset_t offset)
{
    return _chsize_s(_fileno(fp), offset) == 0 ? 0 : -1;
}
#else
typedef off_t monitor_file_offset_t;

monitor_file_offset_t monitor_file_tell(std::FILE* fp)
{
    return ftello(fp);
}

int monitor_file_seek(std::FILE* fp, monitor_file_offset_t offset, int origin)
{
    return fseeko(fp, offset, origin);
}

int monitor_file_truncate(std::FILE* fp, monitor_file_offset_t offset)
{
    return ftruncate(fileno(fp), offset);
}
#endif

struct MonitorCsvRestartPlan {
    bool exists = false;
    bool has_checkpoint_row = false;
    monitor_file_offset_t keep_bytes = 0;
    std::vector<std::pair<int, double> > row_keys;
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

std::string make_csv_header(const Param& param)
{
    std::string header;
    bool first = true;
    csv_append_name(header, "step", first);
    csv_append_name(header, "time_s", first);
    for (int d = 0; d < NDIMS; ++d) {
        char name[32];
        std::snprintf(name, sizeof(name), "query_%s", axis_name(d));
        csv_append_name(header, name, first);
    }
    csv_append_name(header, "matched_node", first);
    csv_append_name(header, "matched_elem", first);

    if (param.monitor.output_coord) {
        for (int d = 0; d < NDIMS; ++d) {
            char name[32];
            std::snprintf(name, sizeof(name), "coord_%s", axis_name(d));
            csv_append_name(header, name, first);
        }
    }
    if (param.monitor.output_velocity) {
        for (int d = 0; d < NDIMS; ++d) {
            char name[32];
            std::snprintf(name, sizeof(name), "velocity_%s", axis_name(d));
            csv_append_name(header, name, first);
        }
    }
    if (param.monitor.output_force) {
        for (int d = 0; d < NDIMS; ++d) {
            char name[32];
            std::snprintf(name, sizeof(name), "force_%s", axis_name(d));
            csv_append_name(header, name, first);
        }
    }
    if (param.monitor.output_temperature) csv_append_name(header, "temperature", first);
    if (param.monitor.output_pore_pressure) csv_append_name(header, "pore_pressure", first);
    if (param.monitor.output_bcflag) csv_append_name(header, "bcflag", first);

    if (param.monitor.output_stress) append_component_header(header, "stress", NSTR, first);
    if (param.monitor.output_strain) append_component_header(header, "strain", NSTR, first);
    if (param.monitor.output_strain_rate) append_component_header(header, "strain_rate", NSTR, first);
    if (param.monitor.output_plastic_strain) csv_append_name(header, "plastic_strain", first);
    if (param.monitor.output_plastic_strain_rate) csv_append_name(header, "plastic_strain_rate", first);
    if (param.monitor.output_radiogenic_source) csv_append_name(header, "radiogenic_source", first);
    if (param.monitor.output_density) csv_append_name(header, "density", first);
    if (param.monitor.output_mesh_quality) csv_append_name(header, "mesh_quality", first);
    if (param.monitor.output_viscosity) csv_append_name(header, "viscosity", first);
    if (param.monitor.output_material) csv_append_name(header, "material", first);
    if (param.monitor.output_dynamic_friction) csv_append_name(header, "dynamic_friction", first);
    if (param.monitor.output_state_variable) csv_append_name(header, "state_variable", first);

    return header;
}

bool write_csv_header(std::FILE* fp, const Param& param)
{
    const std::string header = make_csv_header(param);
    return std::fwrite(header.data(), 1, header.size(), fp) == header.size() &&
           std::fputc('\n', fp) != EOF;
}

bool compute_elem_centroid(const Variables& var, int e, double c[NDIMS])
{
    if (e < 0 || e >= var.nelem) return false;
    ConstArrayIndirectAccessor coord = var.coord->view_const((*var.connectivity)[e]);
    for (int d = 0; d < NDIMS; ++d) c[d] = 0.0;
    for (int k = 0; k < NODES_PER_ELEM; ++k) {
        for (int d = 0; d < NDIMS; ++d) c[d] += coord[k][d];
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

std::string monitor_filename(const Param& param, int point_id)
{
    char suffix[64];
    std::snprintf(suffix, sizeof(suffix), "_point_%d.csv", point_id);
    return param.monitor.output_prefix + suffix;
}

bool read_monitor_csv_line(std::FILE* fp, const std::string& filename,
                           std::string& line,
                           monitor_file_offset_t& end_offset)
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

    end_offset = monitor_file_tell(fp);
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

void validate_monitor_csv_values(const std::string& filename,
                                 const std::string& line)
{
    const char* begin = line.c_str();
    while (true) {
        char* end = nullptr;
        std::strtod(begin, &end);
        if (end == begin || (*end != ',' && *end != '\0')) {
            monitor_csv_restart_error(filename, "has a non-numeric data column");
        }
        if (*end == '\0') return;
        begin = end + 1;
    }
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
    monitor_file_offset_t line_start = 0;
    monitor_file_offset_t line_end = 0;
    std::string line;

    if (param.monitor.write_header) {
        if (!read_monitor_csv_line(fp, filename, line, line_end)) {
            std::fclose(fp);
            monitor_csv_restart_error(filename, "is missing its CSV header");
        }
        if (line != expected_header) {
            std::fclose(fp);
            monitor_csv_restart_error(
                filename, "has a header that does not match the current monitor schema");
        }
        plan.keep_bytes = line_end;
        line_start = line_end;
    }

    bool stale_suffix = false;
    bool have_previous = false;
    int previous_step = -1;
    double previous_time = 0.0;
    std::string previous_line;

    while (read_monitor_csv_line(fp, filename, line, line_end)) {
        if (monitor_csv_column_count(line) != expected_columns) {
            std::fclose(fp);
            monitor_csv_restart_error(
                filename, "has a data row with the wrong number of columns");
        }
        validate_monitor_csv_values(filename, line);

        int step = -1;
        double time = 0.0;
        parse_monitor_csv_key(filename, line, step, time);
        if (step != 0 && step % param.monitor.step_interval != 0) {
            std::fclose(fp);
            monitor_csv_restart_error(
                filename, "does not match the current monitor step interval");
        }

        if (have_previous) {
            if (step < previous_step || time < previous_time) {
                std::fclose(fp);
                monitor_csv_restart_error(
                    filename, "has non-monotonic step/time rows");
            }
            if (step == previous_step) {
                if (time != previous_time || line != previous_line) {
                    std::fclose(fp);
                    monitor_csv_restart_error(
                        filename, "has ambiguous rows for the same step");
                }
                const bool discardable_duplicate = stale_suffix ||
                    (step == var.steps && time == var.time &&
                     plan.has_checkpoint_row);
                if (!discardable_duplicate) {
                    std::fclose(fp);
                    monitor_csv_restart_error(
                        filename, "has a duplicate row before the restart cutoff");
                }
                if (!stale_suffix) {
                    stale_suffix = true;
                    plan.keep_bytes = line_start;
                }
            }
        }

        if (!stale_suffix) {
            if (step < var.steps) {
                if (time > var.time) {
                    std::fclose(fp);
                    monitor_csv_restart_error(
                        filename, "has a pre-checkpoint step after the checkpoint time");
                }
                plan.keep_bytes = line_end;
            }
            else if (step == var.steps) {
                if (time != var.time) {
                    std::fclose(fp);
                    monitor_csv_restart_error(
                        filename, "has a checkpoint step with a different time");
                }
                plan.has_checkpoint_row = true;
                plan.keep_bytes = line_end;
            }
            else {
                stale_suffix = true;
                plan.keep_bytes = line_start;
            }
        }

        plan.row_keys.push_back(std::make_pair(step, time));
        have_previous = true;
        previous_step = step;
        previous_time = time;
        previous_line = line;
        line_start = line_end;
    }

    if (var.steps % param.monitor.step_interval == 0 &&
        !plan.has_checkpoint_row) {
        std::fclose(fp);
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
        if (!p.fp) {
            std::cerr << "Error: cannot open monitor file '" << filename << "'\n";
            continue;
        }
        if (param.monitor.write_header) {
            if (!write_csv_header(p.fp, param)) {
                std::cerr << "Error: cannot write monitor file '" << filename << "'\n";
                std::fclose(p.fp);
                p.fp = nullptr;
                continue;
            }
            std::fflush(p.fp);
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
        plans[i].row_keys.clear();
    }

    // Open every point before truncating any of them. A permission failure on
    // one point therefore leaves every existing history byte-for-byte intact.
    for (std::size_t i = 0; i < manager.points.size(); ++i) {
        MonitorPointState& p = manager.points[i];
        const std::string filename = monitor_filename(param, p.id);
        p.fp = std::fopen(filename.c_str(), plans[i].exists ? "r+b" : "w+b");
        if (p.fp == nullptr) {
            for (std::size_t j = 0; j < i; ++j) {
                std::fclose(manager.points[j].fp);
                manager.points[j].fp = nullptr;
            }
            monitor_csv_io_error(filename, "open for restart", EXIT_IO_OPEN);
        }
    }

    for (std::size_t i = 0; i < manager.points.size(); ++i) {
        MonitorPointState& p = manager.points[i];
        const std::string filename = monitor_filename(param, p.id);
        if (plans[i].exists) {
            if (monitor_file_seek(p.fp, 0, SEEK_END) != 0)
                monitor_csv_io_error(filename, "seek in", EXIT_IO_RW);
            const monitor_file_offset_t file_size = monitor_file_tell(p.fp);
            if (file_size < 0)
                monitor_csv_io_error(filename, "measure", EXIT_IO_RW);
            if (plans[i].keep_bytes < file_size &&
                monitor_file_truncate(p.fp, plans[i].keep_bytes) != 0) {
                monitor_csv_io_error(filename, "trim stale rows from", EXIT_IO_RW);
            }
        }
        else if (param.monitor.write_header && !write_csv_header(p.fp, param)) {
            monitor_csv_io_error(filename, "write", EXIT_IO_RW);
        }

        std::clearerr(p.fp);
        if (monitor_file_seek(p.fp, 0, SEEK_END) != 0)
            monitor_csv_io_error(filename, "seek in", EXIT_IO_RW);
        if (std::fflush(p.fp) != 0)
            monitor_csv_io_error(filename, "flush", EXIT_IO_RW);
    }

    if (!plans.empty() && plans[0].exists &&
        plans[0].has_checkpoint_row) {
        manager.last_written_step = var.steps;
    }
}

void monitor_open_files(const Param& param, const Variables& var,
                        MonitorManager& manager)
{
    if (param.sim.is_restarting)
        open_monitor_files_restart(param, var, manager);
    else
        open_monitor_files_fresh(param, manager);
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
#ifdef THREED
        p.query_coord_initial[1] = param.monitor.points_y[i];
        p.query_coord_initial[2] = param.monitor.points_z[i];
#else
        p.query_coord_initial[1] = param.monitor.points_z[i];
#endif
        for (int d = 0; d < NDIMS; ++d) {
            p.query_coord_rebind[d] = p.query_coord_initial[d];
        }
    }

    monitor_match_all_points(param, var, g_monitor);
    monitor_open_files(param, var, g_monitor);
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
