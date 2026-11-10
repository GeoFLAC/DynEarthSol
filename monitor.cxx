#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdint>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <io.h>
#else
#include <sys/types.h>
#include <unistd.h>
#endif

#include "monitor.hpp"
#include "barycentric-fn.hpp"
#include "geometry.hpp"
#include "matprops.hpp"
#include "utils.hpp"

namespace {

// A point is located by its containing element and barycentric weights, as a marker is: it moves
// with the mesh (follows the material) and is relocated by its position after a remesh.
struct MonitorPointState {
    int id = -1;
    double query[NDIMS] = {0.0};
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
    bool legacy_csv = false;
    int last_written_step = std::numeric_limits<int>::min();
    std::vector<MonitorPointState> points;
    double_vec restart_positions;
    int_vec restart_elements;
    double_vec restart_weights;
};

MonitorManager g_monitor;

#ifdef _WIN32
typedef std::int64_t monitor_file_offset_t;

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
    bool legacy_csv = false;
    bool has_checkpoint_row = false;
    monitor_file_offset_t keep_bytes = 0;
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
    csv_append_name(header, "row", first);

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

int find_nearest_elem_id(const Variables& var, const double q[NDIMS])
{
    double min_d2 = std::numeric_limits<double>::infinity();
    int best = -1;
    double c[NDIMS];
    for (int e = 0; e < var.nelem; ++e) {
        compute_elem_centroid(var, e, c);
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
    ConstArrayIndirectAccessor coord = var.coord->view_const((*var.connectivity)[p.elem]);
    for (int d = 0; d < NDIMS; ++d) {
        x[d] = 0.0;
        for (int k = 0; k < NODES_PER_ELEM; ++k) x[d] += p.weight[k] * coord[k][d];
    }
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
    std::size_t expected_columns = monitor_csv_column_count(expected_header);
    monitor_file_offset_t line_start = 0;
    monitor_file_offset_t line_end = 0;
    std::string line;

    if (param.monitor.write_header) {
        if (!read_monitor_csv_line(fp, filename, line, line_end)) {
            std::fclose(fp);
            monitor_csv_restart_error(filename, "is missing its CSV header");
        }
        if (line == expected_header.substr(0,expected_header.size()-4)) {
            plan.legacy_csv = true;
            --expected_columns;
        } else if (line != expected_header) {
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
    int previous_kind = -1;
    double previous_time = 0.0;
    std::string previous_line;

    while (read_monitor_csv_line(fp, filename, line, line_end)) {
        if (!param.monitor.write_header && !have_previous &&
            monitor_csv_column_count(line) == expected_columns-1) {
            plan.legacy_csv = true;
            --expected_columns;
        }
        if (monitor_csv_column_count(line) != expected_columns) {
            std::fclose(fp);
            monitor_csv_restart_error(
                filename, "has a data row with the wrong number of columns");
        }
        validate_monitor_csv_values(filename, line);

        int step = -1;
        double time = 0.0;
        parse_monitor_csv_key(filename, line, step, time);
        const int kind = plan.legacy_csv ? monitor_row_regular : parse_monitor_csv_kind(filename,line);
        if (kind == monitor_row_regular && step != 0 && step % param.monitor.step_interval != 0) {
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
            if (step == previous_step && time == previous_time && kind != previous_kind) {
                // Rows of one step follow the order they are written in: regular, before, after.
                if (kind < previous_kind) {
                    std::fclose(fp);
                    monitor_csv_restart_error(
                        filename, "has rows of one step out of order");
                }
            }
            else if (step == previous_step) {
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
                // The checkpoint precedes the step's remesh, which the restart repeats:
                // its remesh rows are written again.
                if (kind == monitor_row_regular) {
                    plan.has_checkpoint_row = true;
                    plan.keep_bytes = line_end;
                }
                else {
                    stale_suffix = true;
                    plan.keep_bytes = line_start;
                }
            }
            else {
                stale_suffix = true;
                plan.keep_bytes = line_start;
            }
        }

        plan.row_keys.push_back(std::make_tuple(step, kind, time));
        have_previous = true;
        previous_step = step;
        previous_kind = kind;
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
        if (!existence_matches || !rows_match || plans[i].legacy_csv != plans[0].legacy_csv) {
            monitor_csv_restart_error(
                filename, "does not have the same history as the other monitor points (including CSV schema)");
        }
        std::vector<std::tuple<int, int, double> >().swap(plans[i].row_keys);
    }

    manager.legacy_csv = !plans.empty() && plans[0].legacy_csv;
    if (manager.legacy_csv)
        std::printf("Continuing legacy monitor CSVs with regular rows. Use a new monitor.output_prefix for remesh rows.\n");

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

// A nodal field at the point: interpolated (particle) or at the nearest node (grid).
double nodal(const Param& param, const Variables& var, const MonitorPointState& p, int dim, const array_t* field)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (field == nullptr) return nan;
    if (param.monitor.point_type == monitor_point_grid)
        return p.node_id < 0 ? nan : (*field)[p.node_id][dim];
    if (p.elem < 0) return nan;
    ConstConnAccessor conn = (*var.connectivity)[p.elem];
    double v = 0.0;
    for (int k = 0; k < NODES_PER_ELEM; ++k) v += p.weight[k] * (*field)[conn[k]][dim];
    return v;
}

double nodal_scalar(const Param& param, const Variables& var, const MonitorPointState& p, const double_vec* field)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (field == nullptr) return nan;
    if (param.monitor.point_type == monitor_point_grid)
        return p.node_id < 0 ? nan : (*field)[p.node_id];
    if (p.elem < 0) return nan;
    ConstConnAccessor conn = (*var.connectivity)[p.elem];
    double v = 0.0;
    for (int k = 0; k < NODES_PER_ELEM; ++k) v += p.weight[k] * (*field)[conn[k]];
    return v;
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
                             std::FILE* fp, int kind, bool legacy)
{
    if (fp == nullptr) return;

    bool first = true;
    csv_write_int(fp, var.steps, first);
    csv_write_double(fp, var.time, first);

    double q[NDIMS];
    current_position(var, p, q);
    for (int d = 0; d < NDIMS; ++d) csv_write_double(fp, q[d], first);

    csv_write_int(fp, p.node_id, first);
    csv_write_int(fp, p.elem_id, first);

    if (param.monitor.output_coord) {
        for (int d = 0; d < NDIMS; ++d) {
            csv_write_double(fp, nodal(param, var, p, d, var.coord), first);
        }
    }
    if (param.monitor.output_velocity) {
        for (int d = 0; d < NDIMS; ++d) {
            csv_write_double(fp, nodal(param, var, p, d, var.vel), first);
        }
    }
    if (param.monitor.output_force) {
        for (int d = 0; d < NDIMS; ++d) {
            csv_write_double(fp, nodal(param, var, p, d, var.force), first);
        }
    }
    if (param.monitor.output_temperature) {
        csv_write_double(fp, nodal_scalar(param, var, p, var.temperature), first);
    }
    if (param.monitor.output_pore_pressure) {
        csv_write_double(fp, nodal_scalar(param, var, p, var.ppressure), first);
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
    if (!legacy) csv_write_int(fp,kind,first);

    std::fputc('\n', fp);
    std::fflush(fp);
}

void monitor_write_all_points(const Param& param, const Variables& var, MonitorManager& manager,
                              int kind = monitor_row_regular)
{
    // Legacy files have no row-kind column: preserve their schema and regular
    // series without appending ambiguous off-interval remesh samples.
    if (manager.legacy_csv && kind != monitor_row_regular) return;
    for (std::size_t i = 0; i < manager.points.size(); ++i) {
        bind(param, var, manager.points[i]);
        monitor_write_point_row(param,var,manager.points[i],manager.points[i].fp,kind,manager.legacy_csv);
    }
}

} // namespace

void monitor_initialize(const Param& param, Variables& var)
{
    monitor_finalize(var);
    if (!param.monitor.enabled) return;

    const bool resume = g_monitor.restart_positions.size() == std::size_t(param.monitor.num_points) * NDIMS;
    const bool reference = resume && g_monitor.restart_elements.size() == std::size_t(param.monitor.num_points);
    g_monitor.points.resize(param.monitor.num_points);
    for (int i = 0; i < param.monitor.num_points; ++i) {
        MonitorPointState& p = g_monitor.points[i];
        p.id = i;
        p.query[0] = param.monitor.points_x[i];
#ifdef THREED
        p.query[1] = param.monitor.points_y[i];
#endif
        p.query[NDIMS-1] = param.monitor.points_z[i];
        // A restart continues the material point where the checkpoint left it.
        for (int d = 0; d < NDIMS; ++d)
            p.position[d] = resume ? g_monitor.restart_positions[i*NDIMS + d] : p.query[d];
        if (reference) {
            p.elem = g_monitor.restart_elements[i];
            for (int k = 0; k < NODES_PER_ELEM; ++k)
                p.weight[k] = g_monitor.restart_weights[i*NODES_PER_ELEM+k];
        }
    }
    if (!reference) {
        Barycentric_transformation bary(*var.coord,*var.connectivity,*var.volume);
        #pragma acc wait
        for (MonitorPointState& p : g_monitor.points) locate(var,bary,p.position,p);
    }
    g_monitor.restart_positions.clear();
    g_monitor.restart_elements.clear();
    g_monitor.restart_weights.clear();

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

void monitor_checkpoint(const Param& param, const Variables& var,
                        double_vec& positions, int_vec& elements, double_vec& weights)
{
    positions.clear(); elements.clear(); weights.clear();
    if (!param.monitor.enabled || !g_monitor.initialized) return;
    const std::size_t count = g_monitor.points.size();
    positions.resize(count*NDIMS);
    elements.resize(count);
    weights.assign(count*NODES_PER_ELEM,0);
    for (std::size_t i = 0; i < count; ++i) {
        const MonitorPointState& p = g_monitor.points[i];
        current_position(var,p,&positions[i*NDIMS]);
        elements[i] = p.elem;
        if (p.elem >= 0)
            for (int k = 0; k < NODES_PER_ELEM; ++k) weights[i*NODES_PER_ELEM+k] = p.weight[k];
    }
}

void monitor_restore(const Variables& var, const double_vec& positions,
                     const int_vec& elements, const double_vec& weights)
{
    const std::size_t count = positions.size()/NDIMS;
    if (positions.size()%NDIMS ||
        std::any_of(positions.begin(),positions.end(),[](double v) { return !std::isfinite(v); }))
        die(EXIT_IO_RESTART,"Invalid monitor checkpoint positions.");
    const bool reference = !elements.empty() || !weights.empty();
    if (reference && (elements.size() != count || weights.size() != count*NODES_PER_ELEM))
        die(EXIT_IO_RESTART,"Inconsistent monitor checkpoint reference sizes.");
    if (reference) {
        for (std::size_t i = 0; i < count; ++i) {
            const int e = elements[i];
            if (e < -1 || e >= var.nelem)
                die(EXIT_IO_RESTART,"Invalid monitor checkpoint element.");
            double sum = 0;
            for (int k = 0; k < NODES_PER_ELEM; ++k) {
                const double w = weights[i*NODES_PER_ELEM+k];
                if (!std::isfinite(w) || (e >= 0 && (w < -1e-10 || w > 1+NODES_PER_ELEM*1e-10)))
                    die(EXIT_IO_RESTART,"Invalid monitor checkpoint weight.");
                sum += w;
            }
            if (e < 0) continue;
            if (std::fabs(sum-1) > 1e-10)
                die(EXIT_IO_RESTART,"Monitor checkpoint weights do not sum to one.");
            for (int d = 0; d < NDIMS; ++d) {
                double point = 0, scale = 1;
                for (int k = 0; k < NODES_PER_ELEM; ++k) {
                    const double x = (*var.coord)[(*var.connectivity)[e][k]][d];
                    point += weights[i*NODES_PER_ELEM+k]*x;
                    scale = std::max(scale,std::fabs(x));
                }
                if (std::fabs(point-positions[i*NDIMS+d]) > 128*std::numeric_limits<double>::epsilon()*scale)
                    die(EXIT_IO_RESTART,"Monitor checkpoint position and weights disagree.");
            }
        }
    }
    g_monitor.restart_positions = positions;
    g_monitor.restart_elements = elements;
    g_monitor.restart_weights = weights;
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
    g_monitor.legacy_csv = false;
}
