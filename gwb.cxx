#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "constants.hpp"
#include "parameters.hpp"
#include "utils.hpp"
#include "gwb.hpp"

#ifdef USEGWB

#include "world_builder/world.h"
#include "world_builder/wrapper_c.h"
#include "rapidjson/document.h"

namespace {

void *world = nullptr;
std::vector<int> fault_tags;            // GWB tag indices of the world's fault features
std::vector<std::string> fault_names;   // their names, for the distance to each fault's plane
// A 2D model's cross section: DES x maps to the world's map point cs_origin + x * cs_dir.
double cs_origin[2] = {0, 0}, cs_dir[2] = {1, 0};

// GWB reports a bad world file or query by throwing; turn that into a config error.
[[noreturn]] void gwb_failed(const Param &param, const char *what)
{
    std::cerr << "Error: Geodynamic World Builder, world '" << param.ic.gwb_filename
              << "': " << what << '\n';
    die(EXIT_CONFIG_VALUE);
}

// GWB's tag property is the index of a feature's tag in the list of distinct tags, in feature
// order; a feature's tag is its "tag" entry or, if empty, its model name. Rebuild that list from
// the world file, as GWB does, and keep the indices that a "fault" feature's tag has, the fault
// names, and the cross section.
void read_world_file(const Param &param)
{
    std::ifstream f(param.ic.gwb_filename);
    std::stringstream ss;
    ss << f.rdbuf();
    rapidjson::Document doc;
    // GWB's own flags: comments and nan/inf allowed
    doc.Parse<rapidjson::kParseCommentsFlag | rapidjson::kParseNanAndInfFlag>(ss.str().c_str());
    std::vector<std::string> tags;
    fault_tags.clear();
    fault_names.clear();
    if (doc.HasParseError() || !doc.IsObject() || !doc.HasMember("features") || !doc["features"].IsArray())
        return;   // GWB read the file, so this cannot happen; without a features array there is no fault
    if (doc.HasMember("cross section")) {   // GWB checked it: two points of two numbers
        const auto &cs = doc["cross section"];
        const double dx = cs[1][0].GetDouble() - cs[0][0].GetDouble();
        const double dy = cs[1][1].GetDouble() - cs[0][1].GetDouble();
        const double len = std::sqrt(dx*dx + dy*dy);
        cs_origin[0] = cs[0][0].GetDouble();
        cs_origin[1] = cs[0][1].GetDouble();
        cs_dir[0] = dx / len;
        cs_dir[1] = dy / len;
    }
    for (const auto &feat : doc["features"].GetArray()) {
        if (!feat.IsObject() || !feat.HasMember("model") || !feat["model"].IsString()) continue;
        const std::string model = feat["model"].GetString();
        std::string tag = (feat.HasMember("tag") && feat["tag"].IsString()) ? feat["tag"].GetString() : "";
        if (tag.empty()) tag = model;
        auto it = std::find(tags.begin(), tags.end(), tag);
        const int index = static_cast<int>(it - tags.begin());
        if (it == tags.end()) tags.push_back(tag);
        if (model != "fault") continue;
        if (std::find(fault_tags.begin(), fault_tags.end(), index) == fault_tags.end())
            fault_tags.push_back(index);
        fault_names.push_back((feat.HasMember("name") && feat["name"].IsString()) ? feat["name"].GetString() : "");
    }
}

// The GWB properties (temperature: {1,0,0}, composition c: {2,c,0}) at DES point x.
void query(const Param &param, const double *x, const unsigned int props[][3],
           unsigned int nprops, double *values)
{
    if (!world)
        die(EXIT_INTERNAL_ASSERT, "Geodynamic World Builder queried before gwb_build().");
    const double depth = -x[NDIMS-1];
    // GWB's vertical coordinate is the height above the base; only the depth places features.
    const double z = param.mesh.zlength + x[NDIMS-1];
    try {
#ifdef THREED
        properties_3d(world, x[0], x[1], z, depth, props, nprops, values);
#else
        properties_2d(world, x[0], z, depth, props, nprops, values);
#endif
    }
    catch (const std::exception &e) { gwb_failed(param, e.what()); }
    catch (...) { gwb_failed(param, "unknown error in a query"); }
}

} // anonymous namespace


void gwb_build(const Param &param)
{
    if (param.ic.temperature_option != 936 && param.ic.mattype_option != 936 &&
        param.ic.weakzone_option != 936)
        return;
    std::cout << "  Reading the Geodynamic World Builder world " << param.ic.gwb_filename << '\n';
    const bool has_output_dir = false;
    try {
        create_world(&world, param.ic.gwb_filename.c_str(), &has_output_dir, "",
                     param.markers.random_seed);
    }
    catch (const std::exception &e) { gwb_failed(param, e.what()); }
    catch (...) { gwb_failed(param, "unknown error while reading the world"); }
    read_world_file(param);
    if (param.ic.weakzone_option != 936) return;
    if (fault_tags.empty()) {
        std::cerr << "Error: ic.weakzone_option 936 needs a fault feature in the Geodynamic World "
                     "Builder world '" << param.ic.gwb_filename << "'; it has none.\n";
        die(EXIT_CONFIG_DATA);
    }
    // GWB finds a fault's plane by its name, the first feature of that name.
    if (param.ic.gwb_fault_gaussian)
        for (std::size_t i = 0; i < fault_names.size(); ++i)
            if (std::count(fault_names.begin(), fault_names.end(), fault_names[i]) > 1) {
                std::cerr << "Error: ic.gwb_fault_gaussian needs the faults of the Geodynamic "
                             "World Builder world '" << param.ic.gwb_filename << "' to have distinct "
                             "names; '" << fault_names[i] << "' is used more than once.\n";
                die(EXIT_CONFIG_DATA);
            }
}


double gwb_temperature(const Param &param, const double *x)
{
    const unsigned int props[1][3] = {{1, 0, 0}};
    double t;
    query(param, x, props, 1, &t);
    return t;
}


int gwb_mattype(const Param &param, const double *x)
{
    // GWB composition i is mattype i. The composition holding at least half of the point sets
    // its mattype; the world must give every marker one.
    const int n = param.mat.nmat;
    std::vector<unsigned int> props(3 * n);
    std::vector<double> c(n);
    for (int i = 0; i < n; ++i) {
        props[3*i] = 2; props[3*i+1] = i; props[3*i+2] = 0;
    }
    query(param, x, reinterpret_cast<const unsigned int (*)[3]>(props.data()), n, c.data());
    int best = 0;
    for (int i = 1; i < n; ++i)
        if (c[i] > c[best]) best = i;
    if (c[best] < 0.5) {
        std::cerr << "Error: the Geodynamic World Builder world '" << param.ic.gwb_filename
                  << "' defines no material at x = (";
        for (int d = 0; d < NDIMS; ++d) std::cerr << (d ? ", " : "") << x[d];
        std::cerr << ") m: no composition 0 to " << n-1 << " (mat.num_materials - 1) reaches 0.5 there"
                  << " (fractions";
        for (int i = 0; i < n; ++i) std::cerr << ' ' << c[i];
        std::cerr << "). With ic.mattype_option 936 the world must give every point of the domain a "
                     "composition below mat.num_materials, the mantle included.\n";
        die(EXIT_CONFIG_DATA);
    }
    return best;
}


double gwb_fault_weight(const Param &param, const double *x)
{
    const unsigned int props[1][3] = {{4, 0, 0}};   // tag; -1 outside every feature
    double tag;
    query(param, x, props, 1, &tag);
    if (std::find(fault_tags.begin(), fault_tags.end(), static_cast<int>(tag)) == fault_tags.end())
        return 0;
    if (!param.ic.gwb_fault_gaussian) return 1;
    const double sigma = param.ic.weakzone_standard_deviation;

    // The point in the world's 3D Cartesian frame, as properties_2d/3d place it.
    const double depth = -x[NDIMS-1];
#ifdef THREED
    const std::array<double,3> p = {{x[0], x[1], param.mesh.zlength + x[NDIMS-1]}};
#else
    const std::array<double,3> p = {{cs_origin[0] + x[0] * cs_dir[0], cs_origin[1] + x[0] * cs_dir[1],
                                     param.mesh.zlength + x[NDIMS-1]}};
#endif
    // The distance from each fault's plane; inf for a fault the point is not beside. The tag
    // above, not this distance, decides whether the point is in a fault: the distance ignores
    // the fault's max depth.
    const auto *w = static_cast<const WorldBuilder::World *>(world);
    double weight = 0;
    try {
        for (const auto &name : fault_names) {
            const double d = w->distance_to_plane(p, depth, name).get_distance_from_surface();
            if (std::isfinite(d))
                weight = std::max(weight, std::exp(-d * d / (2 * sigma * sigma)));
        }
    }
    catch (const std::exception &e) { gwb_failed(param, e.what()); }
    catch (...) { gwb_failed(param, "unknown error in a fault distance"); }
    return weight;
}


void gwb_release()
{
    if (world) release_world(world);
    world = nullptr;
    fault_tags.clear();
}

#else  // !USEGWB: validate_parameters rejects option 936 in this build

void gwb_build(const Param &) {}

double gwb_temperature(const Param &, const double *)
{
    die(EXIT_UNSUPPORTED_LIB, "Geodynamic World Builder is not compiled in (build with usegwb=1).");
}

int gwb_mattype(const Param &, const double *)
{
    die(EXIT_UNSUPPORTED_LIB, "Geodynamic World Builder is not compiled in (build with usegwb=1).");
}

double gwb_fault_weight(const Param &, const double *)
{
    die(EXIT_UNSUPPORTED_LIB, "Geodynamic World Builder is not compiled in (build with usegwb=1).");
}

void gwb_release() {}

#endif
