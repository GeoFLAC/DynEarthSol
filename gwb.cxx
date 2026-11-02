#include <algorithm>
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

#include "world_builder/wrapper_c.h"
#include "rapidjson/document.h"

namespace {

void *world = nullptr;
std::vector<int> fault_tags;   // GWB tag indices of the world's fault features

// GWB reports a bad world file or query by throwing; turn that into a config error.
[[noreturn]] void gwb_failed(const Param &param, const char *what)
{
    std::cerr << "Error: Geodynamic World Builder, world '" << param.ic.gwb_filename
              << "': " << what << '\n';
    die(EXIT_CONFIG_VALUE);
}

// GWB's tag property is the index of a feature's tag in the list of distinct tags, in feature
// order; a feature's tag is its "tag" entry or, if empty, its model name. Rebuild that list from
// the world file, as GWB does, and keep the indices that a "fault" feature's tag has.
void find_fault_tags(const Param &param)
{
    std::ifstream f(param.ic.gwb_filename);
    std::stringstream ss;
    ss << f.rdbuf();
    rapidjson::Document doc;
    // GWB's own flags: comments and nan/inf allowed
    doc.Parse<rapidjson::kParseCommentsFlag | rapidjson::kParseNanAndInfFlag>(ss.str().c_str());
    std::vector<std::string> tags;
    fault_tags.clear();
    if (doc.HasParseError() || !doc.IsObject() || !doc.HasMember("features") || !doc["features"].IsArray())
        return;   // GWB read the file, so this cannot happen; without a features array there is no fault
    for (const auto &feat : doc["features"].GetArray()) {
        if (!feat.IsObject() || !feat.HasMember("model") || !feat["model"].IsString()) continue;
        const std::string model = feat["model"].GetString();
        std::string tag = (feat.HasMember("tag") && feat["tag"].IsString()) ? feat["tag"].GetString() : "";
        if (tag.empty()) tag = model;
        auto it = std::find(tags.begin(), tags.end(), tag);
        const int index = static_cast<int>(it - tags.begin());
        if (it == tags.end()) tags.push_back(tag);
        if (model == "fault" && std::find(fault_tags.begin(), fault_tags.end(), index) == fault_tags.end())
            fault_tags.push_back(index);
    }
}

void *open_world(const Param &param)
{
    if (world) return world;
    std::cout << "  Reading the Geodynamic World Builder world " << param.ic.gwb_filename << '\n';
    const bool has_output_dir = false;
    try {
        create_world(&world, param.ic.gwb_filename.c_str(), &has_output_dir, "",
                     param.markers.random_seed);
    }
    catch (const std::exception &e) { gwb_failed(param, e.what()); }
    catch (...) { gwb_failed(param, "unknown error while reading the world"); }
    find_fault_tags(param);
    return world;
}

// The GWB properties (temperature: {1,0,0}, composition c: {2,c,0}) at DES point x.
void query(const Param &param, const double *x, const unsigned int props[][3],
           unsigned int nprops, double *values)
{
    void *w = open_world(param);
    const double depth = -x[NDIMS-1];
    // GWB's vertical coordinate is the height above the base; only the depth places features.
    const double z = param.mesh.zlength + x[NDIMS-1];
    try {
#ifdef THREED
        properties_3d(w, x[0], x[1], z, depth, props, nprops, values);
#else
        properties_2d(w, x[0], z, depth, props, nprops, values);
#endif
    }
    catch (const std::exception &e) { gwb_failed(param, e.what()); }
    catch (...) { gwb_failed(param, "unknown error in a query"); }
}

} // anonymous namespace


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


bool gwb_in_fault(const Param &param, const double *x)
{
    open_world(param);
    if (fault_tags.empty()) {
        std::cerr << "Error: ic.weakzone_option 936 needs a fault feature in the Geodynamic World "
                     "Builder world '" << param.ic.gwb_filename << "'; it has none.\n";
        die(EXIT_CONFIG_DATA);
    }
    const unsigned int props[1][3] = {{4, 0, 0}};   // tag; -1 outside every feature
    double tag;
    query(param, x, props, 1, &tag);
    return std::find(fault_tags.begin(), fault_tags.end(), static_cast<int>(tag)) != fault_tags.end();
}


void gwb_close()
{
    if (world) release_world(world);
    world = nullptr;
    fault_tags.clear();
}

#else  // !USEGWB: validate_parameters rejects option 936 in this build

double gwb_temperature(const Param &, const double *)
{
    die(EXIT_UNSUPPORTED_LIB, "Geodynamic World Builder is not compiled in (build with usegwb=1).");
}

int gwb_mattype(const Param &, const double *)
{
    die(EXIT_UNSUPPORTED_LIB, "Geodynamic World Builder is not compiled in (build with usegwb=1).");
}

bool gwb_in_fault(const Param &, const double *)
{
    die(EXIT_UNSUPPORTED_LIB, "Geodynamic World Builder is not compiled in (build with usegwb=1).");
}

void gwb_close() {}

#endif
