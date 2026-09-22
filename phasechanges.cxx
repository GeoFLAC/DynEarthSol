#include <algorithm>
#include <string>
#include <vector>

#include "constants.hpp"
#include "parameters.hpp"
#include "markerset.hpp"
#include "utils.hpp"
#include "mesh.hpp"
#include "geometry.hpp"   // SurfaceTopo: surface-relative depth for the serpentinization grid

#include "phasechanges.hpp"

namespace {

    void simple_subduction(const Param& param, const Variables& var, const MarkerSet& ms, int m, int& new_mt, int& hyd_inc)
    {
        const int mt_mantle = 0;
        const int mt_serpentinized_mantle = 1;
        const int mt_oceanic_crust = 2;
        const int mt_eclogite = 3;
        const int mt_sediment = 4;
        const int mt_schist = 5;
        const int mt_upper_continental_crust = 6;
        const int mt_lower_continental_crust = 7;

        double Z, P, T;
        ms.get_ZPT(param, var, m, Z, P, T);

        int current_mt = new_mt;

        switch (current_mt) {
        case mt_oceanic_crust:
            {
                // basalt -> eclogite
                // Phase diagram from Hacker, 1996, Subduction: Top to Bottom
                const double min_eclogite_T = 500 + 273;
                const double transition_pressure = -0.3e9 + 2.2e6*T;
                const double dehydration_T = 150 + 273;
                if (T > min_eclogite_T && P > transition_pressure) {
                    new_mt = mt_eclogite;
                }
                else if (T > dehydration_T) {
                    hyd_inc = 1;
                }
            }
            break;
        case mt_sediment:
            {
                // sediment -> schist/gneiss
                // from sediment solidus in Nichols et al, 1994, Nature
                const double min_schist_T = 650 + 273;
                const double min_schist_Z = -20e3;
                const double dehydration_T = 150 + 273;
                if (T > min_schist_T && Z < min_schist_Z) {
                    new_mt = mt_schist;
                }
                else if (T > dehydration_T) {
                    hyd_inc = 1;
                }
            }
            break;
        case mt_serpentinized_mantle:
            {
                // serpentinite -> normal mantle
                // Phase diagram taken from Ulmer and Trommsdorff, 1995, Nature
                // Fixed points (730 C, 2.1 GPa) (500 C, 7.5 GPa)
                const double transition_pressure = 2.1e9 + (7.5e9 - 2.1e9) * (T - (730+273)) / (500 - 730);
                const double min_serpentine_T = 550 + 273;
                if (T > min_serpentine_T && P > transition_pressure) {
                    new_mt = mt_mantle;
                    hyd_inc = 1;
                }
            }
            break;
        case mt_mantle:
            {
                const double min_serpentine_T = 550 + 273;
                const int el = ms.get_elem(m);
                // The hydrous marker set exists only with has_hydration_processes;
                // without it no water is present, so the mantle stays dry.
                if (T <= min_serpentine_T && param.control.has_hydration_processes &&
                    (*var.hydrous_elemmarkers)[el][0]) {
                    new_mt = mt_serpentinized_mantle;
                    //hyd_inc = -1;
                }
            }
            break;
        }

    }

    // A template to phase change function
    void custom_phase_change(const Param& param, const Variables& var, const MarkerSet& ms, int m, int& new_mt)
    {
        // Placeholder for custom phase change logic
        // This function should be replaced with actual phase change logic
        // For now, it just returns the current mattype
        double Z, P, T;
        ms.get_ZPT(param, var, m, Z, P, T);

        int current_mt = new_mt;
        // Set new mattype the same as current mattype for now

        switch (current_mt) {
        case 0:
            break;
        case 1:
            break;
        }
    }

    // ===== phase_change_option 3: subduction phase-change rule engine =======
    //
    // A marker converts to a metamorphic product when it crosses a reaction
    // boundary in pressure and temperature.  Ported from geoflac's
    // change_phase routine (github.com/tan2/geoflac, src/change_phase.f90);
    // each predicate cites the source that routine or its history names, or
    // says that it names none.
    //
    // Arguments passed to every predicate:
    //   depth   metres below the CURRENT local free surface (positive down):
    //           get_ZPT() is called with a SurfaceTopo, so depth = z_surf(x) - z
    //           and a rigid vertical motion of the whole mesh (e.g. the initial
    //           isostatic rise) does not change it. Reduces to depth below z = 0
    //           when there is no relief / no top boundary.
    //   press   pressure (Pa) from get_ZPT() -- the mean stress of the marker's
    //           element, compression positive (in 2D the in-plane mean DES
    //           itself uses), so it carries the model's tectonic over- and
    //           underpressure (the source used a fixed-density lithostatic)
    //   T       temperature in Kelvin (converted to Celsius inside each
    //           predicate, the unit of geoflac's thresholds and the cited
    //           diagrams)
    //
    // TO ADD A TRANSITION: write a predicate (cite its phase diagram) and add
    // one row to build_phase_transition_rules() with its mat.has_* switch.
    // A new source/target phase needs a mat.mattype_* parameter; the rule
    // stays inactive unless its switch is on AND that index is a distinct,
    // in-range material.

    // oceanic crust (basalt) -> eclogite.  Hacker (2003) JGR, Fig. 1/8 (as
    // geoflac names it): low-P boundary trpres and steep high-P boundary
    // trpres2; requires T >= 400 C and depth >= 20 km.
    bool basalt_to_eclogite(double depth, double press, double T)
    {
        const double Tc = T - 273.0;
        const double min_eclogite_temp  = 400.0;   // C
        const double min_eclogite_depth = 20e3;    // m
        const double trpres  = -0.3e9 + 2.2e6 * Tc;
        const double trpres2 = 20.1e9 - 0.0375e9 * Tc;
        return Tc >= min_eclogite_temp && depth >= min_eclogite_depth
               && press >= trpres && press >= trpres2;
    }

    // Antigorite (serpentinite) stability wedge, the boundary of both
    // serpentinite rules.  Ulmer & Trommsdorff (1995) Science; fixed points
    // (730 C, 2.1 GPa)-(500 C, 7.5 GPa) and (730 C, 2.1 GPa)-(650 C, 0.2 GPa).
    // Mantle serpentinizes only inside it (pass 2); serpentinite dehydrates
    // outside it (serpentinite_to_mantle). One field for both rules leaves no
    // P-T where both apply, so a marker cannot flip between mantle and
    // serpentinite at fixed P-T. GeoFLAC master's mantle rule takes
    // (670 C, 0.6 GPa) as the low-P point: between the two low-P lines, at
    // 642-730 C, its mantle serpentinizes and the serpentinite dehydrates
    // again at the next call (here even within one call: pass 1, then pass 2).
    bool antigorite_stable(double press, double T)
    {
        const double Tc = T - 273.0;
        const double trpres  = 2.1e9 + (7.5e9 - 2.1e9) * (Tc - 730.0) / (500.0 - 730.0);
        const double trpres2 = 2.1e9 + (0.2e9 - 2.1e9) * (Tc - 730.0) / (650.0 - 730.0);
        return press < trpres && press > trpres2;
    }

    // serpentinite -> mantle (antigorite dehydration): outside the stability
    // wedge, above 550 C (geoflac's threshold, which it gives no source for).
    bool serpentinite_to_mantle(double depth, double press, double T)
    {
        const double Tc = T - 273.0;
        const double serpentine_temp = 550.0;   // C
        return Tc >= serpentine_temp && !antigorite_stable(press, T);
    }

    // hydrous (chlorite-bearing) mantle -> mantle (chlorite dehydration).
    // Grove et al. (2009) Nature: chlorite-out at Tr(C) = 880 - 3.5e-8*(depth - 62 km)^2.
    bool hydrous_mantle_to_mantle(double depth, double press, double T)
    {
        const double Tc = T - 273.0;
        const double trtmpr = 880.0 - 35e-9 * (depth - 62e3) * (depth - 62e3);
        return Tc >= trtmpr;
    }

    // unconsolidated sediment -> sedimentary rock (compaction/lithification).
    // 80 C and 2 km: geoflac's thresholds, which it gives no source for. Both
    // must hold, as in geoflac master until July 2026; geoflac now converts
    // on either. Heuristic, untuned.
    bool sediment_compaction(double depth, double press, double T)
    {
        const double Tc = T - 273.0;
        return Tc >= 80.0 && depth >= 2e3;
    }

    // sedimentary rock -> metamorphic sediment (low-grade metamorphism).
    // 250 C or 7 km: geoflac's thresholds. The geoflac commit that added them
    // cites Montanez & Crossey (2017), "Diagenesis", Encyclopedia of
    // Engineering Geology (Springer), doi:10.1007/978-3-319-39193-9_35-1, for
    // T > 250 C and depth > 7 km; geoflac and this port convert on either.
    bool sediment_to_metasediment(double depth, double press, double T)
    {
        const double Tc = T - 273.0;
        return Tc >= 250.0 || depth >= 7e3;
    }

    // metamorphic sediment -> schist (dehydration).  Metapelite KFMASH grid,
    // invariant points i5 (710.8 C, 0.880 GPa) & i3 (662.3 C, 0.704 GPa);
    // Wei, Powell & Clarke (2004) J. Metamorph. Geol., Fig. 1.
    bool metasediment_to_schist(double depth, double press, double T)
    {
        const double Tc = T - 273.0;
        const double trpres = 0.88e9 - (0.88e9 - 0.704e9) * (710.8 - Tc) / (710.8 - 662.3);
        return press >= trpres;
    }

    // One pure-P-T transition: a marker of material `from` becomes `to` when
    // applies(depth, press, T) is true.  (Mantle serpentinization is handled
    // separately in the driver because it additionally needs slab-below geometry.)
    struct PhaseTransition {
        int from;
        int to;
        bool (*applies)(double depth, double press, double T);
        const char* source;   // the mat.mattype_* that sets `from`, for check_phase_transition_rules()
        const char* target;   // the mat.mattype_* that sets `to`, for phase_change_keys()
    };

    // Build the active rule set from the process switches and configured
    // material indices.  Add a row to introduce a transition; a rule is dropped
    // unless its mat.has_* switch is on AND its source and target are distinct,
    // in-range materials (so an unset mat.mattype_* also keeps it off).
    std::vector<PhaseTransition> build_phase_transition_rules(const Param& param)
    {
        const int nmat = param.mat.nmat;
        struct Candidate { bool on; int from; int to; bool (*applies)(double, double, double); const char* source; const char* target; };
        const Candidate all[] = {
            { param.mat.has_eclogitization,             param.mat.mattype_oceanic_crust,  param.mat.mattype_eclogite,       basalt_to_eclogite,       "mattype_oceanic_crust",  "mattype_eclogite" },
            { param.mat.has_serpentinization,           param.mat.mattype_serpentinite,   param.mat.mattype_mantle,         serpentinite_to_mantle,   "mattype_serpentinite",   "mattype_mantle" },
            { param.mat.has_hydrous_mantle_dehydration, param.mat.mattype_hydrous_mantle, param.mat.mattype_mantle,         hydrous_mantle_to_mantle, "mattype_hydrous_mantle", "mattype_mantle" },
            { param.mat.has_sediment_metamorphism,      param.mat.mattype_sediment,       param.mat.mattype_sediment_rock,  sediment_compaction,      "mattype_sediment",       "mattype_sediment_rock" },
            { param.mat.has_sediment_metamorphism,      param.mat.mattype_sediment_rock,  param.mat.mattype_metasediment,   sediment_to_metasediment, "mattype_sediment_rock",  "mattype_metasediment" },
            { param.mat.has_sediment_metamorphism,      param.mat.mattype_metasediment,   param.mat.mattype_schist,         metasediment_to_schist,   "mattype_metasediment",   "mattype_schist" },
        };
        std::vector<PhaseTransition> rules;
        for (const Candidate& c : all)
            if (c.on && c.from != c.to && c.from >= 0 && c.from < nmat && c.to >= 0 && c.to < nmat)
                rules.push_back({c.from, c.to, c.applies, c.source, c.target});
        return rules;
    }

    // Whether the serpentinization grid counts a material as crust: oceanic
    // crust, sediment, and sedimentary rock when its index is set (-1, unset,
    // equals no material). GeoFLAC master's mantle rule looks below a marker
    // for oceanic crust and sedimentary rock (kocean1, kocean2, ksed2);
    // unconsolidated sediment counts too, as in tan2's geoflac fork.
    bool counts_as_crust(const Param& param, int mt)
    {
        return mt == param.mat.mattype_oceanic_crust || mt == param.mat.mattype_sediment
               || mt == param.mat.mattype_sediment_rock;
    }

    // True iff a crust-dominant cell (as recorded in vt.dominant_cells_ by
    // build_serpentinization_grid()) exists within `thickness` STRICTLY BELOW
    // `cell0` (the querying marker's build-time label, or its live cell when
    // the label was invalidated -- see the pass-2 call site). cell0 itself
    // is skipped -- a mantle marker sitting inside the crust body must not
    // self-trigger; only crust strictly below it does. dominant_cells_ is
    // precomputed once per remesh and kept consistent between rebuilds by
    // invalidate_dominant_cell() in pass 1 (see below), so is_dominant() is a
    // single O(1) hash lookup per walked cell. margin >= 0 additionally
    // requires the marker's horizontal position p to lie within the walked
    // cell's recorded crust extent (+/- margin) -- rejecting crust that is in
    // the marker's column but horizontally distant (the slab-bend
    // column-quantization false positive); margin < 0 keeps whole-column
    // behavior.
    bool crust_dominant_cell_below(const VirtualTopology& vt, int cell0,
                                    const double* p, double thickness,
                                    double margin)
    {
        const double dx = vt.cell_size();
        int cell = vt.neighbor_cell(cell0, NDIMS-1, -1);   // one cell below, not cell0 itself
        for (double dist = dx; cell >= 0 && dist <= thickness; dist += dx) {
            if (vt.is_dominant(cell) &&
                (margin < 0 || vt.dominant_xext_ok(cell, p, margin)))
                return true;
            cell = vt.neighbor_cell(cell, NDIMS-1, -1);
        }
        return false;
    }

    // phase_change_option 3 driver: geometric mantle serpentinization plus the
    // pure-P-T rule table above, in TWO marker passes over the same #pragma
    // omp parallel region (an OpenMP `for` ends with an implicit barrier
    // unless marked `nowait`, so pass 2 starts only once every thread has
    // finished pass 1). Pass 1 applies at most one P-T rule per marker, as
    // FLAC's select case does. Pass 2 then reads each marker's material after
    // pass 1, so a marker can convert twice in one call, where FLAC converts
    // it once; only hydrous mantle can, in a thin P-T sliver near the 150 km
    // depth cap. Both passes iterate by element with the same #pragma omp
    // for, so a marker's per-element elemmarkers update is race-free in each
    // pass.
    //
    // Two passes, not one, because crust_dominant_cell_below() checks
    // var.vtopo->dominant_cells_, which records which cells are crust-dominant.
    // Pass 1 can turn a crust marker into a material the grid does not count
    // (eclogite, metasediment), so its cell may no longer be crust-dominant --
    // it therefore calls
    // invalidate_dominant_cell() immediately so pass 2 sees the updated answer.
    // Without the barrier, pass 2 could observe a cell as still dominant while
    // pass 1 on another thread is in the middle of converting the marker that
    // made it dominant -- a race on dominant_cells_. The implicit barrier
    // between the two #pragma omp for loops ensures all pass-1 invalidations
    // complete before pass 2 reads.
    //
    // Pass 2 calls crust_dominant_cell_below(), an O(1) is_dominant() hash
    // lookup per walked cell. dominant_cells_ is computed at each
    // build_serpentinization_grid() and, apart from pass 1's invalidations,
    // stays valid across the phase_changes() calls until the next rebuild
    // (between remeshes the slab advances only ~dx per tens of thousands of
    // steps, far slower than the 10-step phase_change interval).
    void subduction_phase_changes(const Param& param, Variables& var, MarkerSet& ms)
    {
        int n = ms.get_nmarkers();
        if (n == 0) return;

        std::vector<PhaseTransition> rules = build_phase_transition_rules(param);

        // Geometric input for serpentinization: var.vtopo's general marker-cell
        // index (see mesh.hpp), built eagerly after init() and at the end of
        // remeshing, and frozen in between except for mark_removed() (see
        // markerset.cxx) and invalidate_dominant_cell() (pass 1 below). It is
        // rebuilt here, before the next remesh, only when a restart did not
        // bring a usable one: absent (the checkpoint was written without it),
        // or without crust extents while the margin feature is on.
        bool serpentinize = uses_serpentinization_grid(param);
        if (serpentinize &&
                (var.vtopo == nullptr ||
                 (param.mat.serpentinization_crust_margin >= 0 &&
                  !var.vtopo->has_xext_data())))
            build_serpentinization_grid(param, var);

        // Depth cap on mantle serpentinization, below the stability wedge's own cold-slab
        // extrapolation (> 200 km at 500 C): slab dehydration, the water source, is
        // essentially complete by ~150 km. Heuristic, untuned.
        double max_basalt_depth = 150e3;  // m

        // Surface-relative depth for BOTH the rules' depth argument (get_ZPT) and the
        // grid query coordinate. Built from the current mesh, so it matches the
        // zeff axis build_serpentinization_grid() binned into (identical under a
        // rigid vertical rise). Built for every phase_change_option-3 marker, not
        // just serpentinizing ones, so eclogite/sediment depths are surface-relative
        // too. Inactive (bit-exact fixed-datum) when there are no top nodes.
        SurfaceTopo topo;
        topo.build(param, var, /*elev_only=*/true);

        array_t points(serpentinize ? n : 0);
        if (serpentinize) {
            ms.calculate_marker_coord(var, points);
            #pragma acc wait
            // transform to the grid's (x, zeff) space for crust_dominant_cell_below.
            // Parallel: it runs over all n markers every call, each i writes only
            // its own row, and topo.zeff() is a pure const lookup.
            #pragma omp parallel for default(none) shared(points, topo, n)
            for (int i=0; i<n; i++) {
                double q[NDIMS];
                for (int d=0; d<NDIMS; d++) q[d] = points[i][d];
                points[i][NDIMS-1] = topo.zeff(q);
            }
        }

        #pragma omp parallel default(none) \
            shared(param, var, ms, rules, points, serpentinize, max_basalt_depth, topo)
        {
            // Pass 1: every transition EXCEPT mantle -> serpentinite (that
            // one reads dominant_cells_, so pass 1's invalidations -- from
            // basalt_to_eclogite -- must be complete first; see the function
            // doc above).
            #pragma omp for
            for (int e=0; e<var.nelem; e++) {
                int nmarkers = (*var.markers_in_elem)[e].size();
                for (int i=0; i<nmarkers; i++) {
                    int m = (*var.markers_in_elem)[e][i];
                    const int mt = ms.get_mattype(m);
                    if (serpentinize && mt == param.mat.mattype_mantle)
                        continue;   // resolved in pass 2

                    double Z, press, T;
                    ms.get_ZPT(var, m, topo, Z, press, T);  // Z below the current surface
                    const double depth = std::max(0.0, -Z);   // m below the local free surface

                    int to = -1;
                    for (const PhaseTransition& r : rules)
                        if (r.from == mt) {
                            if (r.applies(depth, press, T)) to = r.to;
                            break;   // at most one rule per source phase
                        }

                    if (to >= 0 && to != mt) {
                        ms.set_mattype(m, to);
                        --(*var.elemmarkers)[e][mt];
                        ++(*var.elemmarkers)[e][to];
                        // When a crust marker turns into a material the grid
                        // does not count, its cell may no longer be
                        // crust-dominant. Invalidate it immediately
                        // so pass 2's is_dominant() queries see the updated
                        // answer. invalidate_dominant_cell() takes mutex_ and
                        // is an O(1) erase -- cheap even under contention.
                        if (serpentinize && counts_as_crust(param, mt) && !counts_as_crust(param, to)) {
                            // Invalidate the cell whose dominance this marker
                            // contributed to -- its BUILD-time cell (the cell
                            // the dominance scan actually counted it in; the
                            // live cell can differ once the marker has
                            // drifted). Live cell_of() only as fallback for
                            // a label mark_removed() dropped (see label_cell()).
                            int cell = var.vtopo->label_cell(m);
                            if (cell < 0) {
                                double pcell[NDIMS];
                                for (int d = 0; d < NDIMS; d++) pcell[d] = points[m][d];
                                cell = var.vtopo->cell_of(pcell);
                            }
                            var.vtopo->invalidate_dominant_cell(cell);
                        }
                    }
                }
            }
            // implicit barrier: all pass-1 writes (set_mattype + invalidations)
            // are visible to every thread before pass 2 starts.

            // Pass 2: mantle -> serpentinite. crust_dominant_cell_below() uses
            // is_dominant(), an O(1) hash lookup into dominant_cells_
            // precomputed at build time.
            if (serpentinize) {
                #pragma omp for
                for (int e=0; e<var.nelem; e++) {
                    int nmarkers = (*var.markers_in_elem)[e].size();
                    for (int i=0; i<nmarkers; i++) {
                        int m = (*var.markers_in_elem)[e][i];
                        if (ms.get_mattype(m) != param.mat.mattype_mantle)
                            continue;

                        double Z, press, T;
                        ms.get_ZPT(var, m, topo, Z, press, T);  // below current surface
                        const double depth = std::max(0.0, -Z);
                        if (depth > max_basalt_depth || !antigorite_stable(press, T))
                            continue;

                        double p[NDIMS];
                        for (int d=0; d<NDIMS; d++)
                            p[d] = points[m][d];
                        // The walk starts from the marker's BUILD-time cell
                        // (its label from the last grid rebuild), so the
                        // marker-grid relationship is evaluated in the
                        // rebuild-time configuration and stays constant
                        // between rebuilds -- a drifting marker cannot flip
                        // the answer mid-interval by crossing a lattice line
                        // against the frozen dominance data. A label
                        // mark_removed() dropped falls back to the live cell
                        // until the next rebuild. The
                        // margin test still uses the live horizontal position
                        // -- drift between rebuilds is far smaller than any
                        // useful margin.
                        int cell0 = var.vtopo->label_cell(m);
                        if (cell0 < 0)
                            cell0 = var.vtopo->cell_of(p);
                        if (crust_dominant_cell_below(*var.vtopo, cell0, p,
                                param.mat.serpentinization_thickness,
                                param.mat.serpentinization_crust_margin)) {
                            ms.set_mattype(m, param.mat.mattype_serpentinite);
                            --(*var.elemmarkers)[e][param.mat.mattype_mantle];
                            ++(*var.elemmarkers)[e][param.mat.mattype_serpentinite];
                        }
                    }
                }
            }
        }
    }

} // anonymous namespace


void check_phase_transition_rules(const Param& param)
{
    // Pass 1 applies the first rule of a marker's material and stops, so a rule
    // that shares its source with an earlier rule would never fire. No rule is
    // meant to start from the mantle: with serpentinization on, pass 1 leaves
    // mantle markers to pass 2 and such a rule never fires; with it off, the
    // rule would convert the whole mantle.
    const std::vector<PhaseTransition> rules = build_phase_transition_rules(param);
    for (std::size_t i = 0; i < rules.size(); ++i) {
        for (std::size_t j = i + 1; j < rules.size(); ++j)
            if (rules[i].from == rules[j].from) {
                const std::string msg = std::string("mat.phase_change_option 3: mat.") + rules[i].source
                    + " and mat." + rules[j].source + " are both " + std::to_string(rules[i].from)
                    + ", so only the first rule would apply; give them different indices.";
                die(EXIT_CONFIG_VALUE, msg.c_str());
            }
        if (rules[i].from == param.mat.mattype_mantle) {
            const std::string msg = std::string("mat.phase_change_option 3: mat.") + rules[i].source
                + " equals mat.mattype_mantle, and no rule may start from the mantle (it would never "
                "apply with serpentinization on, or convert the whole mantle with it off); give it a "
                "different index.";
            die(EXIT_CONFIG_VALUE, msg.c_str());
        }
    }
}


std::vector<const char*> phase_change_keys(const Param& param)
{
    std::vector<const char*> keys;
    const std::vector<PhaseTransition> rules = build_phase_transition_rules(param);
    for (const PhaseTransition& r : rules) {
        keys.push_back(r.source);
        keys.push_back(r.target);
    }
    if (!rules.empty())
        keys.push_back("mattype_mantle");   // check_phase_transition_rules() compares every source with it
    if (uses_serpentinization_grid(param)) {
        keys.push_back("mattype_mantle");
        keys.push_back("mattype_serpentinite");
        keys.push_back("mattype_oceanic_crust");
        keys.push_back("mattype_sediment");   // counts as crust (counts_as_crust())
    }
    return keys;
}


bool uses_serpentinization_grid(const Param& param)
{
    return param.mat.phase_change_option == 3 && param.mat.has_serpentinization;
}


void build_serpentinization_grid(const Param& param, Variables& var)
{
    if (!uses_serpentinization_grid(param)) return;

    delete var.vtopo;
    var.vtopo = nullptr;

    MarkerSet& ms = *(var.markersets[0]);
    const int n = ms.get_nmarkers();
    if (n == 0) return;

    array_t points(n);
    ms.calculate_marker_coord(var, points);   // waits for its own kernel

    // Express the grid's vertical axis as depth below the CURRENT free surface
    // (zeff = z - z_surf(x)) instead of absolute z. The "crust-below" geometry
    // the grid encodes is then a RELATIVE relationship: a rigid vertical motion
    // of the whole mesh (e.g. the initial isostatic rise) shifts every marker and
    // the surface by the same amount, so every marker's zeff -- and hence its
    // cell -- is unchanged, and the crust-dominant cells stay registered to the
    // material that actually sits above/below. Built from the current top-boundary
    // nodes, so it moves with the mesh with no persistent state (SurfaceTopo doc).
    SurfaceTopo topo;
    topo.build(param, var, /*elev_only=*/true);
    for (int i=0; i<n; i++) {
        double q[NDIMS];
        for (int d=0; d<NDIMS; d++) q[d] = points[i][d];
        points[i][NDIMS-1] = topo.zeff(q);   // read horizontal q first, then overwrite z
    }

    // bounding box in the zeff-transformed space (must enclose every binned point)
    double lo[NDIMS], hi[NDIMS];
    for (int d=0; d<NDIMS; d++)
        lo[d] = hi[d] = points[0][d];
    for (int i=1; i<n; i++)
        for (int d=0; d<NDIMS; d++) {
            lo[d] = std::min(lo[d], points[i][d]);
            hi[d] = std::max(hi[d], points[i][d]);
        }

    var.vtopo = new VirtualTopology(lo, hi, param.mesh.resolution * param.mesh.vt_size);

    // Bin every marker, not just crust: the dominance scan below counts each
    // cell's crust share of ITS OWN population, and pass
    // 2 starts each mantle marker's walk from its label.
    var.vtopo->build_index(points, n);

    // Precompute which cells are crust-dominant (more than half of their
    // markers count as crust, counts_as_crust())
    // so subduction_phase_changes() pass 2 can do O(1) is_dominant() lookups
    // instead of recomputing cell composition at every phase_changes() call
    // (every 10 steps).
    // removed_since_rebuild_ is empty right after build_index(), so
    // cell_markers() here is an exact, zero-overhead bucket read.
    // 0.50, not GeoFLAC master's 0.8 (phase_ratio > 0.8 in an element below
    // a mantle marker, which it gives no source for): on these half-element
    // cells 0.8 leaves out the cell straddling the slab top, which cuts the
    // serpentinite band of examples/subduction.cfg by more than half. A cell
    // holds markers_per_element times (cell / element measure) markers:
    // nominally ~2.3 in 2D and ~4.2 in 3D (vt_size doc), so a two-marker cell
    // needs both. Untuned.
    const double crust_fraction_threshold = 0.50;
    const int ncells = var.vtopo->num_cells();
    // Record each dominant cell's crust horizontal extent only when the
    // margin feature is on: the extents are dead weight otherwise, in memory
    // and in the checkpoint (they ride in vtopo_dbl -- see
    // VirtualTopology::serialize).
    const bool want_xext = param.mat.serpentinization_crust_margin >= 0;
    var.vtopo->clear_dominant();
    int_vec ids;
    for (int c = 0; c < ncells; ++c) {
        var.vtopo->cell_markers(c, ids);
        int crust = 0;
        double xlo[NDIMS-1], xhi[NDIMS-1];
        for (int id : ids)
            if (counts_as_crust(param, ms.get_mattype(id))) {
                if (want_xext) {
                    for (int d = 0; d < NDIMS-1; d++) {
                        const double x = points[id][d];
                        if (crust == 0 || x < xlo[d]) xlo[d] = x;
                        if (crust == 0 || x > xhi[d]) xhi[d] = x;
                    }
                }
                ++crust;
            }
        if (!ids.empty() && static_cast<double>(crust) / ids.size() > crust_fraction_threshold) {
            var.vtopo->set_dominant(c);
            if (want_xext)
                var.vtopo->set_dominant_xext(c, xlo, xhi);
        }
    }
    if (want_xext)
        var.vtopo->mark_xext_built();
}


void phase_changes(const Param& param, Variables& var)
{
    if (param.mat.nmat == 1 || param.mat.phase_change_option == 0) return; // no phase change        

#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif

    MarkerSet& ms = *(var.markersets[0]);

    if (param.mat.phase_change_option == 3) {
        // Host-side path: it runs its own marker loop, and geometric
        // serpentinization needs the VirtualTopology grid.  Kept out of the
        // OpenACC-annotated per-marker loop below. Timed on its own (rather
        // than folded into a generic phase_changes timer), so its cost shows
        // in the run summary.
        const int64_t time_tmp = get_nanoseconds();
        subduction_phase_changes(param, var, ms);
        var.func_time.phase_change_time += get_nanoseconds() - time_tmp;
        ++var.nphasechange;
#ifdef NPROF
        nvtxRangePop();
#endif
        return;
    }

    int is_error = 0;

    #pragma acc parallel loop gang vector reduction(+:is_error) async
    for (int e=0; e<var.nelem; ++e) {
        int nmarkers = (*var.markers_in_elem)[e].size();
        #pragma acc loop seq
        for (int i=0; i<nmarkers; ++i) {
            int m = (*var.markers_in_elem)[e][i];
            int current_mt = ms.get_mattype(m);
            int new_mt = current_mt;
            int hyd_inc = 0;


            switch (param.mat.phase_change_option) {
            case 1:
                simple_subduction(param, var, ms, m, new_mt, hyd_inc);
                ms.set_tmp(m, hyd_inc);
                // marker_loop_wrapper(param, var, ms, *var.elemmarkers, simple_subduction);
                break;
            case 101:
                custom_phase_change(param, var, ms, m, new_mt);
                // marker_loop_wrapper(param, var, ms, *var.elemmarkers, custom_phase_change);
                break;
            default:
                is_error++;
            }            

            if (new_mt != current_mt) {
                ms.set_mattype(m, new_mt);
                --(*var.elemmarkers)[e][current_mt];
                ++(*var.elemmarkers)[e][new_mt];
            }
        }
    }

    #pragma acc wait

    if (is_error) {
        // validate_parameters() admits only the options this dispatch implements.
        std::cerr << "Error: unknown phase_change_option: " << param.mat.phase_change_option << '\n';
        die(EXIT_INTERNAL_UNREACHABLE);
    }

    if (param.control.has_hydration_processes) {
        MarkerSet& hydms = *var.markersets[var.hydrous_marker_index];
        auto& hydem = *var.hydrous_elemmarkers;

        for (int e=0; e<var.nelem; ++e) {
            int nmarkers = (*var.markers_in_elem)[e].size();
            for (int i=0; i<nmarkers; ++i) {
                int m = (*var.markers_in_elem)[e][i];

                if( ms.get_tmp(m) == 0) continue; 

                // Dehydration metamorphism, hydrous marker is released.
                const int el = ms.get_elem(m);
                ConstShapefnAccessor eta = ms.get_eta(m);

                // #pragma omp critical(phase_change_simple_subduction)
                {
                    // Add new marker, which has the same coordinate as the dehydrated marker
                    hydms.append_marker(eta, el, 0, 0., 0., 0., 0., 0.);
                    ++hydem[el][0];
                }

                /*** Disable hyd marker deletion
                else if (param.control.has_hydration_processes && hyd_inc == -1) {
                    const int el = ms.get_elem(m);
                    // Find the hydrous marker belong to el
                    int mh;
                    for (mh=0; mh<hydms.get_nmarkers(); mh++) {
                        if (hydms.get_elem(mh) == el) break;
                    }
                    if (mh >= hydms.get_nmarkers()) {
                        die(EXIT_RUNTIME_LOOKUP, "hydrous marker phase change");
                    }
                    #pragma omp critical(phase_change_simple_subduction)
                    {
                        // delete the marker
                        hydms.remove_marker(mh);
                        --hydem[el][0];
                    }
                }
                */
            }
        }
    }


#ifdef NPROF
    nvtxRangePop();
#endif
}
