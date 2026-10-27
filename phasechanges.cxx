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
    // TO ADD A TRANSITION: write a predicate (cite its phase diagram; mark it
    // #pragma acc routine seq, as pass 1 calls it on the device), give it a
    // TransitionKind and its case in transition_rule_applies(), and add one
    // row to phase_transition_candidates() with its mat.has_* switch.  A new
    // source/target phase needs a mat.mattype_* parameter; while the switch
    // is on, check_phase_transition_rules() requires its source and target to
    // be distinct, in-range materials.

    // oceanic crust (basalt) -> eclogite.  Hacker (2003) JGR, Fig. 1/8 (as
    // geoflac names it): low-P boundary trpres and steep high-P boundary
    // trpres2; requires T >= 400 C and depth >= 20 km.
    #pragma acc routine seq
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
    #pragma acc routine seq
    bool antigorite_stable(double press, double T)
    {
        const double Tc = T - 273.0;
        const double trpres  = 2.1e9 + (7.5e9 - 2.1e9) * (Tc - 730.0) / (500.0 - 730.0);
        const double trpres2 = 2.1e9 + (0.2e9 - 2.1e9) * (Tc - 730.0) / (650.0 - 730.0);
        return press < trpres && press > trpres2;
    }

    // serpentinite -> mantle (antigorite dehydration): outside the stability
    // wedge, above 550 C (geoflac's threshold, which it gives no source for).
    #pragma acc routine seq
    bool serpentinite_to_mantle(double depth, double press, double T)
    {
        const double Tc = T - 273.0;
        const double serpentine_temp = 550.0;   // C
        return Tc >= serpentine_temp && !antigorite_stable(press, T);
    }

    // hydrous (chlorite-bearing) mantle -> mantle (chlorite dehydration).
    // Grove et al. (2009) Nature: chlorite-out at Tr(C) = 880 - 3.5e-8*(depth - 62 km)^2.
    #pragma acc routine seq
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
    #pragma acc routine seq
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
    #pragma acc routine seq
    bool sediment_to_metasediment(double depth, double press, double T)
    {
        const double Tc = T - 273.0;
        return Tc >= 250.0 || depth >= 7e3;
    }

    // metamorphic sediment -> schist (dehydration).  Metapelite KFMASH grid,
    // invariant points i5 (710.8 C, 0.880 GPa) & i3 (662.3 C, 0.704 GPa);
    // Wei, Powell & Clarke (2004) J. Metamorph. Geol., Fig. 1.
    #pragma acc routine seq
    bool metasediment_to_schist(double depth, double press, double T)
    {
        const double Tc = T - 273.0;
        const double trpres = 0.88e9 - (0.88e9 - 0.704e9) * (710.8 - Tc) / (710.8 - 662.3);
        return press >= trpres;
    }

    // The predicates by index: device code cannot call through a function
    // pointer, so pass 1 dispatches on this.
    enum class TransitionKind { basalt_to_eclogite, serpentinite_to_mantle, hydrous_mantle_to_mantle,
                                sediment_compaction, sediment_to_metasediment, metasediment_to_schist };

    #pragma acc routine seq
    bool transition_rule_applies(int kind, double depth, double press, double T)
    {
        switch (static_cast<TransitionKind>(kind)) {   // -Wswitch reports a kind without its case
        case TransitionKind::basalt_to_eclogite:       return basalt_to_eclogite(depth, press, T);
        case TransitionKind::serpentinite_to_mantle:   return serpentinite_to_mantle(depth, press, T);
        case TransitionKind::hydrous_mantle_to_mantle: return hydrous_mantle_to_mantle(depth, press, T);
        case TransitionKind::sediment_compaction:      return sediment_compaction(depth, press, T);
        case TransitionKind::sediment_to_metasediment: return sediment_to_metasediment(depth, press, T);
        case TransitionKind::metasediment_to_schist:   return metasediment_to_schist(depth, press, T);
        }
        return false;   // unreachable: every kind is a case
    }

    // One pure-P-T transition: a marker of material `from` becomes `to` when
    // its predicate holds.  (Mantle serpentinization is handled separately in
    // the driver because it additionally needs slab-below geometry.)
    struct PhaseTransition {
        int from;
        int to;
        TransitionKind kind;  // which predicate applies
        const char* source;   // the mat.mattype_* that sets `from`, for check_phase_transition_rules()
        const char* target;   // the mat.mattype_* that sets `to`, for phase_change_keys()
    };

    // Every transition the rule engine knows, whether its mat.has_* switch
    // (`flag`) is on or not. Add a row to introduce a transition.
    struct Candidate { bool on; int from; int to; TransitionKind kind; const char* flag; const char* source; const char* target; };
    std::vector<Candidate> phase_transition_candidates(const Param& param)
    {
        return {
            { param.mat.has_eclogitization,             param.mat.mattype_oceanic_crust,  param.mat.mattype_eclogite,       TransitionKind::basalt_to_eclogite,       "has_eclogitization",             "mattype_oceanic_crust",  "mattype_eclogite" },
            { param.mat.has_serpentinization,           param.mat.mattype_serpentinite,   param.mat.mattype_mantle,         TransitionKind::serpentinite_to_mantle,   "has_serpentinization",           "mattype_serpentinite",   "mattype_mantle" },
            { param.mat.has_hydrous_mantle_dehydration, param.mat.mattype_hydrous_mantle, param.mat.mattype_mantle,         TransitionKind::hydrous_mantle_to_mantle, "has_hydrous_mantle_dehydration", "mattype_hydrous_mantle", "mattype_mantle" },
            { param.mat.has_sediment_metamorphism,      param.mat.mattype_sediment,       param.mat.mattype_sediment_rock,  TransitionKind::sediment_compaction,      "has_sediment_metamorphism",      "mattype_sediment",       "mattype_sediment_rock" },
            { param.mat.has_sediment_metamorphism,      param.mat.mattype_sediment_rock,  param.mat.mattype_metasediment,   TransitionKind::sediment_to_metasediment, "has_sediment_metamorphism",      "mattype_sediment_rock",  "mattype_metasediment" },
            { param.mat.has_sediment_metamorphism,      param.mat.mattype_metasediment,   param.mat.mattype_schist,         TransitionKind::metasediment_to_schist,   "has_sediment_metamorphism",      "mattype_metasediment",   "mattype_schist" },
        };
    }

    // The active rule set: each candidate whose switch is on. Its indices are
    // safe to read per marker only after check_phase_transition_rules()
    // (input), which dies on an active switch whose source or target is unset,
    // out of range, or the same material.
    std::vector<PhaseTransition> build_phase_transition_rules(const Param& param)
    {
        std::vector<PhaseTransition> rules;
        for (const Candidate& c : phase_transition_candidates(param))
            if (c.on)
                rules.push_back({c.from, c.to, c.kind, c.source, c.target});
        return rules;
    }

    // Whether the serpentinization grid counts a material as wet crust: the
    // water-bearing slab material whose dehydration hydrates the mantle above
    // it -- oceanic crust, sediment, and sedimentary rock when its index is set
    // (-1, unset, equals no material). Continental crust, eclogite,
    // metasediment and schist are crust but not wet crust. GeoFLAC master's
    // mantle rule looks below a marker for oceanic crust and sedimentary rock
    // (kocean1, kocean2, ksed2); unconsolidated sediment counts too, as in
    // tan2's geoflac fork.
    #pragma acc routine seq
    bool counts_as_wet_crust(const Param& param, int mt)
    {
        return mt == param.mat.mattype_oceanic_crust || mt == param.mat.mattype_sediment
               || mt == param.mat.mattype_sediment_rock;
    }

    // How many cells pass 2 looks down from a mantle marker's cell: those
    // one, two, ... cell sizes below it, up to `thickness`. Counted by
    // repeated addition, not by division, which can round a boundary case
    // the other way and change which markers convert; capped at num_cells(),
    // more than any column holds, so a huge thickness cannot run away.
    int wet_crust_walk_steps(const VirtualTopology& vt, double thickness)
    {
        int nsteps = 0;
        for (double dist = vt.cell_size(); dist <= thickness && nsteps < vt.num_cells();
                dist += vt.cell_size())
            ++nsteps;
        return nsteps;
    }

    // The wet-crust-margin test (margin >= 0): true iff, at most nsteps cells
    // STRICTLY BELOW `cell0` (the querying marker's build-time label, or its
    // live cell when the label was invalidated -- see the pass-2 call site),
    // a wet-crust-dominant cell (vt.is_dominant(), after the drops between the
    // passes) holds the marker's horizontal position p within its recorded
    // wet-crust extent (+/- margin) -- rejecting wet crust that is in the
    // marker's column but horizontally distant (the slab-bend
    // column-quantization false positive). cell0 itself is skipped -- a mantle marker sitting inside
    // the wet-crust body must not self-trigger; only wet crust strictly below it
    // does. With the margin off, pass 2 reads the whole-column answer from
    // VirtualTopology::cells_with_dominant_below() instead.
    #pragma acc routine seq
    bool wet_crust_dominant_cell_below(const VirtualTopology& vt, int cell0,
                                       const double* p, int nsteps,
                                       double margin)
    {
        int cell = vt.neighbor_cell(cell0, NDIMS-1, -1);   // one cell below, not cell0 itself
        for (int k = 0; cell >= 0 && k < nsteps; ++k) {
            if (vt.is_dominant(cell) && vt.dominant_xext_ok(cell, p, margin))
                return true;
            cell = vt.neighbor_cell(cell, NDIMS-1, -1);
        }
        return false;
    }

    // phase_change_option 3 driver: geometric mantle serpentinization plus the
    // pure-P-T rule table above, in TWO marker passes, each an element loop
    // (OpenMP threads on the CPU, a device kernel on the GPU). Pass 1 applies
    // at most one P-T rule per marker, as FLAC's select case does. Pass 2 then
    // reads each marker's material after pass 1, so a marker can convert twice
    // in one call, where FLAC converts it once; only hydrous mantle can, in a
    // thin P-T sliver near the 150 km depth cap. Each pass iterates by element,
    // so a marker's per-element elemmarkers update is race-free in it.
    //
    // Two passes, not one, because pass 2 reads the grid's dominant-cell
    // flags. Pass 1 can turn a wet-crust marker into a material the grid does not
    // count (eclogite, metasediment), so its cell may no longer be
    // wet-crust-dominant: pass 1 flags the cell, and the host step between the
    // passes drops each flagged cell from the dominant set. A dropped cell is
    // not recounted; it stays non-dominant until the next rebuild. Pass 1's
    // end (the OpenMP join, or the acc wait) makes its flags visible to the
    // host step, and pass 2 starts after that step.
    //
    // The host step then turns the dominant set into a per-cell flag:
    // wet-crust-dominant cell within serpentinization_thickness below. With the
    // margin off, pass 2 reads one flag per mantle marker; with it on, the
    // flag prefilters wet_crust_dominant_cell_below(). The dominant set is
    // computed at each build_serpentinization_grid() and, apart from those
    // drops, stays valid across the phase_changes() calls until the next
    // rebuild, at the next quality check or remesh.
    void subduction_phase_changes(const Param& param, Variables& var, MarkerSet& ms)
    {
        int n = ms.get_nmarkers();
        if (n == 0) return;

        // Pass 1's rule by source material, a table device code can read: one
        // rule per source material (check_phase_transition_rules() refuses two).
        int_vec rule_kind(param.mat.nmat, -1), rule_to(param.mat.nmat, -1);
        for (const PhaseTransition& r : build_phase_transition_rules(param)) {
            rule_kind[r.from] = static_cast<int>(r.kind);
            rule_to[r.from] = r.to;
        }

        // Geometric input for serpentinization: var.vtopo's general marker-cell
        // index (see mesh.hpp), built eagerly after init(), at each quality
        // check and at the end of remeshing, and frozen in between except for
        // invalidate_dominant_cell() (between the passes below). It is rebuilt
        // here, before the next quality check, only when a restart did not
        // bring a usable one: absent (the checkpoint was written without it,
        // or restart() dropped it for another cell size), or without
        // wet-crust extents while the margin feature is on.
        bool serpentinize = uses_serpentinization_grid(param);
        if (serpentinize &&
                (var.vtopo == nullptr ||
                 (param.mat.serpentinization_wet_crust_margin >= 0 &&
                  !var.vtopo->has_xext_data())))
            build_serpentinization_grid(param, var);
        const VirtualTopology* vt = var.vtopo;   // null unless serpentinize
        // dominant_xext_ok() reads the extents unchecked (device code cannot die).
        if (serpentinize && param.mat.serpentinization_wet_crust_margin >= 0 && !vt->has_xext_data())
            die(EXIT_INTERNAL_ASSERT, "phase_change_option 3: the wet-crust margin is on, but the grid has no wet-crust extents");
        // label_cell() reads unchecked: every marker must have its label (see the
        // VirtualTopology doc in mesh.hpp for why no marker is added or removed).
        if (serpentinize && (vt->num_labels() != n || vt->built_last_id() != ms.get_last_id()))
            die(EXIT_INTERNAL_ASSERT, "phase_change_option 3: a marker was added or removed since the serpentinization grid was built");

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

        // Pass 2's look-down distance in cells, its per-cell answer (filled
        // between the passes), and pass 1's per-cell flags.
        int nsteps = serpentinize ? wet_crust_walk_steps(*vt, param.mat.serpentinization_thickness) : 0;
        std::vector<char> wet_crust_below;
        // 1: a wet-crust marker there turned into a material that is not wet crust
        int_vec wet_crust_lost(serpentinize ? vt->num_cells() : 0, 0);
        int* lost = wet_crust_lost.data();   // an atomic write needs a plain scalar target

        // Pass 1: every transition EXCEPT mantle -> serpentinite (that one
        // reads the dominant set, which must first drop the cells whose wet crust
        // this pass converts; see the function doc above).
#ifndef ACC
        #pragma omp parallel for default(none) \
            shared(param, var, ms, rule_kind, rule_to, serpentinize, vt, topo, lost)
#endif
        #pragma acc parallel loop gang vector async
        for (int e=0; e<var.nelem; e++) {
            int nmarkers = (*var.markers_in_elem)[e].size();
            #pragma acc loop seq
            for (int i=0; i<nmarkers; i++) {
                int m = (*var.markers_in_elem)[e][i];
                const int mt = ms.get_mattype(m);
                if (serpentinize && mt == param.mat.mattype_mantle)
                    continue;   // resolved in pass 2

                // The marker's rule first: a material without one needs no
                // P-T evaluation at all.
                const int kind = rule_kind[mt];
                if (kind < 0)
                    continue;

                double Z, press, T, xz[NDIMS];
                ms.get_ZPT(var, m, topo, Z, press, T, xz);  // Z below the current surface
                const double depth = std::max(0.0, -Z);   // m below the local free surface
                if (!transition_rule_applies(kind, depth, press, T))
                    continue;

                const int to = rule_to[mt];   // never mt: build_phase_transition_rules() drops from == to
                ms.set_mattype(m, to);
                --(*var.elemmarkers)[e][mt];
                ++(*var.elemmarkers)[e][to];
                // When a wet-crust marker turns into a material the grid does not
                // count, its cell may no longer be wet-crust-dominant: flag it for
                // the drop between the passes.
                if (serpentinize && counts_as_wet_crust(param, mt) && !counts_as_wet_crust(param, to)) {
                    // The cell whose dominance this marker contributed to --
                    // its BUILD-time cell (the cell the dominance scan
                    // actually counted it in; the live cell can differ once
                    // the marker has drifted).
                    const int cell = vt->label_cell(m);
#ifndef ACC
                    #pragma omp atomic write
#endif
                    #pragma acc atomic write
                    lost[cell] = 1;
                }
            }
        }
        #pragma acc wait

        if (!serpentinize) return;

        // Between the passes, on the host: drop the flagged cells, then build
        // the wet-crust-below map pass 2 reads.
        for (int cell = 0; cell < vt->num_cells(); ++cell)
            if (wet_crust_lost[cell])
                var.vtopo->invalidate_dominant_cell(cell);
        vt->cells_with_dominant_below(nsteps, wet_crust_below);

        // Pass 2: mantle -> serpentinite.
#ifndef ACC
        #pragma omp parallel for default(none) \
            shared(param, var, ms, vt, max_basalt_depth, topo, nsteps, wet_crust_below)
#endif
        #pragma acc parallel loop gang vector async
        for (int e=0; e<var.nelem; e++) {
            int nmarkers = (*var.markers_in_elem)[e].size();
            #pragma acc loop seq
            for (int i=0; i<nmarkers; i++) {
                int m = (*var.markers_in_elem)[e][i];
                if (ms.get_mattype(m) != param.mat.mattype_mantle)
                    continue;

                // The wet-crust-below test reads the marker's BUILD-time cell (its
                // label from the last grid rebuild), so the marker-grid
                // relationship is evaluated in the rebuild-time configuration
                // and stays constant between rebuilds -- a drifting marker
                // cannot flip the answer mid-interval by crossing a lattice
                // line against the frozen dominance data. The margin test still
                // uses the live horizontal position -- drift between rebuilds is
                // far smaller than any useful margin. A marker with no wet crust
                // below skips the P-T evaluation; the tests are pure, so their
                // order cannot change the result.
                const int cell0 = vt->label_cell(m);
                if (!wet_crust_below[cell0])
                    continue;

                double Z, press, T, p[NDIMS];
                ms.get_ZPT(var, m, topo, Z, press, T, p);  // below current surface
                const double depth = std::max(0.0, -Z);
                if (depth > max_basalt_depth || !antigorite_stable(press, T))
                    continue;

                if (param.mat.serpentinization_wet_crust_margin >= 0 &&
                        !wet_crust_dominant_cell_below(*vt, cell0, p, nsteps,
                            param.mat.serpentinization_wet_crust_margin))
                    continue;
                ms.set_mattype(m, param.mat.mattype_serpentinite);
                --(*var.elemmarkers)[e][param.mat.mattype_mantle];
                ++(*var.elemmarkers)[e][param.mat.mattype_serpentinite];
            }
        }
        #pragma acc wait
    }

} // anonymous namespace


void check_phase_transition_rules(const Param& param)
{
    // A switch that is on names materials the rules index per marker, and
    // it must not turn itself off because one of them was left unset or
    // names the other.
    const int nmat = param.mat.nmat;
    for (const Candidate& c : phase_transition_candidates(param)) {
        if (!c.on) continue;
        if (c.from < 0 || c.from >= nmat || c.to < 0 || c.to >= nmat) {
            const std::string msg = std::string("mat.phase_change_option 3: mat.") + c.flag + " is on, but mat."
                + c.source + " (" + std::to_string(c.from) + ") or mat." + c.target + " ("
                + std::to_string(c.to) + ") is unset or not in [0, mat.num_materials); set both, "
                "or turn mat." + c.flag + " off.";
            die(EXIT_CONFIG_VALUE, msg.c_str());
        }
        if (c.from == c.to) {
            const std::string msg = std::string("mat.phase_change_option 3: mat.") + c.flag + " is on, but mat."
                + c.source + " and mat." + c.target + " are both " + std::to_string(c.from)
                + "; give them different indices, or turn mat." + c.flag + " off.";
            die(EXIT_CONFIG_VALUE, msg.c_str());
        }
    }
    // Pass 1 reads one rule per source material, so two rules that share a
    // source cannot both apply. No rule is meant to start from the mantle:
    // with serpentinization on, pass 1 leaves mantle markers to pass 2 and
    // such a rule never fires; with it off, the rule would convert the whole
    // mantle.
    const std::vector<PhaseTransition> rules = build_phase_transition_rules(param);
    for (std::size_t i = 0; i < rules.size(); ++i) {
        for (std::size_t j = i + 1; j < rules.size(); ++j)
            if (rules[i].from == rules[j].from) {
                const std::string msg = std::string("mat.phase_change_option 3: mat.") + rules[i].source
                    + " and mat." + rules[j].source + " are both " + std::to_string(rules[i].from)
                    + ", so only one of the two rules would apply; give them different indices.";
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
    // Surface processes deposit mattype_sediment; were it oceanic crust, the
    // deposits would turn into eclogite.
    if (param.mat.has_eclogitization && param.mat.mattype_sediment == param.mat.mattype_oceanic_crust) {
        const std::string msg = "mat.phase_change_option 3: mat.mattype_sediment equals mat.mattype_oceanic_crust ("
            + std::to_string(param.mat.mattype_oceanic_crust) + ") while mat.has_eclogitization is on, so "
            "sediment would turn into eclogite; give mat.mattype_sediment (default 0) its own index.";
        die(EXIT_CONFIG_VALUE, msg.c_str());
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
    // check_phase_transition_rules() compares every rule source with the
    // mantle, and the grid converts it.
    if (!rules.empty() || uses_serpentinization_grid(param))
        keys.push_back("mattype_mantle");
    if (uses_serpentinization_grid(param)) {
        keys.push_back("mattype_serpentinite");
        keys.push_back("mattype_oceanic_crust");
        keys.push_back("mattype_sediment");   // counts as wet crust (counts_as_wet_crust())
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
    // (zeff = z - z_surf(x)) instead of absolute z. The "wet-crust-below" geometry
    // the grid encodes is then a RELATIVE relationship: a rigid vertical motion
    // of the whole mesh (e.g. the initial isostatic rise) shifts every marker and
    // the surface by the same amount, so every marker's zeff -- and hence its
    // cell -- is unchanged, and the wet-crust-dominant cells stay registered to the
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

    // Bin every marker, not just wet crust: the dominance scan below counts each
    // cell's wet-crust share of ITS OWN population, and pass 2 reads each mantle
    // marker's wet-crust-below flag at its label.
    var.vtopo->build_index(points, n, ms.get_last_id());

    // Precompute which cells are wet-crust-dominant (more than half of their
    // markers count as wet crust, counts_as_wet_crust()) so subduction_phase_changes()
    // pass 2 builds its wet-crust-below map from them instead of recomputing cell
    // composition at every phase_changes() call (every 10 steps).
    // 0.50, not GeoFLAC master's 0.8 (phase_ratio > 0.8 in an element below
    // a mantle marker, which it gives no source for): on these half-element
    // cells 0.8 leaves out the cell straddling the slab top, which cuts the
    // serpentinite band of examples/subduction.cfg by more than half. A cell
    // holds markers_per_element times (cell / element measure) markers:
    // nominally ~2.3 in 2D and ~4.2 in 3D (vt_size doc), so a two-marker cell
    // needs both. Untuned.
    const double wet_crust_fraction_threshold = 0.50;
    const int ncells = var.vtopo->num_cells();
    // Record each dominant cell's wet-crust horizontal extent only when the
    // margin feature is on: the extents are dead weight otherwise, in memory
    // and in the checkpoint (they ride in vtopo_dbl -- see
    // VirtualTopology::serialize).
    const bool want_xext = param.mat.serpentinization_wet_crust_margin >= 0;
    var.vtopo->clear_dominant();
    for (int c = 0; c < ncells; ++c) {
        const VirtualTopology::IdRange ids = var.vtopo->cell_markers(c);
        int wet_crust = 0;
        double xlo[NDIMS-1], xhi[NDIMS-1];
        for (int id : ids)
            if (counts_as_wet_crust(param, ms.get_mattype(id))) {
                if (want_xext) {
                    for (int d = 0; d < NDIMS-1; d++) {
                        const double x = points[id][d];
                        if (wet_crust == 0 || x < xlo[d]) xlo[d] = x;
                        if (wet_crust == 0 || x > xhi[d]) xhi[d] = x;
                    }
                }
                ++wet_crust;
            }
        if (!ids.empty() && static_cast<double>(wet_crust) / ids.size() > wet_crust_fraction_threshold) {
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
        // Its own two marker loops (device kernels in an ACC build), and
        // geometric serpentinization needs the VirtualTopology grid. Timed on
        // its own (rather than folded into a generic phase_changes timer), so
        // its cost shows in the run summary.
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
