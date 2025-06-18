#include "constants.hpp"
#include "parameters.hpp"
#include "markerset.hpp"
#include "utils.hpp"

#include "phasechanges.hpp"

namespace {

    void simple_subduction(const Param& param, const Variables& var, const MarkerSet& ms, int m, int& new_mt, int& hyd_inc)
    {
        auto& hydem = *var.hydrous_elemmarkers;

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
                if (T <= min_serpentine_T && hydem[el][0]) {
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

} // anonymous namespace


void phase_changes(const Param& param, Variables& var)
{
    if (param.mat.nmat == 1 || param.mat.phase_change_option == 0) return; // no phase change        

#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    MarkerSet& ms = *(var.markersets[0]);

    #pragma acc parallel loop async
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
                // std::cerr << "Error: unknown phase_change_option: " << param.mat.phase_change_option << '\n';
                // std::exit(1);
            }            

            if (new_mt != current_mt) {
                ms.set_mattype(m, new_mt);
                --(*var.elemmarkers)[e][current_mt];
                ++(*var.elemmarkers)[e][new_mt];
            }
        }
    }

    #pragma acc wait

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
                const double *eta = ms.get_eta(m);

                // #pragma omp critical(phase_change_simple_subduction)
                {
                    // Add new marker, which has the same coordinate as the dehydrated marker
                    hydms.append_marker(eta, el, 0, 0., 0., 0., 0.);
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
                        std::cerr << "Error: hydrous marker phase change\n";
                        std::exit(12);
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


#ifdef USE_NPROF
    nvtxRangePop();
#endif
}
