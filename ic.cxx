#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"

#include "ic.hpp"


double get_prem_pressure(double depth)
{
    // reference pressure profile from isotropic PREM model
    const int nlayers = 46;
    const double ref_depth[] = { 0e3,    3e3,    15e3,   24.4e3, 40e3,
                                 60e3,   80e3,   115e3,  150e3,  185e3,
                                 220e3,  265e3,  310e3,  355e3,  400e3,
                                 450e3,  500e3,  550e3,  600e3,  635e3,
                                 670e3,  721e3,  771e3,  871e3,  971e3,
                                 1071e3, 1171e3, 1271e3, 1371e3, 1471e3,
                                 1571e3, 1671e3, 1771e3, 1871e3, 1971e3,
                                 2071e3, 2171e3, 2271e3, 2371e3, 2471e3,
                                 2571e3, 2671e3, 2741e3, 2771e3, 2871e3,
                                 2891e3 };

    // pressure in PREM table is given in kilobar, convert to 10^8 Pa
    const double ref_pressure[] = { 0e8,      0.3e8,    3.3e8,    6.0e8,    11.2e8,
                                    17.8e8,   24.5e8,   36.1e8,   47.8e8,   59.4e8,
                                    71.1e8,   86.4e8,   102.0e8,  117.7e8,  133.5e8,
                                    152.2e8,  171.3e8,  190.7e8,  210.4e8,  224.3e8,
                                    238.3e8,  260.7e8,  282.9e8,  327.6e8,  372.8e8,
                                    418.6e8,  464.8e8,  511.6e8,  558.9e8,  606.8e8,
                                    655.2e8,  704.1e8,  753.5e8,  803.6e8,  854.3e8,
                                    905.6e8,  957.6e8,  1010.3e8, 1063.8e8, 1118.2e8,
                                    1173.4e8, 1229.7e8, 1269.7e8, 1287.0e8, 1345.6e8,
                                    1357.5e8 };

    // PREM model doesn't work if depth is above sea level, always returns 0 pressure
    if (depth <= 0) return 0;

    int n;
    for (n=1; n<nlayers; n++) {
        if (depth <= ref_depth[n]) break;
    }

    // linear interpolation
    double pressure = ref_pressure[n-1] + depth *
        (ref_pressure[n] - ref_pressure[n-1]) / (ref_depth[n] - ref_depth[n-1]);

    return pressure;
}


void initial_stress_state(const Param &param, const Variables &var,
                          tensor_t &stress, tensor_t &strain,
                          double &compensation_pressure)
{
    if (param.control.gravity == 0) {
        compensation_pressure = 0;
        return;
    }

    // lithostatic condition for stress and strain
    double rho = var.mat->rho(0);
    double ks = var.mat->bulkm(0);

    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        double zcenter = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            zcenter += (*var.coord)[conn[i]][NDIMS-1];
        }
        zcenter /= NODES_PER_ELEM;

        if (param.control.ref_pressure_option == 1) {
            ks = var.mat->bulkm(e);
            double p = get_prem_pressure(-zcenter);
            for (int i=0; i<NDIMS; ++i) {
                stress[e][i] = -p;
                strain[e][i] = -p / ks / NDIMS;
            }
        }
        else if (param.control.ref_pressure_option == 0) {
            for (int i=0; i<NDIMS; ++i) {
                stress[e][i] = rho * param.control.gravity * zcenter;
                strain[e][i] = rho * param.control.gravity * zcenter / ks / NDIMS;
            }
        }
    }

    switch (param.control.ref_pressure_option) {
    case 0:
        compensation_pressure = rho * param.control.gravity * param.mesh.zlength;
        break;
    case 1:
        compensation_pressure = get_prem_pressure(param.mesh.zlength);
        break;
    default:
        std::cerr << "Error: unknown option for control.ref_pressure_option: " << param.control.ref_pressure_option << '\n';
        std::exit(1);
    }
}


void initial_weak_zone(const Param &param, const Variables &var,
                       double_vec &plstrain)
{
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        // the coordinate of the center of this element
        double center[3] = {0,0,0};
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            for (int d=0; d<NDIMS; ++d) {
                center[d] += (*var.coord)[conn[i]][d];
            }
        }
        for (int d=0; d<NDIMS; ++d) {
            center[d] /= NODES_PER_ELEM;
        }
        // TODO: adding different types of weak zone
        const double d = param.mesh.resolution;
        if (std::fabs(center[0] - param.mesh.xlength * 0.5) < 2*d &&
            std::fabs(center[NDIMS-1] + param.mesh.zlength) < 1.5*d)
            plstrain[e] = 0.1;
    }
}


void initial_temperature(const Param &param, const Variables &var,
                         double_vec &temperature)
{
    const double oceanic_plate_age = 60e6 * YEAR2SEC;
    const double diffusivity = 1e-6;

    for (int i=0; i<var.nnode; ++i) {
        double w = -(*var.coord)[i][NDIMS-1] / std::sqrt(4 * diffusivity * oceanic_plate_age);
        temperature[i] = param.bc.surface_temperature +
            (param.bc.mantle_temperature - param.bc.surface_temperature) * std::erf(w);
    }
}


