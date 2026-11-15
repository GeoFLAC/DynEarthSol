#include <fstream>
#include <iterator>
#include <stdexcept>
#include "earthquake_state.cxx"

static void check(double actual, double expected)
{
    if (std::fabs(actual-expected) > 1e-12*std::max(1.0,std::fabs(expected)))
        throw std::runtime_error("seismic moment differs from G times fault area times slip");
}

int main()
{
    Param param{};
    Variables var;
    constexpr int elements = 2;
    const int nodes = NODES_PER_ELEM;
    array_t coordinate(nodes,0.0), velocity(nodes,0.0);
    conn_t connectivity(elements);
    tensor_t stress(elements,0.0), strain_rate(elements,0.0);
    double_vec temperature(nodes,273.0), pressure(nodes,0.0), lookup;
    // In 2D these are fault lengths and element areas per unit thickness.
    const double fault_area[elements] = {2.0,3.0};
    const double thickness[elements] = {0.25,0.5};
    const double slip[elements] = {0.02,0.03};
    double_vec volume(elements), plastic_increment(elements);
    int_vec2D markers = {{1,0},{0,1}};
    for (int e=0;e<elements;++e) {
        volume[e] = fault_area[e]*thickness[e];
        plastic_increment[e] = slip[e]/(2*thickness[e]);
        for (int j=0;j<NODES_PER_ELEM;++j) connectivity[e][j]=j;
    }
    var.nnode=nodes; var.nelem=elements;
    var.coord=&coordinate; var.vel=&velocity; var.connectivity=&connectivity;
    var.stress=&stress; var.strain_rate=&strain_rate; var.temperature=&temperature;
    var.ppressure=&pressure; var.dppressure=&pressure;
    var.log_table=var.tan_table=var.sin_table=&lookup;
    var.volume=&volume; var.delta_plstrain=&plastic_increment; var.elemmarkers=&markers;
    var.max_vbc_val=1e-9; var.steps=1; var.time=0;
    auto& mat=param.mat;
    mat.nmat=elements; mat.rheol_type=MatProps::rh_ep_rsf;
    mat.shear_modulus={3e10,4e10}; mat.bulk_modulus={5e10};
    mat.visc_exponent={1,1}; mat.visc_coefficient={1,1};
    mat.porosity={0}; mat.biot_coeff={0}; mat.heat_capacity={1};
    mat.therm_cond={0}; mat.fluid_bulk_modulus={1};
    MatProps properties(param,var); var.mat=&properties;
    param.sim.seismic_moment_calculate_output=true;
    param.sim.earthquake_start_factor=2; param.sim.earthquake_end_factor=1;
    param.sim.earthquake_output_step_interval=1;
    const double expected[elements] = {mat.shear_modulus[0]*fault_area[0]*slip[0],
                                      mat.shear_modulus[1]*fault_area[1]*slip[1]};
    EarthquakeState state;
    // Uniform translation and the physical interval cannot change the moment
    // of a prescribed simple-shear plastic slip.
    for (double dt : {0.01,10.0}) for (double speed : {0.001,7.0}) {
        var.dt=dt;
        for (int n=0;n<nodes;++n) velocity[n][0]=speed;
        init_earthquake_state(param,state);
        update_earthquake_tracking(param,var,state);
        for (int e=0;e<elements;++e) check(state.cumulative_moment_by_mat[e],expected[e]);
    }
    init_earthquake_state(param,state);
    for (double& increment : plastic_increment) increment/=4;
    for (int step=0;step<4;++step) {
        ++var.steps;
        update_earthquake_tracking(param,var,state);
    }
    for (int e=0;e<elements;++e) check(state.cumulative_moment_by_mat[e],expected[e]);
    // An unyielding element in the same moving event releases no moment.
    plastic_increment[1]=0;
    init_earthquake_state(param,state);
    update_earthquake_tracking(param,var,state);
    check(state.cumulative_moment_by_mat[0],expected[0]/4);
    check(state.cumulative_moment_by_mat[1],0);
    std::fill(plastic_increment.begin(),plastic_increment.end(),0.0);
    update_earthquake_tracking(param,var,state);
    std::ifstream stream("seismic_moment_magnitude.txt");
    const std::string output((std::istreambuf_iterator<char>(stream)),{});
#ifdef THREED
    if (output.find("M0, N m") == std::string::npos ||
        output.find("Moment magnitude (Mw)") == std::string::npos)
        throw std::runtime_error("3D seismic moment units or magnitude are missing");
#else
    if (output.find("M0, N m/m") == std::string::npos || output.find("Mw") != std::string::npos)
        throw std::runtime_error("2D seismic moment must be per unit length without Mw");
#endif
    std::cout << "PASS simple-shear slip, material sums, subdivision, translation and output units\n";
}
