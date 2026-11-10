#include <cassert>
#include <cmath>
#include <cstdio>
#include "../../remeshing.cxx"

int main()
{
    Param param{};
    Variables var;
    var.nnode = 6;
    array_t coordinate(6), reference(6), normals(nbdrytypes, 0.0);
    uint_vec flags = {BOUNDX0, BOUNDZ0, BOUNDX0|BOUNDZ0, BOUNDZ1, BOUNDN0, BOUNDX1};
    var.coord = &coordinate;
    var.coord0 = &reference;
    var.bnormals = &normals;
    var.bcflag = &flags;
    param.mesh.remeshing_reset_boundary_mask = BOUNDX0|BOUNDZ0|BOUNDN0;
    const double original[6][2] = {{-.1,.4},{.3,-1.1},{-.1,-1.1},{.4,.05},{.9,.6},{1.1,.2}};
    const double initial[6][2] = {{0,.2},{.2,-1},{0,-1},{.2,0},{.5,.5},{1,0}};
    const double restored[6][2] = {{0,.4},{.3,-1},{0,-1},{.4,.05},{.65,.35},{1.1,.2}};
    for (int n = 0; n < 6; ++n)
        for (int d = 0; d < NDIMS; ++d) {
            coordinate[n][d] = original[n][d];
            reference[n][d] = initial[n][d];
        }
    normals[iboundn0][0] = normals[iboundn0][1] = 1/std::sqrt(2.0);
    assert(std::abs(selected_boundary_reference_distortion(param, var, 0)-.1) < 1e-14);
    assert(selected_boundary_reference_distortion(param, var, 3) == 0);
    double_vec packed;
    coordinate.pack_to(packed);
    restore_selected_boundary_reference(param, var, flags, packed.data());
    for (int n = 0; n < 6; ++n)
        for (int d = 0; d < NDIMS; ++d)
            assert(std::abs(packed[n*NDIMS+d]-restored[n][d]) < 1e-14);
    reset_selected_boundary_coord0(param, var);
    #pragma acc wait
    for (int n = 0; n < 6; ++n)
        assert(selected_boundary_reference_distortion(param, var, n) < 1e-14);
    assert(reference[0][1] == initial[0][1]);
    assert(reference[1][0] == initial[1][0]);
    assert(reference[3][0] == initial[3][0] && reference[3][1] == initial[3][1]);
    assert(reference[5][0] == initial[5][0] && reference[5][1] == initial[5][1]);
    std::puts("PASS: selected axes, oblique normals, corners and tangential reference components");
}
