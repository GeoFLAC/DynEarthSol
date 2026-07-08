#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "constants.hpp"
#include "parameters.hpp"
#include "utils.hpp"
#include "mmg_utils.hpp"

#ifdef USEMMG

#ifdef THREED
#include "mmg/mmg3d/libmmg3d.h"
#else
#include "mmg/mmg2d/libmmg2d.h"
#endif

// Die with a message if an MMG call does not return success (1).
#define MMG_OK(call) do { \
    if ((call) != 1) { \
        std::fprintf(stderr, "Error: MMG call failed: %s\n", #call); \
        die(EXIT_MESH_MMG); \
    } } while (0)


void mmg_adapt(const Mesh &mesh, const MMGInput &in, MMGOutput &out)
{
    // Element-size band from the unified largest_size/smallest_size convention:
    // an equilateral element of the target volume has edge = resolution * size^(1/NDIMS).
    const double res   = mesh.resolution;
    const double hmax  = std::pow(mesh.largest_size,  1.0 / NDIMS) * res;
    const double hmin  = std::pow(mesh.smallest_size, 1.0 / NDIMS) * res;
    const double hausd = mesh.mmg_hausd_factor * res;

    // 1-indexed copies for MMG (never mutate the caller's buffers). segref is a mutable
    // copy because MMG's setters take non-const int*.
    std::vector<int> conn1(in.conn, in.conn + (std::size_t)in.nelem * NODES_PER_ELEM);
    for (int &v : conn1) ++v;
    std::vector<int> seg1(in.seg, in.seg + (std::size_t)in.nseg * NODES_PER_FACET);
    for (int &v : seg1) ++v;
    std::vector<int> segref(in.segflag, in.segflag + in.nseg);

    MMG5_pMesh mmgMesh = NULL;
    MMG5_pSol  mmgSol  = NULL;

#ifdef THREED
    MMG3D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh,
                    MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end);

    MMG_OK(MMG3D_Set_meshSize(mmgMesh, in.nnode, in.nelem, 0, in.nseg, 0, 0));
    MMG_OK(MMG3D_Set_vertices(mmgMesh, const_cast<double*>(in.coord), NULL));
    MMG_OK(MMG3D_Set_tetrahedra(mmgMesh, conn1.data(), NULL));
    MMG_OK(MMG3D_Set_triangles(mmgMesh, seg1.data(), segref.data()));

    MMG_OK(MMG3D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, in.nnode, MMG5_Scalar));
    MMG_OK(MMG3D_Set_scalarSols(mmgSol, const_cast<double*>(in.metric)));

    MMG_OK(MMG3D_Chk_meshData(mmgMesh, mmgSol));

    if (in.required_node)
        for (int n = 0; n < in.nnode; ++n)
            if (in.required_node[n]) MMG_OK(MMG3D_Set_requiredVertex(mmgMesh, n + 1));
    if (in.required_elem)
        for (int e = 0; e < in.nelem; ++e)
            if (in.required_elem[e]) MMG_OK(MMG3D_Set_requiredTetrahedron(mmgMesh, e + 1));

    MMG_OK(MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_optim,   0));
    MMG_OK(MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_verbose, mesh.mmg_verbose));
    MMG_OK(MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_debug,   mesh.mmg_debug));
    MMG_OK(MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmax,  hmax));
    MMG_OK(MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmin,  hmin));
    MMG_OK(MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hausd, hausd));

    const int ier = MMG3D_mmg3dlib(mmgMesh, mmgSol);
    if (ier == MMG5_STRONGFAILURE) {
        std::fprintf(stdout, "BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
        die(EXIT_MESH_MMG);
    } else if (ier == MMG5_LOWFAILURE) {
        std::fprintf(stdout, "BAD ENDING OF MMG3DLIB\n");
        // Init tolerates a low failure (mesh is saved but imperfect); remesh treats it as fatal.
        if (!in.tolerate_low_failure) die(EXIT_MESH_MMG);
    }

    int na;
    MMG_OK(MMG3D_Get_meshSize(mmgMesh, &out.nnode, &out.nelem, NULL, &out.nseg, NULL, &na));
    out.coord.resize((std::size_t)out.nnode * NDIMS);
    out.conn.resize((std::size_t)out.nelem * NODES_PER_ELEM);
    out.seg.resize((std::size_t)out.nseg * NODES_PER_FACET);
    out.segflag.resize(out.nseg);

    for (int i = 0; i < out.nnode; ++i)
        MMG_OK(MMG3D_Get_vertex(mmgMesh, &out.coord[i*NDIMS+0], &out.coord[i*NDIMS+1],
                                &out.coord[i*NDIMS+2], NULL, NULL, NULL));
    for (int i = 0; i < out.nelem; ++i) {
        int *c = &out.conn[(std::size_t)i*NODES_PER_ELEM];
        MMG_OK(MMG3D_Get_tetrahedron(mmgMesh, &c[0], &c[1], &c[2], &c[3], NULL, NULL));
        for (int j = 0; j < NODES_PER_ELEM; ++j) c[j] -= 1;
    }
    for (int i = 0; i < out.nseg; ++i) {
        int *s = &out.seg[(std::size_t)i*NODES_PER_FACET];
        MMG_OK(MMG3D_Get_triangle(mmgMesh, &s[0], &s[1], &s[2], &out.segflag[i], NULL));
        for (int j = 0; j < NODES_PER_FACET; ++j) s[j] -= 1;
    }

    MMG3D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh,
                   MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end);

#else // 2D

    MMG2D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh,
                    MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end);

    MMG_OK(MMG2D_Set_meshSize(mmgMesh, in.nnode, in.nelem, 0, in.nseg));
    MMG_OK(MMG2D_Set_vertices(mmgMesh, const_cast<double*>(in.coord), NULL));
    MMG_OK(MMG2D_Set_triangles(mmgMesh, conn1.data(), NULL));
    MMG_OK(MMG2D_Set_edges(mmgMesh, seg1.data(), segref.data()));

    MMG_OK(MMG2D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, in.nnode, MMG5_Scalar));
    MMG_OK(MMG2D_Set_scalarSols(mmgSol, const_cast<double*>(in.metric)));

    MMG_OK(MMG2D_Chk_meshData(mmgMesh, mmgSol));

    if (in.required_node)
        for (int n = 0; n < in.nnode; ++n)
            if (in.required_node[n]) MMG_OK(MMG2D_Set_requiredVertex(mmgMesh, n + 1));
    if (in.required_elem)
        for (int e = 0; e < in.nelem; ++e)
            if (in.required_elem[e]) MMG_OK(MMG2D_Set_requiredTriangle(mmgMesh, e + 1));

    MMG_OK(MMG2D_Set_iparameter(mmgMesh, mmgSol, MMG2D_IPARAM_optim,   0));
    MMG_OK(MMG2D_Set_iparameter(mmgMesh, mmgSol, MMG2D_IPARAM_verbose, mesh.mmg_verbose));
    MMG_OK(MMG2D_Set_iparameter(mmgMesh, mmgSol, MMG2D_IPARAM_debug,   mesh.mmg_debug));
    MMG_OK(MMG2D_Set_dparameter(mmgMesh, mmgSol, MMG2D_DPARAM_hmax,  hmax));
    MMG_OK(MMG2D_Set_dparameter(mmgMesh, mmgSol, MMG2D_DPARAM_hmin,  hmin));
    MMG_OK(MMG2D_Set_dparameter(mmgMesh, mmgSol, MMG2D_DPARAM_hausd, hausd));

    const int ier = MMG2D_mmg2dlib(mmgMesh, mmgSol);
    if (ier == MMG5_STRONGFAILURE) {
        std::fprintf(stdout, "BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH\n");
        die(EXIT_MESH_MMG);
    } else if (ier == MMG5_LOWFAILURE) {
        std::fprintf(stdout, "BAD ENDING OF MMG2DLIB\n");
        // Init tolerates a low failure (mesh is saved but imperfect); remesh treats it as fatal.
        if (!in.tolerate_low_failure) die(EXIT_MESH_MMG);
    }

    MMG_OK(MMG2D_Get_meshSize(mmgMesh, &out.nnode, &out.nelem, NULL, &out.nseg));
    out.coord.resize((std::size_t)out.nnode * NDIMS);
    out.conn.resize((std::size_t)out.nelem * NODES_PER_ELEM);
    out.seg.resize((std::size_t)out.nseg * NODES_PER_FACET);
    out.segflag.resize(out.nseg);

    for (int i = 0; i < out.nnode; ++i)
        MMG_OK(MMG2D_Get_vertex(mmgMesh, &out.coord[i*NDIMS+0], &out.coord[i*NDIMS+1],
                                NULL, NULL, NULL));
    for (int i = 0; i < out.nelem; ++i) {
        int *c = &out.conn[(std::size_t)i*NODES_PER_ELEM];
        MMG_OK(MMG2D_Get_triangle(mmgMesh, &c[0], &c[1], &c[2], NULL, NULL));
        for (int j = 0; j < NODES_PER_ELEM; ++j) c[j] -= 1;
    }
    for (int i = 0; i < out.nseg; ++i) {
        int *s = &out.seg[(std::size_t)i*NODES_PER_FACET];
        MMG_OK(MMG2D_Get_edge(mmgMesh, &s[0], &s[1], &out.segflag[i], NULL, NULL));
        for (int j = 0; j < NODES_PER_FACET; ++j) s[j] -= 1;
    }

    MMG2D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh,
                   MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end);

#endif // THREED
}

#undef MMG_OK

#endif // USEMMG
