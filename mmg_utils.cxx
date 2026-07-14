#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
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

namespace {

// Merge duplicate facets (same node set) into one entry, OR-ing their refs. MMG emits each facet
// of a REQUIRED element as an extra ref=0 copy; fed back, MMG keeps one copy at random, and when
// the ref=0 copy wins, nodes inserted on that facet get bcflag 0 and flatten_* never restores them
// (sagging bottom). Refs are BOUND* bitmasks and distinct facets never share a node set, so the
// OR-merge keeps the exact constraint set with a deterministic ref.
void dedup_facets(int &nseg, std::vector<int> &seg, std::vector<int> &segflag)
{
    std::map<std::array<int, NODES_PER_FACET>, int> first_at;   // facet key -> surviving index
    int ns = 0;
    for (int s = 0; s < nseg; ++s) {
        std::array<int, NODES_PER_FACET> key;
        for (int k = 0; k < NODES_PER_FACET; ++k) key[k] = seg[(std::size_t)s*NODES_PER_FACET + k];
        std::sort(key.begin(), key.end());
        auto it = first_at.find(key);
        if (it != first_at.end()) {
            segflag[it->second] |= segflag[s];   // duplicate: merge ref into the first copy
            continue;
        }
        first_at[key] = ns;
        for (int k = 0; k < NODES_PER_FACET; ++k)
            seg[(std::size_t)ns*NODES_PER_FACET + k] = seg[(std::size_t)s*NODES_PER_FACET + k];
        segflag[ns] = segflag[s];
        ++ns;
    }
    nseg = ns;
    seg.resize((std::size_t)ns * NODES_PER_FACET);
    segflag.resize(ns);
}

// Drop ref=0 facets. Every boundary facet carries a BOUND* ref (the OR-merge above guarantees it),
// so a ref=0 facet is interior: on output, an ECHO (MMG re-emits each facet of a required element
// with ref 0); on input at mesh init, an internal segment of the initial mesh (refined-zone box,
// .poly facet). MMG ignores both as constraints (it drops a ref=0 facet between elements of equal
// ref), so the filter never changes what MMG builds. It keeps the INVARIANT that var.segment holds
// boundary facets only, which collapse_short_boundary_segments_2d, check_mesh_tangle and the 2D
// interface-edge append below rely on: they treat every entry as boundary. A Triangle/TetGen-built
// initial mesh gets the same from discard_internal_segments, so validate_parameters rejects
// mesh.is_discarding_internal_segments = no under MMG.
void drop_unflagged_facets(int &nseg, std::vector<int> &seg, std::vector<int> &segflag)
{
    int ns = 0;
    for (int s = 0; s < nseg; ++s) {
        if (segflag[s] == 0) continue;
        for (int k = 0; k < NODES_PER_FACET; ++k)
            seg[(std::size_t)ns*NODES_PER_FACET + k] = seg[(std::size_t)s*NODES_PER_FACET + k];
        segflag[ns] = segflag[s];
        ++ns;
    }
    nseg = ns;
    seg.resize((std::size_t)ns * NODES_PER_FACET);
    segflag.resize(ns);
}

} // anonymous namespace

void mmg_adapt(const Mesh &mesh, const MMGInput &in, MMGOutput &out)
{
    // The size band of the contract (mmg_utils.hpp); validate_parameters caps the ratio so that
    // hmin < hmax.
    const double res   = mesh.resolution;
    const double hmax  = std::pow(mesh.largest_size,  1.0 / NDIMS) * res;
    const double hmin  = std::pow(mesh.smallest_size, 1.0 / NDIMS) * res
                         * std::pow(mesh.mmg_remesh_defensive_tiny_ratio, 1.0 / NDIMS);
    const double hausd = mesh.mmg_hausd_factor * res;

    // A required element is carried verbatim, so each of its nodes must be required too.
    for (int e = 0; in.required_elem && e < in.nelem; ++e) {
        if (!in.required_elem[e]) continue;
        const int *c = in.conn + (std::size_t)e * NODES_PER_ELEM;
        for (int k = 0; k < NODES_PER_ELEM; ++k)
            if (!in.required_node || !in.required_node[c[k]])
                die(EXIT_INTERNAL_ASSERT, "mmg_adapt: a required element has a non-required node");
    }

    // 1-indexed copies for MMG (never mutate the caller's buffers). segref is a mutable
    // copy because MMG's setters take non-const int*.
    std::vector<int> conn1(in.conn, in.conn + (std::size_t)in.nelem * NODES_PER_ELEM);
    for (int &v : conn1) ++v;
    std::vector<int> seg1(in.seg, in.seg + (std::size_t)in.nseg * NODES_PER_FACET);
    std::vector<int> segref(in.segflag, in.segflag + in.nseg);
    // An Exodus mesh (meshing_option 95, 3D only) lists a face once per side set, so a face in two
    // side sets can arrive twice with different refs (a collapse before MMG keeps only the first:
    // rebuild_collapsed_segments); Triangle/TetGen never repeat a facet.
    int nseg_in = in.nseg;
    dedup_facets(nseg_in, seg1, segref);
    drop_unflagged_facets(nseg_in, seg1, segref);   // initial-mesh internal segments
    for (int &v : seg1) ++v;

    MMG5_pMesh mmgMesh = NULL;
    MMG5_pSol  mmgSol  = NULL;

#ifdef THREED
    // MMG3D drops an interior triangle between tetrahedra of equal reference (none carry a ref
    // here), so a required interface triangle would be discarded: front tracking is 2D only.
    if (in.req_edges) die(EXIT_INTERNAL_ASSERT, "mmg_adapt: interface edges are 2D only");
    MMG3D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh,
                    MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end);

    MMG_OK(MMG3D_Set_meshSize(mmgMesh, in.nnode, in.nelem, 0, nseg_in, 0, 0));
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
    // MMG's ref=0 echoes of required-element facets (dedup_facets, drop_unflagged_facets)
    dedup_facets(out.nseg, out.seg, out.segflag);
    drop_unflagged_facets(out.nseg, out.seg, out.segflag);

    MMG3D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh,
                   MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end);

#else // 2D

    MMG2D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh,
                    MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end);

    // req_edges go after the nseg_in boundary edges with a sentinel ref (see MMGInput).
    const int IFACE_EDGE_REF = 1 << 24;   // distinct from every BOUND* flag combination
    const int n_iface = (in.req_edges && in.n_req_edges > 0) ? in.n_req_edges : 0;
    for (int k = 0; k < n_iface; ++k) {
        seg1.push_back(in.req_edges[2*k]     + 1);
        seg1.push_back(in.req_edges[2*k + 1] + 1);
        segref.push_back(IFACE_EDGE_REF);
    }
    const int nseg_all = nseg_in + n_iface;

    MMG_OK(MMG2D_Set_meshSize(mmgMesh, in.nnode, in.nelem, 0, nseg_all));
    MMG_OK(MMG2D_Set_vertices(mmgMesh, const_cast<double*>(in.coord), NULL));
    MMG_OK(MMG2D_Set_triangles(mmgMesh, conn1.data(), NULL));
    MMG_OK(MMG2D_Set_edges(mmgMesh, seg1.data(), segref.data()));
    for (int k = 0; k < n_iface; ++k)   // 1-based; interface edges are the last n_iface entries
        MMG_OK(MMG2D_Set_requiredEdge(mmgMesh, nseg_in + k + 1));

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
    if (n_iface) {   // drop the interface edges (internal constraints) from the boundary segments
        int w = 0;
        for (int i = 0; i < out.nseg; ++i) {
            if (out.segflag[i] == IFACE_EDGE_REF) continue;
            if (w != i) {
                for (int j = 0; j < NODES_PER_FACET; ++j)
                    out.seg[(std::size_t)w*NODES_PER_FACET + j] = out.seg[(std::size_t)i*NODES_PER_FACET + j];
                out.segflag[w] = out.segflag[i];
            }
            ++w;
        }
        out.nseg = w;
        out.seg.resize((std::size_t)out.nseg * NODES_PER_FACET);
        out.segflag.resize(out.nseg);
    }
    // MMG's ref=0 echoes of required-element facets (dedup_facets, drop_unflagged_facets)
    dedup_facets(out.nseg, out.seg, out.segflag);
    drop_unflagged_facets(out.nseg, out.seg, out.segflag);

    MMG2D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh,
                   MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end);

#endif // THREED
}

#undef MMG_OK

#endif // USEMMG
