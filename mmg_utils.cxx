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

// Merge duplicate boundary facets (same node set) into one entry, OR-ing their refs. MMG emits
// the facets of REQUIRED elements as extra ref=0 copies of the input segments; fed back, MMG's
// facet hash keeps one copy at random, and when the ref=0 copy wins the boundary ref is lost:
// nodes inserted on that facet get bcflag=0 and flatten_* never restores them (sagging bottom).
// Refs are BOUND* bitmasks and distinct facets never share a node set, so OR-merging preserves
// the exact constraint set with a deterministic ref.
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

// Drop ref=0 facets. With mesh.is_discarding_internal_segments every REAL boundary facet carries
// a BOUND* ref (the OR-merge above guarantees it), so a remaining ref=0 facet is an interior ECHO
// of a required element's edges. Echoes add no protection within an adapt (requiredness does
// that) but constrain the NEXT remeshes' interior and ratchet remesh over remesh (nseg grew ~35x
// over a long run, output quality degraded, then two echoes crossed: tangled mesh). Applied to
// MMG's output, so var.segment is the real boundary only.
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

// Drop the facets carrying `ref` (the interface facets appended on input come back with it).
static void drop_facets_with_ref(MMGOutput &out, int ref)
{
    int w = 0;
    for (int i = 0; i < out.nseg; ++i) {
        if (out.segflag[i] == ref) continue;
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

} // anonymous namespace

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
    std::vector<int> segref(in.segflag, in.segflag + in.nseg);
    // Defensive input dedup: a checkpoint or an output written before the readback dedup below
    // may still carry duplicate facets with conflicting refs; never hand those to MMG.
    int nseg_in = in.nseg;
    dedup_facets(nseg_in, seg1, segref);
    if (mesh.is_discarding_internal_segments)
        drop_unflagged_facets(nseg_in, seg1, segref);   // stale interior echoes (see above)
    for (int &v : seg1) ++v;
    // Material-interface front-tracking (option 13): append the interface facets after the boundary
    // facets with a sentinel ref, mark them MMG-required below, and drop them from the returned
    // segments (they are internal, not a boundary).
    const int IFACE_REF = 1 << 24;   // distinct from every BOUND* flag combination
    const int n_iface = (in.req_facets && in.n_req_facets > 0) ? in.n_req_facets : 0;
    for (int k = 0; k < n_iface; ++k) {
        for (int j = 0; j < NODES_PER_FACET; ++j)
            seg1.push_back(in.req_facets[NODES_PER_FACET*k + j] + 1);
        segref.push_back(IFACE_REF);
    }
    const int nseg_all = nseg_in + n_iface;

    MMG5_pMesh mmgMesh = NULL;
    MMG5_pSol  mmgSol  = NULL;

#ifdef THREED
    MMG3D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh,
                    MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end);

    MMG_OK(MMG3D_Set_meshSize(mmgMesh, in.nnode, in.nelem, 0, nseg_all, 0, 0));
    MMG_OK(MMG3D_Set_vertices(mmgMesh, const_cast<double*>(in.coord), NULL));
    MMG_OK(MMG3D_Set_tetrahedra(mmgMesh, conn1.data(), NULL));
    MMG_OK(MMG3D_Set_triangles(mmgMesh, seg1.data(), segref.data()));
    for (int k = 0; k < n_iface; ++k)   // 1-based; interface triangles are the last n_iface entries
        MMG_OK(MMG3D_Set_requiredTriangle(mmgMesh, nseg_in + k + 1));

    if (in.metric_aniso) {
        // metric_aniso is nnode*6 row-major (m11,m12,m13,m22,m23,m33) -- exactly the layout
        // MMG3D_Set_tensorSols expects, so set it in one call (mirrors the scalar branch).
        MMG_OK(MMG3D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, in.nnode, MMG5_Tensor));
        MMG_OK(MMG3D_Set_tensorSols(mmgSol, const_cast<double*>(in.metric_aniso)));
    } else {
        MMG_OK(MMG3D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, in.nnode, MMG5_Scalar));
        MMG_OK(MMG3D_Set_scalarSols(mmgSol, const_cast<double*>(in.metric)));
    }

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
        std::fprintf(stderr, "    [mmg] ERROR: bad ending of MMG3DLIB (strong failure, mesh unusable).\n");
        die(EXIT_MESH_MMG);
    } else if (ier == MMG5_LOWFAILURE) {
        std::fprintf(stderr, "    [mmg] ERROR: bad ending of MMG3DLIB (low failure, mesh saved but imperfect).\n");
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
    if (n_iface) drop_facets_with_ref(out, IFACE_REF);   // interface facets are internal constraints
    // MMG returns each required-element facet as an extra ref=0 copy of the same node set;
    // merge duplicates so var.segment stays duplicate-free (see dedup_facets), then drop the
    // surviving ref=0 interior echoes so they cannot ratchet across remeshes (see
    // drop_unflagged_facets).
    dedup_facets(out.nseg, out.seg, out.segflag);
    if (mesh.is_discarding_internal_segments)
        drop_unflagged_facets(out.nseg, out.seg, out.segflag);

    MMG3D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh,
                   MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end);

#else // 2D

    MMG2D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh,
                    MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end);


    MMG_OK(MMG2D_Set_meshSize(mmgMesh, in.nnode, in.nelem, 0, nseg_all));
    MMG_OK(MMG2D_Set_vertices(mmgMesh, const_cast<double*>(in.coord), NULL));
    MMG_OK(MMG2D_Set_triangles(mmgMesh, conn1.data(), NULL));
    MMG_OK(MMG2D_Set_edges(mmgMesh, seg1.data(), segref.data()));
    for (int k = 0; k < n_iface; ++k)   // 1-based; interface edges are the last n_iface entries
        MMG_OK(MMG2D_Set_requiredEdge(mmgMesh, nseg_in + k + 1));

    if (in.metric_aniso) {
        // metric_aniso is nnode*3 row-major (m11,m12,m22) -- exactly the layout
        // MMG2D_Set_tensorSols expects, so set it in one call (mirrors the scalar branch).
        MMG_OK(MMG2D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, in.nnode, MMG5_Tensor));
        MMG_OK(MMG2D_Set_tensorSols(mmgSol, const_cast<double*>(in.metric_aniso)));
    } else {
        MMG_OK(MMG2D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, in.nnode, MMG5_Scalar));
        MMG_OK(MMG2D_Set_scalarSols(mmgSol, const_cast<double*>(in.metric)));
    }

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
        std::fprintf(stderr, "    [mmg] ERROR: bad ending of MMG2DLIB (strong failure, mesh unusable).\n");
        die(EXIT_MESH_MMG);
    } else if (ier == MMG5_LOWFAILURE) {
        std::fprintf(stderr, "    [mmg] ERROR: bad ending of MMG2DLIB (low failure, mesh saved but imperfect).\n");
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
    if (n_iface) drop_facets_with_ref(out, IFACE_REF);   // interface facets are internal constraints
    // MMG returns each required-element facet as an extra ref=0 copy of the same node set;
    // merge duplicates so var.segment stays duplicate-free (see dedup_facets), then drop the
    // surviving ref=0 interior echoes so they cannot ratchet across remeshes (see
    // drop_unflagged_facets).
    dedup_facets(out.nseg, out.seg, out.segflag);
    if (mesh.is_discarding_internal_segments)
        drop_unflagged_facets(out.nseg, out.seg, out.segflag);

    MMG2D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh,
                   MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end);

#endif // THREED
}

#undef MMG_OK

#endif // USEMMG
