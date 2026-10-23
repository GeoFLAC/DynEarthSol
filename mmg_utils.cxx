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

// Drop ref=0 facets. After dedup_facets every boundary facet carries a BOUND* ref, so a ref=0 facet
// is interior: on output, a required element's echo or, with elem_ref (2D), an edge between two
// references; on input, an internal segment of the initial mesh (refined-zone box, .poly facet).
// MMG needs neither as a constraint (it drops a ref=0 facet between equal references and derives
// one between two references from elem_ref), so this never changes what MMG builds. It keeps
// var.segment boundary-only, which collapse_short_boundary_segments_2d, check_mesh_tangle and
// mmg_adapt's 2D req_edges append rely on (hence validate_parameters rejects
// is_discarding_internal_segments = no under MMG).
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

// MMG's own size at a node on a required edge (nosizreq = 0, its default): the mean length of the
// required edges through it (the required elements' edges and req_edges, each once), as an
// isotropic tensor for metric_aniso. MMG applies it at all such nodes or none, so when a
// keep_metric node lies on a required edge this applies it to met (msize values per node) at the
// other such nodes and returns the number of kept nodes (mmg_adapt then sets nosizreq = 1).
// Returns 0 with met untouched otherwise.
int apply_required_edge_size(const MMGInput &in, int msize, std::vector<double> &met)
{
    if (!in.keep_metric) return 0;
    std::vector<std::array<int, 2>> edges;
    for (int e = 0; in.required_elem && e < in.nelem; ++e) {
        if (!in.required_elem[e]) continue;
        const int *c = in.conn + (std::size_t)e * NODES_PER_ELEM;
        for (int i = 0; i < NODES_PER_ELEM; ++i)
            for (int j = i + 1; j < NODES_PER_ELEM; ++j)
                edges.push_back({{std::min(c[i], c[j]), std::max(c[i], c[j])}});
    }
    for (int k = 0; in.req_edges && k < in.n_req_edges; ++k) {
        const int a = in.req_edges[2*k], b = in.req_edges[2*k + 1];
        edges.push_back({{std::min(a, b), std::max(a, b)}});
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

    std::vector<double> len_sum(in.nnode, 0.0);
    std::vector<int> n_edge(in.nnode, 0);
    for (const auto &ed : edges) {
        double l2 = 0.0;
        for (int d = 0; d < NDIMS; ++d) {
            const double dx = in.coord[(std::size_t)ed[1]*NDIMS + d]
                              - in.coord[(std::size_t)ed[0]*NDIMS + d];
            l2 += dx * dx;
        }
        const double len = std::sqrt(l2);
        for (int k = 0; k < 2; ++k) { len_sum[ed[k]] += len; ++n_edge[ed[k]]; }
    }

    int n_kept = 0;
    for (int n = 0; n < in.nnode; ++n)
        if (n_edge[n] && in.keep_metric[n]) ++n_kept;
    if (n_kept == 0) return 0;
    for (int n = 0; n < in.nnode; ++n) {
        if (!n_edge[n] || in.keep_metric[n]) continue;
        double *m = &met[(std::size_t)n * msize];
        if (msize == 1) { m[0] = len_sum[n] / n_edge[n]; continue; }
        const double lm = n_edge[n] / len_sum[n];   // MMG's form: 1/h, squared on the diagonal
        for (int i = 0, k = 0; i < NDIMS; ++i)
            for (int j = i; j < NDIMS; ++j, ++k) m[k] = (i == j) ? lm * lm : 0.0;
    }
    return n_kept;
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

    // The contract: a required element's nodes are required.
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

    // The metric handed to MMG: the caller's, adjusted by apply_required_edge_size.
    const int msize = in.metric_aniso ? NDIMS * (NDIMS + 1) / 2 : 1;
    const double *msrc = in.metric_aniso ? in.metric_aniso : in.metric;
    std::vector<double> met(msrc, msrc + (std::size_t)in.nnode * msize);
    const int n_kept = apply_required_edge_size(in, msize, met);
    if (n_kept)
        std::printf("      [mmg] kept the metric at %d node(s) on a required edge (keep_metric); "
                    "MMG's required-edge size at the others.\n", n_kept);

    MMG5_pMesh mmgMesh = NULL;
    MMG5_pSol  mmgSol  = NULL;

#ifdef THREED
    // 2D only: MMG3D drops an interior triangle between tetrahedra of equal reference, so a
    // required interface triangle would be lost; references and corners are not wired in 3D.
    if (in.req_edges || in.elem_ref || in.corner_node)
        die(EXIT_INTERNAL_ASSERT, "mmg_adapt: interface edges, element references and corners "
                                  "are 2D only");
    MMG3D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh,
                    MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end);

    MMG_OK(MMG3D_Set_meshSize(mmgMesh, in.nnode, in.nelem, 0, nseg_in, 0, 0));
    MMG_OK(MMG3D_Set_vertices(mmgMesh, const_cast<double*>(in.coord), NULL));
    MMG_OK(MMG3D_Set_tetrahedra(mmgMesh, conn1.data(), NULL));
    MMG_OK(MMG3D_Set_triangles(mmgMesh, seg1.data(), segref.data()));

    if (in.metric_aniso) {
        // metric_aniso's layout is MMG3D_Set_tensorSols' own: one call
        MMG_OK(MMG3D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, in.nnode, MMG5_Tensor));
        MMG_OK(MMG3D_Set_tensorSols(mmgSol, met.data()));
    } else {
        MMG_OK(MMG3D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, in.nnode, MMG5_Scalar));
        MMG_OK(MMG3D_Set_scalarSols(mmgSol, met.data()));
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
    // nosizreq 1 only with a kept node (apply_required_edge_size). MMG never splits or swaps a
    // required edge either way.
    MMG_OK(MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_nosizreq, n_kept > 0 ? 1 : 0));
    MMG_OK(MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmax,  hmax));
    MMG_OK(MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmin,  hmin));
    MMG_OK(MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hausd, hausd));

    const int ier = MMG3D_mmg3dlib(mmgMesh, mmgSol);
    if (ier == MMG5_STRONGFAILURE) {
        std::fprintf(stderr, "    [mmg] ERROR: bad ending of MMG3DLIB (strong failure, mesh unusable).\n");
        die(EXIT_MESH_MMG);
    } else if (ier == MMG5_LOWFAILURE) {
        std::fprintf(stderr, "    [mmg] ERROR: bad ending of MMG3DLIB (low failure, mesh saved but imperfect).\n");
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
    std::vector<int> eref;   // mutable copy: MMG's setters take non-const int*
    if (in.elem_ref) eref.assign(in.elem_ref, in.elem_ref + in.nelem);
    MMG_OK(MMG2D_Set_triangles(mmgMesh, conn1.data(), in.elem_ref ? eref.data() : NULL));
    MMG_OK(MMG2D_Set_edges(mmgMesh, seg1.data(), segref.data()));
    for (int k = 0; k < n_iface; ++k)   // 1-based; interface edges are the last n_iface entries
        MMG_OK(MMG2D_Set_requiredEdge(mmgMesh, nseg_in + k + 1));

    if (in.metric_aniso) {
        // metric_aniso's layout is MMG2D_Set_tensorSols' own: one call
        MMG_OK(MMG2D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, in.nnode, MMG5_Tensor));
        MMG_OK(MMG2D_Set_tensorSols(mmgSol, met.data()));
    } else {
        MMG_OK(MMG2D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, in.nnode, MMG5_Scalar));
        MMG_OK(MMG2D_Set_scalarSols(mmgSol, met.data()));
    }

    MMG_OK(MMG2D_Chk_meshData(mmgMesh, mmgSol));

    if (in.required_node)
        for (int n = 0; n < in.nnode; ++n)
            if (in.required_node[n]) MMG_OK(MMG2D_Set_requiredVertex(mmgMesh, n + 1));
    if (in.required_elem)
        for (int e = 0; e < in.nelem; ++e)
            if (in.required_elem[e]) MMG_OK(MMG2D_Set_requiredTriangle(mmgMesh, e + 1));
    if (in.corner_node)
        for (int n = 0; n < in.nnode; ++n)
            if (in.corner_node[n]) MMG_OK(MMG2D_Set_corner(mmgMesh, n + 1));

    MMG_OK(MMG2D_Set_iparameter(mmgMesh, mmgSol, MMG2D_IPARAM_optim,   0));
    MMG_OK(MMG2D_Set_iparameter(mmgMesh, mmgSol, MMG2D_IPARAM_verbose, mesh.mmg_verbose));
    MMG_OK(MMG2D_Set_iparameter(mmgMesh, mmgSol, MMG2D_IPARAM_debug,   mesh.mmg_debug));
    // nosizreq: see the 3D branch.
    MMG_OK(MMG2D_Set_iparameter(mmgMesh, mmgSol, MMG2D_IPARAM_nosizreq, n_kept > 0 ? 1 : 0));
    MMG_OK(MMG2D_Set_dparameter(mmgMesh, mmgSol, MMG2D_DPARAM_hmax,  hmax));
    MMG_OK(MMG2D_Set_dparameter(mmgMesh, mmgSol, MMG2D_DPARAM_hmin,  hmin));
    MMG_OK(MMG2D_Set_dparameter(mmgMesh, mmgSol, MMG2D_DPARAM_hausd, hausd));

    const int ier = MMG2D_mmg2dlib(mmgMesh, mmgSol);
    if (ier == MMG5_STRONGFAILURE) {
        std::fprintf(stderr, "    [mmg] ERROR: bad ending of MMG2DLIB (strong failure, mesh unusable).\n");
        die(EXIT_MESH_MMG);
    } else if (ier == MMG5_LOWFAILURE) {
        std::fprintf(stderr, "    [mmg] ERROR: bad ending of MMG2DLIB (low failure, mesh saved but imperfect).\n");
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
    out.elem_ref.assign(in.elem_ref ? out.nelem : 0, 0);
    for (int i = 0; i < out.nelem; ++i) {
        int *c = &out.conn[(std::size_t)i*NODES_PER_ELEM];
        MMG_OK(MMG2D_Get_triangle(mmgMesh, &c[0], &c[1], &c[2],
                                  in.elem_ref ? &out.elem_ref[i] : NULL, NULL));
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
    // MMG's ref=0 echoes and reference edges (dedup_facets, drop_unflagged_facets)
    dedup_facets(out.nseg, out.seg, out.segflag);
    drop_unflagged_facets(out.nseg, out.seg, out.segflag);

    MMG2D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh,
                   MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end);

#endif // THREED
}

#undef MMG_OK

#endif // USEMMG
