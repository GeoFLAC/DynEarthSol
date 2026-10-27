#ifndef DYNEARTHSOL3D_MESH_HPP
#define DYNEARTHSOL3D_MESH_HPP

void points_to_new_mesh(const Mesh &mesh, int npoints, const double *points,
                        int n_init_segments, const int *init_segments, const int *init_segflags,
                        int n_regions, const double *regattr,
                        double max_elem_size, int vertex_per_polygon,
                        int &nnode, int &nelem, int &nseg,
                        double *&pcoord, int *&pconnectivity,
                        int *&psegment, int *&psegflag, double *&pregattr);
void points_to_new_surface(const Mesh &mesh, int npoints, const double *points,
                           int n_init_segments, const int *init_segments, const int *init_segflags,
                           int n_regions, const double *regattr,
                           double max_elem_size, int vertex_per_polygon,
                           int &nnode, int &nelem, int &nseg,
                           double *&pcoord, int *&pconnectivity,
                           int *&psegment, int *&psegflag, double *&pregattr);
void renumbering_mesh(const Param& param, array_t &coord, conn_t &connectivity,
                      segment_t &segment, regattr_t *regattr, double_vec *init_elem_size_hint = nullptr);
void create_boundary_flags2(uint_vec &bcflag, int nseg,
                            const segment_t& segment, const segflag_t& segflag);
void create_boundary_flags(Variables& var);
void create_boundary_nodes(Variables& var);
void create_boundary_facets(Variables& var);
void create_top_elems(Variables& var);
void create_surface_info(const Param& param, const Variables& var, SurfaceInfo& surfinfo);
void update_surface_info(const Variables& var, SurfaceInfo& surfinfo);
void create_support(Variables& var);
void create_neighbor(Variables& var);
void create_elemmarkers(const Param&, Variables&);
void create_markers(const Param&, Variables&);
void create_new_mesh(const Param&, Variables&);
void create_equilateral_elem(const Variables& var, int *&connectivity);
void create_equilateral_segments(const Variables& var, int *&segments, int *&segflags);
template <typename T_in, typename T_out>
void average_nodal_to_elem(T_in nodal, const conn_t &connectivity, int nelem, T_out elem, bool is_surface=false);


/* VirtualTopology -- a background Cartesian grid that bins markers (or any
 * point cloud) into per-cell marker-id buckets, for cheap spatial queries
 * that would otherwise require searching the unstructured Lagrangian mesh:
 * any indexed marker can look up who shares its cell, and a caller can step
 * to the axis-adjacent cell (neighbor_cell()). O(1) bucket lookup, not a
 * nearest-point search: a cell holds whatever is binned there, which may be
 * empty even when the true nearest neighbor is only slightly further away
 * (quantized by cell size).
 *
 * Used by serpentinization: a mantle marker converts when a cell within
 * mat.serpentinization_thickness below it is more than some threshold
 * fraction wet crust (counts_as_wet_crust(); see build_serpentinization_grid()
 * in phasechanges.cxx, which computes the fraction from cell_markers() +
 * MarkerSet::get_mattype() at each build, and subduction_phase_changes(),
 * which reads the resulting wet-crust-dominant cells).
 *
 * build_index() runs at init, at each quality check
 * (mesh.quality_check_step_interval) and at the end of each remesh: markers
 * move continuously in between, so absolute positions drift, but their
 * cell-adjacency (who is above/below/left/right) changes far more slowly.
 * Between rebuilds no marker of the indexed set is created or deleted:
 * its count changes only in the surface processes of a quality-check step
 * (rebuilt in the same step) and in a remesh (rebuilt at its end). So every
 * marker id has a bucket and a label. MarkerSet::remove_markers() reuses the
 * freed ids for other markers, so a stale index would read the wrong marker;
 * subduction_phase_changes() therefore dies when the marker count is not
 * num_labels() or the set's last id (MarkerSet::get_last_id(), which every
 * append increments) is not built_last_id().
 *
 * Thread safety: no member takes a lock. The mutators run in serial code
 * only: build_index() and the dominance and extent setters in
 * build_serpentinization_grid(), and invalidate_dominant_cell() in the host step between
 * the passes of subduction_phase_changes(). serialize() runs in
 * Output::write_checkpoint(), and restart() replaces the object through
 * deserialize(). Only label_cell(), neighbor_cell(),
 * is_dominant() and dominant_xext_ok() run from the passes (OpenMP threads,
 * or device kernels in an ACC build), and they never meet a writer: the
 * mutators run on the host between those loops, in program order, and each
 * loop ends (a thread join or an acc wait) before the next mutator runs. A
 * change that mutates the index while such a loop runs needs its own
 * synchronization.
 *
 * column_index()/vertical_index()/cell_of()/neighbor_cell() are pure
 * functions of dims_/lo_/dx_, fixed at construction
 * and never mutated after.
 *
 * Because that lagged index (bucket_ and its labels) drives a live
 * decision at every phase_changes() call, it IS part of the checkpoint
 * (serialize()/deserialize()): a restart must carry the identical state the
 * continuous run had at that step, not a fresher rebuild against current
 * marker positions -- a mismatch can change which markers convert, so the
 * restarted run can diverge.
 *
 * `dominant_` (managed by clear_dominant()/set_dominant()/is_dominant()/
 * invalidate_dominant_cell()) flags the precomputed cells whose indexed-marker
 * composition exceeds some threshold fraction of a designated material (for
 * serpentinization: see build_serpentinization_grid()). The caller computes it
 * once per build_index() call (in build_serpentinization_grid()) rather than
 * recomputing cell composition on every subduction_phase_changes() step:
 * which cells are dominant changes far more slowly than phase_changes() is
 * called: a caller reads the flags instead of running cell_markers() + a
 * fraction count per cell per call.
 *
 * The class is deliberately decoupled from Param/Variables/MarkerSet: the
 * caller passes the points it bins and their bounds. The one caller,
 * build_serpentinization_grid(), takes the marker coordinates
 * (MarkerSet::calculate_marker_coord()), replaces each z by its depth below
 * the local surface (SurfaceTopo::zeff()), and bounds the grid by a min/max
 * scan over those (x, zeff) points.
 *
 * Points outside the domain are clamped to the nearest cell.  The vertical
 * axis is index NDIMS-1, positive up (DES convention); the grid works
 * unchanged in 2D and 3D builds. In 2D, axis 0 is horizontal (left/right);
 * in 3D, axes 0 and 1 are both horizontal and axis NDIMS-1 remains vertical.
 */
class VirtualTopology
{
public:
    // domain_lo/domain_hi: NDIMS-length arrays bounding the domain.
    // resolution: uniform cell size along every axis.
    VirtualTopology(const double* domain_lo, const double* domain_hi,
                    double resolution);

    // Rebuild the per-cell marker-id buckets from scratch, binning every
    // marker index m (into coords, 0..npoints-1), and record last_id, the
    // marker set's last id (see built_last_id()). Replaces whatever
    // build_index() built before.
    void build_index(const array_t& coords, int npoints, int last_id);

    // The marker set's last id at the last build_index(): with num_labels(),
    // it shows whether a marker was added or removed since.
    int built_last_id() const { return built_last_id_; }

    // The flat cell id containing point p, for use with neighbor_cell()/
    // cell_markers(). Pure function.
    int cell_of(const double* p) const { return column_index(p) * nz_ + vertical_index(p); }

    // The cell id one step over from `cell` along `axis` (0..NDIMS-1) in
    // direction `dir` (+1 or -1); -1 if that step leaves the grid. Lets a
    // caller walk multiple steps (e.g. down a column) without re-deriving
    // coordinates each step. Pure function.
    #pragma acc routine seq
    int neighbor_cell(int cell, int axis, int dir) const;

    // Marker ids build_index() binned into `cell`, which must be in
    // [0, num_cells()) (dies otherwise). Empty when nothing is bound there,
    // which includes every cell of a new grid before its first build_index().
    const int_vec& cell_markers(int cell) const;

    // The number of markers the last build_index() binned (or deserialize()
    // restored): label_cell() is valid for ids 0..num_labels()-1.
    int num_labels() const { return static_cast<int>(label_.size()); }

    // Total number of cells in the grid (ncol_ * nz_). Exposed so callers
    // (e.g. build_serpentinization_grid()) can iterate all cells to compute
    // dominance without knowing the internal layout. Pure.
    int num_cells() const { return ncol_ * nz_; }

    // Uniform cell size along every axis (the `resolution` the grid was
    // constructed with), for callers converting a physical distance (e.g.
    // mat.serpentinization_thickness) into a number of cells to walk. Pure.
    double cell_size() const { return dx_; }

    // --- Dominance cache (see the class doc above) ---

    // Clear every dominance flag and wet-crust extent. Called at the start of
    // dominance precomputation in build_serpentinization_grid().
    void clear_dominant();

    // Flag `cell` as dominant. Called from build_serpentinization_grid()
    // after counting wet-crust composition per cell.
    void set_dominant(int cell);

    // O(1) query: true iff `cell` is flagged dominant.
    #pragma acc routine seq
    bool is_dominant(int cell) const { return dominant_[cell] != 0; }

    // Set below[c] = 1 for every cell c that has a dominant cell at most
    // nsteps cells strictly below it in its column, 0 elsewhere (below is
    // resized to num_cells()): the walk down from every cell at once, in
    // num_cells() + (dominant cells) * nsteps steps.
    void cells_with_dominant_below(int nsteps, std::vector<char>& below) const;

    // Unflag `cell` (conservative: treated as non-dominant until the next
    // build, not recounted). Called between the passes of
    // subduction_phase_changes() for each cell where pass 1 turned a wet-crust
    // marker into a material that is not wet crust, so pass 2 does not read a
    // dominance that the conversion may have ended. Its wet-crust extent is
    // then never read again.
    void invalidate_dominant_cell(int cell);

    // --- Per-dominant-cell wet-crust horizontal extent (optional refinement) ---
    //
    // A dominant cell answers "this cell holds enough wet crust" for its WHOLE
    // horizontal span (one cell per axis), so a marker at one edge of a
    // column inherits wet crust sitting up to a full cell size away -- the
    // horizontal-quantization false positive at a slab bend. When the caller
    // additionally records the horizontal min/max of the cell's qualifying
    // (wet-crust) markers at build time, dominant_xext_ok() narrows the answer
    // to "wet crust within `margin` of this horizontal position".
    //
    // Storage is caller-driven exactly like dominance: populated in
    // build_serpentinization_grid() (only when the margin feature is
    // enabled -- absent extents cost nothing in memory or checkpoints) and
    // cleared by clear_dominant(). xlo/xhi are the NDIMS-1 horizontal
    // coordinates (axis order 0..NDIMS-2).
    void set_dominant_xext(int cell, const double* xlo, const double* xhi);

    // Signal that the extent scan is complete / query it: a checkpoint
    // written with the margin feature off must trigger a rebuild when a
    // restart turns it on, not silently pass every margin test.
    void mark_xext_built();
    bool has_xext_data() const { return xext_built_; }

    // True iff every horizontal component of p lies within the dominant
    // `cell`'s wet-crust extent widened by margin, [lo - margin, hi + margin].
    // Only for has_xext_data() grids, where every dominant cell has an extent
    // (serialize() asserts it, deserialize() reads one per dominant id).
    #pragma acc routine seq
    bool dominant_xext_ok(int cell, const double* p, double margin) const;

    // --- Per-marker build-time cell label ---
    //
    // label_cell(m) is the cell marker m was binned into at the last
    // build_index() -- the inverse of the bucket mapping.
    // subduction_phase_changes() always queries by this label so marker-grid
    // relationships are evaluated in the build-time configuration and stay
    // constant between rebuilds (live-position queries let a drifting marker
    // flip the answer mid-interval against the frozen dominance data).
    // m must be in [0, num_labels()). Derived data: reconstructed from
    // bucket_ on deserialize, so it adds NO checkpoint fields.
    #pragma acc routine seq
    int label_cell(int m) const { return label_[m]; }

    // Checkpoint round-trip. Kept decoupled from binaryio (per this class's
    // design): the owner does the file I/O. Both vectors are variable length
    // (buckets in ints, optional wet-crust extents in dbls), so the owner
    // must also stash each one's size (length-1 arrays vtopo_len /
    // vtopo_dbl_len) to pre-size the buffers before the read_array() that
    // fills them -- see output.cxx/dynearthsol.cxx for the pattern.
    // ints = [layout_tag, dims_[0..NDIMS-1],
    //         bucket_[c].size() for c in 0..ncells), bucket_ contents concatenated,
    //         built_last_id_,
    //         xext_built_ (0 or 1), the number of dominant cells, their ids
    //         (serialize() writes them ascending; deserialize() accepts any
    //         order, and a reader must not rely on it)].
    // deserialize() dies (EXIT_IO_RESTART) on any other first int.
    // dbls = [lo_[0..NDIMS-1], dx_,
    //         per dominant cell (same order as its ids in ints): its wet crust
    //         extent lo[0..NDIMS-2], hi[0..NDIMS-2] -- present exactly when
    //         xext_built_ (margin feature on at build; otherwise dbls keeps
    //         the plain NDIMS+1 length). Like ints, dbls is therefore
    //         variable length and the owner stashes its size too
    //         (vtopo_dbl_len).]
    // nz_/ncol_ are re-derived from dims_ on load, so they are not stored.
    void serialize(int_vec& ints, double_vec& dbls) const;
    static VirtualTopology* deserialize(const int_vec& ints, const double_vec& dbls);

private:
    // First int of serialize()'s ints; change it whenever the layout changes.
    // A negative tag never equals a cell count (dims_[0] >= 1), so ints that
    // start without a tag never match; nor do ints with another layout's tag.
    enum { layout_tag = -4 };

    VirtualTopology() : dx_(0), nz_(0), ncol_(0)
    { for (int d=0; d<NDIMS; d++) { lo_[d]=0; dims_[d]=0; } }   // for deserialize()

    double lo_[NDIMS];
    double dx_;
    int dims_[NDIMS];   // # of cells per axis
    int nz_;            // dims_[NDIMS-1], vertical
    int ncol_;          // product of the NDIMS-1 horizontal dims_

    // Per-cell marker-id buckets, flat cell numbering (see cell_of()): one
    // per cell, all empty from construction until build_index() fills them;
    // deserialize() restores them as checkpointed.
    int_vec2D bucket_;
    // One flag per cell: 1 when its indexed-marker population exceeds the
    // caller's dominance threshold (build_serpentinization_grid() for
    // serpentinization). Set once per build_index() by the caller via
    // clear_dominant()/set_dominant(); partially cleared between rebuilds via
    // invalidate_dominant_cell() when a wet-crust marker turns into a material
    // that is not wet crust. See the class doc above.
    std::vector<char> dominant_;

    // Wet-crust horizontal extent per cell (see set_dominant_xext()), xext_stride
    // doubles per cell: lo[0..NDIMS-2], then hi[0..NDIMS-2]. Sized, NaN for a
    // cell without one, at the first set_dominant_xext() or mark_xext_built(),
    // or by deserialize() when xext_built_; empty otherwise.
    enum { xext_stride = 2 * (NDIMS-1) };
    double_vec xext_;
    bool xext_built_ = false;       // set by mark_xext_built(); checkpointed as is

    // Build-time cell of every marker id the last build binned (see
    // label_cell()). Rebuilt by build_index(), reconstructed from bucket_ by
    // deserialize() -- never checkpointed itself.
    int_vec label_;
    int built_last_id_ = 0;   // see built_last_id(); checkpointed

    int column_index(const double* p) const;   // horizontal cell -> [0, ncol_)
    int vertical_index(const double* p) const; // vertical cell -> [0, nz_)
};

#endif
