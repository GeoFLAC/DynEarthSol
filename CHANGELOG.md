# Changelog

All notable changes to DynEarthSol are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Entries are grouped by theme rather than by pull request; each references the
pull request that carries the detail. An entry is one or two sentences on what
a user of the code notices -- a parameter, a result that moved, a build that
now works -- not on how; mechanism and evidence stay in the pull request.
Changes with no effect on users of the code -- CI workflows, repository
metadata, release plumbing -- are not listed.
For a release, the full auto-generated list of merged pull requests is in its
[GitHub release notes](https://github.com/GeoFLAC/DynEarthSol/releases).

## [Unreleased]

### Added

- A Docker image with GoSPL coupling: `GOSPL=1 ./build.sh` ships the conda
  environment, gospl_extensions and a 3D executable ready to run ([#96]).
- Contributor documentation: `CONTRIBUTING.md` (renamed from `DEVELOPING.md`)
  with commit-message guidelines and an AI-usage policy, pull request and issue
  templates, and `AGENTS.md` for coding assistants ([#96]).
- Frames embed their step, counts and time, so a run restarts and is read
  without its `.info` index; `utils/recreate_info.py` rebuilds the index ([#96]).
- Provenance records: the executable embeds its build identity, each run writes
  `<modelname>.manifest` (build, host, device, threads, and how the run ended),
  and each frame and checkpoint names the build and machine that wrote it.
  `doc/provenance.md` describes the fields ([#97]).
- `2vtk.py -g` converts a restarted model's whole series, found through the
  runs' `.manifest`, and links the earlier frames under the model's name; the
  new `-hdf` writes updated `.vtkhdf` copies, and `-u` and `-hdf` list the
  frames in a `.vtkhdf.series` for ParaView ([#97]).
- Conservative MMG remeshing (`usemmg=1`) in 2D and 3D: a remesh keeps most
  quiet elements (not yielding, distorted or tiny) as they are, so it no longer
  diffuses the far field, and a yielding band keeps its plastic strain ([#95]).
- With `remeshing_option = 13`, material entering through a restored side wall
  (all four walls in 3D) takes the wall's temperature, layering and radiogenic
  source, and with MMG the mesh along those walls is kept at about `resolution`
  ([#95]).

### Changed

- A GoSPL-coupled run prints its coupling as one summary line after each DES
  progress line, after the first event, instead of four lines per event,
  GoSPL's `--- step done` among them ([#105]).
- `README.md` matches the current build and tools again: several statements
  were stale, among them exported `BOOST_ROOT_DIR` being ignored, a default
  input file, and MMG built from a hand clone rather than the submodule ([#96]).
- The code is called DynEarthSol and DES in the docs and comments, as the
  repository, `CITATION.cff` and the v2.0 paper already do; the 3D suffix
  predates the 2D executable ([#96]).
- `examples/core-complex.cfg` output is time-driven: 34 frames for the 50 kyr
  run instead of 12047 ([#96]).
- **Breaking.** `2vtk.py -u -c` updates the frames in place, as `-u` does;
  `-hdf -c` writes the updated copies into the current directory instead
  ([#97]).
- **Breaking.** `make usegospl=1` looks for the gospl environment in
  `~/miniforge3/envs/gospl`, where GoSPL's install now puts it; a Miniconda
  install sets `CONDA_ENV_PATH`, which, like `GOSPL_EXT_DIR`, may now also be
  exported ([#104]).
- **Breaking.** The make option `use_gospl` is now `usegospl`, matching
  `usemmg` and `useexo`; `make use_gospl=1` stops with an unknown-variable
  error.
- A GoSPL step that fails stops the run with exit code 53 instead of skipping
  that coupling, and GoSPL's diffusion solvers no longer fail on a round-off
  overshoot of the coupling interval ("bad hmax", PETSc error 77).
- The GoSPL docs follow the des3d coupling tutorial and explain when GoSPL
  writes output; the example's `bc` comments now say north and south are
  fixed outlets, as GoSPL reads them ([#104]).
- **Breaking.** `mmg_hmin_factor` and `mmg_hmax_factor` are retired: MMG takes
  its size band from `smallest_size` and `largest_size`, and
  `examples/defaults.cfg` gives the conversion. With `use_mmg_init`, a `.poly`
  region without a size constraint is no longer refined to `resolution` ([#95]).
- **Breaking.** `meshing_option = 2` meshes are finer; multiplying
  `largest_size` by 3.46 (2D) or 5.93 (3D) keeps the old initial mesh ([#95]).
- **Breaking.** Checkpoints written before this change do not restart ([#95]).
- **Breaking.** More cfgs stop at startup: `remeshing_option = 13` on an
  unstructured mesh without MMG (`usemmg=1`), `remeshing_option = 2` or
  `is_discarding_internal_segments = no` with MMG, a velocity period table
  without exactly `num_vbc_period_*` points in non-decreasing time, and a
  one-value `ic.radiogenic_heat_boundry` with radiogenic layers. A remesh now
  stops the run, instead of continuing, when it leaves a broken mesh or when
  the top surface has necked down to near the restored bottom ([#95]).
- **Breaking.** Results change for remeshing on every mesher, for 3D side walls,
  which now follow their velocity period tables as in 2D, and for a one-value
  `ic.radiogenic_heat_mat_in_layer`, now given to every layer. With MMG, a
  plastically strained zone is meshed less finely and keeps the finest size it
  has reached unless `mmg_metric_refine_recover = yes`, and the initial mesh is
  refined in its weak zone ([#95]).

### Fixed

- The thermal mass was built before the initial temperature was set, so thermal
  diffusion ran 5.4 % slow on hot silicate ([#96]).
- A NaN velocity stops the run with exit code 50 instead of writing frames until
  the mesh dies ([#96]).
- Ctrl-C stops a GoSPL-coupled run; Python's SIGINT handler had swallowed it
  ([#96]).
- `openmp=0` links on glibc ([#96]).
- The MMG build ignores an exported `CFLAGS`, which could change the remeshed
  mesh ([#96]).
- HDF5 file locking is off, so ParaView can open a frame while a run writes it.
  Two runs sharing a `modelname` are no longer refused ([#96]).
- Every `mat.mattype_*` is range-checked against `num_materials` at startup;
  an out-of-range index used to read past the per-material arrays ([#96]).
- Every restart cuts `.info` to the rows before its restart frame, keeping the
  old file as `.info.old`; a repeated restart appended duplicate rows. The
  parent frame a restart overwrites is backed up however the two model names
  are spelled ([#97]).
- MMG remeshing removes material that drifted past a restored boundary instead
  of stopping with "BAD ENDING", and a 3D Triangle/TetGen remesh no longer
  moves a steep free surface by hundreds of metres ([#95]).
- With `remeshing_option = 13`, a remesh no longer skews the displacement output
  beside the restored walls ([#95]).
- With `remeshing_option = 12`, a Triangle remesh deletes only the nodes cut off
  under the slanted x0 wall; it used to delete interior nodes too, coarsening
  the mesh at every remesh ([#113]).
- `vbc_val_division_x1_min/max` are fractions of the x1 wall's height, as
  documented; they were taken on the x0 wall ([#113]).
- `--help` gives `bottom_shear_zone_thickness` in metres, the unit the code has
  always used; it said km ([#113]).

## [2.0.2] - 2026-09-17

DOI: [10.5281/zenodo.20330979](https://doi.org/10.5281/zenodo.20330979).
A patch release on top of 2.0.0; the full list of merged pull requests is in the
[release notes](https://github.com/GeoFLAC/DynEarthSol/releases/tag/v2.0.2).
There is no 2.0.1: its tag was spent on a release that had to be withdrawn, and
a published release's tag name cannot be reused.

### Added

- A total-strain-invariant slip rate for rate-and-state friction, selected by
  `rsf_slip_rate_projection_option`, with an optional aging-law timestep bound
  `rsf_dtheta_max`. The strain-rate-based simple-shear benchmark and the paper
  localization configurations come with it ([#91]).
- Homogeneous absolute initial stress at zero gravity, via
  `initial_stress_option` and `initial_stress` ([#91]).
- `mass_scaling_reference_speed`, selecting the K/rho elastic-speed ceiling for
  global velocity scaling instead of the historical G/rho ([#91]).
- SPR stress recovery for remeshing, with a free-surface pin, a surface-relative
  reference and a Deborah-number blend ([#71], [#74]).
- Run-time reporting of the host, device and thread count, and two-digit exit codes
  ([#86]).
- This changelog, which DEVELOPING.md's release workflow had required without the
  file existing.

### Changed

- Boost, libomp and HDF5 are located automatically on Linux and macOS, so a
  fresh clone builds with a plain `make` ([#78], [#81]).
- `mmg` and `nanoflann` are proper git submodules ([#68]).
- Global shape-derivative arrays eliminated, and OpenMP/OpenACC defect fixes,
  for a measured 1.55x and 2.17x ([#70], [#86]).

### Fixed

- Restart is deterministic, and checkpoints carry the lagged state that was
  missing ([#72]).
- Remeshing robustness, restart metadata, and GPU and build fixes ([#76]).
- Rate-and-state aging-law timestep bounds are recomputed from the current
  velocity, and a 1e-16 rate floor is applied to state evolution ([#91]).
- Plastic strain is zeroed again on remesh gap elements ([#93]).
- GPU build failures under nvc++ 26.5 ([#90]).

### Removed

- `control.mixed_stress_reference_viscosity`, which the SPR stress remap
  replaced. An unknown key is a hard error, so a 2.0.0 configuration file that
  sets it no longer loads ([#71]).

### Documentation

- goSPL coupling limitations clarified ([#92]).

## [2.0.0] - 2026-05-19

DOI: [10.5281/zenodo.20293558](https://doi.org/10.5281/zenodo.20293558).
Forty-two merged pull requests; the full list is in the
[release notes](https://github.com/GeoFLAC/DynEarthSol/releases/tag/v2.0.0).

### Added

- Poroelasticity ([#2], [#3]).
- Rate-and-state friction, with a monitor ([#11], [#51]).
- Rayleigh damping ([#18]).
- GPU acceleration and parallelization across the code, later extended with a
  BVH knn search for irregular grids ([#25], [#49], [#56]).
- MMG mesh optimization and MMG mesh initialization ([#44], [#45], [#57]).
- Regular-mesh support ([#16]).
- Multi-segment weak zones, alongside a unified 3D lateral-face velocity
  boundary condition dispatch ([#65]).
- Two-way DES-goSPL coupling, on the ASPECT-FastScape scheme ([#63]).
- `Array2D` structure-of-arrays layouts ([#55]).
- Functional test configurations, and CI with HDF5 ([#40]).
- A profiling level, and timing and performance statistics ([#36], [#39]).
- 2vtk VTKHDF update mode ([#42]).
- Surface-only interpolation for sediment deposits during remeshing, and 3D
  marker correction for surface processes ([#15], [#29]).

### Changed

- **Breaking.** Output and checkpoint files moved from NetCDF to HDF5, which also
  lets ParaView read model data directly ([#34], [#35]).
- **Breaking.** `info_display_interval` is now `info_display_step_interval`
  ([#46]).
- `MatProps` and `Output` simplified, legacy code dropped, and remeshing memory
  management improved ([#37]).

### Fixed

- Stress in `MatProps::rh_evp`, and a missing `break` for the same branch in
  `update_stress()` ([#12], [#24]).
- CPU and HDF5 correctness in barycentric coordinates, interpolation, OpenMP and
  HDF5 output ([#48]), and an index limit in `interpolate_nd` ([#54]).
- Meshing and NumPy compatibility ([#58]), and `normal_vector_of_facet()`,
  `MarkerSet::read_chkpt_file()` and `find_acm_elem_ratios()` ([#33]).
- macOS OpenMP installation and scripts ([#50]).
- Further rounds of bug fixes, including around rate-and-state friction and the
  3D dome parameters ([#13], [#17], [#66]).

### Removed

- `libadaptivity` ([#30]).

## [1.0.0] - 2020-07-01

First tagged release. It predates the changelog and the GitHub release notes; see
the [`v1.0.0`](https://github.com/GeoFLAC/DynEarthSol/releases/tag/v1.0.0) tag.

[Unreleased]: https://github.com/GeoFLAC/DynEarthSol/compare/v2.0.2...master
[2.0.2]: https://github.com/GeoFLAC/DynEarthSol/compare/v2.0.0...v2.0.2
[2.0.0]: https://github.com/GeoFLAC/DynEarthSol/compare/v1.0.0...v2.0.0
[1.0.0]: https://github.com/GeoFLAC/DynEarthSol/releases/tag/v1.0.0
[#2]: https://github.com/GeoFLAC/DynEarthSol/pull/2
[#3]: https://github.com/GeoFLAC/DynEarthSol/pull/3
[#11]: https://github.com/GeoFLAC/DynEarthSol/pull/11
[#12]: https://github.com/GeoFLAC/DynEarthSol/pull/12
[#13]: https://github.com/GeoFLAC/DynEarthSol/pull/13
[#15]: https://github.com/GeoFLAC/DynEarthSol/pull/15
[#16]: https://github.com/GeoFLAC/DynEarthSol/pull/16
[#17]: https://github.com/GeoFLAC/DynEarthSol/pull/17
[#18]: https://github.com/GeoFLAC/DynEarthSol/pull/18
[#24]: https://github.com/GeoFLAC/DynEarthSol/pull/24
[#25]: https://github.com/GeoFLAC/DynEarthSol/pull/25
[#29]: https://github.com/GeoFLAC/DynEarthSol/pull/29
[#30]: https://github.com/GeoFLAC/DynEarthSol/pull/30
[#33]: https://github.com/GeoFLAC/DynEarthSol/pull/33
[#34]: https://github.com/GeoFLAC/DynEarthSol/pull/34
[#35]: https://github.com/GeoFLAC/DynEarthSol/pull/35
[#36]: https://github.com/GeoFLAC/DynEarthSol/pull/36
[#37]: https://github.com/GeoFLAC/DynEarthSol/pull/37
[#39]: https://github.com/GeoFLAC/DynEarthSol/pull/39
[#40]: https://github.com/GeoFLAC/DynEarthSol/pull/40
[#42]: https://github.com/GeoFLAC/DynEarthSol/pull/42
[#44]: https://github.com/GeoFLAC/DynEarthSol/pull/44
[#45]: https://github.com/GeoFLAC/DynEarthSol/pull/45
[#46]: https://github.com/GeoFLAC/DynEarthSol/pull/46
[#48]: https://github.com/GeoFLAC/DynEarthSol/pull/48
[#49]: https://github.com/GeoFLAC/DynEarthSol/pull/49
[#50]: https://github.com/GeoFLAC/DynEarthSol/pull/50
[#51]: https://github.com/GeoFLAC/DynEarthSol/pull/51
[#54]: https://github.com/GeoFLAC/DynEarthSol/pull/54
[#55]: https://github.com/GeoFLAC/DynEarthSol/pull/55
[#56]: https://github.com/GeoFLAC/DynEarthSol/pull/56
[#57]: https://github.com/GeoFLAC/DynEarthSol/pull/57
[#58]: https://github.com/GeoFLAC/DynEarthSol/pull/58
[#63]: https://github.com/GeoFLAC/DynEarthSol/pull/63
[#65]: https://github.com/GeoFLAC/DynEarthSol/pull/65
[#66]: https://github.com/GeoFLAC/DynEarthSol/pull/66
[#68]: https://github.com/GeoFLAC/DynEarthSol/pull/68
[#70]: https://github.com/GeoFLAC/DynEarthSol/pull/70
[#71]: https://github.com/GeoFLAC/DynEarthSol/pull/71
[#72]: https://github.com/GeoFLAC/DynEarthSol/pull/72
[#74]: https://github.com/GeoFLAC/DynEarthSol/pull/74
[#76]: https://github.com/GeoFLAC/DynEarthSol/pull/76
[#78]: https://github.com/GeoFLAC/DynEarthSol/pull/78
[#81]: https://github.com/GeoFLAC/DynEarthSol/pull/81
[#86]: https://github.com/GeoFLAC/DynEarthSol/pull/86
[#90]: https://github.com/GeoFLAC/DynEarthSol/pull/90
[#91]: https://github.com/GeoFLAC/DynEarthSol/pull/91
[#92]: https://github.com/GeoFLAC/DynEarthSol/pull/92
[#93]: https://github.com/GeoFLAC/DynEarthSol/pull/93
[#95]: https://github.com/GeoFLAC/DynEarthSol/pull/95
[#96]: https://github.com/GeoFLAC/DynEarthSol/pull/96
[#97]: https://github.com/GeoFLAC/DynEarthSol/pull/97
[#104]: https://github.com/GeoFLAC/DynEarthSol/pull/104
[#105]: https://github.com/GeoFLAC/DynEarthSol/pull/105
[#113]: https://github.com/GeoFLAC/DynEarthSol/pull/113
