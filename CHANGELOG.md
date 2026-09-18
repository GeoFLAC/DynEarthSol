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

Nothing yet.

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
