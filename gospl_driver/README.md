# GoSPL Coupling for DynEarthSol

This directory contains the code that couples DynEarthSol (DES) with
[GoSPL](https://github.com/Geodels/gospl) (Global Scalable Paleo Landscape
Evolution), a Python landscape evolution model for river incision, sediment
transport, hillslope diffusion and marine deposition. DES handles tectonics;
GoSPL handles the surface.

For a step-by-step walkthrough with exercises, see the tutorial
[Coupling with GoSPL](https://geoflac.github.io/des3d/docs/tutorial/couplinggospl).
This file is the reference kept with the code.

## Files

- `gospl-driver.hpp` - Header file for the GoSPL C++ interface wrapper
- `gospl-driver.cxx` - Implementation of the GoSPL C++ interface wrapper
- `examples/` - A ready-to-run Gaussian weak zone rift, a DES cfg and its GoSPL YAML

## Prerequisites

To skip the setup below, use the Docker image ([Run with Docker](#run-with-docker)).

1. **GoSPL in a conda environment with Python 3.11**, by default at
   `~/miniforge3/envs/gospl`. Follow
   [the GoSPL installation procedure](https://gospl.readthedocs.io/en/latest/getting_started/installConda.html).

2. **gospl_extensions**, the bridge between GoSPL (Python) and DES (C++). It
   provides a C++ entry point (`libgospl_extensions.so` and its header), the
   `EnhancedModel` GoSPL subclass that can be stepped externally, take imposed
   velocities and return an elevation change, and the interpolation between the
   DES and GoSPL meshes, which are built independently and at different
   resolutions. Build it inside the gospl environment:
   ```bash
   git clone https://github.com/GeoFLAC/gospl_extensions.git ~/opt/gospl_extensions
   cd ~/opt/gospl_extensions/cpp_interface
   conda activate gospl
   make install-local
   ```
   On success it reports
   `✅ Installed locally to gospl_extensions/lib and gospl_extensions/include`,
   having created `lib/libgospl_extensions.so` and `include/gospl_extensions.h`,
   the layout DynEarthSol's Makefile expects.

## Build

GoSPL coupling is 3D only. Build **outside** the gospl environment, which keeps
the compiler off conda's libraries:

```bash
conda deactivate   # if any environment is active
make clean
make ndims=3 use_gospl=1 usemmg=1 -j4
```

- `ndims=3` is required; `usemmg=1` (MMG mesh optimization during remeshing) is
  recommended but optional.
- `GOSPL_EXT_DIR` (default `~/opt/gospl_extensions`) and `CONDA_ENV_PATH`
  (default `~/miniforge3/envs/gospl`) locate the two prerequisites; set them on
  the command line, in the Makefile or as exported variables if yours are
  elsewhere.
  `PYTHON_VERSION`, `PYTHON_INCLUDE_DIR` and `PYTHON_LIB_DIR` cover a
  non-conda Python.
- The build also writes `dynearthsol-gospl`, a wrapper that puts
  `gospl_extensions/cpp_interface` on `PYTHONPATH` and runs `dynearthsol3d`.
  A successful build ends with `✅ DynEarthSol built with GoSPL support!`.

### Verify the build

```bash
./dynearthsol3d --help | grep gospl                # options registered?
conda activate gospl && python -c "import gospl"   # GoSPL importable?
ldd dynearthsol3d | grep python                    # linked to Python?
cat dynearthsol-gospl                              # wrapper script written?
```

If `dynearthsol-gospl` is missing, the build did not complete with
`use_gospl=1`.

## Configuration

### DES parameters

All are in the `[control]` section; `examples/defaults.cfg` lists them and
`./dynearthsol3d --help` describes them.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `surface_process_option` | 0 | Set to **11** to enable GoSPL |
| `surface_process_gospl_config_file` | (empty) | Path to the GoSPL YAML, resolved relative to the working directory |
| `gospl_coupling_mode` | `steps` | `steps` or `time`, selecting which of the next two triggers coupling |
| `gospl_coupling_frequency` | 1 | Couple every N DES steps (`steps` mode) |
| `gospl_coupling_interval_in_yr` | 1000 | Couple every T model years (`time` mode) |
| `gospl_velocity_coupling` | `true` | Pass surface velocities (vx, vy, vz) to GoSPL |
| `gospl_mesh_resolution` | -1 | GoSPL node spacing in metres; -1 sizes it from the DES surface nodes |
| `gospl_mesh_padding` | 0.1 | Fraction by which the GoSPL mesh extends beyond the DES domain on each side |
| `gospl_mesh_perturbation` | 0.3 | Random perturbation of GoSPL node positions, as a fraction of the spacing (0-1) |

```cfg
[control]
surface_process_option = 11
surface_process_gospl_config_file = gospl_config.yml
gospl_coupling_mode = steps
gospl_coupling_frequency = 100
gospl_mesh_resolution = 500
```

For slow erosion, a large `gospl_coupling_frequency` (100 or more) cuts run
time: GoSPL runs less often over the accumulated time. If the result changes
when you lower it, the interval was too coarse.

### GoSPL YAML

The coupling drives the `EnhancedModel` from gospl_extensions, not stock GoSPL,
so a standard GoSPL configuration may not work. Start from
`examples/gospl_config_gaussian_weakzone_3D.yml` and keep its five required
sections: `domain`, `time`, `spl`, `diffusion` and `output`. The keys most
often changed:

| Key | What it sets |
|-----|--------------|
| `spl: K` | Bedrock river incision rate (erodibility) |
| `spl: m`, `spl: n` | Drainage-area and slope exponents of the stream power law (`n` defaults to 1) |
| `diffusion: hillslopeKa` | Hillslope diffusivity, m²/yr |
| `domain: flowdir` | Flow routing (`6` is multi-direction) |
| `domain: bc` | Boundaries in the order N, E, S, W: `o` open, `f` fixed, `w` wall. `'wowo'` opens east and west and closes north and south |
| `domain: seadepo` | Marine deposition on or off |
| `sea: position` | Sea level in metres relative to the initial surface |

`time: dt` and `time: end` are overridden by DES at run time, but they still
bound `tout` (see [Output timing](#output-timing)). The GoSPL mesh is
generated at startup and saved as `gospl_mesh.npz` in the working directory.

### Output timing

GoSPL writes output on its own clock, every `time: tout` years, but only at a
coupling event. Each event runs GoSPL for one step, and GoSPL checks for output
at that step's start and end. Its clock starts at `time: start` and advances by
the coupling interval, so it trails DES time by whatever has accumulated since
the last event.

- **Outputs snap to coupling events.** A file is written at the first coupling
  event at or after each multiple of `tout`. If the coupling interval divides
  `tout`, as `gospl_coupling_mode = time` makes easy, outputs are exactly `tout`
  apart; otherwise the spacing is uneven.
- **A coupling interval longer than `tout` mislabels outputs.** Every event then
  writes, but the time stamped in the XDMF is the nominal output time,
  `start + k·tout`, which falls further behind the model time with each write.
  In `steps` mode the interval in years varies with DES's adaptive `dt`, so this
  can happen unnoticed.
- **`dt` and `end` still bound `tout`.** When GoSPL reads the YAML it raises
  `tout` to `dt` if smaller, and lowers it to `end − start` if
  `start + tout > end`, printing a one-line notice for each. Keep `dt ≤ tout`
  and `end` at least `max_time_in_yr`.

To line GoSPL frames up with DES frames, set `tout` to
`output_time_interval_in_yr` and couple at a divisor of it. GoSPL frames still
fall on coupling events, up to one coupling interval from the DES frame.

## Run

Run from the directory that holds the GoSPL YAML, since the path in
`surface_process_gospl_config_file` is resolved relative to the working
directory, not to the cfg. Or give an absolute path.

```bash
conda activate gospl
cd DynEarthSol/gospl_driver/examples
../../dynearthsol-gospl ./gaussian-weakzone-3d-with-gospl.cfg
```

To manage `PYTHONPATH` yourself, run the executable directly:

```bash
conda activate gospl
export PYTHONPATH="$HOME/opt/gospl_extensions/cpp_interface:${PYTHONPATH}"
cd DynEarthSol/gospl_driver/examples
../../dynearthsol3d ./gaussian-weakzone-3d-with-gospl.cfg
```

DES output goes to the working directory, GoSPL output to the YAML's
`output: dir`.

### Run with Docker

`GOSPL=1 ./build.sh` in the repository root builds `dynearthsol/gcc-11-gospl`
with the gospl conda environment, gospl_extensions and a 3D `dynearthsol3d`
already inside, and the environment activated in every login shell. Mount the
directory holding your cfg and GoSPL YAML and run from it:

```bash
docker run --rm -it -v /path/to/case:/home/human/case dynearthsol/gcc-11-gospl \
  bash -lc 'cd ~/case && ~/DynEarthSol/dynearthsol-gospl your_input.cfg'
```

## Coupling Details

DES and GoSPL exchange data following the ASPECT-FastScape simple coupling
scheme.

1. **Initialization**: GoSPL is initialized once from the YAML config file. At the first coupling event, GoSPL's elevation field (`hGlobal`) is seeded from DES's initial surface via `apply_elevation_data()`.
2. **Each coupling event** (every `gospl_coupling_frequency` DES steps, or every `gospl_coupling_interval_in_yr` years when `gospl_coupling_mode = time`):
   - **Time-averaged tectonic velocity** is computed as `Δcoord/Δt` over the coupling interval (not the instantaneous DES velocity), where `Δcoord` is each surface node's displacement since the previous event and `Δt` the model time elapsed since then. DES uses inertial scaling (quasi-dynamic formulation), so instantaneous velocities contain damped-wave components that are numerical artifacts and would perturb GoSPL's drainage network. Averaging over the interval filters these out. On the first coupling event, instantaneous velocity is used as a fallback.
   - DES surface velocities (vx, vy, vz) are IDW-interpolated onto the GoSPL mesh via `set_surface_velocity()`.
   - `run_and_get_erosion(dt)` advances GoSPL by one step of length `dt`. Internally: horizontal advection (vx, vy), vertical uplift (vz via `upsub`), SPL erosion, and hillslope diffusion are applied to GoSPL's `hGlobal`. The returned `delta_h` contains **only the erosion and diffusion component** — the uplift (`upsub * dt`) is subtracted before returning because DES already applied the same displacement through its Lagrangian mechanical solver. Returning the full delta_h would double-count the tectonic uplift.
   - DES adds `delta_h` (erosion + diffusion only) to its surface node z-coordinates. `delta_h` is not `Δcoord_z`: the latter is DES's tectonic displacement, the former the surface-process change added on top of it.
3. GoSPL **owns** the topography between coupling events and accumulates its drainage network state continuously. DES remeshing does **not** reset GoSPL's elevation or drainage state, so drainage divides survive mesh adaptation.
4. **Padded GoSPL mesh**: the GoSPL mesh extends beyond the DES domain by `gospl_mesh_padding` on each side, keeping DES surface nodes away from GoSPL's boundary artifacts.

```
 DES time    t_prev   ────────── DES steps ──────────    t_now
               │◄────────── Δt = t_now − t_prev ──────────►│
          coord_prev                                   coord_now
               └───── Δcoord = coord_now − coord_prev ─────┘
 At t_now:
   DES ───────── v̄ = Δcoord / Δt  (vx, vy, vz) ─────────►  GoSPL
                                                             │  advances Δt:
                                                             │  advection, uplift,
                                                             │  incision, diffusion
   DES ◄──── Δh (erosion + diffusion, uplift removed) ───────┘
   DES then sets   z_surface ← z_surface + Δh
```

### Known Limitations

- **Coupling interval is a trigger, not a clamp.** In `gospl_coupling_mode = time`, `gospl_coupling_interval_in_yr` only gates *when* coupling fires (`accumulated_dt >= gospl_coupling_interval_in_yr`); the `dt` actually passed to GoSPL is `accumulated_dt` itself. If DES's adaptive `dt` exceeds the configured interval, coupling fires every DES step and GoSPL's Δt silently becomes DES's `dt` instead of the configured value — there is no sub-stepping or truncation back to the nominal interval.

- **Remeshing mid-interval.** The coupling clock (`accumulated_dt`/`step_counter`) is time-based and unaffected by DES remeshing — a remesh occurring between two coupling events does not skip or reset the coupling schedule, and GoSPL's elevation state is not reseeded from DES after a remesh (consistent with "GoSPL owns the topography" above). However, the internal time-averaged-velocity calculation only guards against a *change in surface node count* across the remesh; if a remesh happens to leave the top-boundary node count unchanged (common when only the interior remeshes), node identity/order is not otherwise verified, and the "time-averaged" velocity computed for the next coupling event can silently difference unrelated nodes. Treat the coupling event immediately following a remesh with caution until this is hardened.

- **GoSPL mesh domain is fixed at initialization.** `generate_mesh()` runs once, at startup, sized to the DES model's *initial* top-surface extent plus `gospl_mesh_padding` (default 10%) on each side. It is never regenerated during the run (and on restart, an existing mesh file is reused as-is rather than rebuilt). There is currently no mechanism to re-center or resize the GoSPL domain as the DES model deforms. Consequently, the padding fraction effectively upper-bounds the lateral extension the DES model can accumulate before its surface nodes migrate out of the padded domain's interior and approach the GoSPL mesh boundary, where the truncated-drainage-basin / BC-enforcement artifacts the padding was meant to avoid can reappear. No runtime check warns when this happens, so choose a generous `gospl_mesh_padding` for strongly extensional models.

### Coupling API (`GoSPLDriver` C++ class)

These calls are internal, but their names appear in the log.

| Method | Purpose |
|---|---|
| `set_surface_velocity(coords, vx, vy, vz, n)` | IDW-interpolate DES velocities onto GoSPL mesh |
| `run_and_get_erosion(dt, coords, n, out, k, p)` | Advance GoSPL one step; return `delta_h` at query points |
| `apply_elevation_data(coords, elev, n, k, p)` | Seed GoSPL `hGlobal` from DES surface (called once at init) |
| `interpolate_elevation_to_points(coords, n, out, k, p)` | Query current `hGlobal` at arbitrary coordinates |

## Troubleshooting

Most of these come from a path assumption, fixed by a directory variable in the
Makefile or by changing directory before running.

### Build Issues

| Message | Cause and fix |
|---------|---------------|
| `cannot find -lpython3.11` | The gospl environment is not at `~/miniforge3/envs/gospl` with Python 3.11. Set `CONDA_ENV_PATH` |
| `cannot find -lgospl_extensions` | gospl_extensions is not built at `~/opt/gospl_extensions`. Set `GOSPL_EXT_DIR` |
| `gospl-driver.hpp: No such file` | The `gospl_driver` directory is missing from the source tree |

### Runtime Issues

| Message | Cause and fix |
|---------|---------------|
| `GoSPL not initialized`, or an `ImportError` for `gospl` | The gospl environment is not active: `conda activate gospl` |
| `No module named 'gospl_python_interface'` | `PYTHONPATH` lacks `gospl_extensions/cpp_interface`. Use the `dynearthsol-gospl` wrapper, or export it as in [Run](#run) |
| `The input file is not found`, `Unable to open file` | The YAML path is resolved from the working directory. Run from the YAML's directory or use an absolute path |
| `Error in run_and_get_erosion: error code 77` ... `bad hmax in TSAdaptChoose()` | Intermittent PETSc error in marine deposition; the run recovers. See `examples/README.md` |
| Run is slow | GoSPL is called every step (`gospl_coupling_frequency = 1` by default). Raise it, or use `time` mode |

If GoSPL fails to initialize with a valid path, check that the YAML parses and
has the five required sections (see [GoSPL YAML](#gospl-yaml)).
