# GoSPL Integration for DynEarthSol

This directory contains the integration code for GoSPL (Geomorphological Landscape Evolution) surface processes within DynEarthSol.

## Files

- `gospl-driver.hpp` - Header file for the GoSPL C++ interface wrapper
- `gospl-driver.cxx` - Implementation of the GoSPL C++ interface wrapper
- `examples/` - Example configuration files and documentation

## Prerequisites

**GoSPL Extensions Library**: You need to clone and build the gospl_extensions repository:
   ```bash
   git clone https://github.com/GeoFLAC/gospl_extensions.git
   cd gospl_extensions/cpp_interface
   
   # Activate gospl environment for building
   conda activate gospl
   
   # Build and install locally for DynEarthSol integration
   make install-local
   cd ../..
   ```

   This creates:
   - `gospl_extensions/lib/libgospl_extensions.so` - Shared library for linking
   - `gospl_extensions/include/gospl_extensions.h` - Header file for inclusion

   The `install-local` target sets up the standard library structure that DynEarthSol's Makefile expects.

   **Note**: This step requires the GoSPL conda environment to be set up first: Follow [the goSPL installation procedure](https://gospl.readthedocs.io/en/latest/getting_started/installConda.html)
   

## Build and Runtime Workflow

The GoSPL integration provides two ways to run DynEarthSol with GoSPL support:

### Method 1: Use the Auto-Generated Wrapper (Recommended)

The Makefile automatically creates a wrapper script that handles all environment setup:

1. **Build DynEarthSol with GoSPL support**:
   ```bash
   # Make sure use_gospl = 1 in the Makefile (should already be set)
   make clean
   make -j4
   ```
   
   This automatically creates a `dynearthsol-gospl` wrapper script.

2. **Run using the wrapper** (environment automatically set):
   ```bash
   # Activate gospl environment
   conda activate gospl
   
   # Run using the wrapper (PYTHONPATH is set automatically)
   ./dynearthsol-gospl your_config.cfg
   ```

### Method 2: Manual Environment Setup

If you prefer manual control or the wrapper doesn't work:

1. **Build outside gospl environment** (optional but recommended):
   ```bash
   # Make sure gospl environment is NOT activated for building
   conda deactivate  # if any environment is active
   
   # Set use_gospl = 1 in the Makefile (should already be set)
   make clean
   make -j4
   ```

2. **Run with manual environment setup**:
   ```bash
   # Activate gospl environment before running
   conda activate gospl
   
   # Set PYTHONPATH manually to include gospl_extensions interface
   export PYTHONPATH="$HOME/opt/gospl_extensions/cpp_interface:${PYTHONPATH}"
   
   # Run DynEarthSol with GoSPL-enabled configuration
   ./dynearthsol3d your_config.cfg
   ```

### Configuration
3. In your DynEarthSol configuration file, set:
   ```
   control.surface_process_option = 11
   control.surface_process_gospl_config_file = "/absolute/path/to/gospl_config.yml"
   ```
   **Important**: Use absolute paths to avoid file access issues.

4. Create a GoSPL configuration file compatible with gospl_extensions:
   - See working examples in `/home/echoi2/opt/gospl_extensions/examples/`
   - Note: The config format may differ from standard GoSPL configs
   - Use `examples/gospl_config_example.yml` as a starting template, but verify compatibility

## Integration Details

The GoSPL integration works as follows:

1. **Initialization**: GoSPL is initialized once from the YAML config file. Its elevation field (`hGlobal`) is seeded from DES's initial surface via `apply_elevation_data()`.
2. **Each coupling event** (every `gospl_coupling_frequency` DES steps, or every `gospl_coupling_interval_in_yr` years when `gospl_coupling_mode = time`):
   - **Time-averaged tectonic velocity** is computed as `Δcoord/Δt` over the coupling interval (not the instantaneous DES velocity). DES uses inertial scaling (quasi-dynamic formulation), so instantaneous velocities contain damped-wave components that are numerical artifacts. Averaging over the interval filters these out. On the first coupling event or after remeshing, instantaneous velocity is used as a fallback.
   - DES surface velocities (vx, vy, vz) are IDW-interpolated onto the GoSPL mesh via `set_surface_velocity()`.
   - `run_and_get_erosion(dt)` advances GoSPL by one step of length `dt`. Internally: horizontal advection (vx, vy), vertical uplift (vz via `upsub`), SPL erosion, and hillslope diffusion are applied to GoSPL's `hGlobal`. The returned `delta_h` contains **only the erosion and diffusion component** — the uplift (`upsub * dt`) is subtracted before returning because DES already applied the same displacement through its Lagrangian mechanical solver. Returning the full delta_h would double-count the tectonic uplift.
   - DES adds `delta_h` (erosion + diffusion only) to its surface node z-coordinates.
3. GoSPL **owns** the topography between coupling events and accumulates its drainage network state continuously.

### Coupling API (`EnhancedModel` in `gospl_extensions`)

| Method | Purpose |
|---|---|
| `set_surface_velocity(pts, vx, vy, vz)` | IDW-interpolate DES velocities onto GoSPL mesh |
| `run_and_get_erosion(dt, query_pts)` | Advance GoSPL one step; return `delta_h` at query points |
| `apply_elevation_data(elevdata)` | Reset `hGlobal` from external elevation data (init / post-remesh) |
| `interpolate_elevation_to_points(pts)` | Query current `hGlobal` at arbitrary coordinates |
| `apply_drift_correction(pts, elev, alpha)` | Gently blend `hGlobal` toward DES elevation |

## Known Limitations / Caveats

### GoSPL elevation reset after DES remeshing

**Every time DES remeshes, GoSPL's elevation field is completely reset from the DES surface** (`needs_elevation_reset = true` is set in `dynearthsol.cxx` after each call to `remesh()`). This causes:

- GoSPL's accumulated erosion/deposition patterns are discarded.
- GoSPL's implicit drainage network state (flow accumulation, chi, pit-filling history) is wiped and must rebuild from scratch.
- If DES remeshes frequently (governed by `quality_check_step_interval` and `min_quality`), GoSPL may never accumulate enough time to develop stable stream incision before being reset again.

This reset is necessary because DES node positions change after remeshing, so GoSPL's elevation must re-sync with DES. However, it is a hard discontinuity. A gentler post-remesh correction (e.g., `apply_drift_correction` with `alpha=1.0`) is equivalent in effect but could in principle be tuned. GoSPL's own mesh is fixed and independent of DES's mesh, so only the elevation field — not the connectivity — needs updating after remeshing.

**Practical implication**: To give GoSPL time to develop a stable drainage network, prefer infrequent remeshing (larger `quality_check_step_interval`, looser `min_quality`) during the phase when surface processes are most important.

## Troubleshooting

### Build Issues
- **Error: "cannot find -lpython3.11"**: Make sure the gospl conda environment exists at `~/miniconda3/envs/gospl` and contains Python 3.11. If your conda is elsewhere, update `CONDA_ENV_PATH` in the Makefile.
- **Error: "gospl-driver.hpp: No such file"**: Ensure the gospl_driver directory is in the DynEarthSol source directory.
- **Error: "cannot find -lgospl_extensions"**: Make sure gospl_extensions is built at `~/opt/gospl_extensions`. If elsewhere, update `GOSPL_EXT_DIR` in the Makefile.

### Runtime Issues

#### Using the Wrapper Script (Recommended)
If using `./dynearthsol-gospl`:
- **ImportError for gospl**: Make sure you activated the gospl environment before running (`conda activate gospl`).
- **Wrapper script not found**: Make sure you built with `use_gospl = 1` and the build completed successfully.

#### Manual Environment Setup  
If setting PYTHONPATH manually:
- **ModuleNotFoundError: No module named 'gospl_python_interface'**: Ensure both the gospl environment is activated AND the PYTHONPATH includes the gospl_extensions cpp_interface directory:
  ```bash
  conda activate gospl
  export PYTHONPATH="$HOME/opt/gospl_extensions/cpp_interface:${PYTHONPATH}"
  ./dynearthsol3d your_config.cfg
  ```

#### General GoSPL Issues
- **GoSPL initialization fails**: Check that your GoSPL configuration file path is correct and the file is valid YAML.
- **"The input file is not found" / "Unable to open file" error**: This occurs in the gospl_extensions C++ code when creating the EnhancedModel. **Solutions** (try in order):
  
  1. **Run from the gospl_extensions directory**: The C++ code may expect specific working directory
     ```bash
     cd ~/opt/gospl_extensions/examples
     /path/to/DynEarthSol/dynearthsol-gospl /path/to/your/dynearthsol_config.cfg
     ```
  
  2. **Use a working gospl_extensions config file**: Copy and modify an existing working config
     ```bash
     cp ~/opt/gospl_extensions/examples/input-escarpment.yml ./my_gospl_config.yml
     # Edit my_gospl_config.yml for your simulation parameters
     # Update DynEarthSol config to point to this file
     ```
  
  3. **Check file format compatibility**: The EnhancedModel expects specific YAML structure
     - Standard GoSPL configs may not work with gospl_extensions
     - Required sections may include: `domain`, `time`, `spl`, `diffusion`, `output`
     - See working examples in `~/opt/gospl_extensions/examples/`
  
  4. **Verify file permissions and absolute paths**:
     ```bash
     ls -la /absolute/path/to/gospl_config.yml
     # Ensure file is readable and use absolute paths in DynEarthSol config
     ```

### Verification
- Check if GoSPL options appear in help: `./dynearthsol3d --help | grep gospl`
- Test GoSPL import: `conda activate gospl && python -c "import gospl"`
- Check library linking: `ldd dynearthsol3d | grep python`
- Verify wrapper script: `cat dynearthsol-gospl` (should show PYTHONPATH setup)
