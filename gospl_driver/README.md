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

   **Note**: This step requires the GoSPL conda environment to be set up first:
   ```bash
   conda create -n gospl gospl -c conda-forge
   conda activate gospl
   ```

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

**Note**: The integration automatically uses the `EnhancedModel.run_one_step()` method internally to synchronize GoSPL time steps with DynEarthSol's geodynamic time stepping, ensuring precise coupling between surface and deep processes.

## Integration Details

The GoSPL integration works by:
1. Initializing the GoSPL model with a YAML configuration file
2. At each time step, passing surface velocities from DynEarthSol to GoSPL
3. Running GoSPL landscape evolution for the time step
4. Updating DynEarthSol surface elevations with changes from GoSPL

This allows for two-way coupling between geodynamic processes (DynEarthSol) and surface processes (GoSPL).

## Enhanced Model Features

The integration uses an **EnhancedModel** class that extends the standard GoSPL Model with additional methods for granular control over simulation time steps:

### run_one_step Method

The `EnhancedModel` class includes a `run_one_step` method that allows executing exactly one simulation time step instead of running the full simulation loop. This provides fine-grained control needed for coupling with external models like DynEarthSol.

**Key Features:**
- **Single Time Step Execution**: Runs exactly one simulation time step
- **Parameter Control**: Accepts optional `dt` parameter to override default time step
- **State Management**: Properly saves/restores original time step values
- **Full Process Coverage**: Includes all surface processes (flow accumulation, erosion, deposition, tectonics, etc.)

**Usage Example (Python):**
```python
from gospl_model_ext.enhanced_model import EnhancedModel

# Create enhanced model
model = EnhancedModel('gospl_config.yml')

# Run one step with default dt
model.run_one_step()

# Run one step with custom dt (1000 years)
model.run_one_step(dt=1000)

# Use the higher-level wrapper with timing info
elapsed_time = model.runProcessesForDt(dt=500, verbose=True)
```

**What happens in one time step:**
1. Output and visualization updates
2. Stratigraphy management (new stratal layers when needed)
3. Tectonics (advection and tectonic forces)
4. Surface processes (if not in fast mode):
   - Flow accumulation computation
   - River incision (Stream Power Law)
   - Sediment deposition (inland and marine)
   - Hillslope diffusion
5. Stratigraphic compaction (when needed)
6. Flexural isostasy (if enabled)
7. External force updates
8. Time advancement by `dt`

This granular control enables precise synchronization between DynEarthSol's geodynamic time steps and GoSPL's surface process evolution.

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
