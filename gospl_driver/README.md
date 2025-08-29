# GoSPL Integration for DynEarthSol

This directory contains the integration code for GoSPL (Geomorphological Landscape Evolution) surface processes within DynEarthSol.

## Files

- `gospl-driver.hpp` - Header file for the GoSPL C++ interface wrapper
- `gospl-driver.cxx` - Implementation of the GoSPL C++ interface wrapper
- `examples/` - Example configuration files and documentation

## Prerequisites

1. **GoSPL Extensions Library**: You need to clone and build the gospl_extensions repository:
   ```bash
   git clone https://github.com/GeoFLAC/gospl_extensions.git
   cd gospl_extensions
   make
   cd ..
   ```

2. **GoSPL Python Package**: Install GoSPL in a conda environment:
   ```bash
   conda create -n gospl gospl -c conda-forge
   conda activate gospl
   ```

## Build and Runtime Workflow

The GoSPL integration uses a **build outside, run inside** workflow:

### Building DynEarthSol
1. **Build outside gospl environment** (recommended):
   ```bash
   # Make sure gospl environment is NOT activated
   conda deactivate  # if any environment is active
   
   # Set use_gospl = 1 in the Makefile (should already be set)
   make clean
   make -j4
   ```

   The Makefile automatically detects the gospl environment path and links to the necessary Python and Boost libraries from the gospl environment, even when building outside of it.

### Running DynEarthSol with GoSPL
2. **Run inside gospl environment** (required):
   ```bash
   # Activate gospl environment before running
   conda activate gospl
   
   # Run DynEarthSol with GoSPL-enabled configuration
   ./dynearthsol3d your_config.cfg
   ```

### Configuration
3. In your DynEarthSol configuration file, set:
   ```
   control.surface_process_option = 11
   control.surface_process_gospl_config_file = "/path/to/gospl_config.yml"
   ```

4. Create a GoSPL configuration file (see `examples/gospl_config_example.yml`)

## Integration Details

The GoSPL integration works by:
1. Initializing the GoSPL model with a YAML configuration file
2. At each time step, passing surface velocities from DynEarthSol to GoSPL
3. Running GoSPL landscape evolution for the time step
4. Updating DynEarthSol surface elevations with changes from GoSPL

This allows for two-way coupling between geodynamic processes (DynEarthSol) and surface processes (GoSPL).

## Troubleshooting

### Build Issues
- **Error: "cannot find -lpython3.11"**: Make sure the gospl conda environment exists and contains Python 3.11. Check the CONDA_ENV_PATH in the Makefile points to the correct environment.
- **Error: "gospl-driver.hpp: No such file"**: Ensure you're including `gospl_driver/gospl-driver.hpp` with the correct path.

### Runtime Issues
- **ImportError for gospl**: Make sure you activated the gospl environment before running (`conda activate gospl`).
- **GoSPL initialization fails**: Check that your GoSPL configuration file path is correct and the file is valid YAML.
- **Python path issues**: The gospl environment must be activated to ensure all Python modules are found.

### Verification
- Check if GoSPL options appear in help: `./dynearthsol3d --help | grep gospl`
- Test GoSPL import: `conda activate gospl && python -c "import gospl"`
- Check library linking: `ldd dynearthsol3d | grep python`
