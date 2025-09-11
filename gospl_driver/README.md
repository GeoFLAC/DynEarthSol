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
   
   # Set PYTHONPATH to include gospl_extensions interface
   export PYTHONPATH="/home/echoi2/opt/gospl_extensions/cpp_interface:${PYTHONPATH}"
   
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
- **ModuleNotFoundError: No module named 'gospl_python_interface'**: This error occurs when the Python interface module from gospl_extensions cannot be found. **Solution**: Ensure both the gospl environment is activated AND the PYTHONPATH includes the gospl_extensions cpp_interface directory:
  ```bash
  conda activate gospl
  export PYTHONPATH="/home/echoi2/opt/gospl_extensions/cpp_interface:${PYTHONPATH}"
  ./dynearthsol3d your_config.cfg
  ```
- **GoSPL initialization fails**: Check that your GoSPL configuration file path is correct and the file is valid YAML.
- **"The input file is not found" / "Unable to open file" error**: This occurs in the gospl_extensions C++ code when creating the EnhancedModel. **Solutions** (try in order):
  
  1. **Run from the gospl_extensions directory**: The C++ code may expect specific working directory
     ```bash
     cd /home/echoi2/opt/gospl_extensions/examples
     /home/echoi2/opt/DynEarthSol/dynearthsol3d /path/to/your/dynearthsol_config.cfg
     ```
  
  2. **Use a working gospl_extensions config file**: Copy and modify an existing working config
     ```bash
     cp /home/echoi2/opt/gospl_extensions/examples/input-escarpment.yml ./my_gospl_config.yml
     # Edit my_gospl_config.yml for your simulation parameters
     # Update DynEarthSol config to point to this file
     ```
  
  3. **Check file format compatibility**: The EnhancedModel expects specific YAML structure
     - Standard GoSPL configs may not work with gospl_extensions
     - Required sections may include: `domain`, `time`, `spl`, `diffusion`, `output`
     - See working examples in `/home/echoi2/opt/gospl_extensions/examples/`
  
  4. **Verify file permissions and absolute paths**:
     ```bash
     ls -la /absolute/path/to/gospl_config.yml
     # Ensure file is readable and use absolute paths in DynEarthSol config
     ```
- **Python path issues**: The gospl environment must be activated to ensure all Python modules are found.

### Verification
- Check if GoSPL options appear in help: `./dynearthsol3d --help | grep gospl`
- Test GoSPL import: `conda activate gospl && python -c "import gospl"`
- Check library linking: `ldd dynearthsol3d | grep python`
