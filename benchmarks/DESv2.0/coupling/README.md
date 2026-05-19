# DES–GoSPL Coupled Model Benchmark

To reproduce the demo coupled model, follow the full instructions in
[`gospl_driver/README.md`](../../../gospl_driver/README.md); or follow the [Quick Start](#quick-start)

For more information about the demo, 
[`gospl_driver/examples/README.md`](../../../gospl_driver/examples/README.md).

## Quick Start

### Prerequisites

**GoSPL Extensions Library**: You need to clone and build the `gospl_extensions` repository: e.g.,
   ```bash
   git clone https://github.com/GeoFLAC/gospl_extensions.git ~/opt/gospl_extensions
   cd ~/opt/gospl_extensions/cpp_interface
   conda activate gospl
   make install-local
   ```

   This creates:
   - `gospl_extensions/lib/libgospl_extensions.so` - Shared library for linking
   - `gospl_extensions/include/gospl_extensions.h` - Header file for inclusion

   The `install-local` target sets up the standard library structure that DynEarthSol's Makefile expects.

   **Note**: This step requires the GoSPL conda environment to be set up first: Follow [the goSPL installation procedure](https://gospl.readthedocs.io/en/latest/getting_started/installConda.html)

### Run the coupled model

1. `cd DynEarthSol/gospl_driver/examples`
2. `../../dynearthsol-gospl ./gaussian-weakzone-3d-with-gospl.cfg`

### Run the standalone DES model

1. `cd DynEarthSol/gospl_driver/examples`
2. Change `surface_process_option` to `1` (simple diffusion) in `gaussian-weakzone-3d-with-gospl.cfg`.
3. `../../dynearthsol ./gaussian-weakzone-3d-with-gospl.cfg`

