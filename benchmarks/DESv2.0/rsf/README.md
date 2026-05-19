# RSF Benchmark and Example Reproduction Guide

This README summarizes the rate-and-state friction (RSF) example and benchmark cases in DynEarthSol and provides the basic commands needed to reproduce the corresponding results.

The RSF materials are located in:

```text
examples/rate_and_state_friction
benchmarks/simple_shear_rsf
```

The example cases are intended for quick verification. The distributed example mesh is a low-resolution geometry and should be refined before being used for paper-level production simulations.

---

## Quick Start

From the root directory of the DynEarthSol repository, first make sure that the 2-D executable is available:

```bash
ls dynearthsol2d
```

If the executable is located elsewhere, set its path manually:

```bash
export DYNEXE=/path/to/dynearthsol2d
```

Run the RSF examples:

```bash
cd examples/rate_and_state_friction
python3 run_rsf_examples.py --exe ../../dynearthsol2d
python3 plot_rsf_time_series.py -o rsf_time_series.png
```

Return to the repository root and run the simple-shear RSF benchmark:

```bash
cd benchmarks/simple_shear_rsf
python3 run_simple_shear_benchmark.py --clean --output-step-interval 40
python3 plot_simple_shear_benchmark.py -o simple_shear_benchmark.png
```

The expected figure outputs are:

```text
examples/rate_and_state_friction/rsf_time_series.png
benchmarks/simple_shear_rsf/simple_shear_benchmark.png
```

---

## 1. Rate-and-State Friction Examples

### Directory

```text
examples/rate_and_state_friction
```

### Purpose

This directory contains 2-D RSF example cases for quick verification of the RSF implementation.

### Main files

```text
ep_rsf_creep.cfg
ep_rsf_creep_largeDc.cfg
ep_rsf_eq.cfg
rsf_long_strike_model_low_resolution.poly
run_rsf_examples.py
plot_rsf_time_series.py
```

### Example cases

- `ep_rsf_creep.cfg`: velocity-strengthening creep case
- `ep_rsf_eq.cfg`: velocity-weakening earthquake-cycle case
- `ep_rsf_creep_largeDc.cfg`: velocity-weakening case with large `Dc`, producing creep-like behavior

### Output

The simulations generate monitor/output files according to the `modelname` settings in each `.cfg` file. The plotting script generates:

```text
rsf_time_series.png
```

---

## 2. Optional 100 m Fault-Zone Geometry for Higher-Resolution Runs

The example geometry file

```text
examples/rate_and_state_friction/rsf_long_strike_model_low_resolution.poly
```

contains a low-resolution RSF/seismogenic layer bounded by:

```text
z = -4000 m
z = -6000 m
```

This corresponds to a fault-zone thickness of 2 km. For higher-resolution, paper-level simulations, the fault-zone thickness can be reduced to 100 m while preserving the original center depth of 5 km:

```text
fault-zone top    = -4950 m
fault-zone bottom = -5050 m
```

### Recommended workflow

Create a modified copy of the example `.poly` file:

```bash
cd examples/rate_and_state_friction
cp rsf_long_strike_model_low_resolution.poly rsf_long_strike_model_100m.poly
```

Then edit the copied file as described below.

### Node-coordinate changes

Replace the original top and bottom coordinates of the RSF/seismogenic layer with the following values:

```text
1 0.0 -4950.0 # node 1 (z=-4.95 km)
2 0.0 -5050.0 # node 2 (z=-5.05 km)
5 100000.0 -5050.0 # node 5 (z=-5.05 km)
6 100000.0 -4950.0 # node 6 (z=-4.95 km)
8 45000.0 -4950.0 # node 8 (xL=45 km on z=-4.95 km)
9 45000.0 -5050.0 # node 9 (xL=45 km on z=-5.05 km)
10 55000.0 -4950.0 # node 10 (xR=55 km on z=-4.95 km)
11 55000.0 -5050.0 # node 11 (xR=55 km on z=-5.05 km)
```

The horizontal positions `xL = 45 km` and `xR = 55 km` define the left and right boundaries of the velocity-weakening segment. These horizontal coordinates do not need to be changed when only the fault-zone thickness is modified.

### Segment definitions

The segment connectivity can remain unchanged because the same node IDs are used. The comments may be updated to reflect the new layer depths:

```text
8 1 8 0 # z=-4.95 km: 0 -> xL (internal)
9 8 10 0 # z=-4.95 km: xL -> xR (internal)
10 10 6 0 # z=-4.95 km: xR -> 100 km (internal)
11 2 9 0 # z=-5.05 km: 0 -> xL (internal)
12 9 11 0 # z=-5.05 km: xL -> xR (internal)
13 11 5 0 # z=-5.05 km: xR -> 100 km (internal)
14 8 9 0 # vertical split at xL (internal)
15 10 11 0 # vertical split at xR (internal)
```

### Regional mesh-size settings

For a 100 m-thick fault zone, at least 4-5 elements should span the fault-zone thickness. This corresponds to a target element spacing of approximately 20-25 m inside the RSF/seismogenic layer.

The current RSF example uses:

```ini
resolution = 1e2
smallest_size = 1.e0
```

With this base resolution, a practical starting point is to use a smaller regional size value for the three RSF/seismogenic regions:

```text
0 50000.0 -2000.0 4 4.0  # upper elastic (material 4)
1 22500.0 -5000.0 3 0.02 # seismogenic VS-left (material 3)
2 50000.0 -5000.0 2 0.02 # seismogenic VW-mid (material 2)
3 77500.0 -5000.0 1 0.02 # seismogenic VS-right (material 1)
4 50000.0 -8000.0 0 4.0  # lower elastic (material 0)
```

The region points at `z = -5000 m` remain inside the modified 100 m-thick fault zone. The surrounding elastic regions can remain coarser, unless the target production simulation requires additional refinement outside the fault zone.

After generating the mesh, inspect the element distribution and confirm that the 100 m-thick RSF layer is resolved by at least 4-5 elements. If the layer is under-resolved, further reduce the regional size values for the RSF/seismogenic materials.

### Configuration-file changes

After creating the modified `.poly` file, update the `poly_filename` entry in each RSF example configuration file that should use the 100 m geometry:

```ini
poly_filename = ./rsf_long_strike_model_100m.poly
```

For example, update this entry in:

```text
ep_rsf_creep.cfg
ep_rsf_eq.cfg
ep_rsf_creep_largeDc.cfg
```

Also check the monitor-point coordinates. For the 100 m-thick fault-zone geometry described above, the upper host-rock/fault-zone interface is located at 4.95 km depth. A monitor point on this interface can be specified as:

```ini
points_y = [-4.95]
points_unit = km
```

If the lower host-rock/fault-zone interface is needed instead, use `points_y = [-5.05]`.

---

## 3. Simple-Shear RSF Benchmark

### Directory

```text
benchmarks/simple_shear_rsf
```

### Purpose

This directory contains the simple-shear RSF benchmark suite. The benchmark uses a base configuration template and generates multiple benchmark cases for testing the RSF implementation.

### Main files

```text
simple_shear_base.cfg
run_simple_shear_benchmark.py
plot_simple_shear_benchmark.py
check_simple_shear_benchmark.py
```

### Output

The benchmark run generates:

```text
runs/
simple_shear_benchmark.png
```

The `runs/` directory contains the generated benchmark case directories. The file `simple_shear_benchmark.png` summarizes the benchmark results.

---

## 4. Optional Numerical Checks

A smaller functional test can be run with:

```bash
cd benchmarks/simple_shear_rsf
python3 check_simple_shear_benchmark.py
```

The benchmark result can also be checked numerically without saving a figure:

```bash
cd benchmarks/simple_shear_rsf
python3 plot_simple_shear_benchmark.py --check --no-save --max-relative-error 5e-2
```

---

## 5. Requirements

The examples and benchmarks require:

- A compiled DynEarthSol 2-D executable, typically named `dynearthsol2d`
- Python 3
- `numpy`
- `matplotlib`

If needed, install the Python packages with:

```bash
python3 -m pip install numpy matplotlib
```

---

## 6. Notes

The RSF examples and simple-shear benchmark are self-contained in the directories listed above. No additional external input files are required beyond the compiled DynEarthSol executable and the files included in these directories.

The example `.poly` file is provided for low-resolution verification. For production or paper-level simulations, use an appropriately refined mesh and, if needed, replace the 2 km-thick example fault zone with a thinner geometry such as the 100 m fault-zone option described above.
