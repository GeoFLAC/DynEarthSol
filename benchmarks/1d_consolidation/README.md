# 1D Consolidation Benchmarks

## Maintainer

This benchmark workflow was prepared by Sungho Lee.

- Affiliation: Korea Institute of Geoscience and Mineral Resources (KIGAM)
- Email: `slee91@kigam.re.kr`

This directory contains 1D-style consolidation benchmarks for both the
2D and 3D DynEarthSol executables. The comparison is made against the
Terzaghi-type analytical solution evaluated at the bottom of the column.

## Cases

- `1d-consolidation-des2d_traction.cfg`
  - Surface traction generates the excess pore pressure through poroelastic
    coupling.
- `1d-consolidation-des2d_water_loading.cfg`
  - Initial excess pore pressure is prescribed directly.
- `1d-consolidation-des3d_traction.cfg`
  - 3D version of the traction-driven 1D column benchmark.
- `1d-consolidation-des3d_water_loading.cfg`
  - 3D version of the prescribed initial excess pore pressure benchmark.

Case metadata used by the helper scripts lives in
`benchmark_cases.py`.

## Build

Build both executables from the repository root:

```bash
make ndims=2 -j2
make ndims=3 -j2
```

## Run All Cases

From the repository root:

```bash
python benchmarks/1d_consolidation/run_1d_consolidation.py --case all
```

This will:

- run the selected benchmark cases,
- write outputs under `benchmarks/1d_consolidation/runs/`,
- generate `final_comparison.png` for each case,
- print the first / max / final relative error summary,
- automatically pick `dynearthsol2d` or `dynearthsol3d` from the repository root
  based on the case dimension.

## Run One Case

Examples:

```bash
python benchmarks/1d_consolidation/run_1d_consolidation.py --case traction
python benchmarks/1d_consolidation/run_1d_consolidation.py --case water_loading
python benchmarks/1d_consolidation/run_1d_consolidation.py --case des3d_traction
python benchmarks/1d_consolidation/run_1d_consolidation.py --case des3d_water_loading
```

If your executables live elsewhere, override them explicitly:

```bash
python benchmarks/1d_consolidation/run_1d_consolidation.py \
  --case all \
  --exe2d /path/to/dynearthsol2d \
  --exe3d /path/to/dynearthsol3d
```

To regenerate plots from existing outputs without rerunning the solver:

```bash
python benchmarks/1d_consolidation/run_1d_consolidation.py --case all --plot-only
```

## Plot Only

You can also call the plotting script directly.

Example for the traction case:

```bash
python benchmarks/1d_consolidation/plot_1d_consolidation.py \
  --case traction \
  --run-dir benchmarks/1d_consolidation/runs/traction \
  --output benchmarks/1d_consolidation/runs/traction/final_comparison.png
```

Example for the 3D traction case:

```bash
python benchmarks/1d_consolidation/plot_1d_consolidation.py \
  --case des3d_traction \
  --run-dir benchmarks/1d_consolidation/runs/des3d_traction \
  --output benchmarks/1d_consolidation/runs/des3d_traction/final_comparison.png
```

The plotting script also supports explicit model names if needed:

```bash
python benchmarks/1d_consolidation/plot_1d_consolidation.py \
  --model terzaghi_traction \
  --run-dir benchmarks/1d_consolidation/runs/traction
```

## Output Layout

Each case is written to its own run directory:

- `benchmarks/1d_consolidation/runs/traction/`
- `benchmarks/1d_consolidation/runs/water_loading/`
- `benchmarks/1d_consolidation/runs/des3d_traction/`
- `benchmarks/1d_consolidation/runs/des3d_water_loading/`

Each run directory contains:

- `<model>.info`
- `<model>.save.*`
- `<model>.chkpt.*`
- `run.log`
- `final_comparison.png`

## Recommended Workflow

Use the Python scripts as the canonical benchmark workflow:

1. run `run_1d_consolidation.py`
2. inspect `final_comparison.png`
3. regenerate plots with `--plot-only` or `plot_1d_consolidation.py` when needed
