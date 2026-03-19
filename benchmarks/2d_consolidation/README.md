# 2D Consolidation Benchmark

## Maintainer

This benchmark workflow was prepared by Sungho Lee.

- Affiliation: Korea Institute of Geoscience and Mineral Resources (KIGAM)
- Email: `slee91@kigam.re.kr`

This directory contains the Mandel consolidation benchmark for the 2D
DynEarthSol executable. The numerical solution is compared against the
analytical Mandel solution at the drainage boundary monitor point used by the
Python plotting script.

The setup solved here is a quarter-domain Mandel problem. Drainage is applied
only on the right boundary (`x1`, Dirichlet pore-pressure boundary), while the
remaining boundaries are treated as no-flow in the hydraulic problem.

## Case

- `2d-consolidation.cfg`
  - Mandel consolidation benchmark in the 2D executable.
  - Quarter-domain setup with drainage only on `x1`.

Case metadata used by the helper scripts lives in `benchmark_cases.py`.

## Build

Build the 2D executable from the repository root:

```bash
make ndims=2 -j2
```

## Run

From the repository root:

```bash
python benchmarks/2d_consolidation/run_2d_consolidation.py --case all
```

This will:

- run the selected benchmark cases,
- write outputs under `benchmarks/2d_consolidation/runs/`,
- generate a notebook-style `final_comparison.png`,
- print the first / max / final absolute error summary in bar.

## Plot Only

To regenerate the plot from an existing run without rerunning the solver:

```bash
python benchmarks/2d_consolidation/run_2d_consolidation.py --plot-only
```

You can also call the plotting script directly:

```bash
python benchmarks/2d_consolidation/plot_2d_consolidation.py \
  --case mandel \
  --run-dir benchmarks/2d_consolidation/runs/mandel \
  --output benchmarks/2d_consolidation/runs/mandel/final_comparison.png
```

## Output Layout

Each case is written to its own run directory:

- `benchmarks/2d_consolidation/runs/mandel/`

Each run directory contains:

- `<model>.info`
- `<model>.save.*`
- `<model>.chkpt.*`
- `run.log`
- `final_comparison.png`

## Recommended Workflow

`plot_2d_consolidation.py` is the canonical plotting path. Use the Python
scripts as the reproducible benchmark workflow:

1. run `run_2d_consolidation.py`
2. inspect `final_comparison.png`
3. regenerate plots with `--plot-only` or `plot_2d_consolidation.py` when needed
