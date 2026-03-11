# Rate-And-State Friction Examples (Replacement Set)

This directory contains the replacement RSF example set for paper-style verification.

Included cases:
- `ep_rsf_creep.cfg`: velocity-strengthening creep
- `ep_rsf_eq.cfg`: velocity-weakening earthquake cycles
- `ep_rsf_creep_largeDc.cfg`: velocity-weakening with large `Dc` (creep-like)

Utility scripts:
- `run_rsf_examples.py`: run the three CFGs sequentially
- `plot_rsf_time_series.py`: two-panel verification plot from monitor CSV output

## Quick Start

```bash
cd examples/rate_and_state_friction
python3 run_rsf_examples.py
python3 plot_rsf_time_series.py -o rsf_time_series.png
```

## Requirements

- DynEarthSol executable (`dynearthsol2d` or `dynearthsol3d`) in repo root,
  or set `DYNEXE`, or pass `--exe`.
- Python packages for plotting:

```bash
python3 -m pip install numpy matplotlib
```
