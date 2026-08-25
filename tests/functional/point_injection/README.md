# Fluid Point-Injection Schedule Check

This functional check runs the same 2D fluid source with two schedule models:
one prescribes a constant rate and the other the equivalent total amount. Their
monitor pore-pressure histories must agree to roundoff, remain finite, and show
a positive center-pressure effect relative to an otherwise identical
source-disabled baseline.

The fixture uses a `1e-6 s` explicit step and a `20e-6 s` source window. This
keeps the public mechanical update below its elastic stability limit; the older
private fixture's `1 s` step is unstable on this public solver even when the
source is disabled, despite remaining below the hydraulic diffusion limit.

From the repository root:

```bash
python tests/functional/point_injection/run_point_injection.py --clean
```

The comparison runs by default. Use `--skip-check` only to retain independent
case outputs for manual inspection. Outputs are written below
`tests/functional/point_injection/runs/`.
