# Simple-Shear RSF Benchmark

## Maintainer

This benchmark workflow was prepared by Sungho Lee.

- Affiliation: Korea Institute of Geoscience and Mineral Resources (KIGAM)
- Email: `slee91@kigam.re.kr`

## Files

- `simple_shear_base.cfg`: base CFG template
- `run_simple_shear_benchmark.py`: generate and run benchmark cases
- `plot_simple_shear_benchmark.py`: create a figure or run a numerical check
- `check_simple_shear_benchmark.py`: run the small CI subset without plotting

## Full Run

```bash
cd benchmarks/simple_shear_rsf
python3 run_simple_shear_benchmark.py --clean --output-step-interval 40
python3 plot_simple_shear_benchmark.py -o simple_shear_benchmark.png
```

Generated outputs:

- case directories under `runs/`
- figure file `simple_shear_benchmark.png`

## Check Only

Run the small RSF subset used for functional testing:

```bash
cd benchmarks/simple_shear_rsf
python3 check_simple_shear_benchmark.py
```

Default check settings:

- cases: `steady_ab_pos_v0_1e-6`, `aging_ab_neg_dc_1e-6`
- relative error tolerance: `5e-2`
- output directory: temporary and auto-removed

## Manual Check Without PNG

```bash
cd benchmarks/simple_shear_rsf
python3 plot_simple_shear_benchmark.py --check --no-save --max-relative-error 5e-2
```
