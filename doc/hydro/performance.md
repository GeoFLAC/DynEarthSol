# Cost of finite-density transport

The optional finite-density model adds local pressure/temperature state
conversion and extensive mass/heat bookkeeping. Its cost should be compared
with that same model before optimization. Default linear flow and mechanics
keep their existing integration and arithmetic.

The measured optimization reuses the element gradients, density and mobility
in the shared Darcy flux routine. It parallelizes independent nodal state
conversions only for at least 256 nodes. The existing transport workspace
owns and reuses its arrays; no per-node allocation is added. The signed
interior flux, explicit substep update, local roots and acceptance test
retain their mathematical definitions.

## Paired measurements

Each case alternates executable order, discards one warmup pair and reports
the median of three further pairs. CPU cases use eight OpenMP threads and
1,000 steps; H100 cases use one host thread and 200 steps. Values below are
the timed stepping interval in seconds. These are small tests on shared
hosts, so the CPU and GPU columns are not a backend speed comparison.

| Case | CPU before | CPU after | H100 before | H100 after |
| --- | --- | --- | --- | --- |
| 2D mass transport | 4.236 | 2.628 | 2.306 | 1.968 |
| 2D mass and heat | 8.920 | 3.714 | 4.173 | 2.234 |
| 3D mass transport | 8.487 | 4.966 | 4.054 | 2.791 |
| 3D mass and heat | 17.883 | 7.690 | 7.809 | 3.774 |

The optional model is 1.6–2.4 times faster on these eight-thread CPU cases
and 1.2–2.1 times faster on H100. Single-thread CPU mass/heat timings are
roughly unchanged, with a 6% slowdown in the small 2D mass case. There is
no claim that the finite-density model costs the same as linear flow.
Default-model timings are too short to establish differences of a few
percent or general application performance.

CPU physical arrays match bit for bit across executable variants, one/eight
threads and repeats. The 3D H100 arrays also match exactly. The 2D H100 heat
case changes the last bits when local thermodynamic functions execute on
the device: the maximum normalized differences are $1.7\times10^{-15}$
for temperature, $2.5\times10^{-12}$ for pressure and
$2.3\times10^{-11}$ for velocity. Its extra bitwise test fails and is
retained. It passes DES's unchanged normalized maximum-plus-standard-
deviation limit of $10^{-8}$, along with larger-grid physical mass and
sensible-energy residual checks. Default mechanics and linear-flow arrays
remain bit-exact in these paired controls.

Scientific checks remain separate from timing: both dimensions, native and
HDF5 CPU outputs, H100, mass/heat, transport, remapping, restart and undamped
waves must pass before adopting the optimization. Those controls do not
establish performance on arbitrary meshes, retentive laws or complete
earthquake cycles. A100 and A6000 have actual-device kernel and parent-core
controls; no performance speedup is inferred for them from the H100 timing.

## Reproduce

Use the paired input decks and executable/source hashes supplied with the
review's performance evidence. Run each input four times per executable,
alternating order and discarding the first pair:

```bash
OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  dynearthsol2d input.cfg
```

Use the 3D executable for the cube. For the GPU timing, select one device,
set `OMP_NUM_THREADS=1`, and use the 200-step version of the same input.
Read the difference between the final and initial `walltime` in the model's
`.info` file; exclude initialization. Preserve outputs in separate folders
for the physical-array and budget comparisons.
