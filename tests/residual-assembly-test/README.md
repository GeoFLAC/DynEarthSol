# Undamped force assembly

Run from this directory:

```sh
make check ndims=2 backend=asan
make check ndims=3 backend=asan
OMP_NUM_THREADS=4 make check ndims=2 backend=omp
OMP_NUM_THREADS=4 make check ndims=3 backend=omp
make check ndims=2 backend=acc
make check ndims=3 backend=acc
```

OpenACC requires NVHPC, the production build dependencies and a GPU. `GPU_CC`
defaults to the detected device and can be overridden. Accelerator checks use the
root build graph; do not run them concurrently with another root build. The test
entry point replaces only application `main`; all physical routines are real.

Two positively oriented simplices share one vertex. Checks cover cancellation,
Neumann traction, gravity with porous bulk density, and artificial damping.
The optional output is the unprojected assembled physical force before damping.
Supplying it must preserve both existing outputs bit for bit. Analytic assertions
use `abs(x-y) <= 1e-9 + 1e-14*abs(x) + 1e-14*abs(y)` in the fixture's force units.

`mode=BASE` calls the old five-argument API and measures the legacy residual
overwrite defect. The default `mode=B0` exercises the optional undamped output.
These tests do not qualify a corrected PT solver or constraint projection.

The initial fixture was authored by Qwen3.8-27B-FP8 under supervision. Codex added
the real material/body-force fixture, null-output equivalence checks and build
drivers, and independently validated the baseline and changed implementations.
