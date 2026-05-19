# DynEarthSol Performance Benchmark for DESv2.0

This benchmark measures compute throughput (`s/step`) across a range of problem sizes and hardware configurations. Two mesh types are tested (regular and irregular), and a separate CPU core-scaling test is included.

---

## Directory layout

```
performance/
├── test-regular-3d.cfg       base config for regular mesh
├── test-irregular-3d.cfg     base config for irregular mesh
├── setup_benchmarks.sh       generates the three test folders below
├── 01-regular/               regular mesh, CPU + GPU throughput sweep
├── 02-irregular/             irregular mesh, CPU + GPU throughput sweep
└── 03-scaling/               regular mesh, CPU core-scaling test
```

Each test folder contains per-size config files (`test-*-NN-<size>.cfg`), the executables, and ready-to-run shell scripts.

---

## Step 1 — Build the executables

From the DynEarthSol repository root:

### CPU (OpenMP)

```bash
make usemmg=1
```

Produces `dynearthsol3d`. We recommend adding `hdf5=1` to enable HDF5-based vtkhdf output.

### GPU (OpenACC, NVHPC required)

```bash
make openacc=1 usemmg=1
```

Produces `dynearthsol3d.gpu`. We recommend adding `hdf5=1` to enable HDF5-based vtkhdf output.

---

## Step 2 — Generate test folders

Run once from this directory after both executables are built:

```bash
cd benchmarks/DESv2.0/performance
bash setup_benchmarks.sh
```

This creates `01-regular/`, `02-irregular/`, and `03-scaling/`, each pre-populated with numbered config files and copies of both executables.

---

## Step 3 — Run the benchmarks

### 01-regular — regular mesh throughput (8 problem sizes)

```bash
cd 01-regular

# CPU (set thread count to number of physical cores)
OMP_NUM_THREADS=64 ./run-cpu.sh

# GPU (uses full device automatically)
./run-gpu.sh
```

Config files run in ascending problem-size order:
`01-2e3` → `02-4e3` → … → `08-256e3`

### 02-irregular — irregular mesh throughput (8 problem sizes)

```bash
cd 02-irregular
OMP_NUM_THREADS=64 ./run-cpu.sh
./run-gpu.sh
```

### 03-scaling — CPU core-scaling (4 problem sizes × 7 thread counts)

```bash
cd 03-scaling
./run.sh
```

Runs `OMP_NUM_THREADS` = 1, 2, 4, 8, 16, 32, 64 for each of the 4 sizes (`2e3`–`16e3`).

---

## Step 4 — Collect results

Each run script writes a `.log` file per case. The key metric is `s/step` (compute time per time step, excluding remeshing and output), printed at the end of each run:

```
  Compute : hh:mm:ss (XX.XX%)/ <steps> = <value> s/step
```

Extract all results at once:

```bash
# throughput tests
grep 's/step' 01-regular/*.log
grep 's/step' 02-irregular/*.log

# scaling test
grep 's/step' 03-scaling/*.log
```

Log file naming:
- `*-cpu.log` — CPU throughput runs
- `*-gpu.log` — GPU throughput runs
- `*-t<N>.log` — scaling runs at N threads

---

## Hardware reference

| Resource | Count |
|---|---|
| GPU (NVIDIA) | 72 SMs / full device |
| CPU (AMD EPYC) | 64 threads |
