# Provenance

DES records where each build, run and frame came from:

| What | Where |
|---|---|
| The build | in the executable, as a `build.snapshot` block |
| The run | in `<modelname>.manifest`, beside the model's output |
| Each frame | in the frame: a `/provenance` group (HDF5), a `provenance` record (des-binary) |

Unknown values read `"unknown"`, `-1` or `0`.

## The build

The `build.snapshot` block holds the revision, make options, dependency versions and providers,
toolchain, and the compile and link flags. Read it without running the executable:

```bash
strings <exe> | grep '^build\.snapshot\.'
```

A toolchain that also puts the block in `.debug_*` may print it more than once; the copies are
identical.

### Uncommitted changes

`make snapshot_diff=1` (off by default) also embeds the working-tree diff:

```bash
strings <exe> | sed -n '/^build\.code-changes\.begin :$/,/^build\.code-changes\.end   :$/p'
```

* The diff covers `*.c *.h *.cxx *.hpp *.cpp *.cu`, `Makefile` and `3x3-C/Makefile`, tracked in
  `HEAD` only.
* Other changed files are counted, not embedded.
* The `-dirty` suffix of `rev` and the `dirty=` counts use the same file patterns.

## The run: `<modelname>.manifest`

The manifest uses the cfg format, with these sections in this order:

| Section | What it holds |
|---|---|
| `[runtime.model]` | the model name, the start time, and on a restart the model and frame it read |
| `[runtime.host]` | the OS, the user@host running it, the CPU and its cores, the memory |
| `[runtime.device]` | `kernel` (`CPU` or `GPU`); on a GPU run, the device, its driver and memory |
| `[runtime.threads]` | whether OpenMP is on, the team size, and the wait policy |
| `[runtime.env]` | the thread and device variables that are set, DES's own macOS default included |
| `[build.*]` | the executable's `build.snapshot` block, and `exe_mtime_utc`, the file's own mtime |
| `build.code-changes` | under `snapshot_diff=1` the diff, else a line saying it was not embedded |

* `omp_threads` is the team size as measured, not read from `OMP_NUM_THREADS`. `OMP_NUM_THREADS`
  sets the team; unset, OpenMP takes every logical CPU, which on a hybrid CPU includes the
  efficiency cores.
* `omp_wait_policy_src` says who set the wait policy: `env`, `des-default` or `runtime`.
* The run prints the same sections at start, one line each as `[group][topic]`: `[build][...]`
  lines, then `[runtime][...]`. `has_runtime_info_display = no` in the cfg silences them.

### When it is written

The record is written with the first frame, together with the first `.info` row.

* A fresh run starts the file over.
* A restart appends its record, after a `# ---- restart: another record follows ----` line.
* A run that dies before its first frame leaves the file as it was.
* A GPU build that finds no device appends its record at once, after a
  `# ---- a run that could not start follows ----` line.

### Which executable wrote it

A manifest describes the executable that *wrote* it. To see whether that executable was rebuilt
since, compare the manifest's `exe_mtime_utc` and `state_utc` with the executable you have now.
A plain `cp` or `touch` moves the mtime, so `exe_mtime_utc` dates the file, not the build.

## Each frame

Every frame and checkpoint names its own origin:

```bash
h5dump -A -g /provenance <model>.save.000000.vtkhdf    # HDF5 build
strings <model>.save.000000                            # des-binary build
```

### Where it came from

| Fields | What they say |
|---|---|
| `code_rev`, `code_branch`, `code_dirty`, `code_origin`, `code_state_utc` | the source |
| `build_os`, `builder`, `exe_mtime_utc` | the build's OS and user@host, the executable's mtime |
| `os`, `runner`, `cpu_model`, `logical_cores`, `mem_total_gib`, `kernel` | the machine it ran on |
| `omp_threads` | the team size at this write, as measured; `-1` in an `openmp=0` build |
| `restart_from` | `no`, or the `<model>:<frame>` this run restarted from |
| `gpu_model`, `gpu_device`, `gpu_visible_devices`, `gpu_cuda_driver` | the GPU, on a GPU run |
| `gpu_mem_total_gib`, `gpu_mem_free_at_start_gib` | its memory |

`gpu_mem_free_at_start_gib` is genuinely free device memory, unlike the host's MemAvailable
estimate. The static host inventory (physical, performance and efficiency cores, free memory at
start) is in the `.manifest` only.

### Health over the run

Each frame also samples the run and the machine at its write, so a run's frames form a time
series:

| Fields | What they say |
|---|---|
| `mem_rss_gib`, `mem_peak_rss_gib` | this run's memory, now and its peak so far |
| `mem_avail_gib`, `load_avg_1m` | the machine's available memory and 1-minute load |
| `cpu_time_sec` | this run's CPU time |
| `gpu_mem_used_dev_gib` | the device memory in use: the whole device, not this run |
| `write_utc` | when the frame was written |

Read the series with its caveats:

* `mem_rss_gib` and `cpu_time_sec` restart from zero on a restart leg.
* `mem_rss_gib` also falls when the kernel reclaims pages under external pressure (macOS).
* `cpu_time_sec / (walltime_sec x omp_threads)` cannot see a stall while workers spin, which is
  the macOS default, since DES sets `OMP_WAIT_POLICY=active` there.
* `mem_avail_gib` and `load_avg_1m` say whether the machine itself went bad rather than the model,
  and include this run's own load.
