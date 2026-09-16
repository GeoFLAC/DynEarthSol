# Provenance

DES records where each build and run came from:

| What | Where |
|---|---|
| The build | in the executable, as a `build.snapshot` block |
| The run | in `<modelname>.manifest`, beside the model's output |

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
