# Provenance

Every DES executable records where it came from, in a `build.snapshot` block.

## The build

The `build.snapshot` block holds the revision, make options, dependency versions and providers,
toolchain, and the compile and link flags. Read it without running the executable:

```bash
strings <exe> | grep '^build\.snapshot\.'
```

A run prints the same block at start as `[build][...]` lines, ahead of the `[Runtime]` host and
device lines; `has_runtime_info_display = no` in the cfg silences both.

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
