# Unsaturated header backend check

`unsaturated.hpp` is header-only, so this check builds it alone; nothing is
linked from the production tree. Each routine runs serially on the host and
again inside a parallel loop (OpenACC on the device, or OpenMP), and the
two results are compared.

Run from this directory:

```sh
make check backend=serial
OMP_NUM_THREADS=4 make check backend=omp
make check backend=asan
make check backend=acc              # NVHPC and a GPU
make check backend=acc nofma=1      # the production nofma=1 variant
```

`GPU_CC` defaults to the detected device, with a fallback of 80, and can be
overridden. The accelerator flags are the production `openacc=1` ones.

Checks:

| Check | Pass criterion |
| --- | --- |
| Se, dSe/dp, S, kr and Bishop pressure for 5 laws over 11 decades of suction | host/parallel relative difference <= 1e-12, all finite |
| Saturated patches, 2e5 nodes | parallel result bit-identical to `p + dm/hmass`, the existing update |
| Unsaturated single-law and mixed-material patches (one law on the saturated model), 2e5 nodes each | content residual <= 1e-15, no invalid flag, pressure host/parallel difference <= 1e-9 relative |
| NaN pressure and infinite increment | `invalid` set, or a non-finite result reported |

Also reported: throughput in ns/node for each node class. On the
accelerator this includes host-device transfers. Its purpose is to show the
cost spread between saturated nodes (the fast path) and unsaturated nodes
(about 3 Newton steps each).

On a GPU, read the build log first. `-Minfo=accel` must list both parallel
loops as generated for the accelerator, and must show no "routine ...
not found" or host fallback for `advance_content`, `pore_content` or the
`retention_*` routines. If nvc++ rejects `#pragma acc routine seq` placed
before `template <class Patch>`, record the message. That placement is the
one open portability question in the header.

Reference results on the development CPU (g++ 13.3, WSL2):

| Backend | saturated | unsaturated | mixed |
| --- | --- | --- | --- |
| serial | 38 ns/node | 189 ns/node | 377 ns/node |
| OpenMP, 4 threads | 2.8 ns/node | 52 ns/node | 107 ns/node |

All checks pass on serial, OpenMP and ASan/UBSan. OpenACC is not yet run.
