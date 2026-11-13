# Try CPU MPI and multiple GPUs

These are small functional checks for reviewers. On managed clusters, run parallel
commands inside an allocation with enough slots for the requested ranks/threads.
Use mounted storage for `REVIEW_WORK`. On offline clusters, populate dependencies
and submodules on a connected machine before staging the checkout.

## CPU MPI

Prerequisites: C++17 compiler, MPI C++ wrapper, HDF5, Python with NumPy and h5py,
and the repository submodules. MMG also needs CMake and a C compiler.

Native macOS CPU builds use Apple clang and an MPI wrapper backed by clang.
Homebrew can provide the dependencies with
`brew install cmake libomp open-mpi hdf5`; install NumPy and h5py in your Python
environment. The Makefile detects Homebrew's OpenMP and HDF5 paths. The CPU
commands below were checked on Apple Silicon; see the
[verification scope](mpi.md#verification-scope).

Reserve at least six CPU slots for this driver: its empty-partition test uses
six ranks and its hybrid check uses two ranks with two threads each.

From the repository root:

```sh
export DES_ROOT="$PWD"
export REVIEW_WORK=/path/to/mounted/des-review
mkdir -p "$REVIEW_WORK"
git submodule update --init --recursive
for dim in 2 3; do
    make -j4 ndims="$dim" usemmg=1 hdf5=1
    make -j4 ndims="$dim" usempi=1 usemmg=1 hdf5=1 MPI_CXX=mpicxx
done
OMP_NUM_THREADS=1 python3 tests/functional/mpi/run_mpi.py \
    --serial ./dynearthsol2d --mpi ./dynearthsol2d.mpi \
    --serial3d ./dynearthsol3d --mpi3d ./dynearthsol3d.mpi \
    --launcher 'mpiexec --bind-to none' --ranks 1 2 \
    --workdir "$REVIEW_WORK/cpu"
```

This compares serial and MPI physical fields, checkpoints and restarts, including
bounded remeshing cases, and checks that root-side failures terminate the job.
Expected: exit status 0 and `Passed. Results: .../results.json`. The generated
directory retains cfg files, commands, logs and comparisons. These are correctness
checks; do not infer speedup from their runtimes. For parent regression, supply
executables built from `master` as the serial references.

For a direct run of an existing cfg, use
`OMP_NUM_THREADS=1 mpiexec --bind-to none -n 2 /path/to/dynearthsol2d.mpi model.cfg`
from its model directory. No new cfg parameter enables MPI: select the MPI
executable and launcher. `OMP_NUM_THREADS` controls threads per rank; reserve
their combined CPU count.

## One or two NVIDIA GPUs

Use NVHPC and its compatible MPI wrapper/runtime; set the example SDK paths and
`GPU_CC` for your installation/device. Keep the CPU binaries from the preceding
section. The project documents NVHPC 25.7 as its supported minimum.

```sh
export NVHPC_DIR=/path/to/nvidia/hpc_sdk/Linux_x86_64/25.7
export GPU_MPI="$NVHPC_DIR/comm_libs/mpi/bin/mpicxx"
for dim in 2 3; do
    OMPI_CXX="$NVHPC_DIR/compilers/bin/nvc++" \
    make -j4 ndims="$dim" usempi=1 openacc=1 openmp=0 nofma=1 \
        GPU_CC=80 NVHPC_DIR="$NVHPC_DIR" MPI_CXX="$GPU_MPI" \
        usemmg=1 hdf5=1
done
OMP_NUM_THREADS=1 python3 tests/functional/mpi/run_gpu_mpi.py \
    --cpu ./dynearthsol2d --gpu ./dynearthsol2d.mpi.gpu \
    --cpu3d ./dynearthsol3d --gpu3d ./dynearthsol3d.mpi.gpu \
    --launcher "$NVHPC_DIR/comm_libs/mpi/bin/mpiexec --bind-to none" \
    --ranks 1 2 --cases fixed --steps 8 --output-interval 2 \
    --workdir "$REVIEW_WORK/gpu"
```

Set the matching SDK/MPI runtime library paths as required by your installation.
Allocate two visible GPUs for this check. One rank selects one physical GPU;
the report checks distinct GPU UUIDs, actual device execution and all saved
fields/checkpoints against the CPU and one-GPU references. Expected: exit 0,
`all_passed: true` in the generated `results.json`, and distinct UUIDs for the
two-rank cases. For one available GPU use `--ranks 1`. For a direct model run,
use `mpiexec -n 2 /path/to/dynearthsol3d.mpi.gpu model.cfg`.

Halos use pinned host staging. Initialization, MMG, remapping and global output
remain on the root CPU. GPU count does not imply speedup; measure the complete
model's wall time. See [architecture and limitations](mpi.md).
