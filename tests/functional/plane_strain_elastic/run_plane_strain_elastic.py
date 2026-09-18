#!/usr/bin/env python3

"""Check the elastic out-of-plane stress update in 2-D plane strain."""

import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
EXE = Path(sys.argv[1] if len(sys.argv) > 1 else REPO_ROOT / "dynearthsol2d").resolve()
CFG = HERE / "plane_strain_elastic.cfg"
MODEL = "plane-strain-elastic"
LAME_LAMBDA = 2.0e8 - 2.0 * 2.0e8 / 3.0


def read(path: Path, dataset: str):
    with h5py.File(path, "r") as output:
        return output[dataset][...]


if not EXE.is_file():
    raise FileNotFoundError(f"executable not found: {EXE}")

with tempfile.TemporaryDirectory(prefix="des-plane-strain-elastic-") as tmp:
    run_dir = Path(tmp)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    subprocess.run([str(EXE), str(CFG)], cwd=run_dir, env=env, check=True)

    strain0 = read(run_dir / f"{MODEL}.save.000000.vtkhdf", "strain")
    strain1 = read(run_dir / f"{MODEL}.save.000001.vtkhdf", "strain")
    stress0 = read(run_dir / f"{MODEL}.chkpt.000000.vtkhdf", "stressyy")
    stress1 = read(run_dir / f"{MODEL}.chkpt.000001.vtkhdf", "stressyy")

    if strain0.shape != strain1.shape or len(strain1) != len(stress1):
        raise AssertionError("plane-strain output dimensions differ")

    max_trace_increment = 0.0
    for element, (old_strain, new_strain, old_stress, new_stress) in enumerate(
        zip(strain0, strain1, stress0, stress1)
    ):
        trace_increment = float(
            new_strain[0] - old_strain[0] + new_strain[1] - old_strain[1]
        )
        expected = LAME_LAMBDA * trace_increment
        actual = float(new_stress - old_stress)
        max_trace_increment = max(max_trace_increment, abs(trace_increment))
        if not math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-10):
            raise AssertionError(
                f"element {element}: expected stressyy increment {expected:.17e}, "
                f"got {actual:.17e}"
            )

    if max_trace_increment == 0.0:
        raise AssertionError("test produced no in-plane volumetric strain")

print("elastic plane-strain stress update: PASS")
