#!/usr/bin/env python3
"""Check zero-gravity homogeneous absolute initial stress at t=0."""

from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py


sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
DEFAULT_EXE = REPO_ROOT / "dynearthsol2d"
DEFAULT_CFG = HERE / "initial_stress.cfg"
PRESCRIBED_STRESS = (-5.0e6, -3.0e6, 9.0e5)
PRESCRIBED_STRESSYY = -4.0e6


def read_array(path: Path, dataset: str):
    if not path.is_file():
        raise AssertionError(f"missing output file: {path}")
    with h5py.File(path, "r") as output:
        if dataset not in output:
            raise AssertionError(f"{path}: missing dataset {dataset!r}")
        return output[dataset][...]


def render_cfg(
    template: str,
    modelname: str,
    option: int,
    gravity: float,
) -> str:
    replacements = {
        "__MODELNAME__": modelname,
        "__GRAVITY__": f"{gravity:.17g}",
        "__INITIAL_STRESS_OPTION__": str(option),
        "__INITIAL_STRESS__": "["
        + ",".join(f"{value:.17g}" for value in PRESCRIBED_STRESS)
        + "]",
    }
    cfg = template
    for token, value in replacements.items():
        if cfg.count(token) != 1:
            raise AssertionError(f"expected one config token {token}")
        cfg = cfg.replace(token, value)
    if "__" in cfg:
        raise AssertionError("unexpanded config token remains")
    return cfg


def run_model(
    exe: Path,
    template: str,
    root: Path,
    label: str,
    *,
    option: int,
    gravity: float,
    threads: int,
    expect_success: bool = True,
) -> tuple[Path, str, subprocess.CompletedProcess[str]]:
    modelname = f"initial_stress_{label}_omp{threads}"
    run_dir = root / modelname
    run_dir.mkdir()
    cfg_path = run_dir / "input.cfg"
    cfg_path.write_text(
        render_cfg(template, modelname, option, gravity),
        encoding="ascii",
    )
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    result = subprocess.run(
        [str(exe), str(cfg_path)],
        cwd=run_dir,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if expect_success and result.returncode != 0:
        raise RuntimeError(
            f"{label} OMP={threads} failed with {result.returncode}:\n"
            f"{result.stdout}"
        )
    if not expect_success and result.returncode == 0:
        raise AssertionError(f"{label} OMP={threads} unexpectedly succeeded")
    return run_dir, modelname, result


def assert_close(label: str, actual: float, expected: float) -> None:
    if not math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-8):
        raise AssertionError(
            f"{label}: {actual:.17e} != {expected:.17e}"
        )


def check_initial_state(
    run_dir: Path,
    modelname: str,
    expected_stress: tuple[float, float, float],
    expected_stressyy: float,
) -> None:
    stress = read_array(
        run_dir / f"{modelname}.save.000000.vtkhdf", "stress"
    )
    strain = read_array(
        run_dir / f"{modelname}.save.000000.vtkhdf", "strain"
    )
    stressyy = read_array(
        run_dir / f"{modelname}.chkpt.000000.vtkhdf", "stressyy"
    )
    if len(stress) != len(strain) or len(stress) != len(stressyy):
        raise AssertionError("initial stress/strain/stressyy counts differ")
    for element, (stress_e, strain_e, stressyy_e) in enumerate(
        zip(stress, strain, stressyy)
    ):
        for component, expected in enumerate(expected_stress):
            assert_close(
                f"element {element} stress[{component}]",
                float(stress_e[component]),
                expected,
            )
            assert_close(
                f"element {element} strain[{component}]",
                float(strain_e[component]),
                0.0,
            )
        assert_close(
            f"element {element} stressyy",
            float(stressyy_e),
            expected_stressyy,
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", type=Path, default=DEFAULT_EXE)
    parser.add_argument("--cfg", type=Path, default=DEFAULT_CFG)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 4])
    args = parser.parse_args()
    exe = args.exe.expanduser().resolve()
    cfg_path = args.cfg.expanduser().resolve()
    if not exe.is_file() or not os.access(exe, os.X_OK):
        raise FileNotFoundError(f"executable not found: {exe}")
    template = cfg_path.read_text(encoding="ascii")

    with tempfile.TemporaryDirectory(prefix="des-initial-stress-") as temp:
        root = Path(temp)
        for threads in args.threads:
            if threads <= 0:
                raise ValueError("thread counts must be positive")
            baseline_dir, baseline_name, _ = run_model(
                exe,
                template,
                root,
                "legacy",
                option=0,
                gravity=0.0,
                threads=threads,
            )
            check_initial_state(
                baseline_dir,
                baseline_name,
                (0.0, 0.0, 0.0),
                0.0,
            )
            imposed_dir, imposed_name, _ = run_model(
                exe,
                template,
                root,
                "absolute",
                option=1,
                gravity=0.0,
                threads=threads,
            )
            check_initial_state(
                imposed_dir,
                imposed_name,
                PRESCRIBED_STRESS,
                PRESCRIBED_STRESSYY,
            )
            print(f"initial stress t=0 regression OMP={threads}: PASS")

        _, _, invalid = run_model(
            exe,
            template,
            root,
            "gravity_rejected",
            option=1,
            gravity=9.81,
            threads=1,
            expect_success=False,
        )
        if "requires control.gravity=0" not in invalid.stdout:
            raise AssertionError(
                "missing gravity incompatibility diagnostic:\n"
                + invalid.stdout
            )
        print("initial stress gravity compatibility check: PASS")


if __name__ == "__main__":
    main()
