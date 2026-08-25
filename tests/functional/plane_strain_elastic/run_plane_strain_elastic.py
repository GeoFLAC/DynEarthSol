#!/usr/bin/env python3

"""Check the accumulated out-of-plane stress of 2-D elastic plane strain."""

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
DEFAULT_CFG = HERE / "plane_strain_elastic.cfg"
STRESSYY_DATASET = "stressyy"
STRAIN_DATASET = "strain"
BULK_MODULUS = 2.0e8
SHEAR_MODULUS = 2.0e8
LAME_LAMBDA = BULK_MODULUS - 2.0 * SHEAR_MODULUS / 3.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the feature-off 2-D plane-strain elastic stress regression. "
            "The executable must be an HDF5-enabled 2-D build."
        )
    )
    parser.add_argument("--exe", type=Path, default=DEFAULT_EXE)
    parser.add_argument("--cfg", type=Path, default=DEFAULT_CFG)
    parser.add_argument(
        "--run-dir",
        type=Path,
        help="Keep outputs below this directory (default: temporary directory).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        nargs="+",
        default=[1, 4],
        help="OMP thread counts to exercise (default: 1 4).",
    )
    return parser.parse_args()


def read_array(path: Path, dataset: str):
    if not path.is_file():
        raise AssertionError(f"missing output file: {path}")
    with h5py.File(path, "r") as output:
        if dataset not in output:
            raise AssertionError(f"{path}: missing dataset {dataset!r}")
        return output[dataset][...]


def make_case_dir(run_root: Path, threads: int) -> Path:
    case_dir = run_root / f"omp{threads}"
    try:
        case_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError as exc:
        raise FileExistsError(
            f"refusing to overwrite existing regression directory: {case_dir}"
        ) from exc
    return case_dir


def run_model(
    exe: Path, template: str, run_root: Path, threads: int
) -> tuple[Path, str]:
    case_dir = make_case_dir(run_root, threads)
    modelname = f"plane_strain_elastic_omp{threads}"
    if "__MODELNAME__" not in template:
        raise AssertionError("missing config token __MODELNAME__")
    cfg = template.replace("__MODELNAME__", modelname)
    if "__" in cfg:
        raise AssertionError("unexpanded config token remains")

    cfg_path = case_dir / "input.cfg"
    cfg_path.write_text(cfg, encoding="ascii")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    result = subprocess.run(
        [str(exe), str(cfg_path)],
        cwd=case_dir,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    (case_dir / "run.log").write_text(result.stdout, encoding="utf-8")
    if result.returncode != 0:
        print(result.stdout, file=sys.stderr)
        raise RuntimeError(
            f"plane-strain elastic OMP={threads}: "
            f"DynEarthSol exited with {result.returncode}"
        )
    return case_dir, modelname


def check_accumulated_increment(case_dir: Path, modelname: str) -> None:
    initial_strain = read_array(
        case_dir / f"{modelname}.save.000000.vtkhdf", STRAIN_DATASET
    )
    final_strain = read_array(
        case_dir / f"{modelname}.save.000001.vtkhdf", STRAIN_DATASET
    )
    initial_stressyy = read_array(
        case_dir / f"{modelname}.chkpt.000000.vtkhdf", STRESSYY_DATASET
    )
    final_stressyy = read_array(
        case_dir / f"{modelname}.chkpt.000001.vtkhdf", STRESSYY_DATASET
    )

    if initial_strain.shape != final_strain.shape:
        raise AssertionError("strain array shape changed across the two-step run")
    if initial_stressyy.shape != final_stressyy.shape:
        raise AssertionError("stressyy array shape changed across the two-step run")
    if len(final_strain) != len(final_stressyy):
        raise AssertionError("strain/stressyy element count differs")

    max_trace_increment = 0.0
    for element in range(len(final_stressyy)):
        trace_increment = (
            float(final_strain[element][0] - initial_strain[element][0])
            + float(final_strain[element][1] - initial_strain[element][1])
        )
        stress_increment = float(
            final_stressyy[element] - initial_stressyy[element]
        )
        expected = LAME_LAMBDA * trace_increment
        max_trace_increment = max(max_trace_increment, abs(trace_increment))
        if not math.isclose(
            expected, stress_increment, rel_tol=1.0e-12, abs_tol=1.0e-10
        ):
            raise AssertionError(
                "plane-strain elastic stressyy increment mismatch at element "
                f"{element}: expected={expected:.17e}, "
                f"actual={stress_increment:.17e}"
            )

    if not max_trace_increment > 0.0:
        raise AssertionError("fixture produced no in-plane volumetric strain")


def run_all(exe: Path, template: str, run_root: Path, threads: list[int]) -> None:
    seen: set[int] = set()
    for thread_count in threads:
        if thread_count <= 0:
            raise ValueError("--threads values must be positive")
        if thread_count in seen:
            raise ValueError(f"duplicate --threads value: {thread_count}")
        seen.add(thread_count)
        case_dir, modelname = run_model(exe, template, run_root, thread_count)
        check_accumulated_increment(case_dir, modelname)
        print(f"plane-strain elastic regression OMP={thread_count}: PASS")


def main() -> None:
    args = parse_args()
    exe = args.exe.expanduser().resolve()
    cfg_path = args.cfg.expanduser().resolve()
    if not exe.is_file():
        raise FileNotFoundError(f"executable not found: {exe}")
    if not cfg_path.is_file():
        raise FileNotFoundError(f"config not found: {cfg_path}")
    template = cfg_path.read_text(encoding="ascii")

    if args.run_dir:
        run_root = args.run_dir.expanduser().resolve()
        run_root.mkdir(parents=True, exist_ok=True)
        run_all(exe, template, run_root, args.threads)
    else:
        with tempfile.TemporaryDirectory(prefix="des-plane-strain-elastic-") as tmp:
            run_all(exe, template, Path(tmp), args.threads)


if __name__ == "__main__":
    main()
