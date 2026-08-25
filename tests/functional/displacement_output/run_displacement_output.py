#!/usr/bin/env python3

"""Verify vector displacement output in 2-D and, optionally, 3-D."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py


sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
DEFAULT_CFG = HERE / "displacement_output.cfg"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe-2d", type=Path, default=REPO_ROOT / "dynearthsol2d")
    parser.add_argument("--exe-3d", type=Path)
    parser.add_argument("--cfg", type=Path, default=DEFAULT_CFG)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 4])
    parser.add_argument("--expected-format", choices=("hdf5", "binary"))
    parser.add_argument("--run-dir", type=Path)
    return parser.parse_args()


def read_frame(path: Path):
    if not path.is_file():
        raise AssertionError(f"missing output frame: {path}")
    with h5py.File(path, "r") as output:
        names = set(output.keys())
        for required in ("coordinate", "coord0", "displacement"):
            if required not in names:
                raise AssertionError(f"{path}: missing dataset {required!r}")
        if "vertical displacement" in names:
            raise AssertionError(f"{path}: redundant vertical displacement dataset returned")
        return (
            output["coordinate"][...],
            output["coord0"][...],
            output["displacement"][...],
        )


def read_binary_frame(model: Path, frame: int):
    sys.path.insert(0, str(REPO_ROOT))
    try:
        from Dynearthsol import Dynearthsol
    finally:
        sys.path.pop(0)

    output = Dynearthsol(str(model))
    output.read_header(frame)
    names = set(output.field_pos)
    for required in ("coordinate", "coord0", "displacement"):
        if required not in names:
            raise AssertionError(f"{model}: missing binary field {required!r}")
    if "vertical displacement" in names:
        raise AssertionError(f"{model}: redundant vertical displacement field returned")
    return (
        output.read_field(frame, "coordinate"),
        output.read_field(frame, "coord0"),
        output.read_field(frame, "displacement"),
    )


def check_frame(model: Path, frame: int, ndims: int, expect_motion: bool) -> None:
    hdf5_path = Path(f"{model}.save.{frame:06d}.vtkhdf")
    binary_path = Path(f"{model}.save.{frame:06d}")
    if hdf5_path.is_file():
        coordinate, reference, displacement = read_frame(hdf5_path)
        path = hdf5_path
    elif binary_path.is_file():
        coordinate, reference, displacement = read_binary_frame(model, frame)
        path = binary_path
    else:
        raise AssertionError(f"missing output frame for {model}, frame {frame}")
    if coordinate.shape != reference.shape or coordinate.shape != displacement.shape:
        raise AssertionError(
            f"{path}: coordinate/reference/displacement shapes differ: "
            f"{coordinate.shape}, {reference.shape}, {displacement.shape}"
        )
    if len(coordinate.shape) != 2 or coordinate.shape[1] != ndims:
        raise AssertionError(f"{path}: expected an N x {ndims} nodal vector")

    max_error = 0.0
    max_displacement = 0.0
    for node in range(coordinate.shape[0]):
        for component in range(ndims):
            expected = float(coordinate[node][component] - reference[node][component])
            actual = float(displacement[node][component])
            max_error = max(max_error, abs(actual - expected))
            max_displacement = max(max_displacement, abs(actual))
    if max_error != 0.0:
        raise AssertionError(f"{path}: displacement != coordinate - coord0; max error={max_error}")
    if expect_motion and not max_displacement > 0.0:
        raise AssertionError(f"{path}: fixture produced no nodal displacement")
    if not expect_motion and max_displacement != 0.0:
        raise AssertionError(f"{path}: initial displacement is not zero")


def run_case(
    exe: Path,
    template: str,
    run_root: Path,
    ndims: int,
    threads: int,
    expected_format: str | None,
) -> None:
    case_dir = run_root / f"{ndims}d-omp{threads}"
    case_dir.mkdir(parents=True, exist_ok=False)
    model = f"displacement_{ndims}d_omp{threads}"
    cfg = template.replace("__MODELNAME__", model)
    if "__" in cfg:
        raise AssertionError("unexpanded config token remains")
    cfg_path = case_dir / "input.cfg"
    cfg_path.write_text(cfg, encoding="ascii")

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
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
            f"displacement output {ndims}D OMP={threads}: exit {result.returncode}"
        )

    model_path = case_dir / model
    hdf5_output = Path(f"{model_path}.save.000000.vtkhdf").is_file()
    actual_format = "hdf5" if hdf5_output else "binary"
    if expected_format is not None and actual_format != expected_format:
        raise AssertionError(
            f"expected {expected_format} output, found {actual_format} for {model_path}"
        )

    check_frame(model_path, 0, ndims, False)
    check_frame(model_path, 1, ndims, True)
    print(f"displacement output {ndims}D {actual_format} OMP={threads}: PASS")


def run_all(args: argparse.Namespace, template: str, run_root: Path) -> None:
    executables = [(2, args.exe_2d.expanduser().resolve())]
    if args.exe_3d is not None:
        executables.append((3, args.exe_3d.expanduser().resolve()))
    for ndims, exe in executables:
        if not exe.is_file():
            raise FileNotFoundError(exe)
        seen: set[int] = set()
        for threads in args.threads:
            if threads <= 0 or threads in seen:
                raise ValueError("--threads values must be unique positive integers")
            seen.add(threads)
            run_case(exe, template, run_root, ndims, threads, args.expected_format)


def main() -> None:
    args = parse_args()
    cfg = args.cfg.expanduser().resolve()
    template = cfg.read_text(encoding="ascii")
    if args.run_dir is not None:
        run_root = args.run_dir.expanduser().resolve()
        run_root.mkdir(parents=True, exist_ok=True)
        run_all(args, template, run_root)
    else:
        with tempfile.TemporaryDirectory(prefix="des-displacement-output-") as tmp:
            run_all(args, template, Path(tmp))


if __name__ == "__main__":
    main()
