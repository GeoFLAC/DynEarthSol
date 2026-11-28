#!/usr/bin/env python3

"""Verify vector displacement output in 2-D and, optionally, 3-D."""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

import h5py


sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
DEFAULT_CFG = HERE / "displacement_output.cfg"
sys.path.insert(0, str(REPO_ROOT / "tests/functional"))
from des_test_utils import checkpoint_file, run_des


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe-2d", type=Path, default=REPO_ROOT / "dynearthsol2d")
    parser.add_argument("--exe-3d", type=Path)
    parser.add_argument("--cfg", type=Path, default=DEFAULT_CFG)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 4])
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
        association = []
        output.visititems(lambda name, obj: association.append(name)
                          if name.endswith('/PointData/displacement') else None)
        if len(association) != 1:
            raise AssertionError(f"{path}: displacement is not one PointData vector")
        vector = output[association[0]][...]
        if vector.shape != (output['coordinate'].shape[0], 3):
            raise AssertionError(f"{path}: VTK displacement must have three components")
        if output['coordinate'].shape[1] == 2 and vector[:, 2].any():
            raise AssertionError(f"{path}: 2D displacement padding must be zero")
        return (
            output["coordinate"][...],
            output["coord0"][...],
            output["displacement"][...],
        )


def check_frame(model: Path, frame: int, ndims: int, expect_motion: bool) -> None:
    path = checkpoint_file(f"{model}.save.{frame:06d}")
    if path is None:
        raise AssertionError(f"missing output frame for {model}, frame {frame}")
    coordinate, reference, displacement = read_frame(path)
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
) -> None:
    case_dir = run_root / f"{ndims}d-omp{threads}"
    case_dir.mkdir(parents=True, exist_ok=False)
    model = f"displacement_{ndims}d_omp{threads}"
    cfg = template.replace("__MODELNAME__", model)
    if "__" in cfg:
        raise AssertionError("unexpanded config token remains")
    result = run_des(exe, case_dir, "input", cfg, None, threads=threads, check=False)
    if result.returncode != 0:
        print(result.stdout, file=sys.stderr)
        raise RuntimeError(
            f"displacement output {ndims}D OMP={threads}: exit {result.returncode}"
        )

    model_path = case_dir / model
    check_frame(model_path, 0, ndims, False)
    check_frame(model_path, 1, ndims, True)
    print(f"displacement output {ndims}D OMP={threads}: PASS")


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
            run_case(exe, template, run_root, ndims, threads)


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
