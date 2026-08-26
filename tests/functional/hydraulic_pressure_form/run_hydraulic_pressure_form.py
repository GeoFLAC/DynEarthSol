#!/usr/bin/env python3

"""Check pressure-form hydraulics at zero and nonzero gravity in 2-D/3-D."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np


sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
DEFAULT_CFG = HERE / "hydraulic_pressure_form.cfg"
FLUID_DENSITY = 850.0


@dataclass(frozen=True)
class Case:
    name: str
    gravity: float
    source_enabled: bool


CASES = (
    Case("zero_gravity_source", 0.0, True),
    Case("hydrostatic_equilibrium", 10.0, False),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe-2d", type=Path, default=REPO_ROOT / "dynearthsol2d")
    parser.add_argument("--exe-3d", type=Path, required=True)
    parser.add_argument("--cfg", type=Path, default=DEFAULT_CFG)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 4])
    parser.add_argument("--run-dir", type=Path)
    return parser.parse_args()


def injection_section(ndims: int, enabled: bool) -> str:
    if not enabled:
        return "[injection]\nenabled = no\nnum_points = 0"
    points_y = "points_y = [0.47]\n" if ndims == 3 else ""
    return (
        "[injection]\n"
        "enabled = yes\n"
        "num_points = 1\n"
        "points_x = [0.43]\n"
        f"{points_y}"
        "points_z = [-0.46]\n"
        "points_unit = m\n"
        "rate_model = constant_rate\n"
        "rate = [1e-6]\n"
        "total_amount = [0]\n"
        "start_time_in_yr = [0]\n"
        "end_time_in_yr = []\n"
        "duration_in_yr = []"
    )


def render_config(template: str, case: Case, ndims: int, model: str) -> str:
    replacements = {
        "__MODELNAME__": model,
        "__GRAVITY__": f"{case.gravity:g}",
        "__PLANE_STRAIN__": "yes" if ndims == 2 else "no",
        "__SOURCE_SECTION__": injection_section(ndims, case.source_enabled),
    }
    for token, value in replacements.items():
        template = template.replace(token, value)
    if "__" in template:
        raise AssertionError("unexpanded config token remains")
    return template


def read_frames(
    case_dir: Path, model: str, case: Case, ndims: int
) -> list[np.ndarray]:
    frames: list[np.ndarray] = []
    for frame in range(3):
        path = case_dir / f"{model}.save.{frame:06d}.vtkhdf"
        if not path.is_file():
            raise AssertionError(f"missing output frame: {path}")
        with h5py.File(path, "r") as output:
            for name in ("coordinate", "pore pressure"):
                if name not in output:
                    raise AssertionError(f"{path}: missing dataset {name!r}")
            coordinate = output["coordinate"][...]
            pressure = output["pore pressure"][...]
        if coordinate.ndim != 2 or coordinate.shape[1] != ndims:
            raise AssertionError(f"{path}: expected N x {ndims} coordinates")
        if pressure.shape != (coordinate.shape[0],):
            raise AssertionError(f"{path}: pressure/coordinate node count differs")
        if not np.isfinite(pressure).all():
            raise AssertionError(f"{path}: non-finite pore pressure")
        if not np.isfinite(coordinate).all():
            raise AssertionError(f"{path}: non-finite coordinate")
        if not case.source_enabled:
            expected = -FLUID_DENSITY * case.gravity * coordinate[:, ndims - 1]
            tolerance = 2.0e-12 * max(float(np.max(np.abs(expected))), 1.0)
            if not np.allclose(pressure, expected, rtol=0.0, atol=tolerance):
                raise AssertionError(
                    f"{path}: pressure is not the expected hydrostatic profile"
                )
        frames.append(pressure)
    return frames


def check_physics(case: Case, frames: list[np.ndarray]) -> None:
    initial, first, final = frames
    if case.source_enabled:
        if not np.array_equal(initial, np.zeros_like(initial)):
            raise AssertionError("zero-gravity case did not start at zero pressure")
        if not float(first.max()) > 0.0 or not float(final.max()) > 0.0:
            raise AssertionError("zero-gravity source did not raise pore pressure")
        source_nodes = np.count_nonzero(first > 1.0e-12)
        diffused_nodes = np.count_nonzero(final > 1.0e-12)
        if not diffused_nodes > source_nodes:
            raise AssertionError(
                "zero-gravity pressure gradient did not spread beyond its first-step support"
            )
        return

    expected_scale = max(float(np.max(np.abs(initial))), 1.0)
    tolerance = 2.0e-12 * expected_scale
    if not np.allclose(first, initial, rtol=0.0, atol=tolerance):
        raise AssertionError("hydrostatic pressure changed after one step")
    if not np.allclose(final, initial, rtol=0.0, atol=tolerance):
        raise AssertionError("hydrostatic pressure changed after two steps")


def run_case(
    exe: Path,
    template: str,
    run_root: Path,
    case: Case,
    ndims: int,
    threads: int,
) -> list[np.ndarray]:
    case_dir = run_root / f"{case.name}-{ndims}d-omp{threads}"
    case_dir.mkdir(parents=True, exist_ok=False)
    model = f"pressure_form_{case.name}_{ndims}d_omp{threads}"
    cfg_path = case_dir / "input.cfg"
    cfg_path.write_text(render_config(template, case, ndims, model), encoding="ascii")

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
        timeout=90,
        check=False,
    )
    (case_dir / "run.log").write_text(result.stdout, encoding="utf-8")
    if result.returncode != 0:
        print(result.stdout, file=sys.stderr)
        raise RuntimeError(
            f"{case.name} {ndims}D OMP={threads}: exit {result.returncode}; "
            f"see {case_dir / 'run.log'}"
        )

    frames = read_frames(case_dir, model, case, ndims)
    check_physics(case, frames)
    print(f"pressure form {case.name} {ndims}D OMP={threads}: PASS", flush=True)
    return frames


def run_all(args: argparse.Namespace, template: str, run_root: Path) -> None:
    executables = (
        (2, args.exe_2d.expanduser().resolve()),
        (3, args.exe_3d.expanduser().resolve()),
    )
    threads = list(dict.fromkeys(args.threads))
    if len(threads) != len(args.threads) or any(value <= 0 for value in threads):
        raise ValueError("--threads values must be unique positive integers")

    for ndims, exe in executables:
        if not exe.is_file():
            raise FileNotFoundError(exe)
        for case in CASES:
            reference: list[np.ndarray] | None = None
            for thread_count in threads:
                frames = run_case(
                    exe, template, run_root, case, ndims, thread_count
                )
                if reference is None:
                    reference = frames
                    continue
                for frame, (actual, expected) in enumerate(zip(frames, reference)):
                    if not np.allclose(actual, expected, rtol=1.0e-13, atol=1.0e-12):
                        raise AssertionError(
                            f"{case.name} {ndims}D frame {frame}: "
                            "OpenMP results differ"
                        )


def main() -> None:
    args = parse_args()
    template = args.cfg.expanduser().resolve().read_text(encoding="ascii")
    if args.run_dir is not None:
        run_root = args.run_dir.expanduser().resolve()
        run_root.mkdir(parents=True, exist_ok=True)
        run_all(args, template, run_root)
    else:
        with tempfile.TemporaryDirectory(prefix="des-hydraulic-pressure-form-") as tmp:
            run_all(args, template, Path(tmp))


if __name__ == "__main__":
    main()
