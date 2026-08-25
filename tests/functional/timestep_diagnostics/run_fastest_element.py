#!/usr/bin/env python3

"""Check the 3-D fastest-element timestep diagnostic with pure y velocity."""

from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
DEFAULT_EXE = REPO_ROOT / "dynearthsol3d"
DEFAULT_CFG = HERE / "fastest_element_3d.cfg"
EXPECTED_VY = 1.0e-6
REL_TOL = 1.0e-12
ABS_TOL = 1.0e-15
DIAGNOSTIC = re.compile(
    r"max_global_vel_mag=(?P<global_speed>[^ ]+)"
    r".*?\smax_vel_elem=(?P<element>-?\d+)"
    r"\s+max_vel_coord=\((?P<coord>[^)]*)\)"
    r"\s+max_vel=\((?P<velocity>[^)]*)\)"
)


@dataclass(frozen=True)
class FastestElementDiagnostic:
    element: int
    coord: tuple[float, float, float]
    velocity: tuple[float, float, float]
    global_speed: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the 3-D pure-y fastest-element diagnostic regression."
    )
    parser.add_argument("--exe", type=Path, default=DEFAULT_EXE)
    parser.add_argument("--cfg", type=Path, default=DEFAULT_CFG)
    parser.add_argument(
        "--threads",
        type=int,
        nargs="+",
        default=(1, 4),
        help="OpenMP thread counts to test (default: 1 4).",
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        help="Keep outputs below this directory (default: temporary directory).",
    )
    return parser.parse_args()


def parse_triple(text: str, label: str) -> tuple[float, float, float]:
    fields = text.split(",")
    if len(fields) != 3:
        raise AssertionError(
            f"{label}: expected 3 components, found {len(fields)} in ({text})"
        )
    try:
        values = tuple(float(field) for field in fields)
    except ValueError as exc:
        raise AssertionError(f"{label}: non-numeric tuple ({text})") from exc
    if not all(math.isfinite(value) for value in values):
        raise AssertionError(f"{label}: non-finite tuple {values}")
    return values[0], values[1], values[2]


def parse_first_diagnostic(stdout: str, label: str) -> FastestElementDiagnostic:
    dt_lines = [line for line in stdout.splitlines() if "compute_dt:" in line]
    if not dt_lines:
        raise AssertionError(f"{label}: no compute_dt diagnostic found")

    match = DIAGNOSTIC.search(dt_lines[0])
    if not match:
        raise AssertionError(
            f"{label}: first compute_dt record lacks a complete fastest-element "
            f"diagnostic:\n{dt_lines[0]}"
        )

    try:
        element = int(match.group("element"))
        global_speed = float(match.group("global_speed"))
    except ValueError as exc:
        raise AssertionError(f"{label}: invalid scalar diagnostic value") from exc
    if element < 0:
        raise AssertionError(f"{label}: invalid fastest element {element}")
    if not math.isfinite(global_speed):
        raise AssertionError(f"{label}: non-finite global speed {global_speed}")

    return FastestElementDiagnostic(
        element=element,
        coord=parse_triple(match.group("coord"), f"{label}: max_vel_coord"),
        velocity=parse_triple(match.group("velocity"), f"{label}: max_vel"),
        global_speed=global_speed,
    )


def assert_close(actual: float, expected: float, label: str) -> None:
    if not math.isclose(actual, expected, rel_tol=REL_TOL, abs_tol=ABS_TOL):
        raise AssertionError(
            f"{label}: actual={actual:.17e}, expected={expected:.17e}"
        )


def check_diagnostic(record: FastestElementDiagnostic, label: str) -> None:
    vx, vy, vz = record.velocity
    assert_close(vx, 0.0, f"{label}: vx")
    assert_close(vy, EXPECTED_VY, f"{label}: vy")
    assert_close(vz, 0.0, f"{label}: vz")

    tuple_speed = math.sqrt(sum(component * component for component in record.velocity))
    assert_close(
        tuple_speed,
        record.global_speed,
        f"{label}: velocity tuple norm versus max_global_vel_mag",
    )
    assert_close(record.global_speed, EXPECTED_VY, f"{label}: global speed")


def run_case(exe: Path, cfg: Path, run_root: Path, threads: int) -> FastestElementDiagnostic:
    label = f"OMP_NUM_THREADS={threads}"
    case_dir = run_root / f"omp-{threads}"
    case_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["PYTHONDONTWRITEBYTECODE"] = "1"

    result = subprocess.run(
        [str(exe), str(cfg)],
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
        raise RuntimeError(f"{label}: DynEarthSol exited with {result.returncode}")

    record = parse_first_diagnostic(result.stdout, label)
    check_diagnostic(record, label)
    print(
        f"PASS {label}: element={record.element}, coord={record.coord}, "
        f"velocity={record.velocity}"
    )
    return record


def check_thread_invariance(
    baseline: FastestElementDiagnostic,
    candidate: FastestElementDiagnostic,
    threads: int,
) -> None:
    label = f"OMP_NUM_THREADS={threads}: thread invariance"
    if candidate.element != baseline.element:
        raise AssertionError(
            f"{label}: element={candidate.element}, baseline={baseline.element}"
        )
    for name, actual, expected in (
        ("x", candidate.coord[0], baseline.coord[0]),
        ("y", candidate.coord[1], baseline.coord[1]),
        ("z", candidate.coord[2], baseline.coord[2]),
        ("vx", candidate.velocity[0], baseline.velocity[0]),
        ("vy", candidate.velocity[1], baseline.velocity[1]),
        ("vz", candidate.velocity[2], baseline.velocity[2]),
        ("global speed", candidate.global_speed, baseline.global_speed),
    ):
        assert_close(actual, expected, f"{label}: {name}")


def run_all(exe: Path, cfg: Path, run_root: Path, threads: list[int]) -> None:
    baseline: FastestElementDiagnostic | None = None
    for thread_count in threads:
        record = run_case(exe, cfg, run_root, thread_count)
        if baseline is None:
            baseline = record
        else:
            check_thread_invariance(baseline, record, thread_count)


def main() -> int:
    args = parse_args()
    exe = args.exe.resolve()
    cfg = args.cfg.resolve()
    threads = list(dict.fromkeys(args.threads))

    if not exe.is_file():
        raise FileNotFoundError(f"DynEarthSol executable not found: {exe}")
    if not cfg.is_file():
        raise FileNotFoundError(f"regression config not found: {cfg}")
    if not threads or any(value < 1 for value in threads):
        raise ValueError("--threads values must be positive integers")

    if args.run_dir is not None:
        run_root = args.run_dir.resolve()
        run_root.mkdir(parents=True, exist_ok=True)
        run_all(exe, cfg, run_root, threads)
    else:
        with tempfile.TemporaryDirectory(prefix="des-fastest-element-") as tmp:
            run_all(exe, cfg, Path(tmp), threads)

    print("PASS 3-D fastest-element timestep diagnostic regression")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
