#!/usr/bin/env python3

"""Check that remeshing rebuilds geometry-derived boundary projection data."""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py


sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
DEFAULT_EXE = REPO_ROOT / "dynearthsol2d"
DEFAULT_CFG = HERE / "boundary_remesh.cfg"
DEFAULT_POLY = HERE / "boundary_remesh.poly"

BOUNDZ1 = 32
BOUNDN0 = 64


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run a forced-remesh arbitrary-boundary cache regression. "
            "The executable must be an HDF5-enabled 2-D build."
        )
    )
    parser.add_argument("--exe", type=Path, default=DEFAULT_EXE)
    parser.add_argument("--cfg", type=Path, default=DEFAULT_CFG)
    parser.add_argument("--poly", type=Path, default=DEFAULT_POLY)
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
        raise AssertionError(f"missing HDF5 output: {path}")
    with h5py.File(path, "r") as output:
        if dataset not in output:
            raise AssertionError(f"{path}: missing dataset {dataset!r}")
        return output[dataset][...]


def boundary_nodes(flags, mask: int) -> list[int]:
    nodes = [i for i, flag in enumerate(flags) if int(flag) & mask]
    if len(nodes) < 2:
        raise AssertionError(
            f"fixture expected at least two nodes on boundary mask {mask}"
        )
    return nodes


def line_normal(coords, nodes: list[int]) -> tuple[float, float]:
    first = nodes[0]
    second = max(
        nodes,
        key=lambda node: math.dist(coords[first], coords[node]),
    )
    dx = float(coords[second][0] - coords[first][0])
    dz = float(coords[second][1] - coords[first][1])
    length = math.hypot(dx, dz)
    if not length > 0.0:
        raise AssertionError("arbitrary boundary has zero length")
    return dz / length, -dx / length


def normal_alignment(
    lhs: tuple[float, float], rhs: tuple[float, float]
) -> float:
    return abs(lhs[0] * rhs[0] + lhs[1] * rhs[1])


def check_projection(case_dir: Path, model: str) -> None:
    initial_path = case_dir / f"{model}.save.000000.vtkhdf"
    post_remesh_path = case_dir / f"{model}.save.000002.vtkhdf"
    final_path = case_dir / f"{model}.save.000003.vtkhdf"

    initial_flags = read_array(initial_path, "bcflag")
    initial_coord = read_array(initial_path, "coordinate")
    post_flags = read_array(post_remesh_path, "bcflag")
    post_coord = read_array(post_remesh_path, "coordinate")
    final_flags = read_array(final_path, "bcflag")
    final_velocity = read_array(final_path, "velocity")

    initial_normal = line_normal(
        initial_coord, boundary_nodes(initial_flags, BOUNDN0)
    )
    post_nodes = boundary_nodes(post_flags, BOUNDN0)
    post_normal = line_normal(post_coord, post_nodes)
    alignment = normal_alignment(initial_normal, post_normal)
    if alignment > 0.999:
        raise AssertionError(
            "fixture did not rotate the arbitrary boundary enough to distinguish "
            f"old and remeshed caches (|n0 dot n1|={alignment:.17e})"
        )

    checked = 0
    max_speed = 0.0
    max_current_residual = 0.0
    max_stale_residual = 0.0
    for node in boundary_nodes(final_flags, BOUNDN0):
        # The top condition deliberately rotates the boundary during step 1.
        # It remains an intersection in the static boundary-type table, so use
        # the other BOUNDN0 nodes to isolate the arbitrary-normal projection.
        if int(final_flags[node]) & BOUNDZ1:
            continue
        vx = float(final_velocity[node][0])
        vz = float(final_velocity[node][1])
        speed = math.hypot(vx, vz)
        current_residual = abs(vx * post_normal[0] + vz * post_normal[1])
        stale_residual = abs(vx * initial_normal[0] + vz * initial_normal[1])
        max_speed = max(max_speed, speed)
        max_current_residual = max(max_current_residual, current_residual)
        max_stale_residual = max(max_stale_residual, stale_residual)
        checked += 1

    if checked == 0:
        raise AssertionError("fixture has no isolated arbitrary-boundary node")
    if not max_speed > 1.0e-12:
        raise AssertionError("fixture produced no measurable boundary velocity")
    tolerance = max(1.0e-12, max_speed * 1.0e-10)
    if max_current_residual > tolerance:
        raise AssertionError(
            "post-remesh velocity was projected with stale boundary geometry: "
            f"current-normal residual={max_current_residual:.17e}, "
            f"tolerance={tolerance:.17e}, "
            f"old-normal residual={max_stale_residual:.17e}"
        )
    if not max_stale_residual > tolerance * 100.0:
        raise AssertionError(
            "fixture does not discriminate the old cache from the rebuilt cache: "
            f"old-normal residual={max_stale_residual:.17e}"
        )


def run_case(
    exe: Path, template: str, poly: Path, run_root: Path, threads: int
) -> None:
    case_dir = run_root / f"omp{threads}"
    case_dir.mkdir(parents=True, exist_ok=False)
    model = f"boundary_remesh_omp{threads}"
    cfg = template.replace("__MODELNAME__", model)
    if "__" in cfg:
        raise AssertionError("unexpanded config token remains")
    (case_dir / "input.cfg").write_text(cfg, encoding="ascii")
    shutil.copyfile(poly, case_dir / "boundary_remesh.poly")

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    result = subprocess.run(
        [str(exe), str(case_dir / "input.cfg")],
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
            f"boundary-remesh OMP={threads}: DES exited with {result.returncode}"
        )
    # Some public-master remesh schedules perform another quality-triggered
    # remesh after the final checked step.  The regression isolates the first
    # post-remesh projection through save frames 2 and 3, so require that cycle
    # without coupling the assertion to later housekeeping.
    remesh_count = result.stdout.count("Remeshing starts")
    if remesh_count < 1:
        raise AssertionError(
            f"OMP={threads}: expected at least one forced remesh, got "
            f"{remesh_count}"
        )
    check_projection(case_dir, model)
    print(f"boundary-remesh cache regression OMP={threads}: PASS")


def run_all(
    exe: Path,
    template: str,
    poly: Path,
    run_root: Path,
    threads: list[int],
) -> None:
    seen: set[int] = set()
    for thread_count in threads:
        if thread_count <= 0 or thread_count in seen:
            raise ValueError("--threads values must be unique positive integers")
        seen.add(thread_count)
        run_case(exe, template, poly, run_root, thread_count)


def main() -> None:
    args = parse_args()
    exe = args.exe.expanduser().resolve()
    cfg = args.cfg.expanduser().resolve()
    poly = args.poly.expanduser().resolve()
    for path in (exe, cfg, poly):
        if not path.is_file():
            raise FileNotFoundError(path)
    template = cfg.read_text(encoding="ascii")

    if args.run_dir:
        run_root = args.run_dir.expanduser().resolve()
        run_root.mkdir(parents=True, exist_ok=True)
        run_all(exe, template, poly, run_root, args.threads)
    else:
        with tempfile.TemporaryDirectory(prefix="des-boundary-remesh-") as tmp:
            run_all(exe, template, poly, Path(tmp), args.threads)


if __name__ == "__main__":
    main()
