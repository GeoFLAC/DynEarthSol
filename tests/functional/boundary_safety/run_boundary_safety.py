#!/usr/bin/env python3

"""Exercise arbitrary-normal VBC geometry and projection guards."""

from __future__ import annotations

import argparse
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
DEFAULT_CFG = HERE / "boundary_safety_base.cfg"
BOUNDN0 = 64
BOUNDN1 = 128
POINT_DATA = "/VTKHDF/grid/PointData"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe-2d", type=Path, default=REPO_ROOT / "dynearthsol2d")
    parser.add_argument("--exe-3d", type=Path, default=REPO_ROOT / "dynearthsol3d")
    parser.add_argument("--cfg", type=Path, default=DEFAULT_CFG)
    parser.add_argument(
        "--run-dir",
        type=Path,
        help="Keep outputs below this new directory (default: temporary).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        nargs="+",
        default=[1, 4],
        help="OMP thread counts to exercise (default: 1 4).",
    )
    return parser.parse_args()


def render_config(
    template: str,
    *,
    model: str,
    poly_name: str,
    vbc_n0: int,
    vbc_n1: int = 1,
    vbc_val_n0: float = 0.0,
    vbc_z1: int = 0,
    initial_checkpoint: bool = False,
    restarting: bool = False,
    restart_model: str = "unused",
) -> str:
    replacements = {
        "__MODELNAME__": model,
        "__POLY_FILENAME__": poly_name,
        "__HAS_INITIAL_CHECKPOINT__": "yes" if initial_checkpoint else "no",
        "__IS_RESTARTING__": "yes" if restarting else "no",
        "__RESTART_MODEL__": restart_model,
        "__VBC_N0__": str(vbc_n0),
        "__VBC_N1__": str(vbc_n1),
        "__VBC_VAL_N0__": f"{vbc_val_n0:.17e}",
        "__VBC_Z1__": str(vbc_z1),
    }
    result = template
    for token, value in replacements.items():
        if token not in result:
            raise AssertionError(f"missing config token {token}")
        result = result.replace(token, value)
    if "__" in result:
        raise AssertionError("unexpanded config token remains")
    return result


def run_des(
    exe: Path,
    config: str,
    case_dir: Path,
    name: str,
    threads: int,
    *,
    expected_code: int = 0,
    expected_text: str = "",
) -> str:
    cfg_path = case_dir / f"{name}.cfg"
    cfg_path.write_text(config, encoding="ascii")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    completed = subprocess.run(
        [str(exe), str(cfg_path)],
        cwd=case_dir,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=120,
        check=False,
    )
    (case_dir / f"{name}.log").write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != expected_code:
        print(completed.stdout, file=sys.stderr)
        raise AssertionError(
            f"{name} OMP={threads}: exit {completed.returncode}, "
            f"expected {expected_code}"
        )
    if expected_text and expected_text not in completed.stdout:
        raise AssertionError(
            f"{name} OMP={threads}: missing diagnostic {expected_text!r}"
        )
    return completed.stdout


def copy_poly(case_dir: Path, source_name: str) -> str:
    source = HERE / source_name
    target = case_dir / source.name
    shutil.copyfile(source, target)
    return target.name


def set_checkpoint_velocity(save_path: Path, values: tuple[float, ...]) -> None:
    storage = f"{POINT_DATA}/velocity"
    with h5py.File(save_path, "r+") as output:
        if storage not in output:
            raise AssertionError(f"{save_path}: missing storage {storage!r}")
        velocity = output[storage][...]
        if velocity.ndim != 2 or velocity.shape[1] < len(values):
            raise AssertionError(
                f"{save_path}: velocity shape {velocity.shape} does not provide "
                f"{len(values)} components"
            )
        for component, value in enumerate(values):
            velocity[:, component] = value
        output[storage][...] = velocity


def read_root(save_path: Path, name: str):
    if not save_path.is_file():
        raise AssertionError(f"missing output {save_path}")
    with h5py.File(save_path, "r") as output:
        if name not in output:
            raise AssertionError(f"{save_path}: missing root dataset {name!r}")
        return output[name][...]


def check_type11_horizontal_projection(
    exe: Path, template: str, run_root: Path, threads: int
) -> None:
    case_dir = run_root / f"normalization_omp{threads}"
    case_dir.mkdir(parents=True, exist_ok=False)
    poly_name = copy_poly(case_dir, "oblique_2d.poly")

    source_model = f"boundary_source_omp{threads}"
    source = render_config(
        template,
        model=source_model,
        poly_name=poly_name,
        vbc_n0=1,
        initial_checkpoint=True,
    )
    run_des(exe, source, case_dir, "source", threads)
    source_save = case_dir / f"{source_model}.save.000000.vtkhdf"
    set_checkpoint_velocity(source_save, (3.0, -4.0))

    restart_model = f"boundary_type11_omp{threads}"
    restart = render_config(
        template,
        model=restart_model,
        poly_name=poly_name,
        vbc_n0=11,
        vbc_val_n0=2.0,
        restarting=True,
        restart_model=source_model,
    )
    run_des(exe, restart, case_dir, "restart", threads)
    restart_save = case_dir / f"{restart_model}.save.000000.vtkhdf"
    flags = read_root(restart_save, "bcflag")
    velocity = read_root(restart_save, "velocity")
    isolated = [i for i, flag in enumerate(flags) if int(flag) == BOUNDN0]
    if not isolated:
        raise AssertionError("oblique fixture produced no isolated BOUNDN0 node")
    for node in isolated:
        horizontal_speed = abs(float(velocity[node][0]))
        if abs(horizontal_speed - 2.0) > 1.0e-12:
            raise AssertionError(
                "type-11 horizontal-normal projection did not impose the "
                f"configured speed at node {node}: |vx|={horizontal_speed:.17e}"
            )
        if float(velocity[node][1]) != -4.0:
            raise AssertionError(
                f"type-11 projection changed vertical velocity at node {node}"
            )


def check_vertical_normal_rejected(
    exe: Path, template: str, run_root: Path, threads: int, vbc_type: int
) -> None:
    case_dir = run_root / f"vertical_type{vbc_type}_omp{threads}"
    case_dir.mkdir(parents=True, exist_ok=False)
    poly_name = copy_poly(case_dir, "vertical_normal_2d.poly")
    model = f"vertical_type{vbc_type}_omp{threads}"
    config = render_config(
        template,
        model=model,
        poly_name=poly_name,
        vbc_n0=vbc_type,
    )
    run_des(
        exe,
        config,
        case_dir,
        "run",
        threads,
        expected_code=11,
        expected_text="normal has no horizontal projection",
    )


def check_type11_active_intersection_rejected(
    exe: Path, template: str, run_root: Path, threads: int
) -> None:
    case_dir = run_root / f"type11_intersection_omp{threads}"
    case_dir.mkdir(parents=True, exist_ok=False)
    poly_name = copy_poly(case_dir, "oblique_2d.poly")
    model = f"type11_intersection_omp{threads}"
    config = render_config(
        template,
        model=model,
        poly_name=poly_name,
        vbc_n0=11,
        vbc_z1=1,
    )
    run_des(
        exe,
        config,
        case_dir,
        "run",
        threads,
        expected_code=11,
        expected_text="together with active boundary",
    )


def run_2d(exe: Path, template: str, run_root: Path, threads: list[int]) -> None:
    seen: set[int] = set()
    for thread_count in threads:
        if thread_count <= 0 or thread_count in seen:
            raise ValueError("--threads values must be unique positive integers")
        seen.add(thread_count)
        check_type11_horizontal_projection(exe, template, run_root, thread_count)
        check_vertical_normal_rejected(exe, template, run_root, thread_count, 11)
        check_vertical_normal_rejected(exe, template, run_root, thread_count, 13)
        check_type11_active_intersection_rejected(
            exe, template, run_root, thread_count
        )
        print(f"boundary safety 2D OMP={thread_count}: PASS", flush=True)


def check_nonorthogonal_edge_normalized(
    exe: Path, template: str, run_root: Path, threads: int
) -> None:
    case_dir = run_root / f"nonorthogonal_edge_3d_omp{threads}"
    case_dir.mkdir(parents=True, exist_ok=False)
    poly_name = copy_poly(case_dir, "nonorthogonal_edge_3d.poly")

    source_model = f"edge_source_3d_omp{threads}"
    source = render_config(
        template,
        model=source_model,
        poly_name=poly_name,
        vbc_n0=1,
        vbc_n1=1,
        initial_checkpoint=True,
    )
    run_des(exe, source, case_dir, "source", threads)
    source_save = case_dir / f"{source_model}.save.000000.vtkhdf"
    set_checkpoint_velocity(source_save, (2.0, 3.0, 4.0))

    restart_model = f"edge_projection_3d_omp{threads}"
    restart = render_config(
        template,
        model=restart_model,
        poly_name=poly_name,
        vbc_n0=1,
        vbc_n1=1,
        restarting=True,
        restart_model=source_model,
    )
    run_des(exe, restart, case_dir, "restart", threads)
    restart_save = case_dir / f"{restart_model}.save.000000.vtkhdf"
    flags = read_root(restart_save, "bcflag")
    velocity = read_root(restart_save, "velocity")
    intersection_flag = BOUNDN0 | BOUNDN1
    isolated = [
        i for i, flag in enumerate(flags) if int(flag) == intersection_flag
    ]
    if not isolated:
        raise AssertionError(
            "nonorthogonal fixture produced no isolated BOUNDN0/BOUNDN1 edge node"
        )
    expected = (3.0, 0.0, 3.0)
    for node in isolated:
        for component, target in enumerate(expected):
            actual = float(velocity[node][component])
            if abs(actual - target) > 1.0e-12:
                raise AssertionError(
                    "nonorthogonal edge projection was not normalized at "
                    f"node {node}, component {component}: "
                    f"value={actual:.17e}, expected={target:.17e}"
                )


def check_parallel_normals_rejected(
    exe: Path, template: str, run_root: Path, threads: int
) -> None:
    case_dir = run_root / f"parallel_normals_3d_omp{threads}"
    case_dir.mkdir(parents=True, exist_ok=False)
    poly_name = copy_poly(case_dir, "parallel_normals_3d.poly")
    model = f"parallel_normals_3d_omp{threads}"
    config = render_config(
        template,
        model=model,
        poly_name=poly_name,
        vbc_n0=1,
        vbc_n1=1,
    )
    run_des(
        exe,
        config,
        case_dir,
        "run",
        threads,
        expected_code=42,
        expected_text="meet but have no unique edge direction",
    )


def run_3d(exe: Path, template: str, run_root: Path, threads: list[int]) -> None:
    for thread_count in threads:
        check_nonorthogonal_edge_normalized(
            exe, template, run_root, thread_count
        )
        check_parallel_normals_rejected(exe, template, run_root, thread_count)
        print(f"boundary safety 3D OMP={thread_count}: PASS", flush=True)


def main() -> None:
    args = parse_args()
    exe_2d = args.exe_2d.expanduser().resolve()
    exe_3d = args.exe_3d.expanduser().resolve()
    cfg = args.cfg.expanduser().resolve()
    if not exe_2d.is_file() or not os.access(exe_2d, os.X_OK):
        raise FileNotFoundError(f"--exe-2d is not executable: {exe_2d}")
    if not exe_3d.is_file() or not os.access(exe_3d, os.X_OK):
        raise FileNotFoundError(f"--exe-3d is not executable: {exe_3d}")
    template = cfg.read_text(encoding="ascii")

    if args.run_dir is not None:
        run_root = args.run_dir.expanduser().resolve()
        run_root.mkdir(parents=True, exist_ok=False)
        run_2d(exe_2d, template, run_root, args.threads)
        run_3d(exe_3d, template, run_root, args.threads)
    else:
        with tempfile.TemporaryDirectory(prefix="des-boundary-safety-") as tmp:
            run_2d(exe_2d, template, Path(tmp), args.threads)
            run_3d(exe_3d, template, Path(tmp), args.threads)


if __name__ == "__main__":
    main()
