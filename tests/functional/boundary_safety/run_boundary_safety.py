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
POINT_DATA = "/VTKHDF/grid/PointData"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe-2d", type=Path, default=REPO_ROOT / "dynearthsol2d")
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
    vbc_val_n0: float = 0.0,
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
        "__VBC_VAL_N0__": f"{vbc_val_n0:.17e}",
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


def set_checkpoint_velocity(save_path: Path) -> None:
    storage = f"{POINT_DATA}/velocity"
    with h5py.File(save_path, "r+") as output:
        if storage not in output:
            raise AssertionError(f"{save_path}: missing storage {storage!r}")
        velocity = output[storage][...]
        velocity[:, 0] = 3.0
        velocity[:, 1] = -4.0
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
    set_checkpoint_velocity(source_save)

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


def run_2d(exe: Path, template: str, run_root: Path, threads: list[int]) -> None:
    seen: set[int] = set()
    for thread_count in threads:
        if thread_count <= 0 or thread_count in seen:
            raise ValueError("--threads values must be unique positive integers")
        seen.add(thread_count)
        check_type11_horizontal_projection(exe, template, run_root, thread_count)
        check_vertical_normal_rejected(exe, template, run_root, thread_count, 11)
        check_vertical_normal_rejected(exe, template, run_root, thread_count, 13)
        print(f"boundary safety 2D OMP={thread_count}: PASS", flush=True)


def main() -> None:
    args = parse_args()
    exe = args.exe_2d.expanduser().resolve()
    cfg = args.cfg.expanduser().resolve()
    if not exe.is_file() or not os.access(exe, os.X_OK):
        raise FileNotFoundError(f"--exe-2d is not executable: {exe}")
    template = cfg.read_text(encoding="ascii")

    if args.run_dir is not None:
        run_root = args.run_dir.expanduser().resolve()
        run_root.mkdir(parents=True, exist_ok=False)
        run_2d(exe, template, run_root, args.threads)
    else:
        with tempfile.TemporaryDirectory(prefix="des-boundary-safety-") as tmp:
            run_2d(exe, template, Path(tmp), args.threads)


if __name__ == "__main__":
    main()
