#!/usr/bin/env python3

"""Check fixed pore-pressure BCs across initialization, restart, and remesh."""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path


sys.dont_write_bytecode = True

from fixed_hbc_forced_remesh import check_forced_remesh
from fixed_hbc_fresh_restart import check_fresh_and_restart
from fixed_hbc_invalid_reset import check_parser_and_reset_mapping


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
DEFAULT_CFG = HERE / "fixed_hbc_base.cfg"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run fixed pore-pressure boundary lifecycle regressions."
    )
    parser.add_argument("--exe-2d", type=Path, default=REPO_ROOT / "dynearthsol2d")
    parser.add_argument("--exe-3d", type=Path, default=REPO_ROOT / "dynearthsol3d")
    parser.add_argument("--cfg", type=Path, default=DEFAULT_CFG)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 4])
    parser.add_argument(
        "--run-dir",
        type=Path,
        help="Keep outputs below this directory (default: temporary directory).",
    )
    return parser.parse_args()


def _resolve_executable(path: Path, dimensions: int) -> Path:
    exe = path.expanduser().resolve()
    if not exe.is_file():
        raise FileNotFoundError(f"{dimensions}D executable not found: {exe}")
    return exe


def run_dimension(
    exe: Path,
    template: str,
    run_root: Path,
    *,
    dimensions: int,
    threads: list[int],
) -> None:
    dimension_root = run_root / f"{dimensions}d"

    # Parser/reset validation is deterministic and does not exercise shared
    # numerical loops, so one thread is sufficient. Lifecycle cases run on the
    # full requested OpenMP matrix.
    check_parser_and_reset_mapping(
        exe,
        template,
        dimension_root / f"parser_omp{threads[0]}",
        dimensions=dimensions,
        threads=threads[0],
    )
    for thread_count in threads:
        thread_root = dimension_root / f"omp{thread_count}"
        check_fresh_and_restart(
            exe, template, thread_root, threads=thread_count
        )
        check_forced_remesh(exe, template, thread_root, threads=thread_count)
        print(f"fixed HBC: {dimensions}D OMP={thread_count}: PASS")


def run_all(
    exe_2d: Path,
    exe_3d: Path,
    template: str,
    run_root: Path,
    threads: list[int],
) -> None:
    run_dimension(
        exe_2d, template, run_root, dimensions=2, threads=threads
    )
    run_dimension(
        exe_3d, template, run_root, dimensions=3, threads=threads
    )


def main() -> None:
    args = parse_args()
    exe_2d = _resolve_executable(args.exe_2d, 2)
    exe_3d = _resolve_executable(args.exe_3d, 3)
    cfg_path = args.cfg.expanduser().resolve()
    if not cfg_path.is_file():
        raise FileNotFoundError(f"config template not found: {cfg_path}")

    threads = list(dict.fromkeys(args.threads))
    if not threads or any(value < 1 for value in threads):
        raise ValueError("--threads values must be positive integers")
    template = cfg_path.read_text(encoding="ascii")

    if args.run_dir:
        run_root = args.run_dir.expanduser().resolve()
        run_root.mkdir(parents=True, exist_ok=True)
        run_all(exe_2d, exe_3d, template, run_root, threads)
    else:
        with tempfile.TemporaryDirectory(prefix="des-fixed-hbc-") as tmp:
            run_all(exe_2d, exe_3d, template, Path(tmp), threads)

    print("fixed pore-pressure boundary regression: PASS")


if __name__ == "__main__":
    main()
