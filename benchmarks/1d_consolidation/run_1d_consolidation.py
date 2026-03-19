#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]

os.environ.setdefault("MPLBACKEND", "Agg")

if str(SCRIPT_DIR) not in sys.path:
    sys.path.append(str(SCRIPT_DIR))

from benchmark_cases import CASE_CONFIG, CASE_ORDER
from plot_1d_consolidation import plot_results


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run the 1D consolidation benchmarks and generate comparison plots."
    )
    parser.add_argument(
        "--case",
        nargs="+",
        choices=list(CASE_ORDER) + ["all"],
        default=["all"],
        help="Benchmark cases to run. Defaults to all.",
    )
    parser.add_argument(
        "--exe2d",
        type=Path,
        default=REPO_ROOT / "dynearthsol2d",
        help="Path to the 2D DynEarthSol executable. Defaults to repo_root/dynearthsol2d.",
    )
    parser.add_argument(
        "--exe3d",
        type=Path,
        default=REPO_ROOT / "dynearthsol3d",
        help="Path to the 3D DynEarthSol executable. Defaults to repo_root/dynearthsol3d.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=SCRIPT_DIR / "runs",
        help="Directory where benchmark outputs are written.",
    )
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Skip the solver run and regenerate plots/metrics from existing outputs.",
    )
    return parser.parse_args()


def resolve_cases(selected_cases):
    if "all" in selected_cases:
        return list(CASE_ORDER)
    return selected_cases


def select_executable(config, exe2d_path, exe3d_path):
    return exe3d_path if config["dim"] == 3 else exe2d_path


def run_case(case_name, exe2d_path, exe3d_path, output_root, plot_only):
    config = CASE_CONFIG[case_name]
    cfg_path = SCRIPT_DIR / config["cfg"]
    run_dir = output_root / config["run_subdir"]
    run_dir.mkdir(parents=True, exist_ok=True)

    if not plot_only:
        exe_path = select_executable(config, exe2d_path, exe3d_path)
        cmd = [str(exe_path), str(cfg_path)]
        log_path = run_dir / "run.log"
        print(f"[run] {case_name}: {' '.join(cmd)}")
        with log_path.open("w") as log_file:
            subprocess.run(cmd, cwd=run_dir, stdout=log_file, stderr=subprocess.STDOUT, check=True)

    metrics = plot_results(
        [
            {
                "case_name": case_name,
                "model": config["model"],
                "label": config["label"],
                "dim": config["dim"],
                "skip_first_output": config["skip_first_output"],
            }
        ],
        dim=config["dim"],
        output_path=run_dir / "final_comparison.png",
        run_dir=run_dir,
    )
    return run_dir, metrics[0]


def main():
    args = parse_args()
    exe2d_path = args.exe2d.resolve()
    exe3d_path = args.exe3d.resolve()

    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    cases = resolve_cases(args.case)
    if not args.plot_only:
        required_dims = {CASE_CONFIG[case_name]["dim"] for case_name in cases}
        if 2 in required_dims and not exe2d_path.exists():
            raise FileNotFoundError(f"Could not find 2D executable: {exe2d_path}")
        if 3 in required_dims and not exe3d_path.exists():
            raise FileNotFoundError(f"Could not find 3D executable: {exe3d_path}")

    summary = []
    for case_name in cases:
        run_dir, metric = run_case(case_name, exe2d_path, exe3d_path, output_root, args.plot_only)
        summary.append((case_name, run_dir, metric))

    print("\nSummary")
    for case_name, run_dir, metric in summary:
        print(
            f"  {case_name}: first/max/final error = "
            f"{metric['first_error'] * 100:.3f}% / "
            f"{metric['max_error'] * 100:.3f}% / "
            f"{metric['final_error'] * 100:.3f}% "
            f"({run_dir})"
        )


if __name__ == "__main__":
    main()
