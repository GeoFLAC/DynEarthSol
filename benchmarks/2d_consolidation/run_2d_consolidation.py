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
from plot_2d_consolidation import plot_results


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run the 2D consolidation benchmarks and generate comparison plots."
    )
    parser.add_argument(
        "--case",
        nargs="+",
        choices=list(CASE_ORDER) + ["all"],
        default=["all"],
        help="Benchmark cases to run. Defaults to all.",
    )
    parser.add_argument(
        "--exe",
        type=Path,
        default=REPO_ROOT / "dynearthsol2d",
        help="Path to the 2D DynEarthSol executable. Defaults to repo_root/dynearthsol2d.",
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


def run_case(case_name, exe_path, output_root, plot_only):
    config = CASE_CONFIG[case_name]
    cfg_path = SCRIPT_DIR / config["cfg"]
    run_dir = output_root / config["run_subdir"]
    run_dir.mkdir(parents=True, exist_ok=True)

    if not plot_only:
        cmd = [str(exe_path), str(cfg_path)]
        log_path = run_dir / "run.log"
        print(f"[run] {case_name}: {' '.join(cmd)}")
        with log_path.open("w") as log_file:
            subprocess.run(cmd, cwd=run_dir, stdout=log_file, stderr=subprocess.STDOUT, check=True)

    metrics = plot_results(
        case_name,
        run_dir=run_dir,
        model_name=config["model"],
        output_path=run_dir / "final_comparison.png",
    )
    return run_dir, metrics


def main():
    args = parse_args()
    exe_path = args.exe.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    cases = resolve_cases(args.case)
    if not args.plot_only and not exe_path.exists():
        raise FileNotFoundError(f"Could not find 2D executable: {exe_path}")

    summary = []
    for case_name in cases:
        run_dir, metrics = run_case(case_name, exe_path, output_root, args.plot_only)
        summary.append((case_name, run_dir, metrics))

    print("\nSummary")
    for case_name, run_dir, metrics in summary:
        print(
            f"  {case_name}: first/max/final absolute error = "
            f"{metrics['first_error_bar']:.3f} / "
            f"{metrics['max_error_bar']:.3f} / "
            f"{metrics['final_error_bar']:.3f} bar "
            f"({run_dir})"
        )


if __name__ == "__main__":
    main()
