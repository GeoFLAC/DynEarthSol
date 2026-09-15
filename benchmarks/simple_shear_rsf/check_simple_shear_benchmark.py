#!/usr/bin/env python3
"""Run paper benchmark cases, check their references, and clean the outputs."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True

from benchmark_reference import (
    analytical_stress,
    healing_reference,
    load_monitor_case,
    pointwise_error_metrics,
    recovered_rate_errors,
)
from run_simple_shear_benchmark import (
    BENCHMARK_CASES,
    HEALING_CASE,
    BenchmarkCase,
)


DEFAULT_CASES = (
    "steady_ab_neg_v0_1e-6",
    "aging_ab_neg_dc_1e-3",
)
HEALING_RELATIVE_TOLERANCE = 1.0e-12


def cleanup_pycache(root_dir: Path) -> None:
    for path in root_dir.rglob("__pycache__"):
        if path.is_dir():
            shutil.rmtree(path, ignore_errors=True)


def selected_cases(case_names: list[str], use_all: bool) -> list[BenchmarkCase]:
    if use_all:
        return list(BENCHMARK_CASES)
    cases_by_name = {case.name: case for case in BENCHMARK_CASES}
    missing = [name for name in case_names if name not in cases_by_name]
    if missing:
        raise ValueError(f"Unknown benchmark case(s): {', '.join(missing)}")
    return [cases_by_name[name] for name in case_names]


def check_healing(case_dir: Path) -> tuple[float, float]:
    paths = sorted(case_dir.glob("monitor_point_*.csv"))
    if len(paths) != 2:
        raise FileNotFoundError(f"Expected two monitor CSVs in {case_dir}")

    series: list[list[tuple[int, float, float]]] = []
    for path in paths:
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = [
                (
                    int(row["step"]),
                    float(row["time_s"]),
                    float(row["state_variable"]),
                )
                for row in csv.DictReader(handle)
            ]
        rows.sort(key=lambda item: item[0])
        if not rows:
            raise ValueError(f"No monitor rows found in {path}")
        series.append(rows)

    final_step = min(rows[-1][0] for rows in series)
    if final_step != 1_546:
        raise ValueError(f"Unexpected healing final step: {final_step}")
    reference = healing_reference(final_step)

    errors: list[float] = []
    final_states: list[float] = []
    for rows in series:
        for step, time_s, state_s in rows:
            if abs(time_s - step * 259_200.0) > 1.0e-6:
                raise ValueError(f"Healing time does not match step {step}.")
            errors.append(
                abs(state_s - reference[step])
                / max(abs(reference[step]), 1.0e-300)
            )
        final_states.append(rows[-1][2])

    mean_final_state = sum(final_states) / len(final_states)
    theta0_s = HEALING_CASE.characteristic_distance / HEALING_CASE.characteristic_velocity
    final_time_s = final_step * 259_200.0
    linear_deficit = (
        theta0_s + final_time_s - mean_final_state
    ) / final_time_s
    return max(errors), linear_deficit


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    runner = script_dir / "run_simple_shear_benchmark.py"

    parser = argparse.ArgumentParser(
        description="Check the paper's local EP/RSF benchmark references."
    )
    parser.add_argument("--exe", default=None, help="Path to dynearthsol2d. If omitted, auto-detect.")
    parser.add_argument(
        "--cases",
        nargs="+",
        default=list(DEFAULT_CASES),
        help="Representative simple-shear cases to run and verify.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run all nine simple-shear cases instead of the representative subset.",
    )
    parser.add_argument(
        "--max-relative-error",
        type=float,
        default=2.0e-5,
        help="Maximum allowed pointwise stress error as a fraction.",
    )
    parser.add_argument(
        "--max-rate-relative-error",
        type=float,
        default=5.0e-4,
        help="Maximum allowed error in the rate recovered from friction outputs.",
    )
    parser.add_argument(
        "--keep-output",
        action="store_true",
        help="Keep the generated output directory for debugging.",
    )
    args = parser.parse_args()

    cases = selected_cases(args.cases, args.all)
    temp_ctx: tempfile.TemporaryDirectory[str] | None = None
    if args.keep_output:
        output_root = Path(tempfile.mkdtemp(prefix="simple_shear_rsf_check_"))
    else:
        temp_ctx = tempfile.TemporaryDirectory(prefix="simple_shear_rsf_check_")
        output_root = Path(temp_ctx.name)

    try:
        run_cmd = [
            "python3",
            str(runner),
            "--clean",
            "--healing",
            "--output-root",
            str(output_root),
            "--cases",
            *[case.name for case in cases],
        ]
        if args.exe:
            run_cmd.extend(["--exe", str(Path(args.exe).expanduser().resolve())])
        subprocess.run(run_cmd, cwd=script_dir, check=True)

        failures: list[str] = []
        for case in cases:
            data = load_monitor_case(output_root / case.name)
            reference = analytical_stress(case, data.time_s)
            metrics = pointwise_error_metrics(data.mean_abs_stress, reference)
            print(
                f"[stress] {case.name}: "
                f"mean={100.0 * metrics.mean_fraction:.3e}%, "
                f"max={100.0 * metrics.max_fraction:.3e}%"
            )
            if metrics.max_fraction > args.max_relative_error:
                failures.append(
                    f"{case.name} stress={metrics.max_fraction:.3e}"
                )

            rate_errors = recovered_rate_errors(case, data)
            if rate_errors:
                max_rate_error = max(rate_errors)
                print(
                    f"[rate]   {case.name}: "
                    f"max relative error={max_rate_error:.3e}"
                )
                if max_rate_error > args.max_rate_relative_error:
                    failures.append(
                        f"{case.name} rate={max_rate_error:.3e}"
                    )

        healing_error, linear_deficit = check_healing(
            output_root / HEALING_CASE.name
        )
        print(
            "[healing] "
            f"max relative error={healing_error:.3e}, "
            f"linear-growth deficit={linear_deficit:.3e}"
        )
        if healing_error > HEALING_RELATIVE_TOLERANCE:
            failures.append(f"healing={healing_error:.3e}")

        if failures:
            raise SystemExit("Benchmark check failed: " + ", ".join(failures))
    finally:
        if args.keep_output:
            print(f"[kept] {output_root}")
        elif temp_ctx is not None:
            temp_ctx.cleanup()


if __name__ == "__main__":
    _script_dir = Path(__file__).resolve().parent
    try:
        main()
    finally:
        cleanup_pycache(_script_dir)
