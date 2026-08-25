#!/usr/bin/env python3

"""Compare monitor histories from equivalent fluid-source schedules."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("constant_rate_dir", type=Path)
    parser.add_argument("total_amount_dir", type=Path)
    parser.add_argument("source_disabled_dir", type=Path)
    parser.add_argument("--rtol", type=float, default=1.0e-10)
    parser.add_argument("--atol", type=float, default=1.0e-8)
    return parser.parse_args()


def read_monitor_pressures(case_dir: Path) -> list[tuple[list[float], list[float]]]:
    paths = sorted(case_dir.glob("monitor_point_*.csv"))
    if not paths:
        raise FileNotFoundError(f"No monitor CSV files found in {case_dir}")

    series: list[tuple[list[float], list[float]]] = []
    for path in paths:
        with path.open("r", encoding="ascii", newline="") as handle:
            rows = list(csv.DictReader(handle))
        if len(rows) < 2:
            raise ValueError(f"Expected initial and evolved monitor rows in {path}")
        time_s = [float(row["time_s"]) for row in rows]
        pressure = [float(row["pore_pressure"]) for row in rows]
        series.append((time_s, pressure))
    return series


def compare_pressures(
    lhs: list[tuple[list[float], list[float]]],
    rhs: list[tuple[list[float], list[float]]],
    baseline: list[tuple[list[float], list[float]]],
    rtol: float,
    atol: float,
) -> tuple[float, float, float]:
    if len(lhs) != len(rhs) or len(lhs) != len(baseline):
        raise AssertionError(
            f"Monitor count differs: {len(lhs)}, {len(rhs)}, {len(baseline)}"
        )

    max_abs = 0.0
    max_rel = 0.0
    for point_id, (lhs_item, rhs_item, baseline_item) in enumerate(
        zip(lhs, rhs, baseline)
    ):
        lhs_time, lhs_p = lhs_item
        rhs_time, rhs_p = rhs_item
        baseline_time, baseline_p = baseline_item
        if lhs_time != rhs_time or lhs_time != baseline_time:
            raise AssertionError(f"Point {point_id} time axes differ")
        if len(lhs_p) != len(rhs_p) or len(lhs_p) != len(baseline_p):
            raise AssertionError(f"Point {point_id} pressure lengths differ")
        for sample_id, (a, b) in enumerate(zip(lhs_p, rhs_p)):
            baseline_value = baseline_p[sample_id]
            if (not math.isfinite(a) or not math.isfinite(b) or
                    not math.isfinite(baseline_value)):
                raise AssertionError(
                    f"Non-finite pressure at point {point_id}, sample {sample_id}: "
                    f"{a}, {b}, {baseline_value}"
                )
            abs_err = abs(a - b)
            scale = max(abs(a), abs(b), 1.0)
            rel_err = abs_err / scale
            max_abs = max(max_abs, abs_err)
            max_rel = max(max_rel, rel_err)
            if abs_err > atol + rtol * scale:
                raise AssertionError(
                    f"Pressure mismatch at point {point_id}, sample {sample_id}: "
                    f"{a:.16e} vs {b:.16e}; abs={abs_err:.3e}, rel={rel_err:.3e}"
                )

    final_center = lhs[0][1][-1]
    baseline_final_center = baseline[0][1][-1]
    if not math.isfinite(baseline_final_center):
        raise AssertionError("Source-disabled center pore pressure is non-finite")
    source_effect = final_center - baseline_final_center
    effect_tolerance = atol + rtol * max(
        abs(final_center), abs(baseline_final_center), 1.0
    )
    if not source_effect > effect_tolerance:
        raise AssertionError(
            f"Expected injection to exceed source-disabled center pressure: "
            f"{final_center:.16e} vs {baseline_final_center:.16e}"
        )
    return max_abs, max_rel, source_effect


def main() -> None:
    args = parse_args()
    constant = read_monitor_pressures(args.constant_rate_dir)
    total = read_monitor_pressures(args.total_amount_dir)
    baseline = read_monitor_pressures(args.source_disabled_dir)
    max_abs, max_rel, center_source_effect = compare_pressures(
        constant, total, baseline, args.rtol, args.atol
    )
    print(f"max_abs_error={max_abs:.8e}")
    print(f"max_rel_error={max_rel:.8e}")
    print(f"center_pressure_source_effect={center_source_effect:.8e}")


if __name__ == "__main__":
    main()
