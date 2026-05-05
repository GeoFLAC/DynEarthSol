#!/usr/bin/env python3
"""Run a small simple-shear RSF subset, verify against analytics, and clean up outputs."""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True

from run_simple_shear_benchmark import BENCHMARK_CASES, BenchmarkCase


VX_TOP = 1e-5
SHEAR_MODULUS = 200.0e6
COHESION = 1.0e6
DT = 1.0


def cleanup_pycache(root_dir: Path) -> None:
    for path in root_dir.rglob("__pycache__"):
        if path.is_dir():
            shutil.rmtree(path, ignore_errors=True)


def effective_velocity() -> float:
    v1 = VX_TOP / 3.0
    v2 = 2.0 * VX_TOP / 3.0
    return math.sqrt(v1 * v2)


def analytical_ep(phi_deg: float, total_time: float) -> tuple[list[float], list[float]]:
    phi = math.radians(phi_deg)
    sf = math.sin(phi)
    nphi = (1.0 + sf) / (1.0 - sf)
    npsi = 1.0

    nstep = int(total_time / DT) + 1
    time_s = [DT * i for i in range(nstep)]
    stress_xy = [0.0] * nstep

    for i in range(1, nstep):
        d_exy = 0.5 * VX_TOP * DT
        d_stress_el = 2.0 * SHEAR_MODULUS * d_exy
        stress_el = stress_xy[i - 1] + d_stress_el

        s1_el = -stress_el
        s3_el = stress_el
        yield_fn = s1_el - s3_el * nphi + 2.0 * COHESION * math.sqrt(nphi)

        if yield_fn > 0.0:
            stress_xy[i] = stress_el
        else:
            d_beta = yield_fn / (2.0 * SHEAR_MODULUS * (1.0 + nphi * npsi))
            s3_bar = s3_el + 2.0 * SHEAR_MODULUS * d_beta * npsi
            stress_xy[i] = s3_bar

    return time_s, stress_xy


def analytical_rsf_steady(
    phi_deg: float,
    a: float,
    b: float,
    v0: float,
    total_time: float,
) -> tuple[list[float], list[float]]:
    phi0 = math.radians(phi_deg)
    mu0 = math.tan(phi0)
    velocity = effective_velocity()

    nstep = int(total_time / DT) + 1
    time_s = [DT * i for i in range(nstep)]
    stress_xy = [0.0] * nstep

    mu_ss = max(mu0 + (a - b) * math.log(velocity / v0), 1.0e-6)
    phi_eff = math.atan(mu_ss)
    sf_eff = math.sin(phi_eff)
    nphi_eff = (1.0 + sf_eff) / (1.0 - sf_eff)
    npsi = 1.0

    for i in range(1, nstep):
        d_exy = 0.5 * VX_TOP * DT
        d_stress_el = 2.0 * SHEAR_MODULUS * d_exy
        stress_el = stress_xy[i - 1] + d_stress_el

        s1_el = -stress_el
        s3_el = stress_el
        yield_fn = s1_el - nphi_eff * s3_el + 2.0 * COHESION * math.sqrt(nphi_eff)

        if yield_fn > 0.0:
            stress_xy[i] = stress_el
        else:
            d_beta = yield_fn / (2.0 * SHEAR_MODULUS * (1.0 + nphi_eff * npsi))
            s3_bar = s3_el + 2.0 * SHEAR_MODULUS * d_beta * npsi
            stress_xy[i] = s3_bar

    return time_s, stress_xy


def analytical_rsf_aging(
    phi_deg: float,
    a: float,
    b: float,
    dc: float,
    v0: float,
    total_time: float,
) -> tuple[list[float], list[float]]:
    phi0 = math.radians(phi_deg)
    mu0 = math.tan(phi0)
    velocity = max(effective_velocity(), 1.0e-12)
    theta = dc / v0

    nstep = int(total_time / DT) + 1
    time_s = [DT * i for i in range(nstep)]
    stress_xy = [0.0] * nstep
    npsi = 1.0

    for i in range(1, nstep):
        d_exy = 0.5 * VX_TOP * DT
        d_stress_el = 2.0 * SHEAR_MODULUS * d_exy
        stress_el = stress_xy[i - 1] + d_stress_el

        theta = theta + DT * (1.0 - velocity * theta / dc)
        theta = max(theta, 1.0e-30)
        mu = max(mu0 + a * math.log(velocity / v0) + b * math.log((theta * v0) / dc), 1.0e-6)
        phi_eff = math.atan(mu)
        sf_eff = math.sin(phi_eff)
        nphi_eff = (1.0 + sf_eff) / (1.0 - sf_eff)

        s1_el = -stress_el
        s3_el = stress_el
        yield_fn = s1_el - nphi_eff * s3_el + 2.0 * COHESION * math.sqrt(nphi_eff)

        if yield_fn > 0.0:
            stress_xy[i] = stress_el
        else:
            d_beta = yield_fn / (2.0 * SHEAR_MODULUS * (1.0 + nphi_eff * npsi))
            s3_bar = s3_el + 2.0 * SHEAR_MODULUS * d_beta * npsi
            stress_xy[i] = s3_bar

    return time_s, stress_xy


def analytical_series(case: BenchmarkCase, total_time: float) -> tuple[list[float], list[float]]:
    if case.group == "ep":
        return analytical_ep(case.friction_angle_deg, total_time)
    if case.group in ("steady_ab_pos", "steady_ab_neg"):
        return analytical_rsf_steady(
            phi_deg=case.friction_angle_deg,
            a=case.direct_a,
            b=case.evolution_b,
            v0=case.characteristic_velocity,
            total_time=total_time,
        )
    if case.group == "aging_ab_neg":
        return analytical_rsf_aging(
            phi_deg=case.friction_angle_deg,
            a=case.direct_a,
            b=case.evolution_b,
            dc=case.characteristic_distance,
            v0=case.characteristic_velocity,
            total_time=total_time,
        )
    raise ValueError(f"Unsupported group: {case.group}")


def load_monitor_stress(case_dir: Path) -> tuple[list[float], list[float]]:
    csv_paths = sorted(case_dir.glob("monitor_point_*.csv"))
    if not csv_paths:
        raise FileNotFoundError(f"No monitor CSV found in {case_dir}")

    series: list[tuple[list[float], list[float]]] = []
    for csv_path in csv_paths:
        with csv_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            rows = sorted(reader, key=lambda row: float(row["time_s"]))
        time_s = [float(row["time_s"]) for row in rows]
        stress_xy = [abs(float(row["stress_2"])) for row in rows]
        if not time_s:
            raise ValueError(f"No monitor rows found in {csv_path}")
        series.append((time_s, stress_xy))

    ref_time = series[0][0]
    for idx, (time_s, _) in enumerate(series[1:], start=1):
        if len(time_s) != len(ref_time):
            raise ValueError(f"Monitor time axis length mismatch between point 0 and point {idx} in {case_dir}")
        for t0, t1 in zip(ref_time, time_s):
            if abs(t0 - t1) > 1.0e-12:
                raise ValueError(f"Monitor time axis mismatch between point 0 and point {idx} in {case_dir}")

    mean_abs_stress: list[float] = []
    for sample_index in range(len(ref_time)):
        mean_abs_stress.append(
            sum(stress_series[sample_index] for _, stress_series in series) / float(len(series))
        )

    return ref_time, mean_abs_stress


def interpolate_linear(x_values: list[float], y_values: list[float], x: float) -> float:
    if x <= x_values[0]:
        return y_values[0]
    if x >= x_values[-1]:
        return y_values[-1]

    hi = 1
    while x_values[hi] < x:
        hi += 1
    lo = hi - 1
    x0 = x_values[lo]
    x1 = x_values[hi]
    y0 = y_values[lo]
    y1 = y_values[hi]
    if x1 == x0:
        return y0
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def compare_series(
    num_time: list[float],
    num_stress: list[float],
    ana_time: list[float],
    ana_stress: list[float],
) -> float:
    ana_on_num = [interpolate_linear(ana_time, ana_stress, x) for x in num_time]
    scale = max(max(abs(value) for value in ana_on_num), 1.0)
    return max(abs(num - ana) for num, ana in zip(num_stress, ana_on_num)) / scale


def selected_cases(case_names: list[str]) -> list[BenchmarkCase]:
    cases_by_name = {case.name: case for case in BENCHMARK_CASES}
    missing = [name for name in case_names if name not in cases_by_name]
    if missing:
        raise ValueError(f"Unknown benchmark case(s): {', '.join(missing)}")
    return [cases_by_name[name] for name in case_names]


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    runner = script_dir / "run_simple_shear_benchmark.py"

    parser = argparse.ArgumentParser(description="Check a small RSF simple-shear subset and clean outputs.")
    parser.add_argument("--exe", default=None, help="Path to dynearthsol2d. If omitted, auto-detect.")
    parser.add_argument(
        "--cases",
        nargs="+",
        default=["steady_ab_pos_v0_1e-6", "aging_ab_neg_dc_1e-6"],
        help="Representative cases to run and verify.",
    )
    parser.add_argument(
        "--max-relative-error",
        type=float,
        default=5e-2,
        help="Relative error tolerance for CI-style checks.",
    )
    parser.add_argument(
        "--keep-output",
        action="store_true",
        help="Keep the generated output directory for debugging.",
    )
    args = parser.parse_args()

    cases = selected_cases(args.cases)
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
            "--output-step-interval",
            "40",
            "--output-root",
            str(output_root),
            "--cases",
            *args.cases,
        ]
        if args.exe:
            run_cmd.extend(["--exe", args.exe])
        subprocess.run(run_cmd, cwd=script_dir, check=True)

        failures: list[tuple[str, float]] = []
        for case in cases:
            case_dir = output_root / case.name
            num_time, num_stress = load_monitor_stress(case_dir)
            ana_time, ana_stress = analytical_series(case, total_time=num_time[-1])
            rel_error = compare_series(num_time, num_stress, ana_time, ana_stress)
            print(f"[ok] {case.name}: max relative error = {rel_error:.3e}")
            if rel_error > args.max_relative_error:
                failures.append((case.name, rel_error))

        if failures:
            details = ", ".join(f"{name}={err:.3e}" for name, err in failures)
            raise SystemExit(f"Benchmark check failed: {details}")
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
