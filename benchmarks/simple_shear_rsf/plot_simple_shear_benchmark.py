#!/usr/bin/env python3
"""Plot monitor-based simple-shear benchmark results against analytical curves."""

from __future__ import annotations

import argparse
import os
import sys
from math import pi, sqrt
from pathlib import Path
from typing import Iterable

sys.dont_write_bytecode = True

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from run_simple_shear_benchmark import BENCHMARK_CASES, GROUP_ORDER, BenchmarkCase


VX_TOP = 1e-5
SHEAR_MODULUS = 200.0e6
COHESION = 1.0e6
DT = 1.0
STRESS_SCALE = 1e-6
LINE_W = 2.2
MARKER_SIZE = 6.5
MARKER_EDGE_W = 1.5


def set_plot_style() -> None:
    try:
        from matplotlib import font_manager

        nimbus_paths = [
            "/usr/share/fonts/opentype/urw-base35/NimbusSans-Regular.otf",
            "/usr/share/fonts/opentype/urw-base35/NimbusSans-Bold.otf",
            "/usr/share/fonts/opentype/urw-base35/NimbusSans-Italic.otf",
        ]
        for path in nimbus_paths:
            if os.path.isfile(path):
                font_manager.fontManager.addfont(path)
    except Exception:
        pass

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Nimbus Sans", "Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 15,
            "axes.labelsize": 18,
            "axes.titlesize": 16,
            "legend.fontsize": 11,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "axes.linewidth": 2.4,
            "lines.linewidth": LINE_W,
        }
    )


def style_axis(ax: plt.Axes, panel_tag: str, title: str) -> None:
    ax.set_title(title, fontweight="bold")
    ax.set_xlabel("Time (s)", fontweight="bold")
    ax.set_ylabel(r"$|\sigma_{xy}|$ (MPa)", fontweight="bold")
    ax.text(
        0.02,
        0.97,
        panel_tag,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    ax.tick_params(axis="both", which="both", width=2.4, length=6.5)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")
    for spine in ax.spines.values():
        spine.set_linewidth(2.4)
    ax.grid(False)


def style_legend(legend: plt.Legend | None) -> None:
    if legend is None:
        return
    legend.set_frame_on(False)
    for text in legend.get_texts():
        text.set_fontweight("bold")


def pow10_label(x: float) -> str:
    if x == 0:
        return "0"
    exp = int(np.round(np.log10(x)))
    return rf"10^{{{exp}}}"


def analytical_ep(phi_deg: float, total_time: float) -> tuple[np.ndarray, np.ndarray]:
    phi = phi_deg * np.pi / 180.0
    sf = np.sin(phi)
    nphi = (1.0 + sf) / (1.0 - sf)
    npsi = 1.0

    nstep = int(total_time / DT) + 1
    time_s = DT * np.arange(nstep, dtype=float)
    stress_xy = np.zeros(nstep, dtype=float)

    for i in range(1, nstep):
        d_exy = 0.5 * VX_TOP * DT
        d_stress_el = 2.0 * SHEAR_MODULUS * d_exy
        stress_el = stress_xy[i - 1] + d_stress_el

        s1_el = -stress_el
        s3_el = stress_el
        yield_fn = s1_el - s3_el * nphi + 2.0 * COHESION * sqrt(nphi)

        if yield_fn > 0.0:
            stress_xy[i] = stress_el
        else:
            d_beta = yield_fn / (2.0 * SHEAR_MODULUS * (1.0 + nphi * npsi))
            s3_bar = s3_el + 2.0 * SHEAR_MODULUS * d_beta * npsi
            stress_xy[i] = s3_bar

    return time_s, stress_xy


def effective_velocity() -> float:
    v1 = VX_TOP / 3.0
    v2 = 2.0 * VX_TOP / 3.0
    return float(np.sqrt(v1 * v2))


def analytical_rsf_steady(phi_deg: float, a: float, b: float, v0: float, total_time: float) -> tuple[np.ndarray, np.ndarray]:
    phi0 = phi_deg * pi / 180.0
    mu0 = np.tan(phi0)
    velocity = effective_velocity()

    nstep = int(total_time / DT) + 1
    time_s = DT * np.arange(nstep, dtype=float)
    stress_xy = np.zeros(nstep, dtype=float)

    mu_ss = max(mu0 + (a - b) * np.log(velocity / v0), 1.0e-6)
    phi_eff = np.arctan(mu_ss)
    sf_eff = np.sin(phi_eff)
    nphi_eff = (1.0 + sf_eff) / (1.0 - sf_eff)
    npsi = 1.0

    for i in range(1, nstep):
        d_exy = 0.5 * VX_TOP * DT
        d_stress_el = 2.0 * SHEAR_MODULUS * d_exy
        stress_el = stress_xy[i - 1] + d_stress_el

        s1_el = -stress_el
        s3_el = stress_el
        yield_fn = s1_el - nphi_eff * s3_el + 2.0 * COHESION * np.sqrt(nphi_eff)

        if yield_fn > 0.0:
            stress_xy[i] = stress_el
        else:
            d_beta = yield_fn / (2.0 * SHEAR_MODULUS * (1.0 + nphi_eff * npsi))
            s3_bar = s3_el + 2.0 * SHEAR_MODULUS * d_beta * npsi
            stress_xy[i] = s3_bar

    return time_s, stress_xy


def analytical_rsf_aging(phi_deg: float, a: float, b: float, dc: float, v0: float, total_time: float) -> tuple[np.ndarray, np.ndarray]:
    phi0 = phi_deg * pi / 180.0
    mu0 = np.tan(phi0)
    velocity = max(effective_velocity(), 1.0e-12)
    theta = dc / v0

    nstep = int(total_time / DT) + 1
    time_s = DT * np.arange(nstep, dtype=float)
    stress_xy = np.zeros(nstep, dtype=float)
    npsi = 1.0

    for i in range(1, nstep):
        d_exy = 0.5 * VX_TOP * DT
        d_stress_el = 2.0 * SHEAR_MODULUS * d_exy
        stress_el = stress_xy[i - 1] + d_stress_el

        theta = theta + DT * (1.0 - velocity * theta / dc)
        theta = max(theta, 1.0e-30)
        mu = max(mu0 + a * np.log(velocity / v0) + b * np.log((theta * v0) / dc), 1.0e-6)
        phi_eff = np.arctan(mu)
        sf_eff = np.sin(phi_eff)
        nphi_eff = (1.0 + sf_eff) / (1.0 - sf_eff)

        s1_el = -stress_el
        s3_el = stress_el
        yield_fn = s1_el - nphi_eff * s3_el + 2.0 * COHESION * np.sqrt(nphi_eff)

        if yield_fn > 0.0:
            stress_xy[i] = stress_el
        else:
            d_beta = yield_fn / (2.0 * SHEAR_MODULUS * (1.0 + nphi_eff * npsi))
            s3_bar = s3_el + 2.0 * SHEAR_MODULUS * d_beta * npsi
            stress_xy[i] = s3_bar

    return time_s, stress_xy


def analytical_series(case: BenchmarkCase, total_time: float) -> tuple[np.ndarray, np.ndarray]:
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


def load_monitor_stress(case_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    csv_paths = sorted(case_dir.glob("monitor_point_*.csv"))
    if not csv_paths:
        raise FileNotFoundError(f"No monitor CSV found in {case_dir}")

    series: list[tuple[np.ndarray, np.ndarray]] = []
    for csv_path in csv_paths:
        arr = np.genfromtxt(csv_path, delimiter=",", names=True, dtype=None, encoding="utf-8")
        arr = np.atleast_1d(arr)
        names = arr.dtype.names or ()
        if "time_s" not in names or "stress_2" not in names:
            raise ValueError(f"Missing required columns in {csv_path}")

        time_s = np.asarray(arr["time_s"], dtype=float)
        stress_xy = np.abs(np.asarray(arr["stress_2"], dtype=float))
        order = np.argsort(time_s)
        series.append((time_s[order], stress_xy[order]))

    ref_time = series[0][0]
    stacked = [series[0][1]]
    for idx, (time_s, stress_xy) in enumerate(series[1:], start=1):
        if time_s.shape != ref_time.shape or not np.allclose(time_s, ref_time):
            raise ValueError(
                f"Monitor time axis mismatch between point 0 and point {idx} in {case_dir}"
            )
        stacked.append(stress_xy)

    mean_abs_stress = np.mean(np.vstack(stacked), axis=0)
    return ref_time, mean_abs_stress


def case_label(case: BenchmarkCase) -> str:
    if case.group == "ep":
        return rf"$\phi = {case.friction_angle_deg:.0f}^\circ$"
    if case.group in ("steady_ab_pos", "steady_ab_neg"):
        return rf"$V_0 = {pow10_label(case.characteristic_velocity)}$"
    return rf"$D_c = {pow10_label(case.characteristic_distance)}$"


def case_color(case: BenchmarkCase) -> str:
    palette = {
        "ep_phi30": "#1f77b4",
        "ep_phi20": "#ff7f0e",
        "ep_phi10": "#2ca02c",
        "steady_ab_pos_v0_1e-6": "#1f77b4",
        "steady_ab_pos_v0_1e-5": "#ff7f0e",
        "steady_ab_pos_v0_1e-4": "#2ca02c",
        "steady_ab_neg_v0_1e-6": "#1f77b4",
        "steady_ab_neg_v0_1e-5": "#ff7f0e",
        "steady_ab_neg_v0_1e-4": "#2ca02c",
        "aging_ab_neg_dc_1e-6": "#1f77b4",
        "aging_ab_neg_dc_1e-5": "#ff7f0e",
        "aging_ab_neg_dc_1e-4": "#2ca02c",
    }
    return palette[case.name]


def selected_cases(groups: Iterable[str], names: Iterable[str]) -> dict[str, list[BenchmarkCase]]:
    grouped = {group: [] for group in GROUP_ORDER}
    allowed = set(groups)
    allowed_names = set(names)
    for case in BENCHMARK_CASES:
        if case.group not in allowed:
            continue
        if allowed_names and case.name not in allowed_names:
            continue
        grouped[case.group].append(case)
    return grouped


def compare_series(num_time: np.ndarray, num_stress: np.ndarray, ana_time: np.ndarray, ana_stress: np.ndarray) -> float:
    ana_on_num = np.interp(num_time, ana_time, ana_stress)
    scale = max(np.max(np.abs(ana_on_num)), 1.0)
    return float(np.max(np.abs(num_stress - ana_on_num)) / scale)


def cleanup_pycache(root_dir: Path) -> None:
    for path in root_dir.rglob("__pycache__"):
        if path.is_dir():
            import shutil

            shutil.rmtree(path, ignore_errors=True)


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Plot the monitor-based simple-shear benchmark suite.")
    parser.add_argument(
        "--output-root",
        type=Path,
        default=script_dir / "runs",
        help="Directory containing generated case folders.",
    )
    parser.add_argument(
        "--groups",
        nargs="+",
        choices=GROUP_ORDER,
        default=list(GROUP_ORDER),
        help="Benchmark groups to include.",
    )
    parser.add_argument(
        "--cases",
        nargs="+",
        default=[],
        help="Optional explicit case-name filter.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=script_dir / "simple_shear_benchmark.png",
        help="Output figure path.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Fail if any loaded numerical curve exceeds the relative error tolerance.",
    )
    parser.add_argument(
        "--no-save",
        action="store_true",
        help="Do not save a PNG. Useful for CI/functional checks.",
    )
    parser.add_argument(
        "--max-relative-error",
        type=float,
        default=5e-3,
        help="Tolerance used with --check.",
    )
    args = parser.parse_args()

    grouped_cases = selected_cases(args.groups, args.cases)
    output_root = args.output_root.resolve()
    set_plot_style()

    titles = {
        "ep": "EP Benchmark",
        "steady_ab_pos": "Steady RSF, a-b > 0",
        "steady_ab_neg": "Steady RSF, a-b < 0",
        "aging_ab_neg": "Aging RSF, a-b < 0",
    }
    panel_tags = {
        "ep": "(a)",
        "steady_ab_pos": "(b)",
        "steady_ab_neg": "(c)",
        "aging_ab_neg": "(d)",
    }

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.5), constrained_layout=True)
    axes_by_group = dict(zip(GROUP_ORDER, axes.flatten()))

    checked_errors: list[tuple[str, float]] = []

    for group in GROUP_ORDER:
        if group not in grouped_cases or not grouped_cases[group]:
            continue

        ax = axes_by_group[group]
        style_axis(ax, panel_tags[group], titles[group])

        for case in grouped_cases[group]:
            case_dir = output_root / case.name
            if not case_dir.is_dir():
                print(f"[skip] missing case directory: {case_dir}")
                continue

            num_time, num_stress = load_monitor_stress(case_dir)
            ana_time, ana_stress = analytical_series(case, total_time=float(num_time[-1]))
            rel_error = compare_series(num_time, num_stress, ana_time, ana_stress)
            checked_errors.append((case.name, rel_error))

            color = case_color(case)
            label = case_label(case)
            ax.plot(
                ana_time,
                ana_stress * STRESS_SCALE,
                color=color,
                linestyle="-",
                linewidth=LINE_W,
                label=f"Analytic, {label}",
            )
            ax.plot(
                num_time,
                num_stress * STRESS_SCALE,
                linestyle="none",
                marker="o",
                markersize=MARKER_SIZE,
                markerfacecolor="none",
                markeredgewidth=MARKER_EDGE_W,
                color=color,
                label=f"DynEarthSol, {label}",
            )
            print(f"[ok] {case.name}: max relative error = {rel_error:.3e}")

        handles, labels = ax.get_legend_handles_labels()
        if handles:
            legend = ax.legend(handles, labels, loc="best")
            style_legend(legend)

    output_path = args.output.resolve()
    if not args.no_save:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300)
    plt.close(fig)
    if not args.no_save:
        print(f"[saved] {output_path}")

    if args.check:
        failed = [(name, err) for name, err in checked_errors if err > args.max_relative_error]
        if failed:
            details = ", ".join(f"{name}={err:.3e}" for name, err in failed)
            raise SystemExit(f"Benchmark check failed: {details}")


if __name__ == "__main__":
    _script_dir = Path(__file__).resolve().parent
    try:
        main()
    finally:
        cleanup_pycache(_script_dir)
