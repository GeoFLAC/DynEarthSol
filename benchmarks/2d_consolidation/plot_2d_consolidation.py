#!/usr/bin/env python3

import argparse
import os
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np

if str(REPO_ROOT) not in sys.path:
    sys.path.append(str(REPO_ROOT))

from Dynearthsol import Dynearthsol
from benchmark_cases import CASE_CONFIG, CASE_ORDER


PI = np.pi
EPS = 1e-4
N_ROOTS = 500

# Keep these constants aligned with the notebook currently used for Mandel.
ANALYTICAL_GRAVITY = 9.81
NUMERICAL_GRAVITY = 10.0
FLUID_DENSITY = 1000.0
ATMOSPHERIC_PRESSURE = 0.0


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot the 2D Mandel consolidation benchmark against the analytical solution."
    )
    parser.add_argument(
        "--case",
        nargs="+",
        choices=list(CASE_ORDER),
        default=["mandel"],
        help="Named benchmark cases to plot.",
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        default=None,
        help="Directory containing <model>.info and field outputs.",
    )
    parser.add_argument(
        "--model",
        type=str,
        default=None,
        help="Model name prefix to read. Defaults to the case metadata.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output image path. Defaults to <run-dir>/final_comparison.png.",
    )
    return parser.parse_args()


def solve_mandel_roots(eta, n_roots=N_ROOTS):
    roots = np.zeros(n_roots)
    for i in range(1, n_roots + 1):
        e1 = (i - 1) * PI + PI / 4.0
        e2 = e1 + PI / 2.0 - EPS
        em = e1

        for _ in range(500):
            y1 = np.tan(e1) - eta * e1
            y2 = np.tan(e2) - eta * e2
            em = 0.5 * (e1 + e2)
            ym = np.tan(em) - eta * em
            if ym * y1 > 0:
                e1 = em
            else:
                e2 = em
            if abs(y2) < EPS:
                em = e2

        roots[i - 1] = em

    return roots


def analytical_pressure(time_points):
    perm = 1.02e-12
    eta_f = 1.002e-3
    phi = 0.3
    k_bulk = 1e7
    shear_modulus = 1e7
    k_s = 1e38
    force = 1e5
    width = 5.0
    monitor_x = 0.0

    gamma_w = ANALYTICAL_GRAVITY * FLUID_DENSITY
    youngs_modulus = (9 * k_bulk * shear_modulus) / (3 * k_bulk + shear_modulus)
    poisson_ratio = (3 * k_bulk - 2 * shear_modulus) / (2 * (3 * k_bulk + shear_modulus))

    am11 = youngs_modulus * (1 - poisson_ratio) / (1 - poisson_ratio - 2 * poisson_ratio * poisson_ratio)
    am33 = am11
    am13 = youngs_modulus * poisson_ratio / (1 - poisson_ratio - 2 * poisson_ratio * poisson_ratio)

    # The notebook assumes an effectively incompressible solid grain.
    _ = 1 - k_bulk / k_s
    bm = 1 / (4.5e-10 * phi)
    cv = (perm * bm * am11) / (eta_f * (am11 + bm))
    a11 = (am33 - 2 * am13 + am11) / (am11 - am13)
    a1 = a11 + (am11 * am33 - am13 * am13) / (bm * (am11 - am13))
    a2 = (am11 - am13) / am11
    eta = a1 / a2

    roots = solve_mandel_roots(eta)
    pressure = np.zeros_like(time_points, dtype=float)
    for j, time_value in enumerate(time_points):
        value = 0.0
        for root in roots:
            sin_root = np.sin(root)
            cos_root = np.cos(root)
            cos_x = np.cos(monitor_x * root / width)
            decay = np.exp(-root * root * cv * time_value / (width * width))
            value += sin_root * (cos_x - cos_root) * decay / (root - sin_root * cos_root)
        pressure[j] = 2 * force * value / a1

    return pressure


def read_time_steps(info_path):
    data = np.loadtxt(info_path, usecols=(0, 2))
    return data[:, 0].astype(int), data[:, 1]


def numerical_pressure(run_dir, model_name, vertical_axis=1):
    _, times = read_time_steps(run_dir / f"{model_name}.info")
    gamma_w = NUMERICAL_GRAVITY * FLUID_DENSITY
    des = Dynearthsol(str(run_dir / model_name))
    pressure = np.zeros(len(times), dtype=float)
    for step in range(len(times)):
        pore_pressure = des.read_field(step, "pore pressure")
        coordinates = des.read_field(step, "coordinate")
        pressure[step] = (
            pore_pressure[0]
            - ATMOSPHERIC_PRESSURE
            + gamma_w * coordinates[0, vertical_axis]
        )
    return times, pressure


def analytical_plot_time_points(times, n_points=500):
    positive_times = np.asarray(times)[np.asarray(times) > 0]
    if positive_times.size == 0:
        return np.asarray(times, dtype=float)

    t_min = positive_times.min()
    t_max = positive_times.max()
    if np.isclose(t_min, t_max):
        return positive_times

    return np.geomspace(t_min, t_max, n_points)


def style_axis(axis):
    axis.set_xscale("log")
    axis.set_xlabel("Time (sec)", fontsize=18)
    axis.set_ylabel("Pore Pressure (Bar)", fontsize=18)
    axis.grid(True, linestyle="--", linewidth=1.0, alpha=0.55, color="#9aa0a6")
    axis.tick_params(axis="both", labelsize=15, width=1.8, length=6)
    for spine in axis.spines.values():
        spine.set_linewidth(2.2)


def plot_results(case_name, run_dir, model_name, output_path):
    times, numerical = numerical_pressure(run_dir, model_name)
    analytical = analytical_pressure(times)
    absolute_error_bar = np.abs(numerical - analytical) * 1e-5

    plot_times = analytical_plot_time_points(times)
    analytical_smooth = analytical_pressure(plot_times)

    first_error = absolute_error_bar[1] if len(absolute_error_bar) > 1 else absolute_error_bar[0]
    metrics = {
        "first_error_bar": float(first_error),
        "max_error_bar": float(np.nanmax(absolute_error_bar)),
        "final_error_bar": float(absolute_error_bar[-1]),
    }

    fig, axis = plt.subplots(figsize=(12, 8))
    axis.plot(
        plot_times,
        analytical_smooth * 1e-5,
        color="black",
        linewidth=2.4,
        label="Analytical solution",
        zorder=2,
    )
    axis.plot(
        times,
        numerical * 1e-5,
        linestyle="-",
        linewidth=1.3,
        color="#e76f51",
        alpha=0.75,
        zorder=3,
    )
    axis.plot(
        times,
        numerical * 1e-5,
        "o",
        color="#d62828",
        markersize=7.0,
        markeredgewidth=0.0,
        label="DynEarthSol",
        zorder=4,
    )

    style_axis(axis)
    axis.legend(loc="upper right", fontsize=14, frameon=True, fancybox=False, framealpha=0.95)
    axis.set_title(CASE_CONFIG[case_name]["label"], fontsize=18, pad=10)

    metrics_text = (
        f"Absolute error (bar)\n"
        f"first: {metrics['first_error_bar']:.3f}\n"
        f"max:  {metrics['max_error_bar']:.3f}\n"
        f"final:{metrics['final_error_bar']:.3f}"
    )
    axis.text(
        0.035,
        0.08,
        metrics_text,
        transform=axis.transAxes,
        fontsize=12.5,
        va="bottom",
        ha="left",
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "white",
            "edgecolor": "#b0b0b0",
            "alpha": 0.92,
        },
    )

    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return metrics


def main():
    args = parse_args()
    case_name = args.case[0]
    config = CASE_CONFIG[case_name]
    run_dir = args.run_dir.resolve() if args.run_dir is not None else (SCRIPT_DIR / "runs" / config["run_subdir"]).resolve()
    model_name = args.model if args.model is not None else config["model"]
    output_path = args.output.resolve() if args.output is not None else run_dir / "final_comparison.png"

    metrics = plot_results(case_name, run_dir, model_name, output_path)
    print(
        f"{config['label']}: first/max/final absolute error = "
        f"{metrics['first_error_bar']:.3f} / "
        f"{metrics['max_error_bar']:.3f} / "
        f"{metrics['final_error_bar']:.3f} bar"
    )
    print(f"Saved figure to {output_path}")


if __name__ == "__main__":
    main()
