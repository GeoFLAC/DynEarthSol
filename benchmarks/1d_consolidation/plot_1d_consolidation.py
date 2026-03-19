#!/usr/bin/env python3

import argparse
import os
import sys
from pathlib import Path

import matplotlib

if "--output" in sys.argv:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]

if str(REPO_ROOT) not in sys.path:
    sys.path.append(str(REPO_ROOT))
if str(SCRIPT_DIR) not in sys.path:
    sys.path.append(str(SCRIPT_DIR))

from Dynearthsol import Dynearthsol
from benchmark_cases import CASE_CONFIG, CASE_ORDER


ATMOSPHERIC_PRESSURE = 0.0
GRAVITY = 10.0
FLUID_DENSITY = 1000.0
GAMMA_W = GRAVITY * FLUID_DENSITY
NOTEBOOK_NUMERICAL_COLOR = "#c0392b"
NOTEBOOK_ANALYTICAL_COLOR = "#111111"
CASE_COLORS = ("#c0392b", "#1f77b4", "#2d8a4e", "#c47f1f")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot the 1D consolidation benchmark against the analytical solution."
    )
    parser.add_argument(
        "--case",
        nargs="+",
        choices=CASE_ORDER,
        default=None,
        help="Named benchmark cases to plot using the built-in metadata.",
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        default=Path("."),
        help="Directory containing <model>.info and field outputs. Defaults to the current directory.",
    )
    parser.add_argument(
        "--model",
        nargs="+",
        default=None,
        help="Model name prefixes to read (for example terzaghi or output/terzaghi).",
    )
    parser.add_argument(
        "--dim",
        type=int,
        default=None,
        choices=(2, 3),
        help="Problem dimension used to select the vertical coordinate component. If omitted, named cases use their metadata.",
    )
    parser.add_argument(
        "--label",
        nargs="+",
        default=None,
        help="Optional legend labels matching the --model entries.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output image path. If omitted, the figure is shown interactively.",
    )
    parser.add_argument(
        "--skip-first-output",
        type=int,
        default=None,
        help="Override the number of initial outputs excluded from the error panel and metrics.",
    )
    return parser.parse_args()


def analytical_pressure(time_points):
    perm = 1.02e-12
    eta_f = 1.002e-3
    phi = 0.3
    beta_w = 1.0 / 2.17e9
    k_bulk = 1.0e7
    mu = 1.0e7
    alpha_c = 1.0
    height = 10.0
    sigma_zz = -1.0e5
    num_terms = 1000
    z = 0.0

    hydraulic_conductivity = perm * GAMMA_W / eta_f
    lambda_lame = k_bulk - (2.0 / 3.0) * mu
    cv = (
        phi * beta_w + alpha_c * (alpha_c + phi - phi * alpha_c) / (lambda_lame + 2.0 * mu)
    ) ** -1 * (hydraulic_conductivity / GAMMA_W)

    p_excess_time = np.zeros(len(time_points), dtype=float)
    for i, t in enumerate(time_points):
        p_excess = 0.0
        for j in range(num_terms):
            term = (4.0 / np.pi) * ((-1.0) ** j) / (2 * j + 1)
            cos_term = np.cos((2 * j + 1) * np.pi * z / (2.0 * height))
            exp_term = np.exp(-((2 * j + 1) * np.pi / 2.0) ** 2 * cv * t / height**2)
            p_excess += term * cos_term * exp_term
        p_excess_time[i] = (-sigma_zz) * p_excess

    return p_excess_time, height, cv


def read_numerical_pore_pressure(model_name, dim):
    des = Dynearthsol(model_name)
    num_steps = len(des.time)
    pore_pressures = np.zeros(num_steps, dtype=float)
    coord_axis = 2 if dim == 3 else 1
    for time_step in range(num_steps):
        pore_pressure = des.read_field(time_step, "pore pressure")
        coordinates = des.read_field(time_step, "coordinate")
        pore_pressures[time_step] = (
            pore_pressure[0] - ATMOSPHERIC_PRESSURE + GAMMA_W * coordinates[0, coord_axis]
        )
    return np.array(des.time), pore_pressures


def compute_relative_error(numerical_pressure, analytical_pressure_values):
    reference = np.maximum(np.abs(analytical_pressure_values), 1.0e-30)
    return np.abs(numerical_pressure - analytical_pressure_values) / reference


def format_percent(value, _):
    return f"{value * 100:.1f}%"


def resolve_case_specs(args):
    if args.case:
        case_specs = []
        for case_name in args.case:
            config = CASE_CONFIG[case_name]
            case_specs.append(
                {
                    "case_name": case_name,
                    "model": config["model"],
                    "label": config["label"],
                    "dim": config["dim"],
                    "skip_first_output": (
                        args.skip_first_output
                        if args.skip_first_output is not None
                        else config["skip_first_output"]
                    ),
                }
            )
        return case_specs

    if args.model:
        labels = args.label or [Path(model).name.removeprefix("terzaghi_").replace("_", " ") for model in args.model]
        skip = 0 if args.skip_first_output is None else args.skip_first_output
        dim = 2 if args.dim is None else args.dim
        return [
            {
                "case_name": Path(model).name,
                "model": model,
                "label": label,
                "dim": dim,
                "skip_first_output": skip,
            }
            for model, label in zip(args.model, labels)
        ]

    info_files = sorted(args.run_dir.glob("*.info"))
    if not info_files:
        raise FileNotFoundError(
            "Could not infer any models. Provide --case, --model, or run from a directory with *.info files."
        )
    return [
        {
            "case_name": info_file.stem,
            "model": info_file.stem,
            "label": info_file.stem.removeprefix("terzaghi_").replace("_", " "),
            "dim": 2 if args.dim is None else args.dim,
            "skip_first_output": 0 if args.skip_first_output is None else args.skip_first_output,
        }
        for info_file in info_files
    ]


def resolve_plot_dim(case_specs, dim_override):
    if dim_override is not None:
        return dim_override

    dims = {case_spec.get("dim", 2) for case_spec in case_specs}
    if len(dims) == 1:
        return dims.pop()

    raise ValueError("Requested cases mix 2D and 3D data. Plot them separately or pass --dim explicitly.")


def load_case_results(case_specs, run_dir, dim):
    cases = []
    max_time = 0.0
    height = None
    cv = None

    cwd = Path.cwd()
    try:
        os.chdir(run_dir)
        for case_spec in case_specs:
            model_name = case_spec["model"]
            info_path = Path(f"{model_name}.info")
            if not info_path.exists():
                raise FileNotFoundError(f"Could not find info file: {run_dir / info_path}")

            times, pore_pressure = read_numerical_pore_pressure(model_name, dim)
            analytical_at_outputs, case_height, case_cv = analytical_pressure(times)
            if height is None:
                height = case_height
                cv = case_cv
            max_time = max(max_time, float(times[-1]))

            compare_start = min(case_spec["skip_first_output"], len(times))
            relative_error = compute_relative_error(pore_pressure, analytical_at_outputs)
            cases.append(
                {
                    "case_name": case_spec["case_name"],
                    "model_name": model_name,
                    "label": case_spec["label"],
                    "times": times,
                    "pore_pressure": pore_pressure,
                    "analytical_pressure": analytical_at_outputs,
                    "relative_error": relative_error,
                    "compare_start": compare_start,
                }
            )
    finally:
        os.chdir(cwd)

    return cases, max_time, height, cv


def build_metrics(cases):
    metrics = []
    for case in cases:
        compare_start = case["compare_start"]
        compared_error = case["relative_error"][compare_start:]
        first_pressure = np.nan
        first_error = np.nan
        max_error = np.nan
        final_error = np.nan
        if compare_start < len(case["pore_pressure"]):
            first_pressure = case["pore_pressure"][compare_start] * 1.0e-5
            if len(compared_error) > 0:
                first_error = compared_error[0]
                max_error = np.nanmax(compared_error)
                final_error = compared_error[-1]

        metrics.append(
            {
                "case_name": case["case_name"],
                "label": case["label"],
                "first_compared_pressure_bar": first_pressure,
                "first_error": first_error,
                "max_error": max_error,
                "final_error": final_error,
            }
        )
    return metrics


def setup_notebook_style():
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "font.size": 12,
            "axes.labelsize": 15,
            "axes.titlesize": 18,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "legend.fontsize": 12,
        }
    )


def plot_results(case_specs, dim, output_path=None, run_dir=Path(".")):
    setup_notebook_style()
    run_dir = run_dir.resolve()
    cases, max_time, height, cv = load_case_results(case_specs, run_dir, dim)
    metrics = build_metrics(cases)

    analytical_time = np.linspace(
        0.0,
        max_time,
        max(2500, 30 * max(len(case["times"]) for case in cases)),
    )
    analytical_curve, _, _ = analytical_pressure(analytical_time)

    fig, (ax_main, ax_err) = plt.subplots(
        2,
        1,
        figsize=(12, 8),
        sharex=True,
        gridspec_kw={"height_ratios": (4, 1.2), "hspace": 0.06},
        constrained_layout=True,
    )

    analytical_label = "Analytical solution"
    for index, case in enumerate(cases):
        color = NOTEBOOK_NUMERICAL_COLOR if len(cases) == 1 else CASE_COLORS[index % len(CASE_COLORS)]
        markevery = max(1, len(case["times"]) // 18)
        compare_start = case["compare_start"]
        numerical_label = "DynEarthSol" if len(cases) == 1 else case["label"]

        ax_main.plot(
            case["times"],
            case["pore_pressure"] * 1.0e-5,
            linestyle="--",
            linewidth=2.0,
            color=color,
            marker="o",
            markersize=6.0,
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=0.8,
            markevery=markevery,
            label=numerical_label,
            zorder=3,
        )
        ax_err.plot(
            case["times"][compare_start:],
            case["relative_error"][compare_start:],
            color=color,
            linewidth=2.0,
            label=numerical_label,
        )

        ax_main.plot(
            analytical_time,
            analytical_curve * 1.0e-5,
            color=NOTEBOOK_ANALYTICAL_COLOR,
            linestyle="-",
            linewidth=2.6,
            label=analytical_label if index == 0 else None,
            zorder=1,
        )

    for ax in (ax_main, ax_err):
        ax.grid(True, linestyle="--", linewidth=0.9, alpha=0.55)
        for spine in ax.spines.values():
            spine.set_linewidth(2.2)

    ax_main.set_title(f"1D Consolidation at z = {height:.0f} m", pad=12)
    ax_main.set_ylabel("Pore Pressure (bar)")
    ax_err.set_ylabel("Rel. Error")
    ax_err.set_xlabel("Time (s)")
    ax_err.set_yscale("log")
    ax_err.yaxis.set_major_formatter(FuncFormatter(format_percent))
    ax_main.legend(loc="upper right", frameon=True)

    print(f"Consolidation coefficient Cv = {cv:.6e}")
    for metric in metrics:
        print(
            f"{metric['label']}: first compared pore pressure = "
            f"{metric['first_compared_pressure_bar']:.6f} bar, "
            f"first/max/final relative error = "
            f"{metric['first_error'] * 100:.3f}% / "
            f"{metric['max_error'] * 100:.3f}% / "
            f"{metric['final_error'] * 100:.3f}%"
        )

    if output_path is None:
        plt.show()
    else:
        output_path = output_path.resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=220, bbox_inches="tight")
        print(f"Saved figure to {output_path}")

    return metrics


def main():
    args = parse_args()
    case_specs = resolve_case_specs(args)
    plot_dim = resolve_plot_dim(case_specs, args.dim)
    plot_results(case_specs, plot_dim, output_path=args.output, run_dir=args.run_dir)


if __name__ == "__main__":
    main()
