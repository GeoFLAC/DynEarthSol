#!/usr/bin/env python3
"""Plot the paper's local EP/RSF benchmark and its pointwise errors."""

from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import sys
from pathlib import Path
from typing import Iterable

sys.dont_write_bytecode = True

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.ticker import LogLocator, NullFormatter

from benchmark_reference import (
    analytical_state_ratio,
    analytical_stress,
    load_monitor_case,
    pointwise_error_metrics,
    slip_rate,
)
from run_simple_shear_benchmark import (
    BENCHMARK_CASES,
    GROUP_ORDER,
    BenchmarkCase,
)


STRESS_SCALE = 1.0e-6
GROUP_TITLE = {
    "ep": "Elastoplastic",
    "steady_ab_neg": r"Steady-state RSF ($a-b<0$)",
    "aging_ab_neg": r"Aging-law RSF ($a-b<0$)",
}


def set_plot_style() -> None:
    try:
        from matplotlib import font_manager

        for font_path in (
            "/usr/share/fonts/opentype/urw-base35/NimbusSans-Regular.otf",
            "/usr/share/fonts/opentype/urw-base35/NimbusSans-Bold.otf",
            "/usr/share/fonts/opentype/urw-base35/NimbusSans-Italic.otf",
        ):
            if os.path.isfile(font_path):
                font_manager.fontManager.addfont(font_path)
    except Exception:
        print(
            "Warning: failed to configure Nimbus Sans fonts; using defaults.",
            file=sys.stderr,
        )

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Nimbus Sans",
                "Arial",
                "Liberation Sans",
                "DejaVu Sans",
            ],
            "font.size": 9.0,
            "axes.labelsize": 10.0,
            "axes.titlesize": 10.5,
            "legend.fontsize": 7.6,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "axes.linewidth": 1.0,
            "lines.linewidth": 1.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def case_color(index: int) -> str:
    return ("#0072B2", "#D55E00", "#009E73")[index]


def case_linestyle(index: int) -> str:
    return ("-", "--", ":")[index]


def case_marker(index: int) -> str:
    return ("o", "s", "^")[index]


def parameter_text(case: BenchmarkCase) -> str:
    if case.group == "ep":
        return rf"$\phi={case.friction_angle_deg:.0f}^\circ$"
    if case.group == "steady_ab_neg":
        exponent = int(round(math.log10(case.characteristic_velocity)))
        return rf"$V_0=10^{{{exponent}}}$ m s$^{{-1}}$"
    exponent = int(math.floor(math.log10(case.characteristic_distance)))
    coefficient = case.characteristic_distance / 10.0**exponent
    if abs(coefficient - 1.0) < 1.0e-9:
        return rf"$D_c=10^{{{exponent}}}$ m"
    return rf"$D_c={coefficient:g}\times10^{{{exponent}}}$ m"


def selected_cases(
    groups: Iterable[str],
    names: Iterable[str],
) -> dict[str, list[BenchmarkCase]]:
    grouped = {group: [] for group in GROUP_ORDER}
    allowed_groups = set(groups)
    allowed_names = set(names)
    known_names = {case.name for case in BENCHMARK_CASES}
    missing = sorted(allowed_names - known_names)
    if missing:
        raise ValueError(f"Unknown benchmark case(s): {', '.join(missing)}")
    for case in BENCHMARK_CASES:
        if case.group not in allowed_groups:
            continue
        if allowed_names and case.name not in allowed_names:
            continue
        grouped[case.group].append(case)
    return grouped


def draw_schematic(axis: plt.Axes) -> None:
    points = np.array(
        ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
    )
    triangles = np.array(((0, 1, 3), (1, 2, 3)), dtype=int)
    mesh = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
    axis.triplot(mesh, color="black", linewidth=1.2)
    axis.scatter(points[:, 0], points[:, 1], s=22, color="black", zorder=3)

    green = "#009E73"
    axis.annotate(
        "",
        xy=(0.5, 0.5),
        xytext=(0.0, 0.0),
        arrowprops={
            "arrowstyle": "<->",
            "color": green,
            "linewidth": 1.1,
            "shrinkA": 1.5,
            "shrinkB": 0.0,
        },
    )
    foot = np.array((0.5, 0.5))
    along = np.array((1.0, -1.0)) / math.sqrt(2.0)
    inward = np.array((-1.0, -1.0)) / math.sqrt(2.0)
    first = foot + 0.065 * along
    second = foot + 0.065 * inward
    right_angle = np.vstack((first, first + second - foot, second))
    axis.plot(right_angle[:, 0], right_angle[:, 1], color=green, linewidth=0.8)
    axis.text(
        0.32,
        0.19,
        r"$w$",
        color=green,
        fontsize=10.5,
        ha="center",
        va="center",
    )
    axis.text(
        0.5,
        -0.22,
        r"$\dot\varepsilon_{\mathrm{II}}=v_x/2H$ in both elements",
        color="#555555",
        fontsize=8.6,
        ha="center",
        va="center",
    )

    for x_value in (0.0, 1.0):
        axis.plot(
            (x_value, x_value - 0.05, x_value + 0.05, x_value),
            (0.0, -0.10, -0.10, 0.0),
            color="black",
            linewidth=0.9,
        )
        axis.arrow(
            x_value,
            1.0,
            0.19,
            0.0,
            head_width=0.035,
            head_length=0.05,
            length_includes_head=True,
            color="#0072B2",
            linewidth=1.0,
        )

    axis.text(
        0.5,
        1.07,
        r"$v_x=10^{-5}$ m s$^{-1}$",
        ha="center",
        va="bottom",
        fontsize=10,
    )
    axis.set_xlim(-0.12, 1.22)
    axis.set_ylim(-0.30, 1.16)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("x (m)")
    axis.set_ylabel("y (m)")
    axis.set_xticks((0.0, 0.5, 1.0))
    axis.set_yticks((0.0, 0.5, 1.0))
    axis.tick_params(width=1.0, length=4)
    axis.set_title("Simple-shear benchmark", pad=11)


def style_stress_axis(axis: plt.Axes, title: str) -> None:
    axis.set_title(title, pad=5)
    axis.set_ylabel(r"$|\sigma_{xy}|$ (MPa)")
    axis.set_xlim(0.0, 2000.0)
    axis.tick_params(width=1.0, length=4, labelbottom=False)


def style_error_axis(axis: plt.Axes) -> None:
    axis.set_xlim(0.0, 2000.0)
    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Rel. error\n(%)", labelpad=3)
    axis.set_yscale("log")
    axis.yaxis.set_major_locator(LogLocator(base=10, numticks=4))
    axis.yaxis.set_minor_formatter(NullFormatter())
    axis.tick_params(width=0.9, length=3.5)


def add_panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.17,
        1.08,
        label,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=12,
        fontweight="bold",
    )


def cleanup_pycache(root_dir: Path) -> None:
    for path in root_dir.rglob("__pycache__"):
        if path.is_dir():
            shutil.rmtree(path, ignore_errors=True)


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Plot the paper's monitor-based simple-shear benchmark."
    )
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
        help="Fail if any loaded stress curve exceeds the error tolerance.",
    )
    parser.add_argument(
        "--no-save",
        action="store_true",
        help="Do not save a figure or metrics CSV.",
    )
    parser.add_argument(
        "--max-relative-error",
        type=float,
        default=2.0e-5,
        help="Fractional pointwise error tolerance used with --check.",
    )
    args = parser.parse_args()

    grouped_cases = selected_cases(args.groups, args.cases)
    output_root = args.output_root.resolve()
    set_plot_style()

    figure = plt.figure(figsize=(7.4, 7.0))
    outer = figure.add_gridspec(
        2,
        2,
        hspace=0.30,
        wspace=0.30,
        left=0.095,
        right=0.965,
        top=0.945,
        bottom=0.065,
    )
    schematic_axis = figure.add_subplot(outer[0, 0])
    draw_schematic(schematic_axis)
    add_panel_label(schematic_axis, "a")

    slots = {
        "ep": outer[0, 1],
        "steady_ab_neg": outer[1, 0],
        "aging_ab_neg": outer[1, 1],
    }
    panel_labels = {
        "ep": "b",
        "steady_ab_neg": "c",
        "aging_ab_neg": "d",
    }
    metrics_rows: list[dict[str, str | int | float]] = []
    failures: list[str] = []

    for group in GROUP_ORDER:
        inner = slots[group].subgridspec(
            2,
            1,
            height_ratios=[2.15, 1.0],
            hspace=0.10,
        )
        stress_axis = figure.add_subplot(inner[0])
        error_axis = figure.add_subplot(inner[1], sharex=stress_axis)
        add_panel_label(stress_axis, panel_labels[group])

        for index, case in enumerate(grouped_cases[group]):
            case_dir = output_root / case.name
            if not case_dir.is_dir():
                message = f"missing case directory: {case_dir}"
                print(f"[skip] {message}", file=sys.stderr)
                if args.check:
                    failures.append(message)
                continue

            data = load_monitor_case(case_dir)
            numerical = data.mean_abs_stress
            reference = analytical_stress(case, data.time_s)
            metrics = pointwise_error_metrics(numerical, reference)
            errors_pct = [
                100.0 * abs(num - ref) / abs(ref)
                for num, ref in zip(numerical[1:], reference[1:])
            ]

            color = case_color(index)
            linestyle = case_linestyle(index)
            marker = case_marker(index)
            stress_axis.plot(
                data.time_s,
                np.asarray(reference) * STRESS_SCALE,
                color=color,
                linestyle=linestyle,
                label=parameter_text(case),
            )
            marker_slice = slice(index, None, 3)
            stress_axis.plot(
                np.asarray(data.time_s)[marker_slice],
                np.asarray(numerical)[marker_slice] * STRESS_SCALE,
                linestyle="none",
                marker=marker,
                markersize=4.2,
                markerfacecolor="none",
                markeredgewidth=1.0,
                color=color,
            )
            error_axis.plot(
                data.time_s[1:],
                errors_pct,
                color=color,
                linestyle=linestyle,
                linewidth=1.2,
            )

            element_spread = max(
                abs(left - right)
                for left, right in zip(
                    data.stress_by_element[0],
                    data.stress_by_element[1],
                )
            )
            stress_scale = max(
                max(values) for values in data.stress_by_element
            )
            metrics_rows.append(
                {
                    "case": case.name,
                    "group": group,
                    "parameter": parameter_text(case).replace("$", ""),
                    "samples": metrics.samples,
                    "mean_pct": 100.0 * metrics.mean_fraction,
                    "max_pct": 100.0 * metrics.max_fraction,
                    "final_pct": 100.0 * metrics.final_fraction,
                    "element_spread_pct": (
                        100.0 * element_spread / max(stress_scale, 1.0)
                    ),
                }
            )
            print(
                f"[ok] {case.name}: "
                f"mean={100.0 * metrics.mean_fraction:.3e}%, "
                f"max={100.0 * metrics.max_fraction:.3e}%"
            )
            if args.check and metrics.max_fraction > args.max_relative_error:
                failures.append(
                    f"{case.name}={metrics.max_fraction:.3e}"
                )

        style_stress_axis(stress_axis, GROUP_TITLE[group])
        style_error_axis(error_axis)
        handles, labels = stress_axis.get_legend_handles_labels()
        if handles:
            stress_axis.legend(
                handles,
                labels,
                loc="lower right",
                frameon=False,
                handlelength=2.4,
                borderaxespad=0.3,
            )

        if group == "aging_ab_neg":
            inset = stress_axis.inset_axes((0.520, 0.485, 0.445, 0.375))
            for index, case in enumerate(grouped_cases[group]):
                transient_dir = output_root / f"{case.name}_transient"
                if not transient_dir.is_dir():
                    continue
                transient = load_monitor_case(transient_dir)
                state = transient.state_by_element[1]
                if state is None:
                    continue
                theta_ss = [
                    case.characteristic_distance / slip_rate(value)
                    for value in transient.time_s
                ]
                numerical_ratio = [
                    value / steady
                    for value, steady in zip(state, theta_ss)
                ]
                reference_ratio = analytical_state_ratio(
                    case,
                    transient.time_s,
                )
                color = case_color(index)
                linestyle = case_linestyle(index)
                marker = case_marker(index)
                inset.plot(
                    transient.time_s[1:],
                    reference_ratio[1:],
                    color=color,
                    linestyle=linestyle,
                    linewidth=1.3,
                )
                requested_times = np.geomspace(
                    1.0,
                    float(transient.time_s[-1]),
                    11,
                ) * (1.0 + 0.09 * index)
                indices = np.unique(
                    np.clip(
                        np.searchsorted(transient.time_s, requested_times),
                        1,
                        len(transient.time_s) - 1,
                    )
                )
                inset.plot(
                    np.asarray(transient.time_s)[indices],
                    np.asarray(numerical_ratio)[indices],
                    linestyle="none",
                    marker=marker,
                    markersize=3.6,
                    markerfacecolor="none",
                    markeredgewidth=0.9,
                    color=color,
                )
            inset.set_xscale("log")
            inset.set_xlim(1.0, 2000.0)
            inset.set_xlabel("Time (s)", fontsize=7.2, labelpad=1.5)
            inset.set_ylabel(
                r"$\theta/\theta_{\mathrm{ss}}$",
                fontsize=7.6,
                labelpad=1.5,
            )
            inset.set_yticks((0.7, 0.8, 0.9, 1.0))
            inset.tick_params(labelsize=6.6, width=0.8, length=2.8)

    output_path = args.output.resolve()
    if not args.no_save:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output_path, dpi=300)
        metrics_path = output_path.with_name(
            f"{output_path.stem}_metrics.csv"
        )
        if metrics_rows:
            with metrics_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=list(metrics_rows[0].keys()),
                )
                writer.writeheader()
                writer.writerows(metrics_rows)
        print(f"[saved] {output_path}")
        if metrics_rows:
            print(f"[saved] {metrics_path}")
    plt.close(figure)

    if failures:
        raise SystemExit("Benchmark check failed: " + ", ".join(failures))


if __name__ == "__main__":
    _script_dir = Path(__file__).resolve().parent
    try:
        main()
    finally:
        cleanup_pycache(_script_dir)
