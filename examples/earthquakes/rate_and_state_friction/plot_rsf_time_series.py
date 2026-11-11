#!/usr/bin/env python3
"""Create a two-panel RSF verification plot from monitor CSV files.

Panels:
1) Time vs velocity magnitude (log scale on y)
2) Time vs displacement magnitude

Baseline series is intentionally excluded.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np

YEAR_IN_SECONDS = 365.25 * 24.0 * 3600.0
VEL_LOG_MIN = 1e-12


def load_series(csv_path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    arr = np.genfromtxt(
        csv_path,
        delimiter=",",
        names=True,
        dtype=None,
        encoding="utf-8",
    )
    arr = np.atleast_1d(arr)
    names = arr.dtype.names or ()

    if "time_s" not in names:
        raise ValueError(f"Missing column 'time_s' in {csv_path}")

    vel_cols = [c for c in names if c.startswith("velocity_")]
    coord_cols = [c for c in names if c.startswith("coord_")]
    if not vel_cols:
        raise ValueError(f"No velocity_* columns in {csv_path}")
    if not coord_cols:
        raise ValueError(f"No coord_* columns in {csv_path}")

    time_s = np.asarray(arr["time_s"], dtype=float)
    order = np.argsort(time_s)

    time_yr = time_s[order] / YEAR_IN_SECONDS

    vel = np.column_stack([np.asarray(arr[c], dtype=float) for c in vel_cols])[order]
    vel_mag = np.sqrt(np.sum(vel * vel, axis=1))
    vel_mag_plot = np.maximum(vel_mag, VEL_LOG_MIN)

    coord = np.column_stack([np.asarray(arr[c], dtype=float) for c in coord_cols])[order]
    disp_mag = np.sqrt(np.sum((coord - coord[0]) ** 2, axis=1))
    return time_yr, vel_mag_plot, disp_mag


def main() -> None:
    script_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(description="Plot RSF monitor series (velocity/displacement).")
    parser.add_argument(
        "--creep",
        type=Path,
        default=script_dir / "ep_rsf_creep" / "monitor_point_0.csv",
        help="Monitor CSV for creep case.",
    )
    parser.add_argument(
        "--eq",
        type=Path,
        default=script_dir / "ep_rsf_eq" / "monitor_point_0.csv",
        help="Monitor CSV for earthquake-cycle case.",
    )
    parser.add_argument(
        "--creep-large",
        dest="creep_large",
        type=Path,
        default=script_dir / "ep_rsf_creep_largeDc" / "monitor_point_0.csv",
        help="Monitor CSV for large-Dc creep case.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=str(script_dir / "rsf_time_series.png"),
        help="Output figure path.",
    )
    args = parser.parse_args()

    series_paths = {
        "Creep (a-b>0)": Path(args.creep),
        "EQ cycle (a-b<0)": Path(args.eq),
        "Creep, large Dc": Path(args.creep_large),
    }
    styles = {
        "Creep (a-b>0)": dict(color="#1f9d55", linestyle="--", linewidth=2.0),
        "EQ cycle (a-b<0)": dict(color="#d1495b", linestyle="-.", linewidth=2.0),
        "Creep, large Dc": dict(color="#e68a00", linestyle=(0, (6, 2)), linewidth=2.0),
    }

    loaded: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for label, path in series_paths.items():
        if not path.exists():
            print(f"[skip] missing CSV: {path}")
            continue
        loaded[label] = load_series(path)
        print(f"[ok] {label}: {path}")

    if not loaded:
        raise FileNotFoundError("No input monitor CSV found.")

    fig, (ax_vel, ax_disp) = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)

    for label, (time_yr, vel_mag, disp_mag) in loaded.items():
        ax_vel.plot(time_yr, vel_mag, label=label, **styles[label])
        ax_disp.plot(time_yr, disp_mag, label=label, **styles[label])

    ax_vel.set_yscale("log")
    ax_vel.set_xlabel("Time (yr)")
    ax_vel.set_ylabel("Velocity magnitude (m/s)")
    ax_vel.grid(True, linestyle="--", linewidth=0.6, alpha=0.6)
    ax_vel.text(0.02, 0.96, "(a)", transform=ax_vel.transAxes, va="top", ha="left", fontsize=11, fontweight="bold")

    ax_disp.set_xlabel("Time (yr)")
    ax_disp.set_ylabel("Displacement magnitude (m)")
    ax_disp.grid(True, linestyle="--", linewidth=0.6, alpha=0.6)
    ax_disp.text(
        0.02,
        0.96,
        "(b)",
        transform=ax_disp.transAxes,
        va="top",
        ha="left",
        fontsize=11,
        fontweight="bold",
    )

    handles, labels = ax_disp.get_legend_handles_labels()
    if handles:
        ax_disp.legend(handles, labels, loc="best", frameon=False)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    print(f"Saved: {output_path.resolve()}")


if __name__ == "__main__":
    main()
