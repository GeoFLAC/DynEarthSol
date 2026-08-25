#!/usr/bin/env python3

"""Verify that 3-D monitor query coordinates preserve x, y, and z."""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
DEFAULT_EXE = REPO_ROOT / "dynearthsol3d"
DEFAULT_CFG = HERE / "monitor_coordinates_3d.cfg"
EXPECTED_QUERY = {"query_x": 0.25, "query_y": 0.50, "query_z": -0.75}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the 3-D monitor coordinate-binding regression."
    )
    parser.add_argument("--exe", type=Path, default=DEFAULT_EXE)
    parser.add_argument("--cfg", type=Path, default=DEFAULT_CFG)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    exe = args.exe.expanduser().resolve()
    cfg = args.cfg.expanduser().resolve()
    if not exe.is_file():
        raise FileNotFoundError(f"executable not found: {exe}")

    with tempfile.TemporaryDirectory(prefix="des-monitor-coordinates-") as tmp:
        run_dir = Path(tmp)
        env = os.environ.copy()
        env["PYTHONDONTWRITEBYTECODE"] = "1"
        result = subprocess.run(
            [str(exe), str(cfg)],
            cwd=run_dir,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            print(result.stdout, file=sys.stderr)
            raise RuntimeError(f"DynEarthSol exited with {result.returncode}")

        csv_path = run_dir / "monitor_coordinates_point_0.csv"
        if not csv_path.is_file():
            raise AssertionError(f"missing monitor output: {csv_path}")
        with csv_path.open(newline="", encoding="ascii") as stream:
            rows = list(csv.DictReader(stream))
        if not rows:
            raise AssertionError("monitor output has no data rows")

        initial = rows[0]
        for column, expected in EXPECTED_QUERY.items():
            if column not in initial:
                raise AssertionError(f"missing monitor column {column!r}")
            actual = float(initial[column])
            if abs(actual - expected) > 1.0e-15:
                raise AssertionError(
                    f"{column}: expected {expected:.17g}, got {actual:.17g}"
                )

    print("3-D monitor coordinate regression: PASS")


if __name__ == "__main__":
    main()
