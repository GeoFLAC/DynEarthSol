#!/usr/bin/env python3

"""Check the hydraulic CFL against the current domain-wide diffusivity."""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

from hydraulic_cfl_base_checks import run_base_regressions


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
DEFAULT_EXE = REPO_ROOT / "dynearthsol2d"
DEFAULT_CFG = HERE / "hydraulic_cfl_base.cfg"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run feature-on/off hydraulic CFL regression cases."
    )
    parser.add_argument("--exe", type=Path, default=DEFAULT_EXE)
    parser.add_argument("--cfg", type=Path, default=DEFAULT_CFG)
    parser.add_argument(
        "--run-dir",
        type=Path,
        help="Keep case outputs below this directory (default: temporary directory).",
    )
    return parser.parse_args()


def run_all(exe: Path, template: str, run_root: Path) -> None:
    run_base_regressions(exe, template, run_root)
    print("hydraulic CFL regression: PASS")


def main() -> None:
    args = parse_args()
    exe = args.exe.resolve()
    template = args.cfg.resolve().read_text(encoding="ascii")
    if args.run_dir:
        run_root = args.run_dir.resolve()
        run_root.mkdir(parents=True, exist_ok=True)
        run_all(exe, template, run_root)
    else:
        with tempfile.TemporaryDirectory(prefix="des-hydraulic-cfl-") as tmp:
            run_all(exe, template, Path(tmp))


if __name__ == "__main__":
    main()
