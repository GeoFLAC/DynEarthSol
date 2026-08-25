#!/usr/bin/env python3

"""Run the fluid point-injection schedule-equivalence check."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
DEFAULT_EXE = REPO_ROOT / "dynearthsol2d"
DEFAULT_CFG = HERE / "point_injection_base.cfg"
DEFAULT_RUN_DIR = HERE / "runs"
DEFAULT_CHECK = HERE / "check_point_injection.py"
YEAR2SEC = 365.2422 * 86400.0
INJECTION_DURATION_S = 20.0e-6
INJECTION_RATE = 1.0e-6
TOTAL_AMOUNT = INJECTION_RATE * INJECTION_DURATION_S


@dataclass(frozen=True)
class InjectionCase:
    name: str
    enabled: str
    rate_model: str
    rate: float
    total_amount: float


CASES = (
    InjectionCase("source_disabled", "no", "constant_rate", 0.0, 0.0),
    InjectionCase("constant_rate", "yes", "constant_rate", INJECTION_RATE, 0.0),
    InjectionCase("total_amount", "yes", "total_amount", 0.0, TOTAL_AMOUNT),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--exe", type=Path, default=DEFAULT_EXE)
    parser.add_argument("--cfg", type=Path, default=DEFAULT_CFG)
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("--skip-check", action="store_true")
    parser.add_argument("--rtol", type=float, default=1.0e-10)
    parser.add_argument("--atol", type=float, default=1.0e-8)
    return parser.parse_args()


def render_cfg(template: str, case: InjectionCase) -> str:
    replacements = {
        "__MODELNAME__": f"point_injection_{case.name}",
        "__ENABLED__": case.enabled,
        "__RATE_MODEL__": case.rate_model,
        "__RATE__": f"{case.rate:.16e}",
        "__TOTAL_AMOUNT__": f"{case.total_amount:.16e}",
        "__DURATION_YR__": f"{INJECTION_DURATION_S / YEAR2SEC:.16e}",
    }
    for key, value in replacements.items():
        template = template.replace(key, value)
    return template


def run_case(exe: Path, template: str, run_root: Path, case: InjectionCase) -> Path:
    case_dir = run_root / case.name
    case_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = case_dir / "input.cfg"
    cfg_path.write_text(render_cfg(template, case), encoding="ascii")
    log_path = case_dir / "run.log"
    with log_path.open("w", encoding="ascii") as log_file:
        subprocess.run(
            [str(exe), str(cfg_path)],
            cwd=case_dir,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            timeout=60,
            check=True,
        )
    print(f"PASS {case.name}: {case_dir}", flush=True)
    return case_dir


def main() -> None:
    args = parse_args()
    run_root = args.run_dir.resolve()
    if args.clean and run_root.exists():
        shutil.rmtree(run_root)
    run_root.mkdir(parents=True, exist_ok=True)

    exe = args.exe.resolve()
    template = args.cfg.resolve().read_text(encoding="ascii")
    case_dirs = [run_case(exe, template, run_root, case) for case in CASES]

    if not args.skip_check:
        subprocess.run(
            [
                sys.executable,
                str(DEFAULT_CHECK),
                str(case_dirs[1]),
                str(case_dirs[2]),
                str(case_dirs[0]),
                "--rtol",
                str(args.rtol),
                "--atol",
                str(args.atol),
            ],
            check=True,
        )


if __name__ == "__main__":
    main()
