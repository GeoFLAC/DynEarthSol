#!/usr/bin/env python3
"""Generate and run the paper's local EP/RSF benchmark cases."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

sys.dont_write_bytecode = True


SHEAR_MAX_STEPS = 200_000
SHEAR_MONITOR_STEP_INTERVAL = 4_000
TRANSIENT_MONITOR_STEP_INTERVAL = 100
SHEAR_FIXED_DT_S = 0.01
HEALING_MAX_STEPS = 1_546
HEALING_FIXED_DT_S = 259_200.0


@dataclass(frozen=True)
class BenchmarkCase:
    name: str
    group: str
    rheology_type: str
    friction_angle_deg: float
    direct_a: float
    evolution_b: float
    characteristic_distance: float
    characteristic_velocity: float
    state_var_model: int


@dataclass(frozen=True)
class RunSpec:
    name: str
    case: BenchmarkCase
    max_steps: int
    output_step_interval: int
    monitor_step_interval: int
    fixed_dt_s: float
    upper_boundary_x_mode: int
    upper_boundary_x_velocity: float
    characteristic_speed_line: str


BENCHMARK_CASES: tuple[BenchmarkCase, ...] = (
    BenchmarkCase("ep_phi30", "ep", "elasto-plastic", 30.0, 0.0, 0.0, 1e-3, 1e-5, 0),
    BenchmarkCase("ep_phi20", "ep", "elasto-plastic", 20.0, 0.0, 0.0, 1e-3, 1e-5, 0),
    BenchmarkCase("ep_phi10", "ep", "elasto-plastic", 10.0, 0.0, 0.0, 1e-3, 1e-5, 0),
    BenchmarkCase("steady_ab_neg_v0_1e-6", "steady_ab_neg", "elasto-plastic-rate-state-friction", 30.0, 0.2, 0.3, 1e-3, 1e-6, 0),
    BenchmarkCase("steady_ab_neg_v0_1e-5", "steady_ab_neg", "elasto-plastic-rate-state-friction", 30.0, 0.2, 0.3, 1e-3, 1e-5, 0),
    BenchmarkCase("steady_ab_neg_v0_1e-4", "steady_ab_neg", "elasto-plastic-rate-state-friction", 30.0, 0.2, 0.3, 1e-3, 1e-4, 0),
    BenchmarkCase("aging_ab_neg_dc_1e-3", "aging_ab_neg", "elasto-plastic-rate-state-friction", 30.0, 0.2, 0.3, 1e-3, 1e-5, 1),
    BenchmarkCase("aging_ab_neg_dc_3e-3", "aging_ab_neg", "elasto-plastic-rate-state-friction", 30.0, 0.2, 0.3, 3e-3, 1e-5, 1),
    BenchmarkCase("aging_ab_neg_dc_1e-2", "aging_ab_neg", "elasto-plastic-rate-state-friction", 30.0, 0.2, 0.3, 1e-2, 1e-5, 1),
)

HEALING_CASE = BenchmarkCase(
    "healing_zero_velocity",
    "healing",
    "elasto-plastic-rate-state-friction",
    30.0,
    0.011,
    0.017,
    1e-2,
    4e-9,
    1,
)

GROUP_ORDER = ("ep", "steady_ab_neg", "aging_ab_neg")


def candidate_executables(script_dir: Path) -> Iterable[Path]:
    repo_root = script_dir.parent.parent
    yield repo_root / "dynearthsol2d"
    yield repo_root / "binaries" / "dynearthsol2d"


def resolve_executable(script_dir: Path, explicit: str | None) -> Path:
    if explicit:
        path = Path(explicit).expanduser().resolve()
        if path.is_file() and os.access(path, os.X_OK):
            return path
        raise FileNotFoundError(f"--exe is not executable: {path}")

    env_exe = os.environ.get("DYNEXE")
    if env_exe:
        path = Path(env_exe).expanduser().resolve()
        if path.is_file() and os.access(path, os.X_OK):
            return path
        raise FileNotFoundError(f"DYNEXE is set but not executable: {path}")

    for candidate in candidate_executables(script_dir):
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate

    raise FileNotFoundError(
        "DynEarthSol 2-D executable not found. Build dynearthsol2d in the repo root or pass --exe."
    )


def format_float(value: float) -> str:
    return f"{value:.17g}"


def render_cfg(template_text: str, spec: RunSpec) -> str:
    case = spec.case
    replacements = {
        "__MODELNAME__": "result",
        "__MAX_STEPS__": str(spec.max_steps),
        "__OUTPUT_STEP_INTERVAL__": str(spec.output_step_interval),
        "__FIXED_DT__": format_float(spec.fixed_dt_s),
        "__CHARACTERISTIC_SPEED_CONTROL__": spec.characteristic_speed_line,
        "__UPPER_BOUNDARY_X_MODE__": str(spec.upper_boundary_x_mode),
        "__UPPER_BOUNDARY_X_VELOCITY__": format_float(spec.upper_boundary_x_velocity),
        "__RHEOLOGY_TYPE__": case.rheology_type,
        "__FRICTION_ANGLE_DEG__": format_float(case.friction_angle_deg),
        "__DIRECT_A__": format_float(case.direct_a),
        "__EVOLUTION_B__": format_float(case.evolution_b),
        "__CHARACTERISTIC_DISTANCE__": format_float(case.characteristic_distance),
        "__CHARACTERISTIC_VELOCITY__": format_float(case.characteristic_velocity),
        "__STATE_VAR_MODEL__": str(case.state_var_model),
        "__RSF_RATE_CONTROL__": (
            "" if case.group == "ep" else "rsf_slip_rate_projection_option = 1"
        ),
        "__MONITOR_PREFIX__": "monitor",
        "__MONITOR_STEP_INTERVAL__": str(spec.monitor_step_interval),
    }

    rendered = template_text
    for old, new in replacements.items():
        rendered = rendered.replace(old, new)
    if "__" in rendered:
        raise ValueError(f"Unsubstituted CFG placeholder in {spec.name}")
    return rendered


def select_cases(groups: list[str], names: list[str]) -> list[BenchmarkCase]:
    selected = list(BENCHMARK_CASES)
    if groups:
        allowed = set(groups)
        selected = [case for case in selected if case.group in allowed]
    if names:
        known = {case.name for case in BENCHMARK_CASES}
        missing = sorted(set(names) - known)
        if missing:
            raise ValueError(f"Unknown benchmark case(s): {', '.join(missing)}")
        allowed = set(names)
        selected = [case for case in selected if case.name in allowed]
    if not selected:
        raise ValueError("No benchmark cases selected.")
    return selected


def make_run_specs(
    cases: list[BenchmarkCase],
    max_steps: int,
    output_step_interval: int,
    monitor_step_interval: int,
    aging_transient: bool,
    include_healing: bool,
) -> list[RunSpec]:
    specs = [
        RunSpec(
            case.name,
            case,
            max_steps,
            output_step_interval,
            monitor_step_interval,
            SHEAR_FIXED_DT_S,
            4,
            1e-5,
            "",
        )
        for case in cases
    ]
    if aging_transient:
        specs.extend(
            RunSpec(
                f"{case.name}_transient",
                case,
                max_steps,
                output_step_interval,
                TRANSIENT_MONITOR_STEP_INTERVAL,
                SHEAR_FIXED_DT_S,
                4,
                1e-5,
                "",
            )
            for case in cases
            if case.group == "aging_ab_neg"
        )
    if include_healing:
        specs.append(
            RunSpec(
                HEALING_CASE.name,
                HEALING_CASE,
                HEALING_MAX_STEPS,
                HEALING_MAX_STEPS,
                1,
                HEALING_FIXED_DT_S,
                1,
                0.0,
                "characteristic_speed = 4e-9",
            )
        )
    return specs


def cleanup_case_dir(case_dir: Path) -> None:
    for path in case_dir.iterdir():
        if path.name == ".gitkeep":
            continue
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()


def write_case_cfg(case_dir: Path, cfg_text: str) -> Path:
    cfg_path = case_dir / "simple_shear_box.cfg"
    cfg_path.write_text(cfg_text, encoding="utf-8")
    return cfg_path


def cleanup_pycache(root_dir: Path) -> None:
    for path in root_dir.rglob("__pycache__"):
        if path.is_dir():
            shutil.rmtree(path, ignore_errors=True)


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Run the paper's local EP/RSF benchmark suite.")
    parser.add_argument("--exe", default=None, help="Path to dynearthsol2d. If omitted, auto-detect.")
    parser.add_argument(
        "--groups",
        nargs="+",
        choices=GROUP_ORDER,
        default=list(GROUP_ORDER),
        help="Simple-shear benchmark groups to run.",
    )
    parser.add_argument(
        "--cases",
        nargs="+",
        default=[],
        help="Optional explicit simple-shear case-name filter.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=script_dir / "runs",
        help="Directory where generated case folders are written.",
    )
    parser.add_argument(
        "--max-steps",
        type=int,
        default=SHEAR_MAX_STEPS,
        help="Number of 0.01 s steps for each simple-shear run.",
    )
    parser.add_argument(
        "--monitor-step-interval",
        type=int,
        default=SHEAR_MONITOR_STEP_INTERVAL,
        help="Simple-shear monitor interval in solver steps.",
    )
    parser.add_argument(
        "--output-step-interval",
        type=int,
        default=None,
        help="Full-field output interval. Defaults to --max-steps.",
    )
    parser.add_argument(
        "--aging-transient",
        action="store_true",
        help="Also run aging cases with 1 s monitor output for the figure inset.",
    )
    parser.add_argument(
        "--healing",
        action="store_true",
        help="Also run the paper's separate zero-velocity healing test.",
    )
    parser.add_argument(
        "--skip-run",
        action="store_true",
        help="Only generate CFG files and folders.",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Delete existing files inside each selected case directory before regenerating.",
    )
    args = parser.parse_args()

    output_step_interval = args.output_step_interval or args.max_steps
    if args.max_steps < 1:
        raise ValueError("--max-steps must be >= 1")
    if args.monitor_step_interval < 1:
        raise ValueError("--monitor-step-interval must be >= 1")
    if output_step_interval < 1:
        raise ValueError("--output-step-interval must be >= 1")

    selected = select_cases(args.groups, args.cases)
    specs = make_run_specs(
        selected,
        args.max_steps,
        output_step_interval,
        args.monitor_step_interval,
        args.aging_transient,
        args.healing,
    )
    template_text = (script_dir / "simple_shear_base.cfg").read_text(encoding="utf-8")
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    exe = None if args.skip_run else resolve_executable(script_dir, args.exe)
    if exe is not None:
        print(f"[exe] {exe}")

    for index, spec in enumerate(specs, start=1):
        case_dir = output_root / spec.name
        case_dir.mkdir(parents=True, exist_ok=True)
        if args.clean:
            cleanup_case_dir(case_dir)
            case_dir.mkdir(parents=True, exist_ok=True)

        cfg_path = write_case_cfg(case_dir, render_cfg(template_text, spec))
        print(f"[prep] {index:02d}/{len(specs)} {spec.name} -> {cfg_path}")
        if not args.skip_run:
            subprocess.run([str(exe), str(cfg_path)], cwd=case_dir, check=True)

    print("[done] benchmark case generation complete")


if __name__ == "__main__":
    _script_dir = Path(__file__).resolve().parent
    try:
        main()
    finally:
        cleanup_pycache(_script_dir)
