#!/usr/bin/env python3
"""Run the simple-shear benchmark suite from one base CFG template."""

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


BENCHMARK_CASES: tuple[BenchmarkCase, ...] = (
    BenchmarkCase("ep_phi30", "ep", "elasto-plastic", 30.0, 0.0, 0.0, 1e-3, 1e-5, 0),
    BenchmarkCase("ep_phi20", "ep", "elasto-plastic", 20.0, 0.0, 0.0, 1e-3, 1e-5, 0),
    BenchmarkCase("ep_phi10", "ep", "elasto-plastic", 10.0, 0.0, 0.0, 1e-3, 1e-5, 0),
    BenchmarkCase("steady_ab_pos_v0_1e-6", "steady_ab_pos", "elasto-plastic-rate-state-friction", 30.0, 0.2, 0.1, 1e-3, 1e-6, 0),
    BenchmarkCase("steady_ab_pos_v0_1e-5", "steady_ab_pos", "elasto-plastic-rate-state-friction", 30.0, 0.2, 0.1, 1e-3, 1e-5, 0),
    BenchmarkCase("steady_ab_pos_v0_1e-4", "steady_ab_pos", "elasto-plastic-rate-state-friction", 30.0, 0.2, 0.1, 1e-3, 1e-4, 0),
    BenchmarkCase("steady_ab_neg_v0_1e-6", "steady_ab_neg", "elasto-plastic-rate-state-friction", 30.0, 0.2, 0.3, 1e-3, 1e-6, 0),
    BenchmarkCase("steady_ab_neg_v0_1e-5", "steady_ab_neg", "elasto-plastic-rate-state-friction", 30.0, 0.2, 0.3, 1e-3, 1e-5, 0),
    BenchmarkCase("steady_ab_neg_v0_1e-4", "steady_ab_neg", "elasto-plastic-rate-state-friction", 30.0, 0.2, 0.3, 1e-3, 1e-4, 0),
    BenchmarkCase("aging_ab_neg_dc_1e-6", "aging_ab_neg", "elasto-plastic-rate-state-friction", 30.0, 0.2, 0.3, 1e-6, 1e-5, 1),
    BenchmarkCase("aging_ab_neg_dc_1e-5", "aging_ab_neg", "elasto-plastic-rate-state-friction", 30.0, 0.2, 0.3, 1e-5, 1e-5, 1),
    BenchmarkCase("aging_ab_neg_dc_1e-4", "aging_ab_neg", "elasto-plastic-rate-state-friction", 30.0, 0.2, 0.3, 1e-4, 1e-5, 1),
)

GROUP_ORDER = ("ep", "steady_ab_pos", "steady_ab_neg", "aging_ab_neg")


def candidate_executables(script_dir: Path) -> Iterable[Path]:
    repo_root = script_dir.parent.parent
    yield repo_root / "dynearthsol2d"
    yield repo_root / "binaries" / "dynearthsol2d"


def resolve_executable(script_dir: Path, explicit: str | None) -> Path:
    if explicit:
        p = Path(explicit).expanduser().resolve()
        if p.is_file() and os.access(p, os.X_OK):
            return p
        raise FileNotFoundError(f"--exe is not executable: {p}")

    env_exe = os.environ.get("DYNEXE")
    if env_exe:
        p = Path(env_exe).expanduser().resolve()
        if p.is_file() and os.access(p, os.X_OK):
            return p
        raise FileNotFoundError(f"DYNEXE is set but not executable: {p}")

    for candidate in candidate_executables(script_dir):
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate

    raise FileNotFoundError(
        "DynEarthSol 2-D executable not found. Build dynearthsol2d in the repo root or pass --exe."
    )


def format_float(value: float) -> str:
    return f"{value:.17g}"


def render_cfg(template_text: str, case: BenchmarkCase, max_steps: int, output_step_interval: int, monitor_step_interval: int) -> str:
    replacements = {
        "__MODELNAME__": "result",
        "__MAX_STEPS__": str(max_steps),
        "__OUTPUT_STEP_INTERVAL__": str(output_step_interval),
        "__RHEOLOGY_TYPE__": case.rheology_type,
        "__FRICTION_ANGLE_DEG__": format_float(case.friction_angle_deg),
        "__DIRECT_A__": format_float(case.direct_a),
        "__EVOLUTION_B__": format_float(case.evolution_b),
        "__CHARACTERISTIC_DISTANCE__": format_float(case.characteristic_distance),
        "__CHARACTERISTIC_VELOCITY__": format_float(case.characteristic_velocity),
        "__STATE_VAR_MODEL__": str(case.state_var_model),
        "__MONITOR_PREFIX__": "monitor",
        "__MONITOR_STEP_INTERVAL__": str(monitor_step_interval),
    }

    rendered = template_text
    for old, new in replacements.items():
        rendered = rendered.replace(old, new)
    return rendered


def select_cases(groups: list[str], names: list[str]) -> list[BenchmarkCase]:
    selected = list(BENCHMARK_CASES)
    if groups:
        allowed = set(groups)
        selected = [case for case in selected if case.group in allowed]
    if names:
        allowed = set(names)
        selected = [case for case in selected if case.name in allowed]
    if not selected:
        raise ValueError("No benchmark cases selected.")
    return selected


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
    parser = argparse.ArgumentParser(description="Run monitor-based simple-shear benchmark sweeps.")
    parser.add_argument("--exe", default=None, help="Path to dynearthsol2d. If omitted, auto-detect.")
    parser.add_argument(
        "--groups",
        nargs="+",
        choices=GROUP_ORDER,
        default=list(GROUP_ORDER),
        help="Benchmark groups to run.",
    )
    parser.add_argument(
        "--cases",
        nargs="+",
        default=[],
        help="Optional explicit case-name filter.",
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
        default=2000,
        help="Number of time steps for each run.",
    )
    parser.add_argument(
        "--monitor-step-interval",
        type=int,
        default=40,
        help="Monitor output interval in solver steps.",
    )
    parser.add_argument(
        "--output-step-interval",
        type=int,
        default=2000,
        help="Full-field output interval. Keep large because the benchmark reads monitor CSVs.",
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

    if args.max_steps < 1:
        raise ValueError("--max-steps must be >= 1")
    if args.monitor_step_interval < 1:
        raise ValueError("--monitor-step-interval must be >= 1")
    if args.output_step_interval < 1:
        raise ValueError("--output-step-interval must be >= 1")

    selected_cases = select_cases(args.groups, args.cases)
    template_path = script_dir / "simple_shear_base.cfg"
    template_text = template_path.read_text(encoding="utf-8")
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    exe = None if args.skip_run else resolve_executable(script_dir, args.exe)
    if exe is not None:
        print(f"[exe] {exe}")

    for index, case in enumerate(selected_cases, start=1):
        case_dir = output_root / case.name
        case_dir.mkdir(parents=True, exist_ok=True)
        if args.clean:
            cleanup_case_dir(case_dir)
            case_dir.mkdir(parents=True, exist_ok=True)

        cfg_text = render_cfg(
            template_text=template_text,
            case=case,
            max_steps=args.max_steps,
            output_step_interval=args.output_step_interval,
            monitor_step_interval=args.monitor_step_interval,
        )
        cfg_path = write_case_cfg(case_dir, cfg_text)
        print(f"[prep] {index:02d}/{len(selected_cases)} {case.name} -> {cfg_path}")

        if args.skip_run:
            continue

        subprocess.run([str(exe), str(cfg_path)], cwd=case_dir, check=True)

    print("[done] benchmark case generation complete")


if __name__ == "__main__":
    _script_dir = Path(__file__).resolve().parent
    try:
        main()
    finally:
        cleanup_pycache(_script_dir)
