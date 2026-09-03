#!/usr/bin/env python3
"""Focused checks for selectable RSF slip-rate calculations."""

from __future__ import annotations

import argparse
import csv
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
BENCHMARK_DIR = REPO_ROOT / "benchmarks" / "simple_shear_rsf"
sys.path.insert(0, str(BENCHMARK_DIR))

from run_simple_shear_benchmark import BenchmarkCase, render_cfg  # noqa: E402


VX_TOP = 1.0e-5


def replace_once(text: str, old: str, new: str) -> str:
    if text.count(old) != 1:
        raise ValueError(f"Expected exactly one occurrence of {old!r}")
    return text.replace(old, new, 1)


def make_cfg(
    *,
    fixed_dt: float,
    dc: float,
    state_model: int = 1,
    rate_option: int | None = 1,
    max_steps: int = 1,
) -> str:
    case = BenchmarkCase(
        name="rsf_control",
        group="aging_ab_neg",
        rheology_type="elasto-plastic-rate-state-friction",
        friction_angle_deg=30.0,
        direct_a=0.2,
        evolution_b=0.3,
        characteristic_distance=dc,
        characteristic_velocity=1.0e-5,
        state_var_model=state_model,
    )
    template = (BENCHMARK_DIR / "simple_shear_base.cfg").read_text(encoding="utf-8")
    cfg = render_cfg(template, case, max_steps, max_steps, 1)
    cfg = replace_once(cfg, "fixed_dt = 1.0", f"fixed_dt = {fixed_dt:.17e}")
    if rate_option is not None:
        cfg = replace_once(
            cfg,
            "damping_option = 1",
            "damping_option = 1\n"
            f"rsf_slip_rate_projection_option = {rate_option}",
        )
    return cfg


def run_cfg(
    exe: Path,
    root: Path,
    name: str,
    cfg: str,
    *,
    expect_success: bool = True,
):
    run_dir = root / name
    run_dir.mkdir()
    cfg_path = run_dir / "case.cfg"
    cfg_path.write_text(cfg, encoding="utf-8")
    result = subprocess.run(
        [str(exe), str(cfg_path)],
        cwd=run_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if expect_success and result.returncode != 0:
        raise RuntimeError(f"{name} failed with {result.returncode}:\n{result.stdout}")
    if not expect_success and result.returncode == 0:
        raise AssertionError(f"{name} unexpectedly succeeded")
    return run_dir, result


def monitor_rows(run_dir: Path) -> list[list[dict[str, str]]]:
    paths = sorted(run_dir.glob("monitor_point_*.csv"))
    if len(paths) != 2:
        raise AssertionError(f"Expected two monitor files in {run_dir}, found {len(paths)}")
    output = []
    for path in paths:
        with path.open("r", encoding="utf-8", newline="") as handle:
            output.append(list(csv.DictReader(handle)))
    return output


def assert_close(
    label: str,
    actual: float,
    expected: float,
    rel: float = 2.0e-12,
) -> None:
    scale = max(abs(expected), 1.0e-30)
    if not math.isfinite(actual) or abs(actual - expected) > rel * scale:
        raise AssertionError(f"{label}: {actual:.17e} != {expected:.17e}")


def check_invariant_rate(exe: Path, root: Path) -> None:
    dt = 1.0e-2
    dc = 1.0e-3
    run_dir, _ = run_cfg(
        exe,
        root,
        "invariant_rate_option1",
        make_cfg(fixed_dt=dt, dc=dc, rate_option=1),
    )
    rows = monitor_rows(run_dir)
    rate = VX_TOP / math.sqrt(2.0)
    theta0 = dc / 1.0e-5
    theta1 = theta0 + dt * (1.0 - rate * theta0 / dc)
    mu0 = math.tan(math.radians(30.0))
    mu1 = (
        mu0
        + 0.2 * math.log(rate / 1.0e-5)
        + 0.3 * math.log(1.0e-5 * theta1 / dc)
    )
    for index, point_rows in enumerate(rows):
        assert_close(
            f"point {index} theta",
            float(point_rows[1]["state_variable"]),
            theta1,
        )
        assert_close(
            f"point {index} friction",
            float(point_rows[1]["dynamic_friction"]),
            mu1,
        )
    print("[ok] option 1 uses the total-strain invariant in both elements")


def check_invariant_restart_refresh(exe: Path, root: Path) -> None:
    seed_cfg = make_cfg(
        fixed_dt=1.0e-2,
        dc=1.0e-3,
        rate_option=1,
        max_steps=1,
    )
    seed_cfg = replace_once(
        seed_cfg,
        "output_step_interval = 1",
        "output_step_interval = 1\ncheckpoint_frame_interval = 1",
    )
    seed_dir, _ = run_cfg(exe, root, "invariant_restart_seed", seed_cfg)
    seed_save = seed_dir / "result.save.000001.vtkhdf"
    with h5py.File(seed_save, "r") as output:
        expected_time = output["time_sec"][()].item()
        expected_theta = output["friction state variable"][...]
        source = output["/VTKHDF/grid/CellData/strain-rate"]
        source_values = source[...]
        if not any(
            abs(float(value)) > 0.0
            for element in source_values
            for value in element
        ):
            raise AssertionError("restart seed strain rate is already zero")

    restart_outputs = []
    for label in ("clean", "poisoned"):
        if label == "poisoned":
            with h5py.File(seed_save, "r+") as output:
                output["/VTKHDF/grid/CellData/strain-rate"][...] = 0.0
        restart_cfg = make_cfg(
            fixed_dt=1.0e-2,
            dc=1.0e-3,
            rate_option=1,
            max_steps=2,
        )
        restart_cfg = replace_once(
            restart_cfg,
            "modelname = result",
            f"modelname = restarted_{label}\n"
            "is_restarting = yes\n"
            f"restarting_from_modelname = {seed_dir / 'result'}\n"
            "restarting_from_frame = 1",
        )
        restart_dir, _ = run_cfg(
            exe, root, f"invariant_restart_{label}", restart_cfg
        )
        restart_save = (
            restart_dir / f"restarted_{label}.save.000001.vtkhdf"
        )
        with h5py.File(restart_save, "r") as output:
            restart_outputs.append(
                (
                    output["time_sec"][()].item(),
                    output["friction state variable"][...],
                    output["dynamic friction coefficient"][...],
                )
            )

    clean, poisoned = restart_outputs
    if (
        clean[1].shape != expected_theta.shape
        or poisoned[1].shape != clean[1].shape
    ):
        raise AssertionError("restart state-variable shapes changed")
    if clean[2].shape != poisoned[2].shape:
        raise AssertionError("restart friction shapes changed")
    assert_close("clean restart frame time", clean[0], expected_time)
    assert_close("poisoned restart frame time", poisoned[0], expected_time)
    for element in range(len(expected_theta)):
        assert_close(
            f"clean restart element {element} theta",
            float(clean[1][element]),
            float(expected_theta[element]),
        )
        assert_close(
            f"poisoned restart element {element} theta",
            float(poisoned[1][element]),
            float(expected_theta[element]),
        )
        assert_close(
            f"restart element {element} derived friction",
            float(poisoned[2][element]),
            float(clean[2][element]),
        )
    print("[ok] option 1 rebuilds its derived rate on restart")


def check_maximum_shear_default(exe: Path, root: Path) -> None:
    explicit_dir, _ = run_cfg(
        exe,
        root,
        "maximum_shear_explicit",
        make_cfg(fixed_dt=1.0e-2, dc=1.0e-3, state_model=0, rate_option=0),
    )
    default_dir, _ = run_cfg(
        exe,
        root,
        "maximum_shear_default",
        make_cfg(fixed_dt=1.0e-2, dc=1.0e-3, state_model=0, rate_option=None),
    )
    explicit_rows = monitor_rows(explicit_dir)
    default_rows = monitor_rows(default_dir)
    if len(default_rows) != len(explicit_rows):
        raise AssertionError("default and explicit option 0 point counts differ")
    for point, (explicit, default) in enumerate(zip(explicit_rows, default_rows)):
        if len(default) != len(explicit):
            raise AssertionError(
                f"default and explicit option 0 row counts differ at point {point}"
            )
        for row, (explicit_record, default_record) in enumerate(
            zip(explicit, default)
        ):
            for field in (
                "time_s",
                "stress_0",
                "stress_1",
                "stress_2",
                "state_variable",
                "dynamic_friction",
            ):
                assert_close(
                    f"default option 0 point {point} row {row} {field}",
                    float(default_record[field]),
                    float(explicit_record[field]),
                )

    friction = [
        float(point_rows[1]["dynamic_friction"])
        for point_rows in explicit_rows
    ]
    if (
        not all(math.isfinite(value) for value in friction)
        or math.isclose(friction[0], friction[1])
    ):
        raise AssertionError(
            "option 0 did not preserve element-wise maximum-shear projection: "
            f"{friction}"
        )
    mu0 = math.tan(math.radians(30.0))
    expected = sorted(
        mu0 + (0.2 - 0.3) * math.log(rate / 1.0e-5)
        for rate in (
            VX_TOP / (3.0 * math.sqrt(2.0)),
            2.0 * VX_TOP / (3.0 * math.sqrt(2.0)),
        )
    )
    for index, (actual, target) in enumerate(zip(sorted(friction), expected)):
        assert_close(f"option 0 analytic friction {index}", actual, target)
    print("[ok] the default matches option 0 maximum-shear projection")


def check_invalid_option(exe: Path, root: Path) -> None:
    _, result = run_cfg(
        exe,
        root,
        "invalid_projection",
        make_cfg(fixed_dt=1.0e-2, dc=1.0e-3, rate_option=2),
        expect_success=False,
    )
    if "must be 0 or 1" not in result.stdout:
        raise AssertionError(
            "invalid projection did not report its accepted range:\n"
            f"{result.stdout}"
        )
    print("[ok] unsupported projection options are rejected")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", required=True, type=Path)
    args = parser.parse_args()
    exe = args.exe.expanduser().resolve()
    if not exe.is_file() or not os.access(exe, os.X_OK):
        raise FileNotFoundError(f"Executable not found: {exe}")

    with tempfile.TemporaryDirectory(prefix="des_rsf_controls_") as temp:
        root = Path(temp)
        check_invariant_rate(exe, root)
        check_invariant_restart_refresh(exe, root)
        check_maximum_shear_default(exe, root)
        check_invalid_option(exe, root)


if __name__ == "__main__":
    main()
