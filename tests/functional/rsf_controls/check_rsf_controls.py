#!/usr/bin/env python3
"""Focused checks for RSF controls and the GVS elastic-speed ceiling."""

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
    dtheta_max: float | None = None,
    characteristic_velocity: float = 1.0e-5,
    top_velocity: float = VX_TOP,
    moving_mesh: bool | None = None,
    inertial_scaling: float = 1.0e5,
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
        characteristic_velocity=characteristic_velocity,
        state_var_model=state_model,
    )
    template = (BENCHMARK_DIR / "simple_shear_base.cfg").read_text(encoding="utf-8")
    cfg = render_cfg(template, case, max_steps, max_steps, 1)
    cfg = replace_once(cfg, "fixed_dt = 1.0", f"fixed_dt = {fixed_dt:.17e}")
    cfg = replace_once(
        cfg,
        "inertial_scaling = 1e5",
        f"inertial_scaling = {inertial_scaling:.17e}",
    )
    cfg = replace_once(
        cfg,
        "vbc_val_z1 = 1e-5",
        f"vbc_val_z1 = {top_velocity:.17e}",
    )
    controls = ["damping_option = 1", "dt_fraction = 1"]
    if moving_mesh is not None:
        controls.append(
            "has_moving_mesh = " + ("yes" if moving_mesh else "no")
        )
    if rate_option is not None:
        controls.append(f"rsf_slip_rate_projection_option = {rate_option}")
    if dtheta_max is not None:
        controls.append(f"rsf_dtheta_max = {dtheta_max:.17e}")
    cfg = replace_once(cfg, "damping_option = 1", "\n".join(controls))
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


def check_adaptive_state_limit_disabled(exe: Path, root: Path) -> None:
    dc = 1.0e-3
    v0 = 1.0e-5
    expected_dt = 0.5 / math.sqrt(2.0)
    rate = VX_TOP / math.sqrt(2.0)
    theta0 = dc / v0
    expected_theta = theta0 + expected_dt * (
        1.0 - rate * theta0 / dc
    )
    outputs = []
    for label, limit in (("omitted", None), ("explicit_zero", 0.0)):
        run_dir, _ = run_cfg(
            exe,
            root,
            f"state_limit_{label}",
            make_cfg(
                fixed_dt=0.0,
                dc=dc,
                dtheta_max=limit,
                characteristic_velocity=v0,
                moving_mesh=False,
            ),
        )
        rows = monitor_rows(run_dir)
        outputs.append(rows)
        for point, point_rows in enumerate(rows):
            if len(point_rows) != 2:
                raise AssertionError(
                    f"{label} point {point}: expected two rows"
                )
            assert_close(
                f"{label} point {point} adaptive dt",
                float(point_rows[1]["time_s"]),
                expected_dt,
            )
            assert_close(
                f"{label} point {point} theta",
                float(point_rows[1]["state_variable"]),
                expected_theta,
            )

    for point, (omitted, explicit) in enumerate(zip(*outputs)):
        for row, (omitted_record, explicit_record) in enumerate(
            zip(omitted, explicit)
        ):
            for field in ("time_s", "state_variable", "dynamic_friction"):
                assert_close(
                    f"disabled parity point {point} row {row} {field}",
                    float(explicit_record[field]),
                    float(omitted_record[field]),
                )
    print("[ok] a disabled state bound preserves adaptive stepping")


def check_state_rate_limit(exe: Path, root: Path) -> None:
    dc = 1.0e-6
    fraction = 0.2
    v0 = 1.0e-6
    velocity = VX_TOP
    cfg = make_cfg(
        fixed_dt=0.0,
        dc=dc,
        dtheta_max=fraction,
        characteristic_velocity=v0,
        top_velocity=0.0,
        moving_mesh=False,
        max_steps=3,
    )
    cfg = replace_once(cfg, "vbc_z1 = 4", "vbc_z1 = 1")
    cfg = replace_once(
        cfg,
        "vbc_val_x1 = 0",
        f"vbc_val_x1 = {velocity:.17e}\n"
        "num_vbc_period_x1 = 2\n"
        "vbc_period_x1_time_in_yr = [0, 1.0e-12]\n"
        "vbc_period_x1_ratio = [1, 2]",
    )
    run_dir, _ = run_cfg(exe, root, "state_rate_limit", cfg)
    rows = monitor_rows(run_dir)
    rates = (
        velocity / math.sqrt(2.0),
        2.0 * velocity / math.sqrt(2.0),
        2.0 * velocity / math.sqrt(2.0),
    )
    theta = dc / v0
    for index, point_rows in enumerate(rows):
        if len(point_rows) != 4:
            raise AssertionError(
                f"point {index}: expected initial plus three rows, "
                f"got {len(point_rows)}"
            )
        expected_time = 0.0
        expected_theta = theta
        for step, record in enumerate(point_rows):
            assert_close(
                f"point {index} step {step} rate-limited time",
                float(record["time_s"]),
                expected_time,
            )
            assert_close(
                f"point {index} step {step} rate-limited theta",
                float(record["state_variable"]),
                expected_theta,
            )
            if step < len(rates):
                expected_dt = fraction * dc / rates[step]
                expected_time += expected_dt
                expected_theta += expected_dt * (
                    1.0 - rates[step] * expected_theta / dc
                )
    print("[ok] fixed mesh refreshes the bound from the current boundary rate")


def check_state_healing_limit(exe: Path, root: Path) -> None:
    dc = 1.0e-6
    fraction = 0.2
    v0 = 1.0e-5
    run_dir, _ = run_cfg(
        exe,
        root,
        "state_healing_limit",
        make_cfg(
            fixed_dt=0.0,
            dc=dc,
            dtheta_max=fraction,
            characteristic_velocity=v0,
            moving_mesh=True,
        ),
    )
    rows = monitor_rows(run_dir)
    theta0 = dc / v0
    expected_dt = fraction * theta0
    rate = VX_TOP / math.sqrt(2.0)
    theta1 = theta0 + expected_dt * (1.0 - rate * theta0 / dc)
    for index, point_rows in enumerate(rows):
        assert_close(
            f"point {index} healing-limited dt",
            float(point_rows[1]["time_s"]),
            expected_dt,
        )
        assert_close(
            f"point {index} healing-limited theta",
            float(point_rows[1]["state_variable"]),
            theta1,
        )
    print("[ok] the theta arm bounds fractional healing in one step")


def check_moving_mesh_rate_refresh(exe: Path, root: Path) -> None:
    dc = 1.0e-1
    fraction = 0.2
    velocity = 1.0
    cfg = make_cfg(
        fixed_dt=0.0,
        dc=dc,
        dtheta_max=fraction,
        characteristic_velocity=1.0e-1,
        top_velocity=0.0,
        moving_mesh=True,
        inertial_scaling=1.0,
        max_steps=2,
    )
    cfg = replace_once(
        cfg,
        "vbc_val_x1 = 0",
        f"vbc_val_x1 = {velocity:.17e}",
    )
    cfg = replace_once(cfg, "vbc_z1 = 4", "vbc_z1 = 1")
    run_dir, _ = run_cfg(
        exe,
        root,
        "moving_mesh_rate_refresh",
        cfg,
    )
    rows = monitor_rows(run_dir)
    for point, point_rows in enumerate(rows):
        if len(point_rows) != 3:
            raise AssertionError(
                f"moving point {point}: expected three rows"
            )
        length = 1.0
        theta = dc / 1.0e-1
        expected_time = 0.0
        for step, record in enumerate(point_rows):
            assert_close(
                f"moving point {point} step {step} time",
                float(record["time_s"]),
                expected_time,
            )
            assert_close(
                f"moving point {point} step {step} theta",
                float(record["state_variable"]),
                theta,
            )
            if step == 2:
                continue

            rate = velocity / math.sqrt(length * length + 1.0)
            dt = fraction * dc / rate
            if dt >= fraction * theta:
                raise AssertionError("moving test is not rate limited")
            expected_time += dt
            theta += dt * (1.0 - rate * theta / dc)
            length += velocity * dt
    print("[ok] moving mesh uses the analytically updated geometry")


def make_mass_scaling_cfg(
    reference: str | None,
    *,
    dynamic: bool = False,
    global_scaling: bool = True,
) -> str:
    cfg = make_cfg(
        fixed_dt=0.0,
        dc=1.0e-3,
        rate_option=None,
        top_velocity=0.0,
        moving_mesh=False,
        inertial_scaling=1.0e6,
    )
    cfg = replace_once(
        cfg,
        "use_global_velocity_scaling = true",
        "use_global_velocity_scaling = "
        + ("true" if global_scaling else "false"),
    )
    controls = [
        "damping_option = 0",
        "is_quasi_static = " + ("no" if dynamic else "yes"),
    ]
    if reference is not None:
        controls.append(f"mass_scaling_reference_speed = {reference}")
    cfg = replace_once(
        cfg, "damping_option = 1", "\n".join(controls)
    )
    cfg = replace_once(
        cfg,
        "rheology_type = elasto-plastic-rate-state-friction",
        "rheology_type = elastic",
    )
    cfg = replace_once(cfg, "bulk_modulus = [2.0e8]", "bulk_modulus = [5]")
    cfg = replace_once(cfg, "shear_modulus = [2.0e8]", "shear_modulus = [3]")
    cfg = replace_once(cfg, "vbc_val_x0 = 0", "vbc_val_x0 = -1e-5")
    cfg = replace_once(cfg, "vbc_val_x1 = 0", "vbc_val_x1 = 1e-5")
    cfg = replace_once(cfg, "vbc_z0 = 1", "vbc_z0 = 3")
    cfg = replace_once(
        cfg,
        "vbc_z1 = 4",
        "vbc_z1 = 0\nstress_bc_z1 = 3\nstress_val_z1 = -1",
    )
    cfg = replace_once(
        cfg,
        "points_x = [0.3333333333333333, 0.6666666666666666]",
        "points_x = [0, 1]",
    )
    cfg = replace_once(
        cfg,
        "points_y = [-0.6666666666666666, -0.3333333333333333]",
        "points_y = [0, 0]",
    )
    cfg = replace_once(
        cfg,
        "output_velocity = no",
        "output_velocity = yes\noutput_force = yes",
    )
    cfg = replace_once(
        cfg, "output_dynamic_friction = yes", "output_dynamic_friction = no"
    )
    cfg = replace_once(
        cfg, "output_state_variable = yes", "output_state_variable = no"
    )
    return cfg


def check_mass_scaling_reference(exe: Path, root: Path) -> None:
    cases = (
        ("default", None, False),
        ("shear", "shear", False),
        ("bulk", "bulk", False),
        ("dynamic", "bulk", True),
    )
    times = {}
    masses = {}
    for label, reference, dynamic in cases:
        run_dir, _ = run_cfg(
            exe,
            root,
            f"mass_{label}",
            make_mass_scaling_cfg(reference, dynamic=dynamic),
        )
        rows = monitor_rows(run_dir)
        times[label] = float(rows[0][1]["time_s"])
        masses[label] = []
        for point, point_rows in enumerate(rows):
            if len(point_rows) != 2:
                raise AssertionError(
                    f"mass {label} point {point}: expected two rows"
                )
            expected_vx = -1.0e-5 if point == 0 else 1.0e-5
            assert_close(
                f"mass {label} point {point} initial velocity_x",
                float(point_rows[0]["velocity_x"]),
                expected_vx,
            )
            dt = (
                float(point_rows[1]["time_s"])
                - float(point_rows[0]["time_s"])
            )
            dv = (
                float(point_rows[1]["velocity_z"])
                - float(point_rows[0]["velocity_z"])
            )
            force = float(point_rows[1]["force_z"])
            if (
                not all(math.isfinite(value) for value in (dt, dv, force))
                or dt <= 0.0
                or force == 0.0
                or dv == 0.0
                or force * dv <= 0.0
            ):
                raise AssertionError(
                    f"mass {label} point {point}: invalid response "
                    f"dt={dt}, force={force}, dv={dv}"
                )
            masses[label].append(dt * force / dv)

    altitude = 1.0 / math.sqrt(2.0)
    assert_close(
        "default shear-wave dt floor",
        times["default"],
        altitude / (5.0 * math.sqrt(3.0)),
    )
    assert_close("explicit shear dt", times["shear"], times["default"])
    assert_close(
        "bulk-wave dt floor",
        times["bulk"],
        altitude / (5.0 * math.sqrt(5.0)),
    )
    for point in range(2):
        assert_close(
            f"point {point} default pseudo-mass",
            masses["default"][point],
            masses["shear"][point],
        )
        assert_close(
            f"point {point} shear/bulk mass ratio",
            masses["shear"][point] / masses["bulk"][point],
            5.0 / 3.0,
        )
        assert_close(
            f"point {point} bulk physical mass",
            masses["bulk"][point],
            masses["dynamic"][point],
        )
    for index, (actual, expected) in enumerate(
        zip(sorted(masses["bulk"]), (1.0 / 6.0, 1.0 / 3.0))
    ):
        assert_close(f"bulk nodal mass {index}", actual, expected)
    print(
        "[ok] GVS preserves the shear default and reaches physical density "
        "at the bulk-wave ceiling"
    )

    invalid = (
        (
            "invalid_mass_reference",
            make_mass_scaling_cfg("pwave"),
            "must be 'shear' or 'bulk'",
        ),
        (
            "bulk_without_gvs",
            make_mass_scaling_cfg("bulk", global_scaling=False),
            "requires control.use_global_velocity_scaling=true",
        ),
    )
    for name, cfg, message in invalid:
        _, result = run_cfg(
            exe, root, name, cfg, expect_success=False
        )
        if message not in result.stdout:
            raise AssertionError(
                f"{name}: missing diagnostic {message!r}\n{result.stdout}"
            )
    print("[ok] unsupported mass-scaling selections are rejected")


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


def check_invalid_state_limits(exe: Path, root: Path) -> None:
    limit_isostasy = replace_once(
        make_cfg(fixed_dt=0.0, dc=1.0e-6, dtheta_max=0.2),
        "weakzone_option = 0",
        "weakzone_option = 0\nisostasy_adjustment_time_in_yr = 1",
    )
    limit_pt = replace_once(
        make_cfg(fixed_dt=0.0, dc=1.0e-6, dtheta_max=0.2),
        "dt_fraction = 1",
        "dt_fraction = 1\nhas_PT = yes",
    )
    limit_body_force = replace_once(
        make_cfg(fixed_dt=0.0, dc=1.0e-6, dtheta_max=0.2),
        "weakzone_option = 0",
        "weakzone_option = 0\nhas_body_force_adjustment = yes",
    )
    limit_non_rsf = replace_once(
        make_cfg(fixed_dt=0.0, dc=1.0e-6, dtheta_max=0.2),
        "rheology_type = elasto-plastic-rate-state-friction",
        "rheology_type = elastic",
    )
    cases = (
        (
            "limit_option0",
            make_cfg(
                fixed_dt=0.0,
                dc=1.0e-6,
                rate_option=0,
                dtheta_max=0.2,
            ),
            "projection_option=1",
        ),
        (
            "limit_steady",
            make_cfg(
                fixed_dt=0.0,
                dc=1.0e-6,
                state_model=0,
                dtheta_max=0.2,
            ),
            "state_var_model=1",
        ),
        (
            "limit_fixed_dt",
            make_cfg(fixed_dt=1.0e-2, dc=1.0e-6, dtheta_max=0.2),
            "requires control.fixed_dt=0",
        ),
        (
            "limit_two",
            make_cfg(fixed_dt=0.0, dc=1.0e-6, dtheta_max=2.0),
            "finite and in [0, 2)",
        ),
        (
            "limit_nan",
            make_cfg(fixed_dt=0.0, dc=1.0e-6, dtheta_max=math.nan),
            "finite and in [0, 2)",
        ),
        ("limit_pt", limit_pt, "control.has_PT=true"),
        (
            "limit_body_force",
            limit_body_force,
            "ic.has_body_force_adjustment=true",
        ),
        (
            "limit_isostasy",
            limit_isostasy,
            "not supported during isostasy adjustment",
        ),
        ("limit_non_rsf", limit_non_rsf, "requires an RSF rheology"),
    )
    for name, cfg, message in cases:
        _, result = run_cfg(
            exe, root, name, cfg, expect_success=False
        )
        if message not in result.stdout:
            raise AssertionError(
                f"{name}: missing diagnostic {message!r}\n{result.stdout}"
            )
    print("[ok] unsupported aging-law timestep combinations are rejected")


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
        check_adaptive_state_limit_disabled(exe, root)
        check_state_rate_limit(exe, root)
        check_state_healing_limit(exe, root)
        check_moving_mesh_rate_refresh(exe, root)
        check_mass_scaling_reference(exe, root)
        check_invalid_option(exe, root)
        check_invalid_state_limits(exe, root)


if __name__ == "__main__":
    main()
