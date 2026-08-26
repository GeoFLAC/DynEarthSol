"""Domain, feature-off, material-order, mesh, and restart hydraulic CFL checks."""

from __future__ import annotations

import sys
from pathlib import Path

from hydraulic_cfl_support import (
    CflCase,
    assert_close,
    parse_dt_records,
    render_cfg,
    replace_once,
    run_case,
    run_config,
)


CASES = (
    CflCase("enabled_low", True, (1.0e-16, 1.0e-14), 1.0e-4),
    # Put the largest diffusivity in material 0 while mattype_ref is 1. This
    # distinguishes a true domain-wide reduction from sampling only the
    # reference material.
    CflCase("enabled_high", True, (1.0e-8, 1.0e-14), 1.0e2),
    CflCase("disabled_low", False, (1.0e-16, 1.0e-14), None),
    CflCase("disabled_high", False, (1.0e-14, 1.0e-8), None),
    CflCase(
        "reference_order_01",
        True,
        (1.0e-16, 1.0e-8),
        1.0e2,
        quasi_static=True,
        expect_selected_hydraulic=False,
    ),
    CflCase(
        "reference_order_10",
        True,
        (1.0e-16, 1.0e-8),
        1.0e2,
        quasi_static=True,
        layer_mattypes=(1, 0),
        expect_selected_hydraulic=False,
    ),
    CflCase(
        "moving_mesh_refresh",
        True,
        (1.0e-8, 1.0e-14),
        1.0e2,
        moving_mesh=True,
        minimum_dt_records=3,
    ),
)


def check_enabled(case: CflCase, records: list[dict[str, float]]) -> None:
    assert case.expected_diffusivity_max is not None
    if len(records) < case.minimum_dt_records:
        raise AssertionError(
            f"{case.name}: expected at least {case.minimum_dt_records} "
            f"compute_dt records, found {len(records)}"
        )
    for index, record in enumerate(records):
        expected_dt = (
            0.5
            * record["min_length"]
            * record["min_length"]
            / case.expected_diffusivity_max
        )
        assert_close(
            record["hydro_diffusion"],
            expected_dt,
            f"{case.name}: hydraulic CFL record {index}",
        )
        if case.expect_selected_hydraulic:
            assert_close(
                record["selected"],
                expected_dt,
                f"{case.name}: selected timestep record {index}",
            )


def check_disabled(
    low: list[dict[str, float]], high: list[dict[str, float]]
) -> None:
    if len(low) != len(high):
        raise AssertionError("feature-off cases produced different compute_dt counts")
    for index, (low_record, high_record) in enumerate(zip(low, high)):
        if low_record.keys() != high_record.keys():
            raise AssertionError(f"feature-off record {index} has different fields")
        for name in low_record:
            assert_close(
                low_record[name],
                high_record[name],
                f"feature-off permeability leaked into {name}, record {index}",
            )
        if low_record["hydro_diffusion"] != sys.float_info.max:
            raise AssertionError(
                f"feature-off record {index}: hydraulic bound is not disabled"
            )


def check_same_records(
    first_name: str,
    first: list[dict[str, float]],
    second_name: str,
    second: list[dict[str, float]],
) -> None:
    if len(first) != len(second):
        raise AssertionError(f"{first_name}/{second_name}: different compute_dt counts")
    for index, (first_record, second_record) in enumerate(zip(first, second)):
        if first_record.keys() != second_record.keys():
            raise AssertionError(
                f"{first_name}/{second_name}: record {index} fields differ"
            )
        for name in first_record:
            assert_close(
                first_record[name],
                second_record[name],
                f"{first_name}/{second_name}: {name} depends on element ordering, record {index}",
            )


def check_restart_refresh(exe: Path, template: str, run_root: Path) -> None:
    restart_dir = run_root / "restart_refresh"
    restart_dir.mkdir(parents=True, exist_ok=True)
    source_case = CflCase(
        "restart_source", True, (1.0e-14, 1.0e-8), 1.0e2
    )
    source_cfg = render_cfg(template, source_case)
    source_cfg = replace_once(source_cfg, "fixed_dt = 0\n", "fixed_dt = 1.2345\n")
    source_cfg = replace_once(
        source_cfg,
        "has_initial_checkpoint = no\n",
        "has_initial_checkpoint = yes\n",
    )
    source_cfg = replace_once(
        source_cfg, "has_marker_output = no\n", "has_marker_output = yes\n"
    )
    source_path = restart_dir / "source.cfg"
    source_path.write_text(source_cfg, encoding="ascii")
    run_config(exe, source_path, restart_dir, source_case.name)

    target_case = CflCase(
        "restart_target", True, (1.0e-14, 1.0e-8), 1.0e2
    )
    target_cfg = render_cfg(template, target_case)
    target_cfg = replace_once(
        target_cfg,
        "modelname = hydraulic_cfl_restart_target\n",
        "modelname = hydraulic_cfl_restart_target\n"
        "is_restarting = yes\n"
        "restarting_from_modelname = hydraulic_cfl_restart_source\n"
        "restarting_from_frame = 0\n",
    )
    target_path = restart_dir / "restart.cfg"
    target_path.write_text(target_cfg, encoding="ascii")
    restart_stdout = run_config(exe, target_path, restart_dir, target_case.name)
    restart_records = parse_dt_records(restart_stdout, target_case.name)
    check_enabled(target_case, restart_records)


def run_base_regressions(exe: Path, template: str, run_root: Path) -> None:
    results = {
        case.name: run_case(exe, template, run_root, case) for case in CASES
    }
    for case in CASES:
        if case.expected_diffusivity_max is not None:
            check_enabled(case, results[case.name])
    check_disabled(results["disabled_low"], results["disabled_high"])
    check_same_records(
        "reference_order_01",
        results["reference_order_01"],
        "reference_order_10",
        results["reference_order_10"],
    )
    check_restart_refresh(exe, template, run_root)
