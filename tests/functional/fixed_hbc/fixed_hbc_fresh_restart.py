"""Fresh-run and restart fixed pore-pressure boundary checks."""

from __future__ import annotations

from pathlib import Path

from fixed_hbc_support import (
    BCFLAG,
    BOUNDX0,
    BOUNDZ0,
    CORRUPT_DPP,
    FRESH_BOUNDARIES,
    FRESH_MASKS,
    FRESH_X0,
    FRESH_Z0,
    RESTART_BOUNDARIES,
    RESTART_X0,
    RESTART_Z0,
    assert_boundary_intersection_average,
    assert_fixed_boundary_state,
    assert_fixed_dpp,
    assert_fixed_pressure,
    assert_pressure_mapping,
    assert_zero_fixed_dpp,
    corrupt_restart_boundary,
    fixed_indices,
    read_root_array,
    render_cfg,
    run_config,
)


FINITE_TYPES = {"x0": 1, "z0": 1}


def check_fresh_and_restart(
    exe: Path, template: str, run_root: Path, *, threads: int
) -> None:
    case_dir = run_root / "fresh_restart"
    source_model = "fixed_hbc_source"
    fresh_cfg = render_cfg(
        template,
        modelname=source_model,
        has_initial_checkpoint=True,
        hbc_types=FINITE_TYPES,
        hbc_values={"x0": FRESH_X0, "z0": FRESH_Z0},
    )
    run_config(exe, fresh_cfg, case_dir, "fresh", threads=threads)

    fresh_initial = case_dir / f"{source_model}.save.000000.vtkhdf"
    assert_fixed_boundary_state(
        fresh_initial,
        case_dir / f"{source_model}.chkpt.000000.vtkhdf",
        FRESH_BOUNDARIES,
        require_interior_dpp=False,
    )
    assert_boundary_intersection_average(
        fresh_initial, BOUNDX0, BOUNDZ0, (FRESH_X0 + FRESH_Z0) / 2.0
    )
    assert_fixed_boundary_state(
        case_dir / f"{source_model}.save.000001.vtkhdf",
        case_dir / f"{source_model}.chkpt.000001.vtkhdf",
        FRESH_BOUNDARIES,
        require_interior_dpp=True,
    )

    held_pressure = corrupt_restart_boundary(
        fresh_initial,
        case_dir / f"{source_model}.chkpt.000000.vtkhdf",
        FRESH_MASKS,
    )

    restart_model = "fixed_hbc_restart"
    restart_cfg = render_cfg(
        template,
        modelname=restart_model,
        is_restarting=True,
        restarting_from_modelname=source_model,
        hbc_types=FINITE_TYPES,
        hbc_values={"x0": RESTART_X0, "z0": RESTART_Z0},
    )
    run_config(exe, restart_cfg, case_dir, "restart", threads=threads)

    # The initial restart save is emitted before a new physical step. It must
    # already reflect the new finite config values, not the corrupted source.
    restart_initial = case_dir / f"{restart_model}.save.000000.vtkhdf"
    assert_fixed_pressure(restart_initial, RESTART_BOUNDARIES)
    assert_boundary_intersection_average(
        restart_initial,
        BOUNDX0,
        BOUNDZ0,
        (RESTART_X0 + RESTART_Z0) / 2.0,
    )
    assert_fixed_boundary_state(
        case_dir / f"{restart_model}.save.000001.vtkhdf",
        case_dir / f"{restart_model}.chkpt.000001.vtkhdf",
        RESTART_BOUNDARIES,
        require_interior_dpp=True,
    )

    # With every active fixed face left at its NaN default, restart must retain
    # the current pressure trace exactly while still clearing stale dpp state.
    nan_model = "fixed_hbc_nan_hold"
    nan_cfg = render_cfg(
        template,
        modelname=nan_model,
        is_restarting=True,
        restarting_from_modelname=source_model,
        hbc_types=FINITE_TYPES,
        hbc_values={},
    )
    run_config(exe, nan_cfg, case_dir, "nan_hold", threads=threads)

    nan_initial = case_dir / f"{nan_model}.save.000000.vtkhdf"
    nan_step = case_dir / f"{nan_model}.save.000001.vtkhdf"
    nan_checkpoint = case_dir / f"{nan_model}.chkpt.000001.vtkhdf"
    assert_pressure_mapping(nan_initial, held_pressure)
    assert_pressure_mapping(nan_step, held_pressure)
    bcflag = read_root_array(nan_step, BCFLAG)
    fixed = fixed_indices(bcflag, FRESH_MASKS)
    assert_zero_fixed_dpp(nan_step, nan_checkpoint, fixed)

    # HBC configuration was historically inert when neither transport nor
    # effective-stress coupling was active. Finite values must therefore not
    # rewrite restored pressure or dpp in that feature-off mode.
    off_model = "fixed_hbc_coupling_off"
    off_cfg = render_cfg(
        template,
        modelname=off_model,
        is_restarting=True,
        restarting_from_modelname=source_model,
        hydraulic_enabled=False,
        effective_stress_enabled=False,
        hbc_types=FINITE_TYPES,
        hbc_values={"x0": RESTART_X0, "z0": RESTART_Z0},
    )
    run_config(exe, off_cfg, case_dir, "coupling_off", threads=threads)

    off_initial = case_dir / f"{off_model}.save.000000.vtkhdf"
    off_step = case_dir / f"{off_model}.save.000001.vtkhdf"
    off_checkpoint = case_dir / f"{off_model}.chkpt.000001.vtkhdf"
    assert_pressure_mapping(off_initial, held_pressure)
    assert_pressure_mapping(off_step, held_pressure)
    off_flags = read_root_array(off_step, BCFLAG)
    off_fixed = fixed_indices(off_flags, FRESH_MASKS)
    assert_fixed_dpp(off_step, off_checkpoint, off_fixed, CORRUPT_DPP)
