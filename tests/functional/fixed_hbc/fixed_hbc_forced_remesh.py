"""Forced option-1 remesh fixed pore-pressure boundary lifecycle check."""

from __future__ import annotations

from pathlib import Path

from fixed_hbc_support import (
    REMESH_BOTTOM,
    REMESH_BOUNDARIES,
    assert_bottom_was_relocated,
    assert_fixed_pressure,
    render_cfg,
    run_config,
)


def check_forced_remesh(
    exe: Path, template: str, run_root: Path, *, threads: int
) -> None:
    case_dir = run_root / "forced_remesh"
    model = "fixed_hbc_remesh"
    cfg = render_cfg(
        template,
        modelname=model,
        max_steps=1,
        output_step_interval=2,
        has_initial_checkpoint=True,
        moving_mesh=True,
        quality_check_step_interval=1,
        min_quality=0.999999,
        has_output_during_remeshing=True,
        remeshing_option=1,
        # Move the bottom upward on the only physical step. Option 1 then
        # restores it to z=-zlength, so interpolated pressure alone cannot be
        # relied on for the prescribed bottom trace.
        vbc_z0=10.0,
        hbc_types={"z0": 1},
        hbc_values={"z0": REMESH_BOTTOM},
    )
    stdout = run_config(exe, cfg, case_dir, "forced_remesh", threads=threads)
    if stdout.count("Remeshing starts") != 1:
        raise AssertionError(
            "forced-remesh fixture expected exactly one remesh, found "
            f"{stdout.count('Remeshing starts')}"
        )

    # With output_step_interval=2, frame 1 is immediately before remeshing and
    # frame 2 immediately after it. The post-remesh frame is emitted before a
    # later hydraulic update could repair a missing lifecycle call.
    pre_remesh = case_dir / f"{model}.save.000001.vtkhdf"
    post_remesh = case_dir / f"{model}.save.000002.vtkhdf"
    assert_bottom_was_relocated(pre_remesh, post_remesh)
    assert_fixed_pressure(pre_remesh, REMESH_BOUNDARIES)
    assert_fixed_pressure(post_remesh, REMESH_BOUNDARIES)
