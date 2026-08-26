"""Parser and moving-mesh reset validation for fixed pore-pressure BCs."""

from __future__ import annotations

from pathlib import Path

from fixed_hbc_support import (
    FACE_MASKS,
    assert_fixed_boundary_state,
    render_cfg,
    run_config,
)


RESET_ERROR = "requires an explicit finite bc.hbc_val_* value"
INFINITE_ERROR = "bc.hbc_val_* value must be finite or NaN"


def _assert_rejected_before_simulation(
    exe: Path,
    template: str,
    case_dir: Path,
    *,
    name: str,
    threads: int,
    remeshing_option: int,
    face: str,
    value: float | str | None,
    expected_error: str,
    hydraulic_enabled: bool = True,
    effective_stress_enabled: bool = False,
) -> None:
    values = {} if value is None else {face: value}
    cfg = render_cfg(
        template,
        modelname=name,
        moving_mesh=remeshing_option != 0,
        remeshing_option=remeshing_option,
        hydraulic_enabled=hydraulic_enabled,
        effective_stress_enabled=effective_stress_enabled,
        hbc_types={face: 1},
        hbc_values=values,
    )
    stdout = run_config(
        exe,
        cfg,
        case_dir,
        name,
        threads=threads,
        expected_error=expected_error,
    )
    if "Starting simulation..." in stdout:
        raise AssertionError(f"{name}: invalid config reached the simulation loop")
    if any(case_dir.glob("*.vtkhdf")):
        raise AssertionError(f"{name}: invalid config emitted simulation output")


def _check_finite_face_parser_mapping(
    exe: Path,
    template: str,
    run_root: Path,
    *,
    dimensions: int,
    threads: int,
) -> None:
    physical_faces = ["x0", "x1", "z0", "z1"]
    if dimensions == 3:
        physical_faces[2:2] = ["y0", "y1"]
    values = {
        face: float((index + 1) * 100_000)
        for index, face in enumerate(physical_faces)
    }
    boundary_values = tuple(
        (FACE_MASKS[face], values[face]) for face in physical_faces
    )
    types = {face: 1 for face in physical_faces}

    case_dir = run_root / "finite_face_mapping"
    model = f"fixed_hbc_face_mapping_{dimensions}d"
    cfg = render_cfg(
        template,
        modelname=model,
        has_initial_checkpoint=True,
        hbc_types=types,
        hbc_values=values,
    )
    run_config(exe, cfg, case_dir, "finite_face_mapping", threads=threads)
    assert_fixed_boundary_state(
        case_dir / f"{model}.save.000000.vtkhdf",
        case_dir / f"{model}.chkpt.000000.vtkhdf",
        boundary_values,
        require_interior_dpp=False,
    )
    assert_fixed_boundary_state(
        case_dir / f"{model}.save.000001.vtkhdf",
        case_dir / f"{model}.chkpt.000001.vtkhdf",
        boundary_values,
        require_interior_dpp=False,
    )


def _check_inactive_coupling_parser_compatibility(
    exe: Path,
    template: str,
    run_root: Path,
    *,
    threads: int,
) -> None:
    # Option 12 geometrically resets x0. Its fixed HBC value was nonetheless
    # ignored before pore-pressure mechanics existed, so feature-off configs
    # must continue accepting both the NaN sentinel and legacy infinite data.
    accepted_root = run_root / "coupling_off"
    for label, value in (("nan", None), ("infinite", "inf")):
        values = {} if value is None else {"x0": value}
        name = f"option12_x0_{label}_coupling_off"
        cfg = render_cfg(
            template,
            modelname=name,
            moving_mesh=True,
            remeshing_option=12,
            hydraulic_enabled=False,
            effective_stress_enabled=False,
            hbc_types={"x0": 1},
            hbc_values=values,
        )
        run_config(
            exe,
            cfg,
            accepted_root / name,
            name,
            threads=threads,
        )

    # Either hydraulic transport or explicit effective-stress coupling makes
    # the fixed boundary active, so a reset face can no longer use NaN.
    rejected_root = run_root / "coupling_active"
    for label, hydraulic, effective in (
        ("hydraulic", True, False),
        ("effective_stress", False, True),
    ):
        name = f"option12_x0_nan_{label}"
        _assert_rejected_before_simulation(
            exe,
            template,
            rejected_root / name,
            name=name,
            threads=threads,
            remeshing_option=12,
            face="x0",
            value=None,
            expected_error=RESET_ERROR,
            hydraulic_enabled=hydraulic,
            effective_stress_enabled=effective,
        )


def check_parser_and_reset_mapping(
    exe: Path,
    template: str,
    run_root: Path,
    *,
    dimensions: int,
    threads: int,
) -> None:
    if dimensions not in (2, 3):
        raise AssertionError(f"unsupported test dimension {dimensions}")

    _check_finite_face_parser_mapping(
        exe,
        template,
        run_root,
        dimensions=dimensions,
        threads=threads,
    )
    _check_inactive_coupling_parser_compatibility(
        exe,
        template,
        run_root,
        threads=threads,
    )

    invalid_reset_cases = [
        (1, "z0"),
        (11, "z0"),
        (13, "x0"),
        (13, "x1"),
        (13, "z0"),
    ]
    if dimensions == 2:
        invalid_reset_cases.append((2, "z0"))
    else:
        invalid_reset_cases.extend(((13, "y0"), (13, "y1")))

    invalid_root = run_root / "invalid_reset_mapping"
    for option, face in invalid_reset_cases:
        name = f"option{option}_{face}_requires_finite"
        _assert_rejected_before_simulation(
            exe,
            template,
            invalid_root / name,
            name=name,
            threads=threads,
            remeshing_option=option,
            face=face,
            value=None,
            expected_error=RESET_ERROR,
        )

    # These faces are not geometrically reset by the selected options. Their
    # omitted values must retain the current trace rather than being rejected.
    accepted_cases = [
        (0, "z0"),
        (1, "x0"),
        (11, "z1"),
        (12, "z0"),
        (13, "z1"),
    ]
    if dimensions == 2:
        accepted_cases.append((2, "x1"))
        # Option 13 resets y faces only in 3D; a 2D y key has no physical node.
        accepted_cases.append((13, "y0"))

    accepted_root = run_root / "accepted_reset_mapping"
    for option, face in accepted_cases:
        name = f"option{option}_{face}_nan_allowed"
        cfg = render_cfg(
            template,
            modelname=name,
            moving_mesh=option != 0,
            remeshing_option=option,
            hbc_types={face: 1},
            hbc_values={},
        )
        run_config(
            exe,
            cfg,
            accepted_root / name,
            name,
            threads=threads,
        )

    # Infinite values are never the NaN sentinel, even without moving mesh.
    _assert_rejected_before_simulation(
        exe,
        template,
        run_root / "infinite_value",
        name="nonreset_infinite_value",
        threads=threads,
        remeshing_option=0,
        face="x0",
        value="inf",
        expected_error=INFINITE_ERROR,
    )
