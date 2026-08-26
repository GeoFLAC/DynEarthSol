"""Shared configuration, execution, and HDF5 helpers for fixed-HBC tests."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Mapping, Sequence

import h5py


BOUNDX0 = 1
BOUNDX1 = 2
BOUNDY0 = 4
BOUNDY1 = 8
BOUNDZ0 = 16
BOUNDZ1 = 32

FACE_MASKS = {
    "x0": BOUNDX0,
    "x1": BOUNDX1,
    "y0": BOUNDY0,
    "y1": BOUNDY1,
    "z0": BOUNDZ0,
    "z1": BOUNDZ1,
}
FACES = tuple(FACE_MASKS)

PORE_PRESSURE = "pore pressure"
PENDING_DPP = "pore pressure stress increment"
BCFLAG = "bcflag"
COORDINATE = "coordinate"

# Root datasets are VDS views. Mutate their physical PointData sources so the
# restart reader observes the deliberate corruption through the root VDS.
POINT_DATA = "/VTKHDF/grid/PointData"

FRESH_X0 = 100_000.0
FRESH_Z0 = 300_000.0
RESTART_X0 = 400_000.0
RESTART_Z0 = 600_000.0
REMESH_BOTTOM = 500_000.0
CORRUPT_PRESSURE = 9_000_000.0
CORRUPT_DPP = 700_000_000.0

FRESH_BOUNDARIES = ((BOUNDX0, FRESH_X0), (BOUNDZ0, FRESH_Z0))
RESTART_BOUNDARIES = ((BOUNDX0, RESTART_X0), (BOUNDZ0, RESTART_Z0))
FRESH_MASKS = (BOUNDX0, BOUNDZ0)
REMESH_BOUNDARIES = ((BOUNDZ0, REMESH_BOTTOM),)


def _bool(value: bool) -> str:
    return "yes" if value else "no"


def _number(value: float | str) -> str:
    if isinstance(value, str):
        return value
    return f"{value:.17e}"


def render_cfg(
    template: str,
    *,
    modelname: str,
    max_steps: int = 1,
    output_step_interval: int = 1,
    has_initial_checkpoint: bool = False,
    is_restarting: bool = False,
    restarting_from_modelname: str | None = None,
    moving_mesh: bool = False,
    quality_check_step_interval: int = 1000,
    min_quality: float = 0.4,
    has_output_during_remeshing: bool = False,
    remeshing_option: int = 0,
    vbc_z0: float = 0.0,
    hydraulic_enabled: bool = True,
    effective_stress_enabled: bool = False,
    hbc_types: Mapping[str, int] | None = None,
    hbc_values: Mapping[str, float | str] | None = None,
) -> str:
    """Fill the fixed-HBC configuration template."""

    hbc_types = dict(hbc_types or {})
    hbc_values = dict(hbc_values or {})
    unknown_faces = (set(hbc_types) | set(hbc_values)) - set(FACES)
    if unknown_faces:
        raise AssertionError(f"unknown hydraulic boundary faces: {unknown_faces}")

    hbc_lines = [f"hbc_{face} = {hbc_types.get(face, 0)}" for face in FACES]
    hbc_lines.extend(
        f"hbc_val_{face} = {_number(hbc_values[face])}"
        for face in FACES
        if face in hbc_values
    )

    replacements = {
        "__MODELNAME__": modelname,
        "__MAX_STEPS__": str(max_steps),
        "__OUTPUT_STEP_INTERVAL__": str(output_step_interval),
        "__HAS_INITIAL_CHECKPOINT__": _bool(has_initial_checkpoint),
        "__IS_RESTARTING__": _bool(is_restarting),
        "__RESTARTING_FROM_MODELNAME__": (
            restarting_from_modelname or modelname
        ),
        "__HAS_OUTPUT_DURING_REMESHING__": _bool(
            has_output_during_remeshing
        ),
        "__HAS_MOVING_MESH__": _bool(moving_mesh),
        "__QUALITY_CHECK_STEP_INTERVAL__": str(quality_check_step_interval),
        "__MIN_QUALITY__": f"{min_quality:.17e}",
        "__REMESHING_OPTION__": str(remeshing_option),
        "__VBC_VALUE_Z0__": f"{vbc_z0:.17e}",
        "__HYDRAULIC_ENABLED__": _bool(hydraulic_enabled),
        "__EFFECTIVE_STRESS_ENABLED__": _bool(effective_stress_enabled),
        "__HBC_CONFIG__": "\n".join(hbc_lines),
    }
    rendered = template
    for token, value in replacements.items():
        if token not in rendered:
            raise AssertionError(f"missing config token {token}")
        rendered = rendered.replace(token, value)
    if "__" in rendered:
        raise AssertionError("unexpanded config token remains")
    return rendered


def run_config(
    exe: Path,
    cfg: str,
    case_dir: Path,
    name: str,
    *,
    threads: int,
    expected_error: str | None = None,
) -> str:
    case_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = case_dir / f"{name}.cfg"
    cfg_path.write_text(cfg, encoding="ascii")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    result = subprocess.run(
        [str(exe), str(cfg_path)],
        cwd=case_dir,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    (case_dir / f"{name}.log").write_text(result.stdout, encoding="utf-8")
    if expected_error is not None:
        if result.returncode == 0:
            raise AssertionError(f"{name}: invalid config unexpectedly succeeded")
        if expected_error not in result.stdout:
            print(result.stdout, file=sys.stderr)
            raise AssertionError(
                f"{name}: missing expected config error {expected_error!r}"
            )
        return result.stdout
    if result.returncode != 0:
        print(result.stdout, file=sys.stderr)
        raise RuntimeError(f"{name}: DynEarthSol exited with {result.returncode}")
    return result.stdout


def read_root_array(path: Path, name: str):
    if not path.is_file():
        raise AssertionError(f"missing HDF5 output: {path}")
    with h5py.File(path, "r") as handle:
        if name not in handle:
            raise AssertionError(f"{path}: missing root dataset {name!r}")
        return handle[name][...]


def boundary_indices(bcflag, mask: int) -> list[int]:
    indices = [i for i, flag in enumerate(bcflag) if int(flag) & mask]
    if not indices:
        raise AssertionError(f"fixture has no nodes carrying boundary mask {mask}")
    return indices


def fixed_indices(bcflag, masks: Sequence[int]) -> set[int]:
    fixed: set[int] = set()
    for mask in masks:
        fixed.update(boundary_indices(bcflag, mask))
    return fixed


def expected_fixed_pressures(bcflag, boundary_values) -> dict[int, float]:
    for mask, _ in boundary_values:
        boundary_indices(bcflag, mask)

    expected = {}
    for node, flag_value in enumerate(bcflag):
        values = [
            value for mask, value in boundary_values if int(flag_value) & mask
        ]
        if values:
            expected[node] = sum(values) / len(values)
    return expected


def assert_pressure_mapping(save_path: Path, expected: Mapping[int, float]) -> None:
    pressure = read_root_array(save_path, PORE_PRESSURE)
    for node, expected_value in expected.items():
        actual = float(pressure[node])
        if actual != expected_value:
            raise AssertionError(
                f"{save_path}: fixed node {node} pressure={actual:.17e}, "
                f"expected exact {expected_value:.17e}"
            )


def assert_fixed_pressure(save_path: Path, boundary_values) -> set[int]:
    bcflag = read_root_array(save_path, BCFLAG)
    pressure = read_root_array(save_path, PORE_PRESSURE)
    if len(bcflag) != len(pressure):
        raise AssertionError("bcflag and pore-pressure sizes differ")

    expected = expected_fixed_pressures(bcflag, boundary_values)
    assert_pressure_mapping(save_path, expected)
    return set(expected)


def assert_boundary_intersection_average(
    save_path: Path, mask_a: int, mask_b: int, expected: float
) -> None:
    bcflag = read_root_array(save_path, BCFLAG)
    pressure = read_root_array(save_path, PORE_PRESSURE)
    corners = [
        node
        for node, flag in enumerate(bcflag)
        if int(flag) & mask_a and int(flag) & mask_b
    ]
    if not corners:
        raise AssertionError(
            f"{save_path}: fixture has no nodes at mask intersection "
            f"{mask_a}|{mask_b}"
        )
    for node in corners:
        actual = float(pressure[node])
        if actual != expected:
            raise AssertionError(
                f"{save_path}: corner node {node} pressure={actual:.17e}, "
                f"expected average {expected:.17e}"
            )


def assert_fixed_dpp(
    save_path: Path,
    checkpoint_path: Path,
    fixed: set[int],
    expected: float,
) -> None:
    bcflag = read_root_array(save_path, BCFLAG)
    dpp = read_root_array(checkpoint_path, PENDING_DPP)
    if len(bcflag) != len(dpp):
        raise AssertionError("bcflag and dpp sizes differ")
    for node in fixed:
        actual = float(dpp[node])
        if actual != expected:
            raise AssertionError(
                f"{checkpoint_path}: fixed node {node} dpp={actual:.17e}, "
                f"expected {expected:.17e}"
            )


def assert_zero_fixed_dpp(
    save_path: Path, checkpoint_path: Path, fixed: set[int]
) -> None:
    assert_fixed_dpp(save_path, checkpoint_path, fixed, 0.0)


def assert_fixed_boundary_state(
    save_path: Path,
    checkpoint_path: Path,
    boundary_values,
    *,
    require_interior_dpp: bool,
) -> None:
    fixed = assert_fixed_pressure(save_path, boundary_values)
    assert_zero_fixed_dpp(save_path, checkpoint_path, fixed)

    if require_interior_dpp:
        dpp = read_root_array(checkpoint_path, PENDING_DPP)
        interior_max = max(
            (
                abs(float(value))
                for node, value in enumerate(dpp)
                if node not in fixed
            ),
            default=0.0,
        )
        if interior_max == 0.0:
            raise AssertionError(
                f"{checkpoint_path}: fixture generated no interior pressure increment"
            )


def corrupt_restart_boundary(
    save_path: Path, checkpoint_path: Path, masks: Sequence[int]
) -> dict[int, float]:
    bcflag = read_root_array(save_path, BCFLAG)
    fixed = fixed_indices(bcflag, masks)
    held_pressure = {
        node: CORRUPT_PRESSURE + 1024.0 * node for node in fixed
    }

    pressure_storage = f"{POINT_DATA}/{PORE_PRESSURE}"
    with h5py.File(save_path, "r+") as handle:
        if pressure_storage not in handle:
            raise AssertionError(f"{save_path}: missing storage {pressure_storage!r}")
        pressure = handle[pressure_storage][...]
        for node, value in held_pressure.items():
            pressure[node] = value
        handle[pressure_storage][...] = pressure

    dpp_storage = f"{POINT_DATA}/{PENDING_DPP}"
    with h5py.File(checkpoint_path, "r+") as handle:
        if dpp_storage not in handle:
            raise AssertionError(
                f"{checkpoint_path}: missing storage {dpp_storage!r}"
            )
        dpp = handle[dpp_storage][...]
        for node in fixed:
            dpp[node] = CORRUPT_DPP
        handle[dpp_storage][...] = dpp

    # Confirm the root VDS read by DES exposes the intended corruption.
    assert_pressure_mapping(save_path, held_pressure)
    root_dpp = read_root_array(checkpoint_path, PENDING_DPP)
    if any(float(root_dpp[node]) != CORRUPT_DPP for node in fixed):
        raise AssertionError("corrupted dpp is not visible through root VDS")
    return held_pressure


def assert_bottom_was_relocated(pre_remesh: Path, post_remesh: Path) -> None:
    pre_flags = read_root_array(pre_remesh, BCFLAG)
    post_flags = read_root_array(post_remesh, BCFLAG)
    pre_coord = read_root_array(pre_remesh, COORDINATE)
    post_coord = read_root_array(post_remesh, COORDINATE)
    pre_bottom = boundary_indices(pre_flags, BOUNDZ0)
    post_bottom = boundary_indices(post_flags, BOUNDZ0)

    pre_z = [float(pre_coord[node][-1]) for node in pre_bottom]
    post_z = [float(post_coord[node][-1]) for node in post_bottom]
    if min(pre_z) <= -0.95:
        raise AssertionError(
            f"{pre_remesh}: bottom did not move far enough to test interpolation"
        )
    if max(abs(value + 1.0) for value in post_z) > 1.0e-12:
        raise AssertionError(
            f"{post_remesh}: option-1 bottom was not restored to z=-1"
        )
