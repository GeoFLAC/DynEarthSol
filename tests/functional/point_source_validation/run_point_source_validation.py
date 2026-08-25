#!/usr/bin/env python3

"""Exercise fluid point-source input, lookup, and nodal conservation."""

from __future__ import annotations

import argparse
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
DEFAULT_TEMPLATE = HERE / "point_source_validation_base.cfg"
DEFAULT_RUN_DIR = HERE / "runs"


@dataclass(frozen=True)
class Case:
    name: str
    source_section: str
    expected_code: int = 0
    expected_text: str = ""
    meshing_option: int = 1
    fixed_dt: str = "1"
    has_pt: str = "no"
    has_moving_mesh: str = "no"
    hbc_x1: str = "0"
    porosity: str = "0.2"
    permeability: str = "1e-20"
    biot_coeff: str = "1"
    fluid_bulk_modulus: str = "2.17e9"


def injection_section(
    *,
    num_points: int = 1,
    points_x: str = "[0.5]",
    points_y: str = "[]",
    points_z: str = "[-0.5]",
    points_unit: str = "m",
    rate_model: str = "constant_rate",
    rate: str = "[1e-6]",
    total_amount: str = "[0]",
    start: str = "[0]",
    end: str = "[]",
    duration: str = "[]",
) -> str:
    return f"""
[injection]
enabled = yes
num_points = {num_points}
points_x = {points_x}
points_y = {points_y}
points_z = {points_z}
points_unit = {points_unit}
rate_model = {rate_model}
rate = {rate}
total_amount = {total_amount}
start_time_in_yr = {start}
end_time_in_yr = {end}
duration_in_yr = {duration}
""".strip()


CASES_2D = (
    Case(
        "box_boundary_tolerance",
        injection_section(points_x="[1.0000000000005]"),
    ),
    Case(
        "legacy_points_y_alias",
        injection_section(points_y="[-0.5]", points_z="[]"),
    ),
    Case(
        "shared_interior_edge",
        injection_section(points_x="[0.5]", points_z="[-0.5]"),
    ),
    Case(
        "multi_point_conservation",
        injection_section(
            num_points=2,
            points_x="[0.25, 0.75]",
            points_z="[-0.25, -0.75]",
            rate="[1e-6, -2e-7]",
            total_amount="[0, 0]",
            start="[0, 0]",
        ),
    ),
    Case(
        "fixed_pressure_source_rejected",
        injection_section(points_x="[0.999999]", points_z="[-0.5]"),
        51,
        "fixed pore-pressure node",
        hbc_x1="1",
    ),
    Case(
        "nonpositive_hydraulic_mass_rejected",
        injection_section(),
        51,
        "non-positive hydraulic mass",
        porosity="0",
        permeability="0",
        biot_coeff="0",
    ),
    Case(
        "nonfinite_hydraulic_mass_rejected",
        injection_section(),
        50,
        "non-finite hydraulic mass",
        fluid_bulk_modulus="0",
    ),
    Case(
        "pt_moving_mesh_fails_closed",
        injection_section(),
        11,
        "PT-driven moving-mesh remeshing is not yet supported",
        has_pt="yes",
        has_moving_mesh="yes",
    ),
    Case(
        "nonpositive_step_interval_rejected",
        injection_section(),
        50,
        "hydraulic integration interval must be finite and positive",
        fixed_dt="-1",
    ),
    Case(
        "missing_2d_vertical_coordinate",
        injection_section(points_y="[]", points_z="[]"),
        11,
        "requires points_z in 2D",
    ),
    Case(
        "ambiguous_2d_vertical_coordinate",
        injection_section(points_y="[-0.5]", points_z="[-0.5]"),
        11,
        "mutually exclusive in 2D",
    ),
    Case(
        "outside_box_mesh",
        injection_section(points_x="[1.01]"),
        11,
        "outside the box-mesh interval",
    ),
    Case(
        "coordinate_unit_overflow",
        injection_section(points_x="[1e308]", points_unit="km"),
        11,
        "must be finite after unit conversion",
    ),
    Case(
        "nan_rate",
        injection_section(rate="[nan]"),
        11,
        "incorrect format for injection.rate",
    ),
    Case(
        "infinite_total_amount",
        injection_section(total_amount="[inf]"),
        11,
        "incorrect format for injection.total_amount",
    ),
    Case(
        "nonfinite_end_time",
        injection_section(end="[nan]"),
        11,
        "incorrect format for injection.end_time_in_yr",
    ),
    Case(
        "zero_duration",
        injection_section(duration="[0]"),
        11,
        "must be finite and greater than zero",
    ),
    Case(
        "start_time_conversion_overflow",
        injection_section(start="[1e308]"),
        11,
        "overflows when converted to seconds",
    ),
    Case(
        "end_time_conversion_overflow",
        injection_section(end="[1e308]"),
        11,
        "end time must remain finite",
    ),
    # The point is inside the nominal box but outside triangle.poly. Its
    # inactive schedule still has to be mapped against the current mesh.
    Case(
        "inactive_source_outside_actual_mesh",
        injection_section(points_x="[0.9]", points_z="[-0.9]", start="[100]"),
        51,
        "could not be mapped into the current mesh",
        meshing_option=90,
    ),
    Case(
        "actual_mesh_sloping_edge",
        injection_section(points_x="[0.5]", points_z="[-0.5]"),
        meshing_option=90,
    ),
)


CASES_3D = (
    Case(
        "box_boundary_tolerance_3d",
        injection_section(
            points_x="[1.00000000002]",
            points_y="[0.5]",
            points_z="[-0.5]",
        ),
    ),
    Case(
        "shared_internal_face_3d",
        injection_section(
            points_x="[0.5]",
            points_y="[0.5]",
            points_z="[-0.5]",
        ),
    ),
    Case(
        "multi_point_conservation_3d",
        injection_section(
            num_points=2,
            points_x="[0.25, 0.75]",
            points_y="[0.25, 0.75]",
            points_z="[-0.25, -0.75]",
            rate="[1e-6, -2e-7]",
            total_amount="[0, 0]",
            start="[0, 0]",
        ),
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe-2d", type=Path, default=REPO_ROOT / "dynearthsol2d")
    parser.add_argument(
        "--exe-3d",
        type=Path,
        help="Optional 3D executable; when supplied, 3D cases also run.",
    )
    parser.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--clean", action="store_true")
    return parser.parse_args()


def render_config(template: str, case: Case) -> str:
    replacements = {
        "__MODELNAME__": case.name,
        "__MESHING_OPTION__": str(case.meshing_option),
        "__FIXED_DT__": case.fixed_dt,
        "__HAS_PT__": case.has_pt,
        "__HAS_MOVING_MESH__": case.has_moving_mesh,
        "__HBC_X1__": case.hbc_x1,
        "__POROSITY__": case.porosity,
        "__PERMEABILITY__": case.permeability,
        "__BIOT_COEFF__": case.biot_coeff,
        "__FLUID_BULK_MODULUS__": case.fluid_bulk_modulus,
        "__SOURCE_SECTION__": case.source_section,
    }
    for key, value in replacements.items():
        template = template.replace(key, value)
    return template


def run_case(exe: Path, template: str, run_root: Path, case: Case) -> None:
    case_dir = run_root / case.name
    case_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(HERE / "triangle.poly", case_dir / "mesh.poly")
    cfg_path = case_dir / "input.cfg"
    cfg_path.write_text(render_config(template, case), encoding="ascii")

    completed = subprocess.run(
        [str(exe), str(cfg_path)],
        cwd=case_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=60,
        check=False,
    )
    (case_dir / "run.log").write_text(completed.stdout, encoding="utf-8")

    if completed.returncode != case.expected_code:
        raise AssertionError(
            f"{case.name}: exit {completed.returncode}, expected {case.expected_code}; "
            f"see {case_dir / 'run.log'}"
        )
    if case.expected_text and case.expected_text not in completed.stdout:
        raise AssertionError(
            f"{case.name}: expected diagnostic {case.expected_text!r}; "
            f"see {case_dir / 'run.log'}"
        )
    print(f"PASS {case.name} (exit {completed.returncode})", flush=True)


def main() -> None:
    args = parse_args()
    run_root = args.run_dir.resolve()
    if args.clean and run_root.exists():
        shutil.rmtree(run_root)
    run_root.mkdir(parents=True, exist_ok=True)
    template = args.template.resolve().read_text(encoding="ascii")

    exe_2d = args.exe_2d.resolve()
    for case in CASES_2D:
        run_case(exe_2d, template, run_root / "2d", case)

    if args.exe_3d is not None:
        exe_3d = args.exe_3d.resolve()
        for case in CASES_3D:
            run_case(exe_3d, template, run_root / "3d", case)


if __name__ == "__main__":
    main()
