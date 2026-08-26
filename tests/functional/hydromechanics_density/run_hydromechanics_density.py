#!/usr/bin/env python3

"""Regress reference fluid density, hydrostatic pressure, and water loading."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np


sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
TEMPLATE = HERE / "hydromechanics_density_base.cfg"
SOLID_DENSITY = 2000.0
GRAVITY = 10.0


@dataclass(frozen=True)
class Case:
    name: str
    hydraulic: str
    porosity: float
    fluid_density: float
    surf_base_level: float = 0.0
    water_loading: str = "no"
    sea_water_density: float = 1200.0
    stress_bc_z1: int = 3
    stress_val_z1: float = -1000.0
    ref_pressure_option: int = 0
    num_materials: int = 1


@dataclass(frozen=True)
class RunResult:
    initial: dict[str, np.ndarray]
    final: dict[str, np.ndarray]
    case_dir: Path
    modelname: str


CASES = (
    Case(
        "hydrostatic_water",
        hydraulic="yes",
        porosity=0.25,
        fluid_density=500.0,
        surf_base_level=0.25,
        water_loading="yes",
        stress_bc_z1=0,
        stress_val_z1=0.0,
    ),
    Case(
        "signed_pore_datum",
        hydraulic="yes",
        porosity=0.25,
        fluid_density=500.0,
        surf_base_level=-0.25,
        stress_bc_z1=0,
        stress_val_z1=0.0,
    ),
    Case("mass_fluid_500", "yes", 0.25, 500.0),
    Case("mass_fluid_1500", "yes", 0.25, 1500.0),
    Case("feature_off_a", "no", 0.10, 100.0),
    Case("feature_off_b", "no", 0.90, 9000.0),
    Case(
        "prem_without_water",
        hydraulic="yes",
        porosity=0.25,
        fluid_density=500.0,
        surf_base_level=0.25,
        stress_bc_z1=0,
        stress_val_z1=0.0,
        ref_pressure_option=1,
        num_materials=2,
    ),
    Case(
        "prem_with_water",
        hydraulic="yes",
        porosity=0.25,
        fluid_density=500.0,
        surf_base_level=0.25,
        water_loading="yes",
        stress_bc_z1=0,
        stress_val_z1=0.0,
        ref_pressure_option=1,
        num_materials=2,
    ),
    Case(
        "prem_modified_without_water",
        hydraulic="yes",
        porosity=0.25,
        fluid_density=500.0,
        surf_base_level=0.25,
        stress_bc_z1=0,
        stress_val_z1=0.0,
        ref_pressure_option=2,
        num_materials=2,
    ),
    Case(
        "prem_modified_with_water",
        hydraulic="yes",
        porosity=0.25,
        fluid_density=500.0,
        surf_base_level=0.25,
        water_loading="yes",
        stress_bc_z1=0,
        stress_val_z1=0.0,
        ref_pressure_option=2,
        num_materials=2,
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe-2d", type=Path, default=REPO_ROOT / "dynearthsol2d")
    parser.add_argument("--exe-3d", type=Path)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 4])
    parser.add_argument("--run-dir", type=Path)
    return parser.parse_args()


def render_config(
    template: str,
    case: Case,
    modelname: str,
    *,
    restarting: bool = False,
    restart_frame: int = 0,
) -> str:
    replacements = {
        "__MODELNAME__": modelname,
        "__MAX_STEPS__": "1",
        "__IS_RESTARTING__": "yes" if restarting else "no",
        "__RESTART_MODEL__": modelname,
        "__RESTART_FRAME__": str(restart_frame),
        "__HYDRAULIC__": case.hydraulic,
        "__POROSITY__": str(case.porosity),
        "__FLUID_DENSITY__": str(case.fluid_density),
        "__SURF_BASE_LEVEL__": str(case.surf_base_level),
        "__WATER_LOADING__": case.water_loading,
        "__SEA_WATER_DENSITY__": str(case.sea_water_density),
        "__STRESS_BC_Z1__": str(case.stress_bc_z1),
        "__STRESS_VAL_Z1__": str(case.stress_val_z1),
        "__REF_PRESSURE_OPTION__": str(case.ref_pressure_option),
        "__NUM_MATERIALS__": str(case.num_materials),
    }
    for key, value in replacements.items():
        template = template.replace(key, value)
    if "__" in template:
        raise AssertionError("unexpanded configuration token remains")
    return template


def read_frame(model: Path, frame: int) -> dict[str, np.ndarray]:
    path = Path(f"{model}.save.{frame:06d}.vtkhdf")
    if not path.is_file():
        raise AssertionError(f"missing HDF5 output: {path}")
    with h5py.File(path, "r") as output:
        required = ("coordinate", "connectivity", "pore pressure", "stress", "velocity")
        missing = [name for name in required if name not in output]
        if missing:
            raise AssertionError(f"{path}: missing datasets {missing}")
        return {name: output[name][...] for name in required}


def execute(
    exe: Path,
    config: str,
    case_dir: Path,
    *,
    threads: int,
    log_name: str,
) -> None:
    cfg_path = case_dir / "input.cfg"
    cfg_path.write_text(config, encoding="ascii")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    completed = subprocess.run(
        [str(exe), str(cfg_path)],
        cwd=case_dir,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=90,
        check=False,
    )
    (case_dir / log_name).write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        print(completed.stdout, file=sys.stderr)
        raise RuntimeError(
            f"{case_dir.name}: exit {completed.returncode}; "
            f"see {case_dir / log_name}"
        )


def run_case(
    exe: Path,
    template: str,
    run_root: Path,
    ndims: int,
    threads: int,
    case: Case,
) -> RunResult:
    case_dir = run_root / f"{ndims}d" / f"omp{threads}" / case.name
    case_dir.mkdir(parents=True, exist_ok=False)
    modelname = f"hm_density_{ndims}d_omp{threads}_{case.name}"
    execute(
        exe,
        render_config(template, case, modelname),
        case_dir,
        threads=threads,
        log_name="fresh.log",
    )
    model = case_dir / modelname
    return RunResult(
        initial=read_frame(model, 0),
        final=read_frame(model, 1),
        case_dir=case_dir,
        modelname=modelname,
    )


def check_hydrostatic_water(result: RunResult) -> None:
    coord = result.initial["coordinate"]
    connectivity = result.initial["connectivity"].astype(int)
    expected_pressure = 500.0 * GRAVITY * (0.25 - coord[:, -1])
    for frame in (result.initial, result.final):
        np.testing.assert_allclose(
            frame["pore pressure"], expected_pressure, rtol=1e-13, atol=1e-10
        )

    saturated_density = SOLID_DENSITY * 0.75 + 500.0 * 0.25
    center_z = coord[connectivity, -1].mean(axis=1)
    expected_stress_pressure = (
        saturated_density * GRAVITY * (-center_z)
        + 1200.0 * GRAVITY * 0.25
    )
    ndims = coord.shape[1]
    for component in range(ndims):
        np.testing.assert_allclose(
            result.initial["stress"][:, component],
            -expected_stress_pressure,
            rtol=1e-13,
            atol=1e-10,
        )

    max_velocity = float(np.max(np.abs(result.final["velocity"])))
    if max_velocity > 1e-12:
        raise AssertionError(
            "water traction, reference pressure, and saturated body force are not "
            f"in equilibrium; max velocity={max_velocity:.17e}"
        )


def check_signed_pore_datum(result: RunResult) -> None:
    coord = result.initial["coordinate"]
    expected = 500.0 * GRAVITY * (-0.25 - coord[:, -1])
    if not float(expected.min()) < 0.0 or not float(expected.max()) > 0.0:
        raise AssertionError("signed-datum fixture does not cross the pore-pressure datum")
    for frame in (result.initial, result.final):
        np.testing.assert_allclose(
            frame["pore pressure"], expected, rtol=1e-13, atol=1e-10
        )


def check_mass_scaling(low: RunResult, high: RunResult) -> None:
    low_density = SOLID_DENSITY * 0.75 + 500.0 * 0.25
    high_density = SOLID_DENSITY * 0.75 + 1500.0 * 0.25
    low_scaled = low.final["velocity"] * low_density
    high_scaled = high.final["velocity"] * high_density
    if not float(np.max(np.abs(low.final["velocity"]))) > 0.0:
        raise AssertionError("traction fixture produced no dynamic response")
    np.testing.assert_allclose(low_scaled, high_scaled, rtol=1e-11, atol=1e-14)


def check_feature_off(first: RunResult, second: RunResult) -> None:
    for field in ("velocity", "stress", "coordinate"):
        np.testing.assert_allclose(
            first.final[field], second.final[field], rtol=0.0, atol=1e-14
        )


def check_prem_water(dry: RunResult, wet: RunResult) -> None:
    ndims = dry.initial["coordinate"].shape[1]
    sea_pressure = 1200.0 * GRAVITY * 0.25
    expected_delta = np.zeros_like(dry.initial["stress"])
    expected_delta[:, :ndims] = -sea_pressure
    np.testing.assert_allclose(
        wet.initial["stress"] - dry.initial["stress"],
        expected_delta,
        rtol=0.0,
        atol=1e-10,
    )
    np.testing.assert_allclose(
        wet.initial["pore pressure"], dry.initial["pore pressure"],
        rtol=0.0, atol=1e-12,
    )
    np.testing.assert_allclose(
        wet.final["velocity"], dry.final["velocity"],
        rtol=1e-12, atol=1e-14,
    )


def check_restart(
    exe: Path,
    template: str,
    run_root: Path,
    ndims: int,
    threads: int,
    case: Case,
    baseline: RunResult,
) -> None:
    restart_dir = run_root / f"{ndims}d" / f"omp{threads}" / f"{case.name}_restart"
    shutil.copytree(baseline.case_dir, restart_dir)
    execute(
        exe,
        render_config(
            template, case, baseline.modelname, restarting=True, restart_frame=0
        ),
        restart_dir,
        threads=threads,
        log_name="restart.log",
    )
    restarted = read_frame(restart_dir / baseline.modelname, 1)
    for field in ("coordinate", "velocity", "stress", "pore pressure"):
        np.testing.assert_allclose(
            restarted[field], baseline.final[field], rtol=1e-13, atol=1e-14
        )


def run_dimension(
    exe: Path,
    template: str,
    run_root: Path,
    ndims: int,
    threads: list[int],
) -> None:
    reference_by_case: dict[str, dict[str, np.ndarray]] = {}
    case_by_name = {case.name: case for case in CASES}
    for thread_count in threads:
        results = {
            case.name: run_case(exe, template, run_root, ndims, thread_count, case)
            for case in CASES
        }

        check_hydrostatic_water(results["hydrostatic_water"])
        check_signed_pore_datum(results["signed_pore_datum"])
        check_mass_scaling(results["mass_fluid_500"], results["mass_fluid_1500"])
        check_feature_off(results["feature_off_a"], results["feature_off_b"])
        check_prem_water(results["prem_without_water"], results["prem_with_water"])
        check_prem_water(
            results["prem_modified_without_water"],
            results["prem_modified_with_water"],
        )
        check_restart(
            exe,
            template,
            run_root,
            ndims,
            thread_count,
            case_by_name["hydrostatic_water"],
            results["hydrostatic_water"],
        )

        for name, result in results.items():
            if name in reference_by_case:
                for field in ("velocity", "stress", "pore pressure"):
                    np.testing.assert_allclose(
                        result.final[field], reference_by_case[name][field],
                        rtol=1e-13, atol=1e-14,
                    )
            else:
                reference_by_case[name] = result.final
        print(f"hydromechanics density {ndims}D OMP={thread_count}: PASS")


def run_all(args: argparse.Namespace, run_root: Path) -> None:
    template = TEMPLATE.read_text(encoding="ascii")
    executables = [(2, args.exe_2d.expanduser().resolve())]
    if args.exe_3d is not None:
        executables.append((3, args.exe_3d.expanduser().resolve()))
    for ndims, exe in executables:
        if not exe.is_file():
            raise FileNotFoundError(exe)
        run_dimension(exe, template, run_root, ndims, args.threads)


def main() -> None:
    args = parse_args()
    if any(value <= 0 for value in args.threads):
        raise ValueError("--threads values must be positive")
    if len(set(args.threads)) != len(args.threads):
        raise ValueError("--threads values must be unique")

    if args.run_dir is not None:
        run_root = args.run_dir.expanduser().resolve()
        run_root.mkdir(parents=True, exist_ok=False)
        run_all(args, run_root)
    else:
        with tempfile.TemporaryDirectory(prefix="des-hm-density-") as tmp:
            run_all(args, Path(tmp))


if __name__ == "__main__":
    main()
