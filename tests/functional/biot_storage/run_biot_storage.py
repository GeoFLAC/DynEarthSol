#!/usr/bin/env python3

"""Check Biot storage in 2-D/3-D point injection and input validation."""

from __future__ import annotations

import argparse
import os
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
DEFAULT_CFG = HERE / "biot_storage_base.cfg"
BULK_MODULUS = 5.0e9
SHEAR_MODULUS = 3.0e9
POROSITY = 0.2
BIOT_COEFFICIENT = 0.8


@dataclass(frozen=True)
class MaterialCase:
    name: str
    fluid_bulk_modulus: float
    grain_bulk_modulus: float


@dataclass(frozen=True)
class InvalidCase:
    name: str
    replacements: dict[str, str]
    diagnostic: str
    hydraulic: bool = True
    effective_stress: bool = False


MATERIALS = (
    MaterialCase("storage_a", 2.0e9, 2.0e10),
    MaterialCase("storage_b", 1.0e9, 1.0e10),
)

INVALID_CASES = (
    InvalidCase(
        "zero_fluid_bulk_modulus",
        {"__FLUID_BULK_MODULUS__": "0"},
        "mat.fluid_bulk_modulus must be finite and positive",
    ),
    InvalidCase(
        "zero_grain_bulk_modulus",
        {"__GRAIN_BULK_MODULUS__": "0"},
        "mat.bulk_modulus_s must be finite and positive",
    ),
    InvalidCase(
        "zero_fluid_viscosity",
        {"__FLUID_VISCOSITY__": "0"},
        "mat.fluid_visc must be finite and positive",
    ),
    InvalidCase(
        "zero_storage",
        {"__POROSITY__": "0", "__BIOT_COEFFICIENT__": "0"},
        "hydraulic pressure storage must be finite and positive",
    ),
    InvalidCase(
        "zero_diffusivity",
        {"__PERMEABILITY__": "0"},
        "hydraulic diffusivity must be finite and positive",
    ),
    InvalidCase(
        "overflow_diffusivity",
        {
            "__PERMEABILITY__": "1e308",
            "__FLUID_VISCOSITY__": "1e-308",
        },
        "hydraulic diffusivity must be finite and positive",
    ),
    InvalidCase(
        "effective_stress_zero_fluid_bulk_modulus",
        {"__FLUID_BULK_MODULUS__": "0", "__GRAVITY__": "10"},
        "mat.fluid_bulk_modulus must be finite and positive",
        hydraulic=False,
        effective_stress=True,
    ),
    InvalidCase(
        "effective_stress_zero_grain_bulk_modulus",
        {"__GRAIN_BULK_MODULUS__": "0", "__GRAVITY__": "10"},
        "mat.bulk_modulus_s must be finite and positive",
        hydraulic=False,
        effective_stress=True,
    ),
    InvalidCase(
        "effective_stress_zero_storage",
        {
            "__POROSITY__": "0",
            "__BIOT_COEFFICIENT__": "0",
            "__GRAVITY__": "10",
        },
        "initial Skempton pressure storage must be finite and positive",
        hydraulic=False,
        effective_stress=True,
    ),
)

ACTIVE_FIELDS = ("coordinate", "pore pressure")
FEATURE_OFF_FIELDS = (
    "coordinate",
    "displacement",
    "velocity",
    "force",
    "pore pressure",
    "stress",
    "strain",
    "density",
    "temperature",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe-2d", type=Path, default=REPO_ROOT / "dynearthsol2d")
    parser.add_argument("--exe-3d", type=Path, required=True)
    parser.add_argument("--cfg", type=Path, default=DEFAULT_CFG)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 4])
    parser.add_argument("--run-dir", type=Path)
    return parser.parse_args()


def source_section(ndims: int, enabled: bool) -> str:
    if not enabled:
        return "[injection]\nenabled = no\nnum_points = 0"
    y_coordinate = "points_y = [0.5]\n" if ndims == 3 else ""
    return (
        "[injection]\n"
        "enabled = yes\n"
        "num_points = 1\n"
        "points_x = [0.5]\n"
        f"{y_coordinate}"
        "points_z = [-0.5]\n"
        "points_unit = m\n"
        "rate_model = constant_rate\n"
        "rate = [1e-6]\n"
        "total_amount = [0]\n"
        "start_time_in_yr = [0]\n"
        "end_time_in_yr = []\n"
        "duration_in_yr = []"
    )


def render_config(
    template: str,
    model: str,
    ndims: int,
    material: MaterialCase,
    hydraulic: bool,
    feedback: bool,
    effective_stress: bool = False,
    extra_replacements: dict[str, str] | None = None,
) -> str:
    replacements = {
        "__MODELNAME__": model,
        "__GRAVITY__": "0",
        "__HYDRAULIC__": "yes" if hydraulic else "no",
        "__FEEDBACK__": "yes" if feedback else "no",
        "__EFFECTIVE_STRESS__": "yes" if effective_stress else "no",
        "__PLANE_STRAIN__": "yes" if ndims == 2 else "no",
        "__POROSITY__": f"{POROSITY:.17e}",
        "__PERMEABILITY__": "1e-30",
        "__FLUID_BULK_MODULUS__": f"{material.fluid_bulk_modulus:.17e}",
        "__FLUID_VISCOSITY__": "1e-3",
        "__BIOT_COEFFICIENT__": f"{BIOT_COEFFICIENT:.17e}",
        "__GRAIN_BULK_MODULUS__": f"{material.grain_bulk_modulus:.17e}",
        "__SOURCE_SECTION__": source_section(ndims, hydraulic),
    }
    if extra_replacements:
        replacements.update(extra_replacements)
    rendered = template
    for token, value in replacements.items():
        rendered = rendered.replace(token, value)
    if "__" in rendered:
        raise AssertionError(f"{model}: unexpanded configuration token remains")
    return rendered


def read_frame(case_dir: Path, model: str, frame: int, fields: tuple[str, ...]) -> dict[str, np.ndarray]:
    path = case_dir / f"{model}.save.{frame:06d}.vtkhdf"
    if not path.is_file():
        raise AssertionError(f"missing output frame: {path}")
    values: dict[str, np.ndarray] = {}
    with h5py.File(path, "r") as output:
        for field in fields:
            if field not in output:
                raise AssertionError(f"{path}: missing dataset {field!r}")
            values[field] = output[field][...]
    return values


def run_model(
    exe: Path,
    config: str,
    run_root: Path,
    model: str,
    threads: int,
    fields: tuple[str, ...],
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    case_dir = run_root / model
    case_dir.mkdir(parents=True, exist_ok=False)
    cfg_path = case_dir / "input.cfg"
    cfg_path.write_text(config, encoding="ascii")
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
        timeout=90,
        check=False,
    )
    (case_dir / "run.log").write_text(result.stdout, encoding="utf-8")
    if result.returncode != 0:
        print(result.stdout, file=sys.stderr)
        raise RuntimeError(f"{model}: exit {result.returncode}")
    return (
        read_frame(case_dir, model, 0, fields),
        read_frame(case_dir, model, 1, fields),
    )


def storage(material: MaterialCase, ndims: int, feedback: bool) -> float:
    value = (BIOT_COEFFICIENT - POROSITY) / material.grain_bulk_modulus
    value += POROSITY / material.fluid_bulk_modulus
    if feedback:
        constrained_modulus = BULK_MODULUS
        if ndims == 2:
            constrained_modulus += SHEAR_MODULUS / 3.0
        value += BIOT_COEFFICIENT**2 / constrained_modulus
    return value


def assert_arrays_close(
    actual: dict[str, np.ndarray],
    expected: dict[str, np.ndarray],
    context: str,
) -> None:
    if actual.keys() != expected.keys():
        raise AssertionError(f"{context}: output fields differ")
    for field in actual:
        np.testing.assert_allclose(
            actual[field],
            expected[field],
            rtol=2.0e-13,
            atol=2.0e-12,
            err_msg=f"{context}: {field}",
        )


def check_pressure_ratio(
    first: tuple[dict[str, np.ndarray], dict[str, np.ndarray]],
    second: tuple[dict[str, np.ndarray], dict[str, np.ndarray]],
    expected_ratio: float,
    context: str,
) -> None:
    np.testing.assert_array_equal(first[0]["coordinate"], second[0]["coordinate"])
    delta_first = first[1]["pore pressure"] - first[0]["pore pressure"]
    delta_second = second[1]["pore pressure"] - second[0]["pore pressure"]
    amplitude = max(float(np.max(np.abs(delta_first))), float(np.max(np.abs(delta_second))))
    if not np.isfinite(amplitude) or amplitude <= 0.0:
        raise AssertionError(f"{context}: injection did not produce finite pressure")
    support = np.abs(delta_first) + np.abs(delta_second) > amplitude * 1.0e-13
    if not np.any(support):
        raise AssertionError(f"{context}: point-source support is empty")
    np.testing.assert_allclose(
        delta_first[support],
        expected_ratio * delta_second[support],
        rtol=2.0e-12,
        atol=2.0e-12 * amplitude,
        err_msg=f"{context}: pressure increment/storage ratio",
    )


def check_active_cases(
    exe: Path,
    template: str,
    run_root: Path,
    ndims: int,
    threads: list[int],
) -> None:
    for feedback in (False, True):
        thread_results: dict[int, dict[str, tuple[dict[str, np.ndarray], dict[str, np.ndarray]]]] = {}
        for thread_count in threads:
            material_results = {}
            for material in MATERIALS:
                model = (
                    f"biot_{material.name}_{ndims}d_"
                    f"feedback_{int(feedback)}_omp{thread_count}"
                )
                config = render_config(
                    template, model, ndims, material, True, feedback
                )
                material_results[material.name] = run_model(
                    exe, config, run_root, model, thread_count, ACTIVE_FIELDS
                )
            expected_ratio = storage(MATERIALS[1], ndims, feedback) / storage(
                MATERIALS[0], ndims, feedback
            )
            check_pressure_ratio(
                material_results[MATERIALS[0].name],
                material_results[MATERIALS[1].name],
                expected_ratio,
                f"{ndims}D feedback={feedback} OMP={thread_count}",
            )
            thread_results[thread_count] = material_results

        reference = thread_results[threads[0]]
        for thread_count in threads[1:]:
            for material in MATERIALS:
                actual = thread_results[thread_count][material.name]
                expected = reference[material.name]
                assert_arrays_close(
                    actual[0], expected[0], f"{ndims}D {material.name} initial OMP"
                )
                assert_arrays_close(
                    actual[1], expected[1], f"{ndims}D {material.name} final OMP"
                )
        print(
            f"Biot storage point injection {ndims}D feedback={feedback}: PASS",
            flush=True,
        )


def check_feature_off(
    exe: Path,
    template: str,
    run_root: Path,
    ndims: int,
    threads: list[int],
) -> None:
    thread_results = {}
    for thread_count in threads:
        material_results = {}
        for material in MATERIALS:
            model = f"biot_feature_off_{material.name}_{ndims}d_omp{thread_count}"
            config = render_config(
                template, model, ndims, material, False, False
            )
            material_results[material.name] = run_model(
                exe, config, run_root, model, thread_count, FEATURE_OFF_FIELDS
            )
        for frame in (0, 1):
            assert_arrays_close(
                material_results[MATERIALS[0].name][frame],
                material_results[MATERIALS[1].name][frame],
                f"{ndims}D feature-off frame {frame} OMP={thread_count}",
            )
        thread_results[thread_count] = material_results

    reference = thread_results[threads[0]]
    for thread_count in threads[1:]:
        for material in MATERIALS:
            for frame in (0, 1):
                assert_arrays_close(
                    thread_results[thread_count][material.name][frame],
                    reference[material.name][frame],
                    f"{ndims}D feature-off {material.name} frame {frame} OMP",
                )
    print(f"Biot storage feature-off {ndims}D: PASS", flush=True)


def check_invalid_inputs(exe: Path, template: str, run_root: Path) -> None:
    for case in INVALID_CASES:
        model = f"biot_invalid_{case.name}"
        case_dir = run_root / model
        case_dir.mkdir(parents=True, exist_ok=False)
        config = render_config(
            template,
            model,
            2,
            MATERIALS[0],
            case.hydraulic,
            False,
            case.effective_stress,
            case.replacements,
        )
        cfg_path = case_dir / "input.cfg"
        cfg_path.write_text(config, encoding="ascii")
        result = subprocess.run(
            [str(exe), str(cfg_path)],
            cwd=case_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=30,
            check=False,
        )
        (case_dir / "run.log").write_text(result.stdout, encoding="utf-8")
        if result.returncode != 11:
            raise AssertionError(
                f"{case.name}: exit {result.returncode}, expected 11"
            )
        if case.diagnostic not in result.stdout:
            raise AssertionError(
                f"{case.name}: missing diagnostic {case.diagnostic!r}"
            )
        print(f"Biot storage invalid input {case.name}: PASS", flush=True)


def run_all(args: argparse.Namespace, template: str, run_root: Path) -> None:
    threads = list(dict.fromkeys(args.threads))
    if len(threads) != len(args.threads) or any(value <= 0 for value in threads):
        raise ValueError("--threads values must be unique positive integers")
    executables = (
        (2, args.exe_2d.expanduser().resolve()),
        (3, args.exe_3d.expanduser().resolve()),
    )
    for ndims, exe in executables:
        if not exe.is_file():
            raise FileNotFoundError(exe)
        check_active_cases(exe, template, run_root, ndims, threads)
        check_feature_off(exe, template, run_root, ndims, threads)
    check_invalid_inputs(executables[0][1], template, run_root)


def main() -> None:
    args = parse_args()
    template = args.cfg.expanduser().resolve().read_text(encoding="ascii")
    if args.run_dir is not None:
        run_root = args.run_dir.expanduser().resolve()
        run_root.mkdir(parents=True, exist_ok=False)
        run_all(args, template, run_root)
    else:
        with tempfile.TemporaryDirectory(prefix="des-biot-storage-") as tmp:
            run_all(args, template, Path(tmp))


if __name__ == "__main__":
    main()
