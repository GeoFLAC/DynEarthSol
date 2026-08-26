#!/usr/bin/env python3

"""Check total-stress pressure bridging, PT iteration safety, and restart IO."""

from __future__ import annotations

import argparse
import os
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
DEFAULT_CFG = HERE / "total_stress_base.cfg"
HEADER_BYTES = 4096
PRESSURE_NEW = 40.0
DPPRESSURE = -PRESSURE_NEW  # p_old - p_new, with p_old = 0.
BIOT_COEFFICIENT = 0.25
DPP_NAME = "pore pressure stress increment"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe-2d", type=Path)
    parser.add_argument("--exe-3d", type=Path)
    parser.add_argument("--cfg", type=Path, default=DEFAULT_CFG)
    parser.add_argument(
        "--expected-format", choices=("hdf5", "binary"),
        help="Require the executable(s) to emit this output format.",
    )
    parser.add_argument(
        "--threads", type=int, nargs="+", default=[1, 4],
        help="OMP thread counts to exercise (default: 1 4).",
    )
    parser.add_argument(
        "--run-dir", type=Path,
        help="Keep outputs below this new directory (default: temporary).",
    )
    return parser.parse_args()


def render_config(
    template: str,
    *,
    model: str,
    restart_model: str,
    restarting: bool,
    effective_stress: bool,
    dimension: int,
    has_pt: bool = False,
) -> str:
    replacements = {
        "__MODELNAME__": model,
        "__RESTART_MODEL__": restart_model,
        "__HAS_INITIAL_CHECKPOINT__": "no" if restarting else "yes",
        "__IS_RESTARTING__": "yes" if restarting else "no",
        "__PLANE_STRAIN__": "yes" if dimension == 2 else "no",
        "__HAS_PT__": "yes" if has_pt else "no",
        "__EFFECTIVE_STRESS_OPTION__": (
            "has_pore_pressure_effective_stress = yes"
            if effective_stress else ""
        ),
    }
    rendered = template
    for token, value in replacements.items():
        if token not in rendered:
            raise AssertionError(f"missing config token {token}")
        rendered = rendered.replace(token, value)
    if "__" in rendered:
        raise AssertionError("unexpanded config token remains")
    return rendered


def run_des(
    exe: Path,
    config: str,
    case_dir: Path,
    name: str,
    threads: int,
    *,
    expected_code: int = 0,
    expected_text: str = "",
) -> None:
    cfg_path = case_dir / f"{name}.cfg"
    cfg_path.write_text(config, encoding="ascii")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    completed = subprocess.run(
        [str(exe), str(cfg_path)],
        cwd=case_dir,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=120,
        check=False,
    )
    (case_dir / f"{name}.log").write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != expected_code:
        print(completed.stdout, file=sys.stderr)
        raise AssertionError(
            f"{name} OMP={threads}: DynEarthSol exited with {completed.returncode}, "
            f"expected {expected_code}"
        )
    if expected_text and expected_text not in completed.stdout:
        raise AssertionError(
            f"{name} OMP={threads}: missing diagnostic {expected_text!r}"
        )


def binary_offsets(path: Path) -> dict[str, int]:
    with path.open("rb") as stream:
        header = stream.read(HEADER_BYTES)
    if len(header) != HEADER_BYTES:
        raise AssertionError(f"{path}: truncated binary header")
    text = header.split(b"\0", 1)[0].decode("ascii")
    offsets: dict[str, int] = {}
    for line in text.splitlines()[1:]:
        name, location = line.rsplit("\t", 1)
        offsets[name] = int(location)
    return offsets


def read_binary_scalar(path: Path, name: str) -> int:
    offsets = binary_offsets(path)
    if name not in offsets:
        raise AssertionError(f"{path}: missing binary scalar {name!r}")
    with path.open("rb") as stream:
        stream.seek(offsets[name])
        raw = stream.read(struct.calcsize("=i"))
    return int(struct.unpack("=i", raw)[0])


def read_array(
    path: Path, name: str, count: int, components: int = 1
) -> np.ndarray:
    if path.suffix == ".vtkhdf":
        import h5py

        with h5py.File(path, "r") as output:
            if name not in output:
                raise AssertionError(f"{path}: missing HDF5 dataset {name!r}")
            return np.asarray(output[name][...]).copy()

    offsets = binary_offsets(path)
    if name not in offsets:
        raise AssertionError(f"{path}: missing binary array {name!r}")
    with path.open("rb") as stream:
        stream.seek(offsets[name])
        values = np.fromfile(stream, dtype=np.float64, count=count * components)
    if values.size != count * components:
        raise AssertionError(f"{path}: truncated binary array {name!r}")
    if components > 1:
        values = values.reshape((count, components))
    return values


def write_constant(
    path: Path,
    name: str,
    count: int,
    value: float,
    *,
    hdf_storage: str,
) -> None:
    if path.suffix == ".vtkhdf":
        import h5py

        with h5py.File(path, "r+") as output:
            if hdf_storage not in output:
                raise AssertionError(
                    f"{path}: missing HDF5 storage {hdf_storage!r} for {name!r}"
                )
            dataset = output[hdf_storage]
            if dataset.size != count:
                raise AssertionError(
                    f"{path}: {name!r} size {dataset.size}, expected {count}"
                )
            dataset[...] = value
        return

    offsets = binary_offsets(path)
    if name not in offsets:
        raise AssertionError(f"{path}: missing binary array {name!r}")
    values = np.full(count, value, dtype=np.float64)
    with path.open("r+b") as stream:
        stream.seek(offsets[name])
        stream.write(values.tobytes())


def remove_restart_array(path: Path, name: str) -> None:
    if path.suffix == ".vtkhdf":
        import h5py

        with h5py.File(path, "r+") as output:
            if name not in output:
                raise AssertionError(f"{path}: missing HDF5 dataset {name!r}")
            # Restart reads the root compatibility dataset. The underlying VTK
            # storage may remain: removing this link emulates an older checkpoint.
            del output[name]
        return

    encoded = name.encode("ascii")
    with path.open("r+b") as stream:
        header = stream.read(HEADER_BYTES)
        if encoded not in header:
            raise AssertionError(f"{path}: missing binary header entry {name!r}")
        header = header.replace(encoded, b"x" * len(encoded), 1)
        stream.seek(0)
        stream.write(header)


def output_paths(case_dir: Path, model: str, output_format: str, frame: int):
    suffix = ".vtkhdf" if output_format == "hdf5" else ""
    return (
        case_dir / f"{model}.save.{frame:06d}{suffix}",
        case_dir / f"{model}.chkpt.{frame:06d}{suffix}",
    )


def detect_format(case_dir: Path, model: str) -> str:
    hdf = case_dir / f"{model}.save.000000.vtkhdf"
    binary = case_dir / f"{model}.save.000000"
    if hdf.is_file() == binary.is_file():
        raise AssertionError(
            f"{model}: expected exactly one HDF5/raw frame-0 save file"
        )
    return "hdf5" if hdf.is_file() else "binary"


def assert_close(name: str, actual: np.ndarray, expected: np.ndarray) -> None:
    if actual.shape != expected.shape:
        raise AssertionError(
            f"{name}: shape {actual.shape}, expected {expected.shape}"
        )
    if not np.allclose(actual, expected, rtol=1e-13, atol=1e-12):
        index = np.unravel_index(np.argmax(np.abs(actual - expected)), actual.shape)
        raise AssertionError(
            f"{name}: mismatch at {index}: actual={actual[index]:.17e}, "
            f"expected={expected[index]:.17e}"
        )


def check_restart_result(
    case_dir: Path,
    model: str,
    output_format: str,
    dimension: int,
    nnode: int,
    nelem: int,
    *,
    expected_pressure_stress: float,
    expected_dppressure: float,
) -> tuple[np.ndarray, np.ndarray | None]:
    save, checkpoint = output_paths(case_dir, model, output_format, 1)
    nstr = 3 if dimension == 2 else 6
    stress = read_array(save, "stress", nelem, nstr)
    strain = read_array(save, "strain", nelem, nstr)
    expected_stress = np.zeros((nelem, nstr), dtype=np.float64)
    expected_stress[:, :dimension] = expected_pressure_stress
    assert_close(f"{model} stress", stress, expected_stress)
    assert_close(f"{model} zero strain", strain, np.zeros_like(strain))

    dpp = read_array(checkpoint, DPP_NAME, nnode)
    assert_close(
        f"{model} checkpoint dppressure",
        dpp,
        np.full(nnode, expected_dppressure, dtype=np.float64),
    )

    stressyy = None
    if dimension == 2:
        stressyy = read_array(checkpoint, "stressyy", nelem)
        assert_close(
            f"{model} plane-strain stressyy",
            stressyy,
            np.full(nelem, expected_pressure_stress, dtype=np.float64),
        )
    return stress, stressyy


def run_dimension_thread(
    exe: Path,
    template: str,
    run_root: Path,
    dimension: int,
    threads: int,
    expected_format: str | None,
) -> tuple[str, np.ndarray, np.ndarray | None]:
    case_dir = run_root / f"{dimension}d_omp{threads}"
    case_dir.mkdir(parents=True, exist_ok=False)
    source_model = f"total_stress_source_{dimension}d_omp{threads}"
    source_cfg = render_config(
        template,
        model=source_model,
        restart_model="unused",
        restarting=False,
        effective_stress=True,
        dimension=dimension,
    )
    run_des(exe, source_cfg, case_dir, "source", threads)
    output_format = detect_format(case_dir, source_model)
    if expected_format and output_format != expected_format:
        raise AssertionError(
            f"{exe}: emitted {output_format}, expected {expected_format}"
        )

    source_save, source_checkpoint = output_paths(
        case_dir, source_model, output_format, 0
    )
    if output_format == "hdf5":
        nnode = int(read_array(source_save, "pore pressure", 0).size)
        nelem = int(read_array(source_save, "stress", 0).shape[0])
    else:
        nnode = read_binary_scalar(source_save, "nnode")
        nelem = read_binary_scalar(source_save, "nelem")

    initial_dpp = read_array(source_checkpoint, DPP_NAME, nnode)
    assert_close(
        f"{source_model} initial checkpoint dppressure",
        initial_dpp,
        np.zeros(nnode, dtype=np.float64),
    )
    write_constant(
        source_save,
        "pore pressure",
        nnode,
        PRESSURE_NEW,
        hdf_storage="/VTKHDF/grid/PointData/pore pressure",
    )
    write_constant(
        source_checkpoint,
        DPP_NAME,
        nnode,
        DPPRESSURE,
        hdf_storage=f"/VTKHDF/grid/PointData/{DPP_NAME}",
    )

    enabled_model = f"total_stress_enabled_{dimension}d_omp{threads}"
    enabled_cfg = render_config(
        template,
        model=enabled_model,
        restart_model=source_model,
        restarting=True,
        effective_stress=True,
        dimension=dimension,
    )
    run_des(exe, enabled_cfg, case_dir, "enabled_restart", threads)
    enabled_stress, enabled_stressyy = check_restart_result(
        case_dir,
        enabled_model,
        output_format,
        dimension,
        nnode,
        nelem,
        expected_pressure_stress=-BIOT_COEFFICIENT * PRESSURE_NEW,
        expected_dppressure=DPPRESSURE,
    )

    # PT must retain the current pressure level in every constitutive trial,
    # while applying the physical-step pressure increment exactly once.
    pt_model = f"total_stress_pt_{dimension}d_omp{threads}"
    pt_cfg = render_config(
        template,
        model=pt_model,
        restart_model=source_model,
        restarting=True,
        effective_stress=True,
        dimension=dimension,
        has_pt=True,
    )
    run_des(exe, pt_cfg, case_dir, "pt_restart", threads)
    check_restart_result(
        case_dir,
        pt_model,
        output_format,
        dimension,
        nnode,
        nelem,
        expected_pressure_stress=-BIOT_COEFFICIENT * PRESSURE_NEW,
        expected_dppressure=DPPRESSURE,
    )

    # Omit the new option entirely so this also checks its default-off behavior.
    disabled_model = f"total_stress_disabled_{dimension}d_omp{threads}"
    disabled_cfg = render_config(
        template,
        model=disabled_model,
        restart_model=source_model,
        restarting=True,
        effective_stress=False,
        dimension=dimension,
    )
    run_des(exe, disabled_cfg, case_dir, "disabled_restart", threads)
    check_restart_result(
        case_dir,
        disabled_model,
        output_format,
        dimension,
        nnode,
        nelem,
        expected_pressure_stress=0.0,
        expected_dppressure=DPPRESSURE,
    )

    # Emulate a checkpoint written before dppressure became restart state.
    remove_restart_array(source_checkpoint, DPP_NAME)
    fallback_model = f"total_stress_fallback_{dimension}d_omp{threads}"
    fallback_cfg = render_config(
        template,
        model=fallback_model,
        restart_model=source_model,
        restarting=True,
        effective_stress=True,
        dimension=dimension,
    )
    run_des(exe, fallback_cfg, case_dir, "fallback_restart", threads)
    check_restart_result(
        case_dir,
        fallback_model,
        output_format,
        dimension,
        nnode,
        nelem,
        expected_pressure_stress=0.0,
        expected_dppressure=0.0,
    )

    print(
        f"total-stress {dimension}D {output_format} OMP={threads}: PASS",
        flush=True,
    )
    return output_format, enabled_stress, enabled_stressyy


def run_all(
    executables: dict[int, Path],
    template: str,
    run_root: Path,
    threads: list[int],
    expected_format: str | None,
) -> None:
    if not threads or any(value <= 0 for value in threads):
        raise ValueError("--threads values must be positive")
    if len(set(threads)) != len(threads):
        raise ValueError("--threads values must be unique")

    for dimension, exe in sorted(executables.items()):
        baseline: tuple[np.ndarray, np.ndarray | None] | None = None
        seen_format: str | None = None
        for thread_count in threads:
            output_format, stress, stressyy = run_dimension_thread(
                exe,
                template,
                run_root,
                dimension,
                thread_count,
                expected_format,
            )
            if seen_format is not None and output_format != seen_format:
                raise AssertionError(f"{dimension}D output format changed across OMP runs")
            seen_format = output_format
            if baseline is None:
                baseline = stress, stressyy
            else:
                assert_close(
                    f"{dimension}D OMP stress parity", stress, baseline[0]
                )
                if stressyy is not None and baseline[1] is not None:
                    assert_close(
                        f"{dimension}D OMP stressyy parity", stressyy, baseline[1]
                    )


def main() -> None:
    args = parse_args()
    executables: dict[int, Path] = {}
    if args.exe_2d:
        executables[2] = args.exe_2d.expanduser().resolve()
    if args.exe_3d:
        executables[3] = args.exe_3d.expanduser().resolve()
    if not executables:
        executables = {
            2: (REPO_ROOT / "dynearthsol2d").resolve(),
            3: (REPO_ROOT / "dynearthsol3d").resolve(),
        }
    for dimension, exe in executables.items():
        if not exe.is_file():
            raise FileNotFoundError(f"{dimension}D executable not found: {exe}")

    cfg = args.cfg.expanduser().resolve()
    if not cfg.is_file():
        raise FileNotFoundError(f"config not found: {cfg}")
    template = cfg.read_text(encoding="ascii")

    if args.run_dir:
        run_root = args.run_dir.expanduser().resolve()
        run_root.mkdir(parents=True, exist_ok=True)
        run_all(
            executables,
            template,
            run_root,
            args.threads,
            args.expected_format,
        )
    else:
        with tempfile.TemporaryDirectory(prefix="des-total-stress-") as tmp:
            run_all(
                executables,
                template,
                Path(tmp),
                args.threads,
                args.expected_format,
            )


if __name__ == "__main__":
    main()
