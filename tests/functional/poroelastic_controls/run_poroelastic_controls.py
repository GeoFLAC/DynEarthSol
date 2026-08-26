#!/usr/bin/env python3

"""Check independent poroelastic controls and Biot-modulus validation."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]


@dataclass(frozen=True)
class Case:
    name: str
    feedback: str = "no"
    derive_biot: str = "no"
    bulk_modulus: str = "2e8"
    bulk_modulus_s: str = "37e9"
    expected_code: int = 0
    expected_text: str = ""


VALID_CASES = (
    Case("feedback_off"),
    Case("feedback_on_without_hydraulics", feedback="yes"),
    Case(
        "derive_biot_valid",
        feedback="yes",
        derive_biot="yes",
        bulk_modulus_s="4e8",
    ),
    Case(
        "derive_biot_equal_moduli",
        derive_biot="yes",
        bulk_modulus="4e8",
        bulk_modulus_s="4e8",
    ),
)


INVALID_CASES = (
    Case(
        "derive_biot_zero_grain_modulus",
        derive_biot="yes",
        bulk_modulus_s="0",
        expected_code=11,
        expected_text="mat.bulk_modulus_s must be finite and positive",
    ),
    Case(
        "derive_biot_negative_grain_modulus",
        derive_biot="yes",
        bulk_modulus_s="-4e8",
        expected_code=11,
        expected_text="mat.bulk_modulus_s must be finite and positive",
    ),
    Case(
        "derive_biot_nonfinite_grain_modulus",
        derive_biot="yes",
        bulk_modulus_s="nan",
        expected_code=11,
        expected_text="incorrect format for mat.bulk_modulus_s",
    ),
    Case(
        "derive_biot_zero_drained_modulus",
        derive_biot="yes",
        bulk_modulus="0",
        expected_code=11,
        expected_text="mat.bulk_modulus must be finite and positive",
    ),
    Case(
        "derive_biot_negative_drained_modulus",
        derive_biot="yes",
        bulk_modulus="-2e8",
        expected_code=11,
        expected_text="mat.bulk_modulus must be finite and positive",
    ),
    Case(
        "derive_biot_nonfinite_drained_modulus",
        derive_biot="yes",
        bulk_modulus="nan",
        expected_code=11,
        expected_text="incorrect format for mat.bulk_modulus",
    ),
    Case(
        "derive_biot_drained_exceeds_grain",
        derive_biot="yes",
        bulk_modulus="4e8",
        bulk_modulus_s="2e8",
        expected_code=11,
        expected_text="mat.bulk_modulus must not exceed mat.bulk_modulus_s",
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe-2d", type=Path, default=REPO_ROOT / "dynearthsol2d")
    parser.add_argument("--exe-3d", type=Path, default=REPO_ROOT / "dynearthsol3d")
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 4])
    parser.add_argument("--run-dir", type=Path, default=HERE / "runs")
    parser.add_argument("--clean", action="store_true")
    return parser.parse_args()


def render_config(template: str, case: Case, model: str, ndims: int) -> str:
    replacements = {
        "__MODELNAME__": model,
        "__PRESSURE_FEEDBACK__": case.feedback,
        "__DERIVE_BIOT__": case.derive_biot,
        "__BULK_MODULUS__": case.bulk_modulus,
        "__BULK_MODULUS_S__": case.bulk_modulus_s,
        "__PLANE_STRAIN__": "yes" if ndims == 2 else "no",
    }
    rendered = template
    for key, value in replacements.items():
        rendered = rendered.replace(key, value)
    if "__" in rendered:
        raise AssertionError("unexpanded configuration token remains")
    return rendered


def run_case(
    exe: Path,
    template: str,
    run_root: Path,
    case: Case,
    ndims: int,
    threads: int,
) -> None:
    case_dir = run_root / f"{ndims}d" / f"omp{threads}" / case.name
    case_dir.mkdir(parents=True, exist_ok=False)
    model = f"poro_controls_{ndims}d_omp{threads}_{case.name}"
    cfg_path = case_dir / "input.cfg"
    cfg_path.write_text(render_config(template, case, model, ndims), encoding="ascii")

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
    print(
        f"PASS {case.name} ({ndims}D OMP={threads}, exit {completed.returncode})",
        flush=True,
    )


def main() -> None:
    args = parse_args()
    run_root = args.run_dir.expanduser().resolve()
    if args.clean and run_root.exists():
        shutil.rmtree(run_root)
    run_root.mkdir(parents=True, exist_ok=True)

    threads = list(dict.fromkeys(args.threads))
    if len(threads) != len(args.threads) or any(value <= 0 for value in threads):
        raise ValueError("--threads values must be unique positive integers")

    executables = (
        (2, args.exe_2d.expanduser().resolve()),
        (3, args.exe_3d.expanduser().resolve()),
    )
    for _, exe in executables:
        if not exe.is_file():
            raise FileNotFoundError(exe)

    template = (HERE / "poroelastic_controls_base.cfg").read_text(encoding="ascii")

    # Invalid values exercise dimension-independent input validation once.
    for case in INVALID_CASES:
        run_case(executables[0][1], template, run_root, case, 2, threads[0])

    # Valid controls must parse and complete in every supported build/thread mode.
    for ndims, exe in executables:
        for thread_count in threads:
            for case in VALID_CASES:
                run_case(exe, template, run_root, case, ndims, thread_count)


if __name__ == "__main__":
    main()
