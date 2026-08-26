"""Shared configuration, execution, and diagnostic helpers for hydraulic CFL tests."""

from __future__ import annotations

import math
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


DT_LINE = re.compile(r"compute_dt:\s+(?P<fields>.*?)\s+sec$")
DT_FIELD = re.compile(r"(?P<name>[a-z_]+)=(?P<value>[^ ]+)")


@dataclass(frozen=True)
class CflCase:
    name: str
    hydraulic_enabled: bool
    permeability: tuple[float, float]
    expected_diffusivity_max: float | None
    quasi_static: bool = False
    layer_mattypes: tuple[int, int] = (0, 1)
    expect_selected_hydraulic: bool = True
    moving_mesh: bool = False
    minimum_dt_records: int = 1


def parse_numeric_dt_fields(text: str) -> dict[str, float]:
    fields: dict[str, float] = {}
    for field in DT_FIELD.finditer(text):
        try:
            fields[field.group("name")] = float(field.group("value"))
        except ValueError:
            # Debug diagnostics also include tuples such as max_vel_coord.
            continue
    return fields


def render_cfg(template: str, case: CflCase) -> str:
    replacements = {
        "__MODELNAME__": f"hydraulic_cfl_{case.name}",
        "__HYDRAULIC_ENABLED__": "yes" if case.hydraulic_enabled else "no",
        "__QUASI_STATIC__": "yes" if case.quasi_static else "no",
        "__MOVING_MESH__": "yes" if case.moving_mesh else "no",
        "__LAYER_MATTYPES__": ", ".join(
            str(value) for value in case.layer_mattypes
        ),
        "__PERM0__": f"{case.permeability[0]:.17e}",
        "__PERM1__": f"{case.permeability[1]:.17e}",
    }
    rendered = template
    for key, value in replacements.items():
        rendered = rendered.replace(key, value)
    return rendered


def parse_dt_records(stdout: str, case_name: str) -> list[dict[str, float]]:
    records: list[dict[str, float]] = []
    for line in stdout.splitlines():
        match = DT_LINE.search(line)
        if not match:
            continue
        fields = parse_numeric_dt_fields(match.group("fields"))
        records.append(fields)
    if not records:
        raise AssertionError(f"{case_name}: no compute_dt diagnostics found")
    required = {"selected", "hydro_diffusion", "min_length"}
    for index, record in enumerate(records):
        missing = required - record.keys()
        if missing:
            raise AssertionError(
                f"{case_name}: compute_dt record {index} lacks {sorted(missing)}"
            )
    return records


def run_case(
    exe: Path, template: str, run_root: Path, case: CflCase
) -> list[dict[str, float]]:
    case_dir = run_root / case.name
    case_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = case_dir / "input.cfg"
    cfg_path.write_text(render_cfg(template, case), encoding="ascii")
    stdout = run_config(exe, cfg_path, case_dir, case.name)
    return parse_dt_records(stdout, case.name)


def run_config(exe: Path, cfg_path: Path, run_dir: Path, name: str) -> str:
    result = subprocess.run(
        [str(exe), str(cfg_path)],
        cwd=run_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    (run_dir / f"{name}.log").write_text(result.stdout, encoding="utf-8")
    if result.returncode != 0:
        print(result.stdout, file=sys.stderr)
        raise RuntimeError(f"{name}: DynEarthSol exited with {result.returncode}")
    return result.stdout


def assert_close(actual: float, expected: float, message: str) -> None:
    if not math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=0.0):
        raise AssertionError(
            f"{message}: actual={actual:.17e}, expected={expected:.17e}"
        )


def replace_once(text: str, old: str, new: str) -> str:
    if text.count(old) != 1:
        raise AssertionError(f"expected exactly one config token: {old!r}")
    return text.replace(old, new, 1)
