#!/usr/bin/env python3
"""Run representative RSF example CFGs sequentially."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List

sys.dont_write_bytecode = True


def parse_modelname(cfg_path: Path) -> str:
    for raw in cfg_path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            continue
        key, val = line.split("=", 1)
        if key.strip() == "modelname":
            return val.strip()
    return ""


def ensure_output_dir(cfg_path: Path) -> None:
    modelname = parse_modelname(cfg_path)
    if not modelname or "/" not in modelname:
        return
    out_dir = modelname.rsplit("/", 1)[0]
    target = cfg_path.parent / out_dir
    target.mkdir(parents=True, exist_ok=True)
    print(f"[prep] output dir: {target}")


def candidate_executables(script_dir: Path) -> Iterable[Path]:
    repo_root = script_dir.parent.parent
    yield repo_root / "dynearthsol2d"
    yield repo_root / "binaries" / "dynearthsol2d"


def resolve_executable(script_dir: Path, explicit: str | None) -> Path:
    if explicit:
        p = Path(explicit).expanduser().resolve()
        if p.is_file() and os.access(p, os.X_OK):
            return p
        raise FileNotFoundError(f"--exe not executable: {p}")

    env_exe = os.environ.get("DYNEXE")
    if env_exe:
        p = Path(env_exe).expanduser().resolve()
        if p.is_file() and os.access(p, os.X_OK):
            return p
        raise FileNotFoundError(f"DYNEXE is set but not executable: {p}")

    for p in candidate_executables(script_dir):
        if p.is_file() and os.access(p, os.X_OK):
            return p

    raise FileNotFoundError(
        "DynEarthSol 2-D executable not found. Build dynearthsol2d in the repo root or pass --exe /path/to/dynearthsol2d."
    )


def cleanup_pycache(root_dir: Path) -> None:
    removed = 0
    for path in root_dir.rglob("__pycache__"):
        if path.is_dir():
            shutil.rmtree(path, ignore_errors=True)
            removed += 1
    if removed > 0:
        print(f"[cleanup] removed {removed} __pycache__ director{'y' if removed == 1 else 'ies'}")


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Run the three 2-D RSF example CFGs sequentially.")
    parser.add_argument(
        "--exe",
        default=None,
        help="Path to dynearthsol2d. If omitted, use DYNEXE or auto-detect.",
    )
    parser.add_argument(
        "--cfgs",
        nargs="+",
        default=["ep_rsf_creep.cfg", "ep_rsf_eq.cfg", "ep_rsf_creep_largeDc.cfg"],
        help="CFG file list relative to this directory.",
    )
    parser.add_argument(
        "--plot-after",
        action="store_true",
        help="Run plot_rsf_time_series.py after all runs.",
    )
    args = parser.parse_args()

    exe = resolve_executable(script_dir, args.exe)
    print(f"[exe] {exe}")

    cfgs: List[Path] = []
    for cfg in args.cfgs:
        p = (script_dir / cfg).resolve()
        if not p.is_file():
            raise FileNotFoundError(f"CFG not found: {p}")
        cfgs.append(p)

    for idx, cfg in enumerate(cfgs, start=1):
        print("==========================================")
        print(f"[run] {idx}/{len(cfgs)} : {cfg.name}")
        print("==========================================")
        ensure_output_dir(cfg)
        subprocess.run([str(exe), str(cfg)], cwd=script_dir, check=True)

    print("[done] all CFG runs completed")

    if args.plot_after:
        plot_script = script_dir / "plot_rsf_time_series.py"
        if not plot_script.is_file():
            raise FileNotFoundError(f"Plot script not found: {plot_script}")
        subprocess.run(["python3", str(plot_script)], cwd=script_dir, check=True)


if __name__ == "__main__":
    _script_dir = Path(__file__).resolve().parent
    try:
        main()
    finally:
        cleanup_pycache(_script_dir)
