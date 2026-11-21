#!/usr/bin/env python3
"""Reproduce repository examples and compare parent serial / CPU MPI results.

Full mode preserves physical, stopping, mesh and output settings. Regression mode
caps steps and increases output frequency; every override is recorded. Maxwell
uses the separate, documented current-format constitutive benchmark.
"""
import argparse
import configparser
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import signal
import subprocess
import sys
import tempfile
import time

import h5py
import numpy as np

from run_mpi import ROOT, compare_physical_hdf5, command_text
from run_gpu_mpi import compare_fields, device_records

CASES = {
    "core-mmg": ("examples/core-complex-mmg.cfg", 2, 20000),
    "rsf-eq": ("examples/rate_and_state_friction/ep_rsf_eq.cfg", 2, 2000),
    "rsf-creep": ("examples/rate_and_state_friction/ep_rsf_creep.cfg", 2, 2000),
    "rsf-large-dc": ("examples/rate_and_state_friction/ep_rsf_creep_largeDc.cfg", 2, 2000),
    "evp-rifting": ("examples/rifting-2d.cfg", 2, 1000),
    "maxwell-2d": ("benchmarks/maxwell-current.cfg", 2, 5000),
    "maxwell-3d": ("benchmarks/maxwell-current.cfg", 3, 5000),
}


def digest(path):
    with Path(path).open("rb") as stream:
        hasher = hashlib.sha256()
        for chunk in iter(lambda: stream.read(1024*1024), b""):
            hasher.update(chunk)
        return hasher.hexdigest()


def configuration(name, mode, max_steps=None):
    source, ndims, cap = CASES[name]
    path = ROOT / source
    cfg = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    cfg.optionxform = str
    cfg.read(path,encoding='utf-8')
    original = {section: dict(cfg[section]) for section in cfg.sections()}
    for key in ("poly_filename", "exo_filename"):
        if key in cfg["mesh"]:
            asset = (path.parent / cfg["mesh"][key]).resolve()
            if not asset.is_file():
                raise FileNotFoundError(asset)
            cfg["mesh"][key] = str(asset)
    if mode == "regression" or max_steps is not None:
        steps = max_steps if max_steps is not None else cap
        if "max_steps" in cfg["sim"]:
            steps = min(steps, int(cfg["sim"]["max_steps"]))
        cfg["sim"]["max_steps"] = str(steps)
        cfg["sim"]["output_step_interval"] = str(max(1, steps // 4))
        cfg["sim"]["is_outputting_averaged_fields"] = "no"
    changes = {f"{section}.{key}": {"old": original.get(section, {}).get(key), "new": val}
               for section in cfg.sections() for key, val in cfg[section].items()
               if original.get(section, {}).get(key) != val}
    return cfg, dict(source=source, source_sha256=digest(path), ndims=ndims,
                     mode=mode, overrides=changes,
                     stop_conditions_preserved=not any(key.startswith("sim.max_") for key in changes),
                     constitutive_adaptation=name.startswith("maxwell"))


def execute(command, directory, threads, timeout):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE")
    started = time.monotonic()
    with (directory / "run.log").open("w") as log:
        process = subprocess.Popen(command, cwd=directory, env=env, stdout=log,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            status = process.wait(timeout=timeout)
        except BaseException:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                # The group can exit before cleanup reaches killpg; still reap
                # the child with wait below.
                pass
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            raise
    if status:
        raise RuntimeError(f"Exit {status}: {command_text(command)}; {directory / 'run.log'}")
    log = (directory / "run.log").read_text()
    if "Ending simulation." not in log:
        raise RuntimeError(f"No completion marker: {directory}")
    return dict(wall_with_launcher_s=time.monotonic() - started,
                command=command, threads=threads, executable_sha256=digest(command[-2]))


def frames(directory, prefix, kind="save"):
    p = directory / prefix
    return sorted(int(f.name.split(".")[-2]) for f in p.parent.glob(p.name + f".{kind}.*.vtkhdf"))


def compare_monitors(reference, current):
    files = sorted(reference.rglob("*point_*.csv"))
    if [p.relative_to(reference) for p in files] != sorted(
            p.relative_to(current) for p in current.rglob("*point_*.csv")):
        raise AssertionError("Monitor file sets differ")
    reports = []
    for path in files:
        relative = path.relative_to(reference)
        # This CSV contains physical quantities and identities only, with no wall clock.
        with path.open() as left, (current / relative).open() as right:
            a, b = list(csv.DictReader(left)), list(csv.DictReader(right))
        if not a or a != b:
            raise AssertionError(f"Monitor records differ: {relative}")
        for row in a:
            if not all(np.isfinite(float(value)) for value in row.values()):
                raise AssertionError(f"Non-finite monitor record: {relative}")
        reports.append(dict(file=str(relative), records=len(a), bit_exact=True))
    return reports


def physical_summary(directory, prefix, ndims):
    rows = []
    for frame in frames(directory, prefix):
        with h5py.File(directory / f"{prefix}.save.{frame:06d}.vtkhdf") as f:
            grid = f["VTKHDF/grid"]
            cells = grid["CellData"]
            stress = cells["stress"][()]
            rows.append(dict(frame=frame, step=int(grid["FieldData/steps"][()].item()),
                             time_yr=float(grid["FieldData/time_yr"][()].item()),
                             nelem=len(cells["stress"]),
                             mean_stress_zz=float(stress[:, ndims-1].mean()),
                             min_stress_zz=float(stress[:, ndims-1].min()),
                             max_stress_zz=float(stress[:, ndims-1].max()),
                             max_plastic_strain=float(np.max(cells["plastic strain"][()])),
                             max_plastic_strain_rate=float(np.max(np.abs(cells["plastic strain-rate"][()]))),
                             max_speed=float(np.max(np.linalg.norm(grid["PointData/velocity"][()], axis=1)))))
    if not rows:
        raise AssertionError("No saved frames")
    log = (directory / "run.log").read_text()
    steps = re.search(r"Compute\s*:.*?\)/\s*(\d+)\s*=", log)
    return dict(saved=rows, completed="Ending simulation." in log,
                completed_steps=int(steps[1]) if steps else None,
                remeshes=log.count("Remeshing starts"),
                sampled_dynamic_frames=sum(row["max_speed"] > 1e-6 and
                                           row["max_plastic_strain_rate"] > 0 for row in rows))


def check_maxwell(summary):
    rows = summary["saved"][1:]
    # The saved time and .info are formatted; step * fixed_dt is exact here.
    times = np.array([r["step"] for r in rows], dtype=float) * 1e9
    expected = 2e8 * (-np.expm1(-times / 1e12))
    actual = np.array([r["mean_stress_zz"] for r in rows])
    # Check every element via the extrema as well as the spatial mean.
    extrema = np.array([[r["min_stress_zz"], r["max_stress_zz"]] for r in rows])
    error = np.abs((extrema - expected[:, None]) / expected[:, None])
    if not len(error) or not np.isfinite(error).all() or error.max() >= 1e-6:
        raise AssertionError(f"Maxwell analytical error too large: {error}")
    return dict(max_relative_error=float(error.max()), tolerance=1e-6,
                relaxation_times=float(times[-1] / 1e12),
                final_stress_pa=float(actual[-1]), analytical_final_stress_pa=float(expected[-1]))


def analyze(reference, current, cfg, ndims, gpu=False):
    prefix = cfg["sim"]["modelname"]
    checked = {}
    comparisons = []
    for kind in ("save", "chkpt"):
        old, new = frames(reference, prefix, kind), frames(current, prefix, kind)
        if old != new or (kind == "save" and len(old) < 2):
            raise AssertionError(f"Saved frame sets differ: {kind}: {old} / {new}")
        for frame in old:
            if gpu:
                comparisons.append(compare_fields(reference,current,frame,kind,prefix=prefix))
            else:
                compare_physical_hdf5(reference, current, frame, prefix, kind)
        checked[kind] = old
    monitors = compare_monitors(reference, current)
    summary = physical_summary(current, prefix, ndims)
    baseline = physical_summary(reference, prefix, ndims)
    if not baseline["completed"] or not summary["completed"]:
        raise AssertionError("Cannot label an incomplete run reproduced")
    if (baseline["completed_steps"], baseline["remeshes"]) != (
            summary["completed_steps"], summary["remeshes"]):
        raise AssertionError("Completed steps/remesh counts differ")
    # Save the normal repository comparator output as well as full-dataset checks.
    env = dict(os.environ, PYTHONPATH=str(ROOT))
    with (current / "compare-standard.log").open("w") as output:
        subprocess.run([sys.executable, str(ROOT / "benchmarks-cores/compare.py"),
                        str(reference / prefix), str(current / prefix), str(checked["save"][-1])],
                       env=env, stdout=output, stderr=subprocess.STDOUT, check=True)
    return dict(bit_exact=all(row['bit_exact'] for row in comparisons), comparisons=comparisons,
                frames=checked, monitors=monitors, physical=summary,
                reference=str(reference), current=str(current))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--serial", type=Path, required=True)
    parser.add_argument("--mpi", type=Path, required=True)
    parser.add_argument("--serial3d", type=Path)
    parser.add_argument("--mpi3d", type=Path)
    parser.add_argument("--cases", nargs="+", choices=CASES,
                        default=["core-mmg", "rsf-eq", "maxwell-2d"])
    parser.add_argument("--mode", choices=("full", "regression"), default="full")
    parser.add_argument("--max-steps", type=int)
    parser.add_argument("--ranks", type=int, default=2)
    parser.add_argument("--serial-threads", type=int, default=1)
    parser.add_argument("--launcher", default="mpiexec")
    parser.add_argument("--timeout", type=float, default=7200, help="Seconds allowed per executable")
    parser.add_argument("--workdir", type=Path)
    parser.add_argument("--gpu", action="store_true",help="Verify GPU MPI identities and compare physical fields within 1e-8")
    parser.add_argument("--gpu-reference", action="store_true",help="Compare GPU MPI against a one-GPU run; retain CPU comparisons separately")
    args = parser.parse_args()
    if args.gpu_reference and not args.gpu:parser.error('--gpu-reference requires --gpu')
    if args.ranks < 2 or args.serial_threads < 1 or args.timeout <= 0 or (
            args.max_steps is not None and args.max_steps < 4):
        parser.error("ranks >= 2, threads >= 1, timeout > 0 and max-steps >= 4 are required")
    executables = {2: (args.serial.resolve(), args.mpi.resolve())}
    if args.serial3d and args.mpi3d:
        executables[3] = (args.serial3d.resolve(), args.mpi3d.resolve())
    for case in args.cases:
        if CASES[case][1] not in executables:
            parser.error(f"{case} requires both 3D executables")
    for pair in executables.values():
        for path in pair:
            if not path.is_file():
                parser.error(f"Executable missing: {path}")
    if args.workdir:
        args.workdir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="reproduce-", dir=args.workdir))
    print(f"Artifacts: {work}", flush=True)
    reference_kind='gpu1' if args.gpu_reference else 'parent'
    results = dict(cases={}, mode=args.mode, gpu_tested=args.gpu, completed=False,
                   reference_backend='GPU MPI, one rank' if args.gpu_reference else 'parent CPU')
    try:
        for name in args.cases:
            cfg, provenance = configuration(name, args.mode, args.max_steps)
            directory = work / name
            directory.mkdir()
            (directory / "provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
            (directory / "original.cfg").write_text((ROOT / provenance["source"]).read_text(encoding='utf-8'),encoding='utf-8')
            timings = {}
            for kind, exe in zip((reference_kind, "mpi"), executables[provenance["ndims"]]):
                if args.gpu_reference:exe=executables[provenance['ndims']][1]
                target = directory / kind
                target.mkdir()
                for prefix in (cfg["sim"]["modelname"], cfg.get("monitor", "output_prefix", fallback="monitor")):
                    (target / Path(prefix).parent).mkdir(parents=True, exist_ok=True)
                with (target / "model.cfg").open("w") as stream:
                    cfg.write(stream)
                command = [str(exe), "model.cfg"]
                threads = args.serial_threads if kind == "parent" else 1
                if kind == "mpi" or args.gpu_reference:
                    command = shlex.split(args.launcher) + ["-n", str(args.ranks if kind=='mpi' else 1)] + command
                print(f"Running {name}: {kind}", flush=True)
                timings[kind] = execute(command, target, threads, args.timeout)
            report = analyze(directory / reference_kind, directory / "mpi", cfg, provenance["ndims"], args.gpu)
            if args.gpu:
                report['devices']=device_records((directory/'mpi/run.log').read_text(),args.ranks)
                if args.gpu_reference:
                    report['reference_devices']=device_records((directory/'gpu1/run.log').read_text(),1)
            report.update(provenance=provenance, timings=timings)
            if name.startswith("maxwell"):
                report["analytical"] = check_maxwell(report["physical"])
            results["cases"][name] = report
            (directory / "comparison.json").write_text(json.dumps(report, indent=2)+"\n")
            (work / "results.json").write_text(json.dumps(results, indent=2)+"\n")
            print(f"Passed {name}: {len(report['frames']['save'])} frames, bit_exact={report['bit_exact']}", flush=True)
        results["completed"] = True
    except Exception as error:
        results['failure']=str(error)
        raise
    finally:
        (work / "results.json").write_text(json.dumps(results, indent=2)+"\n")
    print(f"Results: {work / 'results.json'}", flush=True)


if __name__ == "__main__":
    main()
