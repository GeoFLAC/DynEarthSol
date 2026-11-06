#!/usr/bin/env python3
"""Compare CPU MPI timestep kernels with a serial executable and measure its overhead.

Pass an executable built on the parent commit as --serial for regression testing.
All generated configs, logs, snapshots and results.json stay under --workdir (a new
temporary directory by default). No existing run or reference is overwritten.
Both executables must use HDF5 output; full physical datasets are checked exactly.
"""
import argparse
import configparser
import json
import os
from pathlib import Path
import re
import shlex
import signal
import statistics
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parents[3]


def configure(ndims, moving, resolution, steps=196, output_interval=49):
    # Checkpoints fall between mesh checks: outputs precede same-step remeshing in DES.
    cfg = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    cfg.optionxform = str
    cfg.read(ROOT / "examples/core-complex-mmg.cfg")
    cfg["sim"].update(modelname="result", max_steps=str(steps), max_time_in_yr="1e20",
                      output_step_interval=str(output_interval),
                      output_time_interval_in_yr="1e20", checkpoint_frame_interval="1")
    cfg["mesh"].update(resolution=str(resolution), quality_check_step_interval="10",
                       max_boundary_distortion="1e-6" if moving else "0.2", mmg_verbose="-1")
    cfg["control"]["has_moving_mesh"] = "yes" if moving else "no"
    if ndims != 2:
        cfg["mesh"].update(xlength="10000", ylength="4000", zlength="4000")
    return cfg


def execute(command, directory, env, logname="run.log", expected_failure=False):
    start = time.perf_counter()
    with (directory / logname).open("w") as log:
        process = subprocess.Popen(command, cwd=directory, env=env, stdout=log,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            returncode = process.wait(timeout=180)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            raise RuntimeError(f"Timed out: {shlex.join(command)}; log: {directory / logname}")
    elapsed = time.perf_counter() - start
    log = (directory / logname).read_text()
    if bool(returncode) != expected_failure:
        raise RuntimeError(f"Command returned {returncode}: {shlex.join(command)}\n"
                           f"{log[-5000:]}")
    return log, elapsed


def compare_physical_hdf5(reference, current, frame, prefix="result", kind="save"):
    import h5py
    import numpy as np

    name = f"{prefix}.{kind}.{frame:06d}.vtkhdf"
    with h5py.File(reference / name) as old, h5py.File(current / name) as new:
        def datasets(file):
            found = {}
            def visit(path, obj):
                if isinstance(obj, h5py.Dataset) and path.split("/")[-1] != "walltime_sec":
                    found[path] = obj
            file["VTKHDF"].visititems(visit)
            return found
        left, right = datasets(old), datasets(new)
        if not left or left.keys() != right.keys():
            raise RuntimeError(f"Physical HDF5 dataset sets differ in frame {frame}")
        for path in left:
            a, b = left[path][()], right[path][()]
            if a.shape != b.shape or a.dtype != b.dtype:
                raise RuntimeError(f"Shape or dtype differs: frame {frame}, {path}")
            if a.dtype.kind in "fc" and (not np.isfinite(a).all() or not np.isfinite(b).all()):
                raise RuntimeError(f"Non-finite physical data: frame {frame}, {path}")
            if a.tobytes() != b.tobytes():
                raise RuntimeError(f"Physical data differ: frame {frame}, {path}")
        report = dict(frame=frame, bit_exact=True, datasets=sorted(left))
        (current / f"full-compare-{kind}-{frame}.json").write_text(json.dumps(report, indent=2) + "\n")


def compare(reference, current, frame, env):
    # The standard comparator checks invariants and selected coordinates; also check
    # complete vectors/tensors, connectivity, material fields, force and marker records.
    compare_physical_hdf5(reference, current, frame)
    log, _ = execute([sys.executable, str(ROOT / "benchmarks-cores/compare.py"),
                      str(reference / "result"), str(current / "result"), str(frame)],
                     current, env, f"compare-{frame}.log")
    return "Status: BIT-EXACT" in log


def metrics(log, elapsed):
    data = {"wall_with_launcher_s": elapsed}
    # The simulation timer excludes MPI initialization and launcher startup.
    match = re.search(r"Execute\s*:\s*(\d+):(\d+):([\d.]+)", log)
    if match:
        data["simulation_s"] = 3600*int(match[1]) + 60*int(match[2]) + float(match[3])
    data["kernels"] = {name: dict(calls=int(calls), seconds=float(seconds))
                       for name, calls, seconds in re.findall(
                           r"\[mpi\]\[kernel\] name=(\w+) calls=(\d+) seconds=([\d.eE+-]+)", log)}
    data["ownership_reports"] = log.count("[mpi][partition]")
    match = re.search(r"\[mpi\]\[communication\] (.*)", log)
    if match:
        data["communication"] = {k: float(v) for k,v in re.findall(r"(\w+)=([\d.eE+-]+)",match[1])}
    data["remeshes"] = log.count("Remeshing starts")
    data["rank_profile"] = []
    for line in re.findall(r"\[mpi\]\[rank-profile\] (.*)",log):
        row=dict(re.findall(r"(\w+)=(\S+)",line))
        data["rank_profile"].append({k:(v if k=="operation" else float(v)) for k,v in row.items()})
    data["overlap"] = [{k:float(v) for k,v in re.findall(r"(\w+)=([\d.eE+-]+)",line)}
                       for line in re.findall(r"\[mpi\]\[overlap\] (.*)",log)]
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--serial", type=Path, required=True)
    parser.add_argument("--mpi", type=Path, required=True)
    parser.add_argument("--serial3d", type=Path)
    parser.add_argument("--mpi3d", type=Path)
    parser.add_argument("--launcher", default="mpiexec", help="Launcher and optional arguments")
    parser.add_argument("--workdir", type=Path)
    parser.add_argument("--benchmark", action="store_true", help="Also time two fixed-mesh sizes")
    parser.add_argument("--previous-mpi", type=Path, help="Also time the preceding MPI implementation")
    parser.add_argument("--benchmark-only", action="store_true", help="Run checked timings only")
    parser.add_argument("--ranks", type=int, nargs="+", default=[1,2,3,4,6,8])
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if bool(args.serial3d) != bool(args.mpi3d):
        parser.error("--serial3d and --mpi3d must be provided together")
    if not args.ranks or min(args.ranks) < 1 or len(set(args.ranks)) != len(args.ranks):
        parser.error("--ranks must contain distinct positive process counts")
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    executables = {"serial": args.serial.resolve(), "mpi": args.mpi.resolve()}
    if args.previous_mpi:
        executables["previous-mpi"] = args.previous_mpi.resolve()
    if args.serial3d:
        executables.update(serial3d=args.serial3d.resolve(), mpi3d=args.mpi3d.resolve())
    for path in executables.values():
        if not path.is_file():
            parser.error(f"Missing executable: {path}")
    work = args.workdir.resolve() if args.workdir else Path(tempfile.mkdtemp(prefix="des-mpi-"))
    work.mkdir(parents=True, exist_ok=True)
    # Reserve a fresh child so even an explicitly reused --workdir cannot replace a run.
    work = Path(tempfile.mkdtemp(prefix="run-", dir=work))
    print(f"Artifacts: {work}", flush=True)
    launcher = shlex.split(args.launcher)
    env = dict(os.environ, OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE")
    env["PYTHONPATH"] = os.pathsep.join(filter(None, (str(ROOT), env.get("PYTHONPATH"))))
    results = {"comparison": "All VTKHDF datasets exact and finite; walltime_sec excluded",
               "executables": {k: str(v) for k, v in executables.items()}, "checks": [],
               "benchmarks": []}

    def run(name, cfg, exe, ranks=None, threads=1):
        directory = work / name
        directory.mkdir()
        with (directory / "model.cfg").open("w") as f:
            cfg.write(f)
        command = [str(exe), "model.cfg"]
        if ranks is not None:
            command = launcher + ["-n", str(ranks)] + command
        run_env = dict(env, OMP_NUM_THREADS=str(threads))
        log, elapsed = execute(command, directory, run_env)
        return directory, metrics(log, elapsed)

    if not args.benchmark_only:
        # Keep example geometry, materials and physics. Only cap run duration and
        # choose deterministic, frequent output; resolve mesh assets against examples/.
        cases = [("defaults.cfg",2,40),("aniso_test.cfg",2,2),
                 ("core-complex.cfg",2,80),("core-complex-mmg.cfg",2,80),
                 ("subduction.cfg",2,40)]
        if args.serial3d:
            cases += [("defaults.cfg",3,40),("aniso_test.cfg",3,2)]
        for example, ndims, steps in cases:
            cfg=configparser.ConfigParser(inline_comment_prefixes=("#",";"))
            cfg.optionxform = str
            cfg.read(ROOT / "examples" / example)
            cfg["sim"].update(modelname="result",max_steps=str(steps),max_time_in_yr="1e20",
                              output_step_interval=str(max(1,steps//4)),
                              output_time_interval_in_yr="1e20",is_outputting_averaged_fields="no")
            if "poly_filename" in cfg["mesh"]:
                cfg["mesh"]["poly_filename"]=str(ROOT / "examples" / cfg["mesh"]["poly_filename"])
            name=f"example-{Path(example).stem}-{ndims}d"
            serial=executables["serial" if ndims==2 else "serial3d"]
            mpi=executables["mpi" if ndims==2 else "mpi3d"]
            reference, baseline=run(name+"-parent",cfg,serial)
            for ranks, threads in [(r,1) for r in args.ranks if r>1]+[(2,2)]:
                current,timing=run(f"{name}-mpi{ranks}-omp{threads}",cfg,mpi,ranks,threads)
                frames=sorted(int(p.name.split(".")[-2]) for p in reference.glob("result.save.*.vtkhdf"))
                exact=all([compare(reference,current,frame,env) for frame in frames])
                required={"Temperature","StrainRate","Stress","Force","Velocity","Mass"}
                if not required.issubset(timing["kernels"]):
                    raise RuntimeError(f"Missing MPI kernels: {required-set(timing['kernels'])}")
                results["checks"].append(dict(case=name,ranks=ranks,threads=threads,frames=frames,
                                             bit_exact=exact,timing=timing,serial_timing=baseline))
                print(f"{name}: ranks={ranks} OMP={threads}: bit-exact ({len(frames)} frames)",flush=True)

        # Exercise optional kernels that the stock examples leave disabled.
        optional = [(name,ndims) for name in ("hydraulic","pt-hydraulic","thermal-off")
                    for ndims in ((2,3) if args.serial3d else (2,))]
        optional.append(("hydration",2))
        for name,ndims in optional:
            if name=="hydration":
                cfg=configparser.ConfigParser(inline_comment_prefixes=("#",";"))
                cfg.optionxform = str
                cfg.read(ROOT / "examples/subduction.cfg")
                cfg["sim"].update(modelname="result",max_steps="40",max_time_in_yr="1e20",
                                  output_step_interval="10",output_time_interval_in_yr="1e20",
                                  is_outputting_averaged_fields="no")
                cfg["mesh"]["poly_filename"]=str(ROOT / "examples/poly-subduction.poly")
                cfg["control"]["has_hydration_processes"]="yes"
                cfg["mat"]["phase_change_option"]="1"
            else:
                cfg=configure(ndims,False,1000,steps=8,output_interval=2)
                if name=="thermal-off":
                    cfg["control"]["has_thermal_diffusion"]="no"
                else:
                    cfg["control"]["has_hydraulic_diffusion"]="yes"
                    cfg["mat"]["porosity"]="[0.01]"
                    cfg["mat"]["hydraulic_perm"]="[1e-20]"
                if name=="pt-hydraulic":
                    cfg["control"].update(has_PT="yes",PT_max_iter="4")
            label=f"optional-{name}-{ndims}d"
            reference,_=run(label+"-parent",cfg,executables["serial" if ndims==2 else "serial3d"])
            current,timing=run(label+"-mpi2",cfg,executables["mpi" if ndims==2 else "mpi3d"],2)
            exact=all([compare(reference,current,frame,env) for frame in range(1,5)])
            results["checks"].append(dict(case=label,ranks=2,bit_exact=exact,timing=timing))
            print(f"{label}: bit-exact",flush=True)

        for moving in (False, True):
            case = "remesh" if moving else "fixed"
            cfg = configure(2, moving, 1000)
            reference, baseline = run(f"{case}-parent", cfg, executables["serial"])
            if moving and not baseline["remeshes"]:
                raise RuntimeError("The remesh test did not exercise remeshing")
            for ranks, threads in [(r,1) for r in args.ranks]+[(2,2)]:
                current, timing = run(f"{case}-mpi{ranks}-omp{threads}", cfg, executables["mpi"],
                                      ranks, threads)
                exact = all([compare(reference, current, frame, env) for frame in range(1, 5)])
                if moving and ranks > 1 and timing["ownership_reports"] < 2*ranks:
                    raise RuntimeError("Ownership was not updated across remeshing")
                results["checks"].append(dict(case=case, ranks=ranks, threads=threads,
                                               bit_exact=exact, timing=timing))
                print(f"{case}: ranks={ranks} OMP={threads}: "
                      f"{'bit-exact' if exact else 'within regression tolerance'}", flush=True)
            cfg["sim"].update(is_restarting="yes", restarting_from_modelname=str(reference / "result"),
                              restarting_from_frame="2")
            current, timing = run(f"{case}-restart", cfg, executables["mpi"], 2)
            exact = compare(reference, current, 4, env)
            results["checks"].append(dict(case=case + "-restart", ranks=2, bit_exact=exact, timing=timing))
            print(f"{case}: restart on two ranks: "
                  f"{'bit-exact' if exact else 'within regression tolerance'}", flush=True)

        # Averaging maintains persistent local accumulators across remeshing.
        # Compare all output frames in both dimensions, including the initial one.
        for ndims in ((2, 3) if args.serial3d else (2,)):
            cfg = configure(ndims, True, 1000, steps=40, output_interval=10)
            cfg["sim"]["is_outputting_averaged_fields"] = "yes"
            suffix = "" if ndims == 2 else "3d"
            label = f"averaged-remesh-{ndims}d"
            reference, _ = run(label + "-parent", cfg, executables["serial" + suffix])
            for ranks in (r for r in args.ranks if r > 1):
                current, timing = run(f"{label}-mpi{ranks}", cfg, executables["mpi" + suffix], ranks)
                exact = all([compare(reference, current, frame, env) for frame in range(5)])
                results["checks"].append(dict(case=label, ranks=ranks, bit_exact=exact,
                                               frames=list(range(5)), timing=timing))
                print(f"{label}: ranks={ranks}: bit-exact", flush=True)

        cfg = configure(2, False, 1000)
        cfg["ic"]["isostasy_adjustment_time_in_yr"] = "10"
        reference, _ = run("isostasy-parent", cfg, executables["serial"])
        current, timing = run("isostasy-mpi2", cfg, executables["mpi"], 2)
        exact = compare(reference, current, 4, env)
        results["checks"].append(dict(case="fixed-isostasy", ranks=2, bit_exact=exact, timing=timing))
        print(f"Fixed mesh after isostasy: {'bit-exact' if exact else 'within tolerance'}", flush=True)

        cfg = configure(2, False, 10000, steps=4, output_interval=1)
        cfg["mesh"].update(meshing_option="1", xlength="1000", zlength="1000")
        cfg["ic"]["weakzone_option"] = "0"
        reference, _ = run("empty-partitions-parent", cfg, executables["serial"])
        current, timing = run("empty-partitions-mpi6", cfg, executables["mpi"], 6)
        exact = compare(reference, current, 4, env)
        results["checks"].append(dict(case="empty-partitions", ranks=6, bit_exact=exact, timing=timing))
        print(f"Empty node partitions: {'bit-exact' if exact else 'within tolerance'}", flush=True)

        if args.serial3d:
            cfg = configure(3, False, 1000)
            reference, _ = run("3d-parent", cfg, executables["serial3d"])
            current, timing = run("3d-mpi2", cfg, executables["mpi3d"], 2)
            exact = all([compare(reference, current, frame, env) for frame in range(1, 5)])
            results["checks"].append(dict(case="3d-fixed", ranks=2, bit_exact=exact, timing=timing))
            print(f"3d: {'bit-exact' if exact else 'within regression tolerance'}", flush=True)

        # Fail after MPI initialization while the worker waits for work; failure must terminate it.
        cfg = configure(2, False, 1000)
        cfg["sim"]["modelname"] = "missing-directory/result"
        directory = work / "root-io-failure"
        directory.mkdir()
        with (directory / "model.cfg").open("w") as f:
            cfg.write(f)
        log, _ = execute(launcher + ["-n", "2", str(executables["mpi"]), "model.cfg"],
                         directory, env, expected_failure=True)
        if "H5Fcreate failed" not in log or "[DES exit 61]" not in log:
            raise RuntimeError("The failure test did not reach the expected HDF5 exception handler")
        print("Root I/O failure terminates the MPI job: passed", flush=True)

        cfg["sim"]["modelname"] = "result"
        cfg["mesh"].update(meshing_option="91", poly_filename="missing.poly")
        with (directory / "model.cfg").open("w") as f:
            cfg.write(f)
        log, _ = execute(launcher + ["-n", "2", str(executables["mpi"]), "model.cfg"],
                         directory, env, "missing-mesh.log", expected_failure=True)
        if "[DES exit 20]" not in log:
            raise RuntimeError("The missing-mesh test did not reach DES exit 20")
        execute(launcher + ["-n", "2", str(executables["mpi"]), "--help"],
                directory, env, "help.log")
        print("DES error and --help exit without hanging: passed", flush=True)

    if args.benchmark or args.benchmark_only:
        for resolution in (500, 125):
            cfg = configure(2, False, resolution, steps=1000, output_interval=1000)
            reference = None
            variants = [("serial", None, 1), ("serial-omp2", None, 2)]
            if args.previous_mpi:
                variants.append(("previous-mpi2", 2, 1))
            variants += [(f"mpi{r}", r, 1) for r in args.ranks]
            for variant, ranks, threads in variants:
                timings = []
                exe = executables["serial"] if ranks is None else executables["mpi"]
                if variant == "previous-mpi2":
                    exe = executables["previous-mpi"]
                for repeat in range(args.repeats):
                    current, timing = run(f"bench-{resolution}-{variant}-{repeat}", cfg, exe,
                                          ranks, threads)
                    if reference is None:
                        reference = current
                    elif repeat == 0:
                        compare(reference, current, 1, env)
                    timings.append(timing)
                median = statistics.median(t["simulation_s"] for t in timings)
                results["benchmarks"].append(dict(resolution=resolution, variant=variant,
                                                   ranks=ranks, threads=threads,
                                                   median_simulation_s=median, samples=timings))
                print(f"resolution={resolution} {variant}: "
                      f"median simulation {median:.6f} s", flush=True)
    (work / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(f"Passed. Results: {work / 'results.json'}")


if __name__ == "__main__":
    main()
