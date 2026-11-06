#!/usr/bin/env python3
"""Compare RCB/METIS on identical pinned CPUs, checking every physical snapshot.

Use the same METIS-enabled MPI executable for both methods. Warmups precede five
paired, randomized rounds by default. Run without other simulations or builds.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import re
import statistics
import sys
import tempfile

from run_mpi import ROOT, compare_physical_hdf5, configure, execute, metrics


CASES = {
    'small-fixed-2d': (2, False, 500, 1000, 1000),
    'large-fixed-2d': (2, False, 125, 1000, 1000),
    'adaptive-2d': (2, True, 1000, 196, 49),
    'small-fixed-3d': (3, False, 1000, 1000, 1000),
    'large-fixed-3d': (3, False, 250, 1000, 1000),
    'adaptive-3d': (3, True, 1000, 80, 19),
}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frame_sets(directory):
    return {kind: sorted(int(p.name.split('.')[-2])
                         for p in directory.glob(f'result.{kind}.*.vtkhdf'))
            for kind in ('save', 'chkpt')}


def partition_metrics(log):
    calls = [dict(re.findall(r'(\w+)=(\S+)', line))
             for line in re.findall(r'\[mpi\]\[partitioner\] (.*)', log)]
    generations = {}
    for line in re.findall(r'\[mpi\]\[partition\] (.*)', log):
        row = {k: int(v) for k, v in re.findall(r'(\w+)=(\d+)', line)}
        generations.setdefault(row['generation'], []).append(row)
    if not calls or len(calls) != len(generations):
        raise RuntimeError('Missing partitioner timings or ownership generations')
    return dict(seconds=sum(float(c['seconds']) for c in calls), calls=calls,
                generations=[sorted(rows, key=lambda r: r['rank'])
                             for _, rows in sorted(generations.items())])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--serial', type=Path, required=True)
    parser.add_argument('--mpi', type=Path, required=True)
    parser.add_argument('--serial3d', type=Path)
    parser.add_argument('--mpi3d', type=Path)
    parser.add_argument('--ranks', nargs='+', type=int, default=[3, 4])
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--cases', nargs='+', choices=CASES, default=list(CASES))
    parser.add_argument('--workdir', type=Path, required=True)
    args = parser.parse_args()
    cpus = sorted(os.sched_getaffinity(0))
    if (args.repeats < 3 or min(args.ranks) < 2 or max(args.ranks) > len(cpus)
            or len(set(args.ranks)) != len(args.ranks)):
        parser.error('Need >=3 repeats and distinct available rank counts >=2')
    executables = {2: (args.serial.resolve(), args.mpi.resolve())}
    if args.serial3d and args.mpi3d:
        executables[3] = (args.serial3d.resolve(), args.mpi3d.resolve())
    if any(CASES[name][0] not in executables for name in args.cases):
        parser.error('3D cases require --serial3d and --mpi3d')
    args.workdir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='partitioners-', dir=args.workdir))
    report = dict(completed=False, cpus=cpus, repeats=args.repeats, cases=[],
                  warmups=[], measurements=[], sources={}, executables={})
    for dim, pair in executables.items():
        for label, path in zip(('parent', 'mpi'), pair):
            report['executables'][f'{label}{dim}d'] = dict(path=str(path), sha256=sha256(path))
    for path in [ROOT/'mpi_partition.cxx', ROOT/'mpi_simulation.cxx', Path(__file__)]:
        report['sources'][str(path.relative_to(ROOT))] = sha256(path)
    print('Artifacts:', work, flush=True)

    def save_report():
        (work/'results.json').write_text(json.dumps(report, indent=2)+'\n')

    env = dict(os.environ, OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE',
               OMP_PROC_BIND='false', OMP_PLACES='cores', DES_MPI_OVERLAP='1')
    env.pop('DES_MPI_PARTITION_TRACE', None)
    try:
        for name in args.cases:
            dim, moving, resolution, steps, interval = CASES[name]
            cfg = configure(dim, moving, resolution, steps=steps, output_interval=interval)
            reference = work/(name+'-parent')
            reference.mkdir()
            with (reference/'model.cfg').open('w') as stream:
                cfg.write(stream)
            print(name, 'parent reference', flush=True)
            log, elapsed = execute([str(executables[dim][0]), 'model.cfg'], reference, env)
            if 'Ending simulation.' not in log:
                raise RuntimeError('Parent did not finish')
            frames = frame_sets(reference)
            if len(frames['save']) < 2:
                raise RuntimeError('Need initial and final physical snapshots')
            initial = re.search(r'initial nodes=(\d+), elements=(\d+)', log)
            baseline = metrics(log, elapsed)
            if moving and baseline['remeshes'] == 0:
                raise RuntimeError('Adaptive case did not remesh')
            report['cases'].append(dict(name=name, ndims=dim, resolution=resolution,
                moving=moving, steps=steps, frames=frames, directory=str(reference),
                initial_nodes=int(initial[1]), initial_elements=int(initial[2]), parent=baseline))

            def run(method, ranks, tag, repeat=None):
                directory = work/f'{name}-{method}{ranks}-{tag}'
                directory.mkdir()
                with (directory/'model.cfg').open('w') as stream:
                    cfg.write(stream)
                command = ['mpiexec', '--bind-to', 'none', '-n', str(ranks), sys.executable,
                    str(Path(__file__).with_name('pin_cpu.py').resolve()), '--mpi',
                    '--cpus', ','.join(map(str, cpus[:ranks])), '--record', str(directory),
                    '--exe', str(executables[dim][1])]
                run_env = dict(env, DES_MPI_PARTITIONER=method)
                log, elapsed = execute(command, directory, run_env)
                if 'Ending simulation.' not in log:
                    raise RuntimeError('MPI simulation did not finish')
                if re.search(r'omp_threads=(\d+)', log)[1] != '1':
                    raise RuntimeError('MPI run used unexpected OpenMP thread count')
                bindings = sorted((json.loads(p.read_text()) for p in directory.glob('binding-*.json')),
                                  key=lambda row: row['rank'])
                if (len(bindings) != ranks or [r['rank'] for r in bindings] != list(range(ranks))
                        or any(row['cpus'] != [cpus[row['rank']]] for row in bindings)):
                    raise RuntimeError('Incorrect rank affinity')
                part = partition_metrics(log)
                if any(c['requested'] != method or c['actual'] != method for c in part['calls']):
                    raise RuntimeError('Requested partitioner was not used')
                if frame_sets(directory) != frames:
                    raise RuntimeError('Snapshot/checkpoint frame sets differ')
                for kind, numbers in frames.items():
                    for frame in numbers:
                        compare_physical_hdf5(reference, directory, frame, kind=kind)
                row = metrics(log, elapsed)
                if row['remeshes'] != baseline['remeshes']:
                    raise RuntimeError('Remeshing counts differ from parent')
                if len(part['calls']) != 1+row['remeshes']:
                    raise RuntimeError('Partitions were not rebuilt after every remesh')
                row.update(case=name, ndims=dim, method=method, ranks=ranks,
                           repeat=repeat, partition=part, command=command, bindings=bindings,
                           bit_exact=True, directory=str(directory))
                return row

            rng = random.Random(20261106+dim+resolution)
            for ranks in args.ranks:
                methods = ['rcb', 'metis']
                rng.shuffle(methods)
                for method in methods:
                    report['warmups'].append(run(method, ranks, 'warmup'))
            for repeat in range(args.repeats):
                order = list(args.ranks)
                rng.shuffle(order)
                for ranks in order:
                    methods = ['rcb', 'metis']
                    rng.shuffle(methods)
                    for method in methods:
                        row = run(method, ranks, f'sample{repeat}', repeat)
                        report['measurements'].append(row)
                        save_report()
                print(name, 'round', repeat+1, 'complete', flush=True)
            for ranks in args.ranks:
                for method in ('rcb', 'metis'):
                    rows = [r for r in report['measurements']
                            if r['case'] == name and r['ranks'] == ranks and r['method'] == method]
                    print(name, ranks, method, 'median seconds:',
                          statistics.median(r['simulation_s'] for r in rows), flush=True)
        report['completed'] = True
    finally:
        save_report()
    print('Results:', work/'results.json', flush=True)


if __name__ == '__main__':
    main()
