#!/usr/bin/env python3
"""Compare fixed/remeshed/restarted 2D/3D runs on scheduler-allocated nodes.

Requires a shared filesystem, parent/MPI HDF5+MMG executables, NumPy and h5py.
Every fresh frame/checkpoint and every saved continuation frame is compared in
full; wall-clock metadata is excluded. This is a correctness check, not a timing
benchmark. Submit this driver through qsub, with at least two allocated nodes.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

from launch_sge import allocation
from run_mpi import configure, execute, compare_physical_hdf5, metrics


ROOT = Path(__file__).resolve().parents[3]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--serial', type=Path, required=True)
    p.add_argument('--mpi', type=Path, required=True)
    p.add_argument('--serial3d', type=Path, required=True)
    p.add_argument('--mpi3d', type=Path, required=True)
    p.add_argument('--ranks', type=int, nargs='+', default=[2, 4])
    p.add_argument('--ndims', type=int, nargs='+', choices=[2,3], default=[2,3])
    p.add_argument('--workdir', type=Path, required=True)
    a = p.parse_args()
    context = allocation()
    if len(context['host_slots']) < 2:
        p.error('Request at least two scheduler-allocated hosts')
    if len(set(a.ranks)) != len(a.ranks) or min(a.ranks) < 2 or max(a.ranks) > context['nslots']:
        p.error('Use distinct rank counts between two and NSLOTS')
    a.workdir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='matrix-', dir=str(a.workdir))).resolve()
    print('Artifacts:', work, flush=True)
    executables = {name:str(getattr(a, name).resolve()) for name in ('serial', 'mpi', 'serial3d', 'mpi3d')}
    report = dict(completed=False, all_passed=False, allocation=context, checks=[], parent_restart_controls=[],
                  executables=executables,
                  executable_sha256={key:hashlib.sha256(Path(path).read_bytes()).hexdigest()
                                     for key, path in executables.items()})
    target = work / 'results.json'
    def save():
        target.write_text(json.dumps(report, indent=2) + '\n')
    save()
    env = dict(os.environ, OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE', OPENBLAS_NUM_THREADS='1')
    env['PYTHONPATH'] = str(ROOT) + os.pathsep + env.get('PYTHONPATH', '')
    def run(name, cfg, executable, ranks=None, overlap=1):
        directory = work / name
        directory.mkdir()
        with (directory / 'model.cfg').open('w') as stream:
            cfg.write(stream)
        command = [executable, 'model.cfg']
        if ranks:
            command = [sys.executable, str(ROOT/'tests/functional/mpi/launch_sge.py'),
                       '--min-nodes', '2', '-n', str(ranks)] + command
        (directory/'command.json').write_text(json.dumps(command, indent=2)+'\n')
        log, elapsed = execute(command, directory, dict(env, DES_MPI_OVERLAP=str(overlap)))
        return directory, metrics(log, elapsed)
    def compare_all(reference, current, frames):
        for kind in ('save', 'chkpt'):
            for frame in frames:
                compare_physical_hdf5(reference, current, frame, kind=kind)
    try:
        for dim in a.ndims:
            suffix = '' if dim == 2 else '3d'
            for moving in (False, True):
                name = ('remesh' if moving else 'fixed') + '-{}d'.format(dim)
                # Checkpoint 2 is at step 38, between mesh checks at 30 and 40.
                cfg = configure(dim, moving, 1000, steps=76, output_interval=19)
                reference, baseline = run(name+'-parent', cfg, executables['serial'+suffix])
                if moving and not baseline['remeshes']:
                    raise RuntimeError('The parent remesh case did not remesh')
                restart_reference = None
                if moving:
                    cfg['sim'].update(is_restarting='yes', restarting_from_modelname=str(reference/'result'),
                                      restarting_from_frame='2')
                    restart_reference, _ = run(name+'-parent-restart', cfg, executables['serial'+suffix])
                    # Restart reinitializes output scheduling in the parent too.
                    # Require uninterrupted physical frames, and compare the full
                    # checkpoint state against the equivalent parent restart.
                    for frame in (3, 4):
                        compare_physical_hdf5(reference, restart_reference, frame)
                    report['parent_restart_controls'].append(dict(case=name, reference=str(restart_reference),
                                                                  uninterrupted_save_frames_exact=[3, 4]))
                    save()
                    cfg['sim']['is_restarting'] = 'no'
                for ranks in a.ranks:
                    for overlap in (0, 1):
                        label = '{}-mpi{}-overlap{}'.format(name, ranks, overlap)
                        current, timing = run(label, cfg, executables['mpi'+suffix], ranks, overlap)
                        compare_all(reference, current, range(5))
                        if timing['remeshes'] != baseline['remeshes']:
                            raise RuntimeError('Remesh counts differ: '+label)
                        report['checks'].append(dict(case=label, bit_exact=True, timing=timing,
                                                     frames=list(range(5)), checkpoint_frames=list(range(5))))
                        save()
                        print(label+': all physical fields/checkpoints exact', flush=True)
                        if moving:
                            cfg['sim'].update(is_restarting='yes', restarting_from_modelname=str(reference/'result'),
                                              restarting_from_frame='2')
                            restart, resumed = run(label+'-restart', cfg, executables['mpi'+suffix], ranks, overlap)
                            compare_all(restart_reference, restart, (3, 4))
                            report['checks'].append(dict(case=label+'-restart', bit_exact=True, timing=resumed,
                                                         frames=[3, 4], checkpoint_frames=[3, 4],
                                                         checkpoint_reference=str(restart_reference)))
                            save()
                            print(label+'-restart: continuation exact', flush=True)
                            cfg['sim']['is_restarting'] = 'no'
    except Exception as error:
        report['failure'] = str(error)
        save()
        raise
    report.update(completed=True, all_passed=True)
    save()
    print('Passed {} comparisons: {}'.format(len(report['checks']), target), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
