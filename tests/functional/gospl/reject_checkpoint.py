#!/usr/bin/env python3
"""Require controlled rejection of damaged/incompatible coupling checkpoints.

Only private copies are modified. The retained continuous fixture is unchanged.
"""
import argparse
import configparser
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[3]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('continuous', 'exe', 'extension', 'gospl-python', 'workdir'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--launcher', required=True)
    parser.add_argument('--gospl-launcher', required=True)
    parser.add_argument('--frame', type=int, default=5)
    parser.add_argument('--ranks', type=int, default=3)
    parser.add_argument('--gospl-ranks', type=int, default=2)
    args = parser.parse_args()
    if args.gospl_ranks < 2:
        parser.error('Use a fixture with at least two GoSPL ranks')
    args.workdir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='checkpoint-rejection-', dir=args.workdir)).resolve()
    print('Artifacts:', work, flush=True)
    checks = []
    for name, expected in (('checksum', 'checksum mismatch'), ('rank-count', 'same MPI rank count')):
        folder = work/name
        folder.mkdir()
        seed = folder/'source'
        seed.mkdir()
        checkpoint = 'result.chkpt.%06d.vtkhdf' % args.frame
        for file in (checkpoint, 'result.save.%06d.vtkhdf' % args.frame, 'result.info'):
            shutil.copy2(args.continuous/file, seed/file)
        shutil.copytree(args.continuous/(checkpoint+'.gospl'), seed/(checkpoint+'.gospl'))
        original = None
        if name == 'checksum':
            target = next((seed/(checkpoint+'.gospl')).glob('*/rank-0000.pkl'))
            original = args.continuous/target.relative_to(seed)
            before = hashlib.sha256(original.read_bytes()).hexdigest()
            payload = bytearray(target.read_bytes())
            payload[-1] ^= 1
            target.write_bytes(payload)
        config = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
        config.optionxform = str
        config.read(args.continuous/'model.cfg')
        config['sim'].update(is_restarting='yes', restarting_from_modelname=str(seed/'result'),
                             restarting_from_frame=str(args.frame))
        with (folder/'model.cfg').open('w') as stream:
            config.write(stream)
        shutil.copy2(args.continuous/'gospl.yml', folder/'gospl.yml')
        count = args.gospl_ranks if name == 'checksum' else args.gospl_ranks-1
        command = [sys.executable, str(ROOT/'gospl_driver/run_coupled.py'),
                   '--gospl-python', str(args.gospl_python.resolve()), '--interface',
                   str(args.extension.resolve()/'cpp_interface/gospl_python_interface.py'),
                   '--gospl-launcher', args.gospl_launcher+' -n '+str(count),
                   '--workdir', str(folder), '--timeout', '60', '--run-timeout', '90', '--']
        command += shlex.split(args.launcher)+['-n', str(args.ranks), str(args.exe.resolve()), 'model.cfg']
        with (folder/'driver.log').open('w') as log:
            result = subprocess.run(command, cwd=folder, stdout=log, stderr=subprocess.STDOUT,
                                    env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'), timeout=120)
        record = json.loads((folder/'coupling-launch.json').read_text())
        if result.returncode == 0 or record['completed'] or record['worker_exit'] != 53:
            raise RuntimeError('Checkpoint was not rejected cleanly: '+name)
        if expected not in (folder/'des.log').read_text()+(folder/'worker.log').read_text():
            raise RuntimeError('Wrong rejection reason: '+name)
        if original and hashlib.sha256(original.read_bytes()).hexdigest() != before:
            raise RuntimeError('Original fixture changed')
        checks.append(dict(case=name, expected=expected, des_exit=record['des_exit'], worker_exit=record['worker_exit']))
        print(name, 'rejection PASS', flush=True)
    (work/'results.json').write_text(json.dumps(dict(completed=True, all_passed=True, checks=checks), indent=2)+'\n')


if __name__ == '__main__':
    main()
