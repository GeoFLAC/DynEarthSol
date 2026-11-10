#!/usr/bin/env python3
"""Launch DES and a persistent GoSPL worker with independent MPI rank counts.

Run inside the scheduler allocation where required. Each launcher must stay
within its reserved slots. The DES root and GoSPL root must run on this host;
the socket carries boundary arrays and is never used for timestep halos.
"""
import argparse
import json
import os
from pathlib import Path
import secrets
import shlex
import signal
import subprocess
import tempfile
import time


def main():
    def terminate(signum, frame):
        raise SystemExit(128 + signum)
    signal.signal(signal.SIGTERM, terminate)
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gospl-python', type=Path, required=True)
    p.add_argument('--interface', type=Path, required=True)
    p.add_argument('--gospl-launcher', default='', help='Complete prefix, e.g. /env/bin/mpiexec -n 2')
    p.add_argument('--worker-env', action='append', default=[], metavar='NAME=VALUE')
    p.add_argument('--audit-dir', type=Path)
    p.add_argument('--workdir', type=Path, required=True)
    p.add_argument('--timeout', type=float, default=300)
    p.add_argument('--run-timeout', type=float, default=0, help='Total DES runtime limit; zero means unlimited')
    p.add_argument('des_command', nargs=argparse.REMAINDER)
    args = p.parse_args()
    command = args.des_command
    if command[:1] == ['--']:
        command = command[1:]
    if not command or args.timeout <= 0 or args.run_timeout < 0:
        p.error('Provide a DES command after -- and a positive timeout')
    folder = args.workdir.resolve()
    folder.mkdir(parents=True, exist_ok=True)
    # Private local endpoint, independent of model storage and pathname length.
    with tempfile.TemporaryDirectory(prefix='des-gospl-') as ipc:
        endpoint = str(Path(ipc) / 'worker.sock')
        env = dict(os.environ, DES_GOSPL_ENDPOINT=endpoint,
                   DES_GOSPL_TOKEN=secrets.token_hex(32), DES_GOSPL_TIMEOUT=str(args.timeout))
        module_dir = str(Path(__file__).resolve().parent)
        env['PYTHONPATH'] = module_dir + os.pathsep + env.get('PYTHONPATH', '')
        worker_env = dict(env)
        for setting in args.worker_env:
            key, sep, value = setting.partition('=')
            if not sep or not key:
                p.error('--worker-env requires NAME=VALUE')
            worker_env[key] = value
        worker_command = shlex.split(args.gospl_launcher) + [str(args.gospl_python.resolve()),
            str(Path(__file__).with_name('mpi_worker.py').resolve()), '--endpoint', endpoint,
            '--interface', str(args.interface.resolve())]
        if args.audit_dir:
            worker_command += ['--audit-dir', str(args.audit_dir.resolve())]
        record = dict(des_command=command, worker_command=worker_command,
                      started=time.time(), completed=False)
        target = folder / 'coupling-launch.json'
        target.write_text(json.dumps(record, indent=2) + '\n')
        with (folder / 'worker.log').open('w') as log:
            worker = subprocess.Popen(worker_command, env=worker_env, stdout=log,
                                      stderr=subprocess.STDOUT, start_new_session=True)
            des = None
            try:
                deadline = time.monotonic() + args.timeout
                while not Path(endpoint + '.ready').exists():
                    if worker.poll() is not None:
                        raise RuntimeError('GoSPL worker exited; see ' + str(folder / 'worker.log'))
                    if time.monotonic() >= deadline:
                        raise TimeoutError('GoSPL worker startup timed out')
                    time.sleep(.05)
                record['worker'] = json.loads(Path(endpoint + '.ready').read_text())
                with (folder / 'des.log').open('w') as des_log:
                    des = subprocess.Popen(command, env=env, stdout=des_log,
                                           stderr=subprocess.STDOUT, start_new_session=True)
                    run_started = time.monotonic()
                    while des.poll() is None:
                        if args.run_timeout and time.monotonic()-run_started > args.run_timeout:
                            raise TimeoutError('Coupled DES run exceeded --run-timeout')
                        if worker.poll() not in (None, 0):
                            raise RuntimeError('GoSPL worker exited before DES completed')
                        time.sleep(.05)
                if des.returncode:
                    raise RuntimeError('DES exited with code %s; see des.log' % des.returncode)
                # Python finalization closes the client and collectively stops
                # the independent worker world without an extra connection.
                worker.wait(timeout=min(args.timeout, 30))
                if worker.returncode:
                    raise RuntimeError('GoSPL worker exited with code %s' % worker.returncode)
                record.update(completed=True, des_exit=des.returncode, worker_exit=worker.returncode)
            finally:
                for process in (des, worker):
                    if process is not None and process.poll() is None:
                        if process is worker:
                            # A reported RPC error is already unwinding all
                            # worker ranks; allow collective cleanup to finish.
                            try:
                                process.wait(timeout=3)
                                continue
                            except subprocess.TimeoutExpired:
                                # The grace period ended; terminate the whole
                                # worker process group below, then reap it.
                                pass
                        os.killpg(process.pid, signal.SIGTERM)
                        try:
                            process.wait(timeout=5)
                        except subprocess.TimeoutExpired:
                            os.killpg(process.pid, signal.SIGKILL)
                            process.wait()
                stats = Path(endpoint + '.jsonl')
                if stats.exists():
                    (folder / 'worker-calls.jsonl').write_bytes(stats.read_bytes())
                record['ended'] = time.time()
                record['des_exit'] = des.returncode if des is not None else None
                record['worker_exit'] = worker.returncode
                target.write_text(json.dumps(record, indent=2) + '\n')


if __name__ == '__main__':
    main()
