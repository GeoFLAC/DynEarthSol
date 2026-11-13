#!/usr/bin/env python3
"""Launch MPI checks inside an existing Grid Engine allocation.

Use as run_mpi.py's --launcher from a qsub job. The source, executable, working
directory and dependencies must be visible at the same paths on allocated hosts.
Open MPI discovers its hosts and starts tasks through its Grid Engine support.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import socket
import subprocess
import sys


def host_key(name):
    return name.split('.')[0].lower()


def allocation():
    if not all(os.environ.get(key) for key in ('JOB_ID', 'NSLOTS', 'PE_HOSTFILE')):
        raise RuntimeError('Submit through qsub: JOB_ID, NSLOTS and PE_HOSTFILE are required')
    slots = {}
    raw = Path(os.environ['PE_HOSTFILE']).read_text()
    for line in raw.splitlines():
        if not line.strip():
            continue
        fields = line.split()
        host, count = host_key(fields[0]), int(fields[1])
        if count < 1:
            raise RuntimeError('Nonpositive slot count in PE_HOSTFILE')
        slots[host] = slots.get(host, 0) + count
    if sum(slots.values()) != int(os.environ['NSLOTS']):
        raise RuntimeError('PE_HOSTFILE slot total does not match NSLOTS')
    return dict(job_id=os.environ['JOB_ID'], pe=os.environ.get('PE'),
                nslots=int(os.environ['NSLOTS']), host_slots=slots, pe_hostfile=raw)


def rank_main(args):
    context = json.loads(Path(args[0]).read_text())
    command = args[1:]
    host = host_key(socket.gethostname())
    if host not in context['host_slots']:
        raise RuntimeError('MPI rank started outside the scheduler allocation: ' + host)
    if os.environ.get('JOB_ID') != context['job_id']:
        raise RuntimeError('MPI worker job ID differs from the submitting job')
    rank = int(os.environ['OMPI_COMM_WORLD_RANK'])
    local_rank = int(os.environ['OMPI_COMM_WORLD_LOCAL_RANK'])
    local_size = int(os.environ['OMPI_COMM_WORLD_LOCAL_SIZE'])
    threads = context['threads']
    if local_size * threads > context['host_slots'][host]:
        raise RuntimeError('MPI/OpenMP workers exceed allocated slots on ' + host)
    digest = hashlib.sha256(Path(command[0]).read_bytes()).hexdigest()
    if digest != context['executable_sha256']:
        raise RuntimeError('MPI worker executable differs from the submitting executable')
    record = dict(rank=rank, host=host, local_rank=local_rank, local_size=local_size,
                  threads=threads, job_id=context['job_id'], cpus=sorted(os.sched_getaffinity(0)))
    os.write(1, ('[sge][rank] ' + json.dumps(record, sort_keys=True) + '\n').encode())
    os.execv(command[0], command)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == '--rank':
        rank_main(sys.argv[2:])
        return 0
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mpiexec', default='mpiexec')
    parser.add_argument('--python', default=sys.executable)
    parser.add_argument('--min-nodes', type=int, default=1)
    parser.add_argument('--local-only', action='store_true',
                        help='Use an allocated single host with an MPI build lacking Grid Engine support')
    parser.add_argument('-n', type=int, required=True)
    parser.add_argument('executable', type=Path)
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    context = allocation()
    if args.local_only and (len(context['host_slots'])!=1 or
                           host_key(socket.gethostname()) not in context['host_slots']):
        parser.error('--local-only requires this host to be the sole scheduler-allocated host')
    threads = int(os.environ.get('OMP_NUM_THREADS', '1'))
    os.environ.setdefault('OMP_NUM_THREADS', str(threads))
    if threads < 1 or args.n < 1 or args.min_nodes < 1 or args.min_nodes > args.n:
        parser.error('Use positive rank/thread counts and a reachable node minimum')
    capacity = sum(count // threads for count in context['host_slots'].values())
    if args.n > capacity:
        parser.error('Requested MPI ranks times OpenMP threads exceed the allocated slot capacity')
    if args.min_nodes > len(context['host_slots']):
        parser.error('The scheduler allocated too few distinct hosts')
    executable = args.executable.resolve()
    context.update(threads=threads, executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest())
    context_path = Path.cwd() / 'sge-allocation.json'
    context_path.write_text(json.dumps(context, indent=2) + '\n')
    command = [args.mpiexec, '--bind-to', 'none']
    if args.local_only:
        command += ['--mca','ras','^gridengine','--host',
                    '{}:{}'.format(socket.gethostname(),context['nslots'])]
    else:
        command += ['--mca','ras','gridengine']
    for name in ('JOB_ID', 'OMP_NUM_THREADS', 'OMP_DYNAMIC', 'OPENBLAS_NUM_THREADS',
                 'PATH', 'LD_LIBRARY_PATH', 'PYTHONHOME', 'PYTHONPATH',
                 'DES_MPI_OVERLAP', 'DES_MPI_PARTITIONER'):
        if name in os.environ:
            command += ['-x', name]
    command += ['-n', str(args.n),
               args.python, str(Path(__file__).resolve()), '--rank', str(context_path),
               str(executable)] + args.arguments
    report = dict(completed=False, passed=False, allocation=context, command=command,
                  min_nodes=args.min_nodes, identities=[])
    report_path = Path.cwd() / 'sge-launch.json'
    report_path.write_text(json.dumps(report, indent=2) + '\n')
    process = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, universal_newlines=True, bufsize=1)
    for line in process.stdout:
        print(line, end='', flush=True)
        if line.startswith('[sge][rank] '):
            report['identities'].append(json.loads(line[len('[sge][rank] '):]))
    code = process.wait()
    rows = report['identities']
    valid = (len(rows) == args.n and {row['rank'] for row in rows} == set(range(args.n)) and
             len({row['host'] for row in rows}) >= args.min_nodes)
    for host in {row['host'] for row in rows}:
        local = [row for row in rows if row['host'] == host]
        valid = valid and host in context['host_slots']
        valid = valid and {row['local_rank'] for row in local} == set(range(len(local)))
        valid = valid and all(row['local_size'] == len(local) for row in local)
        valid = valid and len(local) * threads <= context['host_slots'].get(host, 0)
    report.update(completed=True, passed=code == 0 and valid, returncode=code)
    report_path.write_text(json.dumps(report, indent=2) + '\n')
    if not valid:
        print('FAIL: MPI placement differs from the scheduler allocation', file=sys.stderr)
    return code or (0 if valid else 1)


if __name__ == '__main__':
    raise SystemExit(main())
