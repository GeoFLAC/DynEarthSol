#!/usr/bin/env python3
"""Launch a CPU DES check across explicitly named SSH hosts without a shared cwd.

Install the same executable, dependencies and input assets at identical paths on
each host first. This stages cfg/restart inputs and collects rank-zero outputs.
Use as run_mpi.py's --launcher. Requires passwordless, host-key-verified SSH.
"""
import argparse
import configparser
import hashlib
import json
import os
from pathlib import Path
import shlex
import socket
import subprocess
import sys


def rank_main(command):
    rank = int(os.environ['OMPI_COMM_WORLD_RANK'])
    row = dict(rank=rank, host=socket.gethostname(),
               local_rank=int(os.environ['OMPI_COMM_WORLD_LOCAL_RANK']),
               threads=int(os.environ.get('OMP_NUM_THREADS', '1')),
               cpus=sorted(os.sched_getaffinity(0)))
    print('[cluster][rank] ' + json.dumps(row, sort_keys=True), flush=True)
    os.execv(command[0], command)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == '--rank':
        rank_main(sys.argv[2:])
        return 0
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hosts', required=True, help='Comma-separated SSH host names or IPs')
    p.add_argument('--mpiexec', default='mpiexec')
    p.add_argument('--python', default='python3', help='Python available on every worker')
    p.add_argument('--known-hosts', type=Path, help='Optional existing pinned SSH host-key file')
    p.add_argument('--min-nodes', type=int, default=2)
    p.add_argument('--map-by', choices=['node', 'seq'], default='node',
                   help='seq assigns contiguous rank groups in the supplied host order')
    p.add_argument('-n', type=int, required=True)
    p.add_argument('executable', type=Path)
    p.add_argument('arguments', nargs=argparse.REMAINDER)
    a = p.parse_args()
    hosts = a.hosts.split(',')
    if len(set(hosts)) != len(hosts) or a.min_nodes < 1 or a.n < a.min_nodes:
        p.error('Use distinct hosts and enough ranks for --min-nodes')
    if len(hosts) < a.min_nodes:
        p.error('Too few hosts for --min-nodes')
    if len(a.arguments) != 1:
        p.error('Supply one cfg path or --help after the executable')
    executable = a.executable.resolve()
    cfg = None if a.arguments[0] in ('--help', '-h') else Path(a.arguments[0]).resolve()
    script = Path(__file__).resolve()
    cwd = Path.cwd()
    ssh = ['ssh', '-x', '-o', 'BatchMode=yes', '-o', 'StrictHostKeyChecking=yes',
           '-o', 'UpdateHostKeys=no', '-o', 'ConnectTimeout=8']
    if a.known_hosts:
        ssh += ['-o', 'UserKnownHostsFile=' + str(a.known_hosts.resolve())]
    digest = hashlib.sha256(executable.read_bytes()).hexdigest()
    # Verify workers execute the same solver. Keep remote setup outside MPI so
    # failures are explicit instead of leaving ranks stuck in initialization.
    host_aliases = {}
    for host in hosts:
        remote = 'mkdir -p ' + shlex.quote(str(cwd)) + ' && '
        if cfg:
            remote += 'cat > ' + shlex.quote(str(cfg)) + ' && '
        remote += 'sha256sum ' + shlex.quote(str(executable)) + ' && hostname'
        result = subprocess.run(ssh + [host, remote], input=cfg.read_bytes() if cfg else b'',
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30)
        lines = result.stdout.decode().splitlines()
        if result.returncode or len(lines) != 2 or lines[0].split()[0] != digest:
            raise RuntimeError(f'Failed cfg staging or executable hash check on {host}: '
                               + result.stderr.decode() + result.stdout.decode())
        host_aliases[lines[1]] = host
    restart_files = []
    if cfg:
        settings = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
        settings.read(str(cfg))
        if settings.getboolean('sim', 'is_restarting', fallback=False):
            prefix = Path(settings['sim']['restarting_from_modelname']).resolve()
            restart_files = sorted(prefix.parent.glob(prefix.name + '.*'))
            if not restart_files:
                raise RuntimeError(f'No restart inputs found for {prefix}')
            for host in hosts:
                subprocess.run(ssh + [host, 'mkdir -p ' + shlex.quote(str(prefix.parent))],
                               stdin=subprocess.DEVNULL, check=True, timeout=30)
                subprocess.run(['rsync', '-a', '-e', ' '.join(map(shlex.quote, ssh))] +
                               [str(path) for path in restart_files] +
                               [host + ':' + str(prefix.parent) + '/'], stdin=subprocess.DEVNULL,
                               check=True, timeout=120)
    allocation = ['--host', ','.join(hosts)]
    if a.map_by == 'seq':
        sequence = [host for i, host in enumerate(hosts)
                    for _ in range(a.n // len(hosts) + (i < a.n % len(hosts)))]
        allocation = ['--host', ','.join(sequence)]
    command = [a.mpiexec] + allocation + ['--map-by', a.map_by, '--bind-to', 'none',
               '--mca', 'plm_rsh_agent', ' '.join(map(shlex.quote, ssh)),
               '--wdir', str(cwd)]
    for name in ('PATH', 'LD_LIBRARY_PATH', 'OPAL_PREFIX', 'UCX_TLS', 'OMPI_MCA_pml',
                 'OMPI_MCA_coll_hcoll_enable', 'OMP_NUM_THREADS', 'OMP_DYNAMIC',
                 'DES_MPI_OVERLAP', 'DES_MPI_PARTITIONER'):
        if name in os.environ:
            command += ['-x', name]
    command += ['-n', str(a.n), a.python, str(script), '--rank', str(executable),
                str(cfg) if cfg else a.arguments[0]]
    report = dict(completed=False, passed=False, command=command, hosts=hosts,
                  executable_sha256=digest, min_nodes=a.min_nodes, identities=[],
                  cfg_sha256=hashlib.sha256(cfg.read_bytes()).hexdigest() if cfg else None,
                  restart_inputs={str(path):hashlib.sha256(path.read_bytes()).hexdigest()
                                  for path in restart_files})
    target = cwd / 'cluster-launch.json'
    target.write_text(json.dumps(report, indent=2) + '\n')
    process = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               universal_newlines=True, bufsize=1)
    for line in process.stdout:
        print(line, end='', flush=True)
        if line.startswith('[cluster][rank] '):
            report['identities'].append(json.loads(line[len('[cluster][rank] '):]))
    code = process.wait()
    rows = report['identities']
    report.update(completed=True, returncode=code)
    valid = (len(rows) == a.n and {row['rank'] for row in rows} == set(range(a.n)) and
             len({row['host'] for row in rows}) >= a.min_nodes)
    for host in {row['host'] for row in rows}:
        local = [row for row in rows if row['host'] == host]
        valid = valid and {row['local_rank'] for row in local} == set(range(len(local)))
    roots = [row for row in rows if row['rank'] == 0]
    if len(roots) == 1 and roots[0]['host'] != socket.gethostname() and cfg:
        source_host = host_aliases[roots[0]['host']]
        collect = ['rsync', '-a', '-e', ' '.join(map(shlex.quote, ssh)),
                   '--exclude=cluster-launch.json', '--exclude=*.log',
                   source_host + ':' + str(cwd) + '/', str(cwd) + '/']
        result = subprocess.run(collect, stdin=subprocess.DEVNULL, timeout=120)
        report['output_collection'] = dict(command=collect, returncode=result.returncode)
        if result.returncode:
            code = code or result.returncode
    report['passed'] = code == 0 and valid
    report['returncode'] = code
    target.write_text(json.dumps(report, indent=2) + '\n')
    if not valid:
        print('FAIL: missing/duplicate ranks or insufficient distinct nodes', file=sys.stderr)
    return code or (0 if valid else 1)


if __name__ == '__main__':
    raise SystemExit(main())
