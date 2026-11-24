#!/usr/bin/env python3
"""Build/run shared halo kernels explicitly in OpenACC host or NVIDIA GPU mode.

GPU mode fails if real devices/CUDA-aware MPI are unavailable. Host mode is not
GPU verification. This tests halo infrastructure, not a complete GPU DES solver.
"""
import argparse
import hashlib
import json
import os
import re
import signal
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[3]


def execute(command, logfile, timeout):
    with logfile.open('w') as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        try:
            code = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            raise RuntimeError(f'Timed out: {command}; see {logfile}')
    return code, logfile.read_text()


def identities(log, ranks, mode, min_nodes):
    rows = [dict(re.findall(r'(\w+)=(\S+)', line)) for line in
            re.findall(r'\[accelerator\]\[rank\] (.*)', log)]
    if len(rows) != ranks or {int(row['rank']) for row in rows} != set(range(ranks)):
        raise RuntimeError('Missing or duplicate MPI rank identity records')
    nodes = {row['host'] for row in rows}
    if len(nodes) < min_nodes:
        raise RuntimeError(f'Expected at least {min_nodes} nodes, got {sorted(nodes)}')
    for node in nodes:
        local = [row for row in rows if row['host'] == node]
        if ({int(row['local_rank']) for row in local} != set(range(len(local))) or
                any(int(row['local_size']) != len(local) for row in local)):
            raise RuntimeError(f'Inconsistent local rank mapping on {node}')
        if mode == 'gpu' and (any(not re.fullmatch(r'[0-9a-f]{32}', row['uuid']) for row in local) or
                              len({row['uuid'] for row in local}) != len(local)):
            raise RuntimeError(f'Physical GPU identities missing or duplicated on {node}')
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=['host','gpu'],default='host')
    p.add_argument('--compiler',default='mpicxx')
    p.add_argument('--std',choices=['c++17','c++20'],default='c++17')
    p.add_argument('--acc-flags',help='Default: -fopenacc (host), -acc -cuda (gpu)')
    p.add_argument('--launcher',default='mpiexec --bind-to none')
    p.add_argument('--ranks',nargs='+',type=int,default=[1,2,3,4,6,8])
    p.add_argument('--workdir',type=Path,required=True)
    p.add_argument('--min-nodes',type=int,default=1,help='Require this many distinct processor names per run')
    p.add_argument('--timeout',type=int,default=120)
    p.add_argument('--prebuilt-dir',type=Path,help='Use existing halo-2d/halo-3d binaries without a compiler')
    a=p.parse_args()
    if min(a.ranks)<1 or len(set(a.ranks))!=len(a.ranks):p.error('Ranks must be distinct and positive')
    if a.min_nodes<1 or min(a.ranks)<a.min_nodes:p.error('Every rank count must cover --min-nodes')
    if a.timeout<1:p.error('Timeout must be positive')
    if a.acc_flags is None:a.acc_flags='-acc -cuda' if a.mode=='gpu' else '-fopenacc'
    a.workdir.mkdir(parents=True,exist_ok=True)
    output=Path(tempfile.mkdtemp(prefix='accelerator-',dir=a.workdir));print('Artifacts:',output,flush=True)
    commands={name:shutil.which(name) for name in ['nvc++','nvcc','nvidia-smi','sycl-ls','icpx','mpicxx','ompi_info']}
    report=dict(completed=False,gpu_verified=False,mode=a.mode,cpp_standard=a.std,tools=commands,
                environment={name:os.environ[name] for name in
                             ['PATH','LD_LIBRARY_PATH','OPAL_PREFIX','OMPI_CXX','UCX_TLS',
                              'UCX_MEMTYPE_CACHE','UCX_CUDA_IPC_CACHE','CUDA_VISIBLE_DEVICES',
                              'OMPI_MCA_pml','OMPI_MCA_btl','OMPI_MCA_coll_hcoll_enable']
                             if name in os.environ},
                local_source_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                               for path in [ROOT/'tests/functional/mpi/accelerator_halo.cxx',
                                            ROOT/'mpi_halo_kernels.hpp',ROOT/'mpi_transport.hpp']},
                builds=[], min_nodes=a.min_nodes, device_paths={path:Path(path).exists() for path in ['/dev/dri','/dev/dxg','/dev/nvidia0']},checks=[])
    if commands['nvidia-smi']:
        r=subprocess.run([commands['nvidia-smi'],'topo','-m'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True,timeout=20)
        report['nvidia_topology']=r.stdout;report['nvidia_status']=r.returncode
    if commands['ompi_info']:
        r=subprocess.run([commands['ompi_info'],'--parsable','--all'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True,timeout=20)
        report['mpi_cuda_built']='mpi_built_with_cuda_support:value:true' in r.stdout
    (output/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    try:
        for dim in (2,3):
            exe=output/f'halo-{dim}d'
            command=shlex.split(a.compiler)+['-std='+a.std,'-O2','-DOMPI_SKIP_MPICXX=1','-I'+str(ROOT)]
            command+=shlex.split(a.acc_flags)+(['-DTHREED'] if dim==3 else [])
            if a.mode=='gpu':command+=['-DDES_TEST_CUDA']
            command+=[str(ROOT/'tests/functional/mpi/accelerator_halo.cxx'),'-o',str(exe)]
            report['builds'].append(dict(ndims=dim,command=command if not a.prebuilt_dir else None))
            (output/'results.json').write_text(json.dumps(report,indent=2)+'\n')
            if a.prebuilt_dir:
                source=(a.prebuilt_dir/f'halo-{dim}d').resolve()
                shutil.copy2(source,exe)
                report['builds'][-1]['prebuilt_source']=str(source)
            else:
                code,_=execute(command,output/f'build-{dim}d.log',180)
                if code:raise RuntimeError(f'Build failed: {output}/build-{dim}d.log')
            report['builds'][-1]['executable_sha256']=hashlib.sha256(exe.read_bytes()).hexdigest()
            for ranks in a.ranks:
                command=shlex.split(a.launcher)+['-n',str(ranks),str(exe)]+(['--host'] if a.mode=='host' else [])
                code,log=execute(command,output/f'run-{dim}d-{ranks}.log',a.timeout)
                passed=code==0 and 'PASS:' in log
                row=dict(ndims=dim,ranks=ranks,command=command,passed=passed,returncode=code)
                if passed:
                    try:row['identities']=identities(log,ranks,a.mode,a.min_nodes)
                    except (RuntimeError,KeyError,ValueError) as error:
                        row.update(passed=False,error=str(error));passed=False
                report['checks'].append(row)
                (output/'results.json').write_text(json.dumps(report,indent=2)+'\n')
                if not passed:
                    report['blocked_reason']='GPU prerequisites unavailable or accelerator test failed; inspect retained log'
                    (output/'results.json').write_text(json.dumps(report,indent=2)+'\n')
                    print(report['blocked_reason'],flush=True);return code or 1
                print(dim,ranks,a.mode,'PASS',flush=True)
    except (OSError, RuntimeError, subprocess.SubprocessError) as error:
        report['failure']=str(error)
        (output/'results.json').write_text(json.dumps(report,indent=2)+'\n')
        print(str(error),flush=True)
        return 1
    report['completed']=True;report['gpu_verified']=a.mode=='gpu'
    (output/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    return 0


if __name__=='__main__':raise SystemExit(main())
