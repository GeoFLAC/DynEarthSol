#!/usr/bin/env python3
"""Build/run shared halo kernels explicitly in OpenACC host or NVIDIA GPU mode.

GPU mode fails if real devices/CUDA-aware MPI are unavailable. Host mode is not
GPU verification. This tests halo infrastructure, not a complete GPU DES solver.
"""
import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[3]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=['host','gpu'],default='host')
    p.add_argument('--compiler',default='mpicxx')
    p.add_argument('--acc-flags',default='-fopenacc')
    p.add_argument('--launcher',default='mpiexec --bind-to none')
    p.add_argument('--ranks',nargs='+',type=int,default=[1,2,3,4,6,8])
    p.add_argument('--workdir',type=Path,required=True)
    a=p.parse_args()
    if min(a.ranks)<1:p.error('Ranks must be positive')
    a.workdir.mkdir(parents=True,exist_ok=True)
    output=Path(tempfile.mkdtemp(prefix='accelerator-',dir=a.workdir));print('Artifacts:',output,flush=True)
    commands={name:shutil.which(name) for name in ['nvc++','nvcc','nvidia-smi','sycl-ls','icpx','mpicxx','ompi_info']}
    report=dict(completed=False,gpu_verified=False,mode=a.mode,tools=commands,
                device_paths={path:Path(path).exists() for path in ['/dev/dri','/dev/dxg','/dev/nvidia0']},checks=[])
    if commands['nvidia-smi']:
        r=subprocess.run([commands['nvidia-smi'],'topo','-m'],capture_output=True,text=True,timeout=20)
        report['nvidia_topology']=r.stdout;report['nvidia_status']=r.returncode
    if commands['ompi_info']:
        r=subprocess.run([commands['ompi_info'],'--parsable','--all'],capture_output=True,text=True,timeout=20)
        report['mpi_cuda_built']='mpi_built_with_cuda_support:value:true' in r.stdout
    (output/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    for dim in (2,3):
        exe=output/f'halo-{dim}d'
        command=shlex.split(a.compiler)+['-std=c++11','-O2','-DOMPI_SKIP_MPICXX=1','-I'+str(ROOT)]
        command+=shlex.split(a.acc_flags)+(['-DTHREED'] if dim==3 else [])
        command+=[str(ROOT/'tests/functional/mpi/accelerator_halo.cxx'),'-o',str(exe)]
        with (output/f'build-{dim}d.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
        for ranks in a.ranks:
            command=shlex.split(a.launcher)+['-n',str(ranks),str(exe)]+(['--host'] if a.mode=='host' else [])
            r=subprocess.run(command,env=os.environ,capture_output=True,text=True,timeout=120)
            (output/f'run-{dim}d-{ranks}.log').write_text(r.stdout+r.stderr)
            passed=r.returncode==0 and 'PASS:' in r.stdout
            row=dict(ndims=dim,ranks=ranks,command=command,passed=passed,returncode=r.returncode)
            report['checks'].append(row)
            if not passed:
                report['blocked_reason']='GPU prerequisites unavailable or accelerator test failed; inspect retained log'
                (output/'results.json').write_text(json.dumps(report,indent=2)+'\n')
                print(report['blocked_reason'],flush=True);return r.returncode or 1
            print(dim,ranks,a.mode,'PASS',flush=True)
    report['completed']=True;report['gpu_verified']=a.mode=='gpu'
    (output/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    return 0


if __name__=='__main__':raise SystemExit(main())
