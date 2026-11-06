#!/usr/bin/env python3
"""Bind a benchmark process, record its actual CPU mask, then exec the solver."""
import argparse
import json
import os
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--cpus',required=True)
p.add_argument('--mpi',action='store_true')
p.add_argument('--record',type=Path,required=True)
p.add_argument('--exe',required=True)
p.add_argument('--config',default='model.cfg')
a=p.parse_args();cpus=[int(c) for c in a.cpus.split(',')]
rank=0;local_rank=0
if a.mpi:
    rank=int(os.environ['OMPI_COMM_WORLD_RANK'])
    local_rank=int(os.environ['OMPI_COMM_WORLD_LOCAL_RANK'])
    if local_rank>=len(cpus):raise ValueError('More local ranks than benchmark CPUs')
    selected={cpus[local_rank]}
else:selected=set(cpus)
os.sched_setaffinity(0,selected)
actual=sorted(os.sched_getaffinity(0))
if actual!=sorted(selected):raise RuntimeError('CPU affinity did not take effect')
(a.record/f'binding-{rank}.json').write_text(json.dumps(dict(rank=rank,cpus=actual))+'\n')
os.execv(a.exe,[a.exe,a.config])
