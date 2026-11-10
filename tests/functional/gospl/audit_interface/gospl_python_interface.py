"""Transparent audit of the real GoSPL bridge for the DES coupling smoke test.

This module forwards every operation to the installed bridge. It records the
boundary values and call durations without replacing the surface model.
"""
import importlib.util
import json
import os
from pathlib import Path
import time
import numpy as np

_spec=importlib.util.spec_from_file_location('_coupling_real_interface',os.environ['GOSPL_REAL_INTERFACE'])
_real=importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_real)
for _name,_value in vars(_real).items():
    if not _name.startswith('__'):globals()[_name]=_value
_audit=Path(os.environ['GOSPL_AUDIT_DIR']);_audit.mkdir(parents=True,exist_ok=True)
_events=0
_last_velocity=None

def _record(row):
    with (_audit/'calls.jsonl').open('a') as stream:stream.write(json.dumps(row)+'\n')

def create_enhanced_model(*args,**kwargs):
    start=time.perf_counter();handle=_real.create_enhanced_model(*args,**kwargs)
    from mpi4py import MPI
    _record(dict(operation='create',handle=handle,mpi_size=MPI.COMM_WORLD.size,call_s=time.perf_counter()-start))
    return handle

def apply_elevation_data(*args,**kwargs):
    start=time.perf_counter();result=_real.apply_elevation_data(*args,**kwargs)
    elapsed=time.perf_counter()-start
    audit_start=time.perf_counter()
    if result==0:
        model=_real._models[args[0]]
        np.savez(str(_audit/'seed.npz'),coords=np.asarray(args[1]).reshape(args[3],3),
                 elevation=np.asarray(args[2]),owned_coords=model.mCoords[model.glbIDs],
                 owned_height=model.hGlobal.getArray().copy(),local_coords=model.mCoords[model.locIDs],
                 local_height=model.hLocal.getArray().copy())
    _record(dict(operation='seed',status=int(result),call_s=elapsed,
                 audit_write_s=time.perf_counter()-audit_start))
    return result

def set_surface_velocity(handle,coords,vx,vy,vz,num_points,*args,**kwargs):
    global _last_velocity
    _last_velocity=np.column_stack((vx,vy,vz)).copy()
    start=time.perf_counter()
    result=_real.set_surface_velocity(handle,coords,vx,vy,vz,num_points,*args,**kwargs)
    _record(dict(operation='velocity',status=int(result),points=int(num_points),call_s=time.perf_counter()-start))
    return result

def run_and_get_erosion(handle,dt,coords,num_points,*args,**kwargs):
    global _events
    points=np.asarray(coords).reshape(num_points,3).copy()
    before=float(_real._models[handle].tNow)
    start=time.perf_counter()
    result=_real.run_and_get_erosion(handle,dt,coords,num_points,*args,**kwargs)
    elapsed=time.perf_counter()-start
    after=float(_real._models[handle].tNow)
    _events+=1
    audit_start=time.perf_counter()
    if result is not None:
        np.savez(str(_audit/('event-%04d.npz'%_events)),coords=points,velocity=_last_velocity,
                 dh=np.asarray(result).copy(),dt=float(dt),time_before=before,time_after=after)
    _record(dict(operation='step',event=_events,dt=float(dt),time_before=before,time_after=after,
                 points=int(num_points),call_s=elapsed,returned=result is not None,
                 audit_write_s=time.perf_counter()-audit_start))
    return result

def destroy_model(handle):
    start=time.perf_counter();result=_real.destroy_model(handle)
    _record(dict(operation='destroy',status=int(result),remaining_models=len(_real._models),call_s=time.perf_counter()-start))
    return result
