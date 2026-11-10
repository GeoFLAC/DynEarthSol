"""Lossless, rank-local GoSPL continuation state for the coupling worker.

PETSc vectors/matrices are stored as numeric arrays. Runtime resources are
reconstructed from the same model configuration; a changed MPI partition is
rejected before applying state. Each immutable generation is committed before
its token is recorded in the DES checkpoint.
"""
import hashlib
import json
import os
from pathlib import Path
import pickle
import shutil
import types
import uuid
from collections.abc import Mapping
import numpy as np
from petsc4py import PETSc


def _sha(value):
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()


def _configuration(model):
    def plain(value):
        if isinstance(value,Mapping):return {str(k):plain(v) for k,v in value.items()}
        if isinstance(value,(list,tuple)):return [plain(x) for x in value]
        if isinstance(value,np.generic):return value.item()
        return value
    data=plain(model.input);data.pop('output',None)
    return hashlib.sha256(json.dumps(data,sort_keys=True).encode()).hexdigest()


def _capture(model):
    data, vectors, matrices, solvers, resources = {}, {}, {}, {}, {}
    matrix_names = {v.handle: k for k,v in vars(model).items() if isinstance(v,PETSc.Mat) and bool(v)}
    for name,value in vars(model).items():
        if isinstance(value,PETSc.Vec):
            vectors[name] = None if not value else dict(values=value.getArray().copy(),sizes=value.getSizes(),
                                                        kind=value.getType(),comm_size=value.getComm().size)
        elif isinstance(value,PETSc.Mat):
            if value:
                if value.getType() not in ('aij','seqaij','mpiaij'):
                    raise ValueError('Unsupported checkpoint matrix type: '+name+' '+value.getType())
                matrices[name]=dict(csr=tuple(x.copy() for x in value.getValuesCSR()),sizes=value.getSizes())
            else:matrices[name]=None
        elif isinstance(value,PETSc.KSP):
            if value:
                operators=value.getOperators()
                solvers[name]=dict(kind=value.getType(),tolerances=value.getTolerances(),pc=value.getPC().getType(),
                    prefix=value.getOptionsPrefix(),nonzero=value.getInitialGuessNonzero(),
                    operators=[matrix_names.get(x.handle) for x in operators])
            else:solvers[name]=None
        elif isinstance(value,PETSc.Object) or isinstance(value,types.MethodType) or name in ('profiler','vtkMesh','modelRunTime','outputDir') or name.startswith('_soluteBudget'):
            resources[name]=type(value).__name__
        else:
            # save() serializes the complete snapshot once and rejects unknown
            # nonserializable values before writing this rank's state file.
            data[name]=value
    return dict(data=data,vectors=vectors,matrices=matrices,solvers=solvers,resources=resources,
                owned_ids=np.asarray(model.glbIDs).copy(),local_ids=np.asarray(model.locIDs).copy(),
                mesh_sha=_sha(model.mCoords),configuration=_configuration(model),petsc=PETSc.Sys.getVersion())


def _generation(prefix,token):
    if len(token)!=32 or any(x not in '0123456789abcdef' for x in token):
        raise ValueError('Invalid coupled checkpoint generation')
    return Path(prefix).resolve()/token


def save(model,comm,prefix):
    token=comm.bcast(uuid.uuid4().hex if comm.rank==0 else None,root=0)
    folder=_generation(prefix,token)
    error=None;digest=None
    try:
        folder.mkdir(parents=True,exist_ok=True,mode=0o700)
        payload=pickle.dumps(_capture(model),protocol=5)
        digest=hashlib.sha256(payload).hexdigest()
        tmp=folder/('rank-%04d.tmp'%comm.rank)
        with tmp.open('wb') as stream:stream.write(payload);stream.flush();os.fsync(stream.fileno())
        tmp.rename(folder/('rank-%04d.pkl'%comm.rank))
    except Exception as exc:error=str(exc)
    errors=comm.allgather(error);digests=comm.gather(digest,root=0)
    if any(errors):raise RuntimeError('GoSPL checkpoint write failed: '+str(errors))
    error=None
    if comm.rank==0:
        try:
            shutil.copyfile(model.meshFile,folder/'mesh.npz')
            metadata=dict(format=1,mpi_size=comm.size,files=digests,mesh_sha=_sha(model.mCoords),time=float(model.tNow))
            tmp=folder/'manifest.tmp';tmp.write_text(json.dumps(metadata,indent=2)+'\n');tmp.rename(folder/'manifest.json')
        except Exception as exc:error=str(exc)
    error=comm.bcast(error,root=0)
    if error:raise RuntimeError(error)
    return token


def prepare(comm,prefix,token,destination):
    folder=_generation(prefix,token)
    error=None
    if comm.rank==0:
        try:
            meta=json.loads((folder/'manifest.json').read_text())
            if meta['format']!=1 or meta['mpi_size']!=comm.size:
                raise ValueError('GoSPL checkpoint requires the same MPI rank count')
            shutil.copyfile(folder/'mesh.npz',destination)
        except Exception as exc:error=str(exc)
    error=comm.bcast(error,root=0)
    if error:raise RuntimeError(error)
    return 0


def load(model,comm,prefix,token):
    folder=_generation(prefix,token);state=None;error=None
    try:
        meta=json.loads((folder/'manifest.json').read_text())
        if meta['format']!=1 or meta['mpi_size']!=comm.size:raise ValueError('Incompatible GoSPL checkpoint')
        payload=(folder/('rank-%04d.pkl'%comm.rank)).read_bytes()
        if hashlib.sha256(payload).hexdigest()!=meta['files'][comm.rank]:raise ValueError('GoSPL checkpoint checksum mismatch')
        state=pickle.loads(payload)
        if state['configuration']!=_configuration(model) or state['petsc']!=PETSc.Sys.getVersion():
            raise ValueError('GoSPL checkpoint configuration or PETSc version differs')
        if state['mesh_sha']!=_sha(model.mCoords) or not np.array_equal(state['owned_ids'],model.glbIDs) or not np.array_equal(state['local_ids'],model.locIDs):
            raise ValueError('GoSPL checkpoint mesh or MPI partition differs')
    except Exception as exc:error=str(exc)
    errors=comm.allgather(error)
    if any(errors):raise RuntimeError('GoSPL checkpoint validation failed: '+str(errors))
    for name,value in state['data'].items():setattr(model,name,value)
    for name,record in state['vectors'].items():
        current=getattr(model,name,None)
        if record is None:
            if current is not None and current:current.destroy()
            setattr(model,name,PETSc.Vec());continue
        if current is None or not current:
            if record['kind']=='seq':current=PETSc.Vec().createSeq(record['sizes'][1],comm=PETSc.COMM_SELF)
            elif record['kind']=='mpi':current=PETSc.Vec().createMPI(record['sizes'],comm=comm)
            else:raise ValueError('Unsupported checkpoint vector '+name)
            setattr(model,name,current)
        if current.getSizes()!=record['sizes']:raise ValueError('GoSPL vector size changed: '+name)
        current.getArray()[:]=record['values']
    for name,record in state['matrices'].items():
        current=getattr(model,name,None)
        if record is None:
            if current is not None and current:current.destroy()
            setattr(model,name,PETSc.Mat());continue
        if current is None or not current:
            current=PETSc.Mat().createAIJ(size=record['sizes'],csr=record['csr'],comm=comm)
            current.setLGMap(model.lgmap_row,model.lgmap_col)
            setattr(model,name,current)
        else:
            current.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR,False)
            current.zeroEntries();current.setValuesCSR(*record['csr']);current.assemble()
    for name,record in state['solvers'].items():
        current=getattr(model,name,None)
        if current is not None and current:current.destroy()
        if record is None:setattr(model,name,PETSc.KSP());continue
        current=PETSc.KSP().create(comm);current.setType(record['kind']);current.getPC().setType(record['pc'])
        if record['prefix'] is not None:current.setOptionsPrefix(record['prefix'])
        current.setTolerances(*record['tolerances'])
        current.setInitialGuessNonzero(record['nonzero']);current.setFromOptions()
        names=record['operators']
        if all(names):current.setOperators(*(getattr(model,key) for key in names))
        setattr(model,name,current)
    return 0
