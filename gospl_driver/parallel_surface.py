"""MPI-aware boundary transfers for gospl_extensions' distributed PETSc vectors.

Only coupling interpolation/indexing lives here. GoSPL's native time stepping,
erosion, diffusion and finite-volume advection remain the installed algorithms.
"""
import types
import numpy as np
from mpi4py import MPI
from scipy.spatial import cKDTree


def _weights(source, query, k, power, tree=None):
    k = max(1, min(int(k), len(source)))
    distances, indices = (tree if tree is not None else cKDTree(source)).query(query, k=k)
    if k == 1:
        distances, indices = distances[:, None], indices[:, None]
    weights = 1.0 / np.maximum(distances, 1e-20) ** power
    coincident = distances[:, 0] < 1e-20
    weights[coincident] = 0
    weights[coincident, 0] = 1
    weights /= weights.sum(axis=1, keepdims=True)
    return indices, weights


class SurfaceTransfer:
    def __init__(self, model, comm):
        self.comm = comm
        self.ids = np.asarray(model.glbIDs, dtype=np.int64)
        if len(self.ids) != len(model.hGlobal.getArray()):
            raise ValueError('GoSPL owned coordinate/vector sizes differ')
        ids = comm.gather(self.ids, root=0)
        self.counts = self.offsets = self.all_ids = self.values = self.global_values = None
        error = None
        if comm.rank == 0:
            self.tree = cKDTree(model.mCoords)
            self.counts = np.array([len(x) for x in ids], dtype=np.int32)
            self.offsets = np.r_[0, np.cumsum(self.counts[:-1])].astype(np.int32)
            self.all_ids = np.concatenate(ids)
            if not np.array_equal(np.sort(self.all_ids), np.arange(len(model.mCoords))):
                error = 'GoSPL owned vertex IDs do not form a unique global partition'
            self.values = np.empty(len(self.all_ids), dtype=np.float64)
            self.global_values = np.empty(len(model.mCoords), dtype=np.float64)
        error = comm.bcast(error, root=0)
        if error:
            raise ValueError(error)

    def gather(self, local):
        receive = [self.values, self.counts, self.offsets, MPI.DOUBLE] if self.comm.rank == 0 else None
        self.comm.Gatherv(np.ascontiguousarray(local, dtype=np.float64), receive, root=0)
        if self.comm.rank == 0:
            self.global_values[self.all_ids] = self.values
            return self.global_values
        return None

    def interpolate(self, model, local, points, k, power):
        values = self.gather(local)
        result = None
        if self.comm.rank == 0:
            ids, weights = _weights(model.mCoords, points, k, power, self.tree)
            result = np.sum(weights * values[ids], axis=1)
        return self.comm.bcast(result, root=0)


def bind_parallel_surface(model, comm):
    transfer = SurfaceTransfer(model, comm)

    def seed(self, elevdata, k=3, power=1.0):
        points = np.asarray(elevdata['coords'], dtype=np.float64)
        elevations = np.asarray(elevdata.get('elev', elevdata.get('elevations')), dtype=np.float64)
        ids, weights = _weights(points, self.mCoords, k, power)
        mapped = np.sum(weights * elevations[ids], axis=1)
        self.hGlobal.getArray()[:] = mapped[self.glbIDs]
        self.hLocal.getArray()[:] = mapped[self.locIDs]

    def interpolate(self, points, k=3, power=1.0):
        return transfer.interpolate(self, self.hGlobal.getArray(), points, k, power)

    def step(self, dt, points, k=3, power=1.0):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError('Coupling interval must be finite and positive')
        has_velocity = getattr(self, '_upsub_override', None) is not None
        old_upsub, old_tecdata = getattr(self, 'upsub', None), self.tecdata
        try:
            if has_velocity and getattr(self, '_vx_override', None) is not None:
                if self.advscheme > 0:
                    from gospl._fortran import getfacevelocity
                    velocity = np.zeros((self.lpoints, 3), dtype=np.float64)
                    velocity[:, 0] = self._vx_override[self.locIDs]
                    velocity[:, 1] = self._vy_override[self.locIDs]
                    self.hdisp = velocity
                    getfacevelocity(self.lpoints, velocity)
                    old_dt = self.dt
                    try:
                        self.dt = dt
                        self._varAdvector()
                    finally:
                        self.dt = old_dt
                        self.hdisp = None
                else:
                    values = transfer.gather(self.hGlobal.getArray())
                    mapped = None
                    if comm.rank == 0:
                        displaced = self.mCoords.copy()
                        displaced[:, 0] -= self._vx_override * dt
                        displaced[:, 1] -= self._vy_override * dt
                        ids, weights = _weights(self.mCoords, displaced, k, power, transfer.tree)
                        mapped = np.sum(weights * values[ids], axis=1)
                    mapped = comm.bcast(mapped, root=0)
                    self.hGlobal.getArray()[:] = mapped[self.glbIDs]
                    self.hLocal.getArray()[:] = mapped[self.locIDs]
                self._vx_override = self._vy_override = None
            if has_velocity:
                self.upsub = self._upsub_override[self.locIDs]
                self.tecdata = None
            before = self.hGlobal.getArray().copy()
            self.runProcessesForDt(dt, verbose=False, skip_tectonics=False)
            delta = self.hGlobal.getArray().copy() - before
            if has_velocity:
                delta -= self._upsub_override[self.glbIDs] * dt
            return transfer.interpolate(self, delta, np.asarray(points), k, power)
        finally:
            if has_velocity:
                self.tecdata, self.upsub = old_tecdata, old_upsub
                self._upsub_override = None

    model.apply_elevation_data = types.MethodType(seed, model)
    model.interpolate_elevation_to_points = types.MethodType(interpolate, model)
    model.run_and_get_erosion = types.MethodType(step, model)
