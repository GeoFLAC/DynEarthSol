#!/usr/bin/env python3
"""Compare tetrahedral results at common physical points, independent of IDs.

Reports geometry validity, volume integrals, interior samples and top elevation.
Errors are measurements, not a universal scientific acceptance tolerance.
"""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.spatial import cKDTree


def error(a, b):
    scale = float(np.max(abs(a))) or 1.
    rms = float(np.sqrt(np.mean(a*a))) or 1.
    return dict(absolute_max=float(np.max(abs(a-b))),
                relative_max=float(np.max(abs(a-b)))/scale,
                relative_rms=float(np.sqrt(np.mean((a-b)**2)))/rms)


class Mesh:
    def __init__(self, path):
        self.path = path
        with h5py.File(path) as file:
            grid = file['VTKHDF/grid']
            self.points = grid['Points'][()]
            self.cells = grid['Connectivity'][()].reshape(-1, 4)
            self.nodal = {key: value[()] for key, value in grid['PointData'].items()}
            self.fields = {key: value[()] for key, value in grid['CellData'].items()}
            self.time = float(grid['FieldData/time_yr'][()].item())
        self.vertices = self.points[self.cells]
        edge = self.vertices[:, 1:] - self.vertices[:, :1]
        det = np.linalg.det(edge)
        self.volume = abs(det)/6
        self.inverse = np.linalg.inv(edge.transpose(0, 2, 1))
        self.tree = cKDTree(self.vertices.mean(axis=1))
        self.valid = bool(np.isfinite(self.points).all() and np.isfinite(det).all()
                          and (abs(det) > 0).all() and ((det > 0).all() or (det < 0).all())
                          and all(np.isfinite(v).all() for v in self.fields.values())
                          and all(np.isfinite(v).all() for v in self.nodal.values()))

    def locate(self, points):
        found = np.full(len(points), -1, dtype=int)
        weights = np.zeros((len(points), 4))
        for count in (32, 128, 512):
            remaining = np.flatnonzero(found < 0)
            if not len(remaining):
                break
            count = min(count, len(self.cells))
            _, candidates = self.tree.query(points[remaining], k=count)
            candidates = np.asarray(candidates).reshape(len(remaining), count)
            delta = points[remaining, None, :] - self.vertices[candidates, 0]
            eta = np.einsum('nkij,nkj->nki', self.inverse[candidates], delta)
            bary = np.concatenate((1-eta.sum(axis=2, keepdims=True), eta), axis=2)
            inside = (bary >= -1e-10).all(axis=2) & (bary <= 1+1e-10).all(axis=2)
            rows = np.flatnonzero(inside.any(axis=1))
            columns = inside[rows].argmax(axis=1)
            found[remaining[rows]] = candidates[rows, columns]
            weights[remaining[rows]] = bary[rows, columns]
        return found, weights

    def sample(self, name, cell, weights):
        if name in self.nodal:
            values = self.nodal[name][self.cells[cell]]
            return np.einsum('ni,ni...->n...', weights, values)
        return self.fields[name][cell]

    def integrals(self):
        result = {}
        for name in ('stress', 'strain', 'strain-rate', 'plastic strain', 'density', 'temperature'):
            if name in self.nodal:
                nodal = self.nodal[name][self.cells].reshape(len(self.cells), 4, -1)
                values = nodal.mean(axis=1)
                # Exact integral of the square of a linear tetrahedral field.
                squares = (nodal.sum(axis=1)**2+(nodal**2).sum(axis=1))/20
            else:
                values = self.fields[name].reshape(len(self.cells), -1)
                squares = values**2
            result[name] = dict(mean=(values*self.volume[:, None]).sum(axis=0)/self.volume.sum(),
                                rms=np.sqrt((squares*self.volume[:, None]).sum(axis=0)/self.volume.sum()))
        return result

    def surface(self, xy):
        # Every facet with three top-boundary vertices belongs to the surface.
        top = (self.nodal['bcflag'].astype(int) & 32) != 0
        facets = np.concatenate([self.cells[:, indices] for indices in
                                 ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))])
        facets = facets[top[facets].all(axis=1)]
        triangles = self.points[facets]
        if not len(triangles):
            raise ValueError('No top facets')
        tree = cKDTree(triangles[:, :, :2].mean(axis=1))
        _, ids = tree.query(xy, k=min(32, len(triangles)))
        ids = np.asarray(ids).reshape(len(xy), -1)
        vertices = triangles[ids]
        edge = (vertices[:, :, 1:, :2] - vertices[:, :, :1, :2]).swapaxes(2, 3)
        eta = np.linalg.solve(edge, (xy[:, None, :] - vertices[:, :, 0, :2])[..., None])[..., 0]
        bary = np.concatenate((1-eta.sum(axis=2, keepdims=True), eta), axis=2)
        inside = (bary >= -1e-10).all(axis=2)
        covered = inside.any(axis=1)
        rows = np.flatnonzero(covered)
        col = inside[rows].argmax(axis=1)
        z = np.full(len(xy), np.nan)
        z[rows] = (bary[rows, col]*vertices[rows, col, :, 2]).sum(axis=1)
        return z, covered


def compare(first, second, samples=2048):
    a, b = Mesh(first), Mesh(second)
    if not a.valid or not b.valid:
        raise ValueError('Invalid mesh or non-finite physical field')
    if not np.isclose(a.time, b.time, rtol=1e-12, atol=1e-12):
        raise ValueError('Physical times differ')
    low = np.maximum(a.points.min(axis=0), b.points.min(axis=0))
    high = np.minimum(a.points.max(axis=0), b.points.max(axis=0))
    rng = np.random.default_rng(117)
    points = low + (0.05+0.9*rng.random((samples, 3)))*(high-low)
    ai, aw = a.locate(points)
    bi, bw = b.locate(points)
    both = (ai >= 0) & (bi >= 0)
    if not both.all():
        raise ValueError('Interior sample coverage incomplete: %d/%d' % (both.sum(), samples))
    common = {name: error(a.sample(name, ai, aw), b.sample(name, bi, bw))
              for name in ('temperature', 'velocity', 'stress', 'strain', 'strain-rate', 'plastic strain')}
    x, y = np.meshgrid(np.linspace(low[0]+.05*(high[0]-low[0]), high[0]-.05*(high[0]-low[0]), 23),
                       np.linspace(low[1]+.05*(high[1]-low[1]), high[1]-.05*(high[1]-low[1]), 23))
    xy = np.column_stack((x.ravel(), y.ravel()))
    az, ac = a.surface(xy)
    bz, bc = b.surface(xy)
    if not (ac & bc).all():
        raise ValueError('Top sample coverage incomplete')
    integrals = {}
    other = b.integrals()
    for name, record in a.integrals().items():
        integrals[name] = {kind: error(record[kind], other[name][kind]) for kind in record}
    return dict(geometry_valid=True, time_yr=a.time, common_interior_points=samples,
                common_top_points=len(xy), volume=dict(reference=float(a.volume.sum()),
                candidate=float(b.volume.sum()), relative_difference=float(abs(a.volume.sum()-b.volume.sum())/a.volume.sum())),
                volume_weighted=integrals, interior=common,
                top_elevation=dict(absolute_max_m=float(abs(az-bz).max()), rms_m=float(np.sqrt(np.mean((az-bz)**2)))),
                inputs=[dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                             nodes=len(mesh.points), elements=len(mesh.cells)) for path, mesh in ((first, a), (second, b))])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=Path)
    parser.add_argument('candidate', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = compare(args.reference, args.candidate)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print('Geometry and common-point coverage PASS; measurements:', args.output)
