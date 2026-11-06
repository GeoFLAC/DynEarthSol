#!/usr/bin/env python3
"""Exercise native MMG/ParMmg CLIs on DES snapshots, outside the timestep solver.

Checks geometry, boundary references, frozen cells, metrics, marker relocation,
and diagnostic field transfer. This does not integrate a new DES remeshing backend.
"""
import argparse
import hashlib
from itertools import combinations
import json
import math
import os
from pathlib import Path
import shlex
import signal
import subprocess
import tempfile
import time

import h5py
import numpy as np
from scipy.spatial import cKDTree


def boundary(cells):
    faces = np.concatenate([np.sort(cells[:, keep], axis=1)
                            for keep in combinations(range(cells.shape[1]), cells.shape[1]-1)])
    unique, count = np.unique(faces, axis=0, return_counts=True)
    if np.any(count > 2):
        raise ValueError('Non-manifold facet')
    return unique[count == 1]


def geometry(points, cells):
    vertices = points[cells]
    determinant = np.linalg.det(vertices[:, 1:] - vertices[:, :1])
    edges = np.stack([vertices[:, j] - vertices[:, i]
                      for i, j in combinations(range(cells.shape[1]), 2)], axis=1)
    edge2 = np.sum(edges * edges, axis=2)
    dim = points.shape[1]
    measure = np.abs(determinant) / math.factorial(dim)
    if dim == 2:
        quality = 4 * np.sqrt(3) * measure / edge2.sum(axis=1)
    else:
        quality = 12 * (3 * measure) ** (2/3) / edge2.sum(axis=1)
    return determinant, measure, quality


def read_snapshot(path):
    with h5py.File(path) as f:
        points, cells = f['coordinate'][()], f['connectivity'][()]
        facets = boundary(cells)
        refs = np.bitwise_and.reduce(f['bcflag'][()][facets], axis=1)
        if np.any(refs == 0):
            raise ValueError('Cannot infer boundary labels from source bcflag')
        nodal = {key: f[key][()] for key in ['velocity', 'temperature', 'pore pressure'] if key in f}
        elemental = {key: f[key][()] for key in ['stress', 'strain', 'plastic strain', 'material'] if key in f}
        markers = {key: f[key][()] for key in f if key.startswith('markerset.')
                   and isinstance(f[key], h5py.Dataset) and f[key].shape[0] == len(f['markerset.coord'])}
    return dict(points=points, cells=cells, facets=facets, refs=refs,
                nodal=nodal, elemental=elemental, markers=markers)


def write_mesh(path, state, frozen):
    points, cells, facets = (state[k] for k in ['points', 'cells', 'facets'])
    dim = points.shape[1]
    cells = cells.copy()
    negative = geometry(points, cells)[0] < 0
    cells[negative, :2] = cells[negative, 1::-1]
    with path.open('w') as f:
        f.write(f'MeshVersionFormatted 2\nDimension {dim}\nVertices\n{len(points)}\n')
        np.savetxt(f, np.column_stack([points, np.zeros(len(points), int)]), fmt=['%.17g']*dim+['%d'])
        f.write(f"{'Triangles' if dim == 2 else 'Tetrahedra'}\n{len(cells)}\n")
        # References preserve material regions represented by the original cells.
        refs = np.rint(state['elemental'].get('material', np.zeros(len(cells)))).astype(int)+1
        np.savetxt(f, np.column_stack([cells+1, refs]), fmt='%d')
        f.write(f"{'Edges' if dim == 2 else 'Triangles'}\n{len(facets)}\n")
        np.savetxt(f, np.column_stack([facets+1, state['refs']]), fmt='%d')
        if frozen:
            center = (points.min(axis=0) + points.max(axis=0))/2
            if 'frozen_points' in state:
                _, vertices = cKDTree(points).query(state['frozen_points'])
                cell = next(i for i, c in enumerate(cells) if set(c) == set(vertices))
            else:
                cell = int(np.argmin(np.linalg.norm(points[cells].mean(axis=1)-center, axis=1)))
            required = cells[cell]
            f.write(f'RequiredVertices\n{len(required)}\n')
            np.savetxt(f, required+1, fmt='%d')
            f.write(f"{'RequiredTriangles' if dim == 2 else 'RequiredTetrahedra'}\n1\n{cell+1}\n")
            state['frozen_points'] = points[required].copy()
        f.write('End\n')


def write_metric(path, points, kind, original):
    dim = points.shape[1]
    extent = np.ptp(original['points'], axis=0)
    _, measure, _ = geometry(original['points'], original['cells'])
    base = np.median(measure) ** (1/dim)
    x = (points[:, 0]-original['points'][:, 0].min())/extent[0]
    h = base*(1.3+0.6*x)
    if kind == 'scalar':
        values = h[:, None]
    else:
        # Medit tensors use lower-triangular order: xx,xy,yy,xz,yz,zz.
        values = np.zeros((len(points), dim*(dim+1)//2))
        diagonal = [0, 2] if dim == 2 else [0, 2, 5]
        values[:, diagonal] = 1/h[:, None]**2
        values[:, 0] *= 4
    with path.open('w') as f:
        f.write(f'MeshVersionFormatted 2\nDimension {dim}\nSolAtVertices\n{len(points)}\n1 {1 if kind == "scalar" else 3}\n')
        np.savetxt(f, values, fmt='%.17g')
        f.write('End\n')
    return base


def read_sections(path):
    tokens = ' '.join(line.split('#')[0] for line in path.read_text().splitlines()).split()
    result = {}; cursor = 0; dim = None
    widths = {'Vertices': None, 'Edges': 3, 'Triangles': 4, 'Tetrahedra': 5,
              'Corners': 1, 'RequiredVertices': 1, 'RequiredEdges': 1,
              'RequiredTriangles': 1, 'RequiredTetrahedra': 1, 'Ridges': 1,
              'Normals': 3, 'NormalAtVertices': 2, 'NormalAtTriangleVertices': 3,
              'Tangents': 3, 'TangentAtVertices': 2}
    while cursor < len(tokens):
        key = tokens[cursor]; cursor += 1
        if key == 'End': break
        if key in ('MeshVersionFormatted', 'Dimension'):
            value = int(tokens[cursor]); cursor += 1
            if key == 'Dimension': dim = value
            continue
        if key not in widths: raise ValueError(f'Unsupported Medit section: {key}')
        count = int(tokens[cursor]); cursor += 1
        width = dim+1 if key == 'Vertices' else widths[key]
        result[key] = np.array(tokens[cursor:cursor+count*width], dtype=float).reshape(count, width)
        cursor += count*width
    return dim, result


def read_mesh(path):
    dim, result = read_sections(path)
    points = result['Vertices'][:, :dim]
    cells = result['Triangles' if dim == 2 else 'Tetrahedra'][:, :dim+1].astype(int)-1
    facets = result['Edges' if dim == 2 else 'Triangles'][:, :dim].astype(int)-1
    refs = result['Edges' if dim == 2 else 'Triangles'][:, dim].astype(int)
    return points, cells, facets, refs


def check_metric(path, dim, expected_points):
    tokens = path.read_text().split()
    start = tokens.index('SolAtVertices')
    count, fields, kind = map(int, tokens[start+1:start+4])
    if count != expected_points or fields != 1 or kind not in (1, 3):
        raise ValueError('Invalid exported metric shape')
    width = 1 if kind == 1 else dim*(dim+1)//2
    values = np.array(tokens[start+4:start+4+count*width], float).reshape(count, width)
    if not np.isfinite(values).all(): raise ValueError('Non-finite exported metric')
    if kind == 1:
        eigenvalues = values
    else:
        matrices = np.zeros((count, dim, dim)); column = 0
        for i in range(dim):
            for j in range(i+1):
                matrices[:, i, j] = matrices[:, j, i] = values[:, column]; column += 1
        eigenvalues = np.linalg.eigvalsh(matrices)
    if np.any(eigenvalues <= 0): raise ValueError('Exported metric is not positive definite')
    return dict(metric_finite=True, metric_positive_definite=True)


class Locator:
    def __init__(self, points, cells):
        self.points, self.cells = points, cells
        self.vertices = points[cells]
        self.inverse = np.linalg.inv(np.swapaxes(self.vertices[:, 1:]-self.vertices[:, :1], 1, 2))
        self.tree = cKDTree(self.vertices.mean(axis=1))

    def locate(self, query):
        found = np.full(len(query), -1, int)
        weights = np.zeros((len(query), self.points.shape[1]+1))
        for k in sorted(set([min(n, len(self.cells)) for n in (16, 64, 256, len(self.cells))])):
            remaining = np.where(found < 0)[0]
            if not len(remaining): break
            # Bound temporary memory even if a fallback searches every cell.
            for chunk in np.array_split(remaining, max(1, math.ceil(len(remaining)*k/100000))):
                if not len(chunk): continue
                _, neighbors = self.tree.query(query[chunk], k=k)
                neighbors = np.asarray(neighbors).reshape(len(chunk), k)
                delta = query[chunk, None, :] - self.vertices[neighbors, 0]
                bary = np.einsum('...ij,...j->...i', self.inverse[neighbors], delta)
                bary = np.concatenate([1-bary.sum(axis=2, keepdims=True), bary], axis=2)
                valid = (bary.min(axis=2) >= -1e-10) & (bary.max(axis=2) <= 1+1e-10)
                row = np.where(valid.any(axis=1))[0]; choice = valid[row].argmax(axis=1)
                found[chunk[row]] = neighbors[row, choice]
                weights[chunk[row]] = bary[row, choice]
        if np.any(found < 0): raise ValueError(f'{int((found < 0).sum())} points outside source mesh')
        return found, weights


def transfer(old, points, cells):
    locator = Locator(old['points'], old['cells'])
    element, bary = locator.locate(points)
    nodal = {}
    for key, data in old['nodal'].items():
        nodal[key] = np.einsum('ij,ij...->i...', bary, data[old['cells'][element]])
    # Verify the locator/interpolation with an affine field independent of DES state.
    scale = max(np.ptp(old['points'], axis=0))
    affine = 1 + old['points'].sum(axis=1)/scale
    actual = np.einsum('ij,ij->i', bary, affine[old['cells'][element]])
    error = float(np.max(np.abs(actual-(1+points.sum(axis=1)/scale))))
    if error > 1e-11: raise ValueError('Affine interpolation failed')
    parent, _ = locator.locate(points[cells].mean(axis=1))
    elemental = {key: data[parent] for key, data in old['elemental'].items()}
    # This cell sampling is diagnostic and not conservative DES remapping.
    before = geometry(old['points'], old['cells'])[1]
    after = geometry(points, cells)[1]
    drift = {}
    for key, data in old['elemental'].items():
        if key == 'material': continue
        denominator = np.sum(before.reshape((-1,)+(1,)*(data.ndim-1))*np.abs(data), axis=0)
        change = np.sum(after.reshape((-1,)+(1,)*(data.ndim-1))*elemental[key], axis=0) - np.sum(before.reshape((-1,)+(1,)*(data.ndim-1))*data, axis=0)
        drift[key] = float(np.max(np.abs(change)/np.maximum(denominator, 1e-30)))
    markers = {key: data.copy() for key, data in old['markers'].items()}
    marker_error = 0.
    if markers:
        index, eta = Locator(points, cells).locate(markers['markerset.coord'])
        rebuilt = np.einsum('ij,ijk->ik', eta, points[cells[index]])
        marker_error = float(np.max(np.abs(rebuilt-markers['markerset.coord']))/scale)
        if marker_error > 1e-11: raise ValueError('Marker barycentric reconstruction failed')
        markers['markerset.elem'], markers['markerset.eta'] = index, eta
        for key in old['markers']:
            if key not in ('markerset.elem', 'markerset.eta') and not np.array_equal(markers[key], old['markers'][key]):
                raise ValueError('Marker identity/history changed')
    return nodal, elemental, markers, dict(affine_error=error, marker_reconstruction_error=marker_error,
                                         markers=len(markers.get('markerset.coord', [])), diagnostic_integral_drift=drift)


def run(command, directory, timeout):
    started = time.monotonic()
    with (directory/'run.log').open('w') as log:
        process = subprocess.Popen(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT,
                                   env=dict(os.environ, OMP_NUM_THREADS='1'), start_new_session=True)
        try: status = process.wait(timeout=timeout)
        except BaseException:
            os.killpg(process.pid, signal.SIGKILL); process.wait(); raise
    if status: raise RuntimeError(f'Exit {status}; see {directory}/run.log')
    return dict(command=command, wall_seconds=time.monotonic()-started)


def check_mesh(old, points, cells, facets, refs):
    if not np.isfinite(points).all() or cells.min() < 0 or cells.max() >= len(points):
        raise ValueError('Invalid coordinates/connectivity')
    determinant, volume, quality = geometry(points, cells)
    if np.any(determinant <= 0) or np.any(quality <= 0): raise ValueError('Inverted/degenerate element')
    actual = boundary(cells)
    supplied = {}
    for face, ref in zip(facets, refs):
        key = tuple(sorted(face)); supplied[key] = supplied.get(key, 0) | int(ref)
    labels = np.array([supplied.get(tuple(face), 0) for face in actual])
    if np.any(labels == 0) or set(labels) != set(old['refs']): raise ValueError('Boundary references lost')
    change = abs(volume.sum()/geometry(old['points'], old['cells'])[1].sum()-1)
    if change > 1e-6: raise ValueError(f'Domain area/volume changed by {change:.3g}')
    if 'frozen_points' in old:
        distance, nodes = cKDTree(points).query(old['frozen_points'])
        tolerance = 1e-10*max(np.ptp(points, axis=0))
        if distance.max() > tolerance: raise ValueError('A required vertex moved/disappeared')
        target = set(nodes.tolist())
        if not any(set(cell) == target for cell in cells): raise ValueError('A required element disappeared')
    return actual, labels, dict(nodes=len(points), elements=len(cells), min_quality=float(quality.min()),
                                volume_relative_change=float(change), boundary_references=sorted(set(map(int, labels))))


def restore_boundaries_2d(old, points, cells, facets, refs, tolerance):
    """Explicit prototype adapter: recover labels only within a geometric tolerance."""
    actual = boundary(cells)
    edges = old['points'][old['facets']]
    direction = edges[:, 1]-edges[:, 0]
    def distances(query, subset):
        delta = query[:, None, :]-edges[subset, 0]
        t = np.einsum('ijk,jk->ij', delta, direction[subset])/np.sum(direction[subset]**2, axis=1)
        projected = edges[subset, 0]+np.clip(t, 0, 1)[:, :, None]*direction[subset]
        return np.linalg.norm(query[:, None, :]-projected, axis=2)
    midpoint = points[actual].mean(axis=1)
    distance = distances(midpoint, np.arange(len(edges)))
    nearest = distance.argmin(axis=1); labels = old['refs'][nearest]
    worst = float(distance.min(axis=1).max())
    for label in set(labels):
        queries = points[actual[labels == label]].reshape(-1, 2)
        worst = max(worst, float(distances(queries, np.where(old['refs'] == label)[0]).min(axis=1).max()))
    if worst > tolerance: raise ValueError(f'Boundary reconstruction distance {worst} exceeds {tolerance}')
    native = {tuple(sorted(face)): int(label) for face, label in zip(facets, refs)}
    count = sum(native.get(tuple(face), 0) != label for face, label in zip(actual, labels))
    return actual, labels, dict(boundary_labels_reconstructed=True, changed_boundary_labels=int(count),
                                max_boundary_projection_distance=worst, projection_tolerance=tolerance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for dim in (2, 3):
        parser.add_argument(f'--input{dim}d', type=Path)
        parser.add_argument(f'--parmmg{dim}d', type=Path)
        parser.add_argument(f'--mmg{dim}d', type=Path)
    parser.add_argument('--ranks', type=int, nargs='+', default=[2, 3, 4])
    parser.add_argument('--metrics', nargs='+', choices=['scalar', 'tensor'], default=['scalar', 'tensor'])
    parser.add_argument('--constraints', nargs='+', choices=['free', 'frozen'], default=['free', 'frozen'])
    parser.add_argument('--cycles', type=int, default=2)
    parser.add_argument('--workdir', type=Path, required=True)
    parser.add_argument('--timeout', type=float, default=180)
    parser.add_argument('--launcher', default='mpiexec --bind-to none')
    parser.add_argument('--restore-boundaries2d', action='store_true',
                        help='Explicit geometric-label prototype adapter, not native ParMmg2D label preservation')
    args = parser.parse_args()
    if args.cycles < 1 or min(args.ranks) < 2: parser.error('Need cycles >= 1 and ranks >= 2')
    args.workdir.mkdir(parents=True, exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix='prototype-', dir=args.workdir))
    report = dict(completed=False, scope='External CLI prototype; diagnostic, nonconservative field transfer', cases=[])
    (output/'results.json').write_text(json.dumps(report, indent=2)+'\n')
    print('Artifacts:', output, flush=True)
    for dim in (2, 3):
        source = getattr(args, f'input{dim}d'); parallel = getattr(args, f'parmmg{dim}d'); serial = getattr(args, f'mmg{dim}d')
        if source is None: continue
        if parallel is None: parser.error(f'Need --parmmg{dim}d with an input')
        original = read_snapshot(source)
        if original['points'].shape[1] != dim: raise ValueError('Input dimension mismatch')
        variants = [('parmmg', parallel.resolve(), r) for r in args.ranks]
        if serial: variants.insert(0, ('mmg', serial.resolve(), 1))
        for metric in args.metrics:
            for constraint in args.constraints:
                for backend, executable, ranks in variants:
                    label = f'{dim}d-{metric}-{constraint}-{backend}{ranks}'
                    record = dict(case=label, input=str(source.resolve()), passed=False, cycles=[],
                                  executable=str(executable), executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest())
                    report['cases'].append(record)
                    state = read_snapshot(source)
                    try:
                        for cycle in range(args.cycles):
                            directory = output/label/f'cycle{cycle}'; directory.mkdir(parents=True)
                            write_mesh(directory/'input.mesh', state, constraint == 'frozen')
                            base = write_metric(directory/'input.sol', state['points'], metric, original)
                            command = [str(executable), '-in', 'input.mesh', '-out', 'output.mesh', '-sol', 'input.sol',
                                       '-v', '3', '-hausd', str(base*.01)]
                            if backend == 'parmmg': command = shlex.split(args.launcher)+['-n', str(ranks)]+command+['-centralized-output', '-niter', '3']
                            timing = run(command, directory, args.timeout)
                            points, cells, facets, refs = read_mesh(directory/'output.mesh')
                            metric_checks = check_metric(directory/'output.sol', dim, len(points))
                            restoration = {}
                            if dim == 2 and backend == 'parmmg' and args.restore_boundaries2d:
                                facets, refs, restoration = restore_boundaries_2d(state, points, cells, facets, refs, base*.011)
                            facets, refs, checks = check_mesh(state, points, cells, facets, refs)
                            nodal, elemental, markers, transfer_checks = transfer(state, points, cells)
                            checks.update(restoration); checks.update(metric_checks); checks.update(transfer_checks); checks.update(timing)
                            np.savez(directory/'transferred.npz', coordinate=points, connectivity=cells,
                                     **nodal, **elemental, **markers)
                            pinned = state.get('frozen_points')
                            state = dict(points=points, cells=cells, facets=facets, refs=refs,
                                         nodal=nodal, elemental=elemental, markers=markers)
                            if pinned is not None: state['frozen_points'] = pinned
                            record['cycles'].append(checks)
                        record['passed'] = True
                    except Exception as error:
                        record['error'] = str(error)
                    (output/'results.json').write_text(json.dumps(report, indent=2)+'\n')
                    print(label, 'PASS' if record['passed'] else record.get('error'), flush=True)
    report['completed'] = True
    report['all_passed'] = bool(report['cases']) and all(row['passed'] for row in report['cases'])
    (output/'results.json').write_text(json.dumps(report, indent=2)+'\n')
    print('Results:', output/'results.json', flush=True)
    return 0 if report['all_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
