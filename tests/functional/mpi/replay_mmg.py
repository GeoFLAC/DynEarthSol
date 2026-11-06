#!/usr/bin/env python3
"""Replay exact captured DES/MMG inputs against native serial/parallel remeshers.

No new metric or freezing policy is constructed. Full DES field transfer remains
in the existing solver; this replay checks the external mesh/constraint contract.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import tempfile

import numpy as np
from scipy.spatial import cKDTree
from remesh_prototype import boundary, check_mesh, check_metric, read_mesh, read_sections, run


def constraints(capture, points, cells, facets):
    dim, sections = read_sections(capture/'input.mesh')
    original = sections['Vertices'][:, :dim]
    distance, mapping = cKDTree(points).query(original)
    tolerance = 1e-10*max(np.ptp(original, axis=0))
    pinned = sections.get('RequiredVertices', np.empty((0, 1))).astype(int).ravel()-1
    original_cells = sections['Triangles' if dim == 2 else 'Tetrahedra'][:, :dim+1].astype(int)-1
    required = sections.get('RequiredTriangles' if dim == 2 else 'RequiredTetrahedra', np.empty((0, 1))).astype(int).ravel()-1
    actual = {tuple(sorted(cell)) for cell in cells}
    lost = sum(np.any(distance[original_cells[e]] > tolerance) or
               tuple(sorted(mapping[original_cells[e]])) not in actual for e in required)
    original_facets = sections['Edges' if dim == 2 else 'Triangles'][:, :dim].astype(int)-1
    required_facets = sections.get('RequiredEdges' if dim == 2 else 'RequiredTriangles', np.empty((0, 1))).astype(int).ravel()-1
    actual_facets = {tuple(sorted(face)) for face in facets}
    lost_facets = sum(np.any(distance[original_facets[e]] > tolerance) or
                     tuple(sorted(mapping[original_facets[e]])) not in actual_facets for e in required_facets)
    return dict(required_nodes=len(pinned), required_elements=len(required), required_facets=len(required_facets),
                moved_or_lost_required_nodes=int(np.sum(distance[pinned] > tolerance)),
                lost_required_elements=int(lost), lost_required_facets=int(lost_facets))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--captures', type=Path, nargs='+', required=True)
    for dim in (2, 3):
        p.add_argument(f'--parmmg{dim}d', type=Path)
        p.add_argument(f'--mmg{dim}d', type=Path)
    p.add_argument('--ranks', type=int, nargs='+', default=[2, 3, 4])
    p.add_argument('--workdir', type=Path, required=True)
    p.add_argument('--launcher', default='mpiexec --bind-to none')
    p.add_argument('--timeout', type=float, default=180)
    a = p.parse_args()
    a.workdir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='replay-', dir=a.workdir)); print('Artifacts:', work, flush=True)
    result = dict(completed=False, cases=[], scope='Captured DES MMG boundary; solver policy unchanged')
    for number, capture in enumerate(a.captures):
        dim, _ = read_sections(capture/'input.mesh')
        old_points, old_cells, old_faces, old_refs = read_mesh(capture/'input.mesh')
        refs_by_face = {tuple(sorted(face)): int(ref) for face, ref in zip(old_faces, old_refs)}
        old_boundary = boundary(old_cells)
        old = dict(points=old_points, cells=old_cells, facets=old_boundary,
                   refs=np.array([refs_by_face[tuple(face)] for face in old_boundary]))
        variants = [('mmg', getattr(a, f'mmg{dim}d'), 1)] + [('parmmg', getattr(a, f'parmmg{dim}d'), r) for r in a.ranks]
        for backend, executable, ranks in variants:
            if executable is None: p.error(f'Missing executable for {backend} {dim}d')
            directory = work/f'input{number}-{dim}d-{backend}{ranks}'; directory.mkdir()
            for name in ['input.mesh', 'input.sol']: shutil.copy2(capture/name, directory/name)
            row = dict(capture=str(capture.resolve()), ndims=dim, backend=backend, ranks=ranks,
                       passed=False, executable=str(executable.resolve()),
                       executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest())
            result['cases'].append(row)
            try:
                command = [str(executable.resolve()), '-in', 'input.mesh', '-sol', 'input.sol', '-out', 'output.mesh', '-v', '3']
                command += shlex.split((capture/'parameters.txt').read_text())
                if backend == 'parmmg': command = shlex.split(a.launcher)+['-n', str(ranks)]+command+['-centralized-output', '-niter', '3']
                row.update(run(command, directory, a.timeout))
                points, cells, facets, refs = read_mesh(directory/'output.mesh')
                row.update(constraints(capture, points, cells, facets))
                row.update(check_metric(directory/'output.sol', dim, len(points)))
                _, _, checks = check_mesh(old, points, cells, facets, refs); row.update(checks)
                if any(row[k] for k in ['moved_or_lost_required_nodes', 'lost_required_elements', 'lost_required_facets']):
                    raise ValueError('Required entities were changed')
                row['passed'] = True
            except Exception as error: row['error'] = str(error)
            (work/'results.json').write_text(json.dumps(result, indent=2)+'\n')
            print(directory.name, 'PASS' if row['passed'] else row.get('error'), flush=True)
    result['completed'] = True; result['all_passed'] = all(row['passed'] for row in result['cases'])
    (work/'results.json').write_text(json.dumps(result, indent=2)+'\n')
    print('Results:', work/'results.json', flush=True)
    return 0 if result['all_passed'] else 1


if __name__ == '__main__': raise SystemExit(main())
