#!/usr/bin/env python3
"""Check production simplex moments against Qhull, symmetry and exact limits."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

import numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull, HalfspaceIntersection

ROOT = Path(__file__).resolve().parents[3]


from des_case import run_solver

def calculate(exe, pairs):
    data = ''.join(' '.join(format(v, '.17g') for a in pair for v in a.ravel())+'\n'
                   for pair in pairs)
    run = subprocess.run([str(exe)], input=data, capture_output=True, text=True,
                         check=True, timeout=30)
    values = np.array([[float(v) for v in line.split()] for line in run.stdout.splitlines()])
    n = pairs[0][0].shape[0]
    return values[:, 0], values[:, 1:].reshape(-1, n, n)


def reference_volume(a, b):
    dimension = a.shape[1]
    planes = np.r_[ConvexHull(a).equations, ConvexHull(b).equations]
    interior = linprog([*([0.]*dimension), -1.],
        A_ub=np.c_[planes[:, :dimension], np.ones(len(planes))], b_ub=-planes[:, -1],
        bounds=[(None, None)]*dimension+[(0, None)], method='highs')
    if not interior.success or interior.x[-1] < 1e-10:
        return 0.
    vertices = HalfspaceIntersection(planes, interior.x[:dimension]).intersections
    return ConvexHull(vertices).volume


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(20261109)
    records = []
    for dimension in (2, 3):
        exe = output/f'overlap{dimension}d'
        command = ['g++', '-std=c++17', '-O2', '-fopenmp', '-ffunction-sections',
                   '-fdata-sections', '-Wno-unknown-pragmas', '-DGPP1X', '-DSOA',
                   '-I'+str(ROOT), '-I'+str(ROOT/'nanoflann/include')]
        if dimension == 3:
            command.append('-DTHREED')
        command.extend([str(ROOT/'tests/functional/hydromechanics/overlap_driver.cxx'),
                        str(ROOT/'geometry.cxx'), str(ROOT/'utils.cxx'),
                        '-Wl,--gc-sections', '-o', str(exe)])
        run_solver(command, log=output / f'build-{dimension}d.log', check=True, timeout=120)
        unit = np.r_[np.zeros((1, dimension)), np.eye(dimension)]
        pairs = [(unit, unit), (unit, unit+2), (unit, unit*.5),
                 (unit, unit+np.eye(dimension)[0])]
        expected = [1/math.factorial(dimension), 0,
                    1/math.factorial(dimension)/2**dimension, 0]
        for _ in range(120):
            a, b = rng.normal(size=(2, dimension+1, dimension))
            b += .2
            if min(abs(np.linalg.det(x[1:]-x[0])) for x in (a, b)) < .01:
                continue
            pairs.append((a, b))
            expected.append(reference_volume(a, b))
        if dimension == 3:
            # Shared edges and nearly coincident faces from an actual TetGen remesh.
            a = np.array([[.5000000000022637, 1., -10.],
                          [.4999999999998543, 1., -9.489696420935509],
                          [0., .5000000000007558, -10.],
                          [.4999999999998609, 0., -9.489696420894388]])
            b = np.array([[.5000000000022637, 1., -10.],
                          [0., .5000000000007558, -10.],
                          [.5000000000001494, 0., -10.],
                          [.4999999999999998, .499999999999999, -9.489723737369077]])
            pairs.append((a, b))
            expected.append(reference_volume(a, b))
            a = np.array([[0., .5000000252446245, -7.815959578308086],
                          [0., .49999998222831776, -7.191948591508043],
                          [.5000000534656855, 1., -7.191948802780597],
                          [0., 1., -7.503954127483875]])
            b = np.array([[0., .49999998222831776, -7.191948591508043],
                          [.50000000966533, 1., -7.815959615957586],
                          [0., .5000000252446245, -7.815959578308086],
                          [0., 1., -7.503954127483875]])
            pairs.append((a, b))
            expected.append(reference_volume(a, b))
        volume, product = calculate(exe, pairs)
        reverse, transpose = calculate(exe, [(b, a) for a, b in pairs])
        np.testing.assert_allclose(volume, expected, rtol=2e-10, atol=2e-12)
        np.testing.assert_allclose(volume, reverse, rtol=2e-10, atol=2e-12)
        np.testing.assert_allclose(product.sum(axis=(1, 2)), volume, rtol=1e-12, atol=1e-14)
        np.testing.assert_allclose(product, transpose.swapaxes(1, 2), rtol=1e-9, atol=2e-12)
        assert product.min() > -1e-12
        n = dimension+1
        np.testing.assert_allclose(product[0], (np.ones((n, n))+np.eye(n))/
                                   (n*(n+1)*math.factorial(dimension)), rtol=1e-12, atol=1e-14)
        for shift in (np.zeros(dimension), np.full(dimension, 1e6)):
            for scale in (1e-3, 1., 1e3):
                v, _ = calculate(exe, [(unit*scale+shift, unit*.5*scale+shift)])
                np.testing.assert_allclose(v, scale**dimension/2**dimension/math.factorial(dimension),
                                          rtol=2e-6 if shift.any() else 1e-12)
        records.append(dict(dimension=dimension, pairs=len(pairs),
                            max_volume_error=float(np.max(np.abs(volume-expected))),
                            executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest()))
        print('PASS production overlap', dimension, flush=True)
    (output/'summary.json').write_text(json.dumps(dict(passed=True, cases=records), indent=2)+'\n')


if __name__ == '__main__':
    main()
