#!/usr/bin/env python3
"""Refine the thin-slab pumping comparison of pumping_slab.py on structured meshes.

pumping_slab.py meshes the slab with unstructured tetrahedra, whose size hardly
follows mesh.resolution here (1698 and 1839 nodes at 20 m and 10 m). This
driver keeps its configuration, comparison and figure, and only selects the
structured element shape (mesh.meshing_elem_shape = 1), which does refine.
The remaining difference from the averaged (generalised Bear-Corapcioglu)
model is then the slab-model error O((B/r)^2), not the mesh.

    slab_refinement.py EXECUTABLE OUTPUT [--resolutions 20 10 5]
"""
import argparse
import json
import os
from pathlib import Path

import pumping_slab as ps


import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tests/functional/hydromechanics'))
from des_case import digest, run_solver


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--resolutions', type=float, nargs='+', default=[20, 10, 5])
    a = p.parse_args()
    exe, out = a.executable.resolve(), a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(Path(__file__)),
                  pumping_slab_sha256=digest(Path(ps.__file__)), cases=[])
    for h in a.resolutions:
        folder = out/f'res{h:g}'; case = folder/'slab-3d'; case.mkdir(parents=True)
        cfg = ps.config(h)
        cfg['mesh']['meshing_elem_shape'] = '1'
        with (case/'input.cfg').open('w') as f: cfg.write(f)
        run_solver([str(exe), 'input.cfg'], log=case / 'run.log', cwd=case, check=True, env=dict(os.environ, OMP_NUM_THREADS='1'), timeout=14400)
        aq, rows = ps.compare(case)
        ps.figure(aq, rows, folder)
        model = ps.Dynearthsol(str(case/'slab'))
        summary = [{k: v for k, v in row.items() if k != 'plot'} for row in rows]
        item = dict(resolution_m=h, nodes=len(model.read_field(0, 'coordinate')),
                    elements=len(model.read_field(0, 'connectivity')), rows=summary,
                    pressure_rms=max(r['pressure_rms'] for r in summary),
                    areal_rms=max(r['areal_rms'] for r in summary),
                    vertical_rms=max(r['vertical_rms'] for r in summary),
                    input_sha256=digest(case/'input.cfg'))
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False))
        print(json.dumps({k: v for k, v in item.items() if k != 'rows'}), flush=True)
    assert digest(exe) == result['executable_sha256']


if __name__ == '__main__':
    main()
