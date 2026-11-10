"""Replay the effective-stress material check and draw its selected path."""

import argparse
import json
from pathlib import Path
import shutil
import sys

import numpy as np
import prepare_figures as pf

import diagram_style as ds
from setup_drawings import schematic

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tests/functional/hydromechanics'))
import check_pore_strength as strength


def build(study, output, record_source=None):
    integrated = record_source is not None
    record_source = record_source or (lambda path: None)
    output.mkdir(parents=True, exist_ok=integrated)
    target = output/'figures' if integrated else output
    target.mkdir(exist_ok=True)
    metrics = {}
    for dim in (2, 3):
        group = study/f'v3-candidate-{dim}d'
        report = json.loads((group/'qualification.json').read_text())
        assert report['passed']
        record_source(group/'qualification.json')
        for rheology in ('elasto-plastic', 'elasto-plastic-rate-state-friction'):
            for name, confinement, pressure, biot, retention in strength.CASES:
                label = f'{rheology}-{name}'
                replay = strength.compare(group/label, dim, confinement, pressure, biot, retention)
                assert replay == report['cases'][label]
                record_source(group/label/'input.cfg')
                for path in (group/label).glob('cell.*'):
                    if '.save.' in path.name or path.suffix == '.info': record_source(path)
                metrics[f'{dim}d-{label}'] = replay

    if not integrated: schematic('strength', output)

    folder = study/'v3-candidate-3d/elasto-plastic-pressurized'
    model = strength.cell.Dynearthsol(str(folder/'cell'))
    frames = model.frames
    strain = np.array([model.read_field(f, 'strain')[:, 0].mean() for f in frames])
    load = np.array([(model.read_field(f, 'stress')[:, 0]-model.read_field(f, 'stress')[:, 2]).mean()/2 for f in frames])
    limit = metrics['3d-elasto-plastic-pressurized']['expected_yield_load_Pa']
    grid = np.linspace(0, strength.STRAIN, 301)
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(1); ax = axes[0]
        ds.analytical(ax, grid*1e3, np.minimum(2*strength.G*grid, limit)/1e3)
        ds.des(ax, strain*1e3, load/1e3)
        ds.panel(ax, '', 'Effective-stress yield at constant pore pressure',
                 xlabel='imposed normal strain (mstrain)', ylabel='principal shear load (kPa)')
        pf.save_figure(fig, target, 'strength' if integrated else 'result')
    data = output/'data/strength' if integrated else output
    data.mkdir(parents=True, exist_ok=True)
    if integrated:
        shutil.copy2(folder/'input.cfg', data/'input.cfg')
        shutil.copy2(study/'qualification.json', data/'qualification.json')
        record_source(study/'qualification.json')
        for source in (Path(strength.__file__), Path(strength.cell.__file__), Path(__file__).with_name('unsaturated_reference.py')):
            shutil.copy2(source, data/source.name)
    (data/'replay.json').write_text(json.dumps(dict(passed=True, cases=metrics), indent=2)+'\n')

    return dict(key='strength', theme='mechanics', title='Pore pressure and material strength',
        problem='An affine material cell is extended in one direction and compressed equally in another, at fixed pore pressure. Effective stress must control its shear yield load.',
        result=f"The pressurized 3D case yields at {limit/1e3:.3f} kPa, matching the analytical load. The largest stress/plastic-strain error over 20 load paths (EP, and EP-RSF at constant friction) is {max(c['max_error_over_yield_load'] for c in metrics.values()):.2e} of the yield load.",
        scope='Unit cells in 2D plane strain and 3D; five prescribed pressure/Biot/suction states and 64 load steps. The figure selects 3D EP at total confinement 1 MPa, pore pressure 0.8 MPa and Biot coefficient 1.',
        meaning='This isolates pressure-dependent skeleton strength and its plastic strain. Constant-friction EP-RSF is checked; evolving RSF is checked on the rate-and-state page; flow, dilation, tensile returns and slope/injection applications remain separate requirements.',
        dimension=3, spacing_m=1., metric=max(c['max_error_over_yield_load'] for c in metrics.values()),
        metric_label='stress/plastic-strain error / yield load', evidence=json.dumps(metrics))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', type=Path); parser.add_argument('output', type=Path)
    args = parser.parse_args(); build(args.study.resolve(), args.output.resolve())
