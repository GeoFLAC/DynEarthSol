"""Replay the elastic opening qualification and draw its selected mesh."""

import argparse
import json
from pathlib import Path
import shutil
import sys

import numpy as np

import diagram_style as ds
import prepare_figures as pf
from kirsch_reference import annulus
from setup_drawings import schematic

sys.path.insert(0, str(pf.ROOT/'tests/functional/hydromechanics'))
import check_kirsch as ck


def build(study, output, record_source):
    (output/'figures').mkdir(parents=True, exist_ok=True)
    report = json.loads((study/'qualification.json').read_text())
    assert report['passed']
    record_source(study/'qualification.json')
    for case in report['cases']:
        source = study/case['folder']
        assert ck.digest(source/'input.cfg') == case['input_sha256']
        assert ck.digest(source/'annulus.poly') == case['mesh_input_sha256']
        replay = ck.compare(source)
        for name, value in case['errors'].items():
            np.testing.assert_allclose(replay['errors'][name], value, rtol=0, atol=1e-12)
        for name, digest in case['outputs_sha256'].items():
            assert ck.digest(source/name) == digest
            record_source(source/name)
        for name in ('input.cfg', 'annulus.poly'): record_source(source/name)
    fine = min(report['cases'], key=lambda c:c['spacing_m'])
    folder = study/fine['folder']
    actual = ck.compare(folder)
    for key, value in fine['errors'].items():
        np.testing.assert_allclose(actual['errors'][key], value, rtol=0, atol=1e-12)
    for path in folder.glob('kirsch.*'):
        if '.save.' in path.name or path.suffix == '.info': record_source(path)
    for name in ('input.cfg', 'annulus.poly'): record_source(folder/name)
    data = output/'data/kirsch'; data.mkdir(parents=True, exist_ok=True)
    for name in ('input.cfg', 'annulus.poly'): shutil.copy2(folder/name, data/name)
    shutil.copy2(study/'qualification.json', data/'qualification.json')
    shutil.copy2(Path(ck.__file__), data/'check_kirsch.py')
    shutil.copy2(Path(__file__).with_name('kirsch_reference.py'), data/'kirsch_reference.py')
    model = ck.Dynearthsol(str(folder/'kirsch'))
    x = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    centre = x[conn].mean(axis=1); r = np.linalg.norm(centre, axis=1)
    nx, nz = (centre/r[:, None]).T
    selected = (r >= 1.1*ck.A) & (r <= 5*ck.A)

    def polar(values):
        rr = nx**2*values[:, 0]+nz**2*values[:, 1]+2*nx*nz*values[:, 2]
        hoop = nz**2*values[:, 0]+nx**2*values[:, 1]-2*nx*nz*values[:, 2]
        return rr[selected], hoop[selected]

    stress = polar(model.read_field(model.frames[-1], 'stress'))
    strain = polar(model.read_field(model.frames[-1], 'strain'))
    order = np.argsort(r[selected])
    keep = order[np.linspace(0, len(order)-1, min(90, len(order))).astype(int)]
    grid = np.linspace(ck.A, 5*ck.A, 301)
    rr, hoop = annulus(grid, ck.A, ck.R, ck.P)
    lame, modulus = 1e9/3, 1e9
    ev = (rr+hoop)/(2*(lame+modulus))
    er, et = (rr-lame*ev)/(2*modulus), (hoop-lame*ev)/(2*modulus)
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        for values, exact, label, color in ((stress[0], rr, 'radial', ds.WATER), (stress[1], hoop, 'hoop', ds.HOT)):
            ds.analytical(axes[0], grid, -exact/ck.P, label=label, color=color)
            ds.des(axes[0], r[selected][keep], -values[keep]/ck.P, color=color, size=3, label='_nolegend_')
        ds.panel(axes[0], 'a', 'Radial and hoop stress', xlabel='$r/a$', ylabel='compression / $P$')
        ds.key(axes[0])
        for values, exact, label, color in ((strain[0], er, 'radial', ds.WATER), (strain[1], et, 'hoop', ds.HOT)):
            ds.analytical(axes[1], grid, exact/(ck.P/modulus), label=label, color=color)
            ds.des(axes[1], r[selected][keep], values[keep]/(ck.P/modulus), color=color, size=3, label='_nolegend_')
        ds.panel(axes[1], 'b', 'Elastic strain', xlabel='$r/a$', ylabel=r'strain / $(P/G)$')
        ds.key(axes[1])
        cases = sorted(report['cases'], key=lambda c:c['spacing_m'])
        h = np.array([c['spacing_m'] for c in cases])
        for field, label, color, marker in (('radial_rms_over_load', 'radial', ds.WATER, 'o'),
                                            ('hoop_rms_over_load', 'hoop', ds.HOT, 's')):
            axes[2].loglog(h, [c['errors'][field] for c in cases], marker+'-', color=color, mfc='white', label=label)
        axes[2].loglog(h, cases[0]['errors']['hoop_rms_over_load']*h/h[0], '--', color=ds.MUTED, label='slope 1')
        axes[2].set_xticks(h, [f'{v:g}' for v in h]); axes[2].xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ds.panel(axes[2], 'c', 'Stress convergence', xlabel='inner spacing (m)', ylabel='RMS error / $P$')
        pf.save_figure(fig, output/'figures', 'kirsch')
    return dict(key='kirsch', title='Elastic circular opening: Kirsch limit', theme='tunnels',
                problem='A circular opening concentrates stress in an elastic solid under isotropic compression.',
                result=f"Radial/hoop stress RMS errors: {100*fine['errors']['radial_rms_over_load']:.3f}% / {100*fine['errors']['hoop_rms_over_load']:.3f}% of $P$; largest local error {100*fine['errors']['max_error_over_load']:.2f}%.",
                scope=f"2D plane strain, $a=1$ m and $R=20a$, {fine['spacing_m']:g} m inner spacing. Compared with the finite-annulus solution; its wall hoop stress differs from infinite-domain Kirsch by 0.251%.",
                meaning='This checks dry elastic stress, strain, curved-surface traction and symmetry constraints. Plasticity and excavation history remain separate checks; poroelastic injection is on the borehole page.',
                metric=fine['errors']['hoop_rms_over_load'], metric_label='hoop stress RMS error / P',
                dimension=2, spacing_m=fine['spacing_m'], evidence=json.dumps(report))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', type=Path); parser.add_argument('output', type=Path)
    args = parser.parse_args()
    result = build(args.study.resolve(), args.output.resolve(), lambda p:None)
    target = args.output/'fig/kirsch'; target.mkdir(parents=True, exist_ok=True)
    schematic('kirsch', target)
    for extension in ('png', 'svg', 'pdf'):
        (args.output/'figures'/f'kirsch.{extension}').replace(target/f'result.{extension}')
    if not any((args.output/'figures').iterdir()): (args.output/'figures').rmdir()
    print(json.dumps(result, indent=2))
