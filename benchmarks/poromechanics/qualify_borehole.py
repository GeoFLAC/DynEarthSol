"""Validate the recorded DES borehole-injection study and plot it.

Interface as qualify_theis: build(study, output, record_source, steps)
replays the finest mesh from raw saves, writes figures/borehole.* and
data/borehole/, and returns the review record. steps is the finer-step
study on one mesh.
"""
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np

import prepare_figures as pf
import diagram_style as ds
sys.path.insert(0, str(pf.ROOT/'tests/functional/hydromechanics'))
import check_borehole as cb

KEY = 'borehole'


def _load(study, record_source):
    summary = json.loads((study/'summary.json').read_text())
    assert hashlib.sha256((study/'runner.py').read_bytes()).hexdigest() == summary['runner_sha256']
    record_source(study/'summary.json'); record_source(study/'runner.py')
    return summary


def _binned(r, values, edges):
    idx = np.digitize(r, edges) - 1
    keep = [(edges[i] + edges[i+1])/2 for i in range(len(edges) - 1) if np.any(idx == i)]
    mean = [values[idx == i].mean() for i in range(len(edges) - 1) if np.any(idx == i)]
    return np.array(keep), np.array(mean)


def build(study, output, record_source, steps):
    summary, fine = _load(study, record_source), _load(steps, record_source)
    cases = sorted(summary['cases'], key=lambda c: -c['spacing_m'])
    finest, refined = cases[-1], fine['cases'][0]
    replay_difference = 0.
    for root, case in ((study, finest), (steps, refined)):
        folder = root/f"h{case['spacing_m']:g}"
        assert cb.digest(folder/'hole.cfg') == case['input_sha256']
        assert cb.digest(folder/'annulus.poly') == case['mesh_sha256']
        for path in folder.glob('hole.*'):
            if path.suffix in ('.cfg', '.info') or '.save.' in path.name: record_source(path)
        record_source(folder/'annulus.poly')
        frames, _, _ = cb.compare(folder)
        actual = np.array([f['pressure_error'] for f in frames])
        saved = np.array([f['pressure_error'] for f in case['frames']])
        # Stehfest's alternating sum amplifies math-library roundoff. Bound
        # this summary replay in units of wall pressure, not the small error.
        np.testing.assert_allclose(actual, saved, rtol=1e-10, atol=1e-8)
        replay_difference = max(replay_difference, float(abs(actual-saved).max()))

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    folder = study/f"h{finest['spacing_m']:g}"
    shutil.copy2(folder/'hole.cfg', data/'input.cfg'); shutil.copy2(folder/'annulus.poly', data/'annulus.poly')
    shutil.copy2(study/'summary.json', data/'summary.json'); shutil.copy2(steps/'summary.json', data/'summary-steps.json')
    shutil.copy2(study/'runner.py', data/'check_borehole.py')
    shutil.copy2(Path(__file__).with_name('borehole_reference.py'), data/'borehole_reference.py')

    ref, pw = cb.REF, cb.PW
    model = cb.Dynearthsol(str(folder/'hole')); rows = cb.scan_frames(str(folder/'hole'))
    x = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    rn = np.linalg.norm(x, axis=1); centre = x[conn].mean(axis=1); re = np.linalg.norm(centre, axis=1)
    nx, nz = (centre/re[:, None]).T
    edges = np.linspace(cb.A, 6.0, 41)
    rr = np.linspace(cb.A, 6.0, 400)
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        ax = axes[0]
        for color, row in zip(ds.SERIES, rows[1:]):
            t = row['time']
            ds.analytical(ax, rr, ref.fields(rr, t)['p']/pw, color=color, label=f'$t={t:g}$ s')
            rb, pb = _binned(rn, model.read_field(row['frame'], 'pore pressure')/pw, edges)
            ds.des(ax, rb, pb, color=color, label=None)
        ds.panel(ax, 'a', 'Pore pressure', xlabel='$r/a$', ylabel='$p/p_w$')
        ds.key(ax, loc='upper right')

        ax = axes[1]
        row = rows[-1]; t = row['time']
        s = model.read_field(row['frame'], 'stress')
        srr = nx**2*s[:, 0] + nz**2*s[:, 1] + 2*nx*nz*s[:, 2]
        stt = nz**2*s[:, 0] + nx**2*s[:, 1] - 2*nx*nz*s[:, 2]
        f = ref.fields(rr, t)
        for color, name, exact, des in ((ds.WATER, r'$\sigma_{rr}$', f['srr'], srr), (ds.HOT, r'$\sigma_{\theta\theta}$', f['stt'], stt)):
            ds.analytical(ax, rr, exact/pw, color=color, label=name)
            rb, sb = _binned(re, des/pw, edges)
            ds.des(ax, rb, sb, color=color, label=None)
        ax.axhline(0, color=ds.MUTED, lw=.6)
        ds.panel(ax, 'b', f'Total stress at $t={t:g}$ s', xlabel='$r/a$', ylabel=r'$\sigma/p_w$ (tension positive)')
        ds.key(ax, loc='lower right')

        ax = axes[2]
        h = np.array([c['spacing_m'] for c in cases])
        last = lambda c, k: c['frames'][-1][k]
        ax.loglog(h, [100*last(c, 'pressure_error') for c in cases], 'o-', color=ds.WATER, mfc='white', label='pressure, max')
        ax.loglog(h, [100*last(c, 'pressure_rms') for c in cases], 'o:', color=ds.WATER, mfc='white', label='pressure, rms')
        ax.loglog(h, [100*last(c, 'hoop_error') for c in cases], 's-', color=ds.HOT, mfc='white', label='hoop stress, max')
        ax.loglog(h, [100*last(c, 'hoop_rms') for c in cases], 's:', color=ds.HOT, mfc='white', label='hoop stress, rms')
        ax.loglog([refined['spacing_m']], [100*last(refined, 'pressure_error')], 'D', color=ds.WATER, ms=5,
                  label=f"{fine['steps']} steps")
        ax.set_xticks(h, [f'{v:g}' for v in h]); ax.xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ds.panel(ax, 'c', f'Errors at $t={t:g}$ s by mesh', xlabel='mesh spacing $h/a$', ylabel='error (% of $p_w$)')
        ax.legend(fontsize=10)
        pf.save_figure(fig, output/'figures', KEY)

    worst = lambda c, k: max(f[k] for f in c['frames'])
    return dict(key=KEY, theme='tunnels', title='Fluid injection into a poroelastic borehole',
        problem=('Fluid at pressure $p_w$ enters a borehole wall: it loads the wall and diffuses into the rock, which changes the '
                 'stresses around the hole as it goes. The axisymmetric part of Detournay and Cheng’s borehole problem, on a '
                 'finite annulus.'),
        result=(f"On the {finest['spacing_m']:g} a mesh: pressure within {100*worst(finest, 'pressure_error'):.2f}% of $p_w$ at every "
                f"node for 1.05 a ≤ r ≤ 4 a and all four times; total stress rms {100*worst(finest, 'hoop_rms'):.2f}% "
                f"(hoop) with a largest element error of {100*finest['stress_error']:.2f}%."),
        scope=('2D plane strain, quarter annulus a = 1 m, R = 20 m, sealed traction-free outer arc, PT. Largest pressure error on '
               'meshes ' + ', '.join(f"{c['spacing_m']:g}" for c in cases) + ' a: '
               + ', '.join(f"{100*worst(c, 'pressure_error'):.3f}%" for c in cases)
               + '; largest stress error ' + ', '.join(f"{100*c['stress_error']:.2f}%" for c in cases)
               + f" (constant-stress triangles at the wall, where the hoop stress is steepest). A time-step check with {fine['steps']} steps on the "
               f"{refined['spacing_m']:g} a mesh gives pressure {100*worst(refined, 'pressure_error'):.3f}%. "
               '3D borehole accuracy is outside this study. (a) pressure; (b) stress; (c) errors.'),
        meaning=('This checks the coupled poroelastic response around an opening: diffusion from a curved boundary, the '
                 'pressure-driven hoop stress and the wall load together. It precedes the tunnel-excavation problems.'),
        metric=worst(finest, 'pressure_error'), metric_label='pressure error / p_w', dimension=2, spacing_m=finest['spacing_m'],
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 diffusivity=summary['diffusivity'], cases=cases, refined=refined, refined_steps=fine['steps'],
                                 reference_replay_difference_over_wall_pressure=replay_difference,
                                 reference_replay_absolute_bound=1e-8),
                            indent=1))
