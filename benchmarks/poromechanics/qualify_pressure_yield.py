"""Validate the recorded DES pressure-driven yielding studies and plot them.

Interface as qualify_theis: build(study, output, record_source, steps, extra=())
replays the finest run from raw saves, writes figures/pressure_yield.* and
data/pressure_yield/, and returns the review record. study is the mesh
sequence at 400 steps, steps the 1600-step run on 0.05 m, and extra any
further step/mesh refinements (6400 steps, 0.025 m).
"""
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys

import numpy as np

import prepare_figures as pf
import diagram_style as ds
sys.path.insert(0, str(pf.ROOT/'tests/functional/hydromechanics'))
import check_pressure_yield as cp

KEY = 'pressure_yield'


def _load(study, record_source):
    summary = json.loads((study/'summary.json').read_text())
    assert hashlib.sha256((study/'runner.py').read_bytes()).hexdigest() == summary['runner_sha256']
    record_source(study/'summary.json'); record_source(study/'runner.py')
    return summary


def _history(folder, element_depths=(0.0, 1.0)):
    """Yield-front depth at every saved step, and the effective stresses of two elements (2D)."""
    model = cp.Dynearthsol(str(folder/'layer')); rows = cp.scan_frames(str(folder/'layer'))
    coord = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    depth = -coord[conn][:, :, -1].mean(axis=1)
    p0 = model.read_field(0, 'pore pressure')[conn].mean(axis=1)
    picks = [int(np.argmin(np.abs(depth - d))) for d in element_depths]
    front, paths = [], {i: [] for i in picks}
    for row in rows[1:]:
        p = model.read_field(row['frame'], 'pore pressure')[conn].mean(axis=1)
        dp = p - p0
        plastic = model.read_field(row['frame'], 'plastic strain') > 0
        front.append((float(np.mean(dp)), float(depth[plastic].max()) if plastic.any() else np.nan))
        s = model.read_field(row['frame'], 'stress')                    # 2D: xx, zz, xz (tension positive)
        for i in picks:
            paths[i].append((float(dp[i]), s[i, 1] + cp.REF.alpha*p[i], s[i, 0] + cp.REF.alpha*p[i]))
    return np.array(front), {float(depth[i]): np.array(v) for i, v in paths.items()}


def build(study, output, record_source, steps, extra=(), slabs=()):
    summary, fine = _load(study, record_source), _load(steps, record_source)
    more = [(_load(s, record_source), s) for s in extra]
    runs = [(c, summary['steps'], study) for c in summary['cases']] + [(c, fine['steps'], steps) for c in fine['cases']]
    runs += [(c, s['steps'], root) for s, root in more for c in s['cases']]
    # The same mesh and step count under the other PT option: keep the first run, report the agreement.
    seen, unique, twins = {}, [], []
    for r in runs:
        k = (r[0]['spacing_m'], r[1])
        if k in seen: twins.append((seen[k], r))
        else: seen[k] = r; unique.append(r)
    runs = unique
    twin_text = ''.join(f" With PT option 1 instead of 0 ({a[0]['spacing_m']:g} m, {a[1]} steps) the mean error is "
                        f"{100*b[0]['mean_error']:.3f}% against {100*a[0]['mean_error']:.3f}%." for a, b in twins)
    mesh = sorted(summary['cases'], key=lambda c: -c['spacing_m'])
    best = min(runs, key=lambda r: (r[0]['spacing_m'], -r[1]))         # finest mesh, most steps
    for case, n, root in runs:
        folder = root/f"h{case['spacing_m']:g}"
        for path in folder.glob('layer.*'):
            if path.suffix in ('.cfg', '.info') or '.save.' in path.name: record_source(path)
    best_folder = best[2]/f"h{best[0]['spacing_m']:g}"
    again = cp.compare(best_folder)
    for k in ('yielded', 'max_error', 'mean_error'):
        np.testing.assert_allclose(again[k], best[0][k], rtol=1e-12)

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    shutil.copy2(best_folder/'layer.cfg', data/'input.cfg')
    shutil.copy2(study/'summary.json', data/'summary.json'); shutil.copy2(steps/'summary.json', data/'summary-steps.json')
    for s, r in more: shutil.copy2(r/'summary.json', data/f"summary-{r.name}.json")
    shutil.copy2(study/'runner.py', data/'check_pressure_yield.py')
    shutil.copy2(Path(__file__).with_name('pressure_yield_reference.py'), data/'pressure_yield_reference.py')

    ref, dp_max = cp.REF, summary['dp_max']
    front_best, paths = _history(best_folder)
    front_coarse, _ = _history(study/f"h{mesh[0]['spacing_m']:g}")
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        ax = axes[0]
        depth = np.linspace(0, cp.H, 201); on = ref.onset(depth)
        reach = on <= dp_max
        ds.analytical(ax, on[reach]/1e3, depth[reach], label='exact yield front')
        for front, color, label in ((front_coarse, ds.SERIES[1], f"DES, {mesh[0]['spacing_m']:g} m, {summary['steps']} steps"),
                                    (front_best, ds.HOT, f"DES, {best[0]['spacing_m']:g} m, {best[1]} steps")):
            ok = ~np.isnan(front[:, 1])
            ax.plot(front[ok, 0]/1e3, front[ok, 1], drawstyle='steps-post', color=color, lw=1.1, label=label)
        ax.invert_yaxis(); ax.set_xlim(on[0]/1e3 - 1, dp_max/1e3 + .5)
        ds.panel(ax, 'a', 'Depth of the yield front', xlabel=r'pressure rise $\Delta p$ (kPa)', ylabel='depth (m)')
        ax.legend(fontsize=10, loc='lower left')

        ax = axes[1]
        apex = ref.c/math.tan(ref.phi)/1e3
        sline = np.array([-apex, 30])
        ax.plot(sline, (ref.c*math.cos(ref.phi) + sline*1e3*math.sin(ref.phi))/1e3, color=ds.INK, lw=1.2, label='Mohr–Coulomb')
        for color, (d, path) in zip((ds.SERIES[0], ds.SERIES[2]), sorted(paths.items())):
            grid = np.linspace(0, ref.onset(d), 60)
            sv, sh = ref.stresses(d, grid)
            ax.plot(-(sv + sh)/2e3, np.abs(sv - sh)/2e3, '-', color=color, lw=1.0, label=f'exact path, $d={d:.2f}$ m')
            # After the surcharge has drained: the first steps carry its undrained transient.
            settled = path[path[:, 0] > .02*ref.onset(d)]
            keep = settled[::max(1, len(settled)//40)]
            ax.plot(-(keep[:, 1] + keep[:, 2])/2e3, np.abs(keep[:, 1] - keep[:, 2])/2e3, 'o', ms=3, mfc='white',
                    mec=color, mew=.8)
        ax.plot([], [], 'o', ms=3, mfc='white', mec=ds.INK, label='DES element')
        ax.axvline(0, color=ds.MUTED, lw=.6, ls=':')
        ax.set_xlim(-apex - 2, 30); ax.set_ylim(0, 16)
        ds.panel(ax, 'b', 'Effective stress paths to yield', xlabel=r"$-(\sigma'_1+\sigma'_3)/2$ (kPa)",
                 ylabel=r"$(\sigma'_1-\sigma'_3)/2$ (kPa)")
        ax.legend(fontsize=10, loc='upper left')

        ax = axes[2]
        by_mesh = {}
        for case, n, _ in runs: by_mesh.setdefault(case['spacing_m'], []).append((n, case))
        for color, h in zip(ds.SERIES, sorted(by_mesh, reverse=True)):
            seq = sorted(by_mesh[h], key=lambda r: r[0])
            n = np.array([r[0] for r in seq])
            ax.loglog(n, [100*r[1]['mean_error'] for r in seq], 'o-', color=color, mfc='white', label=f'{h:g} m, mean')
            ax.loglog(n, [100*r[1]['max_error'] for r in seq], 's:', color=color, mfc='white', ms=4)
        ax.plot([], [], 's:', color=ds.MUTED, mfc='white', ms=4, label='largest (dotted)')
        ticks = [.2, .3, .5, 1, 2]
        ax.set_yticks(ticks, [f'{v:g}' for v in ticks]); ax.yaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        steps_all = sorted({r[1] for r in runs})
        ax.set_xticks(steps_all, [f'{v:g}' for v in steps_all]); ax.xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ds.panel(ax, 'c', 'Late yield by step count and mesh', xlabel='steps over the loading',
                 ylabel=r'$(\Delta p-\Delta p^*)/\Delta p^*$ (%)')
        ax.legend(fontsize=10)
        pf.save_figure(fig, output/'figures', KEY)

    b = best[0]
    slab_cases = []
    for slab in slabs:
        saved = _load(slab, record_source)
        case = saved['cases'][0]
        assert saved['ndims'] == 3
        folder = slab/f"h{case['spacing_m']:g}"
        again = cp.compare(folder)
        for name in ('yielded', 'max_error', 'mean_error'):
            np.testing.assert_allclose(again[name], case[name], rtol=1e-10)
        for path in folder.glob('layer.*'):
            if path.suffix in ('.cfg', '.info') or '.save.' in path.name: record_source(path)
        mode = saved.get('mode', 'pt')
        shutil.copy2(slab/'summary.json', data/f'summary-3d-{mode}.json')
        slab_cases.append(dict(mode=mode, steps=saved['steps'], scaling=saved.get('scaling'), **case))
    slab_text = (' 3D checks: ' + '; '.join(
        f"{c['mode']}, {c['spacing_m']:g} m, {c['steps']} steps: {c['yielded']}/{c['elements']} yielded, "
        f"maximum onset error {100*c['max_error']:.3f}%" for c in slab_cases) + '. ') if slab_cases else ''
    return dict(key=KEY, theme='mechanics', title='Pressure-driven yielding of a confined layer',
        problem=('Rising pore pressure lowers the effective stress until the ground yields, the mechanism of injection- and '
                 'rainfall-induced failure. In a laterally confined layer the yield pressure at each depth is exact, so a '
                 'yield front moving down can be compared element by element.'),
        result=(f"On the {b['spacing_m']:g} m mesh with {best[1]} steps every element yields within {100*b['max_error']:.2f}% "
                f"(mean {100*b['mean_error']:.2f}%) of its exact pressure, and {b['yielded']} of {b['elements']} elements "
                'yield, as predicted. DES always yields slightly late; the delay falls with both the step and the mesh.'),
        scope=('2D plane strain, uniaxial strain, elastic–perfectly plastic Mohr–Coulomb without dilation, fully saturated. '
               'Before yield the DES effective stresses follow the exact path within 0.05%, so the remaining error is the '
               'onset itself: it is read at the first saved step with plastic strain (up to one load increment late) and '
               'carries the one-step lag of the sequential coupling; both shrink with the step, and a spatial part with the '
               'mesh. Mean errors: ' + '; '.join(f"{c['spacing_m']:g} m, {n} steps: {100*c['mean_error']:.2f}%"
                                                for c, n, _ in sorted(runs, key=lambda r: (-r[0]['spacing_m'], r[1])))
               + '.' + twin_text + slab_text + '(a) yield front; (b) stress paths; (c) error.'),
        meaning=('This checks effective-stress plasticity in a boundary-value problem driven by pore pressure: the Biot '
                 'coefficient in the effective stress, the Mohr–Coulomb return and the pressure coupling together.'),
        metric=b['max_error'], metric_label='yield-pressure error', dimension=2, spacing_m=b['spacing_m'],
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 dp_max=dp_max, sanity_3d=slab_cases, runs=[dict(steps=n, **{k: v for k, v in c.items() if k != 'profile'})
                                                      for c, n, _ in runs]), indent=1))
