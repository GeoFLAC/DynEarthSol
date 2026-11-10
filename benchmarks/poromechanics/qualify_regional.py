"""Validate the recorded regional-subsidence application and plot its finest mesh.

Interface as qualify_theis: build(study, output, record_source, times, final)
with the mesh sequence, a time-step check (one mesh, four times the steps)
and the final run on the finest mesh with those steps, which supplies the
bowl, the wells and the headline.
"""
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np

import prepare_figures as pf
import diagram_style as ds
import regional_subsidence as rs

KEY = 'regional'


def _check(study, record_source):
    summary = json.loads((study/'summary.json').read_text())
    assert hashlib.sha256((study/'runner.py').read_bytes()).hexdigest() == summary['runner_sha256']
    assert rs.digest(Path(rs.__file__)) == summary['runner_sha256']
    record_source(study/'summary.json'); record_source(study/'runner.py')
    for case in summary['cases']:
        folder = study/f"h{case['spacing_m']:g}"
        assert rs.digest(folder/'basin.cfg') == case['input_sha256']
        frames, nodes = rs.compare(folder)
        assert nodes == case['nodes'] and frames == case['frames']
        assert max(f['pressure_error'] for f in frames) == case['pressure_error']
        assert max(f['subsidence_error'] for f in frames) == case['subsidence_error']
        for path in folder.glob('basin.*'):
            if path.suffix in ('.cfg', '.info') or '.save.' in path.name: record_source(path)
    return summary


def build(study, output, record_source, times, final):
    summary, steps, last = (_check(s, record_source) for s in (study, times, final))
    mesh = sorted(summary['cases'], key=lambda c: -c['spacing_m'])
    coarse_finest = mesh[-1]
    finest = last['cases'][0]; h = finest['spacing_m']
    assert abs(h - coarse_finest['spacing_m']) < 1e-9 and last['steps'] == steps['steps']
    folder = final/f'h{h:g}'
    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    shutil.copy2(folder/'basin.cfg', data/'input.cfg'); shutil.copy2(study/'summary.json', data/'summary.json')
    shutil.copy2(final/'summary.json', data/'summary-final.json'); shutil.copy2(times/'summary.json', data/'summary-steps.json')
    if (final/'NOTE.txt').exists(): shutil.copy2(final/'NOTE.txt', data/'NOTE-final.txt'); record_source(final/'NOTE.txt')
    shutil.copy2(study/'runner.py', data/'regional_subsidence.py')

    model = rs.Dynearthsol(str(folder/'basin')); rows = rs.scan_frames(str(folder/'basin'))
    coord = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    centre = coord[conn].mean(axis=1); re = np.hypot(centre[:, 0], centre[:, 1])
    year = rs.YEAR2SEC
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        ax = axes[0]
        rr = np.linspace(100, 2500, 300)
        ends = [rs.SEASON, 2*rs.SEASON, 3*rs.SEASON, 4*rs.SEASON]
        names = ('end of pumping 1', 'end of rest 1', 'end of pumping 2', 'end of rest 2')
        for color, target, name in zip(ds.SERIES, ends, names):
            row = min(rows[1:], key=lambda r: abs(r['time'] - target))
            t = row['time']
            ds.analytical(ax, rr/1e3, 100*rs.reference(rr, t)['subsidence'], color=color, label=f'{name} ({t/year:.2f} yr)')
            sub = -rs.B*model.read_field(row['frame'], 'strain')[:, 2]
            band = (re > 100) & (re < 2500) & (np.abs(centre[:, 1]) < .5*h)    # elements along the x axis
            order = np.argsort(re[band]); keep = order[::max(1, len(order)//18)]
            ax.plot(re[band][keep]/1e3, 100*sub[band][keep], 'o', ms=3.5, mfc='white', mec=color, mew=.9)
        ds.band(ax, rs.BAND[0]/1e3, rs.BAND[1]/1e3)
        ds.panel(ax, 'a', 'Subsidence bowl', xlabel='distance from the well field (km)', ylabel='subsidence (cm)')
        ds.key(ax)
        ax = axes[1]
        tt = np.linspace(0, rs.DURATION, 400)
        for color, r in zip((ds.SERIES[0], ds.SERIES[2], ds.SERIES[3]), rs.WELLS):
            near = np.argsort(np.abs(re - r))[:8]
            ds.analytical(ax, tt/year, 100*np.array([rs.reference(re[near], t)['subsidence'].mean() for t in tt]),
                          color=color, label=f'$r={r:g}$ m')
            ax.plot([f['time']/year for f in finest['frames']], [100*f['wells'][f'{r:g}']['des'] for f in finest['frames']],
                    'o', ms=3.5, mfc='white', mec=color, mew=.9)
        for k in range(1, 4):
            ax.axvline(k*rs.SEASON/year, color='#9aa5ab', lw=.8, ls=':')
        ds.panel(ax, 'b', 'Observation wells', xlabel='time (yr)', ylabel='subsidence (cm)')
        ds.key(ax)
        ax = axes[2]
        hs = np.array([c['spacing_m'] for c in mesh])
        after = lambda c, k: max(f[k] for f in c['frames'][1:])          # after the first frame
        ax.loglog(hs, [after(c, 'pressure_error') for c in mesh], 'o-', color=ds.WATER, mfc='white', label='pressure')
        ax.loglog(hs, [after(c, 'subsidence_error') for c in mesh], 's-', color=ds.HOT, mfc='white', label='subsidence')
        for t in (steps['cases'][0], finest):
            ax.loglog([t['spacing_m']], [after(t, 'pressure_error')], 'o', color=ds.WATER, ms=7)
            ax.loglog([t['spacing_m']], [after(t, 'subsidence_error')], 's', color=ds.HOT, ms=7)
        ax.plot([], [], 'o', color=ds.INK, ms=6, mfc='white', label=f"{summary['steps']} steps (open)")
        ax.plot([], [], 'o', color=ds.INK, ms=6, label=f"{steps['steps']} steps (solid)")
        ax.set_xticks(hs, [f'{x:g}' for x in hs]); ax.xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ds.panel(ax, 'c', 'Mesh and time step', xlabel='mesh spacing $h$ (m)',
                 ylabel='largest error after the first frame')
        ax.legend(fontsize=10.5)
        pf.save_figure(fig, output/'figures', KEY)

    peak = max(f['peak_subsidence_ref'] for f in finest['frames'])
    late_sub = max(f['subsidence_error'] for f in finest['frames'][1:])
    well = lambda c, r: max(abs(f['wells'][r]['des'] - f['wells'][r]['ref']) for f in c['frames'])
    late_p = max(f['pressure_error'] for f in finest['frames'][1:])
    wells = {f'{r:g}': max(abs(f['wells'][f'{r:g}']['des'] - f['wells'][f'{r:g}']['ref']) for f in finest['frames'])
             for r in rs.WELLS}
    order = float(np.polyfit(np.log(hs), np.log([c['subsidence_error'] for c in mesh]), 1)[0])
    same = next(c for c in mesh if abs(c['spacing_m'] - steps['cases'][0]['spacing_m']) < 1e-9)
    wells_104 = {f'{r:g}': well(coarse_finest, f'{r:g}') for r in rs.WELLS}
    return dict(key=KEY, theme='consolidation', status='preliminary', title='Application: seasonal pumping and land subsidence',
        problem=('A well field pumps a soft confined aquifer for six months, rests for six months and repeats. The ground '
                 'subsides around the wells, partly rebounds in each rest season and keeps a growing subsidence farther out.'),
        result=(f"Largest subsidence error at the observation wells on the {h:g} m mesh (peak subsidence {100*peak:.1f} cm): "
                + ', '.join(f'{100*v:.2f} cm at {r} m' for r, v in wells.items())
                + f". Over {rs.BAND[0]:g}–{rs.BAND[1]:g} m after the first frame, the largest subsidence error is "
                f"{100*late_sub:.1f}% of the largest subsidence in the same frame, and the largest pressure error {100*late_p:.1f}% "
                f"of the largest pressure drop; at the first frame ({finest['frames'][0]['time']/year:.2f} yr, a cone only a few "
                f"elements wide) the subsidence error is {100*finest['frames'][0]['subsidence_error']:.1f}%."),
        scope=(f"3D quarter of a {2*rs.L/1e3:g} km square, {rs.B:g} m thick; Q = {rs.Q:g} m$^3$/s; PT (adaptive dynamic "
               f"relaxation). The run in the result above uses the {h:g} m mesh with {last['steps']} steps; with {summary['steps']} steps its "
               "well errors were " + ', '.join(f'{r} m: {100*v:.2f} cm' for r, v in wells_104.items())
               + f". Meshes at {summary['steps']} steps " + ', '.join(f'{x:g}' for x in hs)
               + " m: subsidence errors after the first frame " + ', '.join(f"{100*after(c, 'subsidence_error'):.1f}%" for c in mesh)
               + f". The step dominates near each rate change: on the {steps['cases'][0]['spacing_m']:g} m mesh, {steps['steps']} steps give pressure "
               f"{100*after(steps['cases'][0], 'pressure_error'):.1f}% and subsidence {100*after(steps['cases'][0], 'subsidence_error'):.1f}% "
               f"after the first frame (against {100*after(same, 'pressure_error'):.1f}% and {100*after(same, 'subsidence_error'):.1f}%), "
               f"and the worst 250 m well error falls from {100*well(same, '250'):.2f} to {100*well(steps['cases'][0], '250'):.2f} cm. "
               "(a) bowl at the ends of the seasons; (b) observation wells; (c) errors by mesh and step: open markers "
               f"{summary['steps']} steps, solid {last['steps']}."),
        meaning=('This applies the verified pumping and averaged-aquifer physics at field scale, with a realistic schedule '
                 'and rebound. The reference is the averaged model of a closed basin: fixed outer walls add a uniform '
                 'areal strain $h(t)$ that deepens the drawdown everywhere, which the infinite-domain line well misses by '
                 'up to about 20% in the rest seasons. Only changes from the initial state are computed; gravity, an overlying aquitard and inelastic '
                 'compaction of clays, which set most real subsidence, are not modelled.'),
        metric=finest['subsidence_error'], metric_label='subsidence error / max', dimension=3, spacing_m=h,
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 order=order, wells_max_abs_m=wells, outer_radius_m=rs.OUTER,
                                 time_check=dict(steps=steps['steps'], **{k: v for k, v in steps['cases'][0].items() if k != 'frames'}),
                                 final=dict(steps=last['steps'], executable_sha256=last['executable_sha256'],
                                            **{k: v for k, v in finest.items() if k != 'frames'}),
                                 cases=[{k: v for k, v in c.items() if k != 'frames'} for c in mesh]), indent=1))
