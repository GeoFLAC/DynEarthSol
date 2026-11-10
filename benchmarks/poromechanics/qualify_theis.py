"""Validate the recorded DES Theis study and plot its finest mesh.

Interface as qualify_convection: build(study, output, record_source) recomputes
errors from the raw pressure fields, writes figures/theis.* and data/theis/,
and returns the review record. SPEC updates review.SPECS['theis'] from an
analytical-only reference to the numerical comparison.
"""
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np

import prepare_figures as pf
import diagram_style as ds
from theis_reference import scheduled_drawdown
sys.path.insert(0, str(pf.ROOT/'tests/functional/hydromechanics'))
import check_theis as ct

KEY = 'theis'


def build(study, output, record_source):
    summary = json.loads((study/'summary.json').read_text())
    assert hashlib.sha256((study/'runner.py').read_bytes()).hexdigest() == summary['runner_sha256']
    record_source(study/'summary.json'); record_source(study/'runner.py')
    ref = summary['reference']
    T, S, scale = ref['transmissivity'], ref['storativity'], ref['pressure_scale_Pa']
    schedule = [(0.0, ref['Q']), (ref['pump_until_s'], 0.0)]
    finest = min(summary['cases'], key=lambda c: c['spacing_m'])
    folder = study/f"h{finest['spacing_m']:g}"
    for path in folder.glob('theis.*'):
        if path.suffix in ('.cfg', '.info') or '.save.' in path.name: record_source(path)
    model = ct.Dynearthsol(str(folder/'theis')); rows = ct.scan_frames(str(folder/'theis'))
    coord = model.read_field(0, 'coordinate'); p0 = model.read_field(0, 'pore pressure')
    conn = model.read_field(0, 'connectivity')
    vertices = coord[conn]
    volume = np.abs(np.linalg.det(vertices[:, 1:]-vertices[:, :1]))/6
    r = np.hypot(coord[:, 0], coord[:, 1]); window = (r >= 10) & (r <= ref['L']/2)
    history, budget = [], []
    assert len(rows)-1 == len(finest['frames'])
    for row, recorded in zip(rows[1:], finest['frames']):
        s = p0 - model.read_field(row['frame'], 'pore pressure')
        exact = scheduled_drawdown(r[window], np.full(window.sum(), row['time']), schedule, T, S)
        err = float(np.max(np.abs(s[window]-exact)))/scale
        np.testing.assert_allclose(err, recorded['max_error_over_scale'], rtol=0, atol=1e-12)
        history.append((row['time'], s))
        stored = -ref['storage_per_Pa']*float(volume @ np.mean(s[conn], axis=1))
        expected_volume = -ref['Q']/4*min(row['time'], ref['pump_until_s'])
        budget.append((row['time'], stored, expected_volume))
        for field in ('pore pressure', 'stress', 'strain', 'velocity'):
            assert np.isfinite(model.read_field(row['frame'], field)).all(), (folder, field)

    budget_error = max(abs(actual-exact) for _, actual, exact in budget)
    budget_scale = ref['Q']/4*ref['pump_until_s']
    assert budget_error/budget_scale < 1e-10

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    shutil.copy2(folder/'theis.cfg', data/'input.cfg'); shutil.copy2(study/'summary.json', data/'summary.json')
    shutil.copy2(study/'runner.py', data/'check_theis.py')
    shutil.copy2(Path(__file__).with_name('theis_reference.py'), data/'theis_reference.py')
    np.savetxt(data/'balance.csv', budget, delimiter=',',
               header='time_s,stored_change_m3,expected_change_m3', comments='')

    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        line = (np.abs(coord[:, 1]) < 1e-9) & (r > 0) & (r <= ref['L']/2) & (np.abs(coord[:, 2]+ref['B']/2) < ref['B']/2+1e-9)
        grid = np.linspace(.5, ref['L']/2, 300)
        colors = ['#1f4e79', '#5b8cbf', '#b8432a', '#7a1f12']
        ax = axes[0]
        picks = [1, 3, 5, 7]
        for k, i in enumerate(picks):
            t, s = history[i]
            ax.plot(grid, scheduled_drawdown(grid, np.full_like(grid, t), schedule, T, S)/scale, '-', color=colors[k], lw=1.2,
                    label=f'$t={t:g}$ s' + (' (recovery)' if t > ref['pump_until_s'] else ''))
            order = np.argsort(r[line]); keep = order[::max(1, len(order)//20)]   # about 20 nodes per profile
            ax.plot(r[line][keep], s[line][keep]/scale, 'o', ms=3.5, mfc='white', mec=colors[k], mew=.9)
        ds.band(ax, 10, ref['L']/2, y=.55)
        ax.set(xlabel='$r$ (m)', ylabel=r'$s_p/C_Q$', title='(a) Radial drawdown profiles', xlim=(0, ref['L']/2))
        ds.key(ax)
        ax = axes[1]
        tt = np.linspace(0, ref['end_s'], 400)
        for k, radius in enumerate((10.0, 20.0)):
            node = np.argmin(np.abs(r-radius) + 1e3*~line)
            ax.plot(tt, scheduled_drawdown(np.full_like(tt, r[node]), tt, schedule, T, S)/scale, '-', color=colors[2*k], lw=1.2,
                    label=f'$r={r[node]:.1f}$ m')
            ax.plot([h[0] for h in history], [h[1][node]/scale for h in history], 'o', ms=4, mfc='white', mec=colors[2*k], mew=1)
        ax.axvline(ref['pump_until_s'], color='#9aa5ab', lw=.8, ls=':')
        ax.set(xlabel='$t$ (s)', ylabel=r'$s_p/C_Q$', title='(b) Pumping, then recovery')
        ds.key(ax)
        ax = axes[2]
        h = np.array([c['spacing_m'] for c in summary['cases']])
        mx = np.array([c['max_error_over_scale'] for c in summary['cases']])
        rms = np.array([c['rms_error_over_scale'] for c in summary['cases']])
        ax.loglog(h, mx, 'o-', color=ds.HOT, mfc='white', label='max error')
        ax.loglog(h, rms, 's-', color=ds.WATER, mfc='white', label='RMS error')
        ax.loglog(h, mx[0]*(h/h[0])**2, '--', color=ds.MUTED, lw=.8, label='slope 2')
        ax.set(xlabel='mesh spacing $h$ (m)', ylabel=r'error / $C_Q$', title='(c) Mesh convergence')
        ax.set_xticks(h, [f'{x:g}' for x in h]); ax.xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ax.legend()
        pf.save_figure(fig, output/'figures', KEY)

    order = float(np.polyfit(np.log(h), np.log(mx), 1)[0])
    worst = finest['max_error_over_scale']
    return dict(key=KEY, theme='groundwater', title='Theis pumping and recovery in a confined aquifer',
        problem='A well pumps a confined aquifer and then shuts in. The drawdown spreads radially and recovers.',
        result=(f'Largest drawdown error: {100*worst:.3f}% of $C_Q$, over $10\\le r\\le30$ m and eight times through pumping and recovery.'),
        scope=(f"3D quarter box, {finest['spacing_m']:g} m mesh, {ref['end_s']/ref['dt_s']:g} steps of {ref['dt_s']:g} s. "
               f"On meshes of {', '.join(f'{x:g}' for x in h)} m the error falls with order {order:.2f}. "
               f"(a) radial profiles; (b) histories at two radii; (c) mesh convergence."),
        meaning=f'The closed-box water budget differs from the prescribed pumping integral by at most {budget_error:.2e} m³. This checks the rate-controlled well, hydraulic storage and radial diffusion, including shut-in superposition. Mechanical pressure feedback is disabled in this flow-isolation test.',
        metric=worst, metric_label='drawdown error / C_Q', dimension=3, spacing_m=finest['spacing_m'],
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 mesh_order=order, maximum_water_budget_error_m3=budget_error, cases=[{k: v for k, v in c.items() if k != 'frames'} for c in summary['cases']]), indent=1))
