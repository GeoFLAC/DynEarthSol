"""Validate the recorded DES Hantush-Jacob study and plot its finest mesh.

Interface as qualify_theis: build(study, output, record_source) recomputes
errors from the raw pressure fields, writes figures/hantush.* and
data/hantush/, and returns the review record.
"""
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np

import prepare_figures as pf
import diagram_style as ds
from hantush_reference import scheduled_drawdown
import theis_reference
sys.path.insert(0, str(pf.ROOT/'tests/functional/hydromechanics'))
import check_hantush as ch

KEY = 'hantush'


def build(study, output, record_source):
    summary = json.loads((study/'summary.json').read_text())
    assert hashlib.sha256((study/'runner.py').read_bytes()).hexdigest() == summary['runner_sha256']
    record_source(study/'summary.json'); record_source(study/'runner.py')
    ref = summary['reference']
    T, S, scale, lam = ref['transmissivity'], ref['storativity'], ref['pressure_scale_Pa'], ref['leakage_factor']
    schedule = [(0.0, ref['Q']), (ref['pump_until_s'], 0.0)]
    finest = min(summary['cases'], key=lambda c: c['spacing_m'])
    h = finest['spacing_m']
    folder = study/f'h{h:g}'
    for path in folder.glob('hantush.*'):
        if path.suffix in ('.cfg', '.info') or '.save.' in path.name: record_source(path)
    model, frames, _ = ch.compare(folder, h)
    for a, b in zip(frames, finest['frames']):
        np.testing.assert_allclose(a['max_error_over_scale'], b['max_error_over_scale'], rtol=0, atol=1e-12)
    rows = ch.scan_frames(str(folder/'hantush'))
    coord = model.read_field(0, 'coordinate'); p0 = model.read_field(0, 'pore pressure')
    r = np.hypot(coord[:, 0], coord[:, 1])
    history = [(row['time'], p0 - model.read_field(row['frame'], 'pore pressure')) for row in rows[1:]]

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    shutil.copy2(folder/'hantush.cfg', data/'input.cfg'); shutil.copy2(study/'summary.json', data/'summary.json')
    shutil.copy2(study/'runner.py', data/'check_hantush.py')
    shutil.copy2(Path(__file__).with_name('hantush_reference.py'), data/'hantush_reference.py')

    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        # Aquifer base, on the x axis: the line where the screen meets the sealed floor.
        line = (np.abs(coord[:, 1]) < 1e-9) & (r > 0) & (r <= 30) & (np.abs(coord[:, 2] + ref['B'] + h) < 1e-9)
        grid = np.linspace(.5, 30, 300)
        ax = axes[0]
        pump = [i for i, (t, _) in enumerate(history) if t <= ref['pump_until_s']]
        picks = sorted(set(pump[k*(len(pump)-1)//3] for k in range(4)))
        for color, i in zip(ds.SERIES, picks):
            t, s = history[i]
            ds.analytical(ax, grid, scheduled_drawdown(grid, np.full_like(grid, t), schedule, T, S, lam)/scale,
                          color=color, label=f'$t={t:g}$ s')
            order = np.argsort(r[line]); keep = order[::max(1, len(order)//16)]   # about 16 nodes per profile
            ax.plot(r[line][keep], s[line][keep]/scale, 'o', ms=3.5, mfc='white', mec=color, mew=.9)
        t = history[picks[-1]][0]
        ax.plot(grid, theis_reference.scheduled_drawdown(grid, np.full_like(grid, t), schedule, T, S)/scale,
                '--', color=ds.MUTED, lw=1, label=f'Theis, $t={t:g}$ s')
        ds.band(ax, *ref['window'], y=.55)
        ds.panel(ax, 'a', 'Drawdown while pumping', xlabel='$r$ (m)', ylabel=r'$s_p/C_Q$')
        ds.key(ax)
        ax.set_xlim(0, 30)
        ax = axes[1]
        tt = np.linspace(0, ref['end_s'], 400)
        for color, radius in zip((ds.SERIES[0], ds.SERIES[2]), (10.0, 20.0)):
            node = np.argmin(np.abs(r - radius) + 1e3*~line)
            ds.analytical(ax, tt, scheduled_drawdown(np.full_like(tt, r[node]), tt, schedule, T, S, lam)/scale,
                          color=color, label=f'$r={r[node]:.1f}$ m')
            ax.plot(tt, theis_reference.scheduled_drawdown(np.full_like(tt, r[node]), tt, schedule, T, S)/scale,
                    '--', color=color, lw=.8, alpha=.6)
            ax.plot([x[0] for x in history], [x[1][node]/scale for x in history], 'o', ms=4, mfc='white', mec=color, mew=1)
        ax.axvline(ref['pump_until_s'], color='#9aa5ab', lw=.8, ls=':')
        ax.plot([], [], '--', color=ds.MUTED, lw=.8, label='Theis, no leakage')
        ds.panel(ax, 'b', 'Leakage limits the drawdown', xlabel='$t$ (s)', ylabel=r'$s_p/C_Q$')
        ds.key(ax)
        ax = axes[2]
        hs = np.array([c['spacing_m'] for c in summary['cases']])
        mx = np.array([c['max_error_over_scale'] for c in summary['cases']])
        rms = np.array([c['rms_error_over_scale'] for c in summary['cases']])
        ax.loglog(hs, mx, 'o-', color=ds.HOT, mfc='white', label='max error')
        ax.loglog(hs, rms, 's-', color=ds.WATER, mfc='white', label='RMS error')
        ax.loglog(hs, mx[0]*(hs/hs[0])**2, '--', color=ds.MUTED, lw=.8, label='slope 2')
        ax.set_xticks(hs, [f'{x:g}' for x in hs]); ax.xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ds.panel(ax, 'c', 'Mesh convergence', xlabel='mesh spacing $h$ (m)', ylabel=r'error / $C_Q$')
        pf.save_figure(fig, output/'figures', KEY)

    order = float(np.polyfit(np.log(hs), np.log(mx), 1)[0])
    worst, leak = finest['max_error_over_scale'], finest['leakage_effect_over_scale']
    return dict(key=KEY, theme='groundwater', title='Hantush–Jacob pumping from a leaky aquifer',
        problem='A well pumps a confined aquifer covered by a leaky aquitard. Water drawn down through the aquitard limits the drawdown, which approaches a steady cone.',
        result=(f'Largest drawdown error: {100*worst:.3f}% of $C_Q$, over ${ref["window"][0]:g}\\le r\\le{ref["window"][1]:g}$ m '
                f'and eight times through pumping and recovery. Leakage changes the drawdown by up to {100*leak:.0f}% of $C_Q$ from Theis.'),
        scope=(f"3D quarter box, {h:g} m mesh, aquitard one element thick, {ref['end_s']/ref['dt_s']:g} steps of {ref['dt_s']:g} s; "
               f"$\\lambda={lam:g}$ m, vertical permeability {ref['anisotropy']:g} times horizontal. "
               f"On meshes of {', '.join(f'{x:g}' for x in hs)} m the error falls with order {order:.2f}. "
               f"(a) radial profiles, with Theis dashed; (b) histories at two radii, Theis dashed; (c) mesh convergence."),
        meaning=('This checks a layered aquifer system: an aquitard of a second material leaking from a fixed-pressure top, '
                 'with the well, storage and radial diffusion of the Theis case. On the finest mesh the aquitard stores '
                 f"{100*finest['aquitard_storage_ratio']:.2f}% of the aquifer's water, close to the no-storage assumption. "
                 'Mechanical pressure feedback is disabled in this flow-isolation test.'),
        metric=worst, metric_label='drawdown error / C_Q', dimension=3, spacing_m=h,
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 mesh_order=order, cases=[{k: v for k, v in c.items() if k != 'frames'} for c in summary['cases']]), indent=1))
