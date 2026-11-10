"""Validate the recorded DES Mandel studies and plot the selected run.

Interface as qualify_convection. build(study, output, record_source, space, times)
takes the selected run (finest time step), the mesh sequence and the time-step
sequence. SPEC is the review entry for 'mandel'; setup(folder) draws its setup.
"""
import json
from pathlib import Path
import shutil
import sys

import numpy as np

import prepare_figures as pf
import diagram_style as ds
sys.path.insert(0, str(pf.ROOT/'tests/functional/hydromechanics'))
import check_mandel as cm

KEY = 'mandel'


def build(study, output, record_source, space, times):
    summary = pf.load_study(study, record_source)
    case = summary['cases'][0]
    folder = study/f"h{case['spacing_m']:g}"
    for path in folder.glob('mandel.*'):
        if path.suffix in ('.cfg', '.info') or '.save.' in path.name: record_source(path)
    _, rows, frames = cm.compare(folder)
    assert len(frames) == len(case['frames'])
    for a, b in zip(frames, case['frames']):
        for key in ('max_error_over_p0', 'centre_des_over_p0', 'settlement_m', 'platen_vz_spread'):
            np.testing.assert_allclose(a[key], b[key], rtol=0, atol=1e-12)
    mesh = pf.load_study(space, record_source)
    steps = [(pf.load_study(t, record_source), t) for t in times]
    tc, p0 = cm.TC, cm.REF.p0

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    shutil.copy2(folder/'mandel.cfg', data/'input.cfg'); shutil.copy2(study/'summary.json', data/'summary.json')
    shutil.copy2(study/'runner.py', data/'check_mandel.py')
    shutil.copy2(Path(__file__).with_name('mandel_reference.py'), data/'mandel_reference.py')

    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        ax = axes[0]
        tt = np.linspace(1e-4, 1, 400)*tc
        for x, color, name in ((0.0, ds.HOT, 'centre $x=0$'), (.5, ds.WATER, '$x=a/2$')):
            ax.plot(tt/tc, cm.REF.pressure(np.full_like(tt, x), tt)/p0, '-', color=color, lw=1.2, label=name)
        model = cm.Dynearthsol(str(folder/'mandel'))
        x = model.read_field(0, 'coordinate')[:, 0]
        mid = np.argmin(np.abs(x-.5))
        ax.plot([f['t_over_tc'] for f in frames], [f['centre_des_over_p0'] for f in frames], 'o', ms=3.5,
                mfc='white', mec=ds.HOT, mew=.9)
        ax.plot([r['time']/tc for r in rows[1:]], [model.read_field(r['frame'], 'pore pressure')[mid]/p0 for r in rows[1:]],
                'o', ms=3.5, mfc='white', mec=ds.WATER, mew=.9)
        ax.axhline(1, color='#9aa5ab', lw=.8, ls=':')
        ax.set(xlabel=r'$ct/a^2$', ylabel='$p/p_0$', title='(a) Pressure histories')
        ds.key(ax)
        ax = axes[1]
        xs = np.linspace(0, 1, 200)
        # One node per distinct x (pressure is uniform in z), about 25 per profile.
        _, first = np.unique(np.round(x, 9), return_index=True)
        pick = first[::max(1, len(first)//25)]
        colors = ['#7a1f12', '#b8432a', '#5b8cbf', '#1f4e79']
        for k, target in enumerate((.02, .1, .3, .7)):
            i = int(np.argmin([abs(f['t_over_tc']-target) for f in frames]))
            row, t = rows[i+1], frames[i]['time']
            p = model.read_field(row['frame'], 'pore pressure')
            ax.plot(xs, cm.REF.pressure(xs, np.full_like(xs, t))/p0, '-', color=colors[k], lw=1.2, label=f'$ct/a^2={t/tc:.2f}$')
            ax.plot(x[pick], p[pick]/p0, 'o', ms=3.5, mfc='white', mec=colors[k], mew=.9)
        ax.set(xlabel='$x/a$', ylabel='$p/p_0$', title='(b) Profiles across the strip')
        ds.key(ax)
        old = ax.get_legend()                            # two columns keep the key below the profiles
        ax.legend(old.legend_handles, [t.get_text() for t in old.texts], ncol=2, fontsize=10.5, loc='lower left')
        ax = axes[2]
        seq = sorted([(s['steps'], s['cases'][0]['max_error_over_p0'], s['cases'][0]['max_settlement_error']) for s, _ in steps]
                     + [(summary['steps'], case['max_error_over_p0'], case['max_settlement_error'])])
        mesh_h = [c for c in mesh['cases'] if abs(c['spacing_m']-case['spacing_m']) < 1e-12]
        if mesh_h:
            seq = sorted(set(seq + [(mesh['steps'], mesh_h[0]['max_error_over_p0'], mesh_h[0]['max_settlement_error'])]))
        dt = np.array([1/n for n, _, _ in seq]); pe = np.array([e for _, e, _ in seq]); se = np.array([e for _, _, e in seq])
        ax.loglog(dt, pe, 'o-', color=ds.HOT, mfc='white', label='pressure, max / $p_0$')
        ax.loglog(dt, se, 's-', color=ds.WATER, mfc='white', label='settlement, max relative')
        ax.loglog(dt, pe[-1]*dt/dt[-1], '--', color=ds.MUTED, lw=.8, label='slope 1')
        ax.set(xlabel=r'$\Delta t\,c/a^2$', ylabel='error', title='(c) Time-step convergence')
        ax.set_xticks(dt, [f'1/{n}' for n, _, _ in seq]); ax.xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ax.legend()
        pf.save_figure(fig, output/'figures', KEY)

    time_order = float(np.polyfit(np.log(dt), np.log(pe), 1)[0])
    peak_des, peak_ref = case['centre_peak_des_over_p0'], case['centre_peak_ref_over_p0']
    spacings = ', '.join(f"{c['spacing_m']:g}" for c in mesh['cases'])
    mesh_errors = ', '.join(f"{100*c['max_error_over_p0']:.2f}%" for c in mesh['cases'])
    return dict(key=KEY, theme='consolidation', title='Mandel consolidation under a rigid plate',
        problem='Rigid frictionless plates squeeze a drained poroelastic strip. The edges drain first, load shifts to the undrained core, and the centre pressure rises before it decays.',
        result=(f'Largest pressure error: {100*case["max_error_over_p0"]:.2f}% of $p_0$, over all nodes and times; largest relative settlement error: '
                f'{100*case["max_settlement_error"]:.2f}%. Sampled centre peak: DES ${peak_des:.4f}\\,p_0$, analytical ${peak_ref:.4f}\\,p_0$.'),
        scope=(f"2D, {case['spacing_m']:g} m mesh, {summary['steps']} PT steps over $a^2/c$. Halving the time step reduces the error with order "
               f"{time_order:.2f}; at {mesh['steps']} steps, meshes of {spacings} m give {mesh_errors}, so the time step dominates. "
               f"(a) pressure histories; (b) profiles; (c) time-step convergence."),
        meaning='This checks two-way poroelastic coupling under a sudden load, the non-monotonic Mandel–Cryer pressure and the rigid-platen boundary. It does not test large strain or plasticity.',
        metric=case['max_error_over_p0'], metric_label='pressure error / p0', dimension=2, spacing_m=case['spacing_m'],
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 time_order=time_order,
                                 time_sequence=[dict(steps=n, pressure=e, settlement=s) for n, e, s in seq],
                                 mesh_sequence=[{k: v for k, v in c.items() if k != 'frames'} for c in mesh['cases']]), indent=1))
