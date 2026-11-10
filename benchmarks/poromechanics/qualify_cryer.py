"""Validate the recorded DES Cryer studies and plot the selected run.

Interface as qualify_convection. build(study, output, record_source, space, times)
takes the selected run, the mesh sequence and the time-step sequence. SPEC is
the review entry for 'cryer'; setup(folder) draws its setup.
"""
import json
from pathlib import Path
import shutil
import sys

import numpy as np

import prepare_figures as pf
import diagram_style as ds
sys.path.insert(0, str(pf.ROOT/'tests/functional/hydromechanics'))
import check_cryer as cc

KEY = 'cryer'


def _late(case):
    """Largest error and RMS after ct/R^2 = 0.05, past the one-step coupling lag."""
    late = [f for f in case['frames'] if f['t_over_tc'] > .05]
    return max(f['max_error_over_p0'] for f in late), max(f['rms_error_over_p0'] for f in late)


def build(study, output, record_source, space, times=()):
    # Every Cryer study beside `study` joins the mesh x time-step table, so
    # runs added later need no change to the caller.
    studies = sorted({study, space, *times, *(d for d in study.parent.glob('cryer-*') if (d/'summary.json').exists())})
    table = []
    for d in studies:
        sm = pf.load_study(d, record_source)
        for c in sm['cases']:
            table.append(dict(study=d, steps=sm['steps'], summary=sm, case=c, late=_late(c)))
    best = min(table, key=lambda e: (e['case']['spacing_m'], -e['steps']))
    study, summary, case = best['study'], best['summary'], best['case']
    folder = study/f"h{case['spacing_m']:g}"
    for path in folder.glob('cryer.*'):
        if path.suffix in ('.cfg', '.info') or '.save.' in path.name: record_source(path)
    record_source(folder/'octant.poly')
    frames, nodes, _ = cc.compare(folder)
    assert len(frames) == len(case['frames'])
    for a, b in zip(frames, case['frames']):
        for key in ('max_error_over_p0', 'rms_error_over_p0', 'centre_des_over_p0'):
            np.testing.assert_allclose(a[key], b[key], rtol=0, atol=1e-12)
    tc, p0 = cc.TC, cc.REF.p0

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    for name in ('cryer.cfg', 'octant.poly'):
        shutil.copy2(folder/name, data/('input.cfg' if name.endswith('cfg') else name))
    shutil.copy2(study/'summary.json', data/'summary.json'); shutil.copy2(study/'runner.py', data/'check_cryer.py')
    shutil.copy2(Path(__file__).with_name('cryer_reference.py'), data/'cryer_reference.py')

    model = cc.Dynearthsol(str(folder/'cryer')); rows = cc.scan_frames(str(folder/'cryer'))
    coord = model.read_field(0, 'coordinate'); r = np.linalg.norm(coord, axis=1)
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        ax = axes[0]
        tt = np.linspace(1e-4, .5, 400)*tc
        ax.plot(tt/tc, cc.REF.pressure(np.zeros_like(tt), tt)/p0, '-', color=ds.HOT, lw=1.2, label='analytical')
        ax.plot([f['t_over_tc'] for f in frames], [f['centre_des_over_p0'] for f in frames], 'o', ms=3.5,
                mfc='white', mec=ds.HOT, mew=.9, label='DES')
        ax.axhline(1, color='#9aa5ab', lw=.8, ls=':')
        ax.set(xlabel=r'$ct/R^2$', ylabel='$p/p_0$', title='(a) Centre pressure')
        ax.legend()
        ax = axes[1]
        rs = np.linspace(0, 1, 200)
        colors = ['#7a1f12', '#b8432a', '#5b8cbf', '#1f4e79']
        for k, target in enumerate((.02, .06, .15, .35)):
            i = int(np.argmin([abs(f['t_over_tc']-target) for f in frames]))
            t = frames[i]['time']; p = model.read_field(rows[i+1]['frame'], 'pore pressure')
            ax.plot(rs, cc.REF.pressure(rs, np.full_like(rs, t))/p0, '-', color=colors[k], lw=1.2, label=f'$ct/R^2={t/tc:.2f}$')
            # Every node, binned by radius: the bin mean as a marker, the spread as a faint band.
            edges = np.linspace(0, 1, 23); which = np.clip(np.digitize(r, edges) - 1, 0, 21)
            centre = np.array([r[which == j].mean() for j in range(22) if (which == j).any()])
            stats = [(p[which == j].mean(), p[which == j].min(), p[which == j].max()) for j in range(22) if (which == j).any()]
            mean, lo, hi = (np.array(x)/p0 for x in zip(*stats))
            ax.fill_between(centre, lo, hi, color=colors[k], alpha=.15, lw=0)
            ax.plot(centre, mean, 'o', ms=3.5, mfc='white', mec=colors[k], mew=.9)
        ax.set(xlabel='$r/R$', ylabel='$p/p_0$', title='(b) Radial profiles, nodes binned by radius')
        ds.key(ax, loc='lower left')
        old = ax.get_legend()
        ax.set_ylim(-.42, None)                          # a clear strip below p = 0 for the key
        ax.set_yticks(np.arange(0, 1.51, .25))
        ax.legend(old.legend_handles, [t.get_text() for t in old.texts], ncol=3, fontsize=10, loc='lower center')
        ax = axes[2]
        counts = sorted({e['steps'] for e in table})
        for k, n in enumerate(counts):
            row = sorted((e['case']['spacing_m'], e['late'][0]) for e in table if e['steps'] == n)
            ax.loglog([x for x, _ in row], [y for _, y in row], 'o-', color=ds.SERIES[k % len(ds.SERIES)], mfc='white', label=f'{n} steps')
        h = np.array(sorted({e['case']['spacing_m'] for e in table}))
        ax.set_xticks(h, [f'{x:g}' for x in h]); ax.xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ds.panel(ax, 'c', 'Error by mesh and step count', xlabel='mesh spacing $h$ (m)', ylabel=r'max error after $ct/R^2=0.05$ / $p_0$')
        pf.save_figure(fig, output/'figures', KEY)

    first = frames[0]
    late, rms = best['late']
    cells = '; '.join(f"{h_:g} m: " + ', '.join(f"{e['steps']} steps {100*e['late'][0]:.2f}%" for e in sorted(
        (e for e in table if e['case']['spacing_m'] == h_), key=lambda e: e['steps'])) for h_ in sorted({e['case']['spacing_m'] for e in table}, reverse=True))
    return dict(key=KEY, theme='consolidation', title='Cryer consolidation of a sphere',
        problem='A drained poroelastic sphere is suddenly loaded on its surface. The shell drains and contracts first, squeezing the core, so the centre pressure rises before it decays.',
        result=(f"Largest pressure error after $ct/R^2=0.05$: {100*late:.2f}% of $p_0$, and {100*case['max_error_over_p0']:.2f}% including the first-step lag. "
                f"Sampled centre peak: DES ${case['centre_peak_des_over_p0']:.4f}\\,p_0$, analytical ${case['centre_peak_ref_over_p0']:.4f}\\,p_0$."),
        scope=(f"3D octant, {case['spacing_m']:g} m mesh, {summary['steps']} PT steps over $0.5\\,R^2/c$. Largest error after $ct/R^2=0.05$ by mesh and step count: "
               f"{cells}. A finer mesh lowers the centre peak and a longer step raises it, so their errors partly cancel and the two must be refined together. "
               "(a) centre pressure; (b) radial profiles; (c) error by mesh and step count."),
        meaning='This checks 3D two-way coupling under a sudden load, the Mandel–Cryer pressure and the normal-stress load on a curved surface. It does not remove the one-step lag of the sequential coupling in the earliest frames, which falls with the time step.',
        metric=float(late), metric_label='pressure error / p0 after ct/R^2 = 0.05', dimension=3, spacing_m=case['spacing_m'],
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 table=[dict(spacing_m=e['case']['spacing_m'], steps=e['steps'], late_max=e['late'][0], late_rms=e['late'][1],
                                             centre_peak=e['case']['centre_peak_des_over_p0']) for e in table]), indent=1))
