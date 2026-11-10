"""Validate the recorded DES Rudnicki studies and plot the finest mesh.

Interface as qualify_convection. build(study, output, record_source, space, times)
takes the selected run, the mesh sequence (one study per spacing) and the
time-step check.
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
import check_rudnicki as cr

KEY = 'rudnicki'


def _check(study, record_source):
    summary = json.loads((study/'summary.json').read_text())
    assert hashlib.sha256((study/'runner.py').read_bytes()).hexdigest() == summary['runner_sha256']
    record_source(study/'summary.json'); record_source(study/'runner.py')
    return summary


def _binned(ax, r, v, color, edges):
    """Mean of every point in each radial bin as a marker, the spread as a faint band."""
    which = np.digitize(r, edges) - 1
    keep = [j for j in range(len(edges) - 1) if (which == j).any()]
    centre = np.array([r[which == j].mean() for j in keep])
    mean = np.array([v[which == j].mean() for j in keep])
    lo = np.array([v[which == j].min() for j in keep]); hi = np.array([v[which == j].max() for j in keep])
    ax.fill_between(centre, lo, hi, color=color, alpha=.15, lw=0)
    ax.plot(centre, mean, 'o', ms=3.5, mfc='white', mec=color, mew=.9)


def arrival_time():
    """When the diffusion front 2 sqrt(c t) reaches the inner edge of the band."""
    return cr.BAND[0]**2/(4*cr.REF.c)


def _arrived(case):
    """Errors once the front is in the band; earlier frames compare tiny pressures."""
    late = [f for f in case['frames'] if f['time'] >= arrival_time()]
    first = case['frames'][0]
    case.update({f'all_{k}': case[k] for k in ('pressure_rms', 'radial_rms', 'hoop_rms')},
                first_time=first['time'], first_pressure_rms=first['pressure_rms'],
                **{k: max(f[k] for f in late) for k in ('pressure_rms', 'radial_rms', 'hoop_rms')})
    return case


def build(study, output, record_source, space, times=()):
    summary = _check(study, record_source)
    case = summary['cases'][0]
    folder = study/f"h{case['spacing_m']:g}"
    for path in folder.glob('rudnicki.*'):
        if path.suffix in ('.cfg', '.info') or '.save.' in path.name: record_source(path)
    frames, nodes, _ = cr.compare(folder)
    assert len(frames) == len(case['frames'])
    for a, b in zip(frames, case['frames']):
        for key in ('pressure_rms', 'radial_rms', 'hoop_rms'):
            np.testing.assert_allclose(a[key], b[key], rtol=0, atol=1e-12)
    mesh = sorted((_check(s, record_source)['cases'][0] for s in space), key=lambda c: -c['spacing_m'])
    for c in mesh: _arrived(c)
    steps = [_check(t, record_source) for t in times]
    for t in steps: _arrived(t['cases'][0])
    _arrived(case)

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    shutil.copy2(folder/'rudnicki.cfg', data/'input.cfg'); shutil.copy2(study/'summary.json', data/'summary.json')
    shutil.copy2(study/'runner.py', data/'check_rudnicki.py')
    shutil.copy2(Path(__file__).with_name('rudnicki_reference.py'), data/'rudnicki_reference.py')

    model = cr.Dynearthsol(str(folder/'rudnicki')); rows = cr.scan_frames(str(folder/'rudnicki'))
    coord = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    rn = np.linalg.norm(coord, axis=1)
    rs = np.linspace(.5, 12, 300)
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        ax = axes[0]
        picks = [int(round(k*(len(rows)-2)/3)) + 1 for k in range(1, 4)] + [len(rows)-1]
        picks = sorted(set(picks))[-4:]
        for color, i in zip(ds.SERIES, picks):
            t = rows[i]['time']; p = model.read_field(rows[i]['frame'], 'pore pressure')
            ds.analytical(ax, rs, cr.REF.pressure(rs, t)/1e3, color=color, label=f'$t={t:g}$ s')
            _binned(ax, rn, p/1e3, color, np.linspace(1, 12, 34))
        ds.band(ax, *cr.BAND)
        ds.panel(ax, 'a', 'Pressure around the source', xlabel='$r$ (m)', ylabel='$p$ (kPa)')
        ds.key(ax)
        ax.set_ylim(0, 1.3*cr.REF.pressure(cr.BAND[0], rows[-1]['time'])/1e3)
        ax = axes[1]
        t = rows[-1]['time']
        re, err, ett = cr.polar_strains(coord, conn, model.read_field(rows[-1]['frame'], 'strain'))
        sel = (re > .5) & (re < 12)
        _, rr, tt = cr.REF.strains(rs, t)
        ds.analytical(ax, rs, 1e6*rr, color=ds.HOT, label=r'$\varepsilon_{rr}$')
        ds.analytical(ax, rs, 1e6*tt, color=ds.WATER, label=r'$\varepsilon_{\theta\theta}$')
        edges = np.linspace(1, 12, 34)
        _binned(ax, re[sel], 1e6*err[sel], ds.HOT, edges)
        _binned(ax, re[sel], 1e6*ett[sel], ds.WATER, edges)
        ax.set_xlim(0, 12)
        ax.set_ylim(1.4*min(1e6*np.min(rr[rs > 1]), -5), 1.15*1e6*np.max(tt[rs > 1]))   # the singular r < 1 m is off scale
        ds.band(ax, *cr.BAND); ax.axhline(0, color='#9aa5ab', lw=.8)
        ds.panel(ax, 'b', f'Strains at $t={t:g}$ s', xlabel='$r$ (m)', ylabel='strain ($10^{-6}$)')
        ds.key(ax)
        ax = axes[2]
        h = np.array([c['spacing_m'] for c in mesh])
        for name, marker, color, label in (('pressure_rms', 'o', ds.INK, 'pressure'), ('radial_rms', 's', ds.HOT, 'radial strain'),
                                           ('hoop_rms', '^', ds.WATER, 'hoop strain')):
            ax.loglog(h, [c[name] for c in mesh], marker+'-', color=color, mfc='white', label=label)
        ax.loglog(h, mesh[0]['pressure_rms']*(h/h[0]), '--', color=ds.MUTED, lw=.8, label='slope 1')
        ax.set_xticks(h, [f'{x:g}' for x in h]); ax.xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ds.panel(ax, 'c', 'Mesh convergence', xlabel='mesh spacing $h$ (m)', ylabel='RMS error after arrival (relative)')
        pf.save_figure(fig, output/'figures', KEY)

    order = float(np.polyfit(np.log(h), np.log([c['pressure_rms'] for c in mesh]), 1)[0]) if len(h) > 1 else float('nan')
    seq = '; '.join(f"{c['spacing_m']:g} m: {100*c['pressure_rms']:.2f}%, {100*c['radial_rms']:.2f}%, {100*c['hoop_rms']:.2f}%" for c in mesh)
    time_note = ''
    if steps:
        s = steps[0]; c = s['cases'][0]
        same = [m for m in mesh if abs(m['spacing_m']-c['spacing_m']) < 1e-12]
        if same:
            time_note = (f" Doubling the steps to {s['steps']} at {c['spacing_m']:g} m changes the pressure RMS error from "
                         f"{100*same[0]['pressure_rms']:.2f}% to {100*c['pressure_rms']:.2f}%.")
    return dict(key=KEY, theme='injection', title='Fluid injection into a poroelastic solid: Rudnicki',
        problem='Water is injected at a constant rate at a point. Pressure diffuses outward and the solid swells, changing the stress that a nearby fault would feel.',
        result=(f"Largest RMS errors once the front reaches the band ($t>{arrival_time():.0f}$ s): {100*case['pressure_rms']:.2f}% in pressure, "
                f"{100*case['radial_rms']:.2f}% in radial strain and {100*case['hoop_rms']:.2f}% in hoop strain, over "
                f"${cr.BAND[0]:g}\\le r\\le{cr.BAND[1]:g}$ m."),
        scope=(f"3D octant of a {cr.L:g} m box refined within {cr.ZONE:g} m, {case['spacing_m']:g} m mesh, {summary['steps']} PT steps to "
               f"{summary['duration_s']:g} s. RMS errors by mesh (pressure, radial, hoop): {seq}" + (f"; pressure order {order:.2f}." if len(h) > 2 else '.')
               + time_note + f" At $t={case['first_time']:g}$ s the front has not reached the band and the pressure RMS error is "
               f"{100*case['first_pressure_rms']:.1f}% of the much smaller pressures there. (a) pressure; (b) strains; (c) mesh convergence."),
        meaning='This checks two-way poroelastic coupling around a fluid source: the coupled diffusivity and the irrotational strain field that sets injection-induced stress on faults. It does not test the fault itself, frictional slip or the singular region next to the source.',
        metric=case['pressure_rms'], metric_label='pressure RMS error', dimension=3, spacing_m=case['spacing_m'],
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 mesh_sequence=[{k: v for k, v in c.items() if k != 'frames'} for c in mesh],
                                 time_check=[{k: v for k, v in s['cases'][0].items() if k != 'frames'} | dict(steps=s['steps']) for s in steps]), indent=1))
