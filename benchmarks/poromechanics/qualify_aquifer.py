"""Validate the recorded thin-slab refinement and plot its finest mesh.

Interface as qualify_convection: build(study, output, record_source) with the
slab_refinement.py output. SPEC_UPDATE turns review.SPECS['averaged_aquifer']
from an analytical-only reference into the numerical comparison.
"""
import json
from pathlib import Path
import shutil

import numpy as np

import prepare_figures as pf
import diagram_style as ds
import pumping_slab as ps

KEY = 'averaged_aquifer'


def _binned(ax, r, v, color, marker='o', label=None):
    """Mean of every point in each log-radius bin as a marker, the spread as a faint band."""
    edges = np.geomspace(15, 120, 22); which = np.digitize(r, edges) - 1
    keep = [j for j in range(len(edges) - 1) if (which == j).any()]
    centre = np.array([np.exp(np.log(r[which == j]).mean()) for j in keep])
    mean = np.array([v[which == j].mean() for j in keep])
    lo = np.array([v[which == j].min() for j in keep]); hi = np.array([v[which == j].max() for j in keep])
    ax.fill_between(centre, lo, hi, color=color, alpha=.15, lw=0)
    ax.plot(centre, mean, marker, ms=3.5, mfc='white', mec=color, mew=.9, ls='none', label=label)


def build(study, output, record_source):
    summary = pf.load_study(study, record_source)
    cases = sorted(summary['cases'], key=lambda c: -c['resolution_m'])
    finest = cases[-1]
    case = study/f"res{finest['resolution_m']:g}"/'slab-3d'
    for path in case.glob('slab.*'):
        if path.suffix in ('.info',) or '.save.' in path.name: record_source(path)
    record_source(case/'input.cfg')
    aq, rows = ps.compare(case)
    assert len(rows) == len(finest['rows'])
    for a, b in zip(rows, finest['rows']):
        for key in ('pressure_rms', 'areal_rms', 'vertical_rms'):
            np.testing.assert_allclose(a[key], b[key], rtol=0, atol=1e-12)

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    shutil.copy2(case/'input.cfg', data/'input.cfg'); shutil.copy2(study/'summary.json', data/'summary.json')
    for name in ('slab_refinement.py', 'pumping_slab.py', 'averaged_aquifer_reference.py'):
        shutil.copy2(Path(__file__).with_name(name), data/name)

    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        r = np.geomspace(ps.BAND[0], ps.BAND[1], 200)
        ax = axes[0]
        for row, color in ((rows[0], '#5b8cbf'), (rows[-1], '#1f4e79')):
            ref = aq.line_well(r, row['time_s'], ps.Q)
            ax.plot(r, -ref['pressure']/1e3, '-', color=color, lw=1.2, label=f"$t={row['time_s']:.0f}$ s")
            _binned(ax, row['plot']['r_node'], -row['plot']['p']/1e3, color)
        ax.set(xscale='log', xlabel='$r$ (m)', ylabel='pressure drop (kPa)', title='(a) Pressure drawdown')
        radii = [15, 20, 30, 50, 80, 120]
        ax.set_xticks(radii, [str(x) for x in radii]); ax.xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ds.key(ax)
        ax = axes[1]
        last = rows[-1]; ref = aq.line_well(r, last['time_s'], ps.Q)
        ax.plot(r, -ref['areal_strain']*1e6, '-', color=ds.INK, lw=1.2, label='analytical (areal = vertical)')
        _binned(ax, last['plot']['r_elem'], -last['plot']['areal']*1e6, ds.WATER, 'o', 'DES areal (binned)')
        _binned(ax, last['plot']['r_elem'], -last['plot']['vertical']*1e6, ds.HOT, 's', 'DES vertical (binned)')
        ax.set(xscale='log', xlabel='$r$ (m)', ylabel='contraction (microstrain)', title=f"(b) Strain split at $t={last['time_s']:.0f}$ s")
        radii = [15, 20, 30, 50, 80, 120]
        ax.set_xticks(radii, [str(x) for x in radii]); ax.xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ax.legend()
        ax = axes[2]
        h = np.array([c['resolution_m'] for c in cases])
        for key, marker, color, name in (('pressure_rms', 'o', ds.WATER, 'pressure'), ('areal_rms', 's', ds.INK, 'areal strain'),
                                         ('vertical_rms', '^', ds.HOT, 'vertical strain')):
            ax.loglog(h, [c[key] for c in cases], marker+'-', color=color, mfc='white', label=name)
        ax.loglog(h, cases[0]['vertical_rms']*(h/h[0]), '--', color=ds.MUTED, lw=.8, label='slope 1')
        ax.set(xlabel='mesh spacing (m)', ylabel='RMS error (relative)', title='(c) Refinement')
        ax.set_xticks(h, [f'{x:g}' for x in h]); ax.xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ax.legend()
        pf.save_figure(fig, output/'figures', KEY)

    seq = '; '.join(f"{c['resolution_m']:g} m: {100*c['pressure_rms']:.2f}%, {100*c['areal_rms']:.2f}%, {100*c['vertical_rms']:.2f}%"
                    for c in cases)
    return dict(key=KEY, theme='consolidation', title='Pumping subsidence: averaged-aquifer (Bear–Corapcioglu) model',
        problem='A well pumps a thin poroelastic aquifer. The drawdown spreads radially and the layer contracts, splitting its volume change between areal and vertical strain.',
        result=(f"RMS errors: {100*finest['pressure_rms']:.2f}% in pressure, {100*finest['areal_rms']:.2f}% in areal strain and "
                f"{100*finest['vertical_rms']:.2f}% in vertical strain, relative to the largest analytical values."),
        scope=(f"3D quarter slab, {finest['resolution_m']:g} m mesh, four times while the front stays inside. RMS errors by mesh "
               f"(pressure, areal, vertical): {seq}. (a) pressure drawdown; (b) strain split; (c) mesh convergence."),
        meaning='This checks two-way coupled pumping deformation against the thickness-averaged model. It does not remove the $O((B/r)^2)$ slab-model difference, which remains at the finest mesh.',
        metric=finest['vertical_rms'], metric_label='vertical strain RMS error', dimension=3, spacing_m=finest['resolution_m'],
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'],
                                 cases=[{k: v for k, v in c.items() if k != 'rows'} for c in cases]), indent=1))
