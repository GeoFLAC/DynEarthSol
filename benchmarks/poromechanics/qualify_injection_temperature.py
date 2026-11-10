"""Validate the recorded injection-temperature studies and plot them.

Interface as qualify_theis: build(avdonin, output, record_source, coarse, steps, doublet)
replays the Avdonin comparison of the finest mesh and the doublet heat
ledger from raw saves, writes figures/injection_temperature.* and
data/injection_temperature/, and returns the review record. avdonin and
coarse are the Avdonin mesh studies, steps a list of step studies on one
mesh, doublet the closed-box ledger study.
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
import check_injection_temperature as hot
import check_injection_energy as ce

KEY = 'injection_temperature'


def _load(study, record_source):
    summary = json.loads((study/'summary.json').read_text())
    assert hashlib.sha256((study/'runner.py').read_bytes()).hexdigest() == summary['runner_sha256']
    record_source(study/'summary.json'); record_source(study/'runner.py')
    return summary


def _record_run(folder, record_source):
    for path in folder.glob('hot.*'):
        if path.suffix in ('.cfg', '.info') or '.save.' in path.name: record_source(path)


def _ledger_series(folder):
    """Stored heat gain and the independent ledger at every saved step of a doublet run."""
    model = ce.Dynearthsol(str(folder/'hot')); rows = ce.scan_frames(str(folder/'hot'))
    coord = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    area = hot.lumped_area(coord, conn)
    nodes, w = ce.shape_weights(coord, conn, ce.WELLS[1])
    t, stored, produced = [], [], []
    for row in rows:
        T = model.read_field(row['frame'], 'temperature')
        t.append(row['time']); stored.append(hot.CAPACITY*float(area @ T)); produced.append(float(w @ T[nodes]))
    t, stored, produced = map(np.array, (t, stored, produced))
    net = np.concatenate([[0.0], np.cumsum(np.diff(t)*.5*(produced[1:] + produced[:-1]))])
    return t, stored - stored[0], hot.FLUID*ce.Q*(hot.TIN*t - net), produced


def build(avdonin, output, record_source, coarse, steps, doublet, sanity=()):
    mesh, low, energy = (_load(s, record_source) for s in (avdonin, coarse, doublet))
    slabs = []
    for slab in sanity:
        saved = _load(slab, record_source)
        for case in saved['cases']:
            folder = slab/f"h{case['spacing_m']:g}"
            _record_run(folder, record_source)
            replay, _ = hot.compare(folder)
            for name in ('temperature_error', 'heat_error'):
                np.testing.assert_allclose([r[name] for r in replay],
                                           [r[name] for r in case['frames']], rtol=1e-8, atol=1e-12)
            slabs.append(case)
    slabs.sort(key=lambda c: -c['spacing_m'])
    timing = [_load(s, record_source) for s in steps]
    cases = sorted(mesh['cases'] + low['cases'], key=lambda c: -c['spacing_m'])
    finest = cases[-1]
    folder = avdonin/f"h{finest['spacing_m']:g}"; _record_run(folder, record_source)
    again, _ = hot.compare(folder)
    np.testing.assert_allclose([f['temperature_error'] for f in again],
                               [f['temperature_error'] for f in finest['frames']], rtol=1e-10)
    pairs = sorted(energy['cases'], key=lambda c: -c['spacing_m'])
    series = {}
    for case in pairs:
        f = doublet/f"h{case['spacing_m']:g}"; _record_run(f, record_source)
        series[case['spacing_m']] = _ledger_series(f)
        t, gained, expected, _ = series[case['spacing_m']]
        drift = np.abs((gained - expected)[2:] - (gained - expected)[1])/(hot.FLUID*ce.Q*(hot.TIN - hot.T0)*t[2:])
        np.testing.assert_allclose(drift.max(), case['ledger_error'], rtol=1e-8)

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    shutil.copy2(folder/'hot.cfg', data/'input.cfg')
    shutil.copy2(doublet/f"h{pairs[-1]['spacing_m']:g}"/'hot.cfg', data/'input-doublet.cfg')
    for name, root in [('summary-avdonin.json', avdonin), ('summary-avdonin-coarse.json', coarse),
                       ('summary-doublet.json', doublet)] + [(f'summary-steps{s["steps"]}.json', p) for s, p in zip(timing, steps)]:
        shutil.copy2(root/'summary.json', data/name)
    shutil.copy2(avdonin/'runner.py', data/'check_injection_temperature.py')
    shutil.copy2(doublet/'runner.py', data/'check_injection_energy.py')
    shutil.copy2(Path(__file__).with_name('avdonin_reference.py'), data/'avdonin_reference.py')

    ref = hot.REF
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        ax = axes[0]
        model = hot.Dynearthsol(str(folder/'hot')); rows = hot.scan_frames(str(folder/'hot'))
        coord = model.read_field(0, 'coordinate')
        edge = np.isclose(coord[:, 1], 0)
        order = np.argsort(coord[edge, 0])
        r = np.linspace(.01, hot.L, 400)
        for color, target in zip(ds.SERIES, hot.TIMES):
            row = min(rows[1:], key=lambda x: abs(x['time'] - target))
            T = model.read_field(row['frame'], 'temperature')[edge][order]
            ds.analytical(ax, r, ref.temperature(r, target) - hot.T0, color=color, label=f'{target/86400:.1f} d')
            ds.des(ax, coord[edge, 0][order][::2], T[::2] - hot.T0, color=color, label=None, size=3)
        ax.set_xlim(0, 4.5)
        ds.panel(ax, 'a', 'Radial injection (Avdonin)', xlabel='$r$ (m)', ylabel=r'$T-T_0$ (K)')
        ds.key(ax, loc='upper right')

        ax = axes[1]
        t, gained, expected, produced = series[pairs[-1]['spacing_m']]
        days = t/86400
        ax.plot(days, expected/1e9, '-', color=ds.INK, lw=1.2, label='ledger: injected − extracted')
        k = max(1, len(t)//25)
        ds.des(ax, days[::k], gained[::k]/1e9, label='DES stored heat')
        ax2 = ax.twinx()
        ax2.plot(days, produced - hot.T0, '--', color=ds.WATER, lw=1.0, label='producer $T-T_0$')
        ax2.set_ylabel(r'producer $T-T_0$ (K)', color=ds.WATER); ax2.tick_params(axis='y', colors=ds.WATER)
        ax2.set_ylim(0, hot.TIN - hot.T0); ax2.grid(False)
        handles = ax.get_legend_handles_labels(); more = ax2.get_legend_handles_labels()
        ds.panel(ax, 'b', 'Closed-box heat ledger, hot doublet', xlabel='$t$ (d)', ylabel='heat gained (GJ per m)')
        ax.legend(handles[0] + more[0], handles[1] + more[1], loc='upper left', fontsize=8)

        ax = axes[2]
        h = np.array([c['spacing_m'] for c in cases])
        ax.loglog(h, [100*c['frames'][-1]['temperature_error'] for c in cases], 'o-', color=ds.HOT, mfc='white',
                  label='Avdonin, max')
        ax.loglog(h, [100*c['frames'][-1]['temperature_rms'] for c in cases], 's:', color=ds.HOT, mfc='white',
                  label='Avdonin, rms')
        ax.loglog(h, 100*cases[0]['frames'][-1]['temperature_error']*h/h[0], '--', color=ds.MUTED, lw=.8, label='slope 1')
        hd = np.array([c['spacing_m'] for c in pairs])
        ax.loglog(hd, [100*c['ledger_error'] for c in pairs], 'D-', color=ds.WATER, mfc='white', label='doublet ledger')
        ax.set_xticks(sorted(set(h) | set(hd)), [f'{x:g}' for x in sorted(set(h) | set(hd))])
        ax.xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ds.panel(ax, 'c', 'Errors by mesh', xlabel='mesh spacing $h$ (m)', ylabel=r'% of $T_{in}-T_0$ or of injected heat')
        ax.legend(fontsize=8)
        pf.save_figure(fig, output/'figures', KEY)

    last = lambda c, k: c['frames'][-1][k]
    return dict(key=KEY, theme='heat', title='Hot-water injection: Avdonin solution and heat ledger',
        problem=('Injected water carries its temperature. DES adds the injected enthalpy relative to the local temperature, '
                 'which with fluid mass conservation makes its advective heat equation conserve energy. Two checks: radial '
                 'hot-water injection against Avdonin’s solution, and a sealed doublet against an independent heat ledger.'),
        result=(f"Avdonin on the {finest['spacing_m']:g} m mesh: largest nodal error {100*last(finest, 'temperature_error'):.2f}% of "
                f"$T_{{in}}-T_0$ at {hot.TIMES[-1]/86400:.1f} d (rms {100*last(finest, 'temperature_rms'):.2f}%). Doublet: stored "
                f"heat follows the ledger within {100*pairs[-1]['ledger_error']:.3f}% after the first step, through breakthrough "
                f"(producer at {100*pairs[-1]['producer_rise']:.0f}% of $T_{{in}}-T_0$)."),
        scope=('2D. Avdonin errors fall at first order with the mesh ('
               + ', '.join(f"{c['spacing_m']:g} m: {100*last(c, 'temperature_error'):.2f}%" for c in cases)
               + '), as expected of the first-order (limited) advection scheme; time steps '
               + ', '.join(f"{s['steps']}: {100*s['cases'][0]['frames'][-1]['temperature_error']:.2f}%" for s in timing)
               + f" on the {timing[0]['cases'][0]['spacing_m']:g} m mesh. "
               'Thermal advection uses the Darcy flux at the start of each step, so a source switched on from rest loses '
               'part of its first step of heat when that step needs several thermal substeps: '
               + ', '.join(f"{c['start_deficit_steps']:.2f} of one step's injected heat on the {c['spacing_m']:g} m doublet" for c in pairs)
               + '; its share of the total falls as $1/t$. '
               + ('3D check, a slab one element thick: largest error '
                  + ', '.join(f"{100*last(c, 'temperature_error'):.2f}% ({c['spacing_m']:g} m)" for c in slabs)
                  + f", maximum stored-heat discrepancy {100*max(c['heat_error'] for c in slabs):.3f}%. "
                  'The finer slab retains a startup deficit; the percentage decreases with elapsed time. ' if slabs else ' ')
               + '(a) Avdonin profiles; (b) doublet ledger; (c) errors.'),
        meaning=('This checks the injected-enthalpy term with open outflow (Avdonin) and with outflow at a changing '
                 'temperature (doublet), measured against energy conservation rather than inferred from the equation.'),
        metric=last(finest, 'temperature_error'), metric_label='temperature error', dimension=2, spacing_m=finest['spacing_m'],
        evidence=json.dumps(dict(executable_sha256=mesh['executable_sha256'], source_commit=mesh['source_commit'],
                                 avdonin=[{k: v for k, v in c.items()} for c in cases],
                                 steps=[s['cases'] for s in timing], doublet=pairs, sanity_3d=slabs), indent=1))
