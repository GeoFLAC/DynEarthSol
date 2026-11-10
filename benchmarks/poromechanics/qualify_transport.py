#!/usr/bin/env python3
"""Draw the multiwell and open-boundary transport qualification results."""

import argparse
import json
from pathlib import Path
import shutil
import sys

import numpy as np

import diagram_style as ds
import prepare_figures as pf
from setup_drawings import schematic

sys.path.insert(0, str(pf.ROOT/'tests/functional/hydromechanics'))
import check_theis_wells as wells
import check_open_flow as pulse
import check_anisotropic_flow as anis


def build(data, output=None, record_source=None):
    record_source = record_source or (lambda path: None)
    qualification = json.loads((data/'qualification.json').read_text())
    assert qualification['numerical_checks_passed']
    record_source(data/'qualification.json')
    multi = json.loads((data/'multiwell/summary.json').read_text())
    opened = json.loads((data/'open-flow/summary.json').read_text())
    assert multi['passed'] and opened['passed']
    for name, summary in (('multiwell', multi), ('open-flow', opened)):
        path = data/name/'summary.json'
        assert wells.digest(path) == qualification['studies'][name]
        record_source(path)
        cases = summary['cases'].items() if name == 'multiwell' else (
            (f"h{c['spacing_m']:g}", c) for c in summary['cases'])
        for key, case in cases:
            folder = data/name/key
            cfg = folder/('input.cfg' if name == 'multiwell' else 'pulse.cfg')
            assert wells.digest(cfg) == case['input_sha256']
            record_source(cfg)
            for filename, expected in case['outputs_sha256'].items():
                path = folder/filename
                assert wells.digest(path) == expected
                record_source(path)

    def destination(key, slug):
        target = output/'figures' if output else data/'fig'/slug
        target.mkdir(parents=True, exist_ok=True)
        if output is None: schematic(key, target)
        return target, key if output else 'result'

    records = []
    fine = min((c for key, c in multi['cases'].items() if key.startswith('joint-')),
               key=lambda c: (c['spacing_m'], c['dt_s']))
    key = f"joint-h{fine['spacing_m']:g}-dt{fine['dt_s']:g}"
    model = pf.Dynearthsol(str(data/'multiwell'/key/'multiwell'))
    x = model.read_field(0, 'coordinate')
    rows = wells.ct.scan_frames(str(data/'multiwell'/key/'multiwell'))
    line = np.isclose(x[:, 1], 80) & np.isclose(x[:, 2], -wells.ct.B/2)
    radii = np.stack([np.hypot(x[:, 0]-wx, x[:, 1]-wy) for wx, wy, _ in wells.WELLS])
    line &= radii.min(axis=0) >= 10
    chosen = min(rows, key=lambda r: abs(r['time']-wells.PUMP))
    pressure = -model.read_field(chosen['frame'], 'pore pressure')
    grid = np.linspace(0, wells.L, 641)
    xy = np.column_stack((grid, np.full(len(grid), 80.)))
    excluded = np.min(np.stack([np.linalg.norm(xy-[wx, wy], axis=1)
                                for wx, wy, _ in wells.WELLS]), axis=0) < 10
    exact = np.full(len(grid), np.nan)
    exact[~excluded] = wells.reference(xy[~excluded], chosen['time'], wells.WELLS)
    target, figure_name = destination('theis_wells', 'theis-wells')
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        ds.analytical(axes[0], grid, exact/1000)
        ds.des(axes[0], x[line, 0], pressure[line]/1000)
        for centre in (60, 100):
            ds.band(axes[0], centre-10, centre+10, 'excluded near the well' if centre == 60 else '',
                    x=centre - 9, y=.5, rotation=90, va='center')
        ds.panel(axes[0], 'a', f"Drawdown at {chosen['time']:g} s", xlabel='$x$ (m), along $y=80$ m', ylabel='drawdown (kPa)')
        for dt, marker in ((1., 'o'), (.125, 's')):
            seq = sorted((c for name, c in multi['cases'].items() if name.startswith('joint-') and c['dt_s'] == dt),
                         key=lambda c: c['spacing_m'])
            axes[1].loglog([c['spacing_m'] for c in seq], [c['max_error_over_scale'] for c in seq],
                           marker+'-', mfc='white', color=ds.WATER if dt == 1 else ds.HOT, label=f'outer $dt={dt:g}$ s')
        ds.panel(axes[1], 'b', 'Mesh and timestep check', xlabel='mesh spacing (m)', ylabel=r'max pressure error / $C_\Sigma$')
        axes[1].set_xticks([2.5, 5], ['2.5', '5'])
        axes[1].xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        history = fine['frames']
        ds.analytical(axes[2], [0]+[r['time_s'] for r in history],
                      [0]+[-1000*r['exact_stored_m3'] for r in history], label='integrated well rate')
        ds.des(axes[2], [r['time_s'] for r in history], [-1000*r['stored_m3'] for r in history])
        ds.panel(axes[2], 'c', 'Closed fluid budget', xlabel='time (s)', ylabel='fluid deficit (L)')
        pf.save_figure(fig, target, figure_name)
    if output:
        packed = output/'data/theis_wells'; packed.mkdir(parents=True, exist_ok=True)
        shutil.copy2(data/'multiwell'/key/'input.cfg', packed/'input.cfg')
        shutil.copy2(data/'multiwell/summary.json', packed/'summary.json')
        for source in (Path(wells.__file__), Path(wells.ct.__file__), Path(__file__).with_name('theis_reference.py')):
            shutil.copy2(source, packed/source.name)
        records.append(dict(key='theis_wells', theme='groundwater', title='Several wells: pumping, injection and recovery',
            problem='Two wells extract water and one injects it into a confined slab; all stop after 40 s and the pressure recovers.',
            result=f"Maximum pressure error: {100*fine['max_error_over_scale']:.3f}% of $C_\\Sigma$; independent fluid-budget error below {max(c['max_volume_error_m3'] for c in multi['cases'].values()):.2e} m³.",
            scope=f"3D, 160 by 160 by 5 m, {fine['spacing_m']:g} m spacing and {fine['dt_s']:g} s outer step, through 80 s. Mesh and time-step studies show slight error cancellation. Separate single-well superposition agrees within {multi['max_discrete_superposition_error_over_scale']:.2e} of the pressure scale.",
            meaning='This checks linear saturated flow, signed pumping/injection, recovery and the closed volume budget. Mechanical feedback and nonlinear property changes are outside this comparison.',
            dimension=3, spacing_m=fine['spacing_m'], metric=fine['max_error_over_scale'], metric_label='pressure error / summed scale',
            evidence=json.dumps(multi)))

    fine = min(opened['cases'], key=lambda c: c['spacing_m'])
    folder = data/'open-flow'/f"h{fine['spacing_m']:g}"
    model = pf.Dynearthsol(str(folder/'pulse'))
    rows = pulse.ct.scan_frames(str(folder/'pulse'))
    x = model.read_field(0, 'coordinate')
    line = np.isclose(x[:, 1], .5) & np.isclose(x[:, 2], -.5)
    target, figure_name = destination('open_flow', 'open-flow')
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        grid = np.linspace(0, pulse.L, 401)
        for time, color in zip((0., 1., 5., 10.), ds.SERIES):
            row = min(rows, key=lambda r: abs(r['time']-time))
            exact = np.sin(np.pi*grid/pulse.L)*np.exp(-pulse.RATE*row['time'])
            ds.analytical(axes[0], grid, exact, color=color, label=f"$t={row['time']:g}$ s")
            ds.des(axes[0], x[line, 0], model.read_field(row['frame'], 'pore pressure')[line]/pulse.P,
                   color=color, label='_nolegend_')
        ds.panel(axes[0], 'a', 'Pressure pulse draining', xlabel='$x$ (m)', ylabel='$p/P$')
        ds.key(axes[0], loc='lower center')
        history = fine['history'][::20]+[fine['history'][-1]]
        t = np.linspace(0, pulse.END, 201)
        initial = 2*pulse.ct.STORAGE*pulse.P*pulse.WIDTH**2*pulse.L/np.pi
        for value, name, color, field, sign in ((np.exp(-pulse.RATE*t), 'stored', ds.WATER, 'stored_fluid_m3', 1),
                                                (1-np.exp(-pulse.RATE*t), 'outflow', ds.HOT, 'cumulative_boundary_m3', -1)):
            ds.analytical(axes[1], t, value, color=color, label=name)
            ds.des(axes[1], [r['time_s'] for r in history], [sign*r[field]/initial for r in history],
                   color=color, label='_nolegend_')
        ds.panel(axes[1], 'b', 'Stored fluid and outflow', xlabel='time (s)', ylabel='volume / initial analytical storage')
        ds.key(axes[1])
        seq = sorted(opened['cases'], key=lambda c: c['spacing_m'])
        for field, label, marker in (('max_pressure_error_over_amplitude', 'pressure / amplitude', 'o'),
                                     ('max_storage_error_over_initial', 'storage / initial', 's')):
            axes[2].loglog([c['spacing_m'] for c in seq], [c[field] for c in seq], marker+'-', mfc='white',
                           color=ds.WATER if marker == 'o' else ds.HOT, label=label)
        ds.panel(axes[2], 'c', 'Mesh convergence', xlabel='mesh spacing (m)', ylabel='maximum error (relative)')
        axes[2].set_xticks([.25, .5, 1], ['0.25', '0.5', '1'])
        axes[2].xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        pf.save_figure(fig, target, figure_name)
    if output:
        packed = output/'data/open_flow'; packed.mkdir(parents=True, exist_ok=True)
        shutil.copy2(folder/'pulse.cfg', packed/'input.cfg')
        shutil.copy2(data/'open-flow/summary.json', packed/'summary.json')
        for source in (Path(pulse.__file__), Path(pulse.ct.__file__)):
            shutil.copy2(source, packed/source.name)
        records.append(dict(key='open_flow', theme='groundwater', title='Open boundaries: a draining pressure pulse',
            problem='A smooth pressure pulse drains through two zero-pressure ends while the side faces remain sealed.',
            result=f"Maximum pressure error: {100*fine['max_pressure_error_over_amplitude']:.4f}% of the initial amplitude; storage error {100*fine['max_storage_error_over_initial']:.4f}%. Stored fluid plus outflow equals the initial volume within {fine['max_global_budget_error_m3']:.2e} m³ (independent integration).",
            scope=f"3D 10 by 1 by 1 m column, {fine['spacing_m']:g} m spacing, 0.05 s steps through 10 s. All 200 steps on each of three meshes enter the pressure, storage and boundary-budget checks.",
            meaning='This checks open-boundary drainage and the discrete global fluid ledger at constant capacity. Unsaturated coupling, changing storage and remapping still require separate budgets.',
            dimension=3, spacing_m=fine['spacing_m'], metric=fine['max_pressure_error_over_amplitude'], metric_label='pressure error / amplitude',
            evidence=json.dumps(opened)))
    return records


def build_anisotropic(data, output, record_source):
    """Replay the two-mode, two-phase flow check and select its finest run."""
    studies = [data/name for name in ('anisotropic-time', 'anisotropic-time-steps40k', 'anisotropic-time-steps80k')]
    summaries = []
    for study in studies:
        summary = json.loads((study/'summary.json').read_text())
        assert anis.digest(study/'runner.py') == summary['runner_sha256']
        record_source(study/'summary.json'); record_source(study/'runner.py')
        for case in summary['cases']:
            folder = study/f"h{case['spacing_m']:g}"
            assert anis.analyze(folder, case['spacing_m'], case['steps'], 2) == case
            for filename, digest in case['files_sha256'].items():
                path = folder/filename
                assert anis.digest(path) == digest; record_source(path)
        summaries.append(summary)
    mesh = sorted(summaries[0]['cases'], key=lambda c: -c['spacing_m'])
    fine = summaries[-1]['cases'][0]; h = fine['spacing_m']; folder = studies[-1]/f'h{h:g}'
    final = [pf.Dynearthsol(str(study/f'h{h:g}/second')).read_field(4, 'pore pressure') for study in studies]
    changes = [float(np.max(np.abs(a-b))) for a, b in zip(final, final[1:])]
    ratio = changes[0]/changes[1]
    orders = np.log([a['maximum_pressure_error']/b['maximum_pressure_error'] for a, b in zip(mesh, mesh[1:])])/np.log(2)
    assert np.all(orders > 1.8) and 1.8 < ratio < 2.2
    target = output/'figures'; target.mkdir(exist_ok=True)
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        model = pf.Dynearthsol(str(folder/'second')); x = model.read_field(2, 'coordinate')
        line = np.isclose(x[:, -1], -.5); order = np.argsort(x[line, 0]); xx = np.linspace(0, 1, 401)
        points = np.column_stack((xx, np.full_like(xx, -.5)))
        history = []
        for name in ('first', 'second'):
            m = pf.Dynearthsol(str(folder/name))
            for row in anis.scan_frames(str(folder/name)):
                if row['frame'] == 0 or (name == 'second' and row['frame'] == 2): continue
                history.append((m, row))
        for color, (m, row) in zip(ds.SERIES, history):
            exact = anis.basis(points)@anis.amplitudes(row['time'], fine['switch_time'])
            ds.analytical(axes[0], xx, exact, color=color, label=f"{row['time']:.2f} s")
            p = m.read_field(row['frame'], 'pore pressure')/anis.P0
            ds.des(axes[0], x[line, 0][order][::2], p[line][order][::2], color=color, label=None)
        ds.panel(axes[0], 'a', 'Pressure evolution', xlabel='$x$ (m)', ylabel='$p/P_0$ at mid-depth')
        ds.key(axes[0])
        tt = np.linspace(0, 1, 401); amplitudes = np.array([anis.amplitudes(t, fine['switch_time']) for t in tt])
        for n, color in enumerate((ds.SERIES[0], ds.SERIES[2])):
            ds.analytical(axes[1], tt, amplitudes[:, n], color=color, label=f'mode {n+1}')
            ds.des(axes[1], [r['time'] for r in fine['frames']], [r['modal_amplitudes'][n] for r in fine['frames']], color=color, label=None)
        axes[1].axvline(fine['switch_time'], color=ds.MUTED, ls=':', lw=.8)
        ds.panel(axes[1], 'b', 'Independent directional decay rates', xlabel='$t$ (s)', ylabel='mode amplitude / $P_0$')
        ds.key(axes[1])
        hs = np.array([c['spacing_m'] for c in mesh]); errors = np.array([c['maximum_pressure_error'] for c in mesh])
        axes[2].loglog(hs, 100*errors, 'o-', color=ds.WATER, mfc='white', label='20,000 steps')
        axes[2].loglog([h], [100*fine['maximum_pressure_error']], 'D', color=ds.HOT, label='80,000 steps')
        axes[2].loglog(hs, 100*errors[0]*(hs/hs[0])**2, '--', color=ds.MUTED, label='slope 2')
        axes[2].set_xticks(hs, [f'{v:g}' for v in hs]); axes[2].xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ds.panel(axes[2], 'c', 'Separate mesh and step studies', xlabel='spacing $h$ (m)', ylabel=r'max $|p-p_{ref}|/P_0$ (%)')
        pf.save_figure(fig, target, 'anisotropic_time')
    packed = output/'data/anisotropic_time'; packed.mkdir(parents=True, exist_ok=True)
    for name, source in [('input.cfg', folder/'first.cfg'), ('continuation.cfg', folder/'second.cfg'),
                          ('summary.json', studies[-1]/'summary.json'), ('summary-mesh.json', studies[0]/'summary.json'),
                          ('summary-steps.json', studies[1]/'summary.json')]: shutil.copy2(source, packed/name)
    replay = dict(spatial_orders=orders.tolist(), time_difference_Pa=changes, temporal_difference_ratio=ratio,
                  finest_pressure_error=fine['maximum_pressure_error'], numerical_checks_passed=True)
    (packed/'replay.json').write_text(json.dumps(replay, indent=2)+'\n')
    return dict(key='anisotropic_time', theme='groundwater', title='Anisotropic flow with changing permeability',
        problem='Two pressure modes drain from a square; scalar permeability and directional anisotropy change together at 0.5 s.',
        result=f"Maximum pressure error {100*fine['maximum_pressure_error']:.4f}% of P₀; spatial orders {orders[0]:.3f} and {orders[1]:.3f} between successive meshes. Halving the time step halves the change in the solution (ratio {ratio:.5f}): first order in time. A restart inside the second phase reproduces the uninterrupted run exactly.",
        scope=f'2D fixed square, h={h:g} m, physical step {1/fine["steps"]:g} s, two prescribed axis-aligned tensors. Storage uses the physical strain-form capacity, but volumetric strain does not feed back into the pressure equation. Pressure still loads mechanics; the physical step respects the elastic limit.',
        meaning='This checks directional diffusion, a prescribed coefficient jump and operator reconstruction at restart. Continuous state-dependent permeability, off-axis time changes and coupled fluid-content budgets remain separate requirements.',
        dimension=2, spacing_m=h, metric=fine['maximum_pressure_error'], metric_label='pressure error / initial scale', evidence=json.dumps(replay))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data', type=Path)
    build(parser.parse_args().data.resolve())
