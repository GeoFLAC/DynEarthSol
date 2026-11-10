"""Validate the recorded DES McTigue studies and plot the selected run.

Interface as qualify_theis. build(study, output, record_source, times, space)
takes the 2D mesh sequence, the finest mesh with twice the steps (the
selected run, since the mesh sequence reaches its time-step floor) and the
3D study.
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
import check_mctigue as cm

KEY = 'mctigue'
METRICS = (('pressure_error', 'o', ds.HOT, 'pressure'), ('temperature_error', 'D', ds.INK, 'temperature'),
           ('strain_error', '^', ds.WATER, 'vertical strain'), ('heave_error', 's', ds.MUTED, 'surface heave'))


def _check(study, record_source):
    summary = json.loads((study/'summary.json').read_text())
    assert hashlib.sha256((study/'runner.py').read_bytes()).hexdigest() == summary['runner_sha256']
    record_source(study/'summary.json'); record_source(study/'runner.py')
    return summary


def heat_ledger(folder):
    """Stored heat per unit area against the exact input through the heated surface, at each saved time.

    The base is insulated, so the exact temperature is an image series; the stored heat
    rho c_p INT (T - T0) dz is integrated from it. DES's lumped nodal volumes integrate its
    piecewise-linear temperature exactly.
    """
    from scipy.special import erfc
    capacity = cm.MAT['rho']*cm.MAT['cp']
    base = cm.Dynearthsol(str(folder/'start')); model = cm.Dynearthsol(str(folder/'mctigue'))
    coord = base.read_field(0, 'coordinate'); conn = base.read_field(0, 'connectivity')
    edges = coord[conn][:, 1:] - coord[conn][:, :1]
    size = np.abs(np.linalg.det(edges))/(2 if coord.shape[1] == 2 else 6)
    volume = np.zeros(len(coord)); np.add.at(volume, conn, size[:, None]/conn.shape[1])
    area = np.ptp(coord[:, 0])*(np.ptp(coord[:, 1]) if coord.shape[1] == 3 else 1.0)
    z = np.linspace(0, cm.L, 20001)
    out = []
    for row in [r for r in cm.scan_frames(str(folder/'mctigue')) if r['time'] > 0]:
        s = 2*np.sqrt(cm.KAPPA*row['time'])
        T = sum((-1)**k*(erfc((2*k*cm.L + z)/s) + erfc((2*(k + 1)*cm.L - z)/s)) for k in range(6))
        exact = capacity*cm.DT*np.trapezoid(T, z)
        stored = capacity*float(volume @ (model.read_field(row['frame'], 'temperature') - cm.T0))/area
        out.append(dict(time=row['time'], stored=stored, exact=exact, error=stored/exact - 1))
    return out


def build(study, output, record_source, times, space):
    summary, steps, other = (_check(s, record_source) for s in (study, times, space))
    mesh = sorted(summary['cases'], key=lambda c: -c['spacing_m'])
    finest = steps['cases'][0]; h = finest['spacing_m']
    assert abs(h - mesh[-1]['spacing_m']) < 1e-12, 'the time-step check must use the finest mesh'
    folder = times/f'h{h:g}'
    for path in folder.iterdir():
        if path.suffix in ('.cfg', '.info') or '.save.' in path.name: record_source(path)
    frames, _ = cm.compare(folder, h)
    for a, b in zip(frames, finest['frames']):
        np.testing.assert_allclose(a['pressure_error'], b['pressure_error'], rtol=0, atol=1e-12)

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    shutil.copy2(folder/'mctigue.cfg', data/'input.cfg'); shutil.copy2(folder/'start.cfg', data/'start.cfg')
    shutil.copy2(times/'summary.json', data/'summary.json'); shutil.copy2(study/'summary.json', data/'summary-mesh.json')
    shutil.copy2(study/'runner.py', data/'check_mctigue.py')
    shutil.copy2(Path(__file__).with_name('mctigue_reference.py'), data/'mctigue_reference.py')

    base = cm.Dynearthsol(str(folder/'start')); model = cm.Dynearthsol(str(folder/'mctigue'))
    rows = [r for r in cm.scan_frames(str(folder/'mctigue')) if r['time'] > 0]
    coord = base.read_field(0, 'coordinate'); depth = -coord[:, -1]
    line = (np.abs(coord[:, 0]) < 1e-9) & (depth <= .6 + 1e-9)    # both profiles vanish below 0.6 m
    order = np.argsort(depth[line])
    zf = np.linspace(0, .6, 400)
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        picks = [rows[i] for i in (0, len(rows)//4, len(rows)//2, len(rows)-1)]
        ax = axes[0]
        for color, row in zip(ds.SERIES, picks):
            t = row['time']
            ds.analytical(ax, zf, cm.REF.pressure(zf, t)/1e3, color=color, label=f'$t={t/1e3:g}$ ks')
            p = model.read_field(row['frame'], 'pore pressure')[line][order]
            keep = slice(None, None, max(1, len(order)//25))
            ax.plot(depth[line][order][keep], p[keep]/1e3, 'o', ms=3.5, mfc='white', mec=color, mew=.9)
        ds.panel(ax, 'a', 'Thermal pressurization', xlabel='depth $z$ (m)', ylabel='$p$ (kPa)')
        ds.key(ax); ax.set_xlim(0, .6)
        ax = axes[1]
        for color, row in zip(ds.SERIES, picks):
            t = row['time']
            ds.analytical(ax, zf, cm.REF.temperature(zf, t), color=color, label=f'$t={t/1e3:g}$ ks')
            T = model.read_field(row['frame'], 'temperature')[line][order] - cm.T0
            keep = slice(None, None, max(1, len(order)//25))
            ax.plot(depth[line][order][keep], T[keep], 'o', ms=3.5, mfc='white', mec=color, mew=.9)
        ds.panel(ax, 'b', 'Temperature', xlabel='depth $z$ (m)', ylabel=r'$T-T_0$ (K)')
        ds.key(ax); ax.set_xlim(0, .6)
        ax = axes[2]
        hs = np.array([c['spacing_m'] for c in mesh])
        for name, marker, color, label in METRICS:
            ax.loglog(hs, [c[name] for c in mesh], marker+'-', color=color, mfc='white', ms=4, label=label)
        # The finest mesh with twice the steps: solid markers at the same spacing.
        for name, marker, color, label in METRICS:
            ax.loglog([h], [finest[name]], marker, color=color, ms=5, mfc=color)
        ax.loglog(hs, mesh[0]['temperature_error']*(hs/hs[0])**2, '--', color=ds.MUTED, lw=.8, label='slope 2')
        ax.plot([], [], 'o', color=ds.INK, ms=5, mfc=ds.INK, label=f"$h={h:g}$ m, {steps['steps']} steps (solid)")
        ax.set_xticks(hs, [f'{x:g}' for x in hs]); ax.xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ds.panel(ax, 'c', f"Mesh convergence, {summary['steps']} steps", xlabel='mesh spacing $h$ (m)', ylabel='largest relative error')
        ax.set_ylim(1e-5, None)    # a clear strip for the legend below every series
        ax.legend(fontsize=7.5, loc='lower center', ncol=2)
        pf.save_figure(fig, output/'figures', KEY)

    orders = {name: float(np.polyfit(np.log(hs), np.log([c[name] for c in mesh]), 1)[0]) for name, *_ in METRICS}
    ledgers = {c['spacing_m']: heat_ledger(study/f"h{c['spacing_m']:g}") for c in mesh}
    first = lambda h: ledgers[h][0]
    worst3 = min(other['cases'], key=lambda c: c['spacing_m'])
    pct = lambda c, k: f'{100*c[k]:.2f}%'
    return dict(key=KEY, theme='heat', title='Thermal pressurization: McTigue surface heating',
        problem=('The drained surface of a fluid-saturated half-space is suddenly heated. Pore water expands more than the '
                 'grains, so pressure rises below the surface while it drains upward, and the solid heaves.'),
        result=(f"Largest errors over all saved times on the {h:g} m mesh with {steps['steps']} steps: pressure {pct(finest, 'pressure_error')} of its peak, "
                f"temperature {pct(finest, 'temperature_error')} of $\\Delta T$, vertical strain {pct(finest, 'strain_error')} "
                f"and surface heave {pct(finest, 'heave_error')}."),
        scope=(f"2D column {cm.L:g} m deep; mesh study at {summary['steps']} PT steps (the result above uses {steps['steps']}) to {cm.DURATION/1e3:g} ks; meshes "
               + ', '.join(f'{x:g}' for x in hs) + ' m. Mesh orders: '
               + ', '.join(f'{label} {orders[name]:.2f}' for name, _, _, label in METRICS)
               + f". With {summary['steps']} steps the finer meshes reach the time-step floor of the sequential coupling "
               f"(pressure {pct(mesh[-1], 'pressure_error')}); doubling the steps halves the pressure, strain and heave errors. "
               f"3D at {worst3['spacing_m']:g} m: pressure {pct(worst3, 'pressure_error')}, heave {pct(worst3, 'heave_error')}. "
               f"Heat ledger: stored heat against the exact input through the surface, at {first(hs[0])['time']/1e3:g} ks on meshes "
               + ', '.join(f'{x:g}' for x in hs) + ' m: ' + ', '.join(f"{100*first(x)['error']:+.3f}%" for x in hs)
               + ' (second order in the mesh; the relative error is smaller at later times because the surface-node surplus is a fixed amount). Pressure work and the '
               'mechanical energy are not in the ledger. (a) pressure; (b) temperature; (c) mesh convergence.'),
        meaning=('This checks the thermal pressurization source, the thermal stress and their two-way coupling with drainage '
                 'and deformation. Heat advection is off.'),
        metric=finest['pressure_error'], metric_label='pressure error / peak', dimension=2, spacing_m=h,
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 orders=orders, heat_ledger={f'{k:g}': v for k, v in ledgers.items()}, mesh=[{k: v for k, v in c.items() if k != 'frames'} for c in mesh],
                                 selected=dict(steps=steps['steps'], **{k: v for k, v in finest.items() if k != 'frames'}),
                                 three_d=[{k: v for k, v in c.items() if k != 'frames'} for c in other['cases']]), indent=1))
