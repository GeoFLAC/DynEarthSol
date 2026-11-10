"""Application page: rainfall builds a perched water table until a slope fails.

Interface as qualify_theis: build(study, output, record_source, pressure_driven)
re-reads the finest rainfall run from raw saves, writes figures/rainfall.* and
data/rainfall/, and returns the review record. pressure_driven is the slope
study with the same soil (phi = 30 degrees, c = 1 kPa) loaded by bedrock
pressure, for comparison.
"""
import json
import math
from pathlib import Path
import shutil

import numpy as np

import prepare_figures as pf
import diagram_style as ds
import slope_application as sa
import rainfall_slope as rs

KEY = 'rainfall'


SLIP = 1e-3                               # mean plastic strain of the basal 0.6 m that marks failure (see qualify_slope)


def _onset(case):
    return next((f for f in case['frames'] if f['basal_plastic_fraction'] > 0), None)


def _at_slip(case, key):
    """Value of key (time or base pressure) when the mean basal plastic strain first reaches SLIP."""
    e = np.array([f['basal_plastic_strain'] for f in case['frames']])
    return float(np.interp(SLIP, e, [f[key] for f in case['frames']])) if e.max() >= SLIP else float('nan')


def build(study, output, record_source, pressure_driven):
    summary, driven = pf.load_study(study, record_source), pf.load_study(pressure_driven, record_source)
    cases = sorted(summary['cases'], key=lambda c: -c['spacing_m'])
    finest = cases[-1]
    folder = study/f"h{finest['spacing_m']:g}"
    for path in folder.iterdir():
        if path.suffix in ('.cfg', '.info', '.poly') or '.save.' in path.name or '.chkpt.' in path.name: record_source(path)
    again, _, _ = sa.analyse(folder)
    np.testing.assert_allclose([f['base_pressure'] for f in again], [f['base_pressure'] for f in finest['frames']], rtol=1e-12)
    assert all(f['surface_pressure_max'] <= 0 for c in cases for f in c['frames']), 'rain was forced into a saturated surface'

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    shutil.copy2(folder/'slope.cfg', data/'input.cfg')
    for name in ('seed.cfg', 'slope.poly'): shutil.copy2(folder/name, data/name)
    shutil.copy2(study/'summary.json', data/'summary.json'); shutil.copy2(pressure_driven/'summary.json', data/'summary-pressure-driven.json')
    shutil.copy2(study/'runner.py', data/'rainfall_slope.py')

    driven_cases = {c['spacing_m']: c for c in driven['cases'] if c['frames']}
    hour = 3600.0
    model = sa.Dynearthsol(str(folder/'slope')); rows = sa.scan_frames(str(folder/'slope'))
    first = rows[0]['frame']
    x = model.read_field(first, 'coordinate'); conn = model.read_field(first, 'connectivity')
    t = math.tan(sa.BETA)
    base = np.abs(x[:, 1] + x[:, 0]*t + sa.D) < 1e-9
    order = np.argsort(x[base, 0])
    u_point, u_davis, u_star = sa.estimates(summary['cohesion'], math.radians(summary['friction_deg']))
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3, width=13.0)
        ax = axes[0]
        p = model.read_field(rows[-1]['frame'], 'pore pressure')
        tri = pf.matplotlib.tri.Triangulation(x[:, 0], x[:, 1], conn)
        art = ax.tripcolor(tri, p/1e3, cmap='RdBu_r', vmin=-12, vmax=12, shading='gouraud')
        ax.tricontour(tri, p, levels=[0.0], colors=ds.INK, linewidths=1.2)
        ax.set_aspect('equal'); ax.set_xlim(14, 26); ax.set_ylim(-26*t - sa.D - .5, -14*t + .5)
        fig.colorbar(art, ax=ax, location='bottom', shrink=.7, pad=.18,
                     label='pore pressure (kPa)')
        ds.panel(ax, 'a', f"Perched water table just before sliding ({rows[-1]['time']/hour:.1f} h)",
                 xlabel='$x$ (m)', ylabel='$z$ (m)')
        ax.text(.03, .05, 'black line: $p=0$', transform=ax.transAxes, fontsize=10, color=ds.INK)

        ax = axes[1]
        for color, case in zip((ds.SERIES[1], ds.HOT), cases):
            h = case['spacing_m']
            for c, ls, what in ((case, '-', 'rain'), (driven_cases.get(h), '--', 'bedrock pressure')):
                if c is None: continue
                pb = np.array([f['base_pressure'] for f in c['frames']])/1e3
                e = np.array([f['basal_plastic_strain'] for f in c['frames']])
                keep = e > 0
                ax.semilogy(pb[keep], e[keep], ls, color=color, lw=1.2, label=f'{what}, {h:g} m')
        for u, label, color, ls in ((u_point, 'pointwise yield', ds.MUTED, ':'), (u_davis, 'Davis, non-associated', ds.INK, '-.')):
            ax.axvline(u/1e3, color=color, ls=ls, lw=1.0, label=label)
        ax.axhline(SLIP, color=ds.MUTED, lw=.6)
        ax.set_xlim(6, 16); ax.set_ylim(1e-6, 3e-2)
        ds.panel(ax, 'b', 'Basal slip: rain against bedrock pressure', xlabel='base pore pressure (kPa)',
                 ylabel='mean basal plastic strain')
        ax.legend(fontsize=11, loc='upper left', bbox_to_anchor=(0, -.25)).set_in_layout(False)

        ax = axes[2]
        picks = np.linspace(1, len(rows) - 1, 5).astype(int)
        for color, i in zip(ds.SERIES + [ds.INK], picks):
            pp = model.read_field(rows[i]['frame'], 'pore pressure')[base][order]/1e3
            ax.plot(x[base, 0][order], pp, '-', color=color, lw=1.1, label=f"{rows[i]['time']/hour:.1f} h")
        ax.axhline(0, color=ds.MUTED, lw=.6); ax.axvspan(0, sa.END, color='#eef1f3', zorder=0)
        ax.axvspan(sa.L - sa.END, sa.L, color='#eef1f3', zorder=0)
        ax.text(sa.END/2, .97, 'strong', transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=9.5, color=ds.MUTED)
        ds.panel(ax, 'c', 'Base pressure along the slope', xlabel='$x$ (m)', ylabel='base pore pressure (kPa)')
        ax.legend(fontsize=11, loc='upper left', bbox_to_anchor=(0, -.25), ncol=2).set_in_layout(False)
        pf.save_figure(fig, output/'figures', KEY)

    on = _onset(finest)
    fail_t = [_at_slip(c, 'time') for c in cases]; fail_u = [_at_slip(c, 'base_pressure') for c in cases]
    driven_u = {h: _at_slip(c, 'base_pressure') for h, c in driven_cases.items()}
    du = [u - driven_u.get(c['spacing_m'], float('nan')) for u, c in zip(fail_u, cases)]
    return dict(key=KEY, theme='mechanics', status='preliminary', title='Application: rainfall triggers a shallow slide',
        problem=('Rain falls on a partly saturated soil layer over impermeable bedrock. The water spreads down to the bedrock, '
                 'a perched water table rises, and the slope fails when the pore pressure on the bedrock is high enough.'),
        result=(f"Basal yield starts after {on['time']/hour:.1f} h at a base pressure of {on['base_pressure']/1e3:.1f} kPa on the "
                f"{finest['spacing_m']:g} m mesh, and the basal slip reaches a plastic strain of {SLIP:g} after "
                + ' and '.join(f'{v/hour:.2f}' for v in fail_t) + ' h on the ' + ' and '.join(f"{c['spacing_m']:g}" for c in cases)
                + f" m meshes, at {fail_u[-1]/1e3:.1f} kPa: {du[-1]/1e3:.1f} kPa more than the same soil needs when the "
                'pressure comes from the bedrock. The surface never saturates.'),
        scope=(f"2D, the slope of the pressure-driven page with $c={summary['cohesion']/1e3:g}$ kPa and "
               f"$\\phi={summary['friction_deg']:g}^\\circ$, "
               f"Gardner retention ($\\alpha={summary['alpha']*1e4:g}\\times10^{{-4}}$ Pa$^{{-1}}$, $S_r={summary['sr']:g}$), Bishop $\\chi=S_e$, initial "
               f"suction {summary['suction']/1e3:g} kPa (Bishop stress {summary['bishop0']/1e3:.1f} kPa), rain "
               f"{summary['rain']*3.6e6:.0f} mm/h; sealed bedrock, draining toe. Meshes "
               + '; '.join(f"{c['spacing_m']:g} m: onset {_onset(c)['base_pressure']/1e3:.2f} kPa, plastic strain {SLIP:g} at "
                           f"{_at_slip(c, 'base_pressure')/1e3:.2f} kPa after {_at_slip(c, 'time')/hour:.2f} h, last PT equilibrium "
                           f"{c['frames'][-1]['base_pressure']/1e3:.2f} kPa after {c['frames'][-1]['time']/hour:.2f} h" for c in cases)
               + '. Loaded from the bedrock (gravity applied elastically first), the same soil reaches that slip at '
               + ', '.join(f'{v/1e3:.2f}' for _, v in sorted(driven_u.items(), reverse=True)) + ' kPa on the '
               + ', '.join(f'{h:g}' for h in sorted(driven_u, reverse=True)) + f" m meshes; estimates: pointwise onset "
               f"{u_point/1e3:.1f} kPa, Davis {u_davis/1e3:.1f} kPa, classical factor of safety {u_star/1e3:.1f} kPa. As on the "
               'pressure-driven page, the slip curve converges and the last PT equilibrium does not, so failure is read at a fixed '
               'slip. DES has no runoff or ponding, so rain that meets a saturated surface would be forced in; the rate is low '
               'enough that the surface stays unsaturated throughout. (a) pressure field; (b) slip against base pressure; '
               '(c) base pressure along the slope.'),
        meaning=('Rain fails the slope at a higher base pressure than pressure from the bedrock does: the perched table '
                 'leaves suction (Bishop strength) above it and seepage parallel to the slope, instead of an upward seepage '
                 'force through the whole layer. The time to failure is the same on both meshes. No exact solution exists; '
                 'this is an application of the verified unsaturated flow, Bishop stress and effective-stress plasticity together.'),
        metric=fail_t[-1]/hour, metric_label=f'time to basal plastic strain {SLIP:g} (h)', dimension=2, spacing_m=finest['spacing_m'],
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 onset=[_onset(c)['base_pressure'] for c in cases], slip=SLIP, failure_time=fail_t,
                                 failure_pressure=fail_u, failure_pressure_driven=driven_u,
                                 last_equilibrium=[c['frames'][-1]['base_pressure'] for c in cases]), indent=1))
