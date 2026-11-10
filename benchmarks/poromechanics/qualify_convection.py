"""Validate the recorded convection study and plot its selected mesh.

The shared builder owns pages, setup drawings, equations and provenance.
This module recomputes growth rates from raw temperature fields.
"""
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np

import prepare_figures as pf
import diagram_style as ds
from convection_reference import growth_rate, critical_rayleigh, storage_growth_rate
sys.path.insert(0, str(pf.ROOT/'tests/functional/hydromechanics'))
import check_convection as cc

KEY = 'convection'
TAU = cc.CAPACITY/cc.LAMBDA


def sci(x):
    """LaTeX scientific notation, e.g. 5.7\\times10^{-15}."""
    m, e = f'{x:.1e}'.split('e')
    return f'{m}\\times10^{{{int(e)}}}'


def rates(summary):
    """Fitted onset and largest rate error (units of lambda/(C H^2)) of one study."""
    seeded = [c for c in summary['cases'] if c['seed_K']]
    ratio = np.array([c['ra_over_critical'] for c in seeded])*critical_rayleigh(1, 1)
    fitted = np.array([c['fitted_rate'] for c in seeded])*TAU
    slope, intercept = np.polyfit(ratio, fitted, 1)
    return -intercept/slope, float(max(abs(c['fitted_rate']-c['exact_rate'])*TAU for c in seeded))


def build(study, output, record_source, coarse=None):
    summary = json.loads((study/'summary.json').read_text())
    assert hashlib.sha256((study/'runner.py').read_bytes()).hexdigest() == summary['runner_sha256']
    assert hashlib.sha256(Path(cc.__file__).read_bytes()).hexdigest() == summary['runner_sha256']
    record_source(study/'summary.json'); record_source(study/'runner.py')
    for case in summary['cases']:
        folder = study/f"ra{case['ra_over_critical']:g}-seed{case['seed_K']:g}"
        actual = cc.amplitudes(folder)
        np.testing.assert_allclose(actual, case['series'], rtol=0, atol=1e-12)
        model = cc.Dynearthsol(str(folder/'box'))
        for row in cc.scan_frames(str(folder/'box')):
            for field in ('pore pressure', 'stress', 'velocity'):
                assert np.isfinite(model.read_field(row['frame'], field)).all(), (folder, field)
        for path in folder.glob('box.*'):
            if path.suffix == '.cfg' or path.suffix == '.info' or '.save.' in path.name:
                record_source(path)
        if case['seed_K']:
            late = slice(len(actual)//2, None)
            times = np.array(actual)[late, 0]
            amplitude = np.array(actual)[late, 1]
            fitted = np.polyfit(times, np.log(np.abs(amplitude)), 1)[0]
            expected = growth_rate(case['ra_over_critical']*critical_rayleigh(1, 1), cc.LAMBDA, cc.CAPACITY, 1, 1)
            np.testing.assert_allclose(fitted, case['fitted_rate'], rtol=0, atol=1e-15)
            np.testing.assert_allclose(expected, case['exact_rate'], rtol=0, atol=1e-15)
        else:
            assert max(s[2] for s in actual)/cc.DT < 1e-10
    seeded = [c for c in summary['cases'] if c['seed_K']]
    rest = [c for c in summary['cases'] if not c['seed_K']]
    assert seeded and rest and all(c['sign_ok'] for c in seeded)
    ra_c = critical_rayleigh(1, 1)
    ratio = np.array([c['ra_over_critical'] for c in seeded])
    fitted = np.array([c['fitted_rate'] for c in seeded])*TAU
    exact = np.array([c['exact_rate'] for c in seeded])*TAU
    slope, intercept = np.polyfit(ratio*ra_c, fitted, 1)
    onset = -intercept/slope
    worst = float(np.max(np.abs(fitted-exact)))
    biggest = max(max(s[2] for s in c['series']) for c in seeded)/cc.DT

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    top = max(seeded, key=lambda c: c['ra_over_critical'])
    folder = study/f"ra{top['ra_over_critical']:g}-seed{top['seed_K']:g}"
    shutil.copy2(folder/'box.cfg', data/'input.cfg')
    shutil.copy2(study/'summary.json', data/'summary.json')

    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        ax = axes[0]
        grid = np.linspace(0, 2.2, 200)
        ax.plot(grid, [growth_rate(r*ra_c, 1, 1, 1, 1) for r in grid], '-', color=ds.INK, lw=1.2, label='analytical')
        ax.plot(ratio, fitted, 'o', mfc='white', mec=ds.HOT, mew=1.2, ms=6, label='DES')
        ax.axhline(0, color='#9aa5ab', lw=.8); ax.axvline(1, color='#9aa5ab', lw=.8, ls=':')
        ax.set(xlabel=r'$Ra/Ra_c$', ylabel=r'$\sigma\,CH^2/\lambda$', title='(a) Growth rate of the first mode')
        ax.legend()
        ax = axes[1]
        for i, c in enumerate(seeded):
            t = np.array([s[0] for s in c['series']]); a = np.abs([s[1] for s in c['series']])
            below = c['ra_over_critical'] < 1  # blue decays below onset, red grows above it
            rank = sum(1 for d in seeded if (d['ra_over_critical'] < 1) == below and
                       abs(d['ra_over_critical']-1) > abs(c['ra_over_critical']-1))
            color = (('#1f4e79', '#5b8cbf') if below else ('#7a1f12', '#b8432a', '#d9825f'))[rank]
            # DES as markers, as the key says; a thin line keeps each series easy to follow.
            ax.semilogy(t/TAU, a/1e-3, '-', color=color, lw=.6, alpha=.6)
            ax.semilogy(t[::2]/TAU, a[::2]/1e-3, 'o', ms=3.5, mfc='white', mec=color, mew=.9,
                        label=f"$Ra/Ra_c={c['ra_over_critical']:g}$")
            late = len(t)//2
            ax.semilogy(t[late:]/TAU, a[late]/1e-3*np.exp(c['exact_rate']*(t[late:]-t[late])), '--', color=ds.INK, lw=.8)
        ax.set(xlabel=r'$t\lambda/(CH^2)$', ylabel='mode amplitude (mK)', title='(b) Mode amplitudes')
        ds.key(ax, dashed=True)
        ax = axes[2]
        model = cc.Dynearthsol(str(folder/'box')); rows = cc.scan_frames(str(folder/'box'))
        coord = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
        theta = model.read_field(rows[-1]['frame'], 'temperature') - (cc.T_TOP - cc.DT*coord[:, 1])
        level = np.abs(theta).max()
        m = ax.tricontourf(coord[:, 0], coord[:, 1]+1, conn, theta, levels=np.linspace(-level, level, 21), cmap='RdBu_r')
        ax.tricontour(coord[:, 0], coord[:, 1]+1, conn, theta, levels=np.linspace(-level, level, 11), colors=ds.INK, linewidths=.35)
        fig.colorbar(m, ax=ax, label=r'$T-T_c$ (K)', format='%.2f')
        ax.set(aspect='equal', xlabel=r'$x$ (m)', ylabel=r'$z$ (m)',
               title=f"(c) $T-T_c$ at $Ra/Ra_c={top['ra_over_critical']:g}$, $t={rows[-1]['time']/TAU:.2g}\\,CH^2/\\lambda$")
        pf.save_figure(fig, output/'figures', KEY)

    nodes = len(coord)
    rest_max = max(c['rest_max_theta_over_dT'] for c in rest)
    rest_text = ('The unseeded box stays exactly at rest.' if rest_max == 0 else
                 f'The unseeded box stays at rest to round-off: $|T-T_c|\\le{sci(rest_max)}\\,\\Delta T$.')
    convergence = ''
    if coarse:
        record_source(coarse/'summary.json')
        coarse_summary = json.loads((coarse/'summary.json').read_text())
        c_onset, c_worst = rates(coarse_summary)
        convergence = (f" On the {coarse_summary['spacing']:g} m mesh the onset error is {100*abs(c_onset-ra_c)/ra_c:.2f}% "
                       f"and the largest rate error {c_worst:.3f} (in units of $\\lambda/(CH^2)$); halving the spacing reduces the rate error {c_worst/worst:.1f}-fold "
                       f"(order {np.log2(c_worst/worst):.1f}).")
    record = dict(key=KEY, theme='heat', title='Onset of hydrothermal convection: Horton–Rogers–Lapwood',
        problem='Water in a porous box heated from below becomes lighter at the base. Above a critical Rayleigh number buoyancy overturns it into a convection cell; below it, a disturbance decays.',
        result=(f'Fitted onset: $Ra={onset:.2f}$ against $4\\pi^2={ra_c:.2f}$, {100*abs(onset-ra_c)/ra_c:.2f}% apart; largest growth-rate error: '
                f'${worst:.3f}\\,\\lambda/(CH^2)$ over {len(seeded)} Rayleigh numbers. ' + rest_text),
        scope=f"2D, {summary['spacing']:g} m mesh, {summary['steps']} steps of {summary['dt']:g} s.{convergence} (a) growth rates; (b) mode amplitudes; (c) the grown cell.",
        meaning=f'This checks Boussinesq buoyancy in Darcy flow together with heat advection, and that the conductive state stays at rest. It does not test finite-amplitude convection (see the Nusselt-number page) or 3D cells; the fastest run reaches a temperature departure of {100*biggest:.1f}% of $\\Delta T$, so its late growth is only approximately linear.',
        metric=abs(onset-ra_c)/ra_c, metric_label='onset Rayleigh number error', dimension=2, spacing_m=summary['spacing'],
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 fitted_onset=onset, max_rate_error=worst, max_departure_over_dT=biggest,
                                 cases=[{k: v for k, v in c.items() if k != 'series'} for c in summary['cases']]), indent=1))
    rows = [(c['ra_over_critical'], *sample) for c in seeded for sample in c['series']]
    np.savetxt(data/'history.csv', rows, delimiter=',',
        header='Ra_over_Rac,time_s,mode_K,max_departure_K', comments='')
    shutil.copy2(study/'runner.py', data/'check_convection.py')
    shutil.copy2(Path(__file__).with_name('convection_reference.py'), data/'convection_reference.py')
    return record


def build_3d(study, output, record_source):
    """Use the independently replayed, same-property 3D refinement study."""
    q = json.loads((study/'qualification.json').read_text())
    assert q['passed'] and len(q['cases']) == 11
    record_source(study/'qualification.json')
    for name, key in (('runner.py', 'runner_sha256'), ('shared_runner.py', 'shared_runner_sha256'),
                      ('reference.py', 'reference_sha256')):
        assert cc.digest(study/name) == q[key]
        record_source(study/name)
    selected = next(c for c in q['cases'] if c['folder'] == q['selected'])
    folder = study/q['selected']
    for name, digest in selected['outputs_sha256'].items():
        if (folder/name).exists():
            assert cc.digest(folder/name) == digest
            record_source(folder/name)
    fine = sorted((c for c in q['cases'] if c['spacing_m'] == .03125 and c['seed_K'] and c['dt_s'] == 500),
                  key=lambda c: c['rayleigh'])
    assert len(fine) == 3 and all(c['correct_sign'] for c in fine)
    rest = next(c for c in q['cases'] if not c['seed_K'])
    assert rest['rest_max_theta_over_dT'] == 0
    data = output/'data/convection_3d'; data.mkdir(parents=True, exist_ok=True)
    for c in q['cases']:
        record_source(study/c['folder']/'summary.json')
        record_source(study/c['folder']/'box.cfg')
    shutil.copy2(folder/'box.cfg', data/'input.cfg')
    shutil.copy2(study/'qualification.json', data/'summary.json')
    shutil.copy2(study/'runner.py', data/'check_convection_3d.py')
    shutil.copy2(study/'reference.py', data/'convection_reference.py')
    model = cc.Dynearthsol(str(folder/'box'))
    x = model.read_field(0, 'coordinate')
    theta = model.read_field(40, 'temperature') - (cc.T_TOP-cc.DT*x[:, 2])
    cut = np.isclose(x[:, 1], .25)
    assert cut.sum() == 33**2
    halves = next(c for c in q['cases'] if c['dt_s'] == 250)
    temporal = abs(halves['fitted_rate_s_inv']/selected['fitted_rate_s_inv']-1)
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        ax = axes[0]
        ratio = np.linspace(.6, 1.5, 101)
        rac = selected['critical_rayleigh']
        reference = [storage_growth_rate(r*rac, selected['hydraulic_to_thermal_diffusivity']*r/1.4,
                     cc.LAMBDA, cc.CAPACITY, 1, 1, transverse=1, width_y=1)*TAU for r in ratio]
        ds.analytical(ax, ratio, reference)
        ds.des(ax, [c['rayleigh']/rac for c in fine], [c['fitted_rate_s_inv']*TAU for c in fine])
        ax.axhline(0, color=ds.MUTED, lw=.6)
        ax.axvline(1, color=ds.MUTED, lw=.6, ls=':')
        ds.panel(ax, 'a', 'Growth of the (1, 1, 1) mode', xlabel=r'$Ra/Ra_c$', ylabel=r'$\sigma H^2/\kappa_T$')
        ax.legend()
        ax = axes[1]
        for i, c in enumerate(fine):
            values = np.array(c['series']); late = len(values)//2
            # Only the late amplitude is anchored; the analytical slope is not fitted.
            t = values[late:, 0]; amplitude = values[late:, 1]
            color = ds.SERIES[i]
            ax.semilogy(t/TAU, amplitude/amplitude[0], 'o', ms=3, mfc='white', color=color,
                        label=f"$Ra/Ra_c={c['rayleigh']/rac:g}$")
            ax.semilogy(t/TAU, np.exp(c['reference_rate_s_inv']*(t-t[0])), '-', color=color, lw=1)
        ds.panel(ax, 'b', 'Late exponential growth or decay', xlabel=r'$t\kappa_T/H^2$', ylabel=r'$A(t)/A(0.1)$')
        ax.legend(fontsize=10)
        ax = axes[2]
        level = float(abs(theta[cut]).max()*1000)
        m = ax.tricontourf(x[cut, 0], x[cut, 2]+1, theta[cut]*1000,
                          levels=np.linspace(-level, level, 21), cmap='RdBu_r')
        ax.set_aspect('equal'); fig.colorbar(m, ax=ax, label=r'$T-T_c$ (mK)', format='%.2f')
        ds.panel(ax, 'c', r'Cell section at $y/H=0.25$', xlabel=r'$x/H$', ylabel=r'$(z+H)/H$')
        pf.save_figure(fig, output/'figures', 'convection_3d')
    rows = [[c['spacing_m'], c['dt_s'], c['rayleigh']/c['critical_rayleigh'],
             c['fitted_rate_s_inv'], c['reference_rate_s_inv'], c['relative_growth_error']]
            for c in q['cases'] if c['seed_K']]
    np.savetxt(data/'rates.csv', rows, delimiter=',',
               header='spacing_m,dt_s,Ra_over_Rac,DES_rate_s_inv,exact_rate_s_inv,relative_error', comments='')
    worst = max(c['relative_growth_error'] for c in fine)
    return dict(key='convection_3d', theme='heat', title='Three-dimensional porous convection cells',
        problem='A unit porous cube heated from below supports a temperature mode that varies in all three directions. Its growth changes sign at a known Rayleigh number.',
        result=(f"At 1/32 m spacing, growth-rate errors are "
                + ', '.join(f"{100*c['relative_growth_error']:.2f}%" for c in fine)
                + f" below and above onset. Halving the physical timestep changes the selected growth rate by {100*temporal:.4f}%. The unseeded cube stays exactly at conduction."),
        scope='3D rigid saturated matrix, explicit Darcy and heat transport, 35,937 nodes and 163,840 tetrahedra. All three spatial meshes use the same 10 Pa fluid bulk modulus; the analytical reference retains its finite pressure storage. No material or amplitude fitting.',
        meaning='This checks three-dimensional buoyancy and heat advection, growth and decay, a rest state, and separate spatial/time controls. The amplitude in panel (b) is normalized at dimensionless time 0.1; only the analytical growth slope is compared. It does not qualify nonlinear steady planforms or realistic water compressibility.',
        metric=worst, metric_label='largest 3D modal growth-rate error', dimension=3, spacing_m=.03125,
        tables=[('Spatial and temporal rates', 'rates.csv'), ('Raw-field replay and source hashes', 'summary.json')],
        evidence=json.dumps(dict(source_commit=selected['source_commit'], executable_sha256=selected['executable_sha256'],
                                 reference_sha256=q['reference_sha256'], temporal_rate_change=temporal,
                                 selected=q['selected'], finite_storage=True), indent=2))


def build_elder(study, output, record_source):
    """Replay the selected Elder fields against the final independent grid."""
    import check_elder as ce

    q = json.loads((study/'qualification.json').read_text())
    assert q['raw_replay_passed']
    record_source(study/'qualification.json')
    for name, digest in q['files_sha256'].items():
        path = study/name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, path
        record_source(path)
    folder = study/q['selected']
    required = {study/'ref-full.npz', study/'ref-full.json', folder.parent/'summary.json',
                folder/'elder.cfg', folder/'elder.poly', folder/'elder.info',
                *folder.glob('elder.save.*')}
    assert {str(p.relative_to(study)) for p in required} <= set(q['files_sha256']), 'unsealed Elder input'
    archived = json.loads((folder.parent/'summary.json').read_text())
    case = min(archived['cases'], key=lambda c: c['spacing_m'])
    assert folder.name == f"h{case['spacing_m']:g}"
    reference_file = study/'ref-full.npz'
    ref = np.load(reference_file)
    n = max(int(k.split('_')[1]) for k in ref.files if k.startswith('theta_'))
    assert n == 256 and 'theta_128' in ref.files
    series, nodes, elements = ce.compare(folder, ref, 1, archived['dt_s'])
    # The archived study used an earlier reference grid. Its DES heat is
    # invariant; field errors below are recomputed against the final grid.
    np.testing.assert_allclose([s['heat_des'] for s in series],
                               [s['heat_des'] for s in case['series']], rtol=1e-12, atol=1e-14)
    for item in series:
        assert np.isfinite(list(item.values())).all()
    assert series == q['series'] and nodes == case['nodes'] and elements == case['elements']
    key = 'elder'; data = output/'data'/key; data.mkdir(parents=True, exist_ok=True)
    shutil.copy2(folder/'elder.cfg', data/'input.cfg')
    shutil.copy2(study/'ref-full.json', data/'reference.json')
    (data/'summary.json').write_text(json.dumps(dict(qualification=q, archived_case=case,
        reference_grid=n, series=series), indent=2)+'\n')
    columns = ('tau', 'max_abs_difference', 'rms_difference', 'heat_des', 'heat_reference')
    np.savetxt(data/'history.csv', [[s[k] for k in columns] for s in series], delimiter=',',
               header=','.join(columns), comments='')
    from plot_elder import fields
    fields(folder.parent, reference_file, output/'figures'/key, [.1], include_heat=True)
    late = series[-1]
    heat_error = abs(late['heat_des']/late['heat_reference']-1)
    early = max(abs(s['heat_des']/s['heat_reference']-1) for s in series)
    return dict(key=key, theme='heat', title="Elder's transient heater problem",
        problem='Heating half the base of a sealed porous layer creates plumes that rise and merge. There is no closed-form solution; an independent stream-function solver provides the reference.',
        result=f"At dimensionless time 0.1, the selected DES mesh differs by {100*heat_error:.2f}% in heat content and {late['max_abs_difference']:.2f} in temperature scaled by the 10 K heater contrast. The early heat discrepancy reaches {100*early:.1f}%.",
        scope=f"2D symmetric half box, Ra = 400, {nodes:,} nodes and {elements:,} triangles, 50 s physical steps. Panels (a,b) show the final time; panel (c) retains all six comparison times. The 128 × 256 and 256 × 512 references differ by up to 0.12 in scaled temperature and less than 0.5% in heat content.",
        meaning='This is a measured nonlinear application limit: early DES plumes differ, and coarse reference grids follow another branch. The incompressible reference and DES with deliberately scaled finite storage are distinct approximations. No material or timestep fitting is used.',
        metric=heat_error, metric_label='final heat-content difference', dimension=2,
        spacing_m=case['spacing_m'], status='preliminary',
        tables=[('Six comparison times', 'history.csv'), ('Field replay and source hashes', 'summary.json'),
                ('Reference grid convergence', 'reference.json')],
        evidence=json.dumps(dict(executable_sha256=archived['executable_sha256'],
            source_commit=q['source_commit'], reference_sha256=hashlib.sha256(reference_file.read_bytes()).hexdigest(),
            archived_reference_sha256=archived['reference_sha256'], limits='Early plume and finite-storage discrepancies retained.'), indent=2))
