"""Selected physical-inertia fault-heating patch and its independent reference."""
import ast
import hashlib
import json
import shutil

import numpy as np

import diagram_style as ds
import fault_heat as fault
import prepare_figures as pf

KEY = 'fault_heat'


def build(study, output, record):
    qualification = json.loads((study/'qualification.json').read_text())
    assert qualification['passed']
    for name, expected in qualification['files_sha256'].items():
        path = study/name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, name
        record(path)
    for name in ('qualification.json', 'onset-checks.json', 'time-convergence.json'):
        record(study/name)
    folder = study/'physical-large-slip'
    summary = json.loads((folder/'summary.json').read_text())
    assert summary['passed'] and summary['mode'] == 'physical'
    assert summary['shear_error_over_initial'] < 2e-5
    assert summary['normal_stress_error_over_load'] < 3e-5
    series = np.loadtxt(folder/'history.csv', delimiter=',', skiprows=1)
    reference = fault.reference(series[:, 0]/fault.HEIGHT)
    np.testing.assert_allclose(series[:, 11:], reference.T, rtol=1e-12, atol=1e-7)
    # Keep the actual executed reference bound to the published source.
    def function(path):
        tree = ast.parse(path.read_text())
        return ast.dump(next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'reference'), include_attributes=False)
    assert function(folder/'runner.py') == function(pf.ROOT/'benchmarks/poromechanics/fault_heat.py')
    model = fault.heat.Dynearthsol(str(folder/'cell'))
    frames = fault.heat.scan_frames(model.modelname)
    np.testing.assert_array_equal(series[:, 0], [r['time'] for r in frames])
    chosen = np.unique(np.r_[np.arange(0, len(series), 12), len(series)-1])
    for i in chosen:
        frame = frames[i]['frame']
        np.testing.assert_allclose(series[i, 1:4], model.read_field(frame, 'stress').mean(axis=0), rtol=1e-13)
        np.testing.assert_allclose(series[i, 4], model.read_field(frame, 'temperature').mean()-300, rtol=1e-13, atol=1e-13)
        np.testing.assert_allclose(series[i, 5], model.read_field(frame, 'pore pressure').mean(), rtol=1e-13, atol=1e-10)
    for path in folder.iterdir():
        if path.suffix in ('.cfg', '.info', '.csv', '.json', '.py') or '.save.' in path.name or '.chkpt.' in path.name:
            record(path)
    target = output/'data'/KEY
    target.mkdir(parents=True, exist_ok=True)
    for source, name in ((folder/'cell.cfg', 'input.cfg'), (folder/'history.csv', 'history.csv'),
                         (folder/'summary.json', 'summary.json'), (study/'qualification.json', 'qualification.json'),
                         (study/'onset-checks.json', 'onset-checks.json'),
                         (study/'time-convergence.json', 'time-convergence.json')):
        shutil.copy2(source, target/name)
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        slip = series[:, 0]*1000
        for ax, actual, exact, scale, title, ylabel, letter in (
            (axes[0], 3, 12, 1e-6, 'Thermal weakening', r'$\tau$ (MPa)', 'a'),
            (axes[1], 4, 14, 1., 'Plastic work becomes heat', r'$T-T_0$ (K)', 'b'),
            (axes[2], 5, 13, 1e-6, 'Sealed fluid pressurization', '$p$ (MPa)', 'c')):
            ds.analytical(ax, slip, series[:, exact]*scale, label='MC / Jaumann reference')
            ds.des(ax, slip[chosen], series[chosen, actual]*scale, size=3, label='DES, physical mass')
            ds.panel(ax, letter, title, xlabel='slip $\delta$ (mm)', ylabel=ylabel)
            ax.legend(fontsize=9.5)
        control = fault.reference(series[:, 0]/fault.HEIGHT, expansion=False)
        axes[0].plot(slip, control[1]/1e6, '--', color=ds.MUTED, lw=1., label='zero-expansion reference')
        axes[0].legend(fontsize=9.5, loc='lower left')
        pf.save_figure(fig, output/'figures', KEY)
    reduction = 1-summary['final_shear_Pa']/fault.TAU0
    return dict(key=KEY, theme='heat', title='Fault-zone heating at fixed normal traction',
        problem='A prestressed, sealed fault patch slips at 1 m/s under 50 MPa normal compression. Plastic dissipation heats the mixture; expanding water raises pore pressure and reduces shear resistance.',
        result=(f"At 100 mm slip, shear resistance falls {100*reduction:.1f}%, temperature rises "
                f"{summary['final_temperature_rise_K']:.2f} K and pore pressure reaches {summary['final_pressure_Pa']/1e6:.2f} MPa. "
                f"Saved-frame shear errors are {100*summary['shear_error_over_initial']:.4f}% of the initial resistance against the differential reference."),
        meaning='Ordinary explicit mechanics with physical mixture inertia and zero damping follows the thermal weakening path. Companion GVS runs and slow-to-fast physical restarts check the mass/timestep transition and retain pending plastic heat separately from kinetic energy.',
        scope='2D homogeneous two-triangle patch with one free normal platen mode; constant-friction Mohr–Coulomb, zero dilation, fixed reference heat capacity, no Darcy flux or conduction. The reference imposes normal equilibrium; physical inertia gives a small measured departure. There is no spatial rupture front, flash heating or complete earthquake cycle.',
        metric=summary['shear_error_over_initial'], metric_label='shear error', dimension=2, spacing_m=fault.HEIGHT,
        tables=[('Every saved state and reference', 'history.csv'), ('Time convergence', 'time-convergence.json'),
                ('First physical intervals after GVS', 'onset-checks.json')],
        evidence=json.dumps(dict(selected_case='physical-large-slip', numerical_source=summary,
                                  qualification='qualification.json'), indent=1))
