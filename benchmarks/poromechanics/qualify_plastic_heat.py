"""Export the finest prescribed EP-RSF heating cell and its continuum reference."""
import json
import shutil

import numpy as np

import diagram_style as ds
import plastic_heat as heat
import prepare_figures as pf

KEY = 'plastic_heat'


def build(study, output, record):
    qualification = json.loads((study/'qualification.json').read_text())
    assert qualification['passed']
    for name in ('qualification.json', 'build.json', 'continuum.json'):
        record(study/name)
    convergence = json.loads((study/'continuum.json').read_text())
    finest = convergence[-1]
    assert finest['steps'] == 1000
    assert finest['temperature_max_relative'] < .002
    assert finest['pressure_max_relative'] < .002
    for field in ('temperature_max_relative', 'pressure_max_relative'):
        errors = np.array([c[field] for c in convergence])
        assert np.all(np.log2(errors[:-1]/errors[1:]) > .98)
    folder = study/'coupled-2d-n1000'
    summary = json.loads((folder/'summary.json').read_text())
    series = np.loadtxt(folder/'continuum.csv', delimiter=',', skiprows=1)
    model = heat.Dynearthsol(str(folder/'cell'))
    frames = heat.scan_frames(model.modelname)
    np.testing.assert_array_equal(series[:, 0], [r['time'] for r in frames])
    theta0 = model.read_field(0, 'friction state variable')[0]
    tau, energy = heat.rsf_continuum(series[:, 0], theta0, summary['element_widths_m'][0])
    np.testing.assert_allclose(series[:, 2], tau, rtol=1e-13, atol=1e-8)
    np.testing.assert_allclose(series[:, 4], energy/summary['heat_capacity_J_m3_K'], rtol=1e-13, atol=1e-14)
    storage = .2/2.2e9+.6/1e11
    np.testing.assert_allclose(series[:, 6], 4e-5/storage*series[:, 4], rtol=1e-13)
    selected = np.unique(np.r_[np.arange(0, len(series), 25), len(series)-1])
    for i in selected:
        stress = model.read_field(frames[i]['frame'], 'stress')
        np.testing.assert_allclose(series[i, 1], ((stress[:, 0]-stress[:, 1])/2).mean(), rtol=1e-14)
        np.testing.assert_allclose(series[i, 3], model.read_field(frames[i]['frame'], 'temperature').mean()-heat.T0, rtol=1e-13)
        np.testing.assert_allclose(series[i, 5], model.read_field(frames[i]['frame'], 'pore pressure').mean(), rtol=1e-13)
    for path in folder.iterdir():
        if path.suffix in ('.cfg', '.info', '.json', '.csv') or '.save.' in path.name or '.chkpt.' in path.name:
            record(path)

    target = output/'data'/KEY
    target.mkdir(parents=True, exist_ok=True)
    shutil.copy2(folder/'cell.cfg', target/'input.cfg')
    for source, name in ((folder/'continuum.csv', 'history.csv'),
                         (folder/'summary.json', 'summary.json'),
                         (study/'continuum.json', 'convergence.json'),
                         (study/'qualification.json', 'qualification.json'),
                         (study/'build.json', 'build.json')):
        shutil.copy2(source, target/name)
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        for ax, actual, reference, factor, title, label, letter in (
            (axes[0], 1, 2, 1e-6, 'Yielding and state evolution', r'$\tau$ (MPa)', 'a'),
            (axes[1], 3, 4, 1., 'Plastic work becomes heat', r'$T-T_0$ (K)', 'b'),
            (axes[2], 5, 6, 1e-6, 'Sealed thermal pressurization', '$p$ (MPa)', 'c')):
            ds.analytical(ax, series[:, 0], series[:, reference]*factor, label='continuum reference')
            ds.des(ax, series[selected, 0], series[selected, actual]*factor, label='DES', size=3)
            ds.panel(ax, letter, title, xlabel='$t$ (s)', ylabel=label)
            ax.legend(loc='best', fontsize=10)
        pf.save_figure(fig, output/'figures', KEY)
    return dict(key=KEY, theme='heat', title='Plastic shear heating and sealed pressurization',
        problem='A porous cell deforms at a prescribed pure-shear rate. Plastic work in aging rate-and-state plasticity heats the cell and expands the sealed fluid.',
        result=(f"At Δt = 0.0002 s, temperature and pressure differ from the continuum reference by at most "
                f"{100*finest['temperature_max_relative']:.3f}% and {100*finest['pressure_max_relative']:.3f}% "
                "of their final reference rises. Both errors converge at first order. Heat plus pending energy matches accepted plastic work to round-off."),
        meaning='This checks plastic earthquake heating, including evolving EP-RSF strength, without introducing inertia or a rupture-front error. The plotted temperature and pressure are the physical fields; the pending numerical reservoirs are accounted for separately.',
        scope='Prescribed spin-free deformation, fixed mesh and heat capacity, unit heat fraction, no heat loss or fluid flux. Constant effective confinement: total normal stress reacts to pore pressure. This does not test thermal fault weakening at fixed total normal traction, a propagating rupture, or a complete earthquake cycle.',
        metric=finest['pressure_max_relative'], metric_label='pressure error', dimension=2,
        spacing_m=1., tables=[('Every-step history and reference', 'history.csv'), ('Time convergence', 'convergence.json')],
        evidence=json.dumps(dict(implementation_commit=qualification['implementation_commit'],
            build='build.json', qualification='qualification.json', convergence=convergence), indent=1))
