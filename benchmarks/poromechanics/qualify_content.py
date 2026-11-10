"""Selected water, thermal and yielding-cell pages from frozen qualifications.

The original raw studies remain immutable. Figures reuse the study renderers;
portable tables retain the fixed-stress pending increment separately from
physical fluid content. Every consumed artifact is checked against its seal.
"""
import json
import shutil

import numpy as np

import diagram_style as ds
import fluid_content as water
import prepare_figures as pf
import thermal_cell as heat
from Dynearthsol import scan_frames


def study(folder, record_source):
    path = folder/'qualification.json'
    q = json.loads(path.read_text())
    record_source(path)

    def consume(relative):
        path = folder/relative
        assert water.digest(path) == q['artifacts_sha256'][relative], path
        record_source(path)
        return path

    return q, consume


def portable(output, key, consume, cfg, history, summary):
    folder = output/'data'/key
    folder.mkdir(parents=True, exist_ok=True)
    for original, name in ((cfg, 'input.cfg'), (history, 'history.csv'), (summary, 'summary.json')):
        shutil.copy2(consume(original), folder/name)
    return folder


def render(output, key, plot):
    folder = output/'figures'/key
    folder.mkdir(parents=True, exist_ok=True)
    plot(folder)
    for ext in ('svg', 'png', 'pdf'):
        (folder/f'result.{ext}').rename(folder.parent/f'{key}.{ext}')
    folder.rmdir()


def build_water(folder, output, record_source):
    q, consume = study(folder, record_source)
    selected = q['selected']
    case = folder/selected
    summary = json.loads(consume('spatial-pt/summary.json').read_text())
    fine = summary['cases'][-1]
    assert fine['elements'] == 80 and fine['late_pressure_error_over_steady'] < 3e-4
    assert q['max_content_budget_error'] < 3e-15
    assert water.reference_checks()['quadrature_balance']
    key = 'content'
    data = portable(output, key, consume, selected+'/input.cfg', selected+'/history.csv',
                    'spatial-pt/summary.json')
    a = np.loadtxt(data/'history.csv', delimiter=',', skiprows=1)
    np.testing.assert_allclose(a[:, 3], a[:, 1]+a[:, 2], rtol=1e-10, atol=1e-20)
    assert np.max(abs(a[:, 3]-a[:, 4])) < 1e-19
    # The profile renderer reads geometry at frame zero and four pressure saves.
    rows = scan_frames(str(case/'column'))
    consume(selected+'/column.info')
    frames = {0, *(min(rows, key=lambda r: abs(r['time']-t))['frame'] for t in (.1, .5, 1., 4.))}
    for frame in frames:
        consume(selected+f'/column.save.{frame:06d}')
    render(output, key, lambda target: water.plot_result(case, target))
    return dict(key=key, tables=[('History CSV', 'history.csv'), ('Control tables', 'summary.json')], theme='consolidation', title='Water balance in a deforming column',
        problem='Water enters a laterally confined elastic column through its base and leaves through a pressure-controlled top. Pressure and deformation independently measure its stored fluid.',
        result=f"Late pressure error is {100*fine['late_pressure_error_over_steady']:.4f}% of the steady scale. The coupled split budget closes within {q['max_content_budget_error']:.2g} of injected volume across 19 cases.",
        scope='Selected 2D PT column: 0.05 m mesh, 0.005 s step, 0–4 s. Lines are the mixed-boundary diffusion series; symbols are DES. (a) pressure profiles; (b) physical stored content and drain outflow. The small pending PT increment is tabulated separately.',
        meaning='Independent mesh-integrated pressure and strain, boundary reaction and Darcy outflow distinguish a physical budget from the pressure-equation residual. Temporal controls measure first-order PT lag; ordinary explicit, adaptive GVS, tiny 3D and restart controls are supporting evidence.',
        metric=fine['late_pressure_error_over_steady'], metric_label='pressure error / steady pressure',
        dimension=2, spacing_m=.05, evidence=json.dumps(dict(qualification='runs/fluid-content-20261108/qualification.json',
            solver_commit=q['solver_commit'], selected=selected, cases=q['cases'], saved_frames=q['saved_frames'],
            restart_suites=q['restart_suites'], exact_restart_payloads=q['exact_restart_payloads'], selected_metrics=fine), indent=2))


def build_thermal(folder, output, record_source):
    q, consume = study(folder, record_source)
    assert q['fixed_capacity_checks_passed'] and not q['full_open_energy_qualified']
    key = 'thermal_cell'
    data = portable(output, key, consume, 'native2d/explicit-fixed/heated.cfg',
                    'native2d/explicit-fixed/history.csv', 'native2d/summary.json')
    summary = json.loads((data/'summary.json').read_text())
    fine = next(c for c in summary['cases'] if c['method'] == 'explicit' and c['boundary'] == 'fixed')
    closed = np.loadtxt(data/'history.csv', delimiter=',', skiprows=1)
    np.testing.assert_allclose(closed[:, 2], heat.BETA*(closed[:, 1]-heat.T0)/heat.STORAGE, atol=2e-8)
    np.testing.assert_allclose(closed[:, 7], closed[:, 8], rtol=1e-12, atol=1e-7)
    for name in ('pt-fixed', 'explicit-drained'):
        shutil.copy2(consume(f'native2d/{name}/history.csv'), data/(name+'.csv'))
    render(output, key, lambda target: heat.plot_result(folder/'native2d', target, controls=False))
    return dict(key=key, tables=[('History CSV', 'history.csv'), ('Control tables', 'summary.json')], theme='heat', title='Thermal pressurization in a sealed cell',
        problem='Uniform supplied heat warms a rigid, sealed porous material. The temperature rise and fluid expansion give independent heat and pressure references.',
        result=f"The 10 K pulse gives {fine['final_pressure_Pa']/1000:.3f} kPa. Sealed pressure and fixed-capacity heat match their references to round-off; general open-system energy remains unqualified.",
        scope='Selected ordinary-explicit 2D unit cell; 40 K/s for 0.25 s, then a hold to 1 s. Solid thermal expansion is zero. (a) sealed pressure against temperature; (b) stored and supplied heat. PT lag and the drained energy control are retained in the tables.',
        meaning='Thirty-one cases cover 2D/3D, native/HDF5 and timestep controls; eight thermal restarts are exact. PT has a first-order pending-pressure lag during heating. The drained retained-fluid ledger has a persistent 10,450 J deficit (0.03446% of supplied heat), so the sealed comparison does not qualify thermal-plume energy.',
        metric=max(fine['maximum_errors'].values()), metric_label='sealed reference error / reference scale',
        dimension=2, spacing_m=.5, evidence=json.dumps(dict(qualification='runs/thermal-cell-20261108/qualification.json',
            solver_commit=q['solver_commit'], selected='native2d/explicit-fixed', selected_metrics=fine,
            cases=q['replayed_cases'], saved_frames=q['saved_frames'], full_open_energy_qualified=False,
            exact_restart_payloads=q['exact_payloads']), indent=2))


def plastic_result(csv, target):
    a = np.loadtxt(csv, delimiter=',', skiprows=1)
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(2)
        ds.analytical(axes[0], a[:, 0], 4*a[:, 0], label='continuous undrained reference')
        ds.des(axes[0], a[::5, 0], a[::5, 4]/1000, label='PT')
        ds.panel(axes[0], 'a', 'Pressure through plastic yield', xlabel='time (s)', ylabel='mean pressure (kPa)')
        ds.analytical(axes[1], a[:, 0], -a[:, 6]*1e6, label='imposed pore-volume loss')
        ds.des(axes[1], a[::5, 0], (a[::5, 1]-a[::5, 6]+a[::5, 2])*1e6,
               label='pressure storage + pending increment')
        ds.panel(axes[1], 'b', 'Independent volume balance', xlabel='time (s)',
                 ylabel='volume (mL per metre thickness)')
        pf.save_figure(fig, target, 'result')


def build_plastic(folder, output, record_source):
    q, consume = study(folder, record_source)
    assert q['primary_cases'] == 48 and q['exact_restart_payloads'] == 5760
    key = 'plastic_content'
    data = portable(output, key, consume, 'qualified-2d/pt-psi10/column.cfg',
                    'qualified-2d/pt-psi10/history.csv', 'qualified-2d/summary.json')
    summary = json.loads((data/'summary.json').read_text())
    fine = next(c for c in summary['cases'] if c['case'] == 'pt-psi10')
    a = np.loadtxt(data/'history.csv', delimiter=',', skiprows=1)
    p, storage, lag = 0., .25/1e6+.75/1e20, 1/(2e7+1e7/3)
    for row in a:
        if row[0]: p = (.001*row[0]+lag*p)/(storage+lag)
        np.testing.assert_allclose(row[4], p, rtol=1e-12, atol=1e-9)
    assert np.max(abs(a[:, 3])) < 1e-17 and a[:, 7].max() > 1e-5
    for dim in ('3d', '2d-hdf', '3d-hdf'):
        shutil.copy2(consume(f'qualified-{dim}/summary.json'), data/(dim+'.json'))
    render(output, key, lambda target: plastic_result(data/'history.csv', target))
    return dict(key=key, tables=[('History CSV', 'history.csv'), ('Control tables', 'summary.json')], theme='mechanics', title='Water balance through plastic yielding',
        problem='A sealed porous column is compressed through Mohr–Coulomb yield. Total volume change fixes its mean pressure independently of the plastic stress return.',
        result='The corrected source conserves fluid content plus the pending PT increment to round-off. The parent drift was 1.24% in the selected 2D case and 4.49% in the 3D dilation control.',
        scope='Selected 2D plane-strain unit cell: 100 Pa cohesion, 20° friction, 10° dilation, 0.001 m/s top compression, 0.01 s step. (a) pressure; (b) pressure storage plus pending increment against imposed volume loss. Porosity and mesh remain fixed.',
        meaning='Forty-eight primary cases cover both dimensions and containers, negative/zero/positive dilation, hydraulic substeps, open drainage, Gardner suction, evolving-pressure EP-RSF and thermal expansion. Old checkpoints preserve earlier offsets. Conservative remapping, finite-grain plastic porosity and general energy remain open.',
        metric=float(np.max(abs(a[:, 3]))/.001), metric_label='split content error / imposed volume',
        dimension=2, spacing_m=.5, evidence=json.dumps(dict(qualification='runs/plastic-source-20261108/qualification.json',
            solver_commit=q['implementation_commit'], selected='qualified-2d/pt-psi10', selected_metrics=fine,
            primary_cases=q['primary_cases'], saved_frames=q['primary_saved_frames'],
            exact_restart_payloads=q['exact_restart_payloads'], core_comparisons=q['core_comparisons']), indent=2))


def build_mass(folder, output, record_source):
    """Sealed finite-density compression, with the physical split defect visible."""
    q, consume = study(folder, record_source)
    assert q['passed'] and q['native_hdf5_dimensions'] == [2, 3]
    key = 'mass'
    selected = '2d-native/mass/lag/4000'
    data = portable(output, key, consume, selected+'/column.cfg', selected+'/history.csv',
                    '2d-native/mass/lag/summary.json')
    summary = json.loads((data/'summary.json').read_text())
    a = np.loadtxt(data/'history.csv', delimiter=',', skiprows=1)
    assert summary['passed'] and summary['cases'][-1]['pressure_error'] < .005
    np.testing.assert_allclose(a[:, 3]+a[:, 4], a[:, 5], rtol=0, atol=2e-9)
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(2)
        ds.analytical(axes[0], a[:, 0], a[:, 2]/1000)
        ds.des(axes[0], a[:, 0], a[:, 1]/1000)
        ds.panel(axes[0], 'a', 'Sealed compression', xlabel='time (s)', ylabel='pressure (kPa)')
        ds.analytical(axes[1], a[:, 0], -a[:, 4]*1000, label='negative numerical buffer')
        ds.des(axes[1], a[:, 0], (a[:, 3]-a[:, 5])*1000, label='physical water-mass defect')
        ds.panel(axes[1], 'b', 'Physical mass and coupling lag', xlabel='time (s)', ylabel='mass change (g)')
        pf.save_figure(fig, output/'figures', key)
    fine = summary['cases'][-1]
    return dict(key=key, theme='consolidation', title='Finite-density water in a sealed cell',
        problem='A sealed saturated cell is compressed at a prescribed rate. The fluid EOS and evolving Biot porosity give an independent pressure and physical water inventory.',
        result=f"The finest step gives {100*fine['pressure_error']:.3f}% pressure error. Halving the step halves both the pressure error and physical mass defect; the separate mass-plus-buffer budget closes within 2 micrograms.",
        scope='2D unit cell, 0.01/s logarithmic compression to 0.1 s, 4,000 physical steps. (a) DES and the implicit sealed-cell solution; (b) physical mass defect and negative numerical buffer. The temporal table retains 1,000 and 2,000 steps.',
        meaning='The stabilised split conserves M+R. R is a numerical lag, contributes no density or heat capacity, and must not be described as physical water. Native/HDF5, 3D, restart, remap, source and transport controls are supporting evidence.',
        metric=fine['pressure_error'], metric_label='sealed pressure error', dimension=2, spacing_m=1.,
        tables=[('History CSV','history.csv'),('Temporal convergence','summary.json')],
        evidence=json.dumps(dict(qualification=f'runs/{folder.name}/qualification.json',
            solver_commit=q['solver_commit'], selected=selected, selected_metrics=fine), indent=2))


def build_drained_heat(folder, output, record_source):
    """Heat supply, retained EOS mass and exported enthalpy in a rigid drain."""
    from mass_reference import Cell
    q, consume = study(folder, record_source)
    assert q['passed']
    key = 'drained_heat';selected = '2d-native/mass/heat/drained-200'
    data = portable(output, key, consume, selected+'/column.cfg', selected+'/history.csv',
                    '2d-native/mass/heat/summary.json')
    a = np.loadtxt(data/'history.csv', delimiter=',', skiprows=1)
    ref = Cell();theta = np.linspace(0, a[-1, 2]-273., 301);temperature = theta+273.
    exact = np.array([ref.drained_heat(t) for t in temperature])
    outflow = np.array([ref.drained_outflow_heat(t) for t in temperature])
    retained = a[:, 9]-a[0, 9]
    np.testing.assert_allclose(retained-a[:, 5]-a[:, 7], 0, rtol=0, atol=2e-5)
    np.testing.assert_allclose(a[:, 3], 250*np.exp(-.0002*(a[:, 2]-273.)), rtol=2e-12)
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(2)
        ds.analytical(axes[0], theta, exact/1e6)
        ds.des(axes[0], a[::10, 2]-273., a[::10, 5]/1e6)
        ds.panel(axes[0], 'a', 'Heat needed to warm a drained cell', xlabel='temperature rise (K)', ylabel='supplied heat (MJ)')
        ds.analytical(axes[1], theta, outflow/1e6, label='analytical exported heat')
        ds.des(axes[1], a[::10, 2]-273., -a[::10, 7]/1e6, label='DES boundary heat')
        ds.panel(axes[1], 'b', 'Water leaving carries heat', xlabel='temperature rise (K)', ylabel='exported sensible heat (MJ)')
        pf.save_figure(fig, output/'figures', key)
    cases = json.loads((data/'summary.json').read_text())['cases']
    errors = [c['reference_error'] for c in cases if c['kind']=='drained']
    return dict(key=key, theme='heat', title='Heating a drained cell with changing water mass',
        problem='A rigid saturated cell receives uniform heat while every face holds zero pore pressure. Thermal expansion expels water, reducing retained heat capacity and exporting sensible heat.',
        result=f"Temperature errors at 100 and 200 steps are {errors[0]:.6f} and {errors[1]:.6f} K. Retained energy plus measured boundary export closes against supplied heat within 20 microjoules.",
        scope='2D unit volume per metre thickness, 30 MW supplied for 5 s, 200 steps, no gravity, solid expansion or conduction. (a) total supplied heat; (b) exported sensible heat. The reference uses the actual exponential water density.',
        meaning='At zero pressure and zero bulk strain this isolates open-system sensible heat without mechanical boundary work. Owned solid heat capacity stays fixed while fluid mass changes. It does not qualify full thermodynamic pressure work or a thermal plume.',
        metric=errors[-1], metric_label='maximum temperature error (K)', dimension=2, spacing_m=1.,
        tables=[('Water and energy history','history.csv'),('Cell controls','summary.json')],
        evidence=json.dumps(dict(qualification=f'runs/{folder.name}/qualification.json',
            solver_commit=q['solver_commit'], selected=selected, temperature_errors_K=errors), indent=2))


def build_mixing(folder, output, record_source):
    """Well-mixed injection with a different inlet density and changing storage."""
    from mass_reference import Cell
    q, consume = study(folder, record_source)
    assert q['passed']
    key = 'mixing';selected = '2d-native/mass/heat/expanding-mixing-200'
    data = portable(output, key, consume, selected+'/column.cfg', selected+'/history.csv',
                    '2d-native/mass/heat/summary.json')
    a = np.loadtxt(data/'history.csv', delimiter=',', skiprows=1)
    ref=Cell();temperature=np.linspace(273.,a[-1,2],301)
    times=np.array([ref.injection_time(t,323.,.01) for t in temperature])
    np.testing.assert_allclose(a[:,9]-a[0,9]-a[:,6]-a[:,7],0,rtol=0,atol=2e-5)
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(2)
        ds.analytical(axes[0],times,temperature-273.)
        ds.des(axes[0],a[::10,0],a[::10,2]-273.)
        ds.panel(axes[0],'a','Warm-water injection',xlabel='time (s)',ylabel='temperature rise (K)')
        axes[1].plot(a[:,0],a[:,6]/1e6,color=ds.SERIES[0],label='injected heat')
        axes[1].plot(a[:,0],-a[:,7]/1e6,color=ds.SERIES[1],label='boundary export')
        ds.des(axes[1],a[::10,0],(a[::10,9]-a[0,9])/1e6,label='retained energy gain')
        ds.panel(axes[1],'b','Open energy budget',xlabel='time (s)',ylabel='sensible heat (MJ)')
        pf.save_figure(fig,output/'figures',key)
    cases=json.loads((data/'summary.json').read_text())['cases']
    errors=[c['reference_error'] for c in cases if c['kind']=='expanding-mixing']
    return dict(key=key,theme='heat',title='Warm-water mixing in a drained cell',
        problem='Warm water enters a pressure-controlled, well-mixed saturated cell. Natural outflow removes water at the cell temperature; inlet and retained water have different EOS densities.',
        result=f"Maximum temperature error is {errors[-1]:.6f} K at 200 steps and falls with timestep refinement. Retained energy equals injected sensible heat minus boundary export within 20 microjoules.",
        scope='2D unit cell, 0.01 m³/s inlet at 323 K, initial 273 K, 0–5 s, 200 steps. (a) temperature and the independent integral solution; (b) injected, exported and retained sensible heat.',
        meaning='This checks inlet-temperature ownership, injection density, outflow enthalpy and changing retained fluid storage together. Sources are distributed uniformly by nodal volume, so it is a mixing control rather than a spatial plume benchmark.',
        metric=errors[-1],metric_label='maximum temperature error (K)',dimension=2,spacing_m=1.,
        tables=[('Water and energy history','history.csv'),('Cell controls','summary.json')],
        evidence=json.dumps(dict(qualification=f'runs/{folder.name}/qualification.json',
            solver_commit=q['solver_commit'],selected=selected,temperature_errors_K=errors),indent=2))
