#!/usr/bin/env python3
"""Build a selective, portable implementation review from verified snapshots.

One predefined finest-grid result per benchmark; no parameter fitting.
Raw experiment paths remain intact behind the review's evidence appendix.
"""

import argparse
import configparser
import html
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

import consolidation_refinement as cr
import prepare_figures as pf
import diagram_style as ds
from unsaturated_reference import Retention, srivastava_yeh, steady_pressure
from qualify_wave import load_case as load_wave_case
from qualify_heat import load_case as load_heat_case
from qualify_convection import build as build_convection, build_3d as build_convection_3d, build_elder, sci
import qualify_theis, qualify_mandel, qualify_cryer, qualify_aquifer, qualify_rudnicki
import qualify_salencon
import qualify_hantush, qualify_bishop, qualify_mctigue, qualify_bishop_cell, qualify_nusselt, qualify_regional, qualify_pressure_yield, qualify_drainage, qualify_injection_temperature, qualify_borehole, qualify_rsf, qualify_slope, qualify_rainfall, qualify_embankment, qualify_drawdown
import qualify_transport, qualify_kirsch, qualify_strength, qualify_content
import qualify_plastic_heat
import qualify_fault_heat
from natural_math import typeset
from review import ORDER, REFERENCE_ONLY, SLUGS, load_catalogues, write_review



def archived_path(name, data, group):
    """Locate a sealed artifact after its results tree has been transferred."""
    path = Path(name)
    if not path.is_absolute():
        return group/path
    if path.is_relative_to(data):
        return path
    parts = path.parts
    return data.joinpath(*parts[parts.index('results')+1:])


def terzaghi_3d(data, record):
    """One sentence on the 3D sanity run of the full-history study, if recorded."""
    path = data/'runs/coupled-3d-20261109/model-terzaghi-3d/results/summary.json'
    if not path.exists(): path = data/'runs/terzaghi-full-20261108/sanity-3d/summary.json'
    if not path.exists(): return ''
    record(path)
    runs = json.loads(path.read_text())['runs']
    finest = min((r for r in runs if r['dt_s'] == min(x['dt_s'] for x in runs)), key=lambda r: r['spacing_m'])
    return (f"3D check: {100*finest['max_pressure_error_over_load']:.3f}% of the load on the {finest['spacing_m']:g} m mesh. ")

def portable_inputs(data, output, consumed, record):
    """Keep the selected input and its mesh/restart dependencies together."""
    originals = [(data/name, digest) for name, digest in consumed.items() if name.endswith('.cfg')]
    changes = {}
    for path in sorted((output/'data').glob('*/input.cfg')):
        digest = cr.digest(path)
        matches = [source for source, expected in originals if expected == digest]
        assert matches, f'Unrecorded input: {path}'
        source = matches[0]
        cfg = configparser.ConfigParser(interpolation=None); cfg.optionxform = str
        cfg.read(path)
        dependencies, rewritten = [], False

        def copy_dependency(original):
            assert original.is_file(), original
            target = path.parent/original.name
            if target.exists():
                assert cr.digest(target) == cr.digest(original), target
            else:
                shutil.copy2(original, target)
            record(original)
            dependencies.append(original.name)

        for section in cfg.sections():
            for name, value in list(cfg[section].items()):
                if not name.lower().endswith('filename') or not value.strip(): continue
                original = Path(value)
                if original.is_absolute():
                    original = archived_path(original, data, source.parent)
                else:
                    original = source.parent/original
                copy_dependency(original)
                if value != original.name:
                    cfg[section][name] = original.name; rewritten = True
        if cfg.getboolean('sim', 'is_restarting', fallback=False):
            prefix = Path(cfg['sim']['restarting_from_modelname'])
            prefix = archived_path(prefix, data, source.parent)
            frame = cfg.getint('sim', 'restarting_from_frame')
            for kind in ('save', 'chkpt'):
                original = prefix.with_name(f'{prefix.name}.{kind}.{frame:06d}')
                if not original.is_file(): original = original.with_suffix(original.suffix+'.vtkhdf')
                copy_dependency(original)
            info = prefix.with_name(prefix.name+'.info')
            if info.is_file(): copy_dependency(info)
            if cfg['sim']['restarting_from_modelname'] != prefix.name:
                cfg['sim']['restarting_from_modelname'] = prefix.name; rewritten = True
        if rewritten:                     # paths made local; provenance.json keeps the original's path and hash
            with path.open('w') as stream: cfg.write(stream)
        changes[path.parent.name] = dict(source=str(source.relative_to(data)), source_sha256=digest,
                                        paths_rewritten=rewritten, dependencies=dependencies)
    (output/'REPRODUCE.md').write_text('''# Running the examples

## Build DES

From the DES source folder:

```bash
make ndims=2      # dynearthsol2d
make ndims=3      # dynearthsol3d
```

Put the executables on your `PATH`, or call them with their full path.
For the Python reader and ParaView converter, set the source-folder path:

```bash
export DES_SOURCE=/path/to/DynEarthSol
export PYTHONPATH="$DES_SOURCE${PYTHONPATH:+:$PYTHONPATH}"
```

## Run an example

Each example's files are in `data/<example>/`: `input.cfg`, and any mesh
(`.poly`) or restart frames (`.save`, `.chkpt`, `.info`) it starts from.
Copy that folder somewhere writable and run the commands given in the
example's "Input and how to run" section, for instance

```bash
mkdir -p ~/runs
cp -r data/theis ~/runs/theis
cd ~/runs/theis
export OMP_NUM_THREADS=1
dynearthsol3d input.cfg
```

`dynearthsol2d -h` lists every input option with its default and meaning.
Each page lists the settings that define its example.

## Look at the output

DES names its output after `sim.modelname`: `<model>.info` lists the saved
frames (frame, step, time in seconds), `<model>.save.NNNNNN` holds the
fields of each frame and `<model>.chkpt.NNNNNN` the restart checkpoints.

```python
import sys; sys.path.append('/path/to/DynEarthSol')
from Dynearthsol import Dynearthsol
model = Dynearthsol('theis')
p = model.read_field(model.frames[-1], 'pore pressure')
```

`python3 "$DES_SOURCE/2vtk.py" <model>` converts the frames for ParaView; add `-g` for a
run that restarted from another, to include the frames before the restart.

## Reproduce a figure

Each page names the script that ran the whole study behind its figure (all
meshes and time steps). The scripts are in `benchmarks/poromechanics/` and
`tests/functional/hydromechanics/`. `provenance.json` records the exact
source and data versions the published figures were made with.
For an exact historical reproduction, use the solver/source version recorded
there; a newer solver can contain corrections absent from the archived run.
''')
    return changes


def build(data, output, refresh=False):
    if output.exists() and refresh:
        old = json.loads((output / 'provenance.json').read_text())
        keys = {r['key'] for r in old['cases']}
        assert {'terzaghi', 'gardner', 'van_genuchten', 'brooks_corey'} <= keys <= set(ORDER)
    output.mkdir(parents=True, exist_ok=refresh)
    for name in ('figures', 'data', 'appendix'):
        (output / name).mkdir(exist_ok=refresh)
    dynamic_doc = pf.ROOT/'doc/hydro/dynamics.md'
    (output/'appendix/physical-inertia.md').write_text(dynamic_doc.read_text().replace(
        '(index.md)', '(../README.md)'))
    mechanics_doc = pf.ROOT/'doc/hydro/mechanics.md'
    (output/'appendix/mechanics.md').write_text(mechanics_doc.read_text().replace(
        '(dynamics.md)', '(physical-inertia.md)'))
    aquifer_doc = pf.ROOT/'doc/hydro/aquifer.md'
    aquifer_text = aquifer_doc.read_text().split('## DES pressure storage:')[0]
    aquifer_text = aquifer_text.replace('../../benchmarks/poromechanics/results/review/appendix/analytical-solutions.html', 'analytical-solutions.html')
    (output/'appendix/averaged-aquifer.md').write_text(aquifer_text)
    supporting_notes = {
        'terzaghi-pulse.md': 'terzaghi-pulse.md',
        'slope-audit.md': 'slope.md',
        'rainfall-audit.md': 'rainfall.md',
        'plastic-heat-implementation.md': 'plastic-heat.md',
        'fault-heat-reference.md': 'fault-heat.md',
        'mass-model.md': 'mass.md',
        'performance.md': 'performance.md',
    }
    for target, source in supporting_notes.items():
        shutil.copy2(pf.ROOT/'doc/hydro'/source, output/'appendix'/target)
    records = []
    consumed = {}

    def record(path):
        consumed[str(path.relative_to(data))] = cr.digest(path)

    performance = data/'runs/performance-20261110'
    if (performance/'qualification.json').exists():
        shutil.copytree(performance, output/'appendix/performance-data', dirs_exist_ok=refresh)
        for path in performance.rglob('*'):
            if path.is_file(): record(path)

    def value(cfg, key):
        return json.loads(cfg['mat'][key])[0]

    def save_case(key, folder):
        target = output / 'data' / key
        target.mkdir(exist_ok=refresh)
        shutil.copy2(folder / ('column.cfg' if (folder / 'column.cfg').exists() else 'input.cfg'), target / 'input.cfg')
        record(folder / ('column.cfg' if (folder / 'column.cfg').exists() else 'input.cfg'))
        for path in folder.glob('*.save.*'):
            record(path)
        for path in folder.glob('*.info'):
            record(path)
        return target

    # Selection is prescribed by the refinement study, not by a fitted curve.
    # The full drainage history, to T_v = 2; the earlier 400 s study covered only T_v < 0.1.
    group = data / 'runs/terzaghi-full-20261108/history-2d'
    source = json.loads((group / 'summary.json').read_text())
    assert source['passed']
    run = next(r for r in source['runs'] if r['key'] == '2d-h0.25-dt0.25')
    folder = group / run['key'] / 'terzaghi-2d'
    assert cr.digest(folder / 'input.cfg') == run['input_sha256']
    for name, digest in run['saves_sha256'].items():
        assert cr.digest(folder / name) == digest
    case = pf.load_case(group / run['key'], 2)
    np.testing.assert_allclose(case['metrics']['max_pressure_error_over_load'],
                               run['max_pressure_error_over_load'], rtol=0, atol=1e-12)
    key = 'terzaghi'
    target = save_case(key, folder)
    record(group / 'summary.json')
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(2)
        last = case['times'][-1]
        depth = np.linspace(0, case['reference'].height, 401)
        cv = case['reference'].height**2/case['reference'].diffusivity
        shown = [i for i, t in enumerate(case['times']) if t > 0][::2] + [len(case['times']) - 1]
        for color, i in zip(ds.SERIES, sorted(set(shown))):
            t = case['times'][i]
            ds.analytical(axes[0], case['reference'].pressure(depth, [t])[0]/1000, depth, color=color,
                          label=f'$T_v={t/cv:.2f}$')
            ds.des(axes[0], case['pressure'][i]/1000, case['depth'], color=color, label=None)
        axes[0].invert_yaxis()
        ds.panel(axes[0], 'a', 'Pressure profiles', xlabel='excess pressure (kPa)', ylabel='depth below the drained top (m)')
        ds.key(axes[0], loc='upper right')
        active = case['times'] > 0
        t = np.linspace(1e-3*last, last, 801)
        ds.analytical(axes[1], t, case['reference'].settlement(t)*1000)
        axes[1].set_xlim(0, 1.02*last); axes[1].set_ylim(0, None)
        ds.des(axes[1], case['times'][active], case['settlement'][active]*1000, size=5)
        axes[1].axhline(case['reference'].drained_settlement*1000, color=ds.MUTED, lw=.8, ls='--', label='drained')
        ds.panel(axes[1], 'b', 'Settlement over the full drainage', xlabel='$t$ (s)', ylabel='settlement (mm)')
        pf.save_figure(fig, output / 'figures', key)
    np.savetxt(target / 'pressure.csv', np.column_stack((case['depth'], case['pressure'][-1],
        case['reference'].pressure(case['depth'], [last])[0])), delimiter=',', header='depth_m,DES_Pa,reference_Pa', comments='')
    pe, se = case['metrics']['max_pressure_error_over_load'], case['metrics']['max_settlement_error_over_drained']
    records.append(dict(key=key, theme='consolidation', title='Terzaghi consolidation of a column',
        problem='A sudden load compresses a saturated column. Water drains through the top, the pressure falls and the column settles.',
        result=f'Largest pressure error: {100*pe:.3f}% of the applied load; largest settlement error: {100*se:.3f}% of the final drained settlement.',
        scope=('2D plane strain, 0.25 m mesh, 0.25 s steps, compared every 1200 s to 8400 s, '
               f"$T_v=c_vt/H^2={last/cv:.2f}$ (about 99% consolidated). Maximum pressure errors on 1, 0.5 and 0.25 m meshes: "
               + ', '.join(f"{100*r['max_pressure_error_over_load']:.3f}%" for r in source['runs'] if r['dt_s'] == .25)
               + ' of the load; halving the step from 1 s changes the result at first order. '
               + terzaghi_3d(data, record) + '(a) pressure profiles; (b) settlement history.'),
        meaning='This checks fluid storage and elastic coupling together, through both pressure and settlement. It does not test a sudden load with GVS or large strain.',
        metric=pe, metric_label='pressure error / load', dimension=2, spacing_m=.25,
        evidence='runs/terzaghi-full-20261108/history-2d/index.html'))

    integration = data / 'unsaturated-integration-20261107'
    check = json.loads((integration / 'integration-summary.json').read_text())
    assert check['checks_passed'] == 76
    record(integration / 'integration-summary.json')
    for dim in (2, 3):
        record(integration / f'unsaturated-{dim}d/summary.json')
    for name, digest in json.loads((integration / 'figures/provenance.json').read_text())['files'].items():
        assert cr.digest(integration / name) == digest

    def load(folder):
        cfg = configparser.ConfigParser(); cfg.read(folder / 'column.cfg')
        model = pf.Dynearthsol(str(folder / 'column'))
        frame = model.frames[-1]
        times = {r['frame']: r['time'] for r in pf.scan_frames(str(folder / 'column'))}
        z = model.read_field(0, 'coordinate')[:, -1] + cfg.getfloat('mesh', 'zlength')
        p = model.read_field(frame, 'pore pressure')
        return cfg, model, frame, times[frame], z, p

    key = 'gardner'
    folder = integration / 'unsaturated-3d/transient-0.025'
    cfg, model, frame, time, z, pressure = load(folder)
    target = save_case(key, folder)
    gamma = value(cfg, 'fluid_rho0')*cfg.getfloat('control', 'gravity')
    alpha = value(cfg, 'retention_alpha')
    kr = np.exp(alpha*pressure)
    opts = dict(length=cfg.getfloat('mesh', 'zlength'), alpha=alpha*gamma,
                ks=value(cfg, 'hydraulic_perm')*gamma/value(cfg, 'fluid_visc'),
                dtheta=value(cfg, 'porosity')*(1-value(cfg, 'residual_saturation')),
                qa=0, qb=-json.loads(cfg['bc']['hbc_flux'])[-1])
    exact = srivastava_yeh(z, time, **opts)
    selected = check['dimensions']['3']['transient'][-1]
    error = float(np.max(np.abs(kr-exact)))
    np.testing.assert_allclose(error, selected['max_error'][-1], rtol=0, atol=1e-9)
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, (ax,) = ds.result_figure(1)
        grid = np.linspace(0, opts['length'], 401)
        ds.analytical(ax, grid, srivastava_yeh(grid, time, **opts))
        ds.des(ax, z, kr, size=3)
        ds.panel(ax, None, f'Conductivity profile at {time:g} s', xlabel='height above the base (m)', ylabel='$K/K_s$ from the pressure')
        pf.save_figure(fig, output / 'figures', key)
    np.savetxt(target / 'profile.csv', np.column_stack((z, pressure, kr, exact)), delimiter=',',
               header='elevation_m,pressure_Pa,DES_K_over_Ks,reference_K_over_Ks', comments='')
    records.append(dict(key=key, theme='unsaturated', title='Transient infiltration: Gardner law',
        problem='Water infiltrates a partly saturated column from the top while the base stays at the water table.',
        result=f'Largest error in $K/K_s$: {100*error:.3f}% at {time:g} s and {100*max(selected["max_error"]):.2f}% over all saved times.',
        scope=f'3D, 0.025 m mesh, all saved times compared; shown at {time:g} s. Halving the mesh reduces the error about fivefold.',
        meaning='This checks conservative unsaturated flow against an independent transient solution. It does not test Bishop effective stress.',
        metric=error, metric_label='final K/Ks error', dimension=3, spacing_m=.025,
        evidence='unsaturated-integration-20261107/figures/index.html'))

    for law, model_id, title in [('van_genuchten', 1, 'van Genuchten'), ('brooks_corey', 2, 'Brooks–Corey')]:
        key = law
        folder = integration / f'unsaturated-2d/steady-{law}-0.025'
        cfg, model, frame, time, z, pressure = load(folder)
        assert value(cfg, 'retention_model') == model_id
        target = save_case(key, folder)
        retention = Retention(law, value(cfg, 'retention_alpha'), value(cfg, 'retention_n'), value(cfg, 'residual_saturation'))
        opts = dict(flux=-json.loads(cfg['bc']['hbc_flux'])[-1], perm=value(cfg, 'hydraulic_perm'),
                    mu=value(cfg, 'fluid_visc'), rho=value(cfg, 'fluid_rho0'), g=cfg.getfloat('control', 'gravity'))
        exact = steady_pressure(retention, z, **opts)
        selected = check['dimensions']['2']['steady'][law][-1]
        error = float(np.max(np.abs(pressure-exact)))
        np.testing.assert_allclose(error, selected['max_error_pa'], rtol=0, atol=1e-7)
        with pf.plt.rc_context(ds.FIGURE_RC):
            fig, (ax,) = ds.result_figure(1)
            grid = np.linspace(0, cfg.getfloat('mesh', 'zlength'), 401)
            ds.analytical(ax, grid, steady_pressure(retention, grid, **opts)/1000)
            ds.des(ax, z, pressure/1000, size=4)
            ds.panel(ax, None, f'Steady pressure profile: {title}', xlabel='height above the base (m)', ylabel='pressure (kPa)')
            pf.save_figure(fig, output / 'figures', key)
        np.savetxt(target / 'profile.csv', np.column_stack((z, pressure, exact)), delimiter=',', header='elevation_m,DES_Pa,reference_Pa', comments='')
        records.append(dict(key=key, theme='unsaturated', title=f'Steady infiltration: {title} law',
            problem=('A constant inflow from the top produces a steady suction profile above a water table at the base. '
                     + {'van Genuchten': 'The van Genuchten law has no air-entry pressure, so the soil desaturates as soon as the pressure falls below zero.',
                        'Brooks–Corey': 'The Brooks–Corey law stays saturated up to its air-entry suction, so the profile has a kink where desaturation begins.'}.get(title, '')),
            result=f'Largest pressure error: {error:.2f} Pa, {100*error/selected["p_range"]:.2f}% of the profile range.',
            scope='2D, 0.025 m mesh, steady state; storage does not enter the steady solution.',
            meaning='This checks the retention-dependent conductivity and the flux boundary. It does not test transient storage or coupled slope behaviour.',
            metric=error, metric_label='maximum pressure error (Pa)', dimension=2, spacing_m=.025,
            evidence='unsaturated-integration-20261107/figures/index.html'))

    group = data / 'runs/physical-pressure-fix-20261107'
    wave_check = json.loads((group/'wave-qualification.json').read_text())
    assert wave_check['passed']
    record(group/'wave-qualification.json')
    for name, sha in wave_check['consumed_files'].items():
        path = archived_path(name, data, group)
        assert cr.digest(path) == sha, path
        record(path)
    chosen = wave_check['selected']
    study = group/chosen['study']
    wave_consumed = {}
    wave = load_wave_case(study, chosen['steps'], wave_consumed)
    for name in wave_consumed:
        record(Path(name))
    key = 'undrained_wave'
    folder = study/f'steps-{chosen["steps"]}'
    target = output/'data'/key
    target.mkdir(exist_ok=refresh)
    shutil.copy2(folder/'wave.cfg', target/'input.cfg')
    shutil.copy2(study/'undrained_wave.py', target/'undrained_wave.py')
    record(study/'undrained_wave.py')
    x, times = wave['coordinate'], wave['times']
    centre = (x.min(axis=0)+x.max(axis=0))/2
    mid = int(np.argmin(np.linalg.norm(x-centre,axis=1)))
    base_point = centre.copy(); base_point[-1] = x[:,-1].min()
    base = int(np.argmin(np.linalg.norm(x-base_point,axis=1)))
    phase = times/wave['period_s']
    v = wave['velocity'][:,mid]/wave['velocity_amplitude']
    p = wave['pressure'][:,base]/wave['pressure_amplitude']
    zphase = np.pi*(x[:,-1]-x[:,-1].min())/np.ptp(x[:,-1])
    exact_v = np.sin(zphase[mid])*np.cos(2*np.pi*phase)
    exact_p = -np.cos(zphase[base])*np.sin(2*np.pi*phase)
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(2)
        t = np.linspace(0,1,401)
        ds.analytical(axes[0],t,np.sin(zphase[mid])*np.cos(2*np.pi*t))
        ds.des(axes[0],phase,v,size=4)
        ds.panel(axes[0],'a','Velocity at mid-height',xlabel='$t/T$',ylabel='$v_z/V$')
        ds.analytical(axes[1],t,-np.cos(zphase[base])*np.sin(2*np.pi*t))
        ds.des(axes[1],phase,p,size=4)
        ds.panel(axes[1],'b','Pressure at the sealed base',xlabel='$t/T$',ylabel='$p$ / pressure amplitude')
        pf.save_figure(fig,output/'figures',key)
    np.savetxt(target/'history.csv',np.column_stack((times,v,exact_v,p,exact_p)),delimiter=',',
               header='time_s,DES_v_over_V,reference_v_over_V,DES_p_over_P,reference_p_over_P',comments='')
    np.savetxt(target/'initial-velocity.csv',np.column_stack((x,wave['velocity'][0])),delimiter=',',
               header='x_m,y_m,z_m,initial_vz_m_s',comments='')
    ve, pe = (wave['maximum'][k] for k in ('velocity_error_over_amplitude','pressure_error_over_amplitude'))
    records.append(dict(key=key,theme='dynamics',title='Undrained wave in a sealed column',
        problem='A small initial velocity sets a sealed saturated column oscillating; compression raises the pressure and expansion lowers it.',
        result=f'Largest velocity error: {100*ve:.3f}% and largest pressure error: {100*pe:.3f}% of their amplitudes, over all nodes and 65 saved times.',
        scope=f'3D, 0.0125 m mesh, 8192 steps over one {1000*wave["period_s"]:.3f} ms period $T$; fixed mesh, no damping, physical inertia. (a) velocity at mid-height; (b) pressure at the sealed base.',
        meaning='This checks fluid storage and elastic wave coupling together; local fluid content is conserved to round-off. It does not test a GVS earthquake cycle.',
        metric=pe,metric_label='pressure error / wave amplitude',dimension=3,spacing_m=.0125,
        evidence='runs/physical-pressure-fix-20261107/wave-qualification.json'))

    heat_group = data/'runs/private-des-integration-20261107'
    for name in ('heat-source-sha256.json', 'heat-source.diff', 'heat-implementation-commit.txt'):
        record(heat_group/name)
    heat_check = json.loads((heat_group/'heat-qualification.json').read_text())
    assert heat_check['passed']
    record(heat_group/'heat-qualification.json')
    for name, sha in heat_check['consumed_files'].items():
        path = archived_path(name, data, heat_group)
        assert cr.digest(path) == sha, path
        record(path)
    chosen = heat_check['selected']
    heat_consumed = {}
    study = heat_group/chosen['study']
    binary = json.loads((study/'summary.json').read_text())['executable']
    heat = load_heat_case(study, chosen['spacing_m'], chosen['steps'], heat_consumed,
                          executable=archived_path(binary, data, heat_group))
    for name in heat_consumed:
        record(Path(name))
    key = 'ogata_banks'
    target = output/'data'/key
    target.mkdir(exist_ok=refresh)
    shutil.copy2(heat['folder']/'fresh.cfg', target/'input.cfg')
    shutil.copy2(study/'runner.py', target/'check_heat_transport.py')
    x, time = heat['coordinate'], heat['times'][-1]
    ref = heat['reference']
    grid = np.linspace(0, 1, 801)
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, (ax,) = ds.result_figure(1)
        ds.analytical(ax, grid, heat['exact_function'](grid, time, ref['velocity'], ref['diffusivity']))
        # One node per distinct x, about 40 along the column, so the analytical curve stays visible.
        _, first = np.unique(np.round(x, 9), return_index=True)
        pick = first[np.argsort(x[first])][::max(1, len(first)//40)]
        ds.des(ax, x[pick], heat['theta'][-1][pick])
        ds.panel(ax, None, f'Temperature front at {time:g} s', xlabel='distance from the hot inlet (m)', ylabel='$(T-T_0)/(T_{in}-T_0)$')
        pf.save_figure(fig, output/'figures', key)
    np.savetxt(target/'profile.csv', np.column_stack((x, heat['theta'][-1], heat['exact'][-1])),
               delimiter=',', header='x_m,DES_normalized_temperature,reference_normalized_temperature', comments='')
    error = float(max(heat['errors']))
    records.append(dict(key=key, theme='heat', title='Groundwater heat transport: Ogata–Banks',
        problem='Steady groundwater carries heat from a hot inlet along a saturated column while conduction spreads the front.',
        result=f'Largest temperature error: {100*error:.3f}% of the imposed 1 K jump, over all nodes and five saved times; the pressure stays exactly steady.',
        scope='3D, 0.003125 m mesh, 8192 steps to 0.47 s, with scaled coefficients that isolate transport. Mesh and time-step refinement studies and an exact-restart check also pass in 2D and 3D.',
        meaning='This checks Darcy heat advection and the mixed solid–fluid heat capacity: the front moves at the bulk-capacity speed, not the pore-water velocity. It does not test enthalpy conservation with evolving storage.',
        metric=error, metric_label='temperature error / imposed jump', dimension=3, spacing_m=.003125,
        evidence='runs/private-des-integration-20261107/heat-qualification.json'))

    key = 'point_sources'
    for name in ('point-source-sha256.json', 'point-source.diff', 'point-implementation-commit.txt'):
        if (heat_group/name).exists():record(heat_group/name)
    style_group = heat_group/'style'
    folder = style_group/'points-3d'
    for name in ('source-sha256.json', 'parent.json', 'verification.json'):
        record(style_group/name)
    point_check = json.loads((folder/'summary.json').read_text())
    assert len(point_check['checks']) == 22
    assert cr.digest(folder/'runner.py') == point_check['runner_sha256']
    record(folder/'summary.json'); record(folder/'runner.py')
    spec = importlib.util.spec_from_file_location('frozen_point_reference', folder/'runner.py')
    sys.path.insert(0, str(pf.ROOT/'tests/functional'))
    point_reference = importlib.util.module_from_spec(spec); spec.loader.exec_module(point_reference)
    selected_point = folder/'mixed-capacity'
    cfg = configparser.ConfigParser(); cfg.read(selected_point/'fresh.cfg')
    point = point_reference.check_balance(selected_point, cfg)
    saved_point = point_check['checks']['mixed-capacity']['balance']
    assert len(point) == len(saved_point)
    for actual, saved in zip(point, saved_point):
        assert actual.keys() == saved.keys()
        # BLAS reduction order can change the final bit of a replayed summary.
        # The frozen check_balance still enforces the physical budget bounds.
        np.testing.assert_array_max_ulp(list(actual.values()),
                                       [saved[key] for key in actual], maxulp=2)
    assert not point_check['checks']['mixed-capacity']['restart']['legacy_output_clock_rounding']
    assert point_check['checks']['remesh']['all_continuation_payloads_exact']
    target = output/'data'/key; target.mkdir(exist_ok=refresh)
    shutil.copy2(selected_point/'fresh.cfg', target/'input.cfg')
    shutil.copy2(folder/'runner.py', target/'check_point_sources.py')
    record(selected_point/'fresh.cfg')
    for path in selected_point.glob('fresh.*'):
        if '.save.' in path.name or '.chkpt.' in path.name or path.suffix == '.info':record(path)
    time = np.array([r['time_s'] for r in point])
    fluid = np.array([r['fluid_m3'] for r in point])
    energy = np.array([r['heat_J'] for r in point])
    exact_fluid = np.array([r['expected_fluid_m3'] for r in point])
    exact_heat = np.array([r['expected_heat_J'] for r in point])
    line_times = np.array([0,.25,.8,1,2,2.1,2.75,3,4])
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(2)
        for ax, letter, name, values, rates, scale, label in (
            (axes[0],'a','Stored water',fluid,[8e-8,-3e-8],1e9,'net supplied volume (µL)'),
            (axes[1],'b','Stored heat',energy,[4000,-1500],1e-3,'net supplied heat (kJ)')):
            ds.analytical(ax,line_times,[scale*point_reference.amount(t,rates) for t in line_times],label='schedule integral')
            ds.des(ax,time,scale*values,size=5)
            ds.panel(ax,letter,name,xlabel='$t$ (s)',ylabel=label)
        pf.save_figure(fig,output/'figures',key)
    np.savetxt(target/'history.csv',np.column_stack((time,fluid,exact_fluid,energy,exact_heat)),
        delimiter=',',header='time_s,DES_fluid_m3,exact_fluid_m3,DES_heat_J,exact_heat_J',comments='')
    ferror = float(np.max(np.abs(fluid-exact_fluid))/max(abs(exact_fluid)))
    herror = float(np.max(np.abs(energy-exact_heat))/max(abs(exact_heat)))
    records.append(dict(key=key,theme='groundwater',title='Scheduled water and heat sources in a closed box',
        problem='One point supplies water and heat while another extracts them, over time windows that end inside time steps.',
        result=f'Largest stored-volume error: ${sci(ferror)}$ and largest stored-heat error: ${sci(herror)}$ of the final net supply.',
        scope='3D closed box, 0.25 m mesh, 0–4 s, fixed mesh and constant capacity. Further checks pass in 2D and 3D: signed schedules, thermal substeps, a source at a well, unit conversion of point coordinates and restart after remeshing. (a) stored water; (b) stored heat.',
        meaning='This checks that every source acts exactly once over each integration interval. It does not test the spatial point-source solution.',
        metric=max(ferror,herror),metric_label='stored amount error / final net supply',dimension=3,spacing_m=.25,
        evidence='runs/private-des-integration-20261107/style/points-3d/summary.json'))

    buoyancy = heat_group/'buoyancy'
    for spacing in ('0.05', '0.025'):
        folder = buoyancy/f'hrl-{spacing}'
        assert folder.is_dir(), folder
    records.append(build_convection(buoyancy/'hrl-0.025', output, record,
                                    coarse=buoyancy/'hrl-0.05'))
    convection_3d = data/'runs/convection-3d-20261110/study'
    if (convection_3d/'qualification.json').exists():
        records.append(build_convection_3d(convection_3d, output, record))
    elder = data/'runs/elder-20261110'
    if (elder/'qualification.json').exists():
        records.append(build_elder(elder, output, record))
    for name in ('incoming-source.json', 'review-validation.json', 'fine-replay.json',
                 'qualification.json'):
        record(buoyancy/name)

    # Poroelastic benchmarks: Theis, Mandel (rigid platen), Cryer (curved load), thin-slab aquifer.
    poro = data/'runs/private-des-integration-20261107/platen/studies'
    records.append(qualify_theis.build(poro/'theis', output, record))
    records.append(qualify_mandel.build(poro/'mandel-t800', output, record, poro/'mandel-h',
                                        [poro/'mandel-t100', poro/'mandel-t400']))
    cryer_times = [poro/'cryer-t400'] if (poro/'cryer-t400'/'summary.json').exists() else []
    records.append(qualify_cryer.build(poro/'cryer-h', output, record, poro/'cryer-h', cryer_times))
    records.append(qualify_aquifer.build(poro/'slabref', output, record))
    rud = sorted(poro.glob('rudnicki-h*-t48'), key=lambda f: float(f.name.split('-h')[1].split('-t')[0]))
    records.append(qualify_rudnicki.build(rud[0], output, record, rud,
                                          sorted(poro.glob('rudnicki-h*-t96'))))
    # Leaky aquifer and Bishop column (U8); built once their archived studies exist.
    bishop = data/'runs/bishop-20261108'
    if (bishop/'hantush'/'runner.py').exists():
        records.append(qualify_hantush.build(bishop/'hantush', output, record))
    if all((bishop/name/'runner.py').exists() for name in ('u8-2d', 'u8-3d')):
        records.append(qualify_bishop.build(bishop/'u8-2d', output, record, bishop/'u8-3d'))
    cell = data/'runs/bishop-transient-20261108'
    if all((cell/name/'runner.py').exists() for name in ('cell-2d', 'cell-3d', 'cell-before-2d')):
        records.append(qualify_bishop_cell.build(cell/'cell-2d', output, record, cell/'cell-3d',
                                                 cell/'cell-before-2d', cell/'restart-2d'))
    slope = data/'runs/slope-application-20261108'
    # Salençon: 16 unloading steps with mixed stress on three graded meshes; accepted for the manual with its limits stated.
    salencon = data/'runs/salencon-audit-20261108/refinement-replay'
    if (salencon/'qualification.json').exists():
        records.append(qualify_salencon.build(salencon, output, record))
    srm_study = data/'runs/strength-reduction-20261108/explicit'
    if (srm_study/'runner.py').exists():
        records.append(qualify_embankment.build(srm_study, output, record))
    drawdown = data/'runs/drawdown-20261108/study'
    fine_drawdown = data/'runs/drawdown-20261110/study'
    if (fine_drawdown/'audit.json').exists():
        drawdown = fine_drawdown
    if (drawdown/'runner.py').exists():
        records.append(qualify_drawdown.build(drawdown, output, record))
    # -6354: rerun on the binary with the Mohr-Coulomb edge-return fix; the 90cfdaa studies stopped PT earlier.
    if all((slope/name/'runner.py').exists() for name in ('mesh-6354', 'steps808-6354')):
        records.append(qualify_slope.build(slope/'mesh-6354', output, record, slope/'steps808-6354'))
    rain = data/'runs/rainfall-slope-20261108/mesh-6354'
    if (rain/'runner.py').exists() and (slope/'mesh-phi30-6354'/'runner.py').exists():
        records.append(qualify_rainfall.build(rain, output, record, slope/'mesh-phi30-6354'))
    rsf = data/'runs/rsf-step-20261108/steps-2d'
    if (rsf/'runner.py').exists():
        records.append(qualify_rsf.build(rsf, output, record))
    borehole = data/'runs/borehole-20261108'
    if all((borehole/name/'runner.py').exists() for name in ('mesh-2d', 'steps576-2d')):
        records.append(qualify_borehole.build(borehole/'mesh-2d', output, record, borehole/'steps576-2d'))
    hot = data/'runs/injection-temperature-20261108'
    if all((hot/name/'runner.py').exists() for name in ('mesh-2d', 'mesh-coarse-2d', 'steps200-2d', 'steps1000-2d', 'doublet-mesh')):
        records.append(qualify_injection_temperature.build(hot/'mesh-2d', output, record, hot/'mesh-coarse-2d',
                                                           [hot/'steps200-2d', hot/'steps1000-2d'], hot/'doublet-mesh',
                                                           [f for f in (hot/'sanity-3d',
                                                            data/'runs/coupled-3d-20261109/model-avdonin-3d/results') if (f/'runner.py').exists()]))
    drainage = data/'runs/drainage-20261108'
    if all((drainage/name/'runner.py').exists() for name in ('mesh-2d', 'steps960-2d')):
        sanity = data/'runs/coupled-3d-20261109/model-u9-3d/results'
        if not (sanity/'runner.py').exists(): sanity = drainage/'sanity-3d'
        explicit = drainage/'explicit-2d'
        records.append(qualify_drainage.build(drainage/'mesh-2d', output, record, drainage/'steps960-2d',
                                              sanity if (sanity/'runner.py').exists() else None,
                                              explicit if (explicit/'runner.py').exists() else None))
    yielding = data/'runs/pressure-yield-20261108'
    if all((yielding/name/'runner.py').exists() for name in ('mesh-2d', 'steps1600-2d')):
        extra = [yielding/name for name in ('steps6400-2d', 'pt1-steps1600-2d', 'pt1-steps6400-2d') if (yielding/name/'runner.py').exists()]
        slabs = [data/f'runs/coupled-3d-20261109/model-yield-{mode}-3d/results' for mode in ('pt', 'explicit', 'gvs')]
        records.append(qualify_pressure_yield.build(yielding/'mesh-2d', output, record, yielding/'steps1600-2d', extra,
                                                   [p for p in slabs if (p/'runner.py').exists()]))
    nusselt = data/'runs/convection-nusselt-20261108/study'
    if (nusselt/'runner.py').exists():
        records.append(qualify_nusselt.build(nusselt, output, record))
    regional = data/'runs/regional-subsidence-20261108'
    if all((regional/name/'runner.py').exists() for name in ('study', 'steps416', 'final-h25-s416')):
        records.append(qualify_regional.build(regional/'study', output, record, regional/'steps416',
                                              regional/'final-h25-s416'))
    mctigue = data/'runs/mctigue-20261108'
    if all((mctigue/name/'runner.py').exists() for name in ('mesh-2d', 'steps400-2d', 'mesh-3d')):
        records.append(qualify_mctigue.build(mctigue/'mesh-2d', output, record, mctigue/'steps400-2d', mctigue/'mesh-3d'))
    transport = data/'runs/transport-qualification-20261108'
    records.extend(qualify_transport.build(transport, output, record))
    if all((transport/name/'summary.json').exists() for name in
           ('anisotropic-time', 'anisotropic-time-steps40k', 'anisotropic-time-steps80k')):
        records.append(qualify_transport.build_anisotropic(transport, output, record))
    records.append(qualify_kirsch.build(data/'runs/application-qualification-20261108/kirsch', output, record))
    records.append(qualify_strength.build(data/'runs/ep-pressure-20261108', output, record))
    for name, build_case in (('fluid-content-20261108', qualify_content.build_water),
                             ('thermal-cell-20261108', qualify_content.build_thermal),
                             ('plastic-source-20261108', qualify_content.build_plastic)):
        folder = data/'runs'/name
        if (folder/'qualification.json').exists():
            records.append(build_case(folder, output, record))
    finite = data/'runs/finite-content-20261109'
    if (data/'runs/finite-content-20261110/qualification.json').exists():
        finite = data/'runs/finite-content-20261110'
    if (finite/'qualification.json').exists():
        for builder in (qualify_content.build_mass, qualify_content.build_drained_heat, qualify_content.build_mixing):
            records.append(builder(finite, output, record))
    plastic_heat = data/'runs/shear-energy-20261109'
    if (plastic_heat/'qualification.json').exists():
        records.append(qualify_plastic_heat.build(plastic_heat, output, record))
    fault_heat = data/'runs/fault-heat-20261109'
    if (fault_heat/'qualification.json').exists():
        records.append(qualify_fault_heat.build(fault_heat, output, record))
    supplements = {
        'terzaghi': [('Partial unloading and material regimes', 'terzaghi-pulse.md')],
        'mass': [('Explicit inventory model and limits', 'mass-model.md')],
        'drained_heat': [('Explicit inventory model and limits', 'mass-model.md')],
        'mixing': [('Explicit inventory model and limits', 'mass-model.md')],
        'undrained_wave': [('Switching from explicit/GVS to physical inertia', 'physical-inertia.md')],
        'slope': [('Current-source coupling and fluid-content checks', 'slope-audit.md')],
        'rainfall': [('Current-source coupling and fluid-content checks', 'rainfall-audit.md')],
        'plastic_heat': [('Implementation and energy ownership', 'plastic-heat-implementation.md')],
        'fault_heat': [('Reference, method and transition checks', 'fault-heat-reference.md')],
    }
    for rec in records:
        rec['supplements'] = supplements.get(rec['key'], [])
        if not (output/'data'/rec['key']/'input.cfg').is_file():
            raise RuntimeError(f"{rec['key']}: the problem builder must export data/{rec['key']}/input.cfg")
    input_changes = portable_inputs(data, output, consumed, record)
    records.sort(key=lambda rec: ORDER.index(rec['key']))
    for note in (output/'data').glob('*/NOTE-*.txt'): note.unlink()   # study notes stay with the study; hashes in provenance
    write_review(records, output)

    evidence = ['<!doctype html><meta charset="utf-8"><title>Supporting verification evidence</title>',
                '<h1>Supporting verification evidence</h1><p><a href="index.html">Return to the selected results</a></p>',
                '<p>The transport integration passes 76 recorded checks, both CPU dimensions and exact parent comparisons of 208 save/checkpoint files and 6,344 physical payloads. The physical-inertia correction adds separate wave mesh/time sequences, exact physical restart and automatic-step checks. Its 94-check regression pass preserves 360 fixed-stress save/checkpoint files and 10,714 physical payloads exactly. These checks do not complete Bishop mechanics, GVS earthquake cycles or applications.</p>',
                '<p>Selection uses the predefined finest-grid benchmark snapshot, without fitting material parameters. Raw trials, including failures, remain in the local evidence archive. The input decks and compact profiles in this review can be used with their recorded source revisions; full data hashes are in provenance.json.</p>']
    for rec in records:
        evidence += [f'<h2 id="{rec["key"]}">{html.escape(rec["title"])}</h2><p>{html.escape(rec["result"])}</p><p>{html.escape(rec["scope"])}</p>',
                     '<p>Local source record: <code>'+html.escape(rec['evidence'])+'</code></p>']
    (output / 'evidence.html').write_text('\n'.join(evidence)+'\n')
    cr.write_json(output / 'provenance.json', dict(selection_rule='Finest qualified grid, predefined target time, no parameter fitting.',
        analyzer_sha256=cr.digest(Path(__file__)),
        dependency_sha256={name: cr.digest(Path(__file__).with_name(name)) for name in
                           ('prepare_figures.py', 'consolidation_refinement.py',
                            'terzaghi_reference.py', 'unsaturated_reference.py', 'qualify_wave.py',
                            'qualify_heat.py', 'qualify_convection.py', 'convection_reference.py', 'elder_reference.py', 'plot_elder.py',
                            'qualify_theis.py', 'theis_reference.py', 'qualify_mandel.py', 'mandel_reference.py',
                            'qualify_cryer.py', 'cryer_reference.py', 'qualify_aquifer.py', 'pumping_slab.py',
                            'slab_refinement.py', 'averaged_aquifer_reference.py',
                            'qualify_rudnicki.py', 'rudnicki_reference.py',
                            'qualify_hantush.py', 'hantush_reference.py', 'qualify_bishop.py',
                            'qualify_mctigue.py', 'mctigue_reference.py', 'qualify_bishop_cell.py', 'qualify_nusselt.py', 'convection_steady.py', 'qualify_pressure_yield.py', 'pressure_yield_reference.py', 'qualify_drainage.py', 'drainage_reference.py', 'qualify_injection_temperature.py', 'avdonin_reference.py', 'qualify_borehole.py', 'borehole_reference.py', 'qualify_rsf.py', 'qualify_salencon.py', 'salencon_reference.py', 'qualify_embankment.py', 'strength_reduction_slope.py', 'qualify_drawdown.py', 'drawdown_slope.py', 'bishop_reference.py', 'qualify_slope.py', 'slope_application.py', 'qualify_rainfall.py', 'rainfall_slope.py', 'qualify_regional.py', 'regional_subsidence.py',
                            'qualify_content.py', 'fluid_content.py', 'thermal_cell.py', 'mass_reference.py',
                            'qualify_plastic_heat.py', 'plastic_heat.py', 'qualify_fault_heat.py', 'fault_heat.py', 'manual_pages.py', 'des_options.py',
                            'qualify_transport.py', 'qualify_kirsch.py', 'kirsch_reference.py', 'qualify_strength.py',
                            'diagram_style.py', 'setup_drawings.py', 'natural_math.py', 'review.py')},
        typesetter_sha256={name:cr.digest(pf.ROOT/'third_party/katex'/name) for name in
                           ('package.json', 'dist/katex.js', 'dist/katex.min.css')},
        reader_sha256=cr.digest(pf.ROOT / 'Dynearthsol.py'),
        portable_inputs=input_changes,
        functional_source_sha256={str(Path(module.__file__).resolve().relative_to(pf.ROOT)): cr.digest(Path(module.__file__))
            for module in tuple(sys.modules.values()) if getattr(module, '__file__', None)
            and Path(module.__file__).resolve().is_relative_to(pf.ROOT/'tests/functional')
            and Path(module.__file__).suffix == '.py'},
        restart_reader_sha256=cr.digest(pf.ROOT/'tests/functional/pt_initial_checkpoint.py'),
        review_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=pf.ROOT, text=True).strip(),
        documentation_sha256={'doc/hydro/dynamics.md':cr.digest(dynamic_doc),
                              'doc/hydro/mechanics.md':cr.digest(mechanics_doc),
                              'doc/hydro/aquifer.md':cr.digest(aquifer_doc),
                              **{'doc/hydro/'+name:cr.digest(pf.ROOT/'doc/hydro'/name)
                                 for name in supporting_notes.values()}},
        physical_inertia_commit=(group/'qualification-commit.txt').read_text().strip(),
        heat_transport_commit=(heat_group/'heat-implementation-commit.txt').read_text().strip(),
        point_sources_commit=((heat_group/'point-implementation-commit.txt').read_text().strip()
                              if (heat_group/'point-implementation-commit.txt').exists() else None),
        reference_only_cases=list(REFERENCE_ONLY),
        reference_pages={key:SLUGS[key]+'.html' for key in REFERENCE_ONLY},
        source_files=consumed, integration_commit=check.get('integration_commit'), cases=records))
    typeset(output, pf.ROOT/'third_party/katex')
    return records


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--refresh', action='store_true', help='Regenerate an existing review made by this script')
    parser.add_argument('--help-dir', type=Path, help='Captured des-help2d.txt and des-help3d.txt from qualified executables')
    args = parser.parse_args()
    if args.help_dir:
        load_catalogues(args.help_dir)
    result = build(args.data.resolve(), args.output.resolve(), args.refresh)
    print(json.dumps(result, indent=2))
