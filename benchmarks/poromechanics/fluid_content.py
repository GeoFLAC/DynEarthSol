#!/usr/bin/env python3
"""Coupled fluid-content ledger with an open drain and prescribed basal inflow.

The isothermal, fixed-mesh column expands as it fills. Pressure and strain
independently measure Biot content; the PT split's pending pressure increment
is reported separately. A mixed-boundary diffusion series supplies pressure,
stored volume and cumulative outflow, without using DES's balance field.
Every accepted step must be saved. --replay verifies and reanalyses hashes.
"""

import argparse
import configparser
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames
from terzaghi_reference import TerzaghiReference

PROPS = dict(height=1., load=1., bulk=2e7, shear=1e7, porosity=.25,
             fluid_bulk=1e8, grain_bulk=1e9, biot=.8,
             permeability=1e-11, viscosity=1e-3)
REF = TerzaghiReference(**PROPS)
INFLOW = 1e-6
MOBILITY = PROPS['permeability']/PROPS['viscosity']
PRESSURE = INFLOW*PROPS['height']/MOBILITY


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def reference(depth, time, terms=512, inflow=INFLOW):
    """Pressure, mean pressure and outward drain flux for unit basal inflow.

    p(0,t)=0; p_d(H,t)=q_in/mobility; p(d,0)=0. Integrating each sine
    gives the mean; differentiating at d=0 gives the drain flux. The
    reference is quasistatic; the ledger also applies to inertial runs.
    """
    depth = np.atleast_1d(depth)
    time = np.atleast_1d(time)
    n = np.arange(terms)
    k = (2*n+1)*np.pi/(2*REF.height)
    decay = np.exp(-time[:, None]*REF.diffusivity*k*k)
    coefficients = 2*(-1.)**n/REF.height
    p = inflow/MOBILITY*(depth - (decay*coefficients/k**2) @ np.sin(k[:, None]*depth))
    mean = inflow/MOBILITY*(REF.height/2 - (decay*coefficients/k**3).sum(axis=1)/REF.height)
    flux = inflow*(1 - (decay*coefficients/k).sum(axis=1))
    p[time == 0] = 0
    mean[time == 0] = 0
    flux[time == 0] = 0
    return p, mean, flux


def reference_checks():
    """Independent PDE, endpoint and integrated volume checks of the series."""
    d, dx, dt = np.array([.2, .5, .8]), 1e-4, 1e-5
    worst = 0.
    for t in (.1, .5, 2.):
        p, mean, _ = reference(d, t)
        rate = (reference(d, t+dt)[0]-reference(d, t-dt)[0])/(2*dt)
        curvature = (reference(d+dx, t)[0]-2*p+reference(d-dx, t)[0])/dx**2
        worst = max(worst, float(np.max(np.abs(rate-REF.diffusivity*curvature))/PRESSURE))
        integral = quad(lambda x: reference(x, t)[0][0, 0], 0, REF.height, epsabs=1e-10)[0]
        np.testing.assert_allclose(integral, mean[0]*REF.height, rtol=1e-12, atol=1e-12)
        drained = quad(lambda s: reference(0, s)[2][0], 0, t, epsabs=1e-15)[0]
        np.testing.assert_allclose(REF.effective_storage*mean[0]*REF.height,
                                   INFLOW*t-drained, rtol=1e-10, atol=1e-16)
        np.testing.assert_allclose(reference(0, t)[0], 0, atol=1e-13)
        gradient = (reference(1+dx, t)[0]-reference(1-dx, t)[0])/(2*dx)
        np.testing.assert_allclose(gradient, INFLOW/MOBILITY, rtol=1e-8)
    assert worst < 1e-6, worst
    return dict(pde_relative_residual=worst, quadrature_balance=True)


def configuration(h, dt, duration, method, scaling, inflow):
    cfg = configparser.ConfigParser(); cfg.optionxform = str
    cfg.read(ROOT/'tests/functional/pt_hydrostatic_equilibrium.cfg')
    cfg['sim'].update(modelname='column', max_steps=str(round(duration/dt) if method == 'pt' else 10**8),
                      max_time_in_yr=repr(duration/YEAR2SEC), output_step_interval='1',
                      output_time_interval_in_yr='1e9', checkpoint_frame_interval='1',
                      has_initial_checkpoint='yes', info_display_step_interval='100000')
    cfg['mesh'].update(meshing_elem_shape='1', xlength=repr(2*h), ylength=repr(2*h),
                       zlength='1', resolution=repr(h), largest_size=repr(h), smallest_size=repr(h))
    cfg['control'].update(has_PT=str(method == 'pt').lower(), fixed_dt=repr(dt if method == 'pt' else 0.),
                          gravity='0', has_thermal_diffusion='false',
                          has_poroelastic_pressure_feedback='true', has_moving_mesh='false',
                          is_using_mixed_stress='true', PT_option='1', PT_max_iter='100000',
                          inertial_scaling=repr(scaling), characteristic_speed='1e-5', dt_fraction='.2',
                          use_global_velocity_scaling=str(method.startswith('gvs')).lower(),
                          mass_scaling_reference_speed='bulk' if method == 'gvs-bulk' else 'shear')
    for side in ('x0', 'x1', 'y0', 'y1', 'z0', 'z1'):
        cfg['bc'].update({f'vbc_{side}': '0' if side == 'z1' else '1',
                          f'vbc_val_{side}': '0', f'hbc_{side}': '0'})
    cfg['bc'].update(hbc_z1='2', hbc_z0='3', hbc_pressure='[0,0,0,0,0,0]',
                     hbc_flux=f'[0,0,0,0,{-inflow!r},0]', stress_val_z1='0')
    cfg['ic'].update(weakzone_option='0', has_body_force_adjustment='false',
                     pore_pressure_option='2', pore_pressure_value='0')
    names = dict(bulk='bulk_modulus', shear='shear_modulus', porosity='porosity',
                 fluid_bulk='fluid_bulk_modulus', grain_bulk='bulk_modulus_s',
                 biot='biot_coeff', permeability='hydraulic_perm', viscosity='fluid_visc')
    for key, name in names.items():
        cfg['mat'][name] = repr([PROPS[key]])
    cfg['mat'].update(rho0='[2650]', fluid_rho0='[1000]', fluid_alpha='[0]', alpha='[0]')
    return cfg


def checkpoint_field(model, frame, name, count):
    model.suffix = 'chkpt'
    model.read_header(frame)
    path = model.get_fn(frame)
    if model.format == 'hdf5':
        import h5py
        with h5py.File(path) as stream:
            result = stream[name][:]
    else:
        with open(path, 'rb') as stream:
            stream.seek(model.field_pos[name])
            result = np.fromfile(stream, dtype='<f8', count=count)
    model.suffix = 'save'
    model.read_header(frame)
    return result


def analyze(folder, method):
    cfg = configparser.ConfigParser(); cfg.read(folder/'input.cfg')
    inflow = -json.loads(cfg['bc']['hbc_flux'])[4]
    pressure_scale = inflow*REF.height/MOBILITY
    model = Dynearthsol(str(folder/'column')); dim = int(model.ndims)
    rows = scan_frames(str(folder/'column'))
    assert [row['steps'] for row in rows] == list(range(len(rows))), 'every physical step must be saved'
    x, conn = (model.read_field(0, name) for name in ('coordinate', 'connectivity'))
    vertices = np.concatenate((np.ones((*conn.shape, 1)), x[conn]), axis=2)
    volume = np.abs(np.linalg.det(vertices))/math.factorial(dim)
    grad = np.linalg.inv(vertices)[:, 1:, :].swapaxes(1, 2)
    area = volume.sum()/REF.height
    # A true exterior facet, rather than a guessed element band, measures drain flux.
    top = np.isclose(x[conn, -1], 0).sum(axis=1) == dim
    facet_area = np.zeros(len(conn))
    for e in np.flatnonzero(top):
        face = x[conn[e]][np.isclose(x[conn[e], -1], 0)]
        facet_area[e] = np.linalg.norm(face[1]-face[0]) if dim == 2 else np.linalg.norm(
            np.cross(face[1]-face[0], face[2]-face[0]))/2
    np.testing.assert_allclose(facet_area.sum(), area, rtol=1e-12)
    lag_capacity = PROPS['biot']**2/(PROPS['bulk']+(PROPS['shear']/3 if dim == 2 else 0))
    history, boundary, darcy_drain = [], 0., 0.
    previous_time, previous_flux = 0., 0.
    peak_step_error, peak_storage_error, peak_source_error = 0., 0., 0.
    old_p_integral, old_strain_integral, old_mean_integral = 0., 0., 0.
    for row in rows:
        frame, time = row['frame'], row['time']
        np.testing.assert_array_equal(model.read_field(frame, 'coordinate'), x)
        p = model.read_field(frame, 'pore pressure')
        strain = model.read_field(frame, 'strain')[:, :dim].sum(axis=1)
        stress = model.read_field(frame, 'stress')[:, :dim].mean(axis=1)
        p_integral = float(volume @ p[conn].mean(axis=1))
        strain_integral = float(volume @ strain)
        mean_integral = float(volume @ stress)
        physical = REF.storage*p_integral + PROPS['biot']*strain_integral
        dp = checkpoint_field(model, frame, 'pending pore pressure increment', len(x))
        pending = -lag_capacity*float(volume @ dp[conn].mean(axis=1)) if method == 'pt' else 0.
        algorithmic = physical + pending
        b = model.read_field(frame, 'fluid balance')
        boundary += b[3]
        flux = float(np.sum(facet_area*(-MOBILITY*np.sum(p[conn]*grad[:, :, -1], axis=1))))
        darcy_drain += .5*(time-previous_time)*(flux+previous_flux)
        exact_p, exact_mean, exact_flux = reference(-x[:, -1], time, inflow=inflow)
        exact_content = float(area*REF.height*REF.effective_storage*exact_mean[0])
        scale = max(area*inflow*time, area*inflow*1e-3)
        if frame:
            # Storage and mechanical terms are reconstructed from independent
            # states, with neither the nodal capacity nor residual copied.
            capacity = REF.storage + (lag_capacity if method == 'pt' else 0.)
            expected_storage = capacity*(p_integral-old_p_integral)
            expected_source = (PROPS['biot']*(mean_integral-old_mean_integral) /
                               (PROPS['bulk']+(PROPS['shear']/3 if dim == 2 else 0))
                               if method == 'pt' else PROPS['biot']*(strain_integral-old_strain_integral))
            peak_storage_error = max(peak_storage_error, abs(b[0]-expected_storage)/scale)
            peak_source_error = max(peak_source_error, abs(b[1]-expected_source)/scale)
            peak_step_error = max(peak_step_error, abs(b[4])/scale)
        history.append([time, physical, pending, algorithmic, boundary, exact_content,
                        area*inflow*time-boundary, darcy_drain,
                        area*inflow*time-exact_content, flux, area*exact_flux[0],
                        float(np.max(np.abs(p-exact_p[0])))/pressure_scale])
        previous_time, previous_flux = time, flux
        old_p_integral, old_strain_integral, old_mean_integral = p_integral, strain_integral, mean_integral
    history = np.asarray(history)
    scale = area*inflow*np.maximum(history[:, 0], 1e-3)
    content_error = np.abs(history[:, 3]-history[:, 4])/scale
    assert peak_step_error < 1e-10, peak_step_error
    assert max(peak_storage_error, peak_source_error) < 1e-9, (peak_storage_error, peak_source_error)
    assert content_error.max() < 1e-8, content_error.max()
    header = ('time_s,physical_content_m3,pending_split_m3,algorithmic_content_m3,boundary_net_in_m3,'
              'reference_content_m3,reaction_drain_m3,darcy_drain_m3,reference_drain_m3,'
              'darcy_drain_m3_s,reference_drain_m3_s,pressure_error_over_steady')
    np.savetxt(folder/'history.csv', history, delimiter=',', header=header, comments='')
    active = history[:, 0] >= history[-1, 0]/4
    report = dict(method=method, dimension=dim, nodes=len(x), elements=len(conn), frames=len(rows),
                  final_time_s=float(history[-1, 0]),
                  max_step_residual=float(peak_step_error), max_storage_reconstruction_error=peak_storage_error,
                  max_mechanics_reconstruction_error=peak_source_error,
                  max_content_budget_error=float(content_error.max()),
                  final_pending_over_injected=float(history[-1, 2]/scale[-1]),
                  final_physical_budget_error=float(abs(history[-1, 1]-history[-1, 4])/scale[-1]),
                  final_content_error_over_injected=float(abs(history[-1, 1]-history[-1, 5])/scale[-1]),
                  final_reaction_drain_error_over_injected=float(abs(history[-1, 6]-history[-1, 8])/scale[-1]),
                  final_darcy_drain_error_over_injected=float(abs(history[-1, 7]-history[-1, 8])/scale[-1]),
                  late_pressure_error_over_steady=float(history[active, 11].max()),
                  smallest_step_s=float(np.min(np.diff(history[:, 0]))),
                  largest_step_s=float(np.max(np.diff(history[:, 0]))),
                  measured_expansion=float(strain_integral/area))
    write_json(folder/'analysis.json', report)
    return report


def figures(output, selected):
    """One selected problem result; numerical controls remain in the tables."""
    target = output/'fig'; target.mkdir(exist_ok=True)
    folder = output/selected
    cfg = configparser.ConfigParser(); cfg.read(folder/'input.cfg')
    inflow = -json.loads(cfg['bc']['hbc_flux'])[4]
    from setup_drawings import content
    content(target, inflow)
    plot_result(folder, target, inflow)
    write_json(output/'selection.json', dict(case=selected, history_sha256=digest(folder/'history.csv'),
                                             explanation='Finest PT mesh; method and timestep controls are kept in separate tables.'))


def plot_result(folder, target, inflow=INFLOW):
    import diagram_style as ds
    import prepare_figures as pf
    model = Dynearthsol(str(folder/'column'))
    rows = scan_frames(str(folder/'column'))
    x = model.read_field(0, 'coordinate')
    history = np.loadtxt(folder/'history.csv', delimiter=',', skiprows=1)
    width = float(np.ptp(x[:, 0]))
    depth = np.linspace(0, REF.height, 301)
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(2)
        for t, color in zip((.1, .5, 1., 4.), ds.SERIES):
            row = min(rows, key=lambda r: abs(r['time']-t))
            pressure = model.read_field(row['frame'], 'pore pressure')
            ds.analytical(axes[0], reference(depth, row['time'], inflow=inflow)[0][0], depth,
                          color=color, label=f"$t={row['time']:g}$ s")
            line = np.isclose(x[:, 0], width/2)
            ds.des(axes[0], pressure[line], -x[line, -1], color=color, label='_nolegend_')
        axes[0].invert_yaxis()
        ds.panel(axes[0], 'a', 'Pressure while the column fills', xlabel='pressure (Pa)', ylabel='depth (m)')
        ds.key(axes[0], loc='upper right')
        sampled = history[::max(1, len(history)//24)]
        for exact, discrete, label, color in ((5, 1, 'stored fluid', ds.WATER), (8, 6, 'drain outflow', ds.HOT)):
            ds.analytical(axes[1], history[:, 0], history[:, exact]/width*1e6, color=color, label=label)
            ds.des(axes[1], sampled[:, 0], sampled[:, discrete]/width*1e6, color=color, label='_nolegend_')
        ds.panel(axes[1], 'b', 'Independent water balance', xlabel='time (s)', ylabel=r'volume / area ($\mu$m)')
        ds.key(axes[1])
        pf.save_figure(fig, target, 'result')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path); parser.add_argument('output', type=Path)
    parser.add_argument('--method', choices=('pt', 'explicit', 'gvs-shear', 'gvs-bulk'), default='pt')
    parser.add_argument('--spacings', type=float, nargs='+', default=[.2, .1, .05])
    parser.add_argument('--dt', type=float, default=.01)
    parser.add_argument('--duration', type=float, default=4.)
    parser.add_argument('--scaling', type=float, default=1e5)
    parser.add_argument('--inflow', type=float, default=INFLOW)
    parser.add_argument('--ndims', type=int, choices=(2, 3), default=2)
    parser.add_argument('--replay', action='store_true')
    parser.add_argument('--figures', action='store_true', help='Draw the selected PT case and setup after analysis')
    args = parser.parse_args()
    exe, out = args.executable.resolve(), args.output.resolve()
    assert args.dt > 0 and args.duration > 0 and args.inflow > 0 and min(args.spacings) > 0
    reports = []
    if args.replay:
        manifest = json.loads((out/'manifest.json').read_text())
        for name, sha in manifest['artifacts_sha256'].items():
            assert digest(out/name) == sha, name
    else:
        out.mkdir(parents=True, exist_ok=False)
        manifest = dict(command=sys.argv, executable=str(exe), executable_sha256=digest(exe),
                        source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                        sources_sha256={str(p.relative_to(ROOT)): digest(p) for p in
                                       (Path(__file__), Path(__file__).with_name('terzaghi_reference.py'), ROOT/'Dynearthsol.py')},
                        method=args.method, settings=vars(args) | dict(executable=str(exe), output=str(out)),
                        properties=PROPS, reference_checks=reference_checks(), cases=[])
        env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
        for h in args.spacings:
            case = out/f'h{h:g}'; case.mkdir()
            with (case/'input.cfg').open('w') as stream:
                configuration(h, args.dt, args.duration, args.method, args.scaling, args.inflow).write(stream)
            with (case/'run.log').open('w') as stream:
                subprocess.run([str(exe), 'input.cfg'], cwd=case, env=env, check=True,
                               stdout=stream, stderr=subprocess.STDOUT, timeout=3600)
            assert digest(exe) == manifest['executable_sha256']
            manifest['cases'].append(case.name)
            manifest['artifacts_sha256'] = {str(p.relative_to(out)): digest(p) for p in sorted(out.glob('h*/*')) if p.is_file()}
            write_json(out/'manifest.json', manifest)
            reports.append(analyze(case, args.method))
            assert reports[-1]['dimension'] == args.ndims, 'executable and requested dimension differ'
            print(json.dumps(reports[-1]), flush=True)
        manifest['artifacts_sha256'] = {str(p.relative_to(out)): digest(p) for p in sorted(out.glob('h*/*'))
                                         if p.is_file() and p.name not in ('analysis.json', 'history.csv')}
        write_json(out/'manifest.json', manifest)
    summary = dict(reference_checks=reference_checks(), cases=reports or
                   [analyze(out/name, manifest['method']) for name in manifest['cases']])
    write_json(out/'summary.json', summary)
    if args.figures:
        assert manifest['method'] == 'pt' and summary['cases'][0]['dimension'] == 2
        selected = min(manifest['cases'], key=lambda name: float(name[1:]))
        figures(out, selected)


if __name__ == '__main__':
    main()
