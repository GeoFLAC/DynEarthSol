#!/usr/bin/env python3
"""Sealed Biot water budget through yielding, subcycling and a changed step.

An imposed top velocity compresses a laterally confined unit column. A
mesh-volume integral of strain and pressure, plus the PT pending increment,
must stay zero. The same integral also supplies an independent mean-pressure
recurrence.
"""

import argparse
import configparser
import copy
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT), str(ROOT/'benchmarks/poromechanics'), str(ROOT/'tests/functional')]
from Dynearthsol import Dynearthsol, scan_frames
import fluid_content as water

SPEED, DURATION, DT = .001, 1., .01
BULK, SHEAR, ALPHA, PHI, FLUID, GRAIN = 2e7, 1e7, 1., .25, 1e6, 1e20
STORAGE = PHI/FLUID+(ALPHA-PHI)/GRAIN


from des_case import digest, file_hashes, hydraulic_substeps, reset_boundaries, run_model, write_json
from des_test_utils import compare_continuation, simplex_gradients, simplex_volumes

def configuration(method, psi, elastic=False, flowing=False):
    cfg = water.configuration(.5, DT, DURATION, method, 1e5, water.INFLOW)
    cfg['sim'].update(max_steps='100', max_time_in_yr='1e9')
    cfg['control']['fixed_dt'] = repr(DT)
    reset_boundaries(cfg['bc'], vbc='1', hbc='0')
    cfg['bc'].update(hbc_flux='[0,0,0,0,0,0]', vbc_val_z1=repr(-SPEED), stress_bc_z1='0')
    cfg['mat'].update(rheology_type='elastic' if elastic else 'elasto-plastic', is_plane_strain='true',
                      biot_coeff=f'[{ALPHA}]', porosity=f'[{PHI}]', bulk_modulus=f'[{BULK}]',
                      shear_modulus=f'[{SHEAR}]', fluid_bulk_modulus=f'[{FLUID}]', bulk_modulus_s=f'[{GRAIN}]',
                      hydraulic_perm='[1e-8]' if flowing else '[0]', cohesion0='[100]', cohesion1='[100]',
                      friction_angle0='[20]', friction_angle1='[20]', dilation_angle0=f'[{psi}]',
                      dilation_angle1=f'[{psi}]', max_tension='1e9')
    return cfg


def analyze(folder, method, yielding, continuation=None, open_boundary=False):
    model = Dynearthsol(str(folder/'column'))
    dim = int(model.ndims)
    records = [(model, row) for row in scan_frames(str(folder/'column'))]
    if continuation:
        resumed, start = continuation
        records = records[:start+1] + [(resumed, row) for row in scan_frames(resumed.modelname)[1:]]
    x, conn = (model.read_field(0, key) for key in ('coordinate', 'connectivity'))
    volume = simplex_volumes(x, conn)
    lag = ALPHA**2/(BULK+(SHEAR/3 if dim == 2 else 0)) if method == 'pt' else 0.
    cfg = configparser.ConfigParser(); cfg.read(folder/'column.cfg')
    mobility = json.loads(cfg['mat']['hydraulic_perm'])[0]/json.loads(cfg['mat']['fluid_visc'])[0]
    substeps = 1
    if mobility:
        substeps = hydraulic_substeps(x, conn, STORAGE+lag, mobility, DT)
        assert substeps > 1, 'the flowing control must exercise hydraulic substeps'
    scale = SPEED*DURATION*volume.sum()
    pressure_scale = ALPHA*SPEED*DURATION/STORAGE
    history, expected_pressure = [], 0.
    cumulative_boundary, largest_step_residual = 0., 0.
    for current, row in records:
        frame, time = row['frame'], row['time']
        np.testing.assert_array_equal(current.read_field(frame, 'coordinate'), x)
        p = current.read_field(frame, 'pore pressure')
        strain = current.read_field(frame, 'strain')[:, :dim].sum(axis=1)
        dp = water.checkpoint_field(current, frame, 'pending pore pressure increment')
        strain_volume = float(volume @ strain)
        pressure_volume = float(volume @ p[conn].mean(axis=1))
        content = STORAGE*pressure_volume+ALPHA*strain_volume
        pending = -lag*float(volume @ dp[conn].mean(axis=1))
        balance = current.read_field(frame, 'fluid balance')
        cumulative_boundary += balance[3]
        largest_step_residual = max(largest_step_residual, abs(balance[4])/scale)
        if frame:
            expected_pressure = (cumulative_boundary/volume.sum()+ALPHA*SPEED*time+
                                 lag*expected_pressure)/(STORAGE+lag)
        else:
            expected_pressure = 0.
        plastic = float(np.max(current.read_field(frame, 'plastic strain')))
        history.append([time, content, pending, content+pending-cumulative_boundary,
                        pressure_volume/volume.sum(), expected_pressure, strain_volume, plastic,
                        cumulative_boundary])
    values = np.asarray(history)
    budget = float(np.max(np.abs(values[:, 3]))/scale)
    mean_pressure_error = float(np.max(np.abs(values[:, 4]-values[:, 5]))/pressure_scale)
    strain_error = float(np.max(np.abs(values[:, 6]+SPEED*values[:, 0]*volume.sum()))/scale)
    assert budget < 1e-9, ('physical content + pending reservoir', budget)
    assert mean_pressure_error < 1e-9, ('independent mean-pressure recurrence', mean_pressure_error)
    assert strain_error < 1e-9, ('imposed boundary-volume change', strain_error)
    assert largest_step_residual < 1e-10, largest_step_residual
    if open_boundary:
        assert cumulative_boundary < 0, 'compression must discharge water through the drain'
        top = np.isclose(x[:, -1], 0)
        np.testing.assert_array_equal(p[top], 0)
    else:
        assert cumulative_boundary == 0, cumulative_boundary
    if yielding:
        assert values[:, 7].max() > 1e-5, 'the budget must be checked after a plastic return'
    else:
        np.testing.assert_array_equal(values[:, 7], 0)
    destination = folder/('changed-step.csv' if continuation else 'history.csv')
    np.savetxt(destination, values, delimiter=',', comments='',
               header='time_s,physical_content_m3,pending_split_m3,budget_m3,mean_pressure_Pa,reference_mean_Pa,strain_volume_m3,max_plastic_strain,boundary_in_m3')
    return dict(dimension=dim, frames=len(records), max_budget_error_over_imposed_volume=budget,
                max_pressure_error_over_reference_scale=mean_pressure_error,
                max_imposed_volume_error=strain_error, largest_step_residual=largest_step_residual,
                max_plastic_strain=float(values[:, 7].max()), hydraulic_substeps=substeps,
                changed_step=continuation is not None, cumulative_boundary_in_m3=float(cumulative_boundary),
                pressure_reference='boundary-reaction budget' if open_boundary else 'closed-column recurrence')


def exact_continuation(exe, folder, cfg, source):
    cfg = copy.deepcopy(cfg)
    cfg['sim'].update(is_restarting='true', restarting_from_modelname=str(source),
                      restarting_from_frame='20', max_steps='24')
    resumed = run_model(exe, folder, cfg, 'resumed')
    return compare_continuation(folder, str(source), resumed)


def retention_check(exe, folder, ndims, replay=False):
    """Independent scalar content inversion for uniform Gardner wetting.

    This checks the implemented strain-driven retention rate and raw pressure
    buffer, not a finite nonlinear porosity/deformation mass law.
    """
    cfg = configuration('pt', 10)
    initial, coefficient = -5000., 2e-4
    cfg['ic']['pore_pressure_value'] = repr(initial)
    cfg['mat'].update(retention_model='[3]', retention_alpha=f'[{coefficient}]', residual_saturation='[0]')
    model = Dynearthsol(str(folder/'column')) if replay else run_model(exe, folder, cfg, 'column')
    assert int(model.ndims) == ndims
    rows = scan_frames(str(folder/'column'))
    x, conn = (model.read_field(0, key) for key in ('coordinate', 'connectivity'))
    volume = simplex_volumes(x, conn)
    lag = ALPHA**2/(BULK+(SHEAR/3 if ndims == 2 else 0))
    scale = SPEED*DURATION*volume.sum()
    reference, previous_increment, previous_time = initial, 0., 0.
    previous = np.full(len(x), initial)
    source, worst_budget, worst_pressure, plastic = 0., 0., 0., 0.
    history = []
    for row in rows:
        frame, time = row['frame'], row['time']
        p = model.read_field(frame, 'pore pressure')
        if frame:
            de = -SPEED*(time-previous_time)
            se = math.exp(coefficient*reference)
            increment = brentq(lambda d: (STORAGE+lag)*d + PHI*se*math.expm1(coefficient*d) +
                               ALPHA*se*de - lag*previous_increment, -1000, 1000, xtol=1e-12)
            reference += increment
            previous_increment = increment
            source += ALPHA*de*float(volume @ np.exp(coefficient*previous[conn]).mean(axis=1))
        dp = water.checkpoint_field(model, frame, 'pending pore pressure increment')
        content = float(volume @ (STORAGE*(p[conn].mean(axis=1)-initial) +
                                  PHI*(np.exp(coefficient*p[conn]).mean(axis=1)-math.exp(coefficient*initial))))
        pending = -lag*float(volume @ dp[conn].mean(axis=1))
        worst_budget = max(worst_budget, abs(content+source+pending)/scale)
        worst_pressure = max(worst_pressure, float(np.max(abs(p-reference)))/abs(initial))
        plastic = max(plastic, float(model.read_field(frame, 'plastic strain').max()))
        history.append([time, reference, float(p.mean()), content, source, pending])
        previous, previous_time = p, time
    assert plastic > 1e-5 and np.max(previous) < 0, (plastic, previous.max())
    assert worst_budget < 1e-8 and worst_pressure < 1e-9, (worst_budget, worst_pressure)
    np.savetxt(folder/'history.csv', history, delimiter=',', comments='',
               header='time_s,reference_pressure_Pa,mean_pressure_Pa,retention_content_m3,strain_source_m3,pending_split_m3')
    return dict(case='pt-gardner', dimension=ndims, frames=len(rows),
                max_rate_budget_error_over_imposed_volume=worst_budget,
                max_pressure_error_over_initial_scale=worst_pressure, max_plastic_strain=plastic,
                scope='Uniform nonlinear retention rate and fixed-stress pressure buffer; no finite porosity-remapping claim.')


def edge_return(axial, lateral, nphi, npsi):
    """Sum of two equal lateral flow multipliers at the compression edge."""
    lame = BULK-2*SHEAR/3
    failure = axial-nphi*lateral+200*np.sqrt(nphi)
    denominator = SHEAR*(2+nphi*npsi)+lame*(nphi-1)*(npsi-1)
    multiplier = np.maximum(-failure/denominator, 0.)
    return (axial+(2*SHEAR-lame*(npsi-1))*multiplier,
            lateral-(SHEAR*npsi+lame*(npsi-1))*multiplier)


def rsf_check(exe, folder, ndims, replay=False):
    """Evolving pressure and RSF state against an independent MC edge return."""
    cfg = configuration('pt', 0)
    cfg['mesh'].update(resolution='1', largest_size='1', smallest_size='1')
    cfg['control']['rsf_slip_rate_projection_option'] = '1'
    a, b, dc, speed0 = .01, .015, .01, .0003
    cfg['mat'].update(rheology_type='elasto-plastic-rate-state-friction', direct_a=f'[{a}]',
                      evolution_b=f'[{b}]', characteristic_distance=f'[{dc}]',
                      characteristic_velocity=f'[{speed0}]', state_var_model='1')
    model = Dynearthsol(str(folder/'column')) if replay else run_model(exe, folder, cfg, 'column')
    assert int(model.ndims) == ndims and len(model.read_field(0, 'coordinate')) == 2**ndims
    report = analyze(folder, 'pt', True)
    x, conn = (model.read_field(0, key) for key in ('coordinate', 'connectivity'))
    gradient = simplex_gradients(x, conn)
    width = 1/np.linalg.norm(gradient, axis=2).max(axis=1)
    velocity = width*SPEED if ndims == 2 else 2*width*SPEED/math.sqrt(3)
    theta = np.full(len(conn), dc/speed0)
    axial, lateral = np.zeros(len(conn)), np.zeros(len(conn))
    lame = BULK-2*SHEAR/3
    previous_pressure, reference_pressure = 0., 0.
    lag = ALPHA**2/(BULK+(SHEAR/3 if ndims == 2 else 0))
    worst_state, worst_mu, worst_stress = 0., 0., 0.
    history = []
    for row in scan_frames(str(folder/'column'))[1:]:
        frame, time = row['frame'], row['time']
        theta += (1-velocity*theta/dc)*DT
        mu = math.tan(math.radians(20))+a*np.log(velocity/speed0)+b*np.log(speed0*theta/dc)
        angle = np.arctan(mu)
        nphi = (1+np.sin(angle))/(1-np.sin(angle))
        axial -= (lame+2*SHEAR)*SPEED*DT
        lateral -= lame*SPEED*DT
        axial, lateral = edge_return(axial, lateral, nphi, 1.)
        expected = np.zeros((len(conn), ndims*(ndims+1)//2))
        expected[:, :ndims] = (lateral-ALPHA*previous_pressure)[:, None]
        expected[:, ndims-1] = axial-ALPHA*previous_pressure
        actual = model.read_field(frame, 'stress')
        worst_stress = max(worst_stress, float(np.max(abs(actual-expected)))/(BULK*SPEED*DURATION))
        if ndims == 2:
            yy = water.checkpoint_field(model, frame, 'stressyy')
            worst_stress = max(worst_stress, float(np.max(abs(yy-lateral+ALPHA*previous_pressure)))/(BULK*SPEED*DURATION))
        actual_theta = model.read_field(frame, 'friction state variable')
        actual_mu = model.read_field(frame, 'dynamic friction coefficient')
        worst_state = max(worst_state, float(np.max(abs(actual_theta/theta-1))))
        worst_mu = max(worst_mu, float(np.max(abs(actual_mu-mu))))
        reference_pressure = (ALPHA*SPEED*time+lag*previous_pressure)/(STORAGE+lag)
        history.append([time, reference_pressure, float(mu.mean()), float(theta.mean()),
                        float(axial.mean()), float(lateral.mean())])
        previous_pressure = reference_pressure
    assert max(worst_state, worst_mu, worst_stress) < 1e-9, (worst_state, worst_mu, worst_stress)
    np.savetxt(folder/'rsf-reference.csv', history, delimiter=',', comments='',
               header='time_s,pressure_reference_Pa,friction_reference,state_reference_s,effective_axial_Pa,effective_lateral_Pa')
    report.update(case='pt-rsf', rsf_state_error=worst_state, rsf_friction_error=worst_mu,
                  effective_stress_reference_error=worst_stress,
                  scope='Strain-controlled sealed material cell with evolving pressure, one aging update per step and an independent two-face stress reference.')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path); parser.add_argument('output', type=Path)
    parser.add_argument('--ndims', type=int, choices=(2, 3), default=2)
    parser.add_argument('--replay', action='store_true', help='Verify all hashes and reconstruct references without running DES')
    args = parser.parse_args()
    exe, output = args.executable.resolve(), args.output.resolve()
    if args.replay:
        report = json.loads((output/'summary.json').read_text())
        assert report['passed']
        for name, sha in report['artifacts_sha256'].items():
            assert digest(output/name) == sha, name
        for index, previous in enumerate(report['cases']):
            name = previous['case']; folder = output/name
            dim = previous['dimension']
            if name == 'pt-gardner':
                current = retention_check(exe, folder, dim, replay=True)
            elif name == 'pt-rsf':
                current = rsf_check(exe, folder, dim, replay=True)
            else:
                method = 'pt' if name.startswith('pt-') else 'explicit'
                current = dict(case=name, **analyze(folder, method, 'elastic' not in name,
                                                   open_boundary=name.endswith('-drained')))
                if isinstance(previous['changed_step'], dict):
                    resumed = Dynearthsol(str(folder/'changed/column'))
                    current['changed_step'] = analyze(folder, method, True, (resumed, 20))
                    same = Dynearthsol(str(folder/'restart/resumed'))
                    current['exact_restart_payloads'] = compare_continuation(folder/'restart', str(folder/'column'), same)
            assert current == previous, (name, current, previous)
            report['cases'][index] = current
        for name, sha in report['artifacts_sha256'].items():
            assert digest(output/name) == sha, ('replay changed an artifact', name)
        report['analyzer_sha256'] = digest(__file__)
        write_json(output/'summary.json', report)
        print('PASS: all hashes, content references and restart payloads replayed')
        return
    output.mkdir(parents=True, exist_ok=False)
    report = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(__file__), cases=[])
    for method, psi, elastic, flowing in [('pt', 0, True, False), ('pt', 0, False, False),
                                         ('pt', 10, False, False), ('pt', -10, False, False),
                                         ('explicit', 10, False, False), ('pt', 10, False, True),
                                         ('explicit', 10, False, True)]:
        name = f'{method}-psi{psi}'+('-elastic' if elastic else '')+('-flow' if flowing else '')
        folder = output/name; folder.mkdir()
        cfg = configuration(method, psi, elastic, flowing)
        model = run_model(exe, folder, cfg, 'column')
        assert int(model.ndims) == args.ndims
        item = dict(case=name, **analyze(folder, method, not elastic))
        if method == 'pt' and psi == 10 and not flowing:
            same = folder/'restart'; same.mkdir()
            item['exact_restart_payloads'] = exact_continuation(exe, same, cfg, folder/'column')
            changed = folder/'changed'; changed.mkdir()
            switched = copy.deepcopy(cfg)
            switched['control']['fixed_dt'] = '.02'
            switched['sim'].update(is_restarting='true', restarting_from_modelname=str(folder/'column'),
                                    restarting_from_frame='20', max_steps='60')
            resumed = run_model(exe, changed, switched, 'column')
            item['changed_step'] = analyze(folder, method, True, (resumed, 20))
        report['cases'].append(item)
        write_json(output/'summary.json', report)
        print(json.dumps(item), flush=True)
    retained = output/'pt-gardner'; retained.mkdir()
    report['cases'].append(retention_check(exe, retained, args.ndims))
    print(json.dumps(report['cases'][-1]), flush=True)
    rsf = output/'pt-rsf'; rsf.mkdir()
    report['cases'].append(rsf_check(exe, rsf, args.ndims))
    print(json.dumps(report['cases'][-1]), flush=True)
    for method in ('pt', 'explicit'):
        drained = output/f'{method}-drained'; drained.mkdir()
        cfg = configuration(method, 10, flowing=True)
        cfg['bc']['hbc_z1'] = '2'
        run_model(exe, drained, cfg, 'column')
        report['cases'].append(dict(case=drained.name, **analyze(drained, method, True, open_boundary=True)))
        print(json.dumps(report['cases'][-1]), flush=True)
    assert digest(exe) == report['executable_sha256']
    report['artifacts_sha256'] = file_hashes((p for p in output.rglob('*') if p.name != 'summary.json'), output)
    report['passed'] = True
    write_json(output/'summary.json', report)


if __name__ == '__main__':
    main()
