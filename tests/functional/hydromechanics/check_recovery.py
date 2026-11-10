#!/usr/bin/env python3
"""Finite-grain plastic pore dilation, compaction and recovery.

Prescribed affine material paths isolate constitutive history from mechanical
equilibrium. The source ledger is the current tangent Biot model; this test
does not qualify finite-density mass or retained-fluid energy conservation.
"""
import argparse
import copy
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

import check_plastic_content as plastic
from check_heat_restart import checkpoint_copy, reject

ROOT = plastic.ROOT
sys.path.insert(0, str(ROOT/'tests/functional'))
from pt_initial_checkpoint import fields

ALPHA, PHI, GRAIN, FLUID = .8, .25, 1e8, 1e8
T0, TAU = 293., .25


from des_case import reset_boundaries

def configuration(method, psi, *, sealed=False, temperature=T0, energy=0.):
    cfg = plastic.configuration(method, psi)
    cfg['sim'].update(max_steps='100', has_marker_output='true')
    cfg['mesh'].update(resolution='1', largest_size='1', smallest_size='1')
    cfg['control'].update(fixed_dt='.01', is_using_mixed_stress='false',
                          use_global_velocity_scaling=str(method == 'gvs'))
    cfg['mat'].update(biot_coeff=repr([ALPHA]), bulk_modulus_s=repr([GRAIN]),
                      fluid_bulk_modulus=repr([FLUID]), porosity_model='1',
                      permeability_model='2', hydraulic_perm='[1e-18]', pore_healing_time=repr([TAU]),
                      pore_healing_temperature=repr([T0]), pore_healing_energy=repr([energy]))
    cfg['ic'].update(temperature_option='4', temp_background=repr(temperature),
                     temp_anomaly_amplitude='0')
    if not sealed:
        reset_boundaries(cfg['bc'], hbc='2')
    return cfg


def check(model, cfg, method):
    dim = int(model.ndims)
    x, conn = (model.read_field(0, k) for k in ('coordinate', 'connectivity'))
    volume = np.abs(np.linalg.det(x[conn][:, 1:]-x[conn][:, :1]))/math.factorial(dim)
    lag = ALPHA**2/(plastic.BULK+(plastic.SHEAR/3 if dim == 2 else 0)) if method == 'pt' else 0
    tau0 = json.loads(cfg['mat']['pore_healing_time'])[0]
    energy = json.loads(cfg['mat']['pore_healing_energy'])[0]
    old = None
    max_dilation = max_compaction = 0.
    errors = dict(plastic_volume=0., flow_rule=0., porosity=0., recoverable=0., source=0., pressure=0.)
    factor = lambda a: (1+math.sin(math.radians(a)))/(1-math.sin(math.radians(a)))
    nphi = factor(20)
    npsi = factor(json.loads(cfg['mat']['dilation_angle0'])[0])
    lame = plastic.BULK-2*plastic.SHEAR/3
    for frame, time in zip(model.frames, model.time):
        data = {k: model.read_field(frame, k) for k in
                ('pore dilation', 'plastic volume increment', 'porosity', 'permeability',
                 'strain', 'stress', 'temperature', 'pore pressure', 'fluid balance')}
        data['syy'] = (plastic.water.checkpoint_field(model, frame, 'stressyy', len(conn))
                       if dim == 2 else np.zeros(len(conn)))
        data['pending'] = plastic.water.checkpoint_field(model, frame, 'pending pore pressure increment', len(x))
        data['time'] = time
        if old:
            dt = time-old['time']
            t = data['temperature'][conn].mean(axis=1)
            rate = (T0/t)/tau0*np.exp(-energy/8.3144*(1/t-1/T0))
            recovered = old['pore dilation']*(-np.expm1(-dt*rate))
            de = (data['strain']-old['strain'])[:, :dim].sum(axis=1)
            dp = (data['pore pressure']-old['pore pressure'])[conn].mean(axis=1)
            used = old['pending'] if method == 'pt' else data['pending']
            # Full three-dimensional Hooke inverse, independently using saved
            # hidden plane-strain stress, excludes consumed pressure stress.
            trace0 = old['stress'][:, :dim].sum(axis=1)+(old['syy'] if dim == 2 else 0)
            trace1 = data['stress'][:, :dim].sum(axis=1)+(data['syy'] if dim == 2 else 0)
            dep = de-(trace1-trace0-3*ALPHA*used[conn].mean(axis=1))/(3*plastic.BULK)
            computed = data['plastic volume increment']
            # Independent two-face Mohr-Coulomb compression-edge flow rule.
            strength_pressure = old['pore pressure'] if method == 'pt' else data['pore pressure']
            shift = ALPHA*(used+strength_pressure)[conn].mean(axis=1)
            axial = old['stress'][:, dim-1]+shift+(lame+2*plastic.SHEAR)*de
            lateral = old['stress'][:, 0]+shift+lame*de
            cohesion = json.loads(cfg['mat']['cohesion0'])[0]
            multiplier = np.maximum(-(axial-nphi*lateral+2*cohesion*math.sqrt(nphi))/(
                plastic.SHEAR*(2+nphi*npsi)+lame*(nphi-1)*(npsi-1)), 0.)
            if float(cfg['bc']['vbc_val_z1']) < 0:
                errors['flow_rule'] = max(errors['flow_rule'], float(abs(computed-(npsi-1)*multiplier).max()))
            errors['plastic_volume'] = max(errors['plastic_volume'], float(abs(dep-computed).max()))
            expected = np.maximum(0., old['pore dilation']-recovered+dep)
            errors['recoverable'] = max(errors['recoverable'], float(abs(expected-data['pore dilation']).max()))
            phi = ALPHA-(ALPHA-old['porosity'])*np.exp(-de-dp/GRAIN+recovered)
            errors['porosity'] = max(errors['porosity'], float(abs(phi-data['porosity']).max()))
            loss = (ALPHA-old['porosity'])*np.expm1(recovered)
            source = float(volume@(ALPHA*de-loss+lag*old['pending'][conn].mean(axis=1)))
            errors['source'] = max(errors['source'], abs(source-data['fluid balance'][1]))
            p = data['pore pressure'][conn].mean(axis=1)
            p0 = old['pore pressure'][conn].mean(axis=1)
            if cfg['bc']['hbc_z1'] == '0':
                storage = old['porosity']/FLUID+(ALPHA-old['porosity'])/GRAIN
                reference = (-ALPHA*de+loss-lag*old['pending'][conn].mean(axis=1))/(storage+lag)
                errors['pressure'] = max(errors['pressure'], float(abs((p-p0)-reference).max()))
            max_dilation = max(max_dilation, float(computed.max()))
            max_compaction = max(max_compaction, float(-computed.min()))
        np.testing.assert_allclose(data['permeability'],
            json.loads(cfg['mat']['hydraulic_perm'])[0]*(data['porosity']/PHI)**3, rtol=3e-14, atol=0)
        old = data
    for key in ('plastic_volume', 'flow_rule', 'recoverable', 'porosity', 'source'):
        assert errors[key] < 3e-13, (key, errors)
    assert errors['pressure'] < 1e-7, errors
    return dict(errors=errors, max_dilative_increment=max_dilation,
                max_compactive_increment=max_compaction, final_recoverable=float(old['pore dilation'].max()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    exe, output = args.executable.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    report = []
    for method, psi, sealed, t, energy in (
        ('pt', 10, False, T0, 0), ('pt', -10, False, T0, 0), ('pt', 0, False, T0, 0),
        ('explicit', 10, False, T0, 0), ('gvs', 10, False, T0, 0),
        ('pt', 10, True, T0, 0), ('explicit', 10, True, T0, 0),
        ('pt', 10, False, 353., 25000.)):
        folder = output/f'{method}-psi{psi:g}-sealed{int(sealed)}-T{t:g}'
        folder.mkdir()
        cfg = configuration(method, psi, sealed=sealed, temperature=t, energy=energy)
        model = plastic.run(exe, folder, cfg)
        result = check(model, cfg, method)
        assert result['max_dilative_increment'] > 1e-9 if psi > 0 else (
            result['max_compactive_increment'] > 1e-9 if psi < 0 else
            result['max_dilative_increment'] < 1e-16 and result['max_compactive_increment'] < 1e-16)
        # Actual output continuation, including all optional history payloads.
        restart = copy.deepcopy(cfg)
        restart['sim'].update(is_restarting='true', restarting_from_modelname=model.modelname,
                               restarting_from_frame='50', max_steps='100')
        resumed = plastic.run(exe, folder, restart, 'resumed')
        count = plastic.compare_continuation(folder, model.modelname, resumed)
        result.update(method=method, dilation_degrees=psi, sealed=sealed,
                      temperature=t, exact_restart_payloads=count)
        # Stop loading, raise the cohesion to exclude further yielding and
        # change dt. The surviving dilation has an exponential reference.
        restart['bc']['vbc_val_z1'] = '0'
        restart['mat'].update(cohesion0='[1e9]', cohesion1='[1e9]')
        restart['control']['fixed_dt'] = '.025'
        hold = plastic.run(exe, folder, restart, 'hold')
        initial = hold.read_field(hold.frames[0], 'pore dilation')
        effective = TAU*(t/T0)*math.exp(energy/8.3144*(1/t-1/T0))
        worst = 0.
        for frame, time in zip(hold.frames, hold.time):
            target = initial*np.exp(-(time-hold.time[0])/effective)
            worst = max(worst, float(abs(hold.read_field(frame, 'pore dilation')-target).max()))
        assert worst < 1e-14, worst
        result['changed_step_hold_dilation_error'] = worst
        report.append(result)
        print('PASS', folder.name, flush=True)
    for elastic in (False, True):
        folder = output/('elastic-extension' if elastic else 'plastic-extension')
        folder.mkdir()
        cfg = configuration('pt', 10)
        cfg['bc']['vbc_val_z1'] = '.001'
        if elastic:
            cfg['mat'].update(cohesion0='[1e9]', cohesion1='[1e9]')
        model = plastic.run(exe, folder, cfg)
        result = check(model, cfg, 'pt')
        eta = model.read_field(100, 'pore dilation')
        assert eta.max() < 1e-16 if elastic else eta.min() > 1e-5
        restart = copy.deepcopy(cfg)
        restart['sim'].update(is_restarting='true', restarting_from_modelname=model.modelname,
                               restarting_from_frame='100', max_steps='150')
        restart['bc']['vbc_val_z1'] = '0'
        restart['mat'].update(cohesion0='[1e9]', cohesion1='[1e9]')
        restart['control']['fixed_dt'] = '.025'
        hold = plastic.run(exe, folder, restart, 'hold')
        phi0, phi1 = (hold.read_field(f, 'porosity') for f in (100, 150))
        k0, k1 = (hold.read_field(f, 'permeability') for f in (100, 150))
        if elastic:
            np.testing.assert_array_equal(phi0, phi1)
        else:
            assert (phi0 > PHI).all() and (k0 > 1e-18).all()
            assert (phi1 < phi0).all() and (k1 < k0).all()
            h = eta*(1-math.exp(-(hold.time[-1]-hold.time[0])/TAU))
            np.testing.assert_allclose(phi1, ALPHA-(ALPHA-phi0)*np.exp(h), atol=3e-14, rtol=0)
        result.update(case=folder.name, initial_hold_porosity=phi0.tolist(),
                      final_hold_porosity=phi1.tolist(), elastic_porosity_preserved=elastic)
        report.append(result)
    # Real changes of connectivity must preserve continuation of all history.
    cfg = configuration('pt', 10)
    cfg['sim']['has_output_during_remeshing'] = 'true'
    cfg['mesh'].update(meshing_elem_shape='0', quality_check_step_interval='1', min_quality='.999',
                        resolution='.5', largest_size='.5', smallest_size='.25')
    cfg['control']['has_moving_mesh'] = 'true'
    source_cfg = output/'remesh.cfg'
    with source_cfg.open('w') as stream:
        cfg.write(stream)
    subprocess.run([sys.executable, str(ROOT/'tests/functional/pt_initial_checkpoint.py'),
        str(exe), str(source_cfg), str(output/'remesh'), '--restart-frame', '1',
        '--steps', '4', '--require-remesh'], check=True, timeout=240,
        env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'))
    # Partial and non-finite history cannot be reconstructed on restart.
    cfg = configuration('pt', 10)
    source = output/'pt-psi10-sealed0-T293'/'column'
    bad = output/'bad-history'
    checkpoint_copy(source, 50, bad, drop=('pore dilation',))
    cfg['sim'].update(is_restarting='true', restarting_from_modelname=str(bad),
                      restarting_from_frame='50', max_steps='1')
    reject(exe, output, cfg, 'missing-history', 'requires saved recoverable dilation')
    count = len(plastic.Dynearthsol(str(source)).read_field(50, 'pore dilation'))
    checkpoint_copy(source, 50, bad, replace={'pore dilation': np.full(count, np.nan).tobytes()})
    reject(exe, output, cfg, 'nan-history', 'Invalid recoverable pore dilation')
    (output/'summary.json').write_text(json.dumps(dict(passed=True, cases=report,
        rejected_checkpoints=2, scope='Accepted plastic history, tangent-Biot closure source and recovery; no finite-density mass/energy claim.'),
        indent=2)+'\n')


if __name__ == '__main__':
    main()
