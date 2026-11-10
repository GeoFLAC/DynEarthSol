#!/usr/bin/env python3
"""Sealed direct-shear fault patch at prescribed total normal traction.

An independent differential Mohr-Coulomb/Jaumann reference includes normal
compliance, plastic work and fluid expansion. This uniform patch is a local
fault-zone model; it has no rupture front, heat loss or fluid drainage.
"""
import argparse
import copy
import json
import math
from pathlib import Path
import shutil

import numpy as np
from scipy.integrate import solve_ivp

import plastic_heat as heat

HEIGHT, NORMAL = .01, 5e7
BULK, SHEAR, COHESION = 2e10, 1e10, 1e6
ANGLE, BIOT, BETA = math.pi/6, .8, 4e-5
STORAGE, CAPACITY = .2/2.2e9+.6/1e11, 2996000.
LAME = BULK-2*SHEAR/3
TAU0 = NORMAL*math.sin(ANGLE)+COHESION*math.cos(ANGLE)


import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tests/functional/hydromechanics'))
from des_case import source_revision

def reference(gamma, heating=True, expansion=True, rtol=1e-10):
    """Integrate in shear strain, starting at yield under normal stress -N.

    Unknowns: sigma_xx, sigma_xz, p, temperature rise, volume strain, plastic
    work per volume. Normal equilibrium and yield consistency determine the
    normal strain rate and plastic multiplier in a two-by-two linear system.
    The out-of-plane principal must remain intermediate on this smooth face.
    """
    gamma = np.asarray(gamma)
    beta = BETA if expansion else 0.

    def rhs(strain, state):
        x, tau, p, temperature, epsilon, work = state
        d = (x+NORMAL)/2
        q = math.hypot(d, tau)
        a, s = d/q, math.sin(ANGLE)
        # p'=(beta*chi*lambda*q/C-alpha*epsilon')/S.
        thermal = BIOT*beta*int(heating)*q/(STORAGE*CAPACITY)
        ax = LAME+BIOT**2/STORAGE
        bx = -SHEAR*a-thermal
        fe = .5*(a+s)*ax-s*BIOT**2/STORAGE
        fl = .5*(a+s)*bx-SHEAR*(tau/q)**2+s*thermal
        fg = .5*(a+s)*tau+tau/q*(SHEAR-d)
        epsilon_rate, plastic = np.linalg.solve(
            [[LAME+2*SHEAR+BIOT**2/STORAGE, SHEAR*a-thermal], [fe, fl]],
            [tau, -fg])
        assert plastic >= -1e-10, 'Reference leaves the assumed loading face'
        dissipation = plastic*q
        pdot = (beta*int(heating)*dissipation/CAPACITY-BIOT*epsilon_rate)/STORAGE
        return [LAME*epsilon_rate-SHEAR*plastic*a+tau-BIOT*pdot,
                SHEAR-SHEAR*plastic*tau/q-d, pdot,
                int(heating)*dissipation/CAPACITY, epsilon_rate, dissipation]

    initial = [-NORMAL, TAU0, 0., 0., 0., 0.]
    if gamma[-1] == 0: return np.tile(initial, (len(gamma), 1)).T
    sol = solve_ivp(rhs, (0., gamma[-1]), initial, dense_output=True,
                    rtol=rtol, atol=[1e-4, 1e-4, 1e-4, 1e-10, 1e-12, 1e-4])
    assert sol.success, sol.message
    result = sol.sol(gamma)
    x, tau, p, temperature, epsilon, work = result
    d, m = (x+NORMAL)/2, (x-NORMAL)/2
    q = np.hypot(d, tau)
    syy = -NORMAL+LAME*epsilon-BIOT*p
    assert np.all((syy >= m-q-1e-3) & (syy <= m+q+1e-3)), 'Reference crosses a principal edge'
    np.testing.assert_allclose(q+math.sin(ANGLE)*(m+BIOT*p)-COHESION*math.cos(ANGLE), 0., atol=.1)
    np.testing.assert_allclose(STORAGE*p+BIOT*epsilon-beta*temperature, 0., atol=1e-11)
    return result


def configuration(mode, dt, strain, velocity, scaling=1e5, heating=True, expansion=True, samples=500, dt_fraction=.2):
    cfg = heat.configuration(2, 'pt1' if mode == 'pt' else mode, dt, 1)
    duration = strain*HEIGHT/velocity
    steps = round(duration/dt) if mode == 'pt' else 100000000
    cfg['sim'].update(max_steps=str(steps), max_time_in_yr=repr(duration/heat.water.YEAR2SEC),
                      output_step_interval='1' if mode == 'pt' else '100000000',
                      output_time_interval_in_yr=repr(1e9 if mode == 'pt' else duration/samples/heat.water.YEAR2SEC),
                      info_display_step_interval='100000')
    cfg['mesh'].update(xlength=repr(HEIGHT), zlength=repr(HEIGHT), resolution=repr(HEIGHT),
                       largest_size=repr(HEIGHT), smallest_size=repr(HEIGHT))
    cfg['control'].update(has_shear_heating=str(heating).lower(), has_hydraulic_diffusion='true',
        has_thermal_diffusion='true', has_fluid_thermal_advection='true', has_thermal_pressurization='true',
        fixed_dt=repr(dt if mode == 'pt' else 0.), is_quasi_static=str(mode != 'physical').lower(),
        characteristic_speed='0', inertial_scaling=repr(scaling), dt_fraction=repr(dt_fraction))
    if mode == 'physical': cfg['control'].update(damping_option='0', damping_factor='0')
    cfg['bc'].update(vbc_z1='8', vbc_z0='1', vbc_val_z0='0', stress_bc_z1='4',
                     stress_val_z1=repr(-NORMAL), vbc_val_x0=repr(velocity), vbc_val_x1=repr(velocity))
    for side in ('x0', 'x1'):
        cfg['bc'].update({f'vbc_val_division_{side}_min':repr(1/3),
                          f'vbc_val_division_{side}_max':repr(2/3)})
        for k, value in enumerate((1., 2/3, 1/3, 0.)):
            cfg['bc'][f'vbc_val_{side}_ratio{k}'] = repr(value)
    cfg['ic']['initial_stress'] = repr([-NORMAL, -NORMAL, TAU0])
    cfg['mat'].update(bulk_modulus=repr([BULK]), shear_modulus=repr([SHEAR]), porosity='[.2]',
        biot_coeff=repr([BIOT]), bulk_modulus_s='[1e11]', fluid_bulk_modulus='[2.2e9]',
        hydraulic_perm='[0]', fluid_rho0='[1000]', fluid_heat_capacity='[4180]',
        fluid_alpha='[2e-4]' if expansion else '[0]', cohesion0=repr([COHESION]),
        cohesion1=repr([COHESION]), friction_angle0='[30]', friction_angle1='[30]')
    return cfg


def analyze(model, mode, velocity, heating=True, expansion=True, initial_generated=0., initial_checkpoint_model=None,
            reference_start_time=0., reference_start_gamma=0.):
    frames = heat.scan_frames(model.modelname)
    times = np.array([r['time'] for r in frames])
    x, conn = (model.read_field(frames[0]['frame'], key) for key in ('coordinate', 'connectivity'))
    volume = abs(np.linalg.det(x[conn][:, 1:]-x[conn][:, :1]))/2
    vtotal = volume.sum()
    values, generated, every_step = [], initial_generated*vtotal, True
    lag = BIOT**2/(BULK+SHEAR/3) if mode == 'pt' else 0.
    for i, row in enumerate(frames):
        frame = row['frame']
        stress = model.read_field(frame, 'stress')
        temperature = model.read_field(frame, 'temperature')[conn].mean(axis=1)-heat.T0
        pressure = model.read_field(frame, 'pore pressure')[conn].mean(axis=1)
        epsilon = model.read_field(frame, 'strain')[:, :2].sum(axis=1)
        pending = model.read_field(frame, 'pending shear heat').sum() if heating else 0.
        power = model.read_field(frame, 'shear heating') @ volume if heating else 0.
        # A non-PT restart writes its initial save frame but reuses the parent
        # checkpoint. Read that explicit source; do not invent a zero increment.
        checkpoint = initial_checkpoint_model if i == 0 and initial_checkpoint_model is not None else model
        dp = heat.water.checkpoint_field(checkpoint, frame, 'pending pore pressure increment', len(x))
        fluid = (STORAGE*pressure+BIOT*epsilon-(BETA if expansion else 0)*temperature-lag*dp[conn].mean(axis=1)) @ volume
        if i:
            consecutive = row['steps'] == frames[i-1]['steps']+1
            every_step &= consecutive
            # A frame stores the next dt. Its value is the generating interval
            # for the following accepted step, without absolute-time rounding.
            interval = frames[i-1]['dt'] if consecutive else row['time']-frames[i-1]['time']
            assert interval is not None and interval > 0
            generated += power*interval
        values.append([row['time'], *(stress.T @ volume/vtotal),
                       temperature @ volume/vtotal, pressure @ volume/vtotal,
                       epsilon @ volume/vtotal, pending/vtotal,
                       CAPACITY*temperature @ volume/vtotal+pending/vtotal,
                       generated/vtotal, fluid/vtotal])
    series = np.asarray(values)
    exact = reference(reference_start_gamma+(times-reference_start_time)*velocity/HEIGHT, heating, expansion)
    normal_error = float(np.max(abs(series[:, 2]+NORMAL))/NORMAL)
    stress_error = float(np.max(abs(series[:, 3]-exact[1]))/TAU0)
    temperature_scale = max(float(exact[3, -1]), 1.)
    temp_error = float(np.max(abs(series[:, 4]-exact[3]))/temperature_scale)
    pressure_error = float(np.max(abs(series[:, 5]-exact[2]))/(NORMAL/BIOT))
    energy_error = float(np.max(abs(series[:, 8]-int(heating)*exact[5]))/max(float(exact[5, -1]), 1.))
    budget_error = float(np.max(abs(series[:, 8]-series[:, 9]))/TAU0) if every_step else None
    assert max(abs(series[:, 10])) < 2e-11, ('fluid budget / patch volume', max(abs(series[:, 10])))
    if every_step: assert budget_error < 2e-10, ('accepted heat ledger', budget_error)
    if not heating: assert np.all(series[:, 4] == 0) and np.all(series[:, 7:10] == 0)
    assert stress_error < .01 and temp_error < .01 and pressure_error < .005 and energy_error < .01
    assert normal_error < .005
    columns = np.c_[series, exact.T]
    prefix = Path(model.modelname)
    table = prefix.parent/('history.csv' if prefix.name == 'cell' else prefix.name+'-history.csv')
    np.savetxt(table, columns, delimiter=',',
        header='time_s,sxx_Pa,szz_Pa,sxz_Pa,temperature_rise_K,pressure_Pa,volume_strain,pending_heat_J_m3,thermal_plus_pending_J_m3,sampled_work_J_m3,fluid_budget_over_volume,reference_sxx,reference_sxz,reference_p,reference_T,reference_volume_strain,reference_plastic_work', comments='')
    return dict(frames=len(frames), steps=frames[-1]['steps'], actual_end_time_s=float(times[-1]),
        normal_stress_error_over_load=normal_error, shear_error_over_initial=stress_error,
        temperature_error_over_final_rise=temp_error, pressure_error_over_N_over_alpha=pressure_error,
        heat_error_over_reference_work=energy_error, accepted_heat_budget_error_over_tau0=budget_error,
        every_step_saved=every_step, max_split_fluid_budget_over_volume=float(max(abs(series[:, 10]))),
        final_temperature_rise_K=float(series[-1, 4]), final_pressure_Pa=float(series[-1, 5]),
        final_shear_Pa=float(series[-1, 3]), reference_final_shear_Pa=float(exact[1, -1]))


def restart_check(exe, folder, cfg, model, mode, velocity, heating, expansion):
    rows = heat.scan_frames(model.modelname)
    start, stop = rows[len(rows)//2:len(rows)//2+2]
    conn = model.read_field(start['frame'], 'connectivity')
    generated = CAPACITY*(model.read_field(start['frame'], 'temperature')[conn].mean()-heat.T0)
    if heating: generated += model.read_field(start['frame'], 'pending shear heat').sum()/HEIGHT**2
    cfg = copy.deepcopy(cfg)
    cfg['sim'].update(is_restarting='true', restarting_from_modelname=str(folder/'cell'),
        restarting_from_frame=str(start['frame']), max_steps=str(stop['steps']), max_time_in_yr='1e9',
        output_step_interval='1', output_time_interval_in_yr='1e9')
    resumed = heat.run(exe, folder, cfg, 'resumed')
    result = analyze(resumed, mode, velocity, heating, expansion, generated, model)
    fields = ('coordinate', 'connectivity', 'velocity', 'strain', 'stress', 'temperature', 'pore pressure')
    if heating: fields += ('shear heating', 'pending shear heat')
    matched = 0
    for row in heat.scan_frames(resumed.modelname):
        original = next((r for r in (start, stop) if r['steps'] == row['steps']), None)
        if original:
            assert original['time'] == row['time']
            for field in fields:
                np.testing.assert_array_equal(model.read_field(original['frame'], field), resumed.read_field(row['frame'], field))
                matched += 1
    assert matched == 2*len(fields)
    result['exact_endpoint_payloads'] = matched
    if mode == 'pt':
        cfg['control']['fixed_dt'] = repr(float(cfg['control']['fixed_dt'])*.4)
        cfg['sim']['max_steps'] = str(start['steps']+20)
        changed = heat.run(exe, folder, cfg, 'changed')
        result['changed_timestep'] = analyze(changed, mode, velocity, heating, expansion, generated, model)
    return result


def physical_transition(exe, folder, cfg, seed, velocity, scaling, target_strain, reuse_phase=False):
    """Accelerate a slowly shearing GVS patch and reselect physical mass/dt."""
    warm = heat.scan_frames(seed.modelname)[-1]
    gamma = velocity*warm['time']/HEIGHT
    assert target_strain > gamma
    rapid = copy.deepcopy(cfg)
    rapid['sim'].update(is_restarting='true', restarting_from_modelname=str(folder/'cell'),
        restarting_from_frame=str(warm['frame']), max_steps='100000000',
        max_time_in_yr=repr((warm['time']+(target_strain-gamma)*HEIGHT)/heat.water.YEAR2SEC),
        output_step_interval='100000000',
        output_time_interval_in_yr=repr((target_strain-gamma)*HEIGHT/500/heat.water.YEAR2SEC))
    rapid['control'].update(is_quasi_static='false', use_global_velocity_scaling='false',
                            fixed_dt='0', damping_option='0', damping_factor='0')
    rapid['bc'].update(vbc_val_x0='1', vbc_val_x1='1')
    quake = heat.Dynearthsol(str(folder/'quake')) if reuse_phase else heat.run(exe, folder, rapid, 'quake')
    first = heat.scan_frames(quake.modelname)[0]
    for field in ('stress', 'strain', 'pore pressure', 'temperature', 'pending shear heat'):
        np.testing.assert_array_equal(seed.read_field(warm['frame'], field), quake.read_field(first['frame'], field))
    conn = seed.read_field(warm['frame'], 'connectivity')
    generated = CAPACITY*(seed.read_field(warm['frame'], 'temperature')[conn].mean()-heat.T0)
    generated += seed.read_field(warm['frame'], 'pending shear heat').sum()/HEIGHT**2
    result = analyze(quake, 'physical', 1., initial_generated=generated,
        initial_checkpoint_model=seed, reference_start_time=warm['time'], reference_start_gamma=gamma)
    nodal_volume = np.bincount(conn.ravel(), weights=np.repeat(np.full(len(conn), HEIGHT**2/len(conn)/3), 3),
                               minlength=len(seed.read_field(warm['frame'], 'velocity')))
    density = BULK/min(velocity*scaling, math.sqrt(SHEAR/2700))**2
    v0, v1 = seed.read_field(warm['frame'], 'velocity'), quake.read_field(first['frame'], 'velocity')
    result.update(first_physical_dt_s=first['dt'], last_gvs_dt_s=warm['dt'],
        initial_state_exact_fields=5, initial_gamma=gamma,
        artificial_kinetic_J_m=float(.5*density*(nodal_volume @ np.sum(v0*v0, axis=1))),
        physical_kinetic_J_m=float(.5*2360*(nodal_volume @ np.sum(v1*v1, axis=1))))
    assert first['dt'] < warm['dt']/100
    onset_cfg = copy.deepcopy(rapid)
    onset_cfg['sim'].update(max_steps=str(warm['steps']+20), max_time_in_yr='1e9',
                             output_step_interval='1', output_time_interval_in_yr='1e9')
    onset = heat.Dynearthsol(str(folder/'onset')) if reuse_phase and (folder/'onset.info').exists() else heat.run(exe, folder, onset_cfg, 'onset')
    result['first_intervals'] = analyze(onset, 'physical', 1., initial_generated=generated,
        initial_checkpoint_model=seed, reference_start_time=warm['time'], reference_start_gamma=gamma)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mode', choices=('pt', 'physical', 'explicit', 'gvs'), default='physical')
    parser.add_argument('--dt', type=float, default=1e-5)
    parser.add_argument('--strain', type=float, default=3.)
    parser.add_argument('--velocity', type=float, default=1.)
    parser.add_argument('--scaling', type=float, default=1e5)
    parser.add_argument('--samples', type=int, default=500)
    parser.add_argument('--dt-fraction', type=float, default=.2)
    parser.add_argument('--restart-check', action='store_true')
    parser.add_argument('--physical-transition', action='store_true', help='Restart slow GVS slip at 1 m/s with physical mass.')
    parser.add_argument('--transition-strain', type=float, default=3.)
    parser.add_argument('--no-heating', action='store_true')
    parser.add_argument('--no-expansion', action='store_true')
    args = parser.parse_args()
    if args.strain <= 0 or args.velocity <= 0 or args.dt <= 0 or args.samples < 2:
        parser.error('Strain, velocity and dt must be positive; samples must be at least two.')
    if args.physical_transition and (args.mode != 'gvs' or args.no_heating or args.no_expansion):
        parser.error('Physical transition requires GVS with heating and fluid expansion.')
    args.output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(__file__, args.output/'runner.py')
    cfg = configuration(args.mode, args.dt, args.strain, args.velocity, args.scaling,
                        not args.no_heating, not args.no_expansion, args.samples, args.dt_fraction)
    model = heat.run(args.exe.resolve(), args.output, cfg)
    report = analyze(model, args.mode, args.velocity, not args.no_heating, not args.no_expansion)
    if args.restart_check:
        report['restart'] = restart_check(args.exe.resolve(), args.output.resolve(), cfg, model,
            args.mode, args.velocity, not args.no_heating, not args.no_expansion)
    if args.physical_transition:
        report['physical_transition'] = physical_transition(args.exe.resolve(), args.output.resolve(), cfg,
            model, args.velocity, args.scaling, args.transition_strain)
    report.update(passed=True, mode=args.mode, heating=not args.no_heating, expansion=not args.no_expansion,
        executable_sha256=heat.water.digest(args.exe), runner_sha256=heat.water.digest(Path(__file__)),
        working_source_commit=source_revision(heat.water.ROOT),
        scope='2D homogeneous direct shear with one normal platen degree of freedom; fixed capacity and zero flux; no rupture front.')
    heat.water.write_json(args.output/'summary.json', report)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
