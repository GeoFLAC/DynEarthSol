#!/usr/bin/env python3
"""Measure undrained waves in a closed saturated column without PT.

The initial state has zero strain/pressure and a small sinusoidal velocity.
Zero permeability makes local fluid content constant. Compare the physical
undrained standing wave, not a quasistatic consolidation reference.
The default uses physical inertia. --content-model 1 instead tests the
optional finite-density model with a declared constant fictitious mass.
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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(executable, folder, cfg, name, timeout=120):
    with (folder / f'{name}.cfg').open('w') as stream:
        cfg.write(stream)
    with (folder / f'{name}.log').open('w') as stream:
        subprocess.run([str(executable), f'{name}.cfg'], cwd=folder, check=True,
                       stdout=stream, stderr=subprocess.STDOUT, timeout=timeout,
                       env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))


def study(executable, output, spacing, steps, amplitude, samples=4, restart_check=False,
          content_model=0, scaling=1e5, gvs=False, cycles=1, moving_mesh=False):
    output.mkdir(parents=True, exist_ok=False)
    (output / 'undrained_wave.py').write_bytes(Path(__file__).read_bytes())
    # All properties also appear in the generated input deck.
    bulk, shear, phi, alpha, fluid_bulk, grain_bulk = 1e7, 1e7, .3, 1., 2.17e9, 37e9
    density = (1-phi)*2650 + phi*1000
    characteristic_speed = 1e-6
    if content_model:
        # Keep the GVS speed at its declared floor, independently of the wave.
        # This gives a constant, known fictitious density, not physical inertia.
        if gvs and (amplitude >= characteristic_speed/10 or scaling > 1e6):
            raise ValueError('The constant-mass GVS reference requires amplitude < 1e-7 and scaling <= 1e6')
        density = bulk/(characteristic_speed*scaling)**2
    storage = phi/fluid_bulk + (alpha-phi)/grain_bulk
    modulus = bulk + 4*shear/3 + alpha*alpha/storage
    speed = math.sqrt(modulus/density)
    period = 2/speed
    summary = dict(executable=str(executable), executable_sha256=digest(executable),
                   analyzer_sha256=digest(Path(__file__)),
                   source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                   content_model=content_model, scaling=scaling if content_model else None,
                   gvs=gvs, cycles=cycles, moving_mesh=moving_mesh,
                   reference=dict(storage_Pa_inv=storage, undrained_modulus_Pa=modulus,
                                  density_kg_m3=density, wave_speed_m_s=speed, period_s=period), cases={})
    patch = subprocess.check_output(['git', 'diff', 'HEAD', '--binary'], cwd=ROOT)
    (output / 'source.diff').write_bytes(patch)
    summary['source_diff_sha256'] = hashlib.sha256(patch).hexdigest()
    snapshots = {}
    for count in steps:
        folder = output / f'steps-{count}'
        folder.mkdir()
        cfg = configparser.ConfigParser(); cfg.optionxform = str
        cfg.read(ROOT / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
        cfg['sim'].update(modelname='seed', max_steps='1', output_step_interval='1',
                          output_time_interval_in_yr='1', checkpoint_frame_interval='1',
                          has_initial_checkpoint='true', max_time_in_yr='1')
        cfg['mesh'].update(xlength=str(2*spacing), ylength=str(2*spacing), zlength='1',
                           resolution=str(spacing), meshing_elem_shape='1', remeshing_option='0')
        cfg['control'].update(has_PT='false', is_quasi_static='false', damping_option='0',
                              has_hydraulic_diffusion='true', has_thermal_diffusion='false',
                              has_moving_mesh='false', gravity='0', fixed_dt=repr(period/count),
                              is_using_mixed_stress='false', use_global_velocity_scaling='false')
        if content_model:
            cfg['control'].update(hydraulic_content_model='1',is_quasi_static='true',
                inertial_scaling=repr(scaling),characteristic_speed=repr(characteristic_speed),
                has_poroelastic_pressure_feedback='true',use_global_velocity_scaling=str(gvs),
                has_moving_mesh=str(moving_mesh))
            cfg['mat']['porosity_model'] = '1'
            if gvs: cfg['control']['fixed_dt'] = '0'
        cfg['bc'].update(vbc_x1='1', vbc_z1='1', stress_bc_z1='0', stress_val_z1='0')
        for side in ('x0', 'x1', 'y0', 'y1', 'z0', 'z1'):
            cfg['bc'][f'hbc_{side}'] = '0'
        cfg['ic']['has_body_force_adjustment'] = 'false'
        cfg['mat'].update(bulk_modulus=f'[{bulk}]', shear_modulus=f'[{shear}]', porosity=f'[{phi}]',
                          biot_coeff=f'[{alpha}]', fluid_bulk_modulus=f'[{fluid_bulk}]',
                          bulk_modulus_s=f'[{grain_bulk}]', rho0='[2650]', fluid_rho0='[1000]',
                          hydraulic_perm='[0]', retention_model='[0]', alpha='[0]')
        run(executable, folder, cfg, 'seed',900 if content_model else 120)
        cfg['control']['fixed_dt'] = repr(period/count)
        seed = Dynearthsol(str(folder / 'seed'))
        if seed.format != 'binary':
            raise RuntimeError('This checkpoint-initialization fixture currently requires binary output')
        x = seed.read_field(0, 'coordinate')
        conn = seed.read_field(0, 'connectivity')
        dim = seed.ndims
        assert not np.any(seed.read_field(0, 'pore pressure'))
        assert not np.any(seed.read_field(0, 'strain'))
        assert not np.any(seed.read_field(0, 'stress'))
        velocity = np.zeros_like(x)
        velocity[:, -1] = amplitude*np.sin(np.pi*(x[:, -1]+1))
        velocity[np.isclose(x[:, -1], -1) | np.isclose(x[:, -1], 0), -1] = 0
        seed.read_header(0)
        with Path(seed.get_fn(0)).open('r+b') as stream:
            stream.seek(seed.field_pos['velocity']); stream.write(velocity.astype(np.float64).tobytes())
        np.testing.assert_array_equal(seed.read_field(0, 'velocity'), velocity)
        cfg['sim'].update(modelname='wave', is_restarting='true', restarting_from_modelname=str(folder/'seed'),
                          restarting_from_frame='0', max_steps=str(count*cycles),
                          max_time_in_yr=repr(10*cycles*period/YEAR2SEC), output_step_interval=str(count//samples))
        run(executable, folder, cfg, 'wave',900 if content_model else 120)
        model = Dynearthsol(str(folder / 'wave'))
        rows = scan_frames(str(folder / 'wave'))
        volume = np.abs(np.linalg.det(x[conn][:, 1:]-x[conn][:, :1]))/math.factorial(dim)
        nodal_volume = np.zeros(len(x))
        np.add.at(nodal_volume, conn, np.broadcast_to((volume/(dim+1))[:, None], conn.shape))
        measured = []
        snapshots[count] = []
        initial_mass = model.read_field(0,'fluid mass') if content_model else None
        for row in rows:
            frame, time = row['frame'], row['time']
            v = model.read_field(frame, 'velocity')[:, -1]
            p = model.read_field(frame, 'pore pressure')
            strain = model.read_field(frame, 'strain')
            snapshots[count].append(np.stack((v/amplitude,p/(alpha*amplitude/(storage*speed)))))
            exact_v = amplitude*np.sin(np.pi*(x[:, -1]+1))*np.cos(2*np.pi*time/period)
            exact_p = -alpha*amplitude/(storage*speed)*np.cos(np.pi*(x[:, -1]+1))*np.sin(2*np.pi*time/period)
            eps = strain[:, :dim].sum(axis=1)
            projection = np.zeros(len(x))
            np.add.at(projection, conn, np.broadcast_to((volume*eps/(dim+1))[:, None], conn.shape))
            content = alpha*projection/nodal_volume + storage*p
            measured.append(dict(time_s=time, velocity_error_over_amplitude=float(np.max(np.abs(v-exact_v))/amplitude),
                                 pressure_error_over_amplitude=float(np.max(np.abs(p-exact_p))/(alpha*amplitude/(storage*speed))),
                                 max_local_fluid_content_error=float(np.max(np.abs(content))),
                                 local_content_error_over_strain_amplitude=float(np.max(np.abs(content))/(alpha*amplitude/speed))))
            if content_model:
                np.testing.assert_allclose(model.read_field(frame,'fluid mass'),initial_mass,rtol=2e-11,atol=1e-12)
                np.testing.assert_array_equal(model.read_field(frame,'hydraulic split buffer'),0)
                assert np.isfinite(v).all() and np.isfinite(p).all()
        summary['cases'][str(count)] = dict(dimension=dim, spacing_m=spacing, steps=count, dt_s=period/count,
                                           nodes=len(x), elements=len(conn), measured=measured,
                                           seed_save_sha256=digest(Path(seed.get_fn(0))),
                                           seed_checkpoint_sha256=digest(folder/'seed.chkpt.000000'),
                                           input_sha256=digest(folder/'wave.cfg'))
        if restart_check:
            sys.path.insert(0, str(ROOT/'tests/functional'))
            from pt_initial_checkpoint import fields
            quarter = next(r for r in rows if r['steps'] == count//4)
            cfg['sim'].update(modelname='resumed', restarting_from_modelname=str(folder/'wave'),
                              restarting_from_frame=str(quarter['frame']))
            run(executable, folder, cfg, 'resume',900 if content_model else 120)
            resumed = Dynearthsol(str(folder/'resumed'))
            original = {r['steps']:r for r in rows}
            for r in scan_frames(str(folder/'resumed')):
                fresh = original[r['steps']]
                assert r['time'] == fresh['time']
                names = ['coordinate', 'velocity', 'stress', 'strain', 'pore pressure', 'temperature']
                if r['steps'] > quarter['steps']:
                    names += ['force', 'strain-rate']
                    a = fields(folder/f'wave.chkpt.{fresh["frame"]:06d}')
                    b = fields(folder/f'resumed.chkpt.{r["frame"]:06d}')
                    keys = ('pending pore pressure increment', 'dt')
                    keys += ('hydraulic pressure consumed','hydraulic coupling version') if content_model else ('physical hydraulic step',)
                    for key in keys:
                        assert a[key] == b[key], key
                for name in names:
                    np.testing.assert_array_equal(model.read_field(fresh['frame'],name),
                                                  resumed.read_field(r['frame'],name))
            cfg['sim']['modelname']='wrong-order'
            if content_model:
                cfg['control']['has_PT']='true'
            else:
                cfg['control']['is_quasi_static']='true'
                cfg['control']['has_poroelastic_pressure_feedback']='false'
            with (folder/'wrong-order.cfg').open('w') as stream:cfg.write(stream)
            with (folder/'wrong-order.log').open('w') as stream:
                rejected=subprocess.run([str(executable),'wrong-order.cfg'],cwd=folder,
                    stdout=stream,stderr=subprocess.STDOUT,timeout=60)
            assert rejected.returncode == 22
            assert ('finite-density' if content_model else 'hydraulic step ordering differs') in (folder/'wrong-order.log').read_text().lower()
            summary['cases'][str(count)].update(restart_physical_fields_exact=True,
                                              incompatible_ordering_rejected=True)
        summary['cases'][str(count)]['outputs_sha256'] = {
            p.name: digest(p) for p in sorted(folder.glob('wave.*'))
            if p.suffix == '.info' or '.save.' in p.name or '.chkpt.' in p.name}
        (output / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False))
        print(count, 'maximum errors:', {
            key: max(r[key] for r in measured) for key in measured[0] if key != 'time_s'}, flush=True)
    if restart_check and not content_model:
        # Verify compute_dt through a fresh, stationary solve, independently
        # measuring the smallest simplex height from its output coordinates.
        cfg['sim'].update(modelname='automatic', is_restarting='false', max_steps='1',
                          output_step_interval='1')
        cfg['control'].update(is_quasi_static='false', fixed_dt='0', dt_fraction='.8')
        run(executable, folder, cfg, 'automatic')
        auto = Dynearthsol(str(folder/'automatic'))
        points = auto.read_field(0,'coordinate')[auto.read_field(0,'connectivity')]
        volumes = np.abs(np.linalg.det(points[:,1:]-points[:,:1]))/math.factorial(dim)
        measures = []
        for omit in range(dim+1):
            face = np.delete(points,omit,axis=1)
            if dim == 2:
                measures.append(np.linalg.norm(face[:,1]-face[:,0],axis=1))
            else:
                measures.append(np.linalg.norm(np.cross(face[:,1]-face[:,0],face[:,2]-face[:,0]),axis=1)/2)
        height = np.min(dim*volumes/np.max(measures,axis=0))
        expected = .8*height/(5*speed)
        actual = scan_frames(str(folder/'automatic'))[1]['time']
        np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=0)
        assert actual < .8*.5*height/math.sqrt(bulk/2650)
        # Emulate an old step-zero checkpoint: no ordering marker and an
        # obsolete next-step dt. Conversion must select the new wave limit.
        (folder/'legacy-seed.save.000000').write_bytes((folder/'seed.save.000000').read_bytes())
        checkpoint = bytearray((folder/'seed.chkpt.000000').read_bytes())
        entries = dict(line.split(b'\t') for line in bytes(checkpoint[:4096]).split(b'\0',1)[0].splitlines()[1:])
        position = int(entries[b'dt'])
        checkpoint[position:position+8] = np.float64(1).tobytes()
        marker = b'physical hydraulic step\t'
        assert marker in checkpoint[:4096]
        checkpoint[:4096] = checkpoint[:4096].replace(marker,b'ignored_ hydraulic step\t',1)
        (folder/'legacy-seed.chkpt.000000').write_bytes(checkpoint)
        cfg['sim'].update(modelname='converted', is_restarting='true',
                          restarting_from_modelname=str(folder/'legacy-seed'), restarting_from_frame='0')
        run(executable,folder,cfg,'converted')
        converted_step = scan_frames(str(folder/'converted'))[1]['time']
        np.testing.assert_allclose(converted_step, expected, rtol=2e-14, atol=0)
        summary['automatic_step'] = dict(expected_s=float(expected), actual_s=actual,
            matches_undrained_wave_limit=True, smaller_than_drained_limit=True,
            initial_ordering_conversion_reselects_step=True,
            input_sha256=digest(folder/'automatic.cfg'))
        (output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
    assert digest(executable) == summary['executable_sha256']
    if content_model:
        counts = sorted(steps)
        differences = [float(np.max(np.abs(np.array(snapshots[a])-np.array(snapshots[b]))))
                       for a,b in zip(counts,counts[1:])]
        if len(differences) > 1 and all(b == 2*a for a,b in zip(counts,counts[1:])):
            assert all(a > 1.5*b for a,b in zip(differences,differences[1:])),differences
        finest = summary['cases'][str(max(steps))]['measured']
        assert max(row['velocity_error_over_amplitude'] for row in finest) < .1
        assert max(row['pressure_error_over_amplitude'] for row in finest) < .1
        summary['passed'] = True
        summary['temporal_successive_differences'] = differences
        summary['scope'] = 'Small-amplitude constant-scaled-mass undrained wave; no damping or parameter fit'
        (output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--spacing', type=float, default=.1)
    parser.add_argument('--steps', type=int, nargs='+', default=[256, 1024, 4096])
    parser.add_argument('--amplitude', type=float, default=1e-6)
    parser.add_argument('--samples', type=int, default=4)
    parser.add_argument('--restart-check', action='store_true')
    parser.add_argument('--content-model', type=int, choices=(0,1), default=0)
    parser.add_argument('--scaling', type=float, default=1e5)
    parser.add_argument('--gvs', action='store_true')
    parser.add_argument('--cycles', type=int, default=1)
    parser.add_argument('--moving-mesh', action='store_true')
    args = parser.parse_args()
    if (not math.isfinite(args.spacing) or args.spacing <= 0 or
            not math.isfinite(args.amplitude) or args.amplitude <= 0 or
            args.samples < 4 or args.samples % 4 or args.cycles < 1 or
            not math.isfinite(args.scaling) or args.scaling <= 0 or (args.gvs and not args.content_model) or
            (args.moving_mesh and not args.content_model) or
            any(n < args.samples or n % args.samples for n in args.steps)):
        parser.error('Use positive finite spacing/amplitude and step counts divisible by the requested samples')
    study(args.executable.resolve(), args.output.resolve(), args.spacing, args.steps, args.amplitude,
          args.samples, args.restart_check, args.content_model, args.scaling, args.gvs, args.cycles,args.moving_mesh)
