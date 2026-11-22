#!/usr/bin/env python3
"""Run one small model per legacy boundary-condition code and compare two executables.

Every velocity BC code on every side, in 2D and 3D, plus the stress BCs, the period
tables, the 2D depth profile and taper, the z1 loading period, the Winkler and water
loading overrides, the oblique n boundaries with their edge rules, and the PT path.
Each case runs with a parent (reference) executable and a current one: every field of
every save and checkpoint must be bit-identical, and a value the parent rejects must
be rejected with the same exit code. This replaces a make set / make cmp pair per case.

Side-wall codes that set no velocity component are tagged "noop". The odd ones still
enter find_max_vbc, and so dt.

    python3 tests/functional/bc/run_bc_cases.py \
        --parent2d <parent build>/dynearthsol2d --current2d ./dynearthsol2d \
        --parent3d <parent build>/dynearthsol3d --current3d ./dynearthsol3d \
        --workdir <workdir>

Build both with hdf5=1. --cases takes glob patterns of case names; --list prints them.
"""
import argparse
import fnmatch
import json
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT/'tests/functional')]
from des_test_utils import read_config, run_des
from pt_initial_checkpoint import fields

# Distinct values per side and component, so a value moved to another side or sign shows.
VALUES = dict(vbc_val_x0='-1e-9', vbc_val_x1='2e-9', vbc_val_y0='3e-10', vbc_val_y1='-4e-10',
              vbc_val_z0='5e-10', vbc_val_z1='-6e-10', vbc_val_x0_l='7e-10', vbc_val_x1_l='-8e-10',
              vbc_val_y0_l='9e-10', vbc_val_y1_l='-1.1e-9', vbc_val_n0='1.2e-9')
SIDE_CODES = {2: list(range(8)) + [8, 11, 13], 3: list(range(8)) + [8, 11, 13]}
NOOP = {2: {5, 7, 8, 11, 13}, 3: {8, 11, 13}}  # side-wall codes that set no component
Z_CODES = {2: range(5), 3: range(4)}
PT_CONTROL = dict(has_PT='yes', PT_max_iter='20000', PT_relative_tolerance='1e-8',
                  PT_absolute_tolerance='1e-12', is_quasi_static='yes', damping_option='1',
                  damping_factor='0.3', inertial_scaling='1', fixed_dt='0',
                  use_global_velocity_scaling='yes', has_moving_mesh='no')


def base(ndims, mesh):
    cfg = read_config(ROOT/'tests/functional/displacement_output/displacement_output.cfg')
    cfg['sim'].update(modelname='bc', max_steps='6', output_step_interval='2')
    cfg['mesh'].update(xlength='1e4', ylength='1e4', zlength='1e4', resolution='2.5e3')
    if mesh == 'cut':
        cfg['mesh'].update(meshing_option='90', meshing_elem_shape='0',
                           poly_filename=str(HERE/f'cut-corner-{ndims}d.poly'))
    # DES's usual dynamics, so that gravity and every constrained component move nodes.
    cfg['control'].update(gravity='10', fixed_dt='0', inertial_scaling='1e5',
                          damping_option='1', damping_factor='0.8')
    cfg['mat'].update(rho0='[3000]', bulk_modulus='[5e10]', shear_modulus='[3e10]')
    bc = cfg['bc']
    for key in list(bc):
        if key.startswith('vbc'):
            del bc[key]
    bc.update(VALUES)
    bc.update(vbc_x0='1', vbc_x1='1', vbc_z0='1', vbc_z1='0', vbc_val_z0='0')
    if ndims == 3:
        bc.update(vbc_y0='1', vbc_y1='1')
    return cfg


def case(name, ndims, bc=None, control=None, mesh='box', noop=False, exit_code=0):
    return dict(name=name, ndims=ndims, bc=bc or {}, control=control or {}, mesh=mesh,
                noop=noop, exit_code=exit_code)


def cases():
    found = []
    for nd in (2, 3):
        d = f'{nd}d'
        walls = [('x', ('vbc_x0', 'vbc_x1'))] + ([('y', ('vbc_y0', 'vbc_y1'))] if nd == 3 else [])
        for axis, keys in walls:
            for code in SIDE_CODES[nd]:
                found.append(case(f'{d}-{axis}-c{code}', nd, {k: str(code) for k in keys},
                                  noop=code in NOOP[nd]))
        for code in Z_CODES[nd]:
            found.append(case(f'{d}-z-c{code}', nd, dict(vbc_z0=str(code), vbc_z1=str(code))))
        rejected = 5 if nd == 2 else 4
        found.append(case(f'{d}-z-c{rejected}', nd, dict(vbc_z1=str(rejected)), exit_code=11))
        found.append(case(f'{d}-n-c2', nd, dict(vbc_n0='2'), mesh='cut', exit_code=11))
        period = dict(num_vbc_period_x0='2', vbc_period_x0_time_in_yr='[0, 0.3]',
                      vbc_period_x0_ratio='[1, -0.5]', num_vbc_period_x1='3',
                      vbc_period_x1_time_in_yr='[0, 0.2, 0.5]', vbc_period_x1_ratio='[1, 2, 0]')
        if nd == 3:
            period.update(num_vbc_period_y0='2', vbc_period_y0_time_in_yr='[0, 0.3]',
                          vbc_period_y0_ratio='[1, 3]', vbc_y1='6', vbc_x0='7')
        found.append(case(f'{d}-period', nd, period))
        found.append(case(f'{d}-z1-loading-period', nd, dict(vbc_z1='1', vbc_val_z1_loading_period='1.3e7')))
        found.append(case(f'{d}-winkler', nd, dict(has_winkler_foundation='yes', vbc_z0='3')))
        found.append(case(f'{d}-water-loading', nd, dict(has_water_loading='yes', vbc_z1='1')))
        stress = dict(vbc_x0='0', vbc_x1='0', stress_bc_x0='1', stress_val_x0='2e4',
                      stress_bc_x1='3', stress_val_x1='-3e4', stress_bc_z1='3', stress_val_z1='-4e4',
                      stress_bc_z0='1', stress_val_z0='5e4')
        if nd == 3:
            stress.update(vbc_y0='0', vbc_y1='0', stress_bc_y0='2', stress_val_y0='6e4',
                          stress_bc_y1='1', stress_val_y1='-7e4')
        found.append(case(f'{d}-stress', nd, stress))
        for code in (1, 3, 11, 13):
            for sides in ((0, 1) if code in (1, 11) else (1,)):
                others = dict(vbc_x1=str(sides), vbc_z1=str(sides))
                if nd == 3:
                    others.update(vbc_y0=str(sides), vbc_y1=str(sides))
                found.append(case(f'{d}-n-c{code}-sides{sides}', nd, dict(vbc_n0=str(code), **others),
                                  mesh='cut'))
        for code in (1, 3, 6):
            found.append(case(f'{d}-pt-x-c{code}', nd, dict(vbc_x0=str(code), vbc_x1=str(code)),
                              control=PT_CONTROL))
        for code in (1, 11):
            found.append(case(f'{d}-pt-n-c{code}', nd, dict(vbc_n0=str(code)), control=PT_CONTROL,
                              mesh='cut'))
    found.append(case('2d-x0-taper', 2, dict(vbc_x0='3', bottom_shear_zone_thickness='0.4')))
    found.append(case('2d-x-depth-profile', 2, dict(
        vbc_val_division_x0_min='0.25', vbc_val_division_x0_max='0.75', vbc_val_x0_ratio0='1',
        vbc_val_x0_ratio1='0.6', vbc_val_x0_ratio2='0.3', vbc_val_x0_ratio3='0.1',
        vbc_val_division_x1_min='0.4', vbc_val_division_x1_max='0.9', vbc_val_x1_ratio0='0.2',
        vbc_val_x1_ratio1='1', vbc_val_x1_ratio2='0.7', vbc_val_x1_ratio3='0.5')))
    return found


def outputs(folder):
    names = {p.name[:-len('.vtkhdf')] if p.suffix == '.vtkhdf' else p.name
             for p in folder.glob('bc.*') if '.save.' in p.name or '.chkpt.' in p.name}
    return sorted(names)


def run(exe, folder, spec):
    cfg = base(spec['ndims'], spec['mesh'])
    cfg['bc'].update(spec['bc'])
    cfg['control'].update(spec['control'])
    return run_des(exe, folder, 'bc', cfg, timeout=600, check=False)


def compare(parent, current, spec):
    if parent.returncode != spec['exit_code']:
        raise RuntimeError('parent exited {} where {} was expected'.format(parent.returncode, spec['exit_code']))
    if current.returncode != parent.returncode:
        raise RuntimeError('exit {} differs from the parent\'s {}'.format(current.returncode, parent.returncode))
    if parent.returncode:
        return 'rejected as before (exit {})'.format(parent.returncode)
    old, new = parent.folder, current.folder
    names = outputs(old)
    if not names or names != outputs(new):
        raise RuntimeError('output files differ: {} vs {}'.format(names, outputs(new)))
    for name in names:
        a, b = fields(old/name), fields(new/name)
        if a.keys() != b.keys():
            raise RuntimeError('{}: field sets differ: {}'.format(name, sorted(a.keys() ^ b.keys())))
        for key in a:
            if key.split('/')[-1] not in ('walltime_sec', 'provenance') and a[key] != b[key]:
                raise RuntimeError('{}: {} differs'.format(name, key))
    return 'bit-exact ({} files)'.format(len(names))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for nd in (2, 3):
        p.add_argument(f'--parent{nd}d', type=Path)
        p.add_argument(f'--current{nd}d', type=Path)
    p.add_argument('--workdir', type=Path)
    p.add_argument('--cases', nargs='+', default=['*'], help='glob patterns of case names')
    p.add_argument('--list', action='store_true', help='print the case names and exit')
    a = p.parse_args()
    selected = [c for c in cases() if any(fnmatch.fnmatch(c['name'], g) for g in a.cases)]
    if a.list:
        for c in selected:
            print(c['name'] + ('  (noop)' if c['noop'] else '') +
                  ('  (rejected, exit {})'.format(c['exit_code']) if c['exit_code'] else ''))
        return 0
    if not a.workdir:
        p.error('--workdir is required')
    exes = {}
    for nd in (2, 3):
        pair = getattr(a, f'parent{nd}d'), getattr(a, f'current{nd}d')
        if bool(pair[0]) != bool(pair[1]):
            p.error(f'give both --parent{nd}d and --current{nd}d')
        if pair[0]:
            exes[nd] = pair
    selected = [c for c in selected if c['ndims'] in exes]
    if not selected:
        p.error('no case selected for the executables given')
    a.workdir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='bc-', dir=str(a.workdir))).resolve()
    results, failures = [], []
    for spec in selected:
        parent_exe, current_exe = exes[spec['ndims']]
        parent = run(parent_exe, work/'parent'/spec['name'], spec)
        parent.folder = work/'parent'/spec['name']
        current = run(current_exe, work/'current'/spec['name'], spec)
        current.folder = work/'current'/spec['name']
        try:
            verdict = compare(parent, current, spec)
            ok = True
        except RuntimeError as error:
            verdict, ok = str(error), False
            failures.append(spec['name'])
        results.append(dict(case=spec['name'], passed=ok, result=verdict, noop=spec['noop']))
        print('{} {}: {}'.format('PASS' if ok else 'FAIL', spec['name'], verdict), flush=True)
    (work/'results.json').write_text(json.dumps(dict(passed=not failures, cases=results), indent=2) + '\n')
    print('{} of {} cases passed; results in {}'.format(len(results) - len(failures), len(results),
                                                         work/'results.json'))
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
