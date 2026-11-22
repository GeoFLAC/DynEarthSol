#!/usr/bin/env python3
"""Run one small model per legacy boundary-condition code and compare two executables.

Every velocity BC code on every side, in 2D and 3D, plus the stress BCs, the period
tables, the 2D depth profile and taper, the z1 loading period, the Winkler and water
loading overrides, the oblique n boundaries with their edge rules, and the PT path.
Each case runs with a parent (reference) executable and a current one: every field of
every save and checkpoint must be bit-identical, and a value the parent rejects must
be rejected with the same exit code. This replaces a make set / make cmp pair per case.
When the current executable prints its per-side table ([bc][side] lines), the table
must match the conversion tables of #115 in the side frames chosen there.

Side-wall codes that set no velocity component ran on the parent; the current executable
refuses them (#115), and these cases expect its exit code 11.

Cases ending in "-letters" give the same conditions as letters in the side frames
(vbc_x0 = vff, vbc_val_x0 = [..]); they run with the current executable only, and must
match the integer-code case they repeat bit for bit. Other letters cases check
conditions the codes cannot give and the input the letters reject.

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
              vbc_val_y0_l='9e-10', vbc_val_y1_l='-1.1e-9', vbc_val_n0='1.2e-9',
              vbc_val_n1='1.3e-9', vbc_val_n2='1.4e-9', vbc_val_n3='1.5e-9')
SIDE_CODES = {2: list(range(8)) + [8, 11, 13], 3: list(range(8)) + [8, 11, 13]}
NOOP = {2: {5, 7, 8, 11, 13}, 3: {8, 11, 13}}  # side-wall codes that set no component, now refused
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
    # Every code and value is set here, so the expected table needs no DES default.
    bc.update(vbc_x0='1', vbc_x1='1', vbc_y0='1' if ndims == 3 else '0', vbc_y1='1' if ndims == 3 else '0',
              vbc_z0='1', vbc_z1='0', vbc_val_z0='0', vbc_n0='1', vbc_n1='1', vbc_n2='1', vbc_n3='1')
    for side in ('x0', 'x1', 'y0', 'y1', 'z0', 'z1'):
        bc.update({'stress_bc_' + side: '0', 'stress_val_' + side: '0'})
    return cfg


def case(name, ndims, bc=None, control=None, mesh='box', noop=False, exit_code=0, letters=False,
         same_as=None):
    # exit_code is the parent's; a code that sets nothing ran there and is refused now.
    return dict(name=name, ndims=ndims, bc=bc or {}, control=control or {}, mesh=mesh, noop=noop,
                exit_code=exit_code, current_exit=11 if noop else exit_code, letters=letters,
                same_as=same_as)


def effective_cfg(spec):
    cfg = base(spec['ndims'], spec['mesh'])
    cfg['bc'].update(spec['bc'])
    cfg['control'].update(spec['control'])
    return cfg


def as_letters(spec):
    """The same conditions as letters, or None where letters cannot say them."""
    nd = spec['ndims']
    cfg = effective_cfg(spec)
    table = expected_table(nd, cfg)
    winkler, water = support_options(cfg)
    bc = {}
    for side, e in table.items():
        if e['stress'] != 'none' or e['frame'] == 'horizontal-normal':
            return None
        if e['support'] != letters_support(side, e['types'], winkler, water):
            return None  # 2D code 4: supported, although its normal is held
        bc['vbc_' + side] = e['types']
        bc['vbc_val_' + side] = '[' + ', '.join(repr(v) for v in e['values']) + ']'
        if side[0] in 'xy':
            bc['vbc_val_{}_l'.format(side)] = '0'
    taper = float(cfg['bc'].get('bottom_shear_zone_thickness', '0')) > 0
    if nd == 2 and taper and table['x0']['types'][0] == 'v' and not table['x0']['taper']:
        return None  # letters taper any held x0 normal; the codes, 3 only
    letters = case(spec['name'] + '-letters', nd, dict(spec['bc'], **bc), spec['control'], spec['mesh'],
                   letters=True, same_as=spec['name'])
    if max_speed(table) != max_speed(expected_table(nd, effective_cfg(letters))):
        return None  # find_max_vbc skips code 6 and its L; letters count every held value
    return letters


def max_speed(table):
    return max([1e-12] + [e['max_speed'] for e in table.values()])


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
    found.append(case('2d-x0-taper', 2, dict(vbc_x0='3', bottom_shear_zone_thickness='4e3')))
    found.append(case('2d-x-depth-profile', 2, dict(
        vbc_val_division_x0_min='0.25', vbc_val_division_x0_max='0.75', vbc_val_x0_ratio0='1',
        vbc_val_x0_ratio1='0.6', vbc_val_x0_ratio2='0.3', vbc_val_x0_ratio3='0.1',
        vbc_val_division_x1_min='0.4', vbc_val_division_x1_max='0.9', vbc_val_x1_ratio0='0.2',
        vbc_val_x1_ratio1='1', vbc_val_x1_ratio2='0.7', vbc_val_x1_ratio3='0.5')))

    found += [c for c in map(as_letters, [c for c in found if not c['current_exit']]) if c]
    for nd in (2, 3):
        d, zeros = f'{nd}d', '[' + ', '.join(['0'] * nd) + ']'
        held, free = 'v' * nd, 'f' * nd
        # The Winkler foundation and water loading free a held z0 and z1, as they do code 1 to 3.
        found.append(case(f'{d}-winkler-held-letters', nd, dict(has_winkler_foundation='yes',
                          vbc_z0=held, vbc_val_z0=zeros), letters=True, same_as=f'{d}-winkler'))
        found.append(case(f'{d}-water-loading-held-letters', nd, dict(has_water_loading='yes',
                          vbc_z1=held, vbc_val_z1=zeros), letters=True, same_as=f'{d}-water-loading'))
        # What the codes cannot give: a free (supported) n boundary.
        found.append(case(f'{d}-n-free-letters', nd, dict(vbc_n0=free, vbc_val_n0=zeros), mesh='cut',
                          letters=True))
        rejected = [('t', dict(vbc_x0='t' + 'f' * (nd - 1), vbc_val_x0=zeros)),
                    ('length', dict(vbc_x0='v' * (nd + 1), vbc_val_x0='[' + ', '.join(['0'] * (nd + 1)) + ']')),
                    ('letter', dict(vbc_x0='x' + 'f' * (nd - 1), vbc_val_x0=zeros)),
                    ('scalar', dict(vbc_x0='v' + 'f' * (nd - 1), vbc_val_x0='1e-9')),
                    ('code-list', dict(vbc_x0='1', vbc_val_x0=zeros)),
                    ('lateral', dict(vbc_x0=held, vbc_val_x0=zeros, vbc_val_x0_l='1e-10')),
                    ('stress', dict(vbc_x0=free, vbc_val_x0=zeros, stress_bc_x0='1', stress_val_x0='1e4')),
                    ('n-tangent', dict(vbc_n0=held, vbc_val_n0='[1e-9, 1e-10' + ', 0' * (nd - 2) + ']')),
                    ('n-partial', dict(vbc_n0='f' + 'v' * (nd - 1), vbc_val_n0=zeros))]
        if nd == 2:
            rejected.append(('y', dict(vbc_y0='ff', vbc_val_y0='[0, 0]')))
        for what, bc in rejected:
            found.append(case(f'{d}-letters-reject-{what}', nd, bc, mesh='cut' if what.startswith('n-') else 'box',
                              exit_code=11, letters=True))
    # What the codes cannot give: a side wall with a vertical velocity, and a 2D x0 held by
    # its normal only that tapers off over the bottom shear zone.
    found.append(case('3d-x0-vertical-letters', 3, dict(vbc_x0='vvv', vbc_val_x0='[1e-9, 2e-10, -3e-10]',
                                                        vbc_val_x0_l='0'), letters=True))
    found.append(case('2d-x0-taper-vf-letters', 2, dict(vbc_x0='vf', vbc_val_x0='[1e-9, 0]', vbc_val_x0_l='0',
                                                        bottom_shear_zone_thickness='4e3'), letters=True))
    return found


# What each code holds today, per global component (the conversion tables of #115):
# N is the normal value before the period table, V the raw value, L the lateral one;
# None leaves a component free. Side walls list (normal, lateral, z) in 3D, (x, z) in 2D.
WALL_ROWS = {3: {0: (None, None, None), 1: ('N', None, None), 2: (None, 0, 0), 3: ('N', 0, 0),
                 4: (None, 'V', 0), 5: (0, 'V', 0), 6: ('N', 'L', None), 7: ('N', 0, None)},
             2: {0: (None, None), 1: ('N', None), 2: (None, 0), 3: ('N', 0), 4: (0, 'V'), 6: ('N', 'L')}}
Z_ROWS = {3: {0: (None, None, None), 1: (None, None, 'V'), 2: (0, 0, None), 3: (0, 0, 'V')},
          2: {0: (None, None), 1: (None, 'V'), 2: (0, None), 3: (0, 'V'), 4: ('V', 0)}}
# The side frames (n, t_h, t_v) chosen in #115, as signed global axes; 2D keeps the in-plane pair.
FRAMES = {3: dict(x0='-x,-y,+z', x1='+x,+y,+z', y0='-y,+x,+z', y1='+y,-x,+z', z0='-z,+x,-y', z1='+z,+x,+y'),
          2: dict(x0='-x,+z', x1='+x,+z', z0='-z,+x', z1='+z,+x')}


def support_options(cfg):
    gravity = float(cfg['control'].get('gravity', '10')) != 0
    winkler = gravity and cfg['bc'].get('has_winkler_foundation', 'yes') in ('yes', 'true', '1')
    water = gravity and cfg['bc'].get('has_water_loading', 'yes') in ('yes', 'true', '1')
    return winkler, water


def letters_support(side, types, winkler, water):
    free = types[0] == 'f'
    return free and (winkler if side == 'z0' else water if side == 'z1' else True)


def letters_entry(ndims, side, bc, winkler, water):
    """A side given as letters: the values are already in its frame."""
    types = bc['vbc_' + side]
    values = [float(x) for x in bc['vbc_val_' + side].strip('[]').split(',')]
    if (side == 'z0' and winkler) or (side == 'z1' and water):
        types = 'f' * ndims
    held = [v for c, v in zip(types, values) if c == 'v']
    return dict(frame=FRAMES[ndims][side] if side[0] != 'n' else 'normal', types=types,
                values=[v if c == 'v' else 0.0 for c, v in zip(types, values)],
                normal_follows_period=side[0] in 'xy' and types[0] == 'v',
                taper=ndims == 2 and side == 'x0' and types[0] == 'v',
                support=letters_support(side, types, winkler, water),
                max_speed=max([abs(v) for v in held] or [0.0]),
                release_time=bc.get('vbc_val_z1_loading_period') if side == 'z1' else None, stress='none')


def expected_table(ndims, cfg):
    bc = cfg['bc']
    winkler, water = support_options(cfg)
    axes = 'xyz' if ndims == 3 else 'xz'
    table = {}
    for side in ('x0', 'x1', 'y0', 'y1', 'z0', 'z1') if ndims == 3 else ('x0', 'x1', 'z0', 'z1'):
        if bc['vbc_' + side][0].isalpha():
            table[side] = letters_entry(ndims, side, bc, winkler, water)
            continue
        code = int(bc['vbc_' + side])
        if (side == 'z0' and winkler) or (side == 'z1' and water):
            code = 0
        known = dict(N=float(bc['vbc_val_' + side]), V=float(bc['vbc_val_' + side]),
                     L=float(bc.get('vbc_val_{}_l'.format(side), '0')))
        if side[0] == 'z':
            rows, components = Z_ROWS[ndims], axes
        else:
            rows = WALL_ROWS[ndims]
            components = (side[0], 'y' if side[0] == 'x' else 'x', 'z') if ndims == 3 else 'xz'
        row = rows.get(code)
        held = {} if row is None else {a: (known[r] if isinstance(r, str) else r)
                                       for a, r in zip(components, row) if r is not None}
        types, values = '', []
        for axis in FRAMES[ndims][side].split(','):
            sign = -1.0 if axis[0] == '-' else 1.0
            types += 'v' if axis[1] in held else 'f'
            values.append(sign * held[axis[1]] if axis[1] in held else 0.0)
        stress = int(bc['stress_bc_' + side])
        stress_axes = {1: 'x', 2: 'y', 3: 'z'} if ndims == 3 else {1: 'x', 3: 'z'}
        period = row is not None and side[0] != 'z' and row[0] == 'N'
        table[side] = dict(
            frame=FRAMES[ndims][side], types=types, values=values, normal_follows_period=period,
            taper=ndims == 2 and side == 'x0' and code == 3,
            support=code in (0, 2, 4) and (side != 'z0' or winkler) and (side != 'z1' or water),
            max_speed=abs(known['V']) if code % 2 == 1 or code == 4 else 0.0,
            release_time=bc.get('vbc_val_z1_loading_period') if side == 'z1' else None,
            stress='{}:{}'.format(stress_axes[stress], float(bc['stress_val_' + side]))
                   if stress in stress_axes else 'none')
    for side in ('n0', 'n1', 'n2', 'n3'):
        if bc['vbc_' + side][0].isalpha():
            table[side] = letters_entry(ndims, side, bc, winkler, water)
            continue
        code, value = int(bc['vbc_' + side]), float(bc['vbc_val_' + side])
        full = code in (3, 13)
        table[side] = dict(frame='horizontal-normal' if code in (11, 13) else 'normal',
                           types='v' + ('v' if full else 'f') * (ndims - 1),
                           values=[value] + [0.0] * (ndims - 1), normal_follows_period=False, taper=False,
                           support=False, max_speed=abs(value), release_time=None, stress='none')
    return table


def check_table(log, ndims, cfg):
    printed = {}
    for line in log.splitlines():
        if line.startswith('[bc][side] '):
            row = dict(item.split('=', 1) for item in line.split()[1:])
            printed[row.pop('name')] = row
    if not printed:
        return False
    expected = expected_table(ndims, cfg)
    if printed.keys() != expected.keys():
        raise RuntimeError('table sides {} differ from {}'.format(sorted(printed), sorted(expected)))
    for side, want in expected.items():
        got = printed[side]
        release = got['release_time']
        have = dict(frame=got['frame'], types=got['types'],
                    values=[float(x) for x in got['values'].split(',')],
                    normal_follows_period=got['normal_follows_period'] == '1', taper=got['taper'] == '1',
                    support=got['support'] == '1',
                    max_speed=float(got['max_speed']),
                    release_time=None if release == 'none' else float(release),
                    stress=got['stress'] if got['stress'] == 'none' else
                    '{}:{}'.format(got['stress'][0], float(got['stress'][2:])))
        if want['release_time'] is not None:
            want = dict(want, release_time=float(want['release_time']))
        if have != want:
            diff = {k: (have[k], want[k]) for k in want if have[k] != want[k]}
            raise RuntimeError('table side {}: printed vs expected {}'.format(side, diff))
    return True


def outputs(folder):
    names = {p.name[:-len('.vtkhdf')] if p.suffix == '.vtkhdf' else p.name
             for p in folder.glob('bc.*') if '.save.' in p.name or '.chkpt.' in p.name}
    return sorted(names)


def run(exe, folder, spec):
    cfg = effective_cfg(spec)
    result = run_des(exe, folder, 'bc', cfg, timeout=600, check=False)
    result.folder, result.cfg = folder, cfg
    return result


def compare(parent, current, spec):
    if parent.returncode != spec['exit_code']:
        raise RuntimeError('parent exited {} where {} was expected'.format(parent.returncode, spec['exit_code']))
    if current.returncode != spec['current_exit']:
        raise RuntimeError('exit {} where {} was expected (the parent\'s: {})'.format(
            current.returncode, spec['current_exit'], parent.returncode))
    if current.returncode:
        return '{} (exit {})'.format('rejected as before' if parent.returncode else 'now rejected',
                                     current.returncode)
    # A build that prints the per-side table must print the one #115 expects.
    table = check_table(current.stdout, spec['ndims'], current.cfg)
    count = same_outputs(parent.folder, current.folder)
    return 'bit-exact ({} files){}'.format(count, ', table as expected' if table else '')


def same_outputs(old, new):
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
    return len(names)


def compare_letters(current, spec, reference):
    if current.returncode != spec['exit_code']:
        raise RuntimeError('exited {} where {} was expected'.format(current.returncode, spec['exit_code']))
    if current.returncode:
        return 'rejected (exit {})'.format(current.returncode)
    if not check_table(current.stdout, spec['ndims'], current.cfg):
        raise RuntimeError('no per-side table printed')
    if reference is None:
        return 'runs, table as expected'
    count = same_outputs(reference.folder, current.folder)
    return 'bit-exact with {} ({} files), table as expected'.format(spec['same_as'], count)


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
            print(c['name'] + ('  (sets nothing; now rejected)' if c['noop'] else '') +
                  ('  (rejected, exit {})'.format(c['exit_code']) if c['exit_code'] else '') +
                  ('  (same as {})'.format(c['same_as']) if c['same_as'] else ''))
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
    results, failures, ran = [], [], {}
    every = {c['name']: c for c in cases()}
    for spec in selected:
        parent_exe, current_exe = exes[spec['ndims']]
        if not spec['letters']:
            parent = run(parent_exe, work/'parent'/spec['name'], spec)
        current = run(current_exe, work/'current'/spec['name'], spec)
        ran[spec['name']] = current
        reference = None
        if spec['same_as']:
            if spec['same_as'] not in ran:
                ran[spec['same_as']] = run(current_exe, work/'current'/spec['same_as'], every[spec['same_as']])
            reference = ran[spec['same_as']]
        try:
            verdict = compare_letters(current, spec, reference) if spec['letters'] else compare(parent, current, spec)
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
