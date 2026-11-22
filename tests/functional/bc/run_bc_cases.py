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

Input that adds nothing ran on the parent; the current executable refuses it (#115), and
these cases expect its exit code 11: side-wall codes that set no velocity component, and
stress_bc_* of an unknown type or on a 2D y wall.

Cases marked restart also restart from checkpoint 1 with the current executable, and
every later frame must match the run's bit for bit: the n boundaries keep the normals of
the starting mesh, which the checkpoints now carry, and the restart keeps the velocity.
Cases with changed_restart restart both executables from checkpoint 1 with those bc
values changed, and every later frame must match the parent's restart bit for bit.
Cases with edge_restart restart from checkpoint 1 with x1 holding, at an x1-n0 edge node,
the y velocity the run's conditions give the checkpointed one there: the restart must match
one from the checkpoint without its record of the conditions. Restarting with only n0's
traction changed, or n1 freed, which has no node, must keep the checkpointed velocity.
Cases marked corner check, at the nodes on n0 and the given sides only, that the velocity
along the normal of each n boundary there is the value it holds and that the components the
sides hold keep their values.

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
import copy
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
    if mesh != 'box':
        poly = dict(cut='cut-corner', cuts='two-cuts', slope='slope-top')[mesh]
        cfg['mesh'].update(meshing_option='90', meshing_elem_shape='0',
                           poly_filename=str(HERE/f'{poly}-{ndims}d.poly'))
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
         same_as=None, refused=False, restart=False, corner=None, quoted=(), changed_restart=None,
         edge_restart=False):
    # exit_code is the parent's; input that set nothing ran there and is refused now.
    # corner is (sides, frame, {global axis: value held there}); quoted lists the bc keys
    # written in quotes; changed_restart holds the bc values a restart changes.
    return dict(name=name, ndims=ndims, bc=bc or {}, control=control or {}, mesh=mesh, noop=noop,
                exit_code=exit_code, current_exit=11 if noop or refused else exit_code, letters=letters,
                same_as=same_as, restart=restart, corner=corner, quoted=quoted,
                changed_restart=changed_restart, edge_restart=edge_restart)


def effective_cfg(spec, quoted=False):
    """The case's configuration; with quoted, as written for DES (spec['quoted'] in quotes)."""
    cfg = base(spec['ndims'], spec['mesh'])
    cfg['bc'].update(spec['bc'])
    cfg['control'].update(spec['control'])
    if spec['restart'] or spec['changed_restart'] or spec['edge_restart']:
        cfg['sim']['has_marker_output'] = 'yes'  # a restart reads the markers
    for k, key in enumerate(spec['quoted'] if quoted else ()):
        mark = '"' if k % 2 == 0 else "'"
        cfg['bc'][key] = mark + cfg['bc'][key] + mark
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
        bc['vbc_' + side] = e['types']
        if side[0] != 'n':
            bc['stress_bc_' + side] = '0'
        bc['vbc_val_' + side] = '[' + ', '.join(repr(v) for v in e['values']) + ']'
        if side[0] in 'xy':
            bc['vbc_val_{}_l'.format(side)] = '0'
    taper = float(cfg['bc'].get('bottom_shear_zone_thickness', '0')) > 0
    if nd == 2 and taper and table['x0']['types'][0] == 'v' and not table['x0']['taper']:
        return None  # letters taper any held x0 normal; the codes, 3 only
    return case(spec['name'] + '-letters', nd, dict(spec['bc'], **bc), spec['control'], spec['mesh'],
                letters=True, same_as=spec['name'], restart=spec['restart'], corner=spec['corner'])


def max_speed(table):
    return max([1e-12] + [e['max_speed'] for e in table.values()])


def parent_max_speed(ndims, cfg):
    """find_max_vbc before #115: |vbc_val| of the side codes that are odd or 4, and of the
    n codes, whatever the code held."""
    bc = cfg['bc']
    winkler, water = support_options(cfg)
    speed = 1e-12
    for side in ('x0', 'x1', 'y0', 'y1', 'z0', 'z1', 'n0', 'n1', 'n2', 'n3'):
        code = int(bc['vbc_' + side])
        if (side == 'z0' and winkler) or (side == 'z1' and water):
            code = 0
        if code % 2 == 1 or (code == 4 and side[0] != 'n'):
            speed = max(speed, abs(float(bc['vbc_val_' + side])))
    return speed


def parent_support(side, code, winkler, water):
    """Before #115 the codes 0, 2 and 4 carried the support, held normal or not."""
    return code in (0, 2, 4) and (side != 'z0' or winkler) and (side != 'z1' or water)


def changes(spec):
    """What #115 changed for this case, compared with the parent, or None."""
    if spec['letters'] or spec['current_exit']:
        return None
    nd, cfg = spec['ndims'], effective_cfg(spec)
    table = expected_table(nd, cfg)
    winkler, water = support_options(cfg)
    found = []
    old, new = parent_max_speed(nd, cfg), max_speed(table)
    if old != new:
        found.append('find_max_vbc {:g} -> {:g}'.format(old, new))
    for side, e in table.items():
        if side[0] == 'n':
            continue
        code = 0 if (side == 'z0' and winkler) or (side == 'z1' and water) else int(cfg['bc']['vbc_' + side])
        if parent_support(side, code, winkler, water) != e['support']:
            found.append('{} support {}'.format(side, 'added' if e['support'] else 'removed'))
    if spec['mesh'] == 'cut' and table['z1']['release_time'] is not None and table['z1']['types'][0] == 'v':
        found.append('n0 sets its normal at the z1 edge once z1 is released')
    if spec['mesh'] == 'cut' and table['n0']['types'] == 'v' + 'f' * (nd - 1):
        normal = dict(zip('xyz' if nd == 3 else 'xz', CUT_NORMAL[nd]))
        if table['n0']['frame'] == 'horizontal-normal':
            normal['z'] = 0.0
        for side in CUT_SIDES[nd]:
            e = table[side]
            held = {a[1] for a, t in zip(FRAMES[nd][side].split(','), e['types']) if t == 'v'}
            if not held or e['types'] == 'v' + 'f' * (nd - 1):
                continue
            if any(normal[a] for a in held) and any(normal[a] for a in normal if a not in held):
                found.append('n0 sets its normal on the {} edge through what {} leaves free'.format(side, side))
    return '; '.join(found) or None


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
        found.append(case(f'{d}-z1-loading-period', nd, dict(vbc_z1='1', vbc_val_z1_loading_period='5e6')))
        found.append(case(f'{d}-winkler', nd, dict(has_winkler_foundation='yes', vbc_z0='3')))
        found.append(case(f'{d}-water-loading', nd, dict(has_water_loading='yes', vbc_z1='1')))
        stress = dict(vbc_x0='0', vbc_x1='0', stress_bc_x0='1', stress_val_x0='2e4',
                      stress_bc_x1='3', stress_val_x1='-3e4', stress_bc_z1='3', stress_val_z1='-4e4',
                      stress_bc_z0='1', stress_val_z0='5e4')
        if nd == 3:
            stress.update(vbc_y0='0', vbc_y1='0', stress_bc_y0='2', stress_val_y0='6e4',
                          stress_bc_y1='1', stress_val_y1='-7e4')
        found.append(case(f'{d}-stress', nd, stress, restart=True))
        # A tangential component loads a tilted facet: x on a top that rises to the east.
        found.append(case(f'{d}-stress-slope', nd, dict(stress_bc_z1='1', stress_val_z1='-4e4'), mesh='slope'))
        # A held component is loaded once its hold ends.
        found.append(case(f'{d}-stress-z1-released', nd, dict(vbc_z1='1', vbc_val_z1_loading_period='5e6',
                                                              stress_bc_z1='3', stress_val_z1='-4e4')))
        # A component the velocity code holds: nothing on a box, but where the side meets an
        # n boundary the edge map passes the load on (x on z1 at the cut corner).
        found.append(case(f'{d}-stress-held-normal', nd, dict(stress_bc_x0='1', stress_val_x0='2e4')))
        found.append(case(f'{d}-stress-held-tangential', nd, dict(vbc_x0='2', stress_bc_x0='3',
                                                                  stress_val_x0='2e4')))
        found.append(case(f'{d}-stress-held-n-corner', nd, dict(vbc_z1='2', stress_bc_z1='1',
                                                                stress_val_z1='-4e4'), mesh='cut',
                          corner=(('z1',), 3, dict(x=0.0, y=0.0) if nd == 3 else dict(x=0.0))))
        # What adds nothing is refused: an unknown type, and (2D) the y walls.
        found.append(case(f'{d}-stress-type4', nd, dict(vbc_x0='0', stress_bc_x0='4', stress_val_x0='2e4'),
                          refused=True))
        if nd == 2:
            found.append(case('2d-stress-y', 2, dict(stress_bc_y0='2', stress_val_y0='2e4'), refused=True))
        for code in (1, 3, 11, 13):
            for sides in ((0, 1) if code in (1, 11) else (1,)):
                others = dict(vbc_x1=str(sides), vbc_z1=str(sides))
                if nd == 3:
                    others.update(vbc_y0=str(sides), vbc_y1=str(sides))
                found.append(case(f'{d}-n-c{code}-sides{sides}', nd, dict(vbc_n0=str(code), **others),
                                  mesh='cut', restart=True))
        # Once z1 is released, n0 sets its normal on their edge, as beside a free z1.
        found.append(case(f'{d}-n-z1-released', nd, dict(vbc_z1='1', vbc_val_z1_loading_period='5e6'),
                          mesh='cut', corner=(('z1',), 3, {})))
        # A restart that changes a held velocity applies it before its first step, as the
        # parent does.
        found.append(case(f'{d}-restart-z1-changed', nd, dict(vbc_z1='1'),
                          changed_restart=dict(vbc_val_z1='-3e-9')))
        # Quoted values, as the codes and scalars took them before they became text and lists.
        found.append(case(f'{d}-quoted', nd, dict(vbc_x0='3', vbc_val_x0='-1e-9'),
                          quoted=('vbc_val_x0', 'vbc_x0')))
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

    found += [c for c in map(as_letters, [c for c in found if not (c['current_exit'] or c['changed_restart'])])
              if c]
    for nd in (2, 3):
        d, zeros = f'{nd}d', '[' + ', '.join(['0'] * nd) + ']'
        held, free = 'v' * nd, 'f' * nd
        # The Winkler foundation and water loading free a held z0 and z1, as they do code 1 to 3.
        found.append(case(f'{d}-winkler-held-letters', nd, dict(has_winkler_foundation='yes',
                          vbc_z0=held, vbc_val_z0=zeros), letters=True, same_as=f'{d}-winkler'))
        found.append(case(f'{d}-water-loading-held-letters', nd, dict(has_water_loading='yes',
                          vbc_z1=held, vbc_val_z1=zeros), letters=True, same_as=f'{d}-water-loading'))
        # What the codes cannot give: a free (supported) n boundary, and one under a normal
        # traction. Both bend; their frame stays that of the starting mesh, also on restart.
        found.append(case(f'{d}-n-free-letters', nd, dict(vbc_n0=free, vbc_val_n0=zeros), mesh='cut',
                          letters=True, restart=True))
        found.append(case(f'{d}-n-normal-traction-letters', nd, dict(vbc_n0='t' + 'f' * (nd - 1),
                          vbc_val_n0='[-2e4' + ', 0' * (nd - 1) + ']'), mesh='cut', letters=True,
                          restart=True))
        # Tractions the codes cannot give: shear on a side wall and on an n boundary, and a
        # traction that stays when the Winkler foundation or water loading frees the velocity.
        found.append(case(f'{d}-x0-shear-letters', nd, dict(vbc_x0='ft' + 'f' * (nd - 2), vbc_val_x0_l='0',
                          vbc_val_x0='[0, 5e4' + ', 0' * (nd - 2) + ']'), letters=True))
        found.append(case(f'{d}-n-traction-letters', nd, dict(vbc_n0='vt' + 'f' * (nd - 2),
                          vbc_val_n0='[1.2e-9, -2e4' + ', 0' * (nd - 2) + ']'), mesh='cut', letters=True,
                          restart=True, edge_restart=nd == 3))
        found.append(case(f'{d}-winkler-traction-letters', nd, dict(has_winkler_foundation='yes',
                          vbc_z0='tv' + 'f' * (nd - 2), vbc_val_z0='[-3e4, 1e-10' + ', 0' * (nd - 2) + ']'),
                          letters=True))
        found.append(case(f'{d}-water-loading-traction-letters', nd, dict(has_water_loading='yes',
                          vbc_z1='tv' + 'f' * (nd - 2), vbc_val_z1='[-2e4, 1e-10' + ', 0' * (nd - 2) + ']'),
                          letters=True))
        rejected = [('length', dict(vbc_x0='v' * (nd + 1), vbc_val_x0='[' + ', '.join(['0'] * (nd + 1)) + ']')),
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
    # A side wall holding only a tangent that n0's normal does not touch: n0 sets its normal
    # on their edge.
    found.append(case('3d-n-x1-tangent-letters', 3, dict(vbc_x1='fvf', vbc_val_x1='[0, 3e-10, 0]',
                                                         vbc_val_x1_l='0'), mesh='cut', letters=True,
                      corner=(('x1',), 3, dict(y=3e-10))))
    # Two side walls holding a tangent each, x1 its y and y0 its x: on their corner, n0 sets
    # its normal through z, which neither holds.
    found.append(case('3d-n-two-tangents-letters', 3, dict(vbc_x1='fvf', vbc_val_x1='[0, 3e-10, 0]',
                                                          vbc_val_x1_l='0', vbc_y0='fvf',
                                                          vbc_val_y0='[0, 0, 0]', vbc_val_y0_l='0'),
                      mesh='cut', letters=True, corner=(('x1', 'y0'), 3, dict(x=0.0, y=3e-10))))
    # x1 holding only its normal and y0 only z, at values n0's normal can keep: on their corner
    # n0 moves neither.
    found.append(case('3d-n-normal-and-tangent-letters', 3, dict(vbc_x1='vff', vbc_val_x1='[0, 0, 0]',
                                                                vbc_val_x1_l='0', vbc_y0='ffv',
                                                                vbc_val_y0=f'[0, 0, {Z_HOLD!r}]',
                                                                vbc_val_y0_l='0'),
                      mesh='cut', letters=True, corner=(('x1', 'y0'), 3, dict(x=0.0, z=Z_HOLD))))
    # Two n boundaries holding their normals where they meet y0 holding only z: on that
    # corner n1 sets its normal keeping z and n0's normal.
    found.append(case('3d-two-n-and-tangent-letters', 3, dict(vbc_y0='ffv', vbc_val_y0='[0, 0, 0]',
                                                             vbc_val_y0_l='0', vbc_n1='vff',
                                                             vbc_val_n1='[7e-10, 0, 0]'),
                      mesh='cuts', letters=True, corner=(('y0', 'n1'), 3, dict(z=0.0))))
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
    free = types[0] != 'v'  # f, or t, which adds to the support
    return free and (winkler if side == 'z0' else water if side == 'z1' else True)


def side_values(bc, side):
    """The values in vbc_val_<side>: a list such as [1e-9, 0, 0], or one number."""
    return [float(x) for x in bc['vbc_val_' + side].strip('[]').split(',')]


def letters_entry(ndims, side, bc, winkler, water):
    """A side given as letters: the values are already in its frame."""
    types = bc['vbc_' + side]
    values = side_values(bc, side)
    if (side == 'z0' and winkler) or (side == 'z1' and water):
        types = types.replace('v', 'f')  # the velocity is freed; a traction stays
    held = [v for c, v in zip(types, values) if c == 'v']
    return dict(frame=FRAMES[ndims][side] if side[0] != 'n' else 'normal', types=types,
                values=[v if c != 'f' else 0.0 for c, v in zip(types, values)],
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
            support=letters_support(side, types, winkler, water),
            max_speed=max([abs(v) for v in held.values()] or [0.0]),
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
        # The log prints 6 significant digits, as '%g' does.
        printed_as = lambda x: float('{:g}'.format(float(x)))
        want = dict(want, values=[printed_as(v) for v in want['values']],
                    max_speed=printed_as(want['max_speed']))
        if want['release_time'] is not None:
            want = dict(want, release_time=printed_as(want['release_time']))
        if have != want:
            diff = {k: (have[k], want[k]) for k in want if have[k] != want[k]}
            raise RuntimeError('table side {}: printed vs expected {}'.format(side, diff))
    return True


def outputs(folder, model='bc'):
    names = {p.name[:-len('.vtkhdf')] if p.suffix == '.vtkhdf' else p.name
             for p in folder.glob(model + '.*') if '.save.' in p.name or '.chkpt.' in p.name}
    return sorted(names)


def run(exe, folder, spec):
    result = run_des(exe, folder, 'bc', effective_cfg(spec, quoted=True), timeout=600, check=False)
    result.folder, result.cfg = folder, effective_cfg(spec)
    return result


# The outward unit normal of n0 in cut-corner-{2,3}d.poly, the sides it meets there, and the
# boundary flag bits.
CUT_NORMAL = {2: (0.5**0.5, 0.5**0.5), 3: (0.5**0.5, 0.0, 0.5**0.5)}
CUT_SIDES = {2: ('x1', 'z1'), 3: ('x1', 'y0', 'y1', 'z1')}
Z_HOLD = 1.2e-9 * 2**0.5  # with x = 0, the z velocity that gives n0's 1.2e-9 in 3D
FLAG_BITS = dict(x0=1, x1=2, y0=4, y1=8, z0=16, z1=32, n0=64, n1=128)


def saved_field(saved, name):
    """The field name of a save or checkpoint read by fields(), as a flat array."""
    import numpy as np
    key = next(k for k in saved if not k.startswith('@') and k.split('/')[-1] == name)
    return np.frombuffer(saved[key], eval(saved['@' + key + '.layout'].decode())[1])


def corner_normal_velocity(result, spec):
    """At the nodes on n0 and the spec['corner'] sides only, the velocity along the normal of
    each n boundary there is the value it holds, and the components the sides hold keep their
    values."""
    import numpy as np
    sides, frame, held_there = spec['corner']
    side = '+'.join(sides)
    nd = spec['ndims']
    saved = fields(result.folder/'bc.save.{:06d}'.format(frame))
    flag = saved_field(saved, 'bcflag')
    velocity = saved_field(saved, 'velocity')
    velocity = velocity.reshape(flag.size, -1)[:, :nd]  # VTKHDF keeps 3 components in 2D
    on = flag == FLAG_BITS['n0'] + sum(FLAG_BITS[s] for s in sides)
    if not on.any():
        raise RuntimeError('no node lies on n0 and {} only'.format(side))
    chkpt = fields(result.folder/'bc.chkpt.{:06d}'.format(frame))
    normals = saved_field(chkpt, 'boundary normals').reshape(-1, nd)
    checked = []
    for b in ['n0'] + [s for s in sides if s[0] == 'n']:
        held = side_values(result.cfg['bc'], b)[0]
        checked.append('{} {:g}'.format(b, held))
        normal = velocity[on] @ normals[FLAG_BITS[b].bit_length() - 1]
        worst = float(np.max(np.abs(normal - held)))
        if worst > 1e-9*abs(held):
            raise RuntimeError('velocity along the {} normal at the {} nodes on n0 and {}: {} instead '
                               'of {}'.format(b, int(on.sum()), side, normal.tolist(), held))
    for axis, value in held_there.items():
        got = velocity[on][:, 'xyz'.index(axis) if nd == 3 else 'xz'.index(axis)]
        if (got != value).any():
            raise RuntimeError('{} velocity at the {} nodes on n0 and {}: {} instead of {}'.format(
                axis, int(on.sum()), side, got.tolist(), value))
    return 'normal velocity {} at the {} nodes on n0 and {}{}'.format(
        ', '.join(checked), int(on.sum()), side,
        ''.join(', {} {:g}'.format(a, v) for a, v in held_there.items()))


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
    change = changes(spec)
    if change:
        return 'changed as expected ({}): {}'.format(change, largest_difference(parent.folder, current.folder))
    count = same_outputs(parent.folder, current.folder)
    return 'bit-exact ({} files){}'.format(count, ', table as expected' if table else '')


def largest_difference(old, new):
    """The largest relative difference of a floating field over every save and checkpoint."""
    import numpy as np
    names = outputs(old)
    if names != outputs(new):
        raise RuntimeError('output files differ: {} vs {}'.format(names, outputs(new)))
    worst, where = 0.0, None
    for file in names:
        a, b = fields(old/file), fields(new/file)
        for key in a:
            layout = a.get('@' + key + '.layout')
            if key.startswith('@') or layout is None or key.split('/')[-1] == 'walltime_sec':
                continue
            dtype = eval(layout.decode())[1]
            if not dtype.lstrip('<>|=').startswith('f') or key not in b:
                continue
            x, y = np.frombuffer(a[key], dtype), np.frombuffer(b[key], dtype)
            if x.size == 0 or x.size != y.size:
                continue
            diff = float(np.max(np.abs(x - y))) / (float(np.max(np.abs(x))) or 1.0)
            if diff > worst:
                worst, where = diff, '{} of {}'.format(key.split('/')[-1], file)
    if worst == 0:
        raise RuntimeError('the outputs did not change')
    return 'largest relative difference {:.2g} in {}'.format(worst, where)


# What the current build adds to its checkpoints: the boundary normals that a restart keeps,
# and the record of the input of the velocity conditions that it compares.
RECORD = ('boundary condition control count', 'boundary condition controls')
ADDED_FIELDS = ('boundary normals',) + RECORD


def field_name(key):
    """The field a key of fields() holds, for its data and its layout."""
    return key.split('/')[-1].lstrip('@').replace('.layout', '')


def same_fields(a, b, name, ignored=()):
    added = {k for k in b.keys() - a.keys() if field_name(k) in ADDED_FIELDS}
    if a.keys() != b.keys() - added:
        raise RuntimeError('{}: field sets differ: {}'.format(name, sorted(a.keys() ^ (b.keys() - added))))
    for key in a:
        if field_name(key) not in ('walltime_sec', 'provenance') + ignored and a[key] != b[key]:
            raise RuntimeError('{}: {} differs'.format(name, key))


def same_outputs(old, new, model='bc', ignored=()):
    names = outputs(old, model)
    if not names or names != outputs(new, model):
        raise RuntimeError('output files differ: {} vs {}'.format(names, outputs(new, model)))
    for name in names:
        same_fields(fields(old/name), fields(new/name), name, ignored)
    return len(names)


def restart_from_1(exe, folder, result, model, bc=None):
    """Restart the run of result in folder from checkpoint 1 as model, with bc values changed."""
    cfg = copy.deepcopy(result.cfg)
    cfg['bc'].update(bc or {})
    cfg['sim'].update(modelname=model, is_restarting='true', restarting_from_modelname='bc',
                      restarting_from_frame='1')
    again = run_des(exe, folder, model, cfg, timeout=600, check=False)
    if again.returncode:
        raise RuntimeError('the {} restart from checkpoint 1 exited {}'.format(model, again.returncode))


def same_after_restart(exe, result):
    """Restart the run from checkpoint 1: every later frame must match the run's."""
    restart_from_1(exe, result.folder, result, 'again')
    later = [n for n in outputs(result.folder) if int(n.split('.')[-1]) > 1]
    for name in later:
        twin = result.folder/name.replace('bc.', 'again.', 1)
        if not any(twin.parent.glob(twin.name + '*')):
            raise RuntimeError('the restart wrote no {}'.format(twin.name))
        same_fields(fields(result.folder/name), fields(twin), name)
    if not later:
        raise RuntimeError('no frame after checkpoint 1')
    return 'restart from checkpoint 1 matches {} files'.format(len(later))


def same_changed_restart(parent_exe, current_exe, parent, current, spec):
    """Restart both runs from checkpoint 1 with spec['changed_restart']: every later frame
    must match the parent's restart."""
    for exe, result in ((parent_exe, parent), (current_exe, current)):
        restart_from_1(exe, result.folder, result, 'changed', spec['changed_restart'])
    count = same_outputs(parent.folder, current.folder, 'changed')
    return 'changed restart matches the parent\'s ({} files)'.format(count)


def edge_restart(exe, result):
    """Restart from checkpoint 1 with x1 holding, at the x1-n0 edge node where the edge map
    changes it most, the y velocity the run's conditions give the checkpointed one, and n0
    free but for its traction. The conditions differ, so the restart applies them, as one from
    the checkpoint without its record does; their result at that node does not tell. With
    only n0's traction changed, or n1 freed, which has no node, they are the run's, and the
    restart keeps the velocity."""
    import h5py
    import numpy as np
    import shutil
    folder, bare = result.folder, result.folder/'without-record'
    bare.mkdir()
    for name in ('bc.info', 'bc.save.000001.vtkhdf', 'bc.chkpt.000001.vtkhdf'):
        shutil.copyfile(folder/name, bare/name)
    with h5py.File(bare/'bc.chkpt.000001.vtkhdf', 'r+') as stream:
        found = []
        stream.visit(lambda key: found.append(key) if key.split('/')[-1] in RECORD else None)
        for key in found:
            del stream[key]
    restart_from_1(exe, bare, result, 'old')  # the run's conditions applied to the checkpointed velocity
    run = fields(folder/'bc.save.000001')
    flag = saved_field(run, 'bcflag')
    saved = saved_field(run, 'velocity').reshape(flag.size, -1)
    mapped = saved_field(fields(bare/'old.save.000001'), 'velocity').reshape(flag.size, -1)
    edge = np.flatnonzero(flag == FLAG_BITS['n0'] + FLAG_BITS['x1'])
    node = edge[np.argmax(np.abs(mapped[edge, 1] - saved[edge, 1]))]
    if mapped[node, 1] == saved[node, 1]:
        raise RuntimeError('the edge map leaves the y velocity of every node on x1 and n0')
    if (flag & FLAG_BITS['n1']).any():
        raise RuntimeError('a node lies on n1')
    for model, what, bc in (('traction', 'the n0 traction', dict(vbc_val_n0='[1.2e-9, -4e4, 0]')),
                            ('unused', 'n1, which has no node,', dict(vbc_n1='fff', vbc_val_n1='[0, 0, 0]'))):
        restart_from_1(exe, folder, result, model, bc)
        kept = saved_field(fields(folder/(model + '.save.000001')), 'velocity').reshape(flag.size, -1)
        if (kept != saved).any():
            raise RuntimeError('changing only {} moved the velocity of {} nodes on restart'.format(
                what, int((kept != saved).any(axis=1).sum())))
    change = dict(vbc_x1='vvf', vbc_val_x1='[0, {!r}, 0]'.format(float(mapped[node, 1])),
                  vbc_val_x1_l='0', vbc_n0='ftf', vbc_val_n0='[0, -2e4, 0]')
    restart_from_1(exe, folder, result, 'changed', change)
    restart_from_1(exe, bare, result, 'changed', change)
    count = same_outputs(bare, folder, 'changed')
    return 'x1 holding the y velocity of the edge map at node {}: restart matches the one ' \
           'without the record ({} files); ones changing only the n0 traction or n1 keep the ' \
           'velocity'.format(int(node), count)


def compare_letters(current, spec, reference):
    if current.returncode != spec['exit_code']:
        raise RuntimeError('exited {} where {} was expected'.format(current.returncode, spec['exit_code']))
    if current.returncode:
        return 'rejected (exit {})'.format(current.returncode)
    if not check_table(current.stdout, spec['ndims'], current.cfg):
        raise RuntimeError('no per-side table printed')
    if reference is None:
        return 'runs, table as expected'
    # The record of the input may differ where the two inputs differ without effect, as a
    # normal held at 0 does that follows the period table as letters and not as code 4 or 5.
    count = same_outputs(reference.folder, current.folder, ignored=RECORD)
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
            print(c['name'] + ('  (now rejected)' if c['current_exit'] and not c['exit_code'] else '') +
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
            if spec['restart']:
                verdict += '; ' + same_after_restart(current_exe, current)
            if spec['changed_restart']:
                verdict += '; ' + same_changed_restart(parent_exe, current_exe, parent, current, spec)
            if spec['edge_restart']:
                verdict += '; ' + edge_restart(current_exe, current)
            if spec['corner']:
                verdict += '; ' + corner_normal_velocity(current, spec)
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
