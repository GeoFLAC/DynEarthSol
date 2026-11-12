#!/usr/bin/env python3
"""Do monitor points follow the material through deformation, remeshing and restart?

An elastic box (2D, no gravity) is stretched: x0 fixed, x1 pulled at V,
rollers on the base, free top. Mesh quality is checked every step with a high
threshold, so the mesh is rebuilt several times, and DES saves a frame every
step and just before and after each remesh. The reference is rebuilt from those
frames independently of the monitor: a material point keeps its barycentric
coordinates in its element while the mesh moves, and keeps its position when
the mesh is rebuilt; nodal values are interpolated at the point.

For each point type (monitor.point_type: particle, the default, interpolates
at the point; grid reports the node nearest the point and the element nearest
it), two runs:
  material: the full run.
  restart:  stopped halfway and restarted; the CSV keeps its earlier rows and
            continues the same series.
A particle must report the exact position and velocity; a grid point must
report the node nearest the exact point in that step's mesh, and its velocity.

A third pair checks the rows written at remeshes (monitor.remesh_output = both,
monitored every 4 steps): a row before each remesh between monitor steps, a row
after every remesh with the remapped values, and a restart from a remesh step
that reproduces the uninterrupted CSV row for row.
Errors are in units of the mesh spacing h and of V.

    check_monitor_tracking.py EXECUTABLE OUTPUT [--tolerance 1e-6]
"""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT), str(ROOT/'tests/functional')]
from Dynearthsol import Dynearthsol
from des_test_utils import read_config, run_des

LX, LZ, H_MESH = 10.0, 5.0, 0.5
STEPS, V = 60, 10.0/60                    # m/s: the box doubles its length
POINTS = [(2.3, -1.7), (5.1, -3.3), (7.7, -0.9)]


def config(name, steps, point_type, restart_from=None, interval=1, remesh_output=None):
    cfg = read_config(ROOT/'tests/functional/pt_hydrostatic_equilibrium.cfg')
    cfg['sim'].update(modelname=name, max_steps=str(steps), max_time_in_yr='1', output_step_interval='1',
                      output_time_interval_in_yr='1e9', checkpoint_frame_interval='1', has_initial_checkpoint='yes',
                      has_output_during_remeshing='yes')
    if restart_from:
        cfg['sim'].update(is_restarting='yes', restarting_from_modelname=restart_from[0],
                          restarting_from_frame=str(restart_from[1]))
    cfg['mesh'] = dict(meshing_option='1', xlength=repr(LX), ylength='1', zlength=repr(LZ), resolution=repr(H_MESH),
                       remeshing_option='0', quality_check_step_interval='1', min_quality='0.6', min_angle='32')
    cfg['control'].update(gravity='0', fixed_dt='1', has_PT='yes', has_moving_mesh='yes', has_thermal_diffusion='no',
                          has_hydraulic_diffusion='no', PT_relative_tolerance='1e-9', PT_absolute_tolerance='1e-8',
                          PT_max_iter='200000')
    bc = cfg['bc']
    for side in ('x0', 'x1', 'y0', 'y1', 'z0', 'z1'):
        bc.update({f'vbc_{side}': '0', f'vbc_val_{side}': '0', f'stress_bc_{side}': '0', f'hbc_{side}': '0'})
    bc.update(vbc_x0='1', vbc_x1='1', vbc_val_x1=repr(V), vbc_z0='1', has_water_loading='no', has_winkler_foundation='no')
    cfg['ic'].update(initial_stress_option='1', initial_stress='[0,0,0]', has_body_force_adjustment='false',
                     weakzone_option='0', mattype_option='0')
    cfg['mat'].update(rheology_type='elastic', is_plane_strain='true', num_materials='1', rho0='[2650]',
                      bulk_modulus='[1e10]', shear_modulus='[1e10]')
    cfg['monitor'] = dict(enabled='yes', step_interval=str(interval), num_points=str(len(POINTS)),
                          points_x='[' + ','.join(repr(p[0]) for p in POINTS) + ']',
                          points_z='[' + ','.join(repr(p[1]) for p in POINTS) + ']',
                          output_prefix='monitor', output_coord='yes', output_velocity='yes')
    cfg['monitor']['point_type'] = point_type
    if remesh_output: cfg['monitor']['remesh_output'] = remesh_output
    return cfg


def run(exe, folder, name, cfg):
    return run_des(exe, folder, name, cfg, timeout=900).stdout.lower().count('remesh')


def frames(models):
    """Every saved frame in order: (step, coordinates, connectivity, velocity)."""
    out = []
    for model in models:
        for i, frame in enumerate(model.frames):
            if out and model.steps[i] == out[-1][0] and frame == model.frames[0] and model is not models[0]:
                continue                  # a restart's first frame repeats the checkpoint frame
            out.append((model.steps[i], model.read_field(frame, 'coordinate'), model.read_field(frame, 'connectivity'),
                        model.read_field(frame, 'velocity')))
    return out


def locate(x, conn, q):
    """Containing element and barycentric weights of q (triangles)."""
    a, b, c = x[conn[:, 0]], x[conn[:, 1]], x[conn[:, 2]]
    det = (b[:, 0] - a[:, 0])*(c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0])*(b[:, 1] - a[:, 1])
    w1 = ((q[0] - a[:, 0])*(c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0])*(q[1] - a[:, 1]))/det
    w2 = ((b[:, 0] - a[:, 0])*(q[1] - a[:, 1]) - (q[0] - a[:, 0])*(b[:, 1] - a[:, 1]))/det
    w0 = 1 - w1 - w2
    inside = np.where((w0 >= -1e-12) & (w1 >= -1e-12) & (w2 >= -1e-12))[0]
    e = int(inside[0]); return e, np.array([w0[e], w1[e], w2[e]])


def references(models, query, after=False):
    """Per step: the point's position and interpolated velocity, and the nearest node and its velocity.
    The step's first frame is the state before any remesh; with after, its last frame, the remapped one."""
    out = {}; elem = weights = None; position = np.array(query, float); previous = None
    for step, x, conn, v in frames(models):
        same = previous is not None and previous.shape == conn.shape and np.array_equal(previous, conn)
        if same:
            position = weights @ x[conn[elem]]        # carried by the moving mesh
        elem, weights = locate(x, conn, position)    # relocated in a new or moved mesh at the same position
        previous = conn
        if step not in out or after:
            n = int(np.argmin(np.hypot(x[:, 0] - position[0], x[:, 1] - position[1])))
            out[step] = (position.copy(), weights @ v[conn[elem]], x[n], v[n])
    return out


def remesh_steps(model):
    """Steps with frames written before and after a remesh (three frames at one step)."""
    steps, counts = np.unique(model.steps, return_counts=True)
    return [int(s) for s, c in zip(steps, counts) if c >= 3]


def rows(path):
    with open(path) as f: return list(csv.DictReader(f))


def errors(folder, models, grid):
    result = []
    for i, query in enumerate(POINTS):
        ref = references(models, query); post = references(models, query, after=True)
        worst_x = worst_v = 0.0; steps = []
        for row in rows(folder/f'monitor_point_{i}.csv'):
            step = int(row['step']); steps.append(step)
            kind = int(row['row'])
            position, velocity, node, node_velocity = (post if kind == 2 else ref)[step]
            ex, ev = (node, node_velocity) if grid else (position, velocity)
            got = (float(row['coord_x']), float(row['coord_z'])); vel = (float(row['velocity_x']), float(row['velocity_z']))
            worst_x = max(worst_x, np.hypot(got[0] - ex[0], got[1] - ex[1])/H_MESH)
            worst_v = max(worst_v, np.hypot(vel[0] - ev[0], vel[1] - ev[1])/V)
        result.append(dict(point=i, rows=len(steps), first_step=min(steps), last_step=max(steps),
                           position_error_over_h=worst_x, velocity_error_over_V=worst_v))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--tolerance', type=float, default=1e-6)
    a = p.parse_args()
    exe = a.executable.expanduser().resolve()
    out = a.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    report = dict(executable=str(exe), tolerance=a.tolerance, cases={})
    restarts = []
    for point_type in ('particle', 'grid'):
        folder = out/f'{point_type}-material'; folder.mkdir()
        remeshes = run(exe, folder, 'box', config('box', STEPS, point_type))
        report['cases'][f'{point_type}-material'] = dict(remesh_log_lines=remeshes, points=errors(
            folder, [Dynearthsol(str(folder/'box'))], point_type == 'grid'))
        folder = out/f'{point_type}-restart'; folder.mkdir(); half = STEPS//2
        run(exe, folder, 'first', config('first', half, point_type))
        run(exe, folder, 'second', config('second', STEPS, point_type, restart_from=('first', half)))
        points = errors(folder, [Dynearthsol(str(folder/'first')), Dynearthsol(str(folder/'second'))], point_type == 'grid')
        for q in points: q['kept_rows_before_restart'] = q['first_step'] == 0 and q['rows'] >= STEPS + 1
        report['cases'][f'{point_type}-restart'] = dict(points=points); restarts += points
    passed = all(q['position_error_over_h'] <= a.tolerance and q['velocity_error_over_V'] <= a.tolerance
                 for c in report['cases'].values() for q in c['points']) and \
             all(q['kept_rows_before_restart'] for q in restarts)
    folder = out/'remesh-rows-material'; folder.mkdir()
    run(exe, folder, 'box', config('box', STEPS, 'particle', interval=4, remesh_output='both'))
    model = Dynearthsol(str(folder/'box')); remeshes = remesh_steps(model)
    points = errors(folder, [model], False)
    kinds = {}
    for r in rows(folder/'monitor_point_0.csv'):
        kinds.setdefault(int(r['step']), []).append(int(r['row']))
    regular_ok = all(k != 0 or s % 4 == 0 for s, ks in kinds.items() for k in ks)
    before_ok = all((1 in kinds.get(s, [])) == (s % 4 != 0) for s in remeshes)
    after_ok = all(kinds.get(s, [])[-1:] == [2] for s in remeshes)
    start = next(s for s in remeshes if s >= STEPS//3 and s % 4)
    folder = out/'remesh-rows-restart'; folder.mkdir()
    run(exe, folder, 'first', config('first', start, 'particle', interval=4, remesh_output='both'))
    run(exe, folder, 'second', config('second', STEPS, 'particle', restart_from=('first', start),
                                               interval=4, remesh_output='both'))
    same = True
    for i in range(len(POINTS)):
        whole = rows(out/'remesh-rows-material'/f'monitor_point_{i}.csv')
        resumed = rows(folder/f'monitor_point_{i}.csv')
        keys = lambda rr: [(r['step'], r['row']) for r in rr]
        same &= keys(whole) == keys(resumed) and all(
            abs(float(u[c]) - float(w[c])) <= 1e-9*max(1.0, abs(float(w[c])))
            for u, w in zip(resumed, whole) for c in w if c not in ('matched_node', 'matched_elem'))
    report['remesh_rows'] = dict(remesh_steps=remeshes, restart_step=start, points=points,
                                 regular_rows_on_interval=regular_ok, before_rows=before_ok,
                                 after_rows=after_ok, restart_reproduces_csv=bool(same))
    passed = passed and regular_ok and before_ok and after_ok and same and all(
        q['position_error_over_h'] <= a.tolerance and q['velocity_error_over_V'] <= a.tolerance for q in points)
    report['passed'] = passed
    (out/'report.json').write_text(json.dumps(report, indent=1))
    for name, c in report['cases'].items():
        for q in c['points']:
            print(f"{name:16s} point {q['point']}: rows {q['rows']} (steps {q['first_step']}-{q['last_step']}), "
                  f"position {q['position_error_over_h']:.2e} h, velocity {q['velocity_error_over_V']:.2e} V"
                  + (f", earlier rows kept: {q['kept_rows_before_restart']}" if 'kept_rows_before_restart' in q else ''))
    r = report['remesh_rows']
    print(f"remesh rows: remeshes at {r['remesh_steps']}, restart at {r['restart_step']}; regular rows on interval "
          f"{r['regular_rows_on_interval']}, before rows {r['before_rows']}, after rows {r['after_rows']}, "
          f"restart reproduces CSV {r['restart_reproduces_csv']}")
    for q in r['points']:
        print(f"remesh-rows      point {q['point']}: position {q['position_error_over_h']:.2e} h, "
              f"velocity {q['velocity_error_over_V']:.2e} V")
    print('PASSED' if passed else 'FAILED')
    sys.exit(0 if passed else 1)


if __name__ == '__main__':
    main()
