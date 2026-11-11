#!/usr/bin/env python3
"""Rapid drawdown of a reservoir against the embankment: factor of safety in time.

The embankment of strength_reduction_slope.py, saturated (porosity 0.3, Gardner
retention, Bishop chi = S_e). A reservoir at the crest loads the face while the
slope fills to hydrostatic; at t = 0 it drops to the toe platform and the face
and toe drain (elastic, all explicit). At chosen times strength-reduction trials
restart with the pore pressure held (permeability frozen, compressible pore
fluid). Reference: Bishop's method on DES's pore pressure and unit weight.

    drawdown_slope.py EXECUTABLE OUTPUT [--spacings 1] [--times 250 1000 ...] [--factors 0.8 0.9 ...]
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames
import bishop_reference as br
import strength_reduction_slope as srm

G_ACC, RHO_F, PHI_N = srm.G_ACC, 1000.0, 0.3
RHO_S = (br.GAMMA/G_ACC - RHO_F*PHI_N)/(1 - PHI_N)        # saturated unit weight = br.GAMMA
GAMMA_W = RHO_F*G_ACC
ALPHA_G, SR = 2.4e-4, 0.05                                # Gardner (1/Pa), residual saturation
PERMEABILITY, VISCOSITY, FLUID_BULK = 1e-13, 1e-3, 2.2e9  # drainage ~8000 s, slow against the mechanics (~100 s)
FROZEN = 1e-25                                            # m^2: no flow during a strength-reduction trial
FREE_FLUID = 1e3                                          # Pa: trial fluid bulk modulus; deformation leaves p unchanged
FILL_TIME, DRAIN_TIME, FRAME = 25000.0, 20000.0, 250.0     # s; saved (and checkpointed) every FRAME seconds
SEED_DT = 1e-3                                            # s per metre: seed step, below every explicit limit
TIMES = [250.0, 1000.0, 3000.0, 6000.0, 12000.0, 20000.0]  # s after drawdown
FACTORS = [0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4]


sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tests/functional/hydromechanics'))
from des_case import digest, reset_boundaries, source_revision


def polygon(path, h):
    """The embankment outline with the crest marked y0 (unused otherwise in 2D) so that the crest and the face can take
    different hydraulic conditions; the face and the toe platform stay z1."""
    srm.polygon(path, h)
    lines = path.read_text().splitlines()
    npts = int(lines[0].split()[0])
    pts = {int(l.split()[0]): (float(l.split()[1]), float(l.split()[2])) for l in lines[1:1 + npts]}
    out = lines[:2 + npts]
    for l in lines[2 + npts:2 + 2*npts]:
        i, a, b, mk = (int(v) for v in l.split())
        crest = mk == 32 and pts[a][1] == 0.0 and pts[b][1] == 0.0
        out.append(f'{i} {a} {b} {4 if crest else mk}')
    out += lines[2 + 2*npts:]
    path.write_text('\n'.join(out) + '\n')


def config(h, name, stage, restart=None, factor=None, folder=None):
    """stage: seed, reservoir (gravity + water load, filling), drawdown, trial; every stage after the seed is explicit
    DES with the dry embankment's settings (DES refuses to switch hydraulic formulation, PT or explicit, mid-run).
    restart = (model, frame) in folder; DES's max_steps and max_time are absolute, so they count from that frame."""
    start = dict(steps=0, time=0.0)
    if restart: start = next(r for r in scan_frames(str(folder/restart[0])) if r['frame'] == restart[1])
    cfg = srm.config(h, name, factor if stage == 'trial' else None)
    sim, ctl, bc, mat = cfg['sim'], cfg['control'], cfg['bc'], cfg['mat']
    if restart: sim.update(is_restarting='yes', restarting_from_modelname=restart[0], restarting_from_frame=str(restart[1]))
    duration = dict(seed=0.0, reservoir=FILL_TIME, drawdown=DRAIN_TIME, trial=srm.DURATION)[stage]
    sim.update(max_steps='0' if stage == 'seed' else '1000000000', max_time_in_yr=repr((start['time'] + duration)/YEAR2SEC),
               output_step_interval='1000000000', checkpoint_frame_interval='1' if stage in ('reservoir', 'drawdown') else '100000',
               output_time_interval_in_yr=repr((FRAME if stage != 'trial' else srm.DURATION/srm.FRAMES)/YEAR2SEC))
    ctl.update(has_PT='false', fixed_dt=repr(SEED_DT*h) if stage == 'seed' else '0', inertial_scaling=repr(srm.SCALING),
               characteristic_speed='1e-5', dt_fraction='0.5', damping_factor=repr(srm.DAMPING),
               has_hydraulic_diffusion='true', has_poroelastic_pressure_feedback='true')
    reset_boundaries(bc, hbc='0')
    loaded = stage in ('seed', 'reservoir')
    # Reservoir: fill from p = 0 through the crest (p = 0) and the base (hydrostatic); face and toe loaded, sealed.
    # Drawdown: crest and base sealed; face and toe platform drain at p = 0.
    base = GAMMA_W*-br.BASE
    bc.update(has_water_loading='yes' if loaded else 'no', sea_water_density=repr(RHO_F),
              hbc_pressure=f'[0,0,0,0,{base if loaded else 0!r},0]',
              hbc_y0='2' if loaded else '0', hbc_z0='2' if loaded else '0', hbc_z1='0' if loaded else '2')
    ctl['surf_base_level'] = '0'
    k = FROZEN if stage == 'trial' else PERMEABILITY
    kf = FREE_FLUID if stage == 'trial' else FLUID_BULK
    mat.update(rho0=f'[{RHO_S!r}]', fluid_rho0=f'[{RHO_F!r}]', porosity=f'[{PHI_N!r}]', biot_coeff='[1]',
               fluid_bulk_modulus=f'[{kf!r}]', bulk_modulus_s='[1e20]', hydraulic_perm=f'[{k!r}]',
               fluid_visc=f'[{VISCOSITY!r}]', fluid_alpha='[0]', retention_model='[3]', retention_alpha=f'[{ALPHA_G!r}]',
               residual_saturation=f'[{SR!r}]')
    return cfg


def hydrostatic_misfit(folder, name):
    """Largest |p - rho_f g (0 - z)| at the last frame, relative to the hydrostatic pressure at the base."""
    model = Dynearthsol(str(folder/name)); row = scan_frames(str(folder/name))[-1]
    x = model.read_field(row['frame'], 'coordinate'); p = model.read_field(row['frame'], 'pore pressure')
    return float(np.abs(p + GAMMA_W*x[:, 1]).max()/(GAMMA_W*-br.BASE))


def frame_at(folder, name, t):
    rows = scan_frames(str(folder/name))
    t0 = rows[0]['time']
    return min(rows[1:], key=lambda r: abs(r['time'] - t0 - t))


def bishop(folder, name, frame):
    """Bishop's FS on DES's pore pressure (chi p) and unit weight at one saved frame."""
    from matplotlib.tri import Triangulation, LinearTriInterpolator
    model = Dynearthsol(str(folder/name))
    x = model.read_field(frame, 'coordinate'); conn = model.read_field(frame, 'connectivity')
    p = model.read_field(frame, 'pore pressure')
    se = np.where(p < 0, np.exp(ALPHA_G*np.minimum(p, 0)), 1.0)
    sat = SR + (1 - SR)*se
    tri = Triangulation(x[:, 0], x[:, 1], conn)
    chi_p = LinearTriInterpolator(tri, se*p); s = LinearTriInterpolator(tri, sat)
    def pore(xx, zz):
        return np.ma.filled(chi_p(xx, np.maximum(zz, br.BASE + 1e-9)), 0.0)
    def weight(xx, zz):
        return G_ACC*(RHO_S*(1 - PHI_N) + RHO_F*PHI_N*np.ma.filled(s(xx, zz), 1.0))
    slope = br.Slope(br.PROFILE, br.BASE, None, br.COHESION, br.FRICTION, pore=pore, unit_weight=weight, slices=200)
    fs, v = slope.critical(br.starts())
    return fs, list(slope.circle(*v)), float(p.max()), float(p.min())


def study(exe, out, spacings, times, factors, workers, flow_threads=1, resume=None):
    if workers < 1 or flow_threads < 1:
        raise ValueError('Worker and flow-thread counts must be positive')
    if resume is not None and len(spacings) != 1:
        raise ValueError('A reservoir checkpoint belongs to one mesh spacing')
    slots = os.environ.get('NSLOTS')
    if slots and max(flow_threads, min(workers, len(factors))) > int(slots):
        raise ValueError('Flow threads and simultaneous trials must fit the scheduler allocation')
    out.mkdir(parents=True, exist_ok=False)
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(__file__),
                  reference_sha256=digest(br.__file__), srm_sha256=digest(srm.__file__),
                  source_commit=source_revision(ROOT),
                  rho_s=RHO_S, porosity=PHI_N, frozen=FROZEN, free_fluid=FREE_FLUID, alpha=ALPHA_G, sr=SR, permeability=PERMEABILITY,
                  fill_time=FILL_TIME, drain_time=DRAIN_TIME, frame=FRAME, duration=srm.DURATION, scaling=srm.SCALING,
                  damping=srm.DAMPING, flow_threads=flow_threads, trial_workers=workers, cases=[])
    for h in spacings:
        folder = out/f'h{h:g}'; folder.mkdir()
        polygon(folder/'slope.poly', h)
        for name, stage, restart in (('seed', 'seed', None), ('reservoir', 'reservoir', ('seed', 0)), ('drawdown', 'drawdown', None)):
            if resume is not None and stage == 'seed': continue
            if name == 'drawdown': restart = ('reservoir', scan_frames(str(folder/'reservoir'))[-1]['frame'])
            if resume is not None and stage == 'reservoir':
                prefix, frame = resume
                cfg = config(h, name, stage, (prefix.name, frame), folder=prefix.parent)
                row = next(r for r in scan_frames(str(prefix)) if r['frame'] == frame)
                if not 0 < row['time'] < FILL_TIME:
                    raise ValueError('Reservoir continuation must precede the original filling endpoint')
                cfg['sim'].update(restarting_from_modelname=str(prefix), max_time_in_yr=repr(FILL_TIME/YEAR2SEC))
                checkpoint = prefix.with_name(f'{prefix.name}.chkpt.{frame:06d}')
                if not checkpoint.is_file(): checkpoint = checkpoint.with_suffix(checkpoint.suffix+'.vtkhdf')
                result['reservoir_continuation'] = dict(frame=frame, time_s=row['time'],
                    checkpoint_sha256=digest(checkpoint))
            else:
                cfg = config(h, name, stage, restart, folder=folder)
            code = srm.run(exe, folder, name, cfg, threads=flow_threads)
            if code: raise RuntimeError(f'{folder}/{name}: exit {code}')
        case = dict(spacing_m=h, hydrostatic_misfit=hydrostatic_misfit(folder, 'reservoir'), states=[])
        print(json.dumps(dict(spacing_m=h, hydrostatic_misfit=case['hydrostatic_misfit'])), flush=True)
        for t in times:
            row = frame_at(folder, 'drawdown', t)
            fs_bishop, circle, pmax, pmin = bishop(folder, 'drawdown', row['frame'])
            start = ('drawdown', row['frame'])
            names = [f'srf{t:g}-{f:.3f}' for f in factors]
            with ThreadPoolExecutor(workers) as pool:
                codes = list(pool.map(lambda nf: srm.run(exe, folder, nf[0], config(h, nf[0], 'trial', start, nf[1], folder)),
                                      zip(names, factors)))
            state = dict(time=row['time'] - scan_frames(str(folder/'drawdown'))[0]['time'], frame=row['frame'],
                         fs_bishop=fs_bishop, circle=circle, p_max=pmax, p_min=pmin, trials=[])
            for name, f, code in zip(names, factors, codes):
                if code: raise RuntimeError(f'{folder}/{name}: exit {code}')
                frames, _, _ = srm.analyse(folder, name)
                state['trials'].append(dict(factor=f, at_rest=srm.at_rest(frames), deceleration=srm.deceleration(frames),
                                            frames=frames))
            case['states'].append(state)
            rest = [tr['factor'] for tr in state['trials'] if tr['at_rest']]
            print(json.dumps(dict(spacing_m=h, time=state['time'], fs_bishop=fs_bishop, largest_at_rest=max(rest, default=None),
                                  at_rest=[tr['at_rest'] for tr in state['trials']])), flush=True)
            result['cases'] = [c for c in result['cases'] if c['spacing_m'] != h] + [case]
            (out/'summary.json').write_text(json.dumps(result, indent=1, allow_nan=False))
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--spacings', type=float, nargs='+', default=[1.0])
    p.add_argument('--times', type=float, nargs='+', default=TIMES)
    p.add_argument('--factors', type=float, nargs='+', default=FACTORS)
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--flow-threads', type=int, default=1, help='OpenMP threads during reservoir filling and drawdown; trials use one thread each')
    p.add_argument('--resume-reservoir', type=Path, help='Continue a saved reservoir prefix into a new output folder')
    p.add_argument('--resume-frame', type=int, help='Completed checkpoint frame of --resume-reservoir')
    a = p.parse_args()
    if (a.resume_reservoir is None) != (a.resume_frame is None):
        p.error('--resume-reservoir and --resume-frame must be supplied together')
    resume = (a.resume_reservoir.resolve(), a.resume_frame) if a.resume_reservoir else None
    study(a.executable.resolve(), a.output.resolve(), a.spacings, a.times, a.factors, a.workers, a.flow_threads, resume)
