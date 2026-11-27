#!/usr/bin/env python3
"""DES fluid injection into a poroelastic borehole (plane strain).

This refinement study reports numerical errors against the analytical solution.

A quarter annulus a = 1 m, R = 20 m, refined out to 5 m, is unstressed
at p = 0. At t = 0+ the wall takes fluid at p_w: pore
pressure p_w (a 1 ms ramp of the scheduled boundary pressure) and normal
traction -p_w. The outer arc is sealed and traction-free (DES's input
check allows prescribed pressures on only one of two intersecting box
faces, which the annulus maps to x1 and z1); the symmetry cuts are
rollers and sealed. Gravity and heat are off; pressure feedback
is on. Compared with benchmarks/poromechanics/borehole_reference.py,
which includes the finite-annulus coupling: nodal pressure and element
radial/hoop total stress for 1.05 a <= r <= 4 a at four times, relative to
p_w.

    check_borehole.py EXECUTABLE OUTPUT [--spacings 0.2 0.1 0.05] [--steps 144]
"""
import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT), str(ROOT/'benchmarks'/'poromechanics'), str(HERE)]
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames
from borehole_reference import Borehole

A, R, RING = 1., 20., 5.
PW, DURATION, RAMP = 1e6, 36.0, 1e-3
MAT = dict(bulk=2e8, shear=1e8, biot=1.0, porosity=.2, fluid_bulk=2e8, grain_bulk=1e20,
           permeability=1e-12, viscosity=1e-3)
REF = Borehole(radius=A, outer=R, bulk=MAT['bulk'], shear=MAT['shear'], biot=MAT['biot'],
               storage_modulus=1/(MAT['porosity']/MAT['fluid_bulk']),
               mobility=MAT['permeability']/MAT['viscosity'], wall_pressure=PW, sealed=True)
BAND = (1.05*A, 4*A)


from des_case import digest, material_options, provenance, read_config, reset_boundaries, run_des, write_json, write_poly


def polygon(path, spacing, refined_radius=None):
    outer_n = math.ceil(np.pi*R/(8*spacing))
    inner_n = math.ceil(np.pi*A/(2*spacing))
    outer = np.linspace(0, np.pi/2, outer_n+1)
    inner = np.linspace(np.pi/2, 0, inner_n+1)
    points = [(R*np.cos(t), -R*np.sin(t)) for t in outer]
    if refined_radius is not None:
        if not A < refined_radius < R/3:
            raise ValueError('Refined ring must lie inside the outer region seed')
        points.append((0., -refined_radius))
    points += [(A*np.cos(t), -A*np.sin(t)) for t in inner]
    if refined_radius is not None: points.append((refined_radius, 0.))
    markers = [32]*outer_n + [1]*(1+(refined_radius is not None)) + [2]*inner_n + [16]*(1+(refined_radius is not None))
    assert len(points) == len(markers)
    edges = [(i, (i+1)%len(points), marker) for i, marker in enumerate(markers)]
    if refined_radius is not None:
        # A constrained, unlabelled interior arc grades mesh size without a
        # material or boundary-condition interface. Both regions are material 0.
        angles = np.linspace(np.pi/2, 0, math.ceil(np.pi*refined_radius/(2*spacing))+1)
        ids = [outer_n+1]
        for angle in angles[1:-1]:
            ids.append(len(points)); points.append((refined_radius*np.cos(angle), -refined_radius*np.sin(angle)))
        ids.append(len(markers)-1)
        edges += [(i, j, 0) for i, j in zip(ids[:-1], ids[1:])]
    regions = [(R/3, -R/3, 0, 16*spacing**2)]
    if refined_radius is not None:
        inner_seed = (A+refined_radius)/(2*np.sqrt(2))
        regions.append((inner_seed, -inner_seed, 0, spacing**2/2))
    write_poly(path, points, edges, regions)


def polar(values, normal):
    """Radial, hoop and shear components of (xx, zz, xz) tensors for unit radial directions normal."""
    nx, nz = normal.T
    rr = nx**2*values[:, 0]+nz**2*values[:, 1]+2*nx*nz*values[:, 2]
    hoop = nz**2*values[:, 0]+nx**2*values[:, 1]-2*nx*nz*values[:, 2]
    shear = nx*nz*(values[:, 1]-values[:, 0])+(nx**2-nz**2)*values[:, 2]
    return rr, hoop, shear


def config(h, steps):
    cfg = read_config(ROOT / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
    cfg['sim'].update(modelname='hole', max_steps=str(steps), max_time_in_yr=repr(2*DURATION/YEAR2SEC),
                      output_step_interval=str(steps//4), output_time_interval_in_yr='1e9',
                      checkpoint_frame_interval='100000', has_initial_checkpoint='no')
    cfg['mesh'] = dict(meshing_option='90', poly_filename='annulus.poly', xlength=repr(R), ylength='1', zlength=repr(R),
                       resolution=repr(h), remeshing_option='0', min_angle='28', max_steiner_factor='0')
    cfg['control'].update(gravity='0', has_thermal_diffusion='false', has_moving_mesh='false',
                          fixed_dt=repr(DURATION/steps), has_poroelastic_pressure_feedback='true',
                          is_using_mixed_stress='false', PT_option='1', PT_max_iter='200000',
                          PT_relative_tolerance='1e-9', PT_absolute_tolerance='1e-8')
    bc = cfg['bc']
    reset_boundaries(bc, vbc='0', vbc_val='0', stress_bc='0', stress_val='0', hbc='0')
    bc.update(vbc_x0='1', vbc_z0='1', stress_bc_x1='4', stress_val_x1=repr(-PW), hbc_x1='2',
              hbc_val_x1=repr(PW), num_hbc_period_x1='2',
              hbc_period_x1_time_in_yr=f'[0, {RAMP/YEAR2SEC!r}]', hbc_period_x1_ratio='[0, 1]')
    cfg['ic'].update(pore_pressure_option='0', has_body_force_adjustment='false', weakzone_option='0')
    cfg['mat'].update(material_options(MAT), rheology_type='elastic',
                      is_plane_strain='true', alpha='[0]', fluid_alpha='[0]')
    return cfg


def compare(folder):
    model = Dynearthsol(str(folder/'hole')); rows = scan_frames(str(folder/'hole'))
    x = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    rn = np.linalg.norm(x, axis=1)
    centre = x[conn].mean(axis=1); re = np.linalg.norm(centre, axis=1)
    nodes, elems = (rn >= BAND[0]) & (rn <= BAND[1]), (re >= BAND[0]) & (re <= BAND[1])
    out = []
    for row in rows[1:]:
        t = row['time']
        p = model.read_field(row['frame'], 'pore pressure')
        s = model.read_field(row['frame'], 'stress')
        srr, stt, _ = polar(s, centre/re[:, None])
        fn, fe = REF.fields(rn[nodes], t), REF.fields(re[elems], t)
        out.append(dict(time=t, pressure_error=float(np.max(np.abs(p[nodes] - fn['p'])))/PW,
                        pressure_rms=float(np.sqrt(np.mean((p[nodes] - fn['p'])**2)))/PW,
                        radial_error=float(np.max(np.abs(srr[elems] - fe['srr'])))/PW,
                        hoop_error=float(np.max(np.abs(stt[elems] - fe['stt'])))/PW,
                        hoop_rms=float(np.sqrt(np.mean((stt[elems] - fe['stt'])**2)))/PW))
    return out, len(x), len(conn)


def study(exe, out, spacings, steps):
    out.mkdir(parents=True, exist_ok=False)
    result = dict(provenance(exe, __file__, out),
                  steps=steps, material=MAT, wall_pressure=PW, diffusivity=REF.c, cases=[])
    for h in spacings:
        folder = out/f'h{h:g}'; folder.mkdir()
        polygon(folder/'annulus.poly', h, refined_radius=RING)
        run_des(exe, folder, 'hole', config(h, steps), timeout=40000)
        frames, nodes, elements = compare(folder)
        item = dict(spacing_m=h, nodes=nodes, elements=elements, frames=frames, input_sha256=digest(folder/'hole.cfg'),
                    mesh_sha256=digest(folder/'annulus.poly'),
                    pressure_error=max(f['pressure_error'] for f in frames),
                    stress_error=max(max(f['radial_error'], f['hoop_error']) for f in frames))
        result['cases'].append(item)
        write_json(out/'summary.json', result)
        print(json.dumps({k: v for k, v in item.items() if k != 'frames'}), flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--spacings', type=float, nargs='+', default=[0.2, 0.1, 0.05])
    p.add_argument('--steps', type=int, default=144)
    a = p.parse_args()
    study(a.executable.resolve(), a.output.resolve(), a.spacings, a.steps)
