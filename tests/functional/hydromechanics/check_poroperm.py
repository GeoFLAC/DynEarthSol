#!/usr/bin/env python3
"""Verify evolving porosity/permeability and exact hydraulic restart/remeshing.

Usage: check_poroperm.py EXECUTABLE OUTPUT_DIRECTORY
Reuses the existing DES hydraulic PT fixtures and restart comparator.
"""

import argparse
import configparser
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from Dynearthsol import Dynearthsol


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--permeability-model', type=int, choices=(1, 2), default=1)
    parser.add_argument('--exponent', type=float, default=3)
    args = parser.parse_args()
    exe, output = args.executable.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    summary = {}
    for mode in ('fixed', 'remesh', 'thermal-matrix'):
        cfg = configparser.ConfigParser()
        cfg.optionxform = str
        cfg.read(ROOT / 'tests/functional/pt_remesh_restart.cfg')
        cfg['mat'].update(porosity_model='1', permeability_model=str(args.permeability_model))
        if args.permeability_model == 2:
            cfg['mat']['permeability_exponent'] = repr([args.exponent])
        cfg['control']['PT_absolute_tolerance'] = '1e-8'
        if mode != 'remesh':
            cfg['control']['has_moving_mesh'] = 'false'
        if mode == 'thermal-matrix':
            cfg['control']['has_thermal_pressurization'] = 'true'
            cfg['bc']['surface_temperature'] = '283'
            cfg['mat'].update(alpha='[3e-5]', matrix_storage='[1e-9]', matrix_exchange='[1e-12]')
        source = output / f'{mode}.cfg'
        with source.open('w') as stream:
            cfg.write(stream)
        case = output / mode
        command = [sys.executable, str(ROOT / 'tests/functional/pt_initial_checkpoint.py'),
                   str(exe), str(source), str(case), '--restart-frame', '1', '--steps', '4']
        if mode == 'remesh':
            command.append('--require-remesh')
        subprocess.run(command, env=dict(os.environ, OMP_NUM_THREADS='1'), check=True, timeout=240)
        model = Dynearthsol(str(case / 'fresh'))
        phi0 = json.loads(cfg['mat']['porosity'])[0]
        k0 = json.loads(cfg['mat']['hydraulic_perm'])[0]
        alpha = json.loads(cfg['mat']['biot_coeff'])[0]
        grain = json.loads(cfg['mat']['bulk_modulus_s'])[0]
        cells = model.read_field(0, 'connectivity')
        strain0 = model.read_field(0, 'strain')[:, :model.ndims].sum(axis=1)
        pressure0 = model.read_field(0, 'pore pressure')[cells].mean(axis=1)
        max_phi_error = 0.0
        for frame in model.frames:
            phi = model.read_field(frame, 'porosity')
            perm = model.read_field(frame, 'permeability')
            assert np.isfinite(phi).all() and ((phi > 0) & (phi < 1) & (phi <= alpha)).all()
            if args.permeability_model == 2:
                expected_k = k0 * (phi / phi0)**args.exponent
            else:
                expected_k = k0 * (phi / phi0)**3 * ((1 - phi0) / (1 - phi))**2
            np.testing.assert_allclose(perm, expected_k, atol=0, rtol=2e-14)
            if mode == 'fixed':
                strain = model.read_field(frame, 'strain')[:, :model.ndims].sum(axis=1) - strain0
                dp = model.read_field(frame, 'pore pressure')[cells].mean(axis=1) - pressure0
                expected_phi = alpha - (alpha - phi0) * np.exp(-strain - dp / grain)
                max_phi_error = max(max_phi_error, float(np.max(np.abs(phi - expected_phi))))
                np.testing.assert_allclose(phi, expected_phi, atol=5e-13, rtol=0)
        if mode == 'thermal-matrix':
            assert np.isfinite(model.read_field(model.frames[-1], 'matrix pressure')).all()
        assert np.max(np.abs(phi - phi0)) > 1e-8, 'Loaded porosity did not evolve'
        summary[mode] = dict(ndims=int(model.ndims), frames=len(model.frames),
                             phi_min=float(phi.min()), phi_max=float(phi.max()),
                             max_phi_error=max_phi_error if mode == 'fixed' else None,
                             exact_restart=True, changed_mesh=mode == 'remesh')
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print('PASS: Biot porosity, permeability law', args.permeability_model,
          'exact restart and remeshing')


if __name__ == '__main__':
    main()
