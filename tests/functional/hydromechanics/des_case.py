"""Provenance, boundary, solver and analysis helpers for the hydromechanics checks, on top of
the configuration, solver and output helpers of tests/functional/des_test_utils.py."""
import copy
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'tests/functional'))
from des_test_utils import (digest, lumped_volumes, read_config, run_des,  # noqa: E402,F401  (read_config for callers)
                            simplex_gradients, simplex_volumes, write_config)

SIDES = ('x0', 'x1', 'y0', 'y1', 'z0', 'z1')
# [mat] option of each poroelastic property, named as in the reference modules' constructors.
MATERIAL_KEYS = dict(bulk='bulk_modulus', shear='shear_modulus', porosity='porosity', biot='biot_coeff',
                     fluid_bulk='fluid_bulk_modulus', grain_bulk='bulk_modulus_s',
                     permeability='hydraulic_perm', viscosity='fluid_visc')


def file_hashes(paths, root=None):
    """sha256 of each file in paths, keyed by its name or, with root, by its path relative to root."""
    return {str(p.relative_to(root)) if root else p.name: digest(p) for p in sorted(paths) if p.is_file()}


def check_hashes(hashes, root):
    """Assert that each file named in hashes, relative to root, still has its recorded sha256."""
    for name, sha in hashes.items():
        assert digest(Path(root)/name) == sha, name


def write_json(path, value):
    """Write value to path as indented JSON, refusing NaN and infinity."""
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def material_options(props):
    """[mat] lists of one material for the properties in props that MATERIAL_KEYS names; others are skipped."""
    return {MATERIAL_KEYS[name]: f'[{value!r}]' for name, value in props.items() if name in MATERIAL_KEYS}


def material_properties(cfg):
    """The first material's MATERIAL_KEYS properties in cfg [mat]."""
    return {name: float(json.loads(cfg['mat'][key])[0]) for name, key in MATERIAL_KEYS.items()}


def write_poly(path, points, edges, regions):
    """Triangle .poly file of points (x, z), edges (start, end, marker), no holes and regions
    (x, z, attribute, maximum area)."""
    with Path(path).open('w') as stream:
        stream.write(f'{len(points)} 2 0 0\n')
        for i, (x, z) in enumerate(points):
            stream.write(f'{i} {x:.17g} {z:.17g}\n')
        stream.write(f'{len(edges)} 1\n')
        for i, (start, end, marker) in enumerate(edges):
            stream.write(f'{i} {start} {end} {marker}\n')
        stream.write(f'0\n{len(regions)}\n')
        for i, (x, z, attribute, area) in enumerate(regions):
            stream.write(f'{i} {x} {z} {attribute} {area}\n')


def reset_boundaries(bc, **values):
    for side in SIDES:
        for name, value in values.items():
            bc[f'{name}_{side}'] = value


def source_revision(root=ROOT):
    return subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip()


def source_diff(root=ROOT):
    return subprocess.check_output(['git', 'diff', 'HEAD', '--binary'], cwd=root)


def provenance(executable, runner, output=None):
    """Executable path and sha256, runner sha256, source commit; writes output/source.diff if given."""
    if output is not None:
        (Path(output)/'source.diff').write_bytes(source_diff())
    return dict(executable=str(executable), executable_sha256=digest(executable),
                runner_sha256=digest(runner), source_commit=source_revision())


def check_manifest(manifest, binary_hash):
    """Parsed binary manifest, or None, after asserting it lists binary_hash (string or entry sha256)."""
    if not manifest:
        return None
    solver = json.loads(Path(manifest).read_text())
    hashes = {value['sha256'] if isinstance(value, dict) else value for value in solver['binaries'].values()}
    assert binary_hash in hashes, (manifest, binary_hash)
    return solver


def run_solver(command, *, log, **options):
    with Path(log).open('w') as stream:
        return subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, **options)


def run_model(executable, folder, cfg, name, *, timeout=600):
    """Run cfg as model name in folder, logging to name.log, and open its output."""
    from Dynearthsol import Dynearthsol
    cfg = copy.deepcopy(cfg)
    cfg['sim']['modelname'] = name
    run_des(executable, folder, name, cfg, timeout=timeout)
    return Dynearthsol(str(folder/name))


def expect_exit(executable, folder, name, cfg, code, timeout=120, message=None):
    """Run DES expecting exit status code and message in the log; return the run, log as stdout."""
    result = run_des(executable, folder, name, cfg, timeout, check=False)
    assert result.returncode == code, (str(folder), name, result.returncode, code, result.stdout[-2500:])
    assert message is None or message in result.stdout, (str(folder), name, message, result.stdout[-2500:])
    return result


def restart_check(executable, cfg_or_path, case, steps, restart_frame, require_remesh=False,
                  require_stationary=False, timeout=240):
    """Run pt_initial_checkpoint.py into new folder case; a cfg object goes to case.cfg, the log to case.log."""
    case = Path(case)
    if isinstance(cfg_or_path, (str, Path)):
        source = Path(cfg_or_path)
    else:
        source = case.parent/(case.name+'.cfg')
        write_config(source, cfg_or_path)
    command = [sys.executable, str(ROOT/'tests/functional/pt_initial_checkpoint.py'), str(executable),
               str(source), str(case), '--steps', str(steps), '--restart-frame', str(restart_frame)]
    command += ['--require-remesh']*require_remesh + ['--require-stationary']*require_stationary
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    return run_solver(command, log=case.parent/(case.name+'.log'), env=env, check=True, timeout=timeout)


def retention_options(*laws):
    """[mat] retention lists for one unsaturated_reference.Retention law per material."""
    from unsaturated_reference import MODELS
    return dict(retention_model=json.dumps([MODELS.index(law.model)+1 for law in laws]),
                retention_alpha=json.dumps([law.alpha for law in laws]),
                residual_saturation=json.dumps([law.sr for law in laws]),
                retention_n=json.dumps([law.n if law.n is not None else 2.0 for law in laws]))


def layer_strain(coord, conn, ezz, h, origin, height):
    """Mean element ezz of each layer of thickness h in a structured column, counted from z = origin."""
    distance = np.abs(coord[conn][:, :, -1].mean(axis=1) - origin)
    layer = np.floor(distance/h).astype(int)
    n = int(round(height/h))
    count = np.bincount(layer, minlength=n).astype(float)
    assert len(count) == n and count.min() > 0, 'the column mesh must be structured in layers of height h'
    return np.bincount(layer, weights=ezz, minlength=n)/count


def hydraulic_substeps(x, conn, storage, mobility, dt):
    """Hydraulic substeps of step dt: dt over the least lumped storage over absolute P1 diffusion row sum."""
    volume = simplex_volumes(x, conn)
    gradients = simplex_gradients(x, conn)
    operator = mobility*volume[:, None, None]*(gradients @ gradients.swapaxes(1, 2))
    mass, row_bound = lumped_volumes(x, conn, storage*volume), np.zeros(len(x))
    np.add.at(row_bound, conn, np.abs(operator).sum(axis=2))
    active = row_bound > 0
    return math.ceil(dt/float(np.min(mass[active]/row_bound[active])))
