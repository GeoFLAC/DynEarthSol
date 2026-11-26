"""Provenance and boundary helpers for the hydromechanics checks, on top of the
configuration and solver helpers of tests/functional/des_test_utils.py."""
import copy
import hashlib
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'tests/functional'))
from des_test_utils import read_config, run_des  # noqa: E402,F401  (read_config for callers)

SIDES = ('x0', 'x1', 'y0', 'y1', 'z0', 'z1')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reset_boundaries(bc, **values):
    for side in SIDES:
        for name, value in values.items():
            bc[f'{name}_{side}'] = value


def source_revision(root=ROOT):
    return subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip()


def source_diff(root=ROOT):
    return subprocess.check_output(['git', 'diff', 'HEAD', '--binary'], cwd=root)


def run_solver(command, *, log, **options):
    with Path(log).open('w') as stream:
        return subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, **options)


def run_model(executable, folder, cfg, name, *, timeout=600):
    from Dynearthsol import Dynearthsol
    cfg = copy.deepcopy(cfg)
    cfg['sim']['modelname'] = name
    run_des(executable, folder, name, cfg, timeout=timeout)
    return Dynearthsol(str(folder/name))
