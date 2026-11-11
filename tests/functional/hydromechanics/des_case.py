"""Shared configuration, execution and provenance for DES checks."""
import configparser
import copy
import hashlib
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[3]
SIDES = ('x0', 'x1', 'y0', 'y1', 'z0', 'z1')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_config(path, **options):
    cfg = configparser.ConfigParser(**options)
    cfg.optionxform = str
    cfg.read(path)
    return cfg


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


def execute(executable, folder, cfg, name, *, timeout=600):
    with (folder/f'{name}.cfg').open('w') as stream:
        cfg.write(stream)
    return run_solver([str(executable), f'{name}.cfg'], log=folder/f'{name}.log',
                      cwd=folder, check=True, timeout=timeout,
                      env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))


def run_model(executable, folder, cfg, name, *, timeout=600):
    from Dynearthsol import Dynearthsol
    cfg = copy.deepcopy(cfg)
    cfg['sim']['modelname'] = name
    execute(executable, folder, cfg, name, timeout=timeout)
    return Dynearthsol(str(folder/name))
