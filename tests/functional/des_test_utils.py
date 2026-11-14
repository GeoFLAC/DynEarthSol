"""Shared configuration and solver invocation for functional checks."""
import configparser
import os
from pathlib import Path
import subprocess


def checkpoint_file(path):
    """The native file at path or its .vtkhdf sibling, whichever exists; None if neither."""
    path = Path(path)
    hdf = Path(str(path) + '.vtkhdf')
    return path if path.is_file() else hdf if hdf.is_file() else None


def native_offsets(data):
    """Field name -> byte offset, in file order, from a native save or checkpoint header."""
    header = bytes(data[:4096]).split(b'\0', 1)[0]
    return {name.decode(): int(pos) for name, pos in
            (line.rsplit(b'\t', 1) for line in header.splitlines()[1:])}


def read_config(source=None, **kwargs):
    cfg = configparser.ConfigParser(**kwargs)
    cfg.optionxform = str
    if source is not None:
        with Path(source).open() as stream:
            cfg.read_file(stream)
    return cfg


def write_config(path, cfg):
    with Path(path).open('w') as stream:
        if isinstance(cfg, str):
            stream.write(cfg)
        else:
            cfg.write(stream)


def run_des(exe, folder, name, cfg, timeout=120, *, threads=1, check=True):
    folder = Path(folder).resolve()
    folder.mkdir(parents=True, exist_ok=True)
    write_config(folder/f'{name}.cfg', cfg)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads),
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    log = folder/f'{name}.log'
    with log.open('w') as stream:
        result = subprocess.run([str(Path(exe).resolve()), f'{name}.cfg'], cwd=folder,
                                env=env, stdout=stream, stderr=subprocess.STDOUT,
                                text=True, timeout=timeout)
    result.stdout = log.read_text()
    if check:
        result.check_returncode()
    return result
