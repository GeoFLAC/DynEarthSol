"""Shared configuration, solver invocation, output and mesh helpers for functional checks."""
import configparser
import hashlib
import itertools
import math
import os
from pathlib import Path
import subprocess

# Records of a save or checkpoint that differ between a fresh and a restarted run of the same physics.
RUNTIME_RECORDS = ('walltime_sec', 'provenance')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checkpoint_file(path):
    """The .vtkhdf file of a save or checkpoint path; None if it does not exist."""
    hdf = Path(str(path) + '.vtkhdf')
    return hdf if hdf.is_file() else None


def fields(path):
    """Payload bytes by field name (native) or dataset path plus '@<path>.layout' (HDF5 sibling)."""
    found = checkpoint_file(path)
    if found is not None and found != Path(path):
        import h5py
        import numpy as np
        result = {}
        with h5py.File(found, 'r') as stream:
            def read_dataset(name, obj):
                if isinstance(obj, h5py.Dataset):
                    value = obj[...]
                    if np.issubdtype(value.dtype, np.number):
                        assert np.isfinite(value).all(), (found, name)
                    result[name] = value.tobytes()
                    result['@' + name + '.layout'] = str((value.shape, value.dtype.str)).encode()
            stream.visititems(read_dataset)
        return result
    data = Path(path).read_bytes()
    entries = list(native_offsets(data).items())
    return {name: data[pos:entries[i + 1][1] if i + 1 < len(entries) else len(data)]
            for i, (name, pos) in enumerate(entries)}


def checkpoint_array(path, name, dtype='<f8'):
    """Flat array of the field name in a native or HDF5 save/checkpoint; dtype applies to native files."""
    import ast
    import numpy as np
    data = fields(path)
    key = next(k for k in data if not k.startswith('@') and k.split('/')[-1] == name)
    layout = data.get('@' + key + '.layout')
    return np.frombuffer(data[key], dtype=np.dtype(ast.literal_eval(layout.decode())[1] if layout else dtype))


def record_name(key):
    """Field name of a fields() key, without its HDF5 group path or layout marker."""
    name = key.rsplit('/', 1)[-1].lstrip('@')
    return name[:-len('.layout')] if name.endswith('.layout') else name


def runtime_key(key):
    """True for the wall-time and provenance records, which differ between fresh and restarted runs."""
    return record_name(key) in RUNTIME_RECORDS or '/provenance/' in key


def physical_differences(a, b, rtol=None):
    """Sorted non-runtime keys whose fields() payloads differ; rtol: record name -> float64 tolerance."""
    rtol = rtol or {}
    def same(key):
        if a.get(key) == b.get(key):
            return True
        if record_name(key) not in rtol or key not in a or key not in b:
            return False
        import numpy as np
        return np.allclose(np.frombuffer(a[key]), np.frombuffer(b[key]), rtol=rtol[record_name(key)], atol=0)
    return sorted(key for key in a.keys() | b.keys() if not runtime_key(key) and not same(key))


def compare_continuation(folder, original, resumed):
    """Assert resumed's later saves/checkpoints match original's physics; return the payload count."""
    payloads = 0
    for frame in resumed.frames[1:]:
        for kind in ('save', 'chkpt'):
            a = fields(Path(folder)/f'{original}.{kind}.{frame:06d}')
            b = fields(Path(folder)/f'resumed.{kind}.{frame:06d}')
            different = physical_differences(a, b)
            assert not different, different
            payloads += sum(not runtime_key(key) for key in a)
    return payloads


def simplex_volumes(x, conn):
    """Area (2D) or volume (3D) of each simplex conn over the node coordinates x."""
    import numpy as np
    cells = x[conn]
    return np.abs(np.linalg.det(cells[:, 1:]-cells[:, :1]))/math.factorial(conn.shape[1]-1)


def lumped_volumes(x, conn, volumes=None):
    """Nodal sum of an equal vertex share of each cell's volume (default simplex_volumes)."""
    import numpy as np
    if volumes is None:
        volumes = simplex_volumes(x, conn)
    lumped = np.zeros(len(x))
    np.add.at(lumped, conn, np.broadcast_to((volumes/conn.shape[1])[:, None], conn.shape))
    return lumped


def simplex_gradients(x, conn):
    """Constant P1 shape-function gradients of each simplex, shaped (cells, vertices, dimensions)."""
    import numpy as np
    vertices = np.concatenate((np.ones((*conn.shape, 1)), x[conn]), axis=2)
    return np.linalg.inv(vertices)[:, 1:, :].swapaxes(1, 2)


def facet_areas(points):
    """Length (2D) or area (3D) of facets given by their vertices, points shaped (..., dimension, dimension)."""
    import numpy as np
    edge = points[..., 1, :]-points[..., 0, :]
    if points.shape[-1] == 2:
        return np.linalg.norm(edge, axis=-1)
    return np.linalg.norm(np.cross(edge, points[..., 2, :]-points[..., 0, :]), axis=-1)/2


def simplex_heights(x, conn, volumes=None):
    """Smallest altitude of each simplex: dimension times its volume over its largest facet."""
    import numpy as np
    if volumes is None:
        volumes = simplex_volumes(x, conn)
    cells = x[conn]
    dim = conn.shape[1]-1
    facets = [facet_areas(cells[:, ids]) for ids in itertools.combinations(range(dim+1), dim)]
    return dim*volumes/np.max(facets, axis=0)


def roller_free(x):
    """Velocity components free of the rollers of a unit-width column: sides x_d = 0 and 1, base z = -1."""
    import numpy as np
    free = np.ones_like(x, dtype=bool)
    for d in range(x.shape[1]-1):
        free[:, d] = ~(np.isclose(x[:, d], 0) | np.isclose(x[:, d], 1))
    free[:, -1] = ~np.isclose(x[:, -1], -1)
    return free


def mesh_measures(model):
    """Coordinates, smallest element height and lumped nodal volumes of a model's first frame."""
    x = model.read_field(model.frames[0], 'coordinate')
    conn = model.read_field(model.frames[0], 'connectivity')
    volumes = simplex_volumes(x, conn)
    return x, float(simplex_heights(x, conn, volumes).min()), lumped_volumes(x, conn, volumes)


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
    """Run folder/name.cfg, written from cfg unless it is None, logging to name.log.

    threads=None keeps the caller's OMP_NUM_THREADS."""
    folder = Path(folder).resolve()
    folder.mkdir(parents=True, exist_ok=True)
    if cfg is not None:
        write_config(folder/f'{name}.cfg', cfg)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    if threads is not None:
        env['OMP_NUM_THREADS'] = str(threads)
    log = folder/f'{name}.log'
    with log.open('w') as stream:
        result = subprocess.run([str(Path(exe).resolve()), f'{name}.cfg'], cwd=folder,
                                env=env, stdout=stream, stderr=subprocess.STDOUT, timeout=timeout)
    result.stdout = log.read_text()
    if check:
        result.check_returncode()
    return result
