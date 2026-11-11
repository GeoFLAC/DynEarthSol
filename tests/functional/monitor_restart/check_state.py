#!/usr/bin/env python3
"""Reject malformed tracer checkpoints before changing monitor histories."""
import argparse
import configparser
import csv
import io
import json
from pathlib import Path
import struct
import sys

from run_monitor_restart import (DEFAULT_CFG, REPO_ROOT, clone_baseline,
                                read_histories, render_cfg, run_config)

sys.path.insert(0, str(REPO_ROOT))
from Dynearthsol import Dynearthsol


def checkpoint(folder):
    path = folder/'monitor_restart.chkpt.000001'
    return path if path.is_file() else Path(str(path)+'.vtkhdf')


def change(path, key, value=None, *, integer=False, rename=None):
    if path.suffix == '.vtkhdf':
        import h5py
        import numpy as np
        with h5py.File(path, 'r+') as stream:
            if rename:
                stream.move(key, rename)
            else:
                dataset = stream[key]
                if dataset.is_virtual:
                    dataset = stream[dataset.virtual_sources()[0].dset_name]
                data = np.asarray(dataset[...]).copy()
                first = data.reshape(-1)
                first[0] = value(first[0]) if callable(value) else value
                dataset[...] = data
        return
    data = bytearray(path.read_bytes())
    if rename:
        old, new = key.encode()+b'\t', rename.encode()+b'\t'
        assert len(old) == len(new) and data[:4096].count(old) == 1
        data[:4096] = data[:4096].replace(old, new)
    else:
        entries = dict(line.rsplit(b'\t', 1) for line in
                       bytes(data[:4096]).split(b'\0', 1)[0].splitlines()[1:])
        fmt, offset = '<i' if integer else '<d', int(entries[key.encode()])
        if callable(value):
            value = value(struct.unpack_from(fmt, data, offset)[0])
        struct.pack_into(fmt, data, offset, value)
    path.write_bytes(data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    exe, output = args.executable.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    template = DEFAULT_CFG.read_text()
    baseline = output/'baseline'
    run_config(exe, render_cfg(template, restarting=False, write_header=True),
               baseline, 'fresh', threads=1)
    restart_cfg = render_cfg(template, restarting=True, write_header=True)
    cases = [('position-nan', 'monitor point positions', float('nan'), False, None,
              'Invalid monitor checkpoint positions'),
             ('element-range', 'monitor point elements', -2, True, None,
              'Invalid monitor checkpoint element'),
             ('weight-nan', 'monitor point weights', float('nan'), False, None,
              'Invalid monitor checkpoint weight'),
             ('weight-range', 'monitor point weights', 2., False, None,
              'Invalid monitor checkpoint weight'),
             ('weight-sum', 'monitor point weights', lambda w: w/2 if w > .5 else w+.1, False, None,
              'weights do not sum to one'),
             ('position-mismatch', 'monitor point positions', 1e6, False, None,
              'position and weights disagree'),
             ('missing-weights', 'monitor point weights', None, False, 'monitor other weights',
              'Incomplete monitor checkpoint state'),
             ('negative-count', 'monitor point count', -1, True, None,
              'Invalid monitor checkpoint count')]
    for name, key, value, integer, rename, error in cases:
        folder = output/name
        clone_baseline(baseline, folder)
        change(checkpoint(folder), key, value, integer=integer, rename=rename)
        before = read_histories(folder)
        run_config(exe, restart_cfg, folder, 'restart', threads=1, expected_error=error)
        assert read_histories(folder) == before, name
    # Position-only checkpoints retain their documented relocation fallback.
    legacy = output/'legacy-position-only'
    clone_baseline(baseline, legacy)
    change(checkpoint(legacy), 'monitor point elements', rename='monitor other elements')
    change(checkpoint(legacy), 'monitor point weights', rename='monitor other weights')
    before = read_histories(legacy)
    run_config(exe, restart_cfg, legacy, 'restart', threads=1)
    after = read_histories(legacy)
    for old, new in zip(before, after):
        a = list(csv.DictReader(io.StringIO(old.decode())))
        b = list(csv.DictReader(io.StringIO(new.decode())))
        assert len(a) == len(b) == 5 and a[:3] == b[:3]
        for row_a, row_b in zip(a, b):
            for key in row_a:
                if key.startswith(('query_', 'coord_')):
                    assert abs(float(row_a[key])-float(row_b[key])) < 1e-12
    # The locator accepts tiny negative weights at a boundary. Their balancing
    # corner weight can exceed one by several times that tolerance.
    model = Dynearthsol(str(baseline/'monitor_restart'))
    corner_checked = False
    if model.ndims == 3:
        import numpy as np
        x, cells = (model.read_field(0, key) for key in ('coordinate', 'connectivity'))
        node = int(np.argmin(np.sum(x*x, axis=1)))
        e = int(np.flatnonzero((cells == node).any(axis=1))[0])
        weights = np.full(4, -4.6e-11)
        weights[np.flatnonzero(cells[e] == node)[0]] = 1+3*4.6e-11
        q = weights@x[cells[e]]
        cfg = configparser.ConfigParser(interpolation=None)
        cfg.optionxform = str
        cfg.read_string(template)
        for key, value in zip(('points_x', 'points_y', 'points_z'), q):
            cfg['monitor'][key] = repr([float(value), .5 if key != 'points_z' else -.5])
        stream = io.StringIO()
        cfg.write(stream)
        near_corner = output/'near-corner'
        run_config(exe, render_cfg(stream.getvalue(), restarting=False, write_header=True),
                   near_corner, 'fresh', threads=1)
        path = checkpoint(near_corner)
        if path.suffix == '.vtkhdf':
            import h5py
            with h5py.File(path) as source:
                saved = source['monitor point weights'][:4]
        else:
            data = path.read_bytes()
            entries = dict(line.rsplit(b'\t', 1) for line in
                           data[:4096].split(b'\0', 1)[0].splitlines()[1:])
            saved = struct.unpack_from('<4d', data, int(entries[b'monitor point weights']))
        assert max(saved) > 1+1e-10, saved
        expected = read_histories(near_corner)
        run_config(exe, render_cfg(stream.getvalue(), restarting=True, write_header=True),
                   near_corner, 'restart', threads=1)
        assert read_histories(near_corner) == expected
        corner_checked = True
    (output/'summary.json').write_text(json.dumps(dict(passed=True,
        rejected_checkpoints=len(cases), histories_unchanged_on_failure=True,
        legacy_position_only_supported=True, legacy_csv_byte_exact=before == after,
        boundary_tolerance_checked=corner_checked), indent=2)+'\n')
    print('PASS tracer state guards and legacy fallback')


if __name__ == '__main__':
    main()
