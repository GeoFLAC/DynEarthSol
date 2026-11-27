#!/usr/bin/env python3
"""Reject malformed tracer checkpoints before changing monitor histories."""
import argparse
import json
from pathlib import Path
import struct

from run_monitor_restart import (DEFAULT_CFG, clone_baseline, read_histories, render_cfg,
                                 run_config)
from des_test_utils import checkpoint_file, native_offsets


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
        fmt, offset = '<i' if integer else '<d', native_offsets(data)[key]
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
             ('missing-weights', 'monitor point weights', None, False, 'monitor other weights',
              'Incomplete monitor checkpoint state'),
             ('negative-count', 'monitor point count', -1, True, None,
              'Invalid monitor checkpoint count')]
    for name, key, value, integer, rename, error in cases:
        folder = output/name
        clone_baseline(baseline, folder)
        change(checkpoint_file(folder/'monitor_restart.chkpt.000001'), key, value, integer=integer, rename=rename)
        before = read_histories(folder)
        run_config(exe, restart_cfg, folder, 'restart', threads=1, expected_error=error)
        assert read_histories(folder) == before, name
    (output/'summary.json').write_text(json.dumps(dict(passed=True,
        rejected_checkpoints=len(cases), histories_unchanged_on_failure=True), indent=2)+'\n')
    print('PASS tracer state guards')


if __name__ == '__main__':
    main()
