#!/usr/bin/env python3
"""Reject malformed tracer checkpoints before changing monitor histories."""
import argparse
import json
from pathlib import Path

from run_monitor_restart import (DEFAULT_CFG, clone_baseline, read_histories, render_cfg,
                                 run_config)
from des_test_utils import edit_checkpoint


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
    cases = [('position-nan', 'monitor point positions', float('nan'), None,
              'Invalid monitor checkpoint positions'),
             ('element-range', 'monitor point elements', -2, None,
              'Invalid monitor checkpoint element'),
             ('weight-nan', 'monitor point weights', float('nan'), None,
              'Invalid monitor checkpoint weight'),
             ('missing-weights', 'monitor point weights', None, 'monitor other weights',
              'Incomplete monitor checkpoint state'),
             ('negative-count', 'monitor point count', -1, None,
              'Invalid monitor checkpoint count')]
    for name, key, value, rename, error in cases:
        folder = output/name
        clone_baseline(baseline, folder)
        checkpoint = folder/'monitor_restart.chkpt.000001'
        if rename:
            edit_checkpoint(checkpoint, rename={key: rename})
        else:
            def first(data, value=value):
                data.reshape(-1)[0] = value
                return data
            edit_checkpoint(checkpoint, replace={key: first})
        before = read_histories(folder)
        run_config(exe, restart_cfg, folder, 'restart', threads=1, expected_error=error)
        assert read_histories(folder) == before, name
    (output/'summary.json').write_text(json.dumps(dict(passed=True,
        rejected_checkpoints=len(cases), histories_unchanged_on_failure=True), indent=2)+'\n')
    print('PASS tracer state guards')


if __name__ == '__main__':
    main()
