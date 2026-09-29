#!/usr/bin/env python3
"""Compare PT regular output schedules across a checkpoint with unchanged cadence."""
import argparse
import configparser
import json
from pathlib import Path
import struct
import subprocess

from pt_initial_checkpoint import fields


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('config', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--schedule', choices=('step', 'time', 'mixed', 'catchup'), required=True)
    args = parser.parse_args()
    executable = args.executable.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    cfg = configparser.ConfigParser(strict=False)
    cfg.optionxform = str
    cfg.read(args.config)
    assert cfg.getboolean('control', 'has_PT')
    dt = cfg.getfloat('control', 'fixed_dt')
    assert dt > 0, 'Use a positive fixed physical dt for this schedule test'
    # Use the solver YEAR2SEC constant, not a calendar-year approximation.
    year_seconds = 365.2422 * 86400
    interval = (0.2 if args.schedule == 'catchup' else 2.5) * dt / year_seconds
    cfg['sim'].update(max_steps='18', max_time_in_yr='1', has_initial_checkpoint='yes',
                      checkpoint_frame_interval='1', is_outputting_averaged_fields='no',
                      has_output_during_remeshing='no',
                      output_step_interval='3' if args.schedule == 'step' else
                      ('4' if args.schedule == 'mixed' else '2147483647'),
                      output_time_interval_in_yr='1' if args.schedule == 'step' else repr(interval))
    for phase in ('fresh', 'resumed'):
        cfg['sim'].update(modelname=phase, is_restarting=str(phase == 'resumed').lower())
        if phase == 'resumed':
            cfg['sim'].update(restarting_from_modelname=str(output / 'fresh'), restarting_from_frame='1')
        with (output / f'{phase}.cfg').open('w') as stream:
            cfg.write(stream)
        with (output / f'{phase}.log').open('w') as stream:
            subprocess.run([str(executable), f'{phase}.cfg'], cwd=output, stdout=stream,
                           stderr=subprocess.STDOUT, check=True, timeout=240)
    histories = {}
    for phase in ('fresh', 'resumed'):
        rows = []
        for path in sorted(output.glob(f'{phase}.save.*')):
            frame = int(path.name.split('.')[2])
            if frame < 1:
                continue
            data = fields(output / f'{phase}.save.{frame:06d}')
            rows.append((frame, struct.unpack('i', data['steps'])[0],
                         struct.unpack('d', data['time_sec'])[0]))
        histories[phase] = rows
    differences = {}
    if histories['fresh'] == histories['resumed']:
        for frame, _, _ in histories['fresh'][1:]:
            for suffix in ('save', 'chkpt'):
                a = fields(output / f'fresh.{suffix}.{frame:06d}')
                b = fields(output / f'resumed.{suffix}.{frame:06d}')
                differences[f'{suffix}.{frame}'] = [key for key in a.keys() | b.keys()
                    if key.split('/')[-1] != 'walltime_sec' and a.get(key) != b.get(key)]
    result = {'schedule': args.schedule, 'histories': histories, 'differences': differences}
    (output / 'comparison.json').write_text(json.dumps(result, indent=2) + '\n')
    assert histories['fresh'] == histories['resumed'], result
    assert len(histories['fresh']) > 1, 'No post-restart output compared'
    assert not any(differences.values()), result
    print(f'PASS: {args.schedule} schedule and all subsequent save/checkpoint payloads exact')


if __name__ == '__main__':
    main()
