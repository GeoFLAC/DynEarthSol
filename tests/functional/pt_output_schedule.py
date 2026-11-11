#!/usr/bin/env python3
"""Compare regular output schedules and physical state across a checkpoint."""
import argparse
import configparser
import json
import math
from pathlib import Path
import struct
import subprocess

from pt_initial_checkpoint import fields


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('config', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--schedule', choices=('step', 'step-default', 'time', 'mixed', 'catchup', 'overflow'), required=True)
    parser.add_argument('--method', choices=('pt', 'explicit', 'gvs'), default='pt')
    parser.add_argument('--rsf', action='store_true',
                        help='Retain the RSF fixture (forces GVS); use with pt or gvs')
    parser.add_argument('--step-after-restart', type=int,
                        help='Change the step cadence at restart; requires --schedule step')
    parser.add_argument('--time-after-restart', type=float,
                        help='New time cadence in fixed-dt units (>1); requires --schedule time')
    args = parser.parse_args()
    executable = args.executable.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    cfg = configparser.ConfigParser(strict=False)
    cfg.optionxform = str
    cfg.read(args.config)
    cfg['control'].update(has_PT=str(args.method == 'pt').lower(),
                          use_global_velocity_scaling=str(args.method == 'gvs').lower())
    if args.rsf:
        assert args.method != 'explicit', 'RSF forces GVS and cannot test ordinary explicit scheduling'
    else:
        cfg['mat']['rheology_type'] = 'elastic'
    uses_gvs = args.method == 'gvs' or args.rsf
    if args.step_after_restart is not None:
        assert args.schedule == 'step' and args.step_after_restart > 0 and args.time_after_restart is None
    if args.time_after_restart is not None:
        assert args.schedule == 'time' and args.time_after_restart > 1
    dt = cfg.getfloat('control', 'fixed_dt')
    assert dt > 0, 'Use a positive fixed physical dt for this schedule test'
    # Use the solver YEAR2SEC constant, not a calendar-year approximation.
    year_seconds = 365.2422 * 86400
    interval = (0.2 if args.schedule == 'catchup' else 2.5) * dt / year_seconds
    cfg['sim'].update(max_steps='64' if args.schedule == 'overflow' else '18', max_time_in_yr='1', has_initial_checkpoint='yes',
                      checkpoint_frame_interval='1', is_outputting_averaged_fields='no',
                      has_output_during_remeshing='no',
                      output_step_interval='3' if args.schedule in ('step', 'step-default') else
                      ('4' if args.schedule == 'mixed' else '100000000' if args.schedule == 'overflow' else '2147483647'),
                      output_time_interval_in_yr='1' if args.schedule == 'step' else repr(interval))
    if args.schedule == 'step-default':
        cfg['sim'].pop('output_time_interval_in_yr', None)
    changed = args.step_after_restart is not None or args.time_after_restart is not None
    phases = ('fresh', 'reference', 'resumed') if changed or args.schedule == 'overflow' else ('fresh', 'resumed')
    original_step_interval = cfg['sim']['output_step_interval']
    original_time_interval = cfg['sim'].get('output_time_interval_in_yr')
    for phase in phases:
        cfg['sim'].update(modelname=phase, is_restarting=str(phase == 'resumed').lower())
        cfg['sim']['output_step_interval'] = original_step_interval
        if original_time_interval is not None:
            cfg['sim']['output_time_interval_in_yr'] = original_time_interval
        if phase == 'reference':
            cfg['sim'].update(output_step_interval='1', output_time_interval_in_yr='1')
        if phase == 'resumed':
            cfg['sim'].update(restarting_from_modelname=str(output / 'fresh'), restarting_from_frame='1')
            if args.step_after_restart:
                cfg['sim']['output_step_interval'] = str(args.step_after_restart)
            if args.time_after_restart:
                cfg['sim']['output_time_interval_in_yr'] = repr(args.time_after_restart*dt/year_seconds)
        with (output / f'{phase}.cfg').open('w') as stream:
            cfg.write(stream)
        with (output / f'{phase}.log').open('w') as stream:
            subprocess.run([str(executable), f'{phase}.cfg'], cwd=output, stdout=stream,
                           stderr=subprocess.STDOUT, check=True, timeout=240)
    histories = {}
    for phase in phases:
        rows = []
        for path in sorted(output.glob(f'{phase}.save.*')):
            frame = int(path.name.split('.')[2])
            if frame < 1:
                continue
            data = fields(output / f'{phase}.save.{frame:06d}')
            # The resumed run rewrites its initial save, not its source checkpoint.
            if phase == 'fresh' or frame != 1:
                checkpoint = fields(output / f'{phase}.chkpt.{frame:06d}')
                assert math.isfinite(struct.unpack('d', checkpoint['reference_frame_time'])[0]), 'Nonfinite output clock'
            rows.append((frame, struct.unpack('i', data['steps'])[0],
                         struct.unpack('d', data['time_sec'])[0]))
        histories[phase] = rows
    if args.schedule == 'overflow':
        expected, index = [], 1
        for _, step, time in histories['reference']:
            target = index*interval*year_seconds
            if (time >= target if uses_gvs else time > target):
                expected.append(step)
                index += 1
        assert [r[1] for r in histories['fresh']] == expected, ('Overflow changed time cadence', histories, expected)
    differences = {}
    if changed:
        start = histories['fresh'][0][1]
        if args.step_after_restart:
            expected = list(range(start, 19, args.step_after_restart))
        else:
            expected = [start]; index = 1
            epoch = histories['fresh'][0][2]
            period = float(cfg['sim']['output_time_interval_in_yr'])
            for _, step, time in histories['reference']:
                if step <= start: continue
                elapsed, target = time-epoch, index*period*year_seconds
                due = elapsed >= target if uses_gvs else elapsed > target
                if due:
                    expected.append(step); index += 1
        assert [row[1] for row in histories['resumed']] == expected
        reference = {step: (frame, time) for frame, step, time in histories['reference']}
        clock = {'reference_frame_time', 'output start time', 'output start step', 'next regular frame',
                 'PT output start time', 'PT output start step', 'PT next regular frame',
                 'output step interval', 'output time interval'}
        for frame, step, time in histories['resumed'][1:]:
            other, expected_time = reference[step]
            assert time == expected_time
            for suffix in ('save', 'chkpt'):
                a = fields(output/f'reference.{suffix}.{other:06d}')
                b = fields(output/f'resumed.{suffix}.{frame:06d}')
                differences[f'{suffix}.{frame}'] = [key for key in a.keys() | b.keys()
                    if key.split('/')[-1] not in clock | {'walltime_sec', 'provenance'}
                    and a.get(key) != b.get(key)]
        result = dict(method=args.method, schedule=args.schedule, step_after_restart=args.step_after_restart,
                      time_after_restart=args.time_after_restart,
                      histories=histories, differences=differences)
        (output/'comparison.json').write_text(json.dumps(result, indent=2)+'\n')
        assert differences and not any(differences.values()), result
        print('PASS: changed cadence starts at restart; physical fields match the independent output-every-step run')
        return
    if histories['fresh'] == histories['resumed']:
        for frame, _, _ in histories['fresh'][1:]:
            for suffix in ('save', 'chkpt'):
                a = fields(output / f'fresh.{suffix}.{frame:06d}')
                b = fields(output / f'resumed.{suffix}.{frame:06d}')
                differences[f'{suffix}.{frame}'] = [key for key in a.keys() | b.keys()
                    if key.split('/')[-1] not in ('walltime_sec', 'provenance') and a.get(key) != b.get(key)]
    result = {'method': args.method, 'schedule': args.schedule, 'histories': histories, 'differences': differences}
    (output / 'comparison.json').write_text(json.dumps(result, indent=2) + '\n')
    assert histories['fresh'] == histories['resumed'], result
    assert len(histories['fresh']) > 1, 'No post-restart output compared'
    assert not any(differences.values()), result
    print(f'PASS: {args.schedule} schedule and all subsequent save/checkpoint payloads exact')


if __name__ == '__main__':
    main()
