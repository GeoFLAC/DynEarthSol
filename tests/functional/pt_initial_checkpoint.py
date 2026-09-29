#!/usr/bin/env python3
"""Check PT frame-zero state and exact restart using a short supplied PT config.

Usage: pt_initial_checkpoint.py EXECUTABLE CONFIG OUTPUT_DIRECTORY
The output directory retains configs, logs and results for independent review.
"""
import configparser
import json
from pathlib import Path
import struct
import subprocess
import sys


def fields(path):
    data = path.read_bytes()
    entries = []
    for line in data[:4096].split(b'\0', 1)[0].splitlines()[1:]:
        name, pos = line.split(b'\t')
        entries.append((name.decode(), int(pos)))
    return {name: data[pos:entries[i + 1][1] if i + 1 < len(entries) else len(data)]
            for i, (name, pos) in enumerate(entries)}


def main():
    executable, source, output = map(Path, sys.argv[1:])
    executable, output = executable.resolve(), output.resolve()
    if not executable.is_file() or not source.is_file():
        raise FileNotFoundError('Executable and source config must exist')
    output.mkdir(parents=True, exist_ok=False)
    cfg = configparser.ConfigParser(strict=False)
    cfg.optionxform = str
    cfg.read(source)
    assert cfg.getboolean('control', 'has_PT')
    cfg['sim'].update(modelname='fresh', max_steps='3', output_step_interval='1',
                      checkpoint_frame_interval='1', has_initial_checkpoint='yes',
                      is_restarting='no')
    for phase in ('fresh', 'restart'):
        if phase == 'restart':
            cfg['sim'].update(modelname='resumed', is_restarting='yes',
                              restarting_from_modelname=str(output / 'fresh'),
                              restarting_from_frame='0')
        with (output / (phase + '.cfg')).open('w') as stream:
            cfg.write(stream)
        with (output / (phase + '.log')).open('w') as stream:
            subprocess.run([str(executable), phase + '.cfg'], cwd=output,
                           stdout=stream, stderr=subprocess.STDOUT, check=True, timeout=240)
        if phase == 'fresh':
            initial = fields(output / 'fresh.chkpt.000000')
            assert struct.unpack('i', initial['initial equilibrium done'])[0] == 1
            saved = fields(output / 'fresh.save.000000')
            velocity = struct.unpack('d' * (len(saved['velocity']) // 8), saved['velocity'])
            if all(value == 0 for value in velocity):
                rates = struct.unpack('d' * (len(saved['strain-rate']) // 8), saved['strain-rate'])
                assert all(value == 0 for value in rates), 'Numerical correction rate leaked into initial output'
    differences = {}
    for suffix in ('save', 'chkpt'):
        a = fields(output / ('fresh.' + suffix + '.000003'))
        b = fields(output / ('resumed.' + suffix + '.000003'))
        assert a.keys() == b.keys()
        # Wall-clock duration is observational metadata, not physical state.
        differences[suffix] = [key for key in a if key != 'walltime_sec' and a[key] != b[key]]
    (output / 'comparison.json').write_text(json.dumps(differences, indent=2) + '\n')
    assert not any(differences.values()), differences
    print('PASS: equilibrated frame-zero checkpoint; all final physical save/checkpoint fields exact')


if __name__ == '__main__':
    main()
