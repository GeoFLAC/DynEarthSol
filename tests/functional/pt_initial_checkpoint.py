#!/usr/bin/env python3
"""Check PT initialization and exact binary/HDF5 restart with a short PT config.

Usage: pt_initial_checkpoint.py EXECUTABLE CONFIG OUTPUT_DIRECTORY [--restart-frame N] [--steps N]
The output directory retains configs, logs and results for independent review.
Use --require-remesh with a config that triggers mesh quality control to check
continuation across an actual mesh change, not just checkpoint serialization.
"""
import argparse
import configparser
import json
from pathlib import Path
import struct
import subprocess


def fields(path):
    hdf_path = Path(str(path) + '.vtkhdf')
    if not path.is_file() and hdf_path.is_file():
        import h5py
        import numpy as np
        result = {}
        with h5py.File(hdf_path, 'r') as stream:
            def read_dataset(name, obj):
                if isinstance(obj, h5py.Dataset):
                    value = obj[...]
                    if np.issubdtype(value.dtype, np.number):
                        assert np.isfinite(value).all(), (hdf_path, name)
                    result[name] = value.tobytes()
                    result['@' + name + '.layout'] = str((value.shape, value.dtype.str)).encode()
            stream.visititems(read_dataset)
        return result
    data = path.read_bytes()
    entries = []
    for line in data[:4096].split(b'\0', 1)[0].splitlines()[1:]:
        name, pos = line.split(b'\t')
        entries.append((name.decode(), int(pos)))
    return {name: data[pos:entries[i + 1][1] if i + 1 < len(entries) else len(data)]
            for i, (name, pos) in enumerate(entries)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--restart-frame', type=int, default=0)
    parser.add_argument('--steps', type=int, default=3)
    parser.add_argument('--timeout', type=float, default=240,
                        help='Maximum seconds for each fresh/restart process (default: 240)')
    parser.add_argument('--require-remesh', action='store_true')
    parser.add_argument('--require-stationary', action='store_true',
                        help='Require a no-load hydrostatic fixture to retain stress and pressure')
    args = parser.parse_args()
    if not args.timeout > 0:
        parser.error('--timeout must be positive')
    if not 0 <= args.restart_frame < args.steps:
        parser.error('require 0 <= restart-frame < steps')
    executable, source, output = args.executable, args.source, args.output
    executable, output = executable.resolve(), output.resolve()
    if not executable.is_file() or not source.is_file():
        raise FileNotFoundError('Executable and source config must exist')
    output.mkdir(parents=True, exist_ok=False)
    cfg = configparser.ConfigParser(strict=False)
    cfg.optionxform = str
    cfg.read(source)
    assert cfg.getboolean('control', 'has_PT')
    if cfg.getint('mesh', 'meshing_option', fallback=1) in (90, 91):
        poly = Path(cfg['mesh']['poly_filename'])
        if not poly.is_absolute():
            cfg['mesh']['poly_filename'] = str((source.resolve().parent / poly).resolve())
    cfg['sim'].update(modelname='fresh', max_steps=str(args.steps), output_step_interval='1',
                      checkpoint_frame_interval='1', has_initial_checkpoint='yes',
                      has_output_during_remeshing='no',
                      is_restarting='no')
    for phase in ('fresh', 'restart'):
        if phase == 'restart':
            cfg['sim'].update(modelname='resumed', is_restarting='yes',
                              restarting_from_modelname=str(output / 'fresh'),
                              restarting_from_frame=str(args.restart_frame))
        with (output / (phase + '.cfg')).open('w') as stream:
            cfg.write(stream)
        with (output / (phase + '.log')).open('w') as stream:
            subprocess.run([str(executable), phase + '.cfg'], cwd=output,
                           stdout=stream, stderr=subprocess.STDOUT, check=True, timeout=args.timeout)
        if phase == 'fresh':
            initial = fields(output / 'fresh.chkpt.000000')
            assert struct.unpack('i', initial['initial equilibrium done'])[0] == 1
            saved = fields(output / 'fresh.save.000000')
            velocity = struct.unpack('d' * (len(saved['velocity']) // 8), saved['velocity'])
            if all(value == 0 for value in velocity):
                rates = struct.unpack('d' * (len(saved['strain-rate']) // 8), saved['strain-rate'])
                assert all(value == 0 for value in rates), 'Numerical correction rate leaked into initial output'
    # Check the restored velocity before another solve can hide a boundary-map
    # error. Then compare every continuation frame, including time and saved dt.
    a = fields(output / f'fresh.save.{args.restart_frame:06d}')
    b = fields(output / f'resumed.save.{args.restart_frame:06d}')
    differences = {'restored_velocity': [] if a['velocity'] == b['velocity'] else ['velocity']}
    for frame in range(args.restart_frame + 1, args.steps + 1):
        for suffix in ('save', 'chkpt'):
            a = fields(output / f'fresh.{suffix}.{frame:06d}')
            b = fields(output / f'resumed.{suffix}.{frame:06d}')
            assert a.keys() == b.keys()
            # Wall-clock duration is observational metadata, not physical state.
            differences[f'{suffix}.{frame}'] = [key for key in a if key.split('/')[-1] != 'walltime_sec' and a[key] != b[key]]
    (output / 'comparison.json').write_text(json.dumps(differences, indent=2) + '\n')
    assert not any(differences.values()), differences
    if args.require_stationary:
        import numpy as np
        initial = fields(output / 'fresh.save.000000')
        changes = {}
        for frame in range(1, args.steps + 1):
            current = fields(output / f'fresh.save.{frame:06d}')
            for name in ('stress', 'pore pressure'):
                a = np.frombuffer(initial[name], dtype=np.float64)
                b = np.frombuffer(current[name], dtype=np.float64)
                scale = max(1.0, float(np.max(np.abs(a))))
                changes[f'{name}.{frame}'] = float(np.max(np.abs(b - a))) / scale
            assert initial['coordinate'] == current['coordinate'], 'Stationary mesh moved'
        (output / 'stationarity.json').write_text(json.dumps(changes, indent=2) + '\n')
        assert all(value <= 1e-8 for value in changes.values()), changes
    if args.require_remesh:
        for phase in ('fresh', 'restart'):
            assert 'Remeshing finished.' in (output / (phase + '.log')).read_text(), phase
        first = fields(output / 'fresh.save.000000')
        last = fields(output / f'fresh.save.{args.steps:06d}')
        assert first['connectivity'] != last['connectivity'], 'Mesh connectivity did not change'
    print(f'PASS: equilibrated initial state; restored velocity and all continuation fields exact from frame {args.restart_frame}')


if __name__ == '__main__':
    main()
