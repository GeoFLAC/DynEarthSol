#!/usr/bin/env python3
"""Exercise selected reference resets, inactive behavior and exact continuation."""
import argparse
import configparser
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT), str(ROOT/'tests/functional')]
from Dynearthsol import Dynearthsol, scan_frames
from pt_initial_checkpoint import fields


def configuration(bits):
    cfg = configparser.ConfigParser(); cfg.optionxform = str
    cfg.read(ROOT/'tests/functional/displacement_output/displacement_output.cfg')
    cfg['sim'].update(modelname='reset', max_steps='1', has_output_during_remeshing='true',
                      has_marker_output='true')
    cfg['mesh'].update(meshing_elem_shape='0', resolution='.5', quality_check_step_interval='1',
        min_quality='.01', max_boundary_distortion='2e-7', remeshing_option='14',
        num_remeshing_reset_boundaries=str(len(bits)), remeshing_reset_boundaries=json.dumps(bits))
    cfg['control']['fixed_dt'] = '.001'
    cfg['bc'].update(vbc_val_x0='.001', vbc_val_x1='.001', vbc_val_z0='-.001', vbc_val_z1='-.001')
    return cfg


def run(exe, folder, cfg, name='input', expected=0):
    folder.mkdir(parents=True, exist_ok=True)
    with (folder/(name+'.cfg')).open('w') as stream: cfg.write(stream)
    result = subprocess.run([str(exe), name+'.cfg'], cwd=folder, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90,
        env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'))
    (folder/(name+'.log')).write_text(result.stdout)
    assert result.returncode == expected, (expected, result.returncode, result.stdout[-2000:])


def inspect(folder, bits):
    model = Dynearthsol(str(folder/'reset'))
    rows = scan_frames(model.modelname)
    assert len(rows) >= 4, 'No before/after remesh frames'
    final = rows[-1]['frame']; x = model.read_field(final, 'coordinate')
    flags = model.read_field(final, 'bcflag'); ref = model.read_field(final, 'coord0')
    for bit, dimension, initial, velocity in ((1,0,0.,.001),(2,0,1.,.001),
                                            (16,1,-1.,-.001),(32,1,0.,-.001)):
        nodes = (flags & bit) != 0
        expected = initial if bit in bits else initial+.001*velocity
        np.testing.assert_allclose(x[nodes, dimension], expected, rtol=0, atol=1e-12)
        if bit in bits:
            np.testing.assert_allclose(ref[nodes, dimension], x[nodes, dimension], rtol=0, atol=1e-12)
    return dict(frames=len(rows), nodes=len(x), selected=bits)


def check_continuation(exe, out, cfg):
    out.mkdir(parents=True, exist_ok=False)
    cfg['sim']['max_steps'] = '3'
    fresh = out/'fresh'; run(exe,fresh,cfg)
    seed = Dynearthsol(str(fresh/'reset')); rows = scan_frames(seed.modelname)
    def exists(path):
        return path.is_file() or Path(str(path)+'.vtkhdf').is_file()
    checkpoint = next(r['frame'] for r in rows if r['steps']==1 and
                      exists(fresh/f"reset.chkpt.{r['frame']:06d}"))
    cfg['sim'].update(modelname='continued',is_restarting='true',
        restarting_from_modelname=str(fresh/'reset'),restarting_from_frame=str(checkpoint))
    resumed = out/'resumed'; run(exe,resumed,cfg)
    child = Dynearthsol(str(resumed/'continued')); matched = 0
    for row in scan_frames(child.modelname):
        frame=row['frame']
        assert any(r['frame']==frame and r['time']==row['time'] and r['steps']==row['steps'] for r in rows)
        for kind in ('save','chkpt'):
            original=fresh/f'reset.{kind}.{frame:06d}'
            current=resumed/f'continued.{kind}.{frame:06d}'
            if kind == 'chkpt' and frame == checkpoint:
                continue  # Restart opens the seed but writes no duplicate checkpoint.
            assert exists(original)==exists(current),(kind,frame)
            if not exists(original): continue
            a=fields(original)
            b=fields(current)
            assert a.keys() == b.keys(), (kind, frame, a.keys() ^ b.keys())
            for key in a:
                if key.split('/')[-1] not in ('walltime_sec','provenance'):
                    assert a[key]==b[key],(kind,frame,key);matched+=1
    return matched


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('exe', type=Path); parser.add_argument('output', type=Path)
    parser.add_argument('--exe-3d', type=Path)
    args = parser.parse_args(); exe = args.exe.resolve(); out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    results = []
    for bits in ([1], [1,16], [1,2,16,32]):
        folder = out/('selected-'+'-'.join(map(str,bits)))
        run(exe, folder, configuration(bits)); results.append(inspect(folder,bits))
    invalid = [(-1,'[]'), (1,'[0]'), (1,'[3]'), (1,'[1024]'), (0,'[]'), (2,'[1]')]
    for i,(count,bits) in enumerate(invalid):
        cfg = configuration([1]); cfg['mesh'].update(num_remeshing_reset_boundaries=str(count),remeshing_reset_boundaries=bits)
        run(exe, out/f'invalid-{i}', cfg, expected=11)
    cfg = configuration([1]); cfg['mesh']['meshing_elem_shape']='1'
    run(exe,out/'structured',cfg,expected=11)
    if args.exe_3d:
        run(args.exe_3d.resolve(),out/'3d',configuration([1]),expected=30)
    matched = 0
    for selected in (False, True):
        for method in ('explicit', 'gvs'):
            cfg = configuration([1])
            cfg['control']['use_global_velocity_scaling'] = str(method == 'gvs')
            if not selected:
                cfg['mesh']['remeshing_option'] = '0'
                cfg['mesh']['min_quality'] = '.99'
                # This regression also runs with the parent executable.
                del cfg['mesh']['num_remeshing_reset_boundaries']
                del cfg['mesh']['remeshing_reset_boundaries']
            matched += check_continuation(exe, out/f'restart-{method}-{selected}', cfg)
    (out/'summary.json').write_text(json.dumps(dict(passed=True,cases=results,
        input_rejections=len(invalid)+1+bool(args.exe_3d),exact_restart_payloads=matched),indent=2)+'\n')
    print('PASS selected boundary resets and continuation;',matched,'exact payloads')


if __name__=='__main__': main()
