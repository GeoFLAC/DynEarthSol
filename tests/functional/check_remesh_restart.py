#!/usr/bin/env python3
"""Check exact explicit/GVS continuation across a pending remesh stage."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT/'tests/functional')]
from Dynearthsol import Dynearthsol, scan_frames
from pt_initial_checkpoint import fields
from des_test_utils import checkpoint_file, read_config, run_des


def configuration():
    cfg = read_config(ROOT/'tests/functional/displacement_output/displacement_output.cfg')
    cfg['sim'].update(modelname='reset', max_steps='1', has_output_during_remeshing='true',
                      has_marker_output='true')
    cfg['mesh'].update(meshing_elem_shape='0', resolution='.5', quality_check_step_interval='1',
        min_quality='.99', max_boundary_distortion='2e-7', remeshing_option='0')
    cfg['control']['fixed_dt'] = '.001'
    cfg['bc'].update(vbc_val_x0='.001', vbc_val_x1='.001', vbc_val_z0='-.001', vbc_val_z1='-.001')
    return cfg


def check_continuation(exe, out, cfg):
    out.mkdir(parents=True, exist_ok=False)
    cfg['sim']['max_steps'] = '3'
    fresh = out/'fresh'; run_des(exe, fresh, 'input', cfg, timeout=90)
    seed = Dynearthsol(str(fresh/'reset')); rows = scan_frames(seed.modelname)
    assert any(a['steps'] == b['steps'] and a['time'] == b['time']
               for a, b in zip(rows, rows[1:])), 'No before/after remesh frames'
    def exists(path):
        return checkpoint_file(path) is not None
    checkpoint = next(r['frame'] for r in rows if r['steps']==1 and
                      exists(fresh/f"reset.chkpt.{r['frame']:06d}"))
    cfg['sim'].update(modelname='continued',is_restarting='true',
        restarting_from_modelname=str(fresh/'reset'),restarting_from_frame=str(checkpoint))
    resumed = out/'resumed'; run_des(exe, resumed, 'input', cfg, timeout=90)
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
    args = parser.parse_args(); exe = args.exe.resolve(); out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    matched = 0
    for method in ('explicit', 'gvs'):
        cfg = configuration()
        cfg['control']['use_global_velocity_scaling'] = str(method == 'gvs')
        matched += check_continuation(exe, out/f'restart-{method}', cfg)
    (out/'summary.json').write_text(json.dumps(dict(passed=True,
        exact_restart_payloads=matched),indent=2)+'\n')
    print('PASS explicit/GVS remesh continuation;',matched,'exact payloads')


if __name__=='__main__': main()
