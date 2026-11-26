#!/usr/bin/env python
'''Check that the current tools and DES still use the output of an earlier build.

usage: python tests/backward-compat/check_backward_compat.py OLD_DIR EXE [-w WORK]

OLD_DIR holds the output of compat-2d.cfg (model "compat") from an earlier build, in
either format; EXE is the current dynearthsol2d. The check
  1. converts every frame to .vtu and .vtp: des-binary frames with 2vtk-legacy.py,
     .vtkhdf frames with 2vtk.py, which reads them through Dynearthsol.py;
  2. converts the last checkpoint with utils/convert_checkpoint.py;
  3. restarts EXE from it for one more frame, whose fields must be finite;
  4. converts that frame with 2vtk.py.
It stops at the first step that fails. WORK (default: a new temporary directory)
receives the converted files and the restarted run.
'''
import argparse
import glob
import os
import re
import subprocess
import sys
import tempfile

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
MODEL = 'compat'
CFG = os.path.join(HERE, 'compat-2d.cfg')
TIMEOUT = 900       # seconds per command; the case runs in seconds


def fail(msg):
    sys.exit('FAIL: ' + msg)


def run(cmd, cwd, what):
    r = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       universal_newlines=True, timeout=TIMEOUT)
    if r.returncode != 0:
        fail('%s exited with %d:\n%s' % (what, r.returncode, r.stdout[-3000:]))
    return r.stdout


def frames(old, kind):
    '''The frame numbers of OLD_DIR's save or chkpt files, and whether they are .vtkhdf.'''
    pat = re.compile(r'%s\.%s\.(\d{6})(\.vtkhdf)?$' % (MODEL, kind))
    found = [pat.match(os.path.basename(f)) for f in glob.glob(os.path.join(old, '%s.%s.*' % (MODEL, kind)))]
    found = [m for m in found if m]
    if not found:
        fail('no %s.%s.* files in %s' % (MODEL, kind, old))
    if len({bool(m.group(2)) for m in found}) != 1:
        fail('%s holds both des-binary and .vtkhdf %s files' % (old, kind))
    return sorted(int(m.group(1)) for m in found), bool(found[0].group(2))


def convert(tool, model, outdir, numbers):
    '''Run TOOL (2vtk.py or 2vtk-legacy.py) on MODEL into OUTDIR; each frame of NUMBERS must
    give a .vtu and a .vtp file.'''
    os.makedirs(outdir, exist_ok=True)
    run([sys.executable, os.path.join(ROOT, tool), '-c', '-m', model], outdir, tool)
    name = os.path.basename(model)
    for f in numbers:
        for ext in ('vtu', 'vtp'):
            out = os.path.join(outdir, '%s.%06d.%s' % (name, f, ext))
            if not os.path.isfile(out):
                fail('%s wrote no %s' % (tool, out))


def info_steps(path, frame):
    '''The step of FRAME in a .info file (columns as Output::write_info() in output.cxx).'''
    for line in open(path):
        cols = line.split()
        if int(cols[0]) == frame:
            return int(cols[1])
    fail('frame %d is not in %s' % (frame, path))


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('old', help='directory with the earlier build\'s output of compat-2d.cfg')
    p.add_argument('exe', help='the current dynearthsol2d')
    p.add_argument('-w', dest='work', help='work directory (default: a new temporary one)')
    args = p.parse_args()
    old = os.path.abspath(args.old)
    exe = os.path.abspath(args.exe)
    work = os.path.abspath(args.work or tempfile.mkdtemp(prefix='compat-'))
    for sub in ('restart', 'run'):
        os.makedirs(os.path.join(work, sub), exist_ok=True)

    # 1. every frame converts
    saves, is_vtkhdf = frames(old, 'save')
    tool = '2vtk.py' if is_vtkhdf else '2vtk-legacy.py'
    convert(tool, os.path.join(old, MODEL), os.path.join(work, 'vtk'), saves)
    print('1. %s converted %d %s frames' % (tool, len(saves), '.vtkhdf' if is_vtkhdf else 'des-binary'))

    # 2. the last checkpoint converts
    chkpts, _ = frames(old, 'chkpt')
    last = chkpts[-1]
    msg = run([sys.executable, os.path.join(ROOT, 'utils', 'convert_checkpoint.py'),
               os.path.join(old, MODEL), str(last), '-o', os.path.join(work, 'restart')],
              work, 'convert_checkpoint.py')
    print('2. ' + msg.splitlines()[0])

    # 3. the current DES restarts from it and writes one more finite frame
    step = info_steps(os.path.join(old, MODEL + '.info'), last)
    keys = ('[sim]\nis_restarting = yes\nrestarting_from_modelname = %s\nrestarting_from_frame = %d\n'
            % (os.path.join(work, 'restart', MODEL), last))
    cfg = open(CFG).read()
    cfg = re.sub(r'^max_steps = .*$', 'max_steps = %d' % (step + 100), cfg, flags=re.M)
    cfg = re.sub(r'^modelname = .*$', 'modelname = rst', cfg, flags=re.M)
    if '[sim]\n' not in cfg:
        fail('%s has no "[sim]" line for the restart keys' % CFG)
    cfg = cfg.replace('[sim]\n', keys, 1)
    with open(os.path.join(work, 'run', 'rst.cfg'), 'w') as f:
        f.write(cfg)
    run([exe, 'rst.cfg'], os.path.join(work, 'run'), os.path.basename(exe))
    new = os.path.join(work, 'run', 'rst.save.%06d.vtkhdf' % (last + 1))
    if not os.path.isfile(new):
        fail('the restart wrote no %s' % new)
    if os.path.isfile(os.path.join(work, 'run', 'rst.save.000000.vtkhdf')):
        fail('the run wrote frame 0: it started from step 0, not from the checkpoint')
    with h5py.File(new, 'r') as h:
        new_step = int(h['steps'][0])
        bad = [k for k in h.keys() if isinstance(h.get(k), h5py.Dataset)
               and h[k].dtype.kind == 'f' and not np.all(np.isfinite(h[k][()]))]
    if new_step != step + 100:
        fail('%s is step %d, not %d' % (new, new_step, step + 100))
    if bad:
        fail('%s has non-finite values in %s' % (new, ', '.join(bad)))
    print('3. %s restarted from frame %d (step %d) and wrote frame %d, all fields finite'
          % (os.path.basename(exe), last, step, last + 1))

    # 4. the current 2vtk.py converts the restarted run's new frame
    convert('2vtk.py', os.path.join(work, 'run', 'rst'), os.path.join(work, 'vtk-rst'), [last + 1])
    print('4. 2vtk.py converted the restarted frame %d' % (last + 1))


if __name__ == '__main__':
    main()
