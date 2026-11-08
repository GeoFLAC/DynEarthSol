#!/usr/bin/env python3
"""Check DES initial conditions from a Geodynamic World Builder world against GWB itself.

    check_gwb.py --exe ../../dynearthsol2d --gwb-dat GWB_DAT gwb/gwb-2d.cfg
    check_gwb.py --exe ../../dynearthsol2d --unsupported gwb/gwb-2d.cfg

The first form runs the cfg for 0 steps (temperature_option, mattype_option and
weakzone_option 936) and compares its frame 0 with gwb-dat, GWB's point-query
application, at the same points:
  - every node temperature, to gwb-dat's 6 printed digits;
  - every marker's mattype with the composition gwb-dat gives it (composition i
    is mattype i; each marker must have one of at least 0.5);
  - plastic strain: weakzone_plstrain exactly in the elements whose centroid
    gwb-dat tags as a "fault" feature, 0 elsewhere.
It runs again with ic.gwb_fault_is_gaussian and weakzone_standard_deviation a
quarter of the first fault's thickness: the same elements are weak, with values in (0,
weakzone_plstrain] that vary across the fault.
It then checks that a world without a fault (exit 12), a world leaving a marker
without a material (exit 12), a missing world (exit 20) and a restart with
is_restarting_weakzone (exit 11) stop the run, as do two faults of one name with
ic.gwb_fault_is_gaussian (exit 12).
--unsupported checks that a build without usegwb=1 refuses option 936 (exit 31).

gwb-dat is queried in 3D, at (x, 0, z) for a 2D model, whose world must then
have the cross section [[0,0],[L,0]]: before GWB PR #949 its 2D output prints
the compositions one column off.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '..'))
from Dynearthsol import Dynearthsol

ap = argparse.ArgumentParser()
ap.add_argument('cfg')
ap.add_argument('--exe', required=True)
ap.add_argument('--gwb-dat')
ap.add_argument('--unsupported', action='store_true')
a = ap.parse_args()
exe = os.path.abspath(a.exe)
cfg_text = open(a.cfg).read()


def key(name):
    m = re.search(r'^\s*%s\s*=\s*(\S+)' % re.escape(name), cfg_text, re.M)
    return m.group(1) if m else None


ndims = 3 if exe.endswith('3d') else 2
world = os.path.join(os.path.dirname(os.path.abspath(a.cfg)), key('gwb_filename'))
zlength = float(key('zlength'))
nmat = int(key('num_materials'))
plstrain = float(key('weakzone_plstrain'))
failures = []


def run(tmp, world_path, weakzone_option=None, extra=''):
    """Run the cfg for 0 steps in tmp with gwb_filename = world_path, and extra appended;
    return the exit code."""
    text = re.sub(r'^(\s*gwb_filename\s*=).*$', r'\1 ' + world_path, cfg_text, flags=re.M)
    text = re.sub(r'^\s*max_time_in_yr\s*=.*$', 'max_steps = 0', text, flags=re.M)
    if weakzone_option is not None:
        text = re.sub(r'^(\s*weakzone_option\s*=).*$', r'\g<1> %d' % weakzone_option, text, flags=re.M)
    with open(os.path.join(tmp, 'run.cfg'), 'w') as f:
        f.write(text + extra)
    with open(os.path.join(tmp, 'run.log'), 'w') as log:
        return subprocess.call([exe, 'run.cfg'], cwd=tmp, stdout=log, stderr=subprocess.STDOUT)


def expect(what, ok, detail=''):
    print('%s %s%s' % ('PASS' if ok else 'FAIL', what, (': ' + detail) if detail else ''))
    if not ok:
        failures.append(what)


def gwb_query(points, ncomp):
    """gwb-dat temperature, compositions 0..ncomp-1 and tag at DES points (n x ndims)."""
    with tempfile.NamedTemporaryFile('w', suffix='.dat', delete=False) as f:
        f.write('# dim = 3\n# compositions = %d\n' % ncomp)
        y = points[:, 1] if ndims == 3 else np.zeros(len(points))
        z = points[:, ndims - 1]
        np.savetxt(f, np.column_stack([points[:, 0], y, zlength + z, -z]), fmt='%.17g')
        name = f.name
    out = subprocess.run([a.gwb_dat, world, name], capture_output=True, text=True, check=True).stdout
    os.unlink(name)
    v = np.array([l.split() for l in out.splitlines() if l and not l.startswith('#')], float)
    # columns: x y z d T vx vy vz c0 .. c(ncomp-1) tag
    return v[:, 4], v[:, 8:8 + ncomp], v[:, 8 + ncomp].astype(int)


def fault_tags(world_json):
    """GWB tag indices of the fault features, as gwb.cxx derives them."""
    tags, faults = [], set()
    for feat in world_json['features']:
        tag = feat.get('tag') or feat['model']
        if tag not in tags:
            tags.append(tag)
        if feat['model'] == 'fault':
            faults.add(tags.index(tag))
    return faults


tmp = tempfile.mkdtemp(prefix='gwb-check-')
try:
    if a.unsupported:
        code = run(tmp, world)
        expect('a build without usegwb=1 refuses option 936 (exit 31)', code == 31, 'exit %d' % code)
    else:
        wj = json.load(open(world))
        if ndims == 2:
            expect('2D world cross section is [[0,0],[L,0]]',
                   wj['cross section'][0] == [0, 0] and wj['cross section'][1][1] == 0)
        code = run(tmp, world)
        expect('DES runs with option 936', code == 0, 'exit %d' % code)
        if code == 0:
            d = Dynearthsol(os.path.join(tmp, key('modelname')))
            xyz = d.read_field(0, 'coordinate')
            T = d.read_field(0, 'temperature')
            Tg, _, _ = gwb_query(xyz, nmat)
            err = np.abs(T - Tg) / np.abs(Tg)
            expect('node temperatures match gwb-dat (%d nodes)' % len(T), err.max() < 1e-5,
                   'max relative difference %.2g' % err.max())

            mk = d.read_markers(0, 'markerset')
            mt = mk['markerset.mattype']
            _, comp, _ = gwb_query(mk['markerset.coord'], nmat)
            expect('every marker has a GWB composition of at least 0.5 (%d markers)' % len(mt),
                   (comp.max(axis=1) >= 0.5).all())
            bad = (mt != comp.argmax(axis=1)).sum()
            expect('marker mattypes match gwb-dat', bad == 0,
                   '%d of %d differ; counts %s' % (bad, len(mt), np.bincount(mt).tolist()))

            conn = d.read_field(0, 'connectivity')
            pls = d.read_field(0, 'plastic strain')
            _, _, tag = gwb_query(xyz[conn].mean(axis=1), nmat)
            want = np.where(np.isin(tag, list(fault_tags(wj))), plstrain, 0.0)
            nweak = int((want > 0).sum())
            bad = int((pls != want).sum())
            expect('weak zone = elements in GWB faults (%d of %d)' % (nweak, len(pls)),
                   nweak > 0 and bad == 0, '%d elements differ' % bad)

        fault = [f for f in wj['features'] if f['model'] == 'fault'][0]
        sigma = fault['segments'][0]['thickness'][0] / 4
        gauss = '\n[ic]\ngwb_fault_is_gaussian = yes\nweakzone_standard_deviation = %g\n' % sigma
        code = run(tmp, world, extra=gauss)
        expect('DES runs with gwb_fault_is_gaussian, weakzone_standard_deviation = %g' % sigma, code == 0, 'exit %d' % code)
        if code == 0:
            d = Dynearthsol(os.path.join(tmp, key('modelname')))
            xyz = d.read_field(0, 'coordinate')
            conn = d.read_field(0, 'connectivity')
            pls = d.read_field(0, 'plastic strain')
            _, _, tag = gwb_query(xyz[conn].mean(axis=1), nmat)
            infault = np.isin(tag, list(fault_tags(wj)))
            w = pls[infault]
            expect('Gaussian weak zone: the elements in GWB faults, values in (0, %g] (%d elements)'
                   % (plstrain, infault.sum()),
                   infault.any() and (pls[~infault] == 0).all() and (w > 0).all() and (w <= plstrain).all()
                   and w.min() < 0.5 * w.max(),
                   'values %.3g..%.3g' % (w.min(), w.max()) if infault.any() else 'none')

        # error cases, each a world written next to the run
        nofault = dict(wj, features=[f for f in wj['features'] if f['model'] != 'fault'])
        json.dump(nofault, open(os.path.join(tmp, 'nofault.wb'), 'w'))
        code = run(tmp, os.path.join(tmp, 'nofault.wb'))
        expect('a world without a fault stops with weakzone_option 936 (exit 12)', code == 12, 'exit %d' % code)

        nomat = dict(nofault, features=[dict(f, **{'composition models': []}) for f in nofault['features']])
        json.dump(nomat, open(os.path.join(tmp, 'nomat.wb'), 'w'))
        code = run(tmp, os.path.join(tmp, 'nomat.wb'), weakzone_option=0)
        expect('a world leaving markers without a material stops (exit 12)', code == 12, 'exit %d' % code)

        twice = dict(wj, features=wj['features'] + [fault])
        json.dump(twice, open(os.path.join(tmp, 'twice.wb'), 'w'))
        code = run(tmp, os.path.join(tmp, 'twice.wb'), extra=gauss)
        expect('two faults of one name stop a Gaussian weak zone (exit 12)', code == 12, 'exit %d' % code)

        code = run(tmp, os.path.join(tmp, 'missing.wb'))
        expect('a missing world stops (exit 20)', code == 20, 'exit %d' % code)

        code = run(tmp, world, extra='\n[sim]\nis_restarting = yes\nrestarting_from_modelname = %s\n'
                   'restarting_from_frame = 0\n[ic]\nis_restarting_weakzone = yes\n' % key('modelname'))
        expect('a restart re-seeding the weak zone with weakzone_option 936 stops (exit 11)',
               code == 11, 'exit %d' % code)
finally:
    shutil.rmtree(tmp)

if failures:
    sys.exit('%d check(s) failed' % len(failures))
