#!/usr/bin/env python
"""Compare two DynEarthSol output snapshots field-by-field.

Usage
-----
::

    python compare.py <path/to/old-modelname> <path/to/new-modelname> <frame>

    <path/to/old-modelname>  path including the model name prefix for the reference run
                             (e.g. orig-test-3d-equ-tiny.cfg/result)
    <path/to/new-modelname>  path including the model name prefix for the new run
                             (e.g. result  or  ~/data/jobs/test/result)
    <frame>                  output frame index to compare (integer, e.g. 4)

Exit codes
----------
0  no serious differences (bit-exact or round-off only)
1  at least one field has relative difference >= 1e-8, or a NaN/Inf was found

Restart troubleshooting
-----------------------
If a restarted run produces unexpected results, use compare.py to isolate
which field first diverges from the original run.  A same-name restart
(``modelname == restarting_from_modelname``) backs up the restart frame's
save file to ``<model>.<frame>.save.vtkhdf.old`` before overwriting it.
To compare that backed-up frame against the restarted output, rename or
copy the ``.old`` file to the expected name in a separate directory, then::

    python compare.py <dir-with-original>/result result <restart-frame>

A ``Status: BIT-EXACT`` result rules out restart divergence as the cause;
the first field with a large relative difference is the starting point for
debugging.
"""

from __future__ import print_function
import sys, os
import numpy as np
sys.path.append('../')
from Dynearthsol import Dynearthsol

def first_invariant(t):
    nstr = t.shape[1]
    ndims = 2 if (nstr == 3) else 3
    return np.sum(t[:,:ndims], axis=1) / ndims


def second_invariant(t):
    '''The second invariant of the deviatoric part of a symmetric tensor t,
    where t[:,0:ndims] are the diagonal components;
      and t[:,ndims:] are the off-diagonal components.'''
    nstr = t.shape[1]

    # second invariant: sqrt(0.5 * t_ij**2)
    if nstr == 3:  # 2D
        return np.sqrt(0.25 * (t[:,0] - t[:,1])**2 + t[:,2]**2)
    else:  # 3D
        a = (t[:,0] + t[:,1] + t[:,2]) / 3
        return np.sqrt( 0.5 * ((t[:,0] - a)**2 + (t[:,1] - a)**2 + (t[:,2] - a)**2) +
                        t[:,3]**2 + t[:,4]**2 + t[:,5]**2)


class Stuff():
    pass


def read_data(des, frame):
    stuff = Stuff()

    stuff.T = des.read_field(frame,'temperature')
    coordinate = des.read_field(frame, 'coordinate')
    stuff.x = np.array(coordinate[:,0])
    stuff.z = np.array(coordinate[:,-1])
    velocity = des.read_field(frame, 'velocity')
    stuff.vx = np.array(velocity[:,0])
    stuff.vz = np.array(velocity[:,-1])
    stuff.pls = des.read_field(frame,'plastic strain')

    stress = des.read_field(frame, 'stress')
    stuff.tI = first_invariant(stress)
    stuff.tII = second_invariant(stress)
    strain = des.read_field(frame, 'strain')
    stuff.sI = first_invariant(strain)
    stuff.sII = second_invariant(strain)
    strain_rate = des.read_field(frame, 'strain-rate')
    stuff.srI = first_invariant(strain_rate)
    stuff.srII = second_invariant(strain_rate)

    stuff.visc = des.read_field(frame, 'viscosity')

    marker_data = des.read_markers(frame, markersetname)
    field = marker_data[markersetname + '.coord']
    stuff.m_x = field[:,0]
    stuff.m_z = field[:,1]
    stuff.m_id = marker_data[markersetname + '.id']
    stuff.m_mat = marker_data[markersetname + '.mattype']
    stuff.m_time = marker_data[markersetname + '.time']
    return stuff


def reldiff(oldf, newf):
    m = np.abs(oldf).max()
    diff = np.abs(newf - oldf)
    if m == 0.:
        return diff.max(), diff.std()
    else:
        return diff.max()/m, diff.std()/m


def show_msg(kind, max, sigma):
    """Print one field comparison line; return (fail, nonzero).

    fail=1  if max+sigma > 1e-8  (serious divergence)
    nonzero=1 if max+sigma > 0   (any difference, even round-off)
    Both zero means the field is bit-identical.
    """
    if not np.isfinite(max) or not np.isfinite(sigma):
        print('  %s:\t\t%s %s (NaN/Inf — field corrupt)' % (kind, max, sigma))
        return 1, 1
    if max + sigma > 1.e-8:
        print('  %s:\t\t%.3e %.3e (> 1.e-8)' % (kind, max, sigma))
        return 1, 1
    elif max + sigma > 0.:
        print('  %s:\t\t%.3e %.3e' % (kind, max, sigma))
        return 0, 1
    else:
        print('  %s:\t\t%.3e %.3e' % (kind, max, sigma))
        return 0, 0


def reldiff_and_show_msg(oldf, newf, kind):
    if oldf.size != newf.size:
        print('  %s:\t\t%-.d -> %-d (size mismatch)'%(kind, oldf.size, newf.size))
        return 1, 1
    else:
        max, sigma = reldiff(oldf, newf)
        return show_msg(kind, max, sigma)


def compare(old, new):
    """Compare all fields; return (n_fail, n_nonzero).

    n_fail    — fields with relative diff >= 1e-8 (serious)
    n_nonzero — fields with any nonzero diff (round-off or worse)
    Caller uses these to distinguish: bit-exact / round-off / seriously wrong.
    """
    n_fail = n_nonzero = 0
    for a, b, kind in [
            (old.T,      new.T,      'Temperature'),
            (old.x,      new.x,      'X coordinate'),
            (old.z,      new.z,      'Z coordinate'),
            (old.vx,     new.vx,     'X velocity'),
            (old.vz,     new.vz,     'Z velocity'),
            (old.pls,    new.pls,    'Pl. strain'),
            (old.tI,     new.tI,     'Stress I'),
            (old.tII,    new.tII,    'Stress II'),
            (old.sI,     new.sI,     'Strain I'),
            (old.sII,    new.sII,    'Strain II'),
            (old.srI,    new.srI,    'S. rate I'),
            (old.srII,   new.srII,   'S. rate II'),
            (old.visc,   new.visc,   'Viscosity'),
            (old.m_x,    new.m_x,    'Marker X'),
            (old.m_z,    new.m_z,    'Marker Z'),
            (old.m_mat,  new.m_mat,  'Marker Mat'),
            (old.m_time, new.m_time, 'Marker Time')]:
        f, nz = reldiff_and_show_msg(a, b, kind)
        n_fail += f; n_nonzero += nz
    return n_fail, n_nonzero


if len(sys.argv) == 4:
    oldmodelname = sys.argv[1]
    modelname = sys.argv[2]
    frame = int(sys.argv[3])
else:
    print("Usage:")
    print("  python compare.py <path/to/old-modelname> <path/to/new-modelname> <frame>")
    sys.exit(1)

markersetname = 'markerset'

# read old and new results

try:
    des_old = Dynearthsol(oldmodelname)
except OSError:
    print("Error: Directory of old results doesn't exist:", oldmodelname)
    sys.exit(1)
old = read_data(des_old, frame)
fmt_old = des_old.format

try:
    des_new = Dynearthsol(modelname)
except OSError:
    print("Error: Directory of new results doesn't exist:", modelname)
    sys.exit(1)
new = read_data(des_new, frame)
fmt_new = des_new.format

# compare results
print()
if fmt_old != fmt_new:
    print('Comparison between :', fmt_old, 'and', fmt_new, '(new)')
print('Relative difference (max, stddev) of frame = %d  step = %d' \
            % (frame, int(des_old.steps[frame])))
print('  ---')
n_fail, n_nonzero = compare(old, new)
print('')
if n_fail > 0:
    print('  Status: !!!!!!!!!! SOMETHING WRONG !!!!!!!!!!')
elif n_nonzero > 0:
    print('  Status: Normal round-off error~')
else:
    print('  Status: BIT-EXACT (all fields identical)')
print('  ---')

sys.exit(0 if n_fail == 0 else 1)
