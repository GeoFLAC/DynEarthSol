#!/usr/bin/python3
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


def show_msg(kind,max,sigma):
    if max+sigma > 1.e-8:
        print('  %s:\t\t%.3e %.3e (> 1.e-8)'%(kind,max, sigma))
        inc = 1
    else:
        print('  %s:\t\t%.3e %.3e'%(kind,max, sigma))
        inc = 0
    return inc


def reldiff_and_show_msg(oldf, newf, kind):
    if oldf.size != newf.size:
        print('  %s:\t\t%-.d -> %-d (size mismatch)'%(kind, oldf.size, newf.size))
        return 1
    else:
        max, sigma = reldiff(oldf, newf)
        return show_msg(kind, max, sigma)


def compare(old, new):
    inc = 0
    
    inc += reldiff_and_show_msg(old.T, new.T, 'Temperature')

    inc += reldiff_and_show_msg(old.x, new.x, 'X coordinate')

    inc += reldiff_and_show_msg(old.z, new.z, 'Z coordinate')

    inc += reldiff_and_show_msg(old.vx, new.vx, 'X velocity')

    inc += reldiff_and_show_msg(old.vz, new.vz, 'Z velocity')

    inc += reldiff_and_show_msg(old.pls, new.pls, 'Pl. strain')

    inc += reldiff_and_show_msg(old.tI, new.tI, 'Stress I')

    inc += reldiff_and_show_msg(old.tII, new.tII, 'Stress II')

    inc += reldiff_and_show_msg(old.sI, new.sI, 'Strain I')

    inc += reldiff_and_show_msg(old.sII, new.sII, 'Strain II')

    inc += reldiff_and_show_msg(old.srI, new.srI, 'S. rate I')

    inc += reldiff_and_show_msg(old.srII, new.srII, 'S. rate II')

    inc += reldiff_and_show_msg(old.visc, new.visc, 'Viscosity')

    inc += reldiff_and_show_msg(old.m_x, new.m_x, 'Marker X')

    inc += reldiff_and_show_msg(old.m_z, new.m_z, 'Marker Z')

    inc += reldiff_and_show_msg(old.m_mat, new.m_mat, 'Marker Mat')

    inc += reldiff_and_show_msg(old.m_time, new.m_time, 'Marker Time')

    return inc


olddir = sys.argv[1]
curdir = os.getcwd()

if len(sys.argv) > 3:
    frame = int(sys.argv[3])
    newdir = sys.argv[2]
    modelname = 'result'
else:
    frame = int(sys.argv[2])
    newdir = curdir
    modelname = 'benchmark'

# name holder
old = 0
new = 0

markersetname = 'markerset'

try:
    # read old and new results

    os.chdir(olddir)
    des = Dynearthsol(modelname)
    old = read_data(des, frame)

    os.chdir(newdir)
    des = Dynearthsol(modelname)
    new = read_data(des, frame)

    # compare results
    print()
    print('Relative difference (max, stddev) of frame =', frame,
          ' step =', int(des.steps[frame]))
    print('  ---')
    inc = compare(old, new)
    print('')
    if inc == 0:
        print('  Status: Normal round-off error~')
    else:
        print('  Status: !!!!!!!!!! SOMETHING WRONG !!!!!!!!!!')
    print('  ---')


finally:
    # restort to original directory
    os.chdir(curdir)
