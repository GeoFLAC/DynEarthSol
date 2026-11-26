#!/usr/bin/env python3
"""Write a .poly file for a wedge with two offset decollements.

usage: make_ramp_poly.py xa xb out.poly   (xa, xb in km)

A 40 x 5 km box. The lower decollement (mattype 1) is the bottom 1 km from
the pushing wall (x = 0) to x = xa; the upper decollement (mattype 1) lies
between 1.5 and 2.5 km depth from x = xb to the right wall. Everything else
is the brittle layer (mattype 0). xa > xb: the two overlap by xa - xb km;
xa < xb: a gap of xb - xa km separates them.

   0------------------------------------7      z = 0
   |                10-----------------6       z = -1.5  upper decollement
   |                9------------------5       z = -2.5  (x >= xb)
   1-----------8                       |       z = -4    lower decollement
   2-----------3-----------------------4       z = -5    (x <= xa)
"""
import sys

xa, xb, out = float(sys.argv[1]) * 1e3, float(sys.argv[2]) * 1e3, sys.argv[3]
L = 40e3
nodes = [(0, 0), (0, -4e3), (0, -5e3), (xa, -5e3), (L, -5e3), (L, -2.5e3), (L, -1.5e3), (L, 0),
         (xa, -4e3), (xb, -2.5e3), (xb, -1.5e3)]
segs = [(0, 1, 1), (1, 2, 1),                     # left wall (x0)
        (2, 3, 16), (3, 4, 16),                   # bottom (z0)
        (4, 5, 2), (5, 6, 2), (6, 7, 2),          # right wall (x1)
        (7, 0, 32),                               # top (z1)
        (1, 8, 0), (8, 3, 0),                     # lower decollement outline
        (9, 5, 0), (10, 6, 0), (9, 10, 0)]        # upper decollement outline
area = 1.5 * 250.0**2                             # as meshing_option = 1 at resolution 250
regions = [((xa / 2, -4.5e3), 1), (((xb + L) / 2, -2e3), 1), ((1e3, -0.5e3), 0)]

with open(out, 'w') as f:
    f.write(__doc__.replace('\n', '\n# ').join(['# ', '\n']))
    f.write(f'# xa = {xa/1e3:g} km, xb = {xb/1e3:g} km\n')
    f.write('#### node coordinates ####\n# npoints ndims 0 0\n  %d 2 0 0\n' % len(nodes))
    for i, (x, z) in enumerate(nodes):
        f.write('  %2d %10.1f %10.1f\n' % (i, x, z))
    f.write('#### segments ####\n# nsegments 1\n  %d 1\n' % len(segs))
    for j, (a, b, flag) in enumerate(segs):
        f.write('  %2d %2d %2d %2d\n' % (j, a, b, flag))
    f.write('#### holes, must be 0 ####\n0\n')
    f.write('#### regions ####\n# nregions\n  %d\n# k  xk  zk  mattype  size\n' % len(regions))
    for k, ((x, z), m) in enumerate(regions):
        f.write('  %d %10.1f %10.1f %d %g\n' % (k, x, z, m, area))
