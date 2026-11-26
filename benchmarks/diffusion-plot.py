#!/usr/bin/env python
from __future__ import print_function, unicode_literals
import sys, os
from matplotlib.pylab import *
from scipy.special import erf

sys.path.append(os.getcwd() + '/..')

### 'import 2vtk' will not work because the module name starts with 2,
### use the following instead:
_tmp = __import__('2vtk', globals(), locals(), ['Dynearthsol'], -1)
Dynearthsol = _tmp.Dynearthsol

surfaceT = 273
mantleT = 1600

def readdata(frame):
    des = Dynearthsol('diffusion')
    time = des.time[frame]

    t = des.read_field(frame, 'temperature')
    x = des.read_field(frame, 'coordinate')
    return x, t, time


def myplot(t, x, xx, frame, time):
    diffusivity = 3./3000/1e3
    plot(t, x[:,-1]*1e-3, 'ok',
         surfaceT-(mantleT-surfaceT)*erf(xx/sqrt(4*diffusivity*time)), xx*1e-3, markerfacecolor='w')
    return



clf()
xx = linspace(0, -250e3)
myr2sec = 1e6 * 86400 * 365.2422
frames = np.array((0, 9, 59))  # at 1myr, 10myr, 60myr
initial_age = 1  # in myrs

for frame in frames:
    x, t, time = readdata(frame)
    age = time + initial_age*myr2sec  # in myrs
    myplot(t, x, xx, frame, age)

axis((surfaceT, mantleT+100, -250, 10))
xlabel('Temperature (Kelvin)', fontsize=16)
ylabel('Depth (km)', fontsize=16)
show()

