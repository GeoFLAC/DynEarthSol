#!/usr/bin/python

import sys
from math import pi,sin, sqrt
import numpy
import matplotlib.pyplot as plt

sys.path.append('..')
from Dynearthsol import Dynearthsol


def analytical():
	K = 200e6
	mu = 200e6
	coh = 1e6
	phi = 10.0*pi/180
	psi = 10.0*pi/180
	ten = 5.67e6
	rho = 1.0
	vx = 1e-5

	e1 = K+4.0*mu/3.0
	e2 = K-2.0*mu/3.0

	sf = sin(phi)
	sp = sin(psi)
	nf = (1.0+sf)/(1.0-sf)
	np = (1.0+sp)/(1.0-sp)
	rl = (e1-e2*nf)/((e1+e2)*nf*np-2.0*e2*(nf+np)+2.0*e1)

	step1 = 2.0*coh*sqrt(nf)/((e1-e2*nf)*vx) # when yielding occurs

	displacement = vx * numpy.array(range(2001), dtype=float)
	Sxx = numpy.zeros(2001, dtype=float)
	for i in range(1, 2001):
		de = vx / (1 - displacement[i])
		if i < step1:
			Sxx[i] = Sxx[i-1] + e1*de
		else:
			Sxx[i] = Sxx[i-1] + de*(e1+2.0*rl*(e2*np-e1))

	return displacement, Sxx


def numerical():
	des = Dynearthsol('result')
	displacement = numpy.zeros(51, dtype=float)
	Sxx = numpy.zeros(51, dtype=float)
	for i in range(51):
		des.read_header(i)
		x = des.read_field(i, 'coordinate')
		displacement[i] = 1 - x[3,0]
		s = des.read_field(i, 'stress')
		Sxx[i] = abs(s[1,0])

	return displacement, Sxx


displacement, Sxx = analytical()
plt.clf()
plt.plot(displacement, Sxx * 1e-6)
plt.hold(1)

displacement, Sxx = numerical()
plt.plot(displacement, Sxx * 1e-6, 'o')

plt.xlabel('displacement (meter)')
plt.ylabel('Sxx (MPa)')
plt.axis((0, .02, 0, 8))
plt.legend(['Analytical', 'DynEarthSol'], loc='lower right')
plt.show()


