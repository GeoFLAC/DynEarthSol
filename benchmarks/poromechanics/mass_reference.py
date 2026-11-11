#!/usr/bin/env python3
"""Independent saturated finite-grain mass and sensible-heat cell references.

These references describe a constant-compressibility fluid and the Eulerian
Biot porosity law, with tension-positive logarithmic bulk strain. They do not
establish that DES implements the finite-density balance. Pressure work and
thermoelastic caloric effects are outside this sensible-heat model.
"""
from dataclasses import dataclass
import json
import math

from scipy.integrate import quad
from scipy.optimize import brentq
import numpy as np


class Transport:
    """Independent semidiscrete saturated transport on an accepted strain path.

    Shape gradients come from inverse vertex matrices. An antisymmetric face
    matrix transports mass and upwind sensible heat, and an independent Radau
    integrator supplies the temporal reference. Mechanical strain is an input
    to this transport-stage check, not an independently solved equilibrium.
    """
    def __init__(self, model, permeability=1e-10):
        self.cell = Cell()
        self.x = model.read_field(0,'coordinate')
        self.conn = model.read_field(0,'connectivity')
        dim = self.x.shape[1]
        self.nn, self.ne = len(self.x), len(self.conn)
        vertices = np.concatenate((np.ones((*self.conn.shape,1)),self.x[self.conn]),axis=2)
        self.gradient = np.linalg.inv(vertices)[:,1:,:].swapaxes(1,2)
        self.volume = np.abs(np.linalg.det(self.x[self.conn][:,1:]-self.x[self.conn][:,:1]))/math.factorial(dim)
        eps0 = model.read_field(0,'strain')[:,:dim].sum(axis=1)
        eps1 = model.read_field(model.frames[-1],'strain')[:,:dim].sum(axis=1)
        self.v0 = self.volume*np.exp(eps0)
        self.de = eps1-eps0
        self.p0 = model.read_field(0,'pore pressure')
        self.t0 = model.read_field(0,'temperature')
        self.phi0 = model.read_field(0,'porosity')
        self.cs = model.read_field(0,'solid heat capacity')
        self.fixed = np.isclose(self.x[:,0],0)|np.isclose(self.x[:,0],1)
        self.tensor = np.array([[2,.1,.3],[.1,1.2,.2],[.3,.2,1.5]])
        if dim==2:self.tensor=self.tensor[np.ix_([0,2],[0,2])]
        self.mobility = permeability/1e-3
        self.duration = model.time[-1]
        self.lag = self.scatter(self.v0[:,None]/self.conn.shape[1]*
            self.density(self.p0,self.t0)[self.conn]*self.cell.alpha**2/(2e7+(1e7/3 if dim==2 else 0)))
        self.sources = []
        for xyz, rate, window in (([.25]*(dim-1)+[-.25],.001,(.13,.47)),
                                   ([.75]*(dim-1)+[-.75],-.0005,(0.,1.))):
            for e, matrix in enumerate(vertices):
                weights = np.linalg.solve(matrix.T,np.r_[1.,xyz])
                if weights.min() >= -1e-12:
                    nodal = np.zeros(self.nn)
                    np.add.at(nodal,self.conn[e],weights*rate)
                    self.sources.append((nodal,window))
                    break
            else:raise AssertionError('Reference point lies outside the mesh')

    def scatter(self, values):
        result = np.zeros(self.nn)
        np.add.at(result,self.conn,np.broadcast_to(values,self.conn.shape))
        return result

    def density(self, p, t):
        return self.cell.rho0*np.exp(p/self.cell.fluid_bulk-self.cell.fluid_expansion*(t-273.))

    def state(self, p, t, fraction):
        c = self.cell
        phi = c.alpha-(c.alpha-self.phi0[:,None])*np.exp(-fraction*self.de[:,None]-
            (p[self.conn]-self.p0[self.conn].mean(axis=1)[:,None])/c.grain_bulk)
        rho = self.density(p,t)[self.conn]
        weight = self.v0[:,None]*np.exp(fraction*self.de[:,None])/self.conn.shape[1]
        mass = self.scatter(weight*rho*phi)
        mp = self.scatter(weight*rho*(phi/c.fluid_bulk+(c.alpha-phi)/c.grain_bulk))
        mt = self.scatter(-weight*rho*phi*c.fluid_expansion)
        moving = self.scatter(weight*rho*c.alpha*self.de[:,None])/self.duration
        return mass,mp,mt,moving

    def rates(self, p, t, start, stop, hot=False):
        local_p, local_t = p[self.conn],t[self.conn]
        gradient_p = np.einsum('ejd,ej->ed',self.gradient,local_p)
        q = -self.mobility*np.einsum('ab,eb->ea',self.tensor,gradient_p)
        delta = self.gradient[:,:,None,:]-self.gradient[:,None,:,:]
        volume_flux = -self.volume[:,None,None]/self.conn.shape[1]*np.einsum('eijd,ed->eij',delta,q)
        rho = self.density(p,t)[self.conn]
        upwind_rho = np.where(volume_flux>=0,rho[:,:,None],rho[:,None,:])
        upwind_t = np.where(volume_flux>=0,local_t[:,:,None],local_t[:,None,:])
        mass_flux = volume_flux*upwind_rho
        mrate = self.scatter(-mass_flux.sum(axis=2))
        conduction = 2*self.volume[:,None,None]*np.einsum('eid,ejd->eij',self.gradient,self.gradient)
        erate = self.scatter(-(self.cell.fluid_heat*mass_flux*upwind_t).sum(axis=2)-
                             np.einsum('eij,ej->ei',conduction,local_t))
        for i,(source,window) in enumerate(self.sources):
            a,b = np.array(window)*self.duration
            ratio = float(a<=start<b) if start==stop else max(0.,min(stop,b)-max(start,a))/(stop-start)
            source_t = np.full(self.nn,333.) if hot and i==0 else t
            flow = ratio*source*self.density(p,source_t)
            mrate += flow
            erate += self.cell.fluid_heat*source_t*flow
        return mrate,erate

    def solve(self):
        from scipy.integrate import solve_ivp
        def derivative(time,y):
            p,t = y[:self.nn],y[self.nn:]
            mass,mp,mt,moving = self.state(p,t,time/self.duration)
            rm,re = self.rates(p,t,phase,phase)
            cf = self.cell.fluid_heat
            capacity = self.cs+cf*mass
            a,b,c,d = mp+self.lag,mt,cf*t*mp,capacity+cf*t*mt
            determinant = a*d-b*c
            pdot = ((rm-moving)*d-b*(re-cf*t*moving))/determinant
            tdot = (a*(re-cf*t*moving)-c*(rm-moving))/determinant
            pdot[self.fixed] = 0
            tdot[self.fixed] = (re-cf*t*rm)[self.fixed]/capacity[self.fixed]
            return np.r_[pdot,tdot]
        y = np.r_[self.p0,self.t0]
        for a,b in zip((0.,.13,.47),(.13,.47,1.)):
            phase = .5*(a+b)*self.duration
            solution = solve_ivp(derivative,(a*self.duration,b*self.duration),y,method='Radau',
                                 rtol=2e-10,atol=np.r_[np.full(self.nn,1e-5),np.full(self.nn,1e-10)])
            assert solution.success,solution.message
            y = solution.y[:,-1]
        return y[:self.nn],y[self.nn:]


@dataclass(frozen=True)
class Cell:
    alpha: float = .98
    phi0: float = .25
    grain_bulk: float = 1e9
    fluid_bulk: float = 1e8
    rho0: float = 1000.
    fluid_expansion: float = 2e-4
    solid_expansion: float = 0.
    reference_temperature: float = 273.
    reference_volume: float = 1.
    solid_density: float = 2650.
    solid_heat: float = 1000.
    fluid_heat: float = 4180.

    @property
    def initial_mass(self):
        return self.rho0*self.phi0*self.reference_volume

    @property
    def solid_capacity(self):
        return (1-self.phi0)*self.solid_density*self.solid_heat*self.reference_volume

    def porosity(self, pressure, temperature, strain=0., closure=0.):
        theta = temperature-self.reference_temperature
        return self.alpha-(self.alpha-self.phi0)*math.exp(
            -strain-pressure/self.grain_bulk+self.solid_expansion*theta+closure)

    def density(self, pressure, temperature):
        return self.rho0*math.exp(pressure/self.fluid_bulk-
            self.fluid_expansion*(temperature-self.reference_temperature))

    def mass(self, pressure, temperature, strain=0., closure=0.):
        return self.reference_volume*math.exp(strain)*self.porosity(
            pressure, temperature, strain, closure)*self.density(pressure, temperature)

    def energy(self, pressure, temperature, strain=0., closure=0.):
        return (self.solid_capacity+self.fluid_heat*self.mass(
            pressure, temperature, strain, closure))*temperature

    def sealed_pressure(self, temperature, strain=0., closure=0.):
        theta = temperature-self.reference_temperature
        if self.fluid_bulk == self.grain_bulk:
            return self.fluid_bulk*(math.log((self.phi0*math.exp(self.fluid_expansion*theta)+
                (self.alpha-self.phi0)*math.exp(self.solid_expansion*theta+closure))/self.alpha)-strain)
        lower = self.grain_bulk*(math.log((self.alpha-self.phi0)/self.alpha)-strain+
                                 self.solid_expansion*theta+closure)
        upper = max(0., self.fluid_bulk*.1, lower+self.fluid_bulk)
        while self.mass(upper, temperature, strain, closure) < self.initial_mass:
            upper = lower+2*(upper-lower)
        return brentq(lambda p: self.mass(p, temperature, strain, closure)-self.initial_mass,
                       lower, upper, xtol=1e-7)

    def drained_heat(self, temperature):
        """Supplied heat at p=0, fixed volume and zero solid expansion."""
        assert self.solid_expansion == 0
        theta, a = temperature-self.reference_temperature, self.fluid_expansion
        integral = -math.expm1(-a*theta)/a if a else theta
        return self.solid_capacity*theta+self.fluid_heat*self.initial_mass*integral

    def drained_outflow_heat(self, temperature):
        assert self.solid_expansion == 0
        theta, a = temperature-self.reference_temperature, self.fluid_expansion
        if not a:
            return 0.
        lost = -math.expm1(-a*theta)
        integral = (lost-a*theta*math.exp(-a*theta))/a
        return self.fluid_heat*self.initial_mass*(self.reference_temperature*lost+integral)

    def injection_time(self, temperature, inlet_temperature, volume_rate):
        """Well-mixed drained cell; natural outflow is at cell temperature."""
        inflow_mass = self.density(0., inlet_temperature)*volume_rate
        span = inlet_temperature-self.reference_temperature
        upper = -math.log1p(-(temperature-self.reference_temperature)/span)
        # Logarithmic distance from Tin cancels the endpoint pole. The
        # quadrature stays smooth even when a root bracket approaches Tin.
        return quad(lambda s: (self.solid_capacity+self.fluid_heat*self.mass(0.,
            inlet_temperature-span*math.exp(-s)))/(self.fluid_heat*inflow_mass),
            0., upper, epsabs=1e-10, epsrel=1e-12)[0]


def checks():
    cell = Cell()
    worst = dict(sealed_mass=0., mass_tangent=0., energy_identity=0.)
    for strain in (0., -.01, .01):
        for temperature in (273., 283., 323.):
            p = cell.sealed_pressure(temperature, strain)
            worst['sealed_mass'] = max(worst['sealed_mass'], abs(
                cell.mass(p, temperature, strain)/cell.initial_mass-1))
            phi = cell.porosity(p, temperature, strain)
            tangent = cell.reference_volume*math.exp(strain)*cell.density(p, temperature)*(
                phi/cell.fluid_bulk+(cell.alpha-phi)/cell.grain_bulk)
            numeric = (cell.mass(p+10., temperature, strain)-cell.mass(p-10., temperature, strain))/20
            worst['mass_tangent'] = max(worst['mass_tangent'], abs(numeric/tangent-1))
    for theta in (1., 10., 50.):
        t = cell.reference_temperature+theta
        residual = cell.energy(0., t)-cell.energy(0., cell.reference_temperature)+\
                   cell.drained_outflow_heat(t)-cell.drained_heat(t)
        worst['energy_identity'] = max(worst['energy_identity'], abs(residual)/cell.drained_heat(t))
    equal = Cell(grain_bulk=1e8)
    for eps in (-.01, 0., .01):
        p = equal.sealed_pressure(equal.reference_temperature, eps)
        assert abs(p+equal.fluid_bulk*eps) < 3e-8
    # The incompressible, non-expanding mixing limit has a closed exponential.
    mixing = Cell(fluid_expansion=0.)
    inlet, rate, time = 323., .001, 50.
    capacity = mixing.solid_capacity+mixing.fluid_heat*mixing.initial_mass
    target = inlet-(inlet-mixing.reference_temperature)*math.exp(
        -mixing.fluid_heat*mixing.rho0*rate*time/capacity)
    mixing_error = abs(mixing.injection_time(target, inlet, rate)/time-1)
    assert worst['sealed_mass'] < 2e-14 and worst['mass_tangent'] < 2e-8
    assert worst['energy_identity'] < 1e-12 and mixing_error < 1e-12
    return dict(reference_only=True, passed=True, errors=worst, mixing_time_error=mixing_error)


if __name__ == '__main__':
    print(json.dumps(checks(), indent=2))
