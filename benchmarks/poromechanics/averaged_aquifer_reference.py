"""Thickness-averaged poroelastic aquifer around a pumping well.

This generalises the regional subsidence model of Bear and Corapcioglu
(1981a, doi:10.1029/WR017i004p00937; 1981b, doi:10.1029/WR017i004p00947):

- Biot coefficient alpha and grain modulus K_s instead of Terzaghi stress
  with incompressible grains, and fluid compressibility kept throughout;
- an exact Laplace-domain solution for a finite well radius and a finite
  outer radius, for which the integration function g(t) of their eq. (92)
  is no longer zero.

Their closed forms are the limit of ``line_well`` with alpha=1,
K_f=K_s=inf and no well radius; ``bear_corapcioglu_1981b`` evaluates the
printed equations directly. The derivation and the list of approximations
are in doc/hydro/aquifer.md.

Conventions: tension-positive stress and strain; excess pore pressure p is
positive in compression, so extraction (rate > 0) makes p negative. Rates are
total well discharges in m^3/s. Subsidence is positive downward and the
radial displacement is positive outward. Vertical strain is the change of
aquifer thickness divided by the thickness.

mechanics='plane_stress' is Bear and Corapcioglu (1981b): no change of
vertical or vertical-shear total stress, horizontal displacement uniform over
the thickness. mechanics='uniaxial' is the Jacob/Bear and Corapcioglu (1981a)
vertical-compressibility model, with no horizontal displacement.

Requires NumPy and SciPy for benchmark evaluation only.
"""

import math

import numpy as np
from scipy.special import exp1, ive, kve


def _positive(name, value, allow_inf=False):
    if not (value > 0 and (allow_inf or math.isfinite(value))) or math.isnan(value):
        raise ValueError(f'{name} must be finite and positive')


class AveragedAquifer:
    """Homogeneous, isotropic, confined aquifer of constant thickness."""

    def __init__(self, *, thickness, bulk, shear, biot, porosity, fluid_bulk,
                 grain_bulk, permeability, viscosity, fluid_density=1000.0,
                 gravity=9.81, mechanics='plane_stress'):
        for name, value in (('Thickness', thickness), ('Bulk modulus', bulk),
                            ('Shear modulus', shear), ('Permeability', permeability),
                            ('Viscosity', viscosity), ('Fluid density', fluid_density),
                            ('Gravity', gravity)):
            _positive(name, value)
        # Infinite moduli give the incompressible limits used by the papers.
        _positive('Fluid bulk modulus', fluid_bulk, allow_inf=True)
        _positive('Grain bulk modulus', grain_bulk, allow_inf=True)
        if not 0 <= porosity <= biot <= 1 or biot == 0:
            raise ValueError('Require 0 <= porosity <= Biot coefficient <= 1 and Biot > 0')
        if mechanics not in ('plane_stress', 'uniaxial'):
            raise ValueError("mechanics must be 'plane_stress' or 'uniaxial'")
        self.thickness, self.biot, self.shear = thickness, biot, shear
        self.fluid_density, self.gravity = fluid_density, gravity
        self.mechanics = mechanics
        self.lame = bulk - 2 * shear / 3
        if self.lame + shear <= 0:
            raise ValueError('Require lambda + G > 0')
        # 1/M of Biot; the second term vanishes for incompressible grains.
        self.storage = porosity / fluid_bulk + (biot - porosity) / grain_bulk
        if mechanics == 'plane_stress':
            # alpha'' = 1/(lambda+G), Bear and Corapcioglu (1981b) eq. (94).
            self.matrix_compressibility = 1 / (self.lame + shear)
            # d(alpha*eps_v)/dh for the uniform areal-strain function h(t).
            self.coupling = 2 * shear * biot / (self.lame + 2 * shear)
        else:
            # alpha' = 1/(lambda+2G), Bear and Corapcioglu (1981a).
            self.matrix_compressibility = 1 / (self.lame + 2 * shear)
            self.coupling = 0.0
        self.specific_storage = self.storage + biot**2 * self.matrix_compressibility
        self.mobility = permeability / viscosity
        self.diffusivity = self.mobility / self.specific_storage
        self.transmissivity = self.mobility * fluid_density * gravity * thickness
        self.storativity = fluid_density * gravity * thickness * self.specific_storage
        # Plane-stress Lame constants of the averaged horizontal problem.
        self._lame_ps = 2 * shear * self.lame / (self.lame + 2 * shear)
        self._half_alpha_comp = biot / (2 * (self.lame + shear))

    # ------------------------------------------------------------------ fields
    def _fields(self, radius, pressure, integral, h, inner_constant, well_radius):
        """Strains and displacements from p, J=int_rw^r rho p drho and h(t)."""
        biot, lam, shear, B = self.biot, self.lame, self.shear, self.thickness
        if self.mechanics == 'uniaxial':
            areal = np.zeros_like(pressure)
            vertical = biot * self.matrix_compressibility * pressure
            radial = np.zeros_like(pressure)
        else:
            areal = self._half_alpha_comp * pressure + h
            vertical = (biot * pressure - lam * areal) / (lam + 2 * shear)
            radial = (inner_constant + self._half_alpha_comp * integral +
                      h * (radius**2 - well_radius**2) / 2) / radius
        volumetric = areal + vertical
        return {
            'pressure': pressure,
            'head_drawdown': -pressure / (self.fluid_density * self.gravity),
            'areal_strain': areal,
            'vertical_strain': vertical,
            'volumetric_strain': volumetric,
            'fluid_content': self.storage * pressure + biot * volumetric,
            'subsidence': -B * vertical,
            'radial_displacement': radial,
        }

    # -------------------------------------------------------------- line well
    def line_well(self, radius_m, time_s, rate_m3_s, well_radius=None):
        """Infinite aquifer with a line sink: the paper's closed forms.

        Pressure is the Theis solution with storativity S = rho g B S_s.
        The displacement satisfies U_r(well_radius)=0, eq. (90), exactly;
        well_radius=None is the paper's approximation u_w W(u_w)=0,
        exp(-u_w)=1. The pressure itself ignores the well radius.
        """
        if not math.isfinite(rate_m3_s):
            raise ValueError('Well rate must be finite')
        radius, time = np.broadcast_arrays(np.asarray(radius_m, dtype=float),
                                          np.asarray(time_s, dtype=float))
        if not np.isfinite(radius).all() or (radius <= 0).any():
            raise ValueError('Observation radii must be finite and positive')
        if not np.isfinite(time).all() or (time < 0).any():
            raise ValueError('Observation times must be finite and nonnegative')
        rw = 0.0 if well_radius is None else float(well_radius)
        if well_radius is not None and (not math.isfinite(rw) or rw <= 0 or (radius < rw).any()):
            raise ValueError('Well radius must be positive and inside every observation radius')
        pressure = np.zeros(radius.shape)
        integral = np.zeros(radius.shape)
        active = time > 0
        D, t, r = self.diffusivity, time[active], radius[active]
        amp = -rate_m3_s / (4 * math.pi * self.mobility * self.thickness)
        u = r**2 / (4 * D * t)
        pressure[active] = amp * exp1(u)
        # int r E1(r^2/4Dt) dr = 2Dt (u E1(u) - exp(-u)).
        if well_radius is None:
            lower = -1.0
        else:
            uw = rw**2 / (4 * D * t)
            lower = uw * exp1(uw) - np.exp(-uw)
        integral[active] = amp * 2 * D * t * (u * exp1(u) - np.exp(-u) - lower)
        return self._fields(radius, pressure, integral, 0.0, 0.0, rw)

    # ------------------------------------------------- finite well and domain
    def finite_well(self, radius_m, time_s, rate_m3_s, *, well_radius, outer_radius,
                    outer_hydraulic='fixed_pressure', inner_mechanics='fixed',
                    outer_mechanics='fixed', talbot_terms=32):
        """Exact averaged solution on well_radius <= r <= outer_radius.

        The well imposes the total rate on its face. outer_hydraulic is
        'fixed_pressure' (p=0) or 'no_flow'; inner/outer mechanics are
        'fixed' (U_r=0) or 'free' (sigma_rr=0). The uniform areal strain h(t)
        is solved with the pressure, which is the g(t) that the infinite-domain
        papers set to zero. Fixed Talbot inversion (Abate and Valko 2004).
        Returns the fields of ``line_well`` plus 'h' for each time row.
        """
        if outer_hydraulic not in ('fixed_pressure', 'no_flow'):
            raise ValueError("outer_hydraulic must be 'fixed_pressure' or 'no_flow'")
        for option in (inner_mechanics, outer_mechanics):
            if option not in ('fixed', 'free'):
                raise ValueError("mechanical boundaries must be 'fixed' or 'free'")
        if self.mechanics == 'uniaxial' and (inner_mechanics != 'fixed' or outer_mechanics != 'fixed'):
            raise ValueError('Uniaxial mechanics require fixed horizontal displacement at both boundaries')
        if not math.isfinite(rate_m3_s):
            raise ValueError('Well rate must be finite')
        rw, R = float(well_radius), float(outer_radius)
        if not (math.isfinite(rw) and math.isfinite(R) and 0 < rw < R):
            raise ValueError('Require 0 < well_radius < outer_radius')
        radius = np.atleast_1d(np.asarray(radius_m, dtype=float))
        time = np.atleast_1d(np.asarray(time_s, dtype=float))
        if (radius.ndim != 1 or not np.isfinite(radius).all() or (radius < rw).any()
                or (radius > R).any()):
            raise ValueError('Radii must be a 1D array inside the annulus')
        if time.ndim != 1 or not np.isfinite(time).all() or (time <= 0).any():
            raise ValueError('Times must be a positive 1D array')
        if not isinstance(talbot_terms, int) or talbot_terms < 8:
            raise ValueError('talbot_terms must be an integer >= 8')

        lam_ps, G = self._lame_ps, self.shear
        inner_c = 0.0 if inner_mechanics == 'fixed' else (lam_ps + 2 * G) * rw**2 / (2 * G)
        nr = radius.size
        values = np.empty((time.size, 2 * nr + 1))
        for i, t in enumerate(time):
            values[i] = _talbot(lambda s: self._laplace(
                s, radius, rate_m3_s, rw, R, outer_hydraulic, inner_c,
                outer_mechanics), t, talbot_terms)
        pressure, integral, h = values[:, :nr], values[:, nr:2 * nr], values[:, -1:]
        fields = self._fields(radius[None, :], pressure, integral, h, inner_c * h, rw)
        fields['h'] = h[:, 0]
        return fields

    def _laplace(self, s, radius, rate, rw, R, outer_hydraulic, inner_c, outer_mechanics):
        """Transforms of p(r), J(r) and h, for complex s of shape (n,)."""
        m, Ss, B = self.mobility, self.specific_storage, self.thickness
        beta = self.coupling / Ss
        a2 = self._half_alpha_comp
        lam_ps, G = self._lame_ps, self.shear
        q = np.sqrt(s * Ss / m)[:, None]
        rho = q.real

        # Scaled bases A'*I0s + B'*K0s keep every exponential factor <= 1.
        def i_s(order, r):
            return ive(order, q * r) * np.exp(rho * (r - R))

        def k_s(order, r):
            return kve(order, q * r) * np.exp(-q * (r - rw))

        def j_i(r):  # int_rw^r rho I0 / scaled
            return (r * i_s(1, r) - rw * i_s(1, rw)) / q

        def j_k(r):
            return -(r * k_s(1, r) - rw * k_s(1, rw)) / q

        n = s.size
        mat = np.zeros((n, 3, 3), dtype=complex)
        rhs = np.zeros((n, 3), dtype=complex)
        # Well face: 2 pi rw B m dp/dr = rate.
        mat[:, 0, 0] = m * q[:, 0] * i_s(1, rw)[:, 0]
        mat[:, 0, 1] = -m * q[:, 0] * k_s(1, rw)[:, 0]
        rhs[:, 0] = rate / (2 * math.pi * rw * B * s)
        if outer_hydraulic == 'fixed_pressure':
            mat[:, 1, 0] = i_s(0, R)[:, 0]
            mat[:, 1, 1] = k_s(0, R)[:, 0]
            mat[:, 1, 2] = -beta
        else:
            mat[:, 1, 0] = q[:, 0] * i_s(1, R)[:, 0]
            mat[:, 1, 1] = -q[:, 0] * k_s(1, R)[:, 0]
        area2 = (R**2 - rw**2) / 2
        if self.coupling == 0:
            mat[:, 2, 2] = 1.0  # uniaxial: h = 0
        else:
            # r U_r(R) = inner_c h + a2 J(R) + h (R^2-rw^2)/2, J(R) linear in A',B',h.
            ru = np.stack([a2 * j_i(R)[:, 0], a2 * j_k(R)[:, 0],
                           np.full(n, inner_c + area2 - a2 * beta * area2)], axis=1)
            if outer_mechanics == 'fixed':
                mat[:, 2, :] = ru
            else:  # (lambda*+2G) h R^2 = 2G r U_r(R)
                mat[:, 2, :] = -2 * G * ru
                mat[:, 2, 2] += (lam_ps + 2 * G) * R**2
        x = np.linalg.solve(mat, rhs[..., None])[..., 0]
        a, b, h = x[:, :1], x[:, 1:2], x[:, 2:]
        r = radius[None, :]
        p = a * i_s(0, r) + b * k_s(0, r) - beta * h
        J = a * j_i(r) + b * j_k(r) - beta * h * (r**2 - rw**2) / 2
        return np.concatenate([p, J, h], axis=1)


def _talbot(transform, t, terms):
    """Fixed Talbot inversion; transform maps complex s (n,) to (n, k)."""
    theta = np.arange(1, terms) * math.pi / terms
    r = 2 * terms / (5 * t)
    cot = 1 / np.tan(theta)
    s = np.concatenate([[r + 0j], r * theta * (cot + 1j)])
    sigma = theta + (theta * cot - 1) * cot
    weights = np.concatenate([[0.5 + 0j], np.exp(r * theta * (cot + 1j) * t - r * t) * (1 + 1j * sigma)])
    values = transform(s)
    return (r / terms) * np.exp(r * t) * np.real(weights[:, None] * values).sum(axis=0)


def bear_corapcioglu_1981b(radius_m, time_s, rate_m3_s, transmissivity, consolidation_coefficient,
                           thickness, well_radius=None):
    """Printed single-well solution, eqs. (99), (102), (104), (107).

    Terzaghi stress, incompressible fluid and grains (alpha'' >> n beta).
    Returns head drawdown, volumetric strain (negative in contraction),
    subsidence (positive down) and radial displacement (negative inward).
    well_radius=None uses the paper's approximation for eq. (107).
    """
    radius = np.asarray(radius_m, dtype=float)
    time = np.asarray(time_s, dtype=float)
    Q, T, Cv, B = rate_m3_s, transmissivity, consolidation_coefficient, thickness
    u = radius**2 / (4 * Cv * time)
    W = exp1(u)
    if well_radius is None:
        tail = 1.0
    else:
        uw = well_radius**2 / (4 * Cv * time)
        tail = np.exp(-uw) - uw * exp1(uw)
    return {
        'head_drawdown': Q / (4 * math.pi * T) * W,
        'volumetric_strain': -Q / (4 * math.pi * Cv * B) * W,
        'subsidence': Q / (8 * math.pi * Cv) * W,
        'radial_displacement': -Q * time / (4 * math.pi * B * radius) * (u * W - np.exp(-u) + tail),
    }
