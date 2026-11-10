"""Unsaturated (Richards) flow references for the explicit DES formulation.

Contents:

- ``Retention``: van Genuchten-Mualem, Brooks-Corey-Burdine and Gardner
  exponential saturation and relative-permeability laws;
- ``content`` / ``content_slope`` / ``advance_content``: the nodal fluid
  content F(p) = phi*S(p) + s*p that the explicit update conserves, and the
  pointwise pressure for a given content increment;
- ``bishop_pressure``: the equivalent pore pressure chi*p for the mechanics;
- ``Column``: a 1D vertical prototype of the planned DES update (lumped P1
  storage, element mobility, conservative content update, local inversion);
- ``srivastava_yeh``: exact transient infiltration for the Gardner model
  (Srivastava and Yeh, 1991, Water Resour. Res. 27(5), 753-762);
- ``steady_pressure``: exact steady infiltration profile above a water table
  for any retention law;
- ``bishop_column``: exact drained settlement of an elastic, laterally
  confined column when its water table moves (verification level U8).

Conventions follow DES: z upward, pore pressure p positive in compression
(p = rho*g*psi for pressure head psi), suction is -p, and pressure enters
Darcy as q = -(k*kr/mu)*(grad p + rho*g*e_z). Retention ``alpha`` is in 1/Pa.
Infiltration flux is positive into the domain. The derivation and the
explicit stability argument are in doc/hydro/unsaturated.md.

Requires NumPy and SciPy for benchmark evaluation only.
"""

import math

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

MODELS = ('van_genuchten', 'brooks_corey', 'gardner')


class Retention:
    """Effective saturation Se(p), saturation S(p) and kr(p) for one material.

    van_genuchten: Se = (1 + (alpha*s)^n)^-m, m = 1 - 1/n, Mualem kr.
    brooks_corey:  Se = (alpha*s)^-n for alpha*s > 1 (air entry 1/alpha), else
                   1; Burdine kr = Se^(3 + 2/n), n being the pore-size index.
    gardner:       Se = kr = exp(-alpha*s).
    s = max(-p, 0) is the suction.
    """

    def __init__(self, model, alpha, n=None, sr=0.0):
        if model not in MODELS:
            raise ValueError(f'model must be one of {MODELS}')
        if not (alpha > 0 and math.isfinite(alpha)):
            raise ValueError('alpha must be finite and positive')
        if model == 'van_genuchten' and not (n is not None and n > 1):
            raise ValueError('van Genuchten requires n > 1')
        if model == 'brooks_corey' and not (n is not None and n > 0):
            raise ValueError('Brooks-Corey requires a positive pore-size index n')
        if not 0 <= sr < 1:
            raise ValueError('Require 0 <= residual saturation < 1')
        self.model, self.alpha, self.n, self.sr = model, alpha, n, sr
        self.m = 1 - 1 / n if model == 'van_genuchten' else None
        # Pressure above which the material is saturated.
        self.p_sat = -1 / alpha if model == 'brooks_corey' else 0.0

    def se(self, p):
        x = self.alpha * np.maximum(-np.asarray(p, float), 0.0)
        if self.model == 'van_genuchten':
            return (1 + x**self.n) ** -self.m
        if self.model == 'brooks_corey':
            return np.where(x > 1, np.maximum(x, 1.0) ** -self.n, 1.0)
        return np.exp(-x)

    def dse(self, p):
        """dSe/dp, nonnegative."""
        x = self.alpha * np.maximum(-np.asarray(p, float), 0.0)
        if self.model == 'van_genuchten':
            n, m = self.n, self.m
            return m * n * self.alpha * x**(n - 1) * (1 + x**n) ** (-m - 1)
        if self.model == 'brooks_corey':
            xs = np.maximum(x, 1.0)
            return np.where(x > 1, self.n * self.alpha * xs ** (-self.n - 1), 0.0)
        return np.where(x > 0, self.alpha * np.exp(-x), 0.0)

    def sat(self, p):
        return self.sr + (1 - self.sr) * self.se(p)

    def kr(self, p):
        se = self.se(p)
        if self.model == 'van_genuchten':
            # Se^(1/m) = 1/(1+y), y = (alpha*s)^n, so 1 - (1 - Se^(1/m))^m is
            # -expm1(-m*log1p(1/y)): accurate near saturation and when dry.
            y = (self.alpha * np.maximum(-np.asarray(p, float), 0.0)) ** self.n
            with np.errstate(divide='ignore'):
                inner = -np.expm1(-self.m * np.log1p(1 / y))
            return np.where(y > 0, np.sqrt(se) * inner**2, 1.0)
        if self.model == 'brooks_corey':
            return se ** (3 + 2 / self.n)
        return se

    def pressure_of_se(self, se):
        """Inverse of Se(p) on the unsaturated branch, for 0 < Se < 1."""
        se = np.asarray(se, float)
        if self.model == 'van_genuchten':
            return -np.expm1(-np.log(se) / self.m) ** (1 / self.n) / self.alpha
        if self.model == 'brooks_corey':
            return -se ** (-1 / self.n) / self.alpha
        return np.log(se) / self.alpha


def content(ret, phi, s, p):
    """Fluid volume per bulk volume, F(p) = phi*S(p) + s*p.

    s is the saturated fixed-stress storage (1/Pa). For p >= ret.p_sat this is
    phi + s*p, so the saturated update is unchanged.
    """
    return phi * ret.sat(p) + s * np.asarray(p, float)


def content_slope(ret, phi, s, p):
    """dF/dp = phi*(1 - Sr)*dSe/dp + s >= s."""
    return phi * (1 - ret.sr) * ret.dse(p) + s


def advance_content(ret, phi, s, p, dm, iterations=100):
    """Pointwise pressure p1 with F(p1) - F(p) = dm, to round-off.

    The increment form keeps the saturated update exact: when p and the trial
    p + dm/s are both saturated, the result is that trial, as in the existing
    DES update. Otherwise a safeguarded Newton iteration keeps a bracket on
    g(x) = phi*(S(x) - S(p)) + s*(x - p) - dm, which increases with x. Since S
    is monotone, x = p + dm/s and x = p bracket the root. With s = 0 the
    unsaturated inverse is closed-form. This is a local scalar root per node,
    not a global solve.
    """
    p = np.atleast_1d(np.asarray(p, float))
    dm = np.broadcast_to(np.asarray(dm, float), p.shape)
    s_old = ret.sat(p)
    target = (s_old + dm / phi - ret.sr) / (1 - ret.sr)  # Se reached when s = 0
    if s == 0:
        if np.any(target <= 0) or np.any(target > 1 + 1e-15):
            raise ValueError('Without storage, saturation must stay in (Sr, 1]')
        return np.where(target >= 1, ret.p_sat,
                        ret.pressure_of_se(np.minimum(target, 1.0)))
    trial = p + dm / s
    out = trial.copy()
    todo = ~((p >= ret.p_sat) & (trial >= ret.p_sat))
    if not np.any(todo):
        return out
    p0, d, s0 = p[todo], dm[todo], s_old[todo]
    lo, hi = np.minimum(p0, trial[todo]), np.maximum(p0, trial[todo])
    guess = ret.pressure_of_se(np.clip(target[todo], 1e-300, 1 - 1e-16))
    x = np.where((target[todo] > 0) & (target[todo] < 1), guess, 0.5 * (lo + hi))
    x = np.clip(x, lo, hi)
    for _ in range(iterations):
        g = phi * (ret.sat(x) - s0) + s * (x - p0) - d
        lo = np.where(g <= 0, x, lo)
        hi = np.where(g >= 0, x, hi)
        step = x - g / content_slope(ret, phi, s, x)
        new = np.where((step > lo) & (step < hi), step, 0.5 * (lo + hi))
        done = (new == x) | (hi - lo <= 4 * np.spacing(np.maximum(np.abs(lo), np.abs(hi))))
        x = new
        if np.all(done):
            break
    out[todo] = x
    return out


def bishop_pressure(ret, p):
    """Equivalent pore pressure pi = chi*p with Bishop chi = Se (Lu et al., 2010,
    doi:10.1029/2009WR008646). Effective stress is sigma' = sigma + alpha*pi,
    tension positive; pi = p when saturated. Mechanics and porosity take the
    increment of pi, a function of the nodal state, so the coupling cannot drift.
    """
    p = np.asarray(p, float)
    return np.where(p >= ret.p_sat, p, ret.se(p) * p)


class Column:
    """Vertical 1D prototype of the planned explicit DES unsaturated update.

    Nodes z[0] < ... < z[-1] with lumped volumes V_i (unit area). Element
    mobility is (k/mu) times the arithmetic mean of the nodal kr. Bottom
    pressure is prescribed; the top receives an infiltration flux. Each step

        dm_i = dt*(sum of Darcy fluxes into i + source_i)/V_i
        p_i <- p_i' with F(p_i') - F(p_i) = dm_i

    so the content change equals the boundary volume to round-off.
    """

    def __init__(self, z, ret, *, phi, storage, perm, mu=1e-3, rho=1000.0, g=9.81):
        self.z = np.asarray(z, float)
        if np.any(np.diff(self.z) <= 0):
            raise ValueError('Nodes must increase upward')
        self.ret, self.phi, self.s = ret, phi, storage
        self.mob, self.gamma = perm / mu, rho * g
        self.h = np.diff(self.z)
        self.vol = np.zeros(len(self.z))
        self.vol[:-1] += self.h / 2
        self.vol[1:] += self.h / 2

    def fluxes(self, p):
        """Net Darcy inflow at each node, element mobilities and upward fluxes."""
        kr = self.ret.kr(p)
        lam = self.mob * 0.5 * (kr[:-1] + kr[1:])
        q = -lam * (np.diff(p) + self.gamma * self.h) / self.h
        net = np.zeros_like(p)
        net[:-1] -= q
        net[1:] += q
        return net, lam, q

    def stable_dt(self, p, mode='saturated'):
        """Largest dt for which every free node's update is a convex combination.

        'saturated' uses kr <= 1 and dF/dp >= s, so it holds for any state and
        equals 0.5*h^2/D of the existing DES bound on a uniform mesh.
        'local' uses the current kr and slope and needs the chord check below.
        """
        if mode == 'saturated':
            lam, slope = np.full(len(self.h), self.mob), self.s
        else:
            lam, slope = self.fluxes(p)[1], content_slope(self.ret, self.phi, self.s, p)
        with np.errstate(divide='ignore'):
            return float(np.min((self.vol * slope / self.conductance(lam))[1:]))

    def conductance(self, lam):
        """Sum of the element conductances lam/h meeting at each node."""
        cond = np.zeros(len(self.z))
        cond[:-1] += lam / self.h
        cond[1:] += lam / self.h
        return cond

    def step(self, p, dt, q_top, p_bottom):
        net, lam, q = self.fluxes(p)
        net[-1] += q_top
        new = p.copy()
        new[1:] = advance_content(self.ret, self.phi, self.s, p[1:], dt * net[1:] / self.vol[1:])
        new[0] = p_bottom
        # Chord slopes dF/dp over the step; the update is monotone when the
        # local Courant number dt*conductance/(V*chord) stays at most 1.
        # Written as s + phi*dS/dp so that saturated nodes give exactly s.
        dp = new - p
        moved = dp != 0
        chord = np.where(moved, self.phi * (self.ret.sat(new) - self.ret.sat(p))
                         / np.where(moved, dp, 1), self.phi * (1 - self.ret.sr) * self.ret.dse(p))
        chord = self.s + np.maximum(chord, 0.0)
        with np.errstate(divide='ignore'):
            courant = float(np.max((dt * self.conductance(lam) / (self.vol * chord))[1:]))
        # The free nodes receive the first element's upward flux from the fixed node.
        return new, courant, q[0]

    def run(self, p, t_end, q_top, p_bottom, mode='saturated', safety=0.5, max_steps=10**7):
        """Advance to t_end. 'local' halves rejected steps (an explicit retry)."""
        p = np.asarray(p, float).copy()
        t, boundary, steps, retries, worst = 0.0, 0.0, 0, 0, 0.0
        m0 = float(np.dot(self.vol[1:], content(self.ret, self.phi, self.s, p[1:])))
        while t < t_end:
            dt = min(t_end - t, self.stable_dt(p, mode) * (1.0 if mode == 'saturated' else safety))
            while True:
                new, courant, bottom = self.step(p, dt, q_top, p_bottom)
                if courant <= 1 + 1e-12 or mode == 'saturated':
                    break
                dt *= 0.5
                retries += 1
            boundary += dt * (q_top + bottom)
            worst = max(worst, courant)
            p, t, steps = new, t + dt, steps + 1
            if steps > max_steps:
                raise RuntimeError('Too many steps')
        m1 = float(np.dot(self.vol[1:], content(self.ret, self.phi, self.s, p[1:])))
        return p, dict(steps=steps, retries=retries, courant=worst,
                       balance=m1 - m0 - boundary, stored=m1 - m0)


def srivastava_yeh(z, t, *, length, alpha, ks, dtheta, qa, qb, terms=200):
    """Exact Gardner-model infiltration above a water table at z = 0.

    alpha in 1/m of head, ks and fluxes in m/s (infiltration positive),
    dtheta = theta_s - theta_r. The column starts at steady flux qa and the top
    flux becomes qb at t = 0. Returns K/Ks; head is log(K/Ks)/alpha.
    Srivastava and Yeh (1991), homogeneous column. The series is checked
    against an independent matrix-exponential solution in the tests.
    """
    L = alpha * length
    a, b = qa / ks, qb / ks
    Z = alpha * np.asarray(z, float)
    T = alpha * ks * t / dtheta
    k = b - (b - 1) * np.exp(-Z)
    if t == 0:
        return a - (a - 1) * np.exp(-Z)
    total = np.zeros_like(Z)
    for lam in eigenvalues(L, terms):
        total += (np.sin(lam * Z) * math.sin(lam * L) * math.exp(-lam * lam * T)
                  / (1 + L / 2 + 2 * lam * lam * L))
    return k - 4 * (b - a) * np.exp((L - Z) / 2 - T / 4) * total


def eigenvalues(L, terms):
    """Positive roots of tan(lam*L) + 2*lam = 0, one in each ((n-1/2)pi, n*pi)/L."""
    f = lambda x: math.sin(x * L) + 2 * x * math.cos(x * L)
    return [brentq(f, (n - 0.5) * math.pi / L, n * math.pi / L, xtol=1e-15, rtol=1e-15)
            for n in range(1, terms + 1)]


def steady_pressure(ret, z, *, flux, perm, mu=1e-3, rho=1000.0, g=9.81):
    """Exact steady pressure at heights z >= 0 above a water table (p = 0).

    A downward flux q gives dp/dz = gamma*(q/K(p) - 1), K = Ks*kr(p). For
    0 < q < Ks, p decreases monotonically towards p* with K(p*) = q, so the
    right-hand side stays bounded; it is integrated to round-off tolerance.
    """
    gamma = rho * g
    ks = perm * gamma / mu
    if not 0 < flux < ks:
        raise ValueError('Require 0 < flux < Ks')
    z = np.asarray(z, float)
    if np.any(z < 0):
        raise ValueError('Heights must be above the water table')
    rhs = lambda _, p: [gamma * (flux / (ks * float(ret.kr(p[0]))) - 1)]
    top = float(np.max(z))
    if top == 0:
        return np.zeros_like(z)
    sol = solve_ivp(rhs, (0.0, top), [0.0], method='DOP853', rtol=1e-12, atol=1e-9,
                    dense_output=True, max_step=top / 200)
    return sol.sol(z)[0]


def bishop_column(ret, z, *, height, table0, table1, phi, biot, bulk, shear,
                  rho=1000.0, g=9.81, density=True):
    """Vertical displacement u(z) (up positive, u(0) = 0) of a laterally
    confined elastic column when a hydrostatic water table moves from height
    table0 to table1, both states drained and at rest.

    Uniaxial strain with Bishop effective stress sigma' = sigma + biot*pi,
    pi = Se*p, gives eps_zz = (d sigma_zz + biot*d pi)/M with M = K + 4G/3
    (plane strain or 3D). If density is true, the bulk density
    rho_s*(1-phi) + rho*phi*S changes with saturation, so that
    d sigma_zz(z) = -g*rho*phi*(1-Sr)*int_z^H dSe dz'. Exchanging the order of
    integration gives
        u(z) = (biot*int_0^z d pi dz' - g*rho*phi*(1-Sr)*int_0^H min(z', z)*dSe dz') / M.
    Pressure is p = -rho*g*(z - table); the column spans 0 <= z <= height.
    """
    gamma = rho * g
    modulus = bulk + 4 * shear / 3
    p = lambda x, t: -gamma * (x - t)
    dpi = lambda x: float(bishop_pressure(ret, p(x, table1)) - bishop_pressure(ret, p(x, table0)))
    dse = lambda x: float(ret.se(p(x, table1)) - ret.se(p(x, table0)))
    # Kinks: each water table, and the Brooks-Corey air-entry height above it.
    kinks = [t + d for t in (table0, table1) for d in (0.0, -ret.p_sat / gamma)]
    out = []
    for top in np.atleast_1d(np.asarray(z, float)):
        if not 0 <= top <= height:
            raise ValueError('Heights must lie inside the column')
        inside = lambda a, b: sorted({k for k in kinks if a < k < b}) or None
        kw = dict(epsabs=1e-14, epsrel=1e-12, limit=200)
        u = biot * quad(dpi, 0, top, points=inside(0, top), **kw)[0] if top > 0 else 0.0
        if density:
            weight = lambda x: min(x, top) * dse(x)
            pts = sorted(set(inside(0, height) or []) | ({top} if 0 < top < height else set())) or None
            u -= gamma * phi * (1 - ret.sr) * quad(weight, 0, height, points=pts, **kw)[0]
        out.append(u / modulus)
    return np.array(out)
