"""Small-strain Biot consolidation of a laterally confined elastic column.

Depth is measured down from the drained top. The bottom is impermeable and
fixed vertically. A positive compressive load is applied instantaneously.
The reference describes t=0+; an unloaded numerical frame at t=0 is different.
See the PyLith Terzaghi solution cited in this directory's README.md.
"""

import numpy as np


class TerzaghiReference:
    def __init__(self, *, height, load, bulk, shear, porosity, fluid_bulk,
                 grain_bulk, biot, permeability, viscosity):
        values = (height, load, bulk, shear, fluid_bulk, grain_bulk,
                  permeability, viscosity)
        if not all(np.isfinite(x) and x > 0 for x in values):
            raise ValueError('Geometry, load, moduli and transport properties must be positive')
        if not 0 < porosity <= biot <= 1:
            raise ValueError('Require 0 < porosity <= Biot coefficient <= 1')
        self.height, self.load, self.biot = height, load, biot
        self.constrained_modulus = bulk + 4 * shear / 3
        self.storage = porosity / fluid_bulk + (biot - porosity) / grain_bulk
        self.effective_storage = self.storage + biot**2 / self.constrained_modulus
        self.diffusivity = permeability / viscosity / self.effective_storage
        self.initial_pressure = biot * load / self.constrained_modulus / self.effective_storage
        self.drained_settlement = load * height / self.constrained_modulus

    def _decay(self, time_s, terms):
        time = np.atleast_1d(np.asarray(time_s, dtype=float))
        if time.ndim != 1 or not np.isfinite(time).all() or (time < 0).any():
            raise ValueError('Times must be a finite nonnegative scalar or 1D array')
        if not isinstance(terms, int) or terms <= 0:
            raise ValueError('Series term count must be a positive integer')
        odd = 2 * np.arange(terms) + 1
        decay = np.exp(-time[:, None] * self.diffusivity *
                       (odd[None, :] * np.pi / (2 * self.height))**2)
        return time, odd, decay

    def pressure(self, depth_m, time_s, terms=512):
        """Return excess pressure in Pa, with shape (number of times, depths)."""
        depth = np.atleast_1d(np.asarray(depth_m, dtype=float))
        if depth.ndim != 1 or not np.isfinite(depth).all() or (depth < 0).any() or (depth > self.height).any():
            raise ValueError('Depths must be a finite scalar or 1D array within the column')
        time, odd, decay = self._decay(time_s, terms)
        shape = np.sin(depth[:, None] * odd[None, :] * np.pi / (2 * self.height))
        result = self.initial_pressure * (decay * (4 / (odd * np.pi))) @ shape.T
        # Use the exact undrained initial condition rather than a truncated
        # Fourier approximation of the discontinuity at the drained boundary.
        result[time == 0] = np.where(depth == 0, 0.0, self.initial_pressure)
        return result

    def settlement(self, time_s, terms=512):
        """Return positive downward top settlement in metres, including t=0+."""
        time, odd, decay = self._decay(time_s, terms)
        mean_pressure = self.initial_pressure * np.sum(decay * 8 / (odd * np.pi)**2, axis=1)
        mean_pressure[time == 0] = self.initial_pressure
        return (self.drained_settlement - self.biot * self.height * mean_pressure /
                self.constrained_modulus)

    def load_history(self, depth_m, time_s, events, terms=512):
        """Pressure and settlement for signed (time, load increment) pairs.

        A load increment is compression-positive and measured in Pa. Its
        response starts at t=event+, including the undrained jump. A saved
        checkpoint just before a load change belongs to the earlier history.
        """
        time = np.atleast_1d(np.asarray(time_s, dtype=float))
        # Reuse the single-step validation, including depth and series order.
        pressure = np.zeros_like(self.pressure(depth_m, time, terms))
        settlement = np.zeros(len(time))
        changes = np.asarray(events, dtype=float)
        if not changes.size:
            return pressure, settlement
        if changes.ndim != 2 or changes.shape[1] != 2 or not np.isfinite(changes).all():
            raise ValueError('Load events must be finite (time, increment) pairs')
        if (changes[:, 0] < 0).any() or (np.diff(changes[:, 0]) < 0).any():
            raise ValueError('Load-event times must be nonnegative and nondecreasing')
        for start, increment in changes:
            active = time >= start
            if not active.any():
                continue
            elapsed = time[active]-start
            pressure[active] += increment/self.load*self.pressure(depth_m, elapsed, terms)
            settlement[active] += increment/self.load*self.settlement(elapsed, terms)
        return pressure, settlement
