"""Analytical confined-aquifer drawdown, independent of the DES solver.

Theis (1935), doi:10.1029/TR016i002p00519; a USGS reprint is at
https://water.usgs.gov/ogw/pubs/Theis-1935.pdf .

Requires NumPy and SciPy for benchmark evaluation only. Rates are total well
discharges in m^3/s, positive for extraction and negative for injection.
"""

import math

import numpy as np
from scipy.special import exp1


def drawdown(radius_m, time_s, rate_m3_s, transmissivity_m2_s, storativity):
    """Return head drawdown in metres for a constant rate starting at t=0.

    Radius and time broadcast as NumPy arrays. The aquifer is infinite,
    homogeneous, confined and isotropic, with constant properties and a
    fully penetrating line well without wellbore storage. The well axis is
    singular and must be excluded; the initial drawdown is zero for r>0.
    """
    if not math.isfinite(rate_m3_s):
        raise ValueError('Well rate must be finite')
    for name, value in (('Transmissivity', transmissivity_m2_s),
                        ('Storativity', storativity)):
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f'{name} must be finite and positive')
    radius, time = np.broadcast_arrays(np.asarray(radius_m, dtype=float),
                                      np.asarray(time_s, dtype=float))
    if not np.isfinite(radius).all() or (radius <= 0).any():
        raise ValueError('Observation radii must be finite and positive')
    if not np.isfinite(time).all() or (time < 0).any():
        raise ValueError('Observation times must be finite and nonnegative')
    result = np.zeros(radius.shape, dtype=float)
    if rate_m3_s == 0:
        return result
    active = time > 0
    u = radius[active]**2 * storativity / (4 * transmissivity_m2_s * time[active])
    result[active] = rate_m3_s / (4 * math.pi * transmissivity_m2_s) * exp1(u)
    return result


def scheduled_drawdown(radius_m, time_s, schedule, transmissivity_m2_s,
                       storativity):
    """Superpose absolute rates given by ordered (start_seconds, rate) pairs.

    The rate before the first entry is zero. Entries remain active until the
    next entry; a zero rate represents shut-in. This linear reference applies
    only while material properties and the governing operator stay constant.
    """
    # Evaluate a zero-rate response to validate observations and properties,
    # including when the schedule is empty or starts after all observations.
    result = drawdown(radius_m, time_s, 0.0, transmissivity_m2_s, storativity)
    previous_time, previous_rate = -1.0, 0.0
    time = np.asarray(time_s, dtype=float)
    for start, rate in schedule:
        if not math.isfinite(start) or start < 0 or start <= previous_time:
            raise ValueError('Schedule times must be finite, nonnegative and strictly increasing')
        if not math.isfinite(rate):
            raise ValueError('Scheduled rates must be finite')
        lag = np.maximum(time - start, 0.0)
        result += drawdown(radius_m, lag, rate - previous_rate,
                           transmissivity_m2_s, storativity)
        previous_time, previous_rate = start, rate
    return result
