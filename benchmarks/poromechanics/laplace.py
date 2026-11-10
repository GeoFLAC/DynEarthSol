"""Gaver–Stehfest weights shared by the independent poroelastic references."""
import math

import numpy as np


def stehfest_weights(n=14):
    if n < 2 or n % 2:
        raise ValueError('Stehfest requires an even number of terms >= 2')
    weights = np.zeros(n)
    half = n//2
    for k in range(1, n + 1):
        total = 0.0
        for j in range((k + 1)//2, min(k, half) + 1):
            total += (j**half*math.factorial(2*j)/(math.factorial(half - j)*math.factorial(j)
                      * math.factorial(j - 1)*math.factorial(k - j)*math.factorial(2*j - k)))
        weights[k - 1] = (-1)**(k + half)*total
    return weights


def stehfest(transform, t, n=14):
    """Invert a smooth, non-oscillating transform, retaining term order."""
    total = 0.0
    for k, weight in enumerate(stehfest_weights(n), 1):
        total += weight*transform(k*math.log(2)/t)
    return total*math.log(2)/t
