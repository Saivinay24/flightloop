"""Turbulence with MIL-F-8785C low-altitude Dryden intensities and length scales (h <= 1000 ft).

sigma_w = 0.1 W20, sigma_u = sigma_w / (0.177 + 0.000823 h)^0.4, L_w = h, L_u = h / (0.177 + 0.000823 h)^1.2
(h in ft, 'light' W20 = 15 kt). Both components are realised as first-order Gauss-Markov processes with time
constant L/V (exact discretisation, frozen field). This matches the Dryden u spectrum exactly and approximates
the second-order Dryden w spectrum with a first-order one (same variance and integral scale); stated wherever
results are shown.
"""
import numpy as np

FT = 0.3048


def params(h_m, level="light"):
    W20 = {"light": 15.0, "moderate": 30.0, "severe": 45.0}[level] * 0.514444
    h = np.clip(np.asarray(h_m) / FT, 10.0, 1000.0)
    sw = 0.1 * W20
    su = sw / (0.177 + 0.000823 * h) ** 0.4
    return su, sw, h / (0.177 + 0.000823 * h) ** 1.2 * FT, h * FT


class Dryden:
    def __init__(self, n, dt, rng, level="light", vertical=True, on=True):
        self.n, self.dt, self.rng, self.level, self.vertical, self.on = n, dt, rng, level, vertical, on
        self.u = np.zeros(n); self.w = np.zeros(n)

    def step(self, V, h):
        if not self.on:
            return self.u, self.w
        su, sw, Lu, Lw = params(h, self.level)
        a = np.exp(-V * self.dt / Lu)
        self.u = a * self.u + su * np.sqrt(1 - a * a) * self.rng.standard_normal(self.n)
        if self.vertical:
            b = np.exp(-V * self.dt / Lw)
            self.w = b * self.w + sw * np.sqrt(1 - b * b) * self.rng.standard_normal(self.n)
        return self.u, self.w
