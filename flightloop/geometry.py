"""Airfoil geometry: Selig .dat reader, arc-length cubic spline, sharp trailing edge, unit chord."""
from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline


def read_selig(path):
    rows = []
    for line in open(path).read().splitlines()[1:]:
        parts = line.split()
        if len(parts) >= 2:
            try:
                rows.append((float(parts[0]), float(parts[1])))
            except ValueError:
                pass
    xy = np.array(rows)
    return xy


def with_le_point(xy):
    """Insert the design leading edge (0, 0) between upper and lower surfaces when the file omits it
    (UIUC e387.dat jumps from (0.00044, 0.00234) to (0.00091, -0.00286); TM-4062 Table I lists (0, 0))."""
    if np.any(np.all(np.isclose(xy, 0.0), axis=1)):
        return xy
    i = int(np.argmin(xy[:, 0]))
    k = i + 1 if xy[i, 1] > 0 else i          # insert after the last upper-surface point
    return np.vstack([xy[:k], [[0.0, 0.0]], xy[k:]])


class Airfoil:
    """Closed-contour spline from TE (upper) -> LE -> TE (lower), parametrised by normalised arc length."""

    def __init__(self, xy):
        xy = np.asarray(xy, float)
        # normalise chord to 1 with LE at origin
        i_le = np.argmin(xy[:, 0])
        xle, yle = xy[i_le]
        xte = 0.5 * (xy[0, 0] + xy[-1, 0])
        xy = (xy - [xle, yle]) / (xte - xle)
        # sharp TE: close the contour at the mean TE point
        te = 0.5 * (xy[0] + xy[-1])
        xy[0], xy[-1] = te, te
        d = np.r_[0, np.cumsum(np.hypot(*np.diff(xy, axis=0).T))]
        self.L = d[-1]
        self.t = d / d[-1]
        self.sx = CubicSpline(self.t, xy[:, 0])
        self.sy = CubicSpline(self.t, xy[:, 1])
        tt = np.linspace(0, 1, 20001)
        self.t_le = tt[np.argmin(self.sx(tt))]
        self.xy_raw = xy

    def point(self, t):
        return np.stack([self.sx(t), self.sy(t)], -1)

    def upper_t(self, x):
        """parameter t on the upper surface (t < t_le) for chordwise x."""
        tt = np.linspace(0, self.t_le, 20001)
        xx = self.sx(tt)
        return np.interp(x, xx[::-1], tt[::-1])

    def lower_t(self, x):
        tt = np.linspace(self.t_le, 1, 20001)
        return np.interp(x, self.sx(tt), tt)

    def thickness_camber(self, n=400):
        x = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n)))
        yu = self.sy(self.upper_t(x))
        yl = self.sy(self.lower_t(x))
        t = yu - yl
        c = 0.5 * (yu + yl)
        return {"t_max": float(t.max()), "x_t": float(x[t.argmax()]), "camber_max": float(c.max()), "x_c": float(x[c.argmax()])}
