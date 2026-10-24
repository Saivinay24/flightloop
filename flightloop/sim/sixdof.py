"""Synthetic-flight simulator: StratoGlide's own 6-DOF equations and aero (stratoglide.release.sim6dof.deriv,
imported read-only from the sister repo) driven by FlightLoop's RK4 loop, controller, turbulence and thrust.

Vectorised over N independent flights (seeds). State: pos NED, body velocity wrt ground, quaternion, body rates.
Thrust (optional, for the prop rule in module F) acts along body x and is added to StratoGlide's vdot.
"""
from __future__ import annotations

import sys
import numpy as np

from ..config import cfg
import pathlib

_SG = pathlib.Path(cfg()["stratoglide_root"]).expanduser() / "src"
if str(_SG) not in sys.path:
    sys.path.insert(0, str(_SG))
from stratoglide.release import sim6dof as SG  # noqa: E402
from stratoglide.release.aero_tables import coeffs  # noqa: E402,F401

G = 9.80665
q2R, euler2q, q2euler = SG.q2R, SG.euler2q, SG.q2euler


class GustWind:
    """Wind callable for StratoGlide's deriv: steady NED wind plus the current gust sample (held over a step)."""

    def __init__(self, n, steady_ned=(0.0, 0.0, 0.0)):
        self.steady = np.broadcast_to(np.asarray(steady_ned, float), (n, 3)).copy()
        self.gust = np.zeros((n, 3))

    def __call__(self, z):
        return self.steady + self.gust


def deriv(X, p, ctl, wind, thrust=None):
    k, info = SG.deriv(X, p, ctl, wind)
    if thrust is not None:
        k[:, 3] += thrust / np.broadcast_to(np.asarray(p["mass"], float), (len(X),))
    return k, info


def rk4(X, p, ctl, wind, dt, thrust=None):
    k1, info = deriv(X, p, ctl, wind, thrust)
    k2, _ = deriv(X + 0.5 * dt * k1, p, ctl, wind, thrust)
    k3, _ = deriv(X + 0.5 * dt * k2, p, ctl, wind, thrust)
    k4, _ = deriv(X + dt * k3, p, ctl, wind, thrust)
    info["qdot"] = k1[:, 11]; info["qdot_avg"] = (k1[:, 11] + 2 * k2[:, 11] + 2 * k3[:, 11] + k4[:, 11]) / 6
    Xn = X + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    Xn[:, 6:10] /= np.linalg.norm(Xn[:, 6:10], axis=1, keepdims=True)
    return Xn, info


def trim_glide(p, V, rho, thrust=0.0):
    """Wings-level trim at airspeed V: returns alpha, de_deg, gamma, theta (rad). Linear attached model."""
    DEG = np.pi / 180
    W = p["mass"] * G
    qS = 0.5 * rho * V * V * p["S"]
    gam = -0.08
    for _ in range(50):
        CL = (W * np.cos(gam) - thrust * np.sin(0.0)) / qS
        # Cm0 + Cma a + Cmde de = 0 ; CL0 + CLa a + CLde de = CL
        A = np.array([[p["CLa"], p["CLde"]], [p["Cma"], p["Cmde"]]])
        a, de = np.linalg.solve(A, [CL - p["CL0"], -p["Cm0"]])
        CD = p["CD0"] + p["k_ind"] * CL * CL
        gam = np.arcsin(np.clip((thrust - CD * qS) / W, -1, 1))
    return float(a), float(de), float(gam), float(a + gam)


def initial_state(p, V, alt, rho, n, thrust=0.0):
    a, de, gam, th = trim_glide(p, V, rho, thrust)
    X = np.zeros((n, 13))
    X[:, 2] = -alt
    X[:, 3] = V * np.cos(a)
    X[:, 5] = V * np.sin(a)
    X[:, 6:10] = euler2q(np.zeros(n), np.full(n, th), np.zeros(n))
    return X, de
