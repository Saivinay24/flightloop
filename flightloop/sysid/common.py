"""Shared helpers: air data, aero coefficients from logs, smoothing derivatives, quasi-steady windows."""
import numpy as np
from scipy.signal import savgol_filter
from ..atmos import isa

G = 9.80665


def deriv(x, dt, win=11, poly=3):
    return savgol_filter(x, win, poly, deriv=1, delta=dt)


def smooth(x, win=11, poly=3):
    return savgol_filter(x, win, poly)


def air(df):
    rho = isa(df["h_baro"].to_numpy())["rho"]
    V = df["tas"].to_numpy()
    return rho, V, 0.5 * rho * V * V


def coeffs_from_accel(df, ac, alpha):
    """Aerodynamic force coefficients from specific force (zero thrust) and alpha (rad)."""
    rho, V, qd = air(df)
    qS = qd * ac["S"]
    CX = ac["mass"] * df["ax"].to_numpy() / qS
    CZ = ac["mass"] * df["az"].to_numpy() / qS
    ca, sa = np.cos(alpha), np.sin(alpha)
    CD = -(CX * ca + CZ * sa)
    CL = CX * sa - CZ * ca
    return CL, CD, qd


def windows(df, dVdt_max=0.15, q_max=0.03, min_len_s=4.0, rate=50.0, mask=None):
    """Contiguous quasi-steady windows: |dV/dt| and |q| below thresholds for at least min_len_s."""
    dt = 1.0 / rate
    V = smooth(df["tas"].to_numpy(), 51, 2)
    dV = deriv(V, dt, 51, 2)
    q = smooth(df["q"].to_numpy(), 25, 2)
    ok = (np.abs(dV) < dVdt_max) & (np.abs(q) < q_max)
    if mask is not None:
        ok &= mask
    out, start = [], None
    for i, o in enumerate(np.r_[ok, False]):
        if o and start is None:
            start = i
        elif not o and start is not None:
            if (i - start) * dt >= min_len_s:
                out.append((start, i))
            start = None
    return out


def no_vane_alpha(df, wind_ned=(0.0, 0.0, 0.0)):
    """alpha from EKF attitude + inertial velocity minus the (steady) wind estimate, rotated to body axes."""
    phi, th, psi = (df[c].to_numpy() for c in ("phi", "theta", "psi"))
    v = np.stack([df["vn"] - wind_ned[0], df["ve"] - wind_ned[1], df["vd"] - wind_ned[2]], 1)
    cph, sph, cth, sth, cps, sps = np.cos(phi), np.sin(phi), np.cos(th), np.sin(th), np.cos(psi), np.sin(psi)
    # R_bn (NED -> body)
    u = cth * cps * v[:, 0] + cth * sps * v[:, 1] - sth * v[:, 2]
    w = (cph * sth * cps + sph * sps) * v[:, 0] + (cph * sth * sps - sph * cps) * v[:, 1] + cph * cth * v[:, 2]
    return np.arctan2(w, u)
