"""Textbook supersonic theory (gamma = 1.4): oblique shock (theta-beta-M, weak branch), Prandtl-Meyer expansion,
shock-expansion and Ackeret (linear) for a symmetric diamond airfoil. Academic verification geometry only."""
import numpy as np
from scipy.optimize import brentq

GAM = 1.4


def theta_from_beta(M, b):
    return np.arctan(2 / np.tan(b) * (M * M * np.sin(b) ** 2 - 1) / (M * M * (GAM + np.cos(2 * b)) + 2))


def oblique(M, theta):
    """weak oblique shock: beta, p2/p1, M2 (theta in rad)."""
    mu = np.arcsin(1 / M)
    bs = np.linspace(mu + 1e-9, np.pi / 2, 20001)
    th = theta_from_beta(M, bs)
    i_max = int(np.argmax(th))
    if theta > th[i_max]:
        raise ValueError("detached")
    b = brentq(lambda x: theta_from_beta(M, x) - theta, mu + 1e-12, bs[i_max])
    Mn1 = M * np.sin(b)
    p21 = 1 + 2 * GAM / (GAM + 1) * (Mn1 ** 2 - 1)
    Mn2 = np.sqrt((1 + 0.5 * (GAM - 1) * Mn1 ** 2) / (GAM * Mn1 ** 2 - 0.5 * (GAM - 1)))
    M2 = Mn2 / np.sin(b - theta)
    return b, p21, M2


def nu(M):
    g = GAM
    return np.sqrt((g + 1) / (g - 1)) * np.arctan(np.sqrt((g - 1) / (g + 1) * (M * M - 1))) - np.arctan(np.sqrt(M * M - 1))


def M_from_nu(v):
    return brentq(lambda M: nu(M) - v, 1.0 + 1e-12, 50.0)


def p_ratio_isentropic(M1, M2):
    g = GAM
    return ((1 + 0.5 * (g - 1) * M1 ** 2) / (1 + 0.5 * (g - 1) * M2 ** 2)) ** (g / (g - 1))


def surface_p(M, p1, turn):
    """pressure after turning the flow by `turn` rad (positive = compression)."""
    if turn > 0:
        b, p21, M2 = oblique(M, turn)
        return p1 * p21, M2
    if turn < 0:
        M2 = M_from_nu(nu(M) - turn)
        return p1 * p_ratio_isentropic(M, M2), M2
    return p1, M


def diamond_shock_expansion(M, tc, alpha):
    """Symmetric diamond, max thickness at mid-chord, unit chord, freestream p = 1. Returns Cl, Cd."""
    eps = np.arctan(tc)          # half-angle: tan(eps) = (t/2)/(c/2)
    q = 0.5 * GAM * M * M        # with p_inf = 1
    # upper front: turn = eps - alpha ; upper rear: expansion by 2 eps from upper front
    pu1, Mu1 = surface_p(M, 1.0, eps - alpha)
    pu2, _ = surface_p(Mu1, pu1, -2 * eps)
    pl1, Ml1 = surface_p(M, 1.0, eps + alpha)
    pl2, _ = surface_p(Ml1, pl1, -2 * eps)
    L = 0.5 / np.cos(eps)        # panel length
    # panel outward normals (body axes, chord along x): upper front (-sin e, cos e), upper rear (sin e, cos e)
    Fx = L * ((pu1 - pu2) * np.sin(eps) + (pl1 - pl2) * np.sin(eps))
    Fy = L * (-(pu1 + pu2) * np.cos(eps) + (pl1 + pl2) * np.cos(eps))
    ca, sa = np.cos(alpha), np.sin(alpha)
    D = Fx * ca + Fy * sa
    Lf = Fy * ca - Fx * sa
    return Lf / q, D / q


def diamond_ackeret(M, tc, alpha):
    B = np.sqrt(M * M - 1)
    return 4 * alpha / B, 4 * (alpha ** 2 + tc ** 2) / B   # t/c diamond: Cd_thickness = 4 (t/c)^2 / B


def table():
    d = np.radians
    out = {"wedge": {}, "diamond": {}}
    for M in (2.0, 3.0):
        b, p21, M2 = oblique(M, d(10))
        out["wedge"][f"M{M:g}"] = {"beta_deg": float(np.degrees(b)), "p2_p1": float(p21), "M2": float(M2)}
        for a in (0, 2, 4):
            cl, cd = diamond_shock_expansion(M, 0.05, d(a))
            cla, cda = diamond_ackeret(M, 0.05, d(a))
            out["diamond"][f"M{M:g}_a{a}"] = {"cl": float(cl), "cd": float(cd), "cl_ackeret": float(cla), "cd_ackeret": float(cda)}
    return out


if __name__ == "__main__":
    import json
    print(json.dumps(table(), indent=1))
