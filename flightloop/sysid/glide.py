"""Glide drag polar CD = CD0 + k CL^2 from quasi-steady unpowered windows.

Two estimators on the same windows:
  energy : E = h + V^2/(2g), D = -W dE/dt / V, CL = W n_z / (q S)   (plan method; biased by vertical air motion)
  force  : CL, CD from body specific force rotated by alpha          (wind-insensitive; needs alpha)
Each window contributes its mean (CL^2, CD); OLS on window means; SE from the window regression.
"""
import numpy as np
from .common import air, deriv, smooth, windows, coeffs_from_accel, G


def _fit(cl2, cd):
    X = np.c_[np.ones_like(cl2), cl2]
    beta, *_ = np.linalg.lstsq(X, cd, rcond=None)
    r = cd - X @ beta
    dof = max(len(cd) - 2, 1)
    s2 = r @ r / dof
    cov = s2 * np.linalg.inv(X.T @ X)
    return {"CD0": float(beta[0]), "k": float(beta[1]), "se_CD0": float(np.sqrt(cov[0, 0])), "se_k": float(np.sqrt(cov[1, 1])),
            "resid_rms": float(np.sqrt(np.mean(r * r))), "n_windows": int(len(cd)), "dof": dof}


def estimate(df, ac, alpha=None, method="force", wins=None, rate=50.0):
    dt = 1.0 / rate
    wins = wins if wins is not None else windows(df, rate=rate)
    W = ac["mass"] * G
    rho, V, qd = air(df)
    if method == "energy":
        Vs = smooth(df["tas"].to_numpy(), 51, 2)
        hs = smooth(df["h_baro"].to_numpy(), 51, 2)
        E = hs + Vs * Vs / (2 * G)
        Ed = deriv(E, dt, 51, 2)
        D = -W * Ed / Vs
        CD = D / (qd * ac["S"])
        CL = W * (-df["az"].to_numpy() / G) / (qd * ac["S"])
    else:
        CL, CD, _ = coeffs_from_accel(df, ac, alpha)
    cl2 = np.array([np.mean(CL[a:b] ** 2) for a, b in wins])
    cd = np.array([np.mean(CD[a:b]) for a, b in wins])
    out = _fit(cl2, cd) if len(wins) >= 3 else {"CD0": np.nan, "k": np.nan, "n_windows": len(wins)}
    out["windows"] = wins
    out["points"] = {"cl2": cl2.tolist(), "cd": cd.tolist()}
    return out
