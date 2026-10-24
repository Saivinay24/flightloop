"""Output-error short-period identification on excitation segments (3-2-1-1, 2-1-1, doublet).

Model per segment (perturbations about the segment start, q_bar and V held at segment means):
  alpha_dot = Za alpha + q + Zd de,   q_dot = Ma alpha + Mq q + Md de
  Za = -qS CLa/(mV), Zd = -qS CLde/(mV), Ma = qSc Cma/Iyy, Mq = qSc^2 Cmq/(2 V Iyy), Md = qSc Cmde/Iyy
Shared aero parameters [CLa, Cma, Cmq, Cmde] (+CLde fixed at the equation-error value) and per-segment initial
states/offsets. Fit by scipy least_squares on measured (alpha, q); Cramer-Rao bounds = s^2 (J'J)^-1.
Ignores process noise (turbulence), so it is biased when gusts are strong: stated, not hidden.
"""
import numpy as np
from scipy.optimize import least_squares
from .common import air, G


def simulate(th, seg, ac, CLde, dt):
    """Exact ZOH discretisation; input response by IIR filtering (lfilter), initial-state response by matrix powers."""
    from scipy.linalg import expm
    from scipy.signal import lfilter, ss2tf
    CLa, Cma, Cmq, Cmde = th[:4]
    a0, q0, ab, qb = th[4:8]
    qS, V = seg["qS"], seg["V"]
    m, Iyy, c = ac["mass"], ac["Iyy"], ac["c"]
    Za, Zd = -qS * CLa / (m * V), -qS * CLde / (m * V)
    Ma, Mq, Md = qS * c * Cma / Iyy, qS * c * c * Cmq / (2 * V * Iyy), qS * c * Cmde / Iyy
    A = np.array([[Za, 1.0], [Ma, Mq]]); B = np.array([[Zd], [Md]])
    Mx = expm(np.block([[A, B], [np.zeros((1, 3))]]) * dt)
    Ad, Bd = Mx[:2, :2], Mx[:2, 2:]
    de = seg["de"] - seg["de"][0]
    n = len(de)
    num, den = ss2tf(Ad, Bd, np.eye(2), np.zeros((2, 1)))
    ya = lfilter(num[0], den, de)
    yq = lfilter(num[1], den, de)
    w, Vv = np.linalg.eig(Ad)
    k = np.arange(n)
    coef = np.linalg.solve(Vv, np.array([a0, q0], complex))
    free = (Vv[:, :, None] * (coef[:, None] * w[:, None] ** k[None, :])[None, :, :]).sum(1).real
    return ya + free[0] + ab, yq + free[1] + qb


def estimate(df, ac, alpha, segs, CLde, theta0, rate=50.0):
    dt = 1.0 / rate
    rho, V, qd = air(df)
    S = []
    for a, b in segs:
        S.append({"qS": float(np.mean(qd[a:b])) * ac["S"], "V": float(np.mean(V[a:b])), "de": df["de"].to_numpy()[a:b],
                  "alpha": alpha[a:b] - alpha[a], "q": df["q"].to_numpy()[a:b] - df["q"].to_numpy()[a]})
    sa, sq = 0.01, 0.01

    def unpack(th):
        return [np.r_[th[:4], th[4 + 4 * i: 8 + 4 * i]] for i in range(len(S))]

    def res(th):
        out = []
        for p, s in zip(unpack(th), S):
            am, qm = simulate(p, s, ac, CLde, dt)
            out += [(am - s["alpha"]) / sa, (qm - s["q"]) / sq]
        return np.concatenate(out)

    x0 = np.r_[theta0, np.zeros(4 * len(S))]
    sol = least_squares(res, x0, method="lm", max_nfev=4000)
    r = sol.fun
    J = sol.jac
    dof = max(len(r) - len(x0), 1)
    s2 = r @ r / dof
    cov = s2 * np.linalg.pinv(J.T @ J)
    se = np.sqrt(np.diag(cov))[:4]
    d = np.sqrt(np.diag(cov))
    corr = cov[:4, :4] / np.outer(d[:4], d[:4])
    names = ["CLa", "Cma", "Cmq", "Cmde"]
    out = {n: float(sol.x[i]) for i, n in enumerate(names)}
    out.update({f"crb_{n}": float(se[i]) for i, n in enumerate(names)})
    out["corr"] = corr.tolist()
    out["resid_rms_alpha_rad"] = float(np.sqrt(np.mean((r[: len(r) // 2] * sa) ** 2)))
    out["n_segments"] = len(S)
    out["success"] = bool(sol.success)
    return out
