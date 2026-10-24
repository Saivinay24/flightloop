"""Equation-error regression for CL and Cm on [1, alpha, q c/(2V), de].

Cm = Iyy q_dot / (q_bar S c). Consistent filtering (Klein and Morelli 2006, sec. 5.6): q, alpha, q_hat and de are
passed through the same zero-phase 4th-order Butterworth low-pass (lp_hz, default 5 Hz = about 2x the highest
closed-loop pitch pole), then q_dot is the 4th-order central difference of the filtered q (a 2nd-order difference reads 6.5% low at 5 Hz). Because the filter is linear and the
model is linear, filtering both sides keeps the regression unbiased while removing the noise that differentiation
amplifies. (A Savitzky-Golay q_dot on unfiltered regressors under-read Cm by 15 to 23% on noise-free data because the
elevator steps put power above the derivative filter's passband: scripts/diag/sysid_clean.py.)
CL from specific force and alpha. Standard errors with the Morelli-Klein correction for coloured residuals:
Cov = (X'X)^-1 [sum_k R_rr(k) sum_i x_i x_{i+k}'] (X'X)^-1, lags |k| <= L.
"""
import numpy as np
from .common import air, deriv, coeffs_from_accel


def ols(X, y, L=150, Sigma=None):
    """Sigma: regressor measurement-noise covariance (p x p). When given, bias-compensated least squares
    b = (X'X - n Sigma)^-1 X'y (errors-in-variables correction, e.g. Soderstrom 2018, ch. 5)."""
    A = X.T @ X - (len(y) * Sigma if Sigma is not None else 0.0)
    XtX_inv = np.linalg.inv(A)
    b = XtX_inv @ X.T @ y
    r = y - X @ b
    n = len(y)
    R = np.array([np.dot(r[: n - k], r[k:]) / n for k in range(L + 1)])
    S = R[0] * (X.T @ X)
    for k in range(1, L + 1):
        A = X[: n - k].T @ X[k:]
        S += R[k] * (A + A.T)
    cov_mk = XtX_inv @ S @ XtX_inv
    cov_ols = (r @ r / (n - X.shape[1])) * XtX_inv
    return b, np.sqrt(np.diag(cov_mk)), np.sqrt(np.diag(cov_ols)), float(np.sqrt(np.mean(r * r))), r


def regressors(df, ac, alpha, rate=50.0):
    rho, V, qd = air(df)
    qhat = df["q"].to_numpy() * ac["c"] / (2 * V)
    return np.c_[np.ones(len(df)), alpha, qhat, df["de"].to_numpy()], qd


def _lp(x, rate, hz):
    from scipy.signal import butter, filtfilt
    b, a = butter(4, hz / (0.5 * rate))
    return filtfilt(b, a, x, axis=0, padlen=min(len(x) - 1, 3 * max(len(a), len(b))))


def vane_noise_var(alpha):
    """White-noise variance of the alpha vane from the data: var(second difference) = 6 sigma^2 for white noise;
    robust (MAD) estimate. The true alpha (gusts and manoeuvres) has negligible second-difference power at 50 Hz."""
    d2 = np.diff(alpha, 2)
    s = 1.4826 * np.median(np.abs(d2 - np.median(d2)))
    return float(s * s / 6.0)


def lp_noise_gain(n, rate, hz, seed=0):
    """Variance gain of the zero-phase low-pass on unit white noise (measured on a long synthetic sequence)."""
    w = np.random.default_rng(seed).standard_normal(max(n, 20000))
    return float(np.var(_lp(w, rate, hz)[500:-500]))


def d4(x, dt):
    """4th-order central difference (gain error (w dt)^4/30, 0.5% at 5 Hz and 50 Hz) with 2nd-order ends."""
    y = np.gradient(x, dt)
    y[2:-2] = (-x[4:] + 8 * x[3:-1] - 8 * x[1:-3] + x[:-4]) / (12 * dt)
    return y


def estimate(df, ac, alpha, mask=None, rate=50.0, L=150, lp_hz=5.0, eiv=True):
    dt = 1.0 / rate
    X, qd = regressors(df, ac, alpha, rate)
    CL, CD, _ = coeffs_from_accel(df, ac, alpha)
    CLf = _lp(CL, rate, lp_hz)
    qf = _lp(df["q"].to_numpy(), rate, lp_hz)
    qdot = d4(qf, dt)
    Cm = ac["Iyy"] * qdot / (qd * ac["S"] * ac["c"])
    Xf = X.copy(); Xf[:, 1:] = _lp(X[:, 1:], rate, lp_hz)
    m = np.ones(len(df), bool) if mask is None else mask
    X, Xf, CL, CLf, Cm = X[m], Xf[m], CL[m], CLf[m], Cm[m]
    names = ["0", "a", "q", "de"]
    out = {}
    s2 = vane_noise_var(alpha) if eiv else 0.0
    Sig = np.zeros((4, 4)); Sig[1, 1] = s2 * lp_noise_gain(len(alpha), rate, lp_hz)
    out["vane_noise_std_est_rad"] = float(np.sqrt(s2))
    for lab, y, XX in (("CL", CLf, Xf), ("Cm", Cm, Xf)):
        b, se_mk, se_ols, rms, _ = ols(XX, y, L, Sig if eiv else None)
        for i, nmv in enumerate(names):
            out[f"{lab}{nmv}"] = float(b[i]); out[f"se_{lab}{nmv}"] = float(se_mk[i]); out[f"se_ols_{lab}{nmv}"] = float(se_ols[i])
        out[f"resid_rms_{lab}"] = rms
    Xs = X / X.std(0, where=np.r_[False, True, True, True]) if False else X
    c = np.corrcoef(X[:, 1:].T)
    out["regressor_corr"] = c.tolist()
    out["cond"] = float(np.linalg.cond(X[:, 1:] - X[:, 1:].mean(0)))
    return out
