"""US Standard Atmosphere 1976, 0-47 km, geometric altitude input (converted to geopotential), Sutherland mu."""
import numpy as np

R0, G0, R, GAMMA = 6356766.0, 9.80665, 287.05287, 1.4
LAYERS = [(0.0, 288.15, -0.0065), (11000.0, 216.65, 0.0), (20000.0, 216.65, 0.001), (32000.0, 228.65, 0.0028), (47000.0, 270.65, 0.0)]


def _pbase():
    p = [101325.0]
    for i in range(len(LAYERS) - 1):
        hb, tb, L = LAYERS[i]
        ht = LAYERS[i + 1][0]
        p.append(p[-1] * (np.exp(-G0 * (ht - hb) / (R * tb)) if L == 0 else (tb / (tb + L * (ht - hb))) ** (G0 / (R * L))))
    return p


PB = _pbase()


def isa(z):
    z = np.asarray(z, float)
    H = R0 * z / (R0 + z)
    T = np.empty_like(H); p = np.empty_like(H)
    for i, (hb, tb, L) in enumerate(LAYERS[:-1]):
        m = (H >= hb) & (H < LAYERS[i + 1][0]) if i > 0 else (H < LAYERS[1][0])
        T[m] = tb + L * (H[m] - hb)
        p[m] = PB[i] * (np.exp(-G0 * (H[m] - hb) / (R * tb)) if L == 0 else (tb / T[m]) ** (G0 / (R * L)))
    rho = p / (R * T)
    mu = 1.458e-6 * T ** 1.5 / (T + 110.4)
    return {"T": T, "p": p, "rho": rho, "mu": mu, "nu": mu / rho, "a": np.sqrt(GAMMA * R * T)}
