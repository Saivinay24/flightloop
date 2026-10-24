"""Grid Convergence Index, Celik et al. 2008 (J. Fluids Eng. 130, 078001), procedure steps 1-5.

h_i = (1/N_i)^(1/2) for 2D. Apparent order from the fixed-point iteration
p = |ln|e32/e21| + q(p)| / ln r21,  q(p) = ln((r21^p - s)/(r32^p - s)), s = sign(e32/e21).
Oscillatory convergence (e32/e21 < 0) is reported as such with no GCI (never forced).
"""
import math


def celik(N1, N2, N3, f1, f2, f3, Fs=1.25):
    """1 = fine, 3 = coarse."""
    h1, h2, h3 = (1 / N1) ** 0.5, (1 / N2) ** 0.5, (1 / N3) ** 0.5
    r21, r32 = h2 / h1, h3 / h2
    e21, e32 = f2 - f1, f3 - f2
    out = {"r21": r21, "r32": r32, "f": [f1, f2, f3], "e21": e21, "e32": e32}
    if e21 == 0 or e32 == 0:
        out.update(type="exact-repeat", p=None, gci_fine_pct=None)
        return out
    ratio = e32 / e21
    if ratio < 0:
        out.update(type="oscillatory", p=None, gci_fine_pct=None,
                   oscillation_band_pct=100 * (max(f1, f2, f3) - min(f1, f2, f3)) / abs(f1))
        return out
    s = 1.0
    p = abs(math.log(abs(ratio))) / math.log(r21)
    for _ in range(200):
        try:
            q = math.log((r21 ** p - s) / (r32 ** p - s))
        except ValueError:
            break
        pn = abs(math.log(abs(ratio)) + q) / math.log(r21)
        if abs(pn - p) < 1e-10:
            p = pn
            break
        p = pn
    f_ext = (r21 ** p * f1 - f2) / (r21 ** p - 1)
    ea21 = abs((f1 - f2) / f1)
    gci = Fs * ea21 / (r21 ** p - 1)
    out.update(type="monotone", p=p, f_ext=f_ext, ea21_pct=100 * ea21, e_ext21_pct=100 * abs((f_ext - f1) / f_ext),
               gci_fine_pct=100 * gci)
    return out
