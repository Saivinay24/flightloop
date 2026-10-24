"""Diagnostic: is the measured-vs-model pitch acceleration gap a logging-rate effect?"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT))
import numpy as np
from flightloop.sim import truth as TR, maneuvers as M, sixdof as S
from flightloop.atmos import isa
p = TR.load()
for lr in (200, 50):
    rec = M.fly(p, 1, 1, [17.0], alt0=900, turb=False, log_rate=lr)
    dt = 1 / lr
    a = rec["alpha"][:, 0]; q = rec["q"][:, 0]; V = rec["tas"][:, 0]; de = rec["de"][:, 0]; h = rec["h"][:, 0]
    rho = isa(h)["rho"]; qd = 0.5 * rho * V * V
    z = 0 * a
    Cm = S.coeffs(p, a, z, z, q * p["c"] / (2 * V), z, np.degrees(de), z, z, V / 340)[4]
    qdot_model = qd * p["S"] * p["c"] * Cm / p["Iyy"]
    qdot_num = np.gradient(q, dt)
    m = slice(10, -10)
    print(lr, "Hz: slope gradient(q) vs model", np.polyfit(qdot_model[m], qdot_num[m], 1)[0],
          "rms diff", np.sqrt(np.mean((qdot_model - qdot_num)[m] ** 2)), "rms model", np.sqrt(np.mean(qdot_model[m] ** 2)))
    if lr == 200:
        fwd = (q[1:] - q[:-1]) / dt
        avg = 0.5 * (qdot_model[1:] + qdot_model[:-1])
        print("forward diff vs trapezoid model slope", np.polyfit(avg[10:-10], fwd[10:-10], 1)[0])
        i0 = 5396
        for i in range(i0, i0 + 12):
            print(i, round(rec["t"][i, 0], 3), "de", round(np.degrees(de[i]), 3), "q", round(q[i], 5), "fwd", round(fwd[i], 4), "model", round(qdot_model[i], 4), "alpha", round(np.degrees(a[i]), 3))
