"""P0 smoke tests summary: SU2 INC NACA0012 (tutorial) and SU2 LM T3A flat plate (tutorial).

T3A: Re_x from the solution (rho U x / mu, freestream values at the inlet edge). ERCOFTAC data is not shipped with
the SU2 tutorial, so the plot compares with Blasius laminar 0.664/sqrt(Re_x) and turbulent 0.0592 Re_x^-0.2 only;
transition onset = Re_x of minimum Cf.
"""
import glob, json, sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from flightloop.cfd.vtu import read_vtu

out = {}
d = sorted(glob.glob(str(ROOT / "runs/p0-smoke-t3a-*")))[-1]
s = read_vtu(pathlib.Path(d) / "surface.vtu")
P = s["Points"]; o = np.argsort(P[:, 0]); x = P[o, 0]
cf = s["Skin_Friction_Coefficient"][o, 0]
rho = s["Density"][o]; mu = s["Laminar_Viscosity"][o]
v = pathlib.Path(d) / "vol_solution.vtu"
V = read_vtu(v)
Uinf = float(np.max(np.linalg.norm(V["Velocity"][:, :2], axis=1)))
rex = rho * Uinf * x / mu
m = (x > 0.02) & (x < 19)
i_min = np.argmin(np.where(m & (rex > 2e4), cf, 1))
h = pd.read_csv(pathlib.Path(d) / "history.csv"); h.columns = [c.strip().strip('"') for c in h.columns]
out["t3a"] = {"run": pathlib.Path(d).name, "iters": int(len(h)), "U_inf": Uinf, "onset_rex_min_cf": float(rex[i_min]),
              "cf_min": float(cf[i_min]), "cf_at_rex_5e5": float(np.interp(5e5, rex[m], cf[m])),
              "turb_corr_at_rex_5e5": 0.0592 * 5e5 ** -0.2, "lam_ratio_at_rex_5e4": float(np.interp(5e4, rex[m], cf[m]) / (0.664 / np.sqrt(5e4)))}
fig, ax = plt.subplots(figsize=(6, 4))
ax.loglog(rex[m], cf[m], "k-", label="SU2 SST-LM (T3A tutorial grid)")
r = np.logspace(3.5, 6.5, 100)
ax.loglog(r, 0.664 / np.sqrt(r), "b--", label="Blasius laminar")
ax.loglog(r, 0.0592 * r ** -0.2, "r--", label="turbulent 0.0592 Re_x^-0.2")
ax.set_xlabel("Re_x"); ax.set_ylabel("Cf"); ax.legend(); ax.set_title("P0 smoke: T3A flat plate (CFD vs correlations)")
fig.savefig(ROOT / "results/figures/p0_t3a_cf.png", dpi=130, bbox_inches="tight")
d2 = sorted(glob.glob(str(ROOT / "runs/p0-smoke-naca0012_inc-*")))[-1]
# the tutorial cfg writes aero coefficients only to the screen: parse the log table (iter, t_iter, wall, rmsP, rmsNu, CL, CD)
rowsl = [l for l in (pathlib.Path(d2) / "log.txt").read_text().splitlines() if l.startswith("|") and l.count("|") == 8]
vals = [[float(x) for x in l.strip("|").split("|")] for l in rowsl if l.strip("|").split("|")[0].strip().isdigit()]
vals = np.array(vals)
out["naca0012_inc"] = {"run": pathlib.Path(d2).name, "iters": int(vals[-1, 0]) + 1, "alpha_deg": 10.0, "re": 6e6,
                       "cl": float(vals[-1, 5]), "cd": float(vals[-1, 6]),
                       "rms_p_drop_orders": float(vals[:20, 3].max() - vals[-1, 3]),
                       "note": "tutorial case run to its 4000-iteration cap in the smoke test"}
(ROOT / "results/p0_smoke.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
