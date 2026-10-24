"""Diagnostic: one clean flight (no turbulence, no sensor noise), estimators on true signals."""
import sys, pathlib, time
ROOT = pathlib.Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT))
import numpy as np, pandas as pd
from flightloop.sim import truth as TR, maneuvers as M
from flightloop.sysid import glide, equation_error as EE
from flightloop.sysid.common import windows, coeffs_from_accel
p = TR.load(); T = TR.truth_params(p)
AC = {k: p[k] for k in ("S", "c", "b", "mass", "Iyy")}
turb = len(sys.argv) > 1
t0 = time.time()
rec = M.fly(p, 1, 1, [14.5, 17.0], alt0=900, turb=turb)
print("fly s", time.time() - t0, "T", rec["t"].shape)
df = pd.DataFrame({k: rec[k][:, 0] for k in ("t", "ax", "az", "q", "tas", "de", "alpha")})
df["h_baro"] = rec["h"][:, 0]
a = df["alpha"].to_numpy()
CL, CD, qd = coeffs_from_accel(df, AC, a)
# reconstruct truth CL from the model: CL0 + CLa a + CLq qhat + CLde de
qh = df["q"] * AC["c"] / (2 * df["tas"])
CLm = p["CL0"] + p["CLa"] * a + p["CLq"] * qh + T["CLde"] * df["de"]
print("CL from accel vs model: max abs diff", np.abs(CL - CLm).max(), "CL range", CL.min(), CL.max())
print("alpha deg range", np.degrees(a.min()), np.degrees(a.max()), "de deg range", np.degrees(df.de.min()), np.degrees(df.de.max()))
r = EE.estimate(df, AC, a)
for k in ("CL0", "CLa", "CLq", "CLde", "Cm0", "Cma", "Cmq", "Cmde"):
    print(k, round(r[k], 5), "truth", round(T[k], 5))
print("corr", np.round(r["regressor_corr"], 3))
w = windows(df); print("windows", w)
from flightloop.sysid.common import deriv
qdot = deriv(df["q"].to_numpy(), 0.02, 11, 3)
Cm_meas = AC["Iyy"] * qdot / (qd * AC["S"] * AC["c"])
Cm_mod = p["Cm0"] + p["Cma"] * a + p["Cmq"] * qh + T["Cmde"] * df["de"]
print("Cm meas range", Cm_meas.min(), Cm_meas.max(), "model range", Cm_mod.min(), Cm_mod.max())
print("lstsq ratio meas/model", np.polyfit(Cm_mod, Cm_meas, 1))
print("dt in log", np.diff(df.t)[:3], "Iyy", p["Iyy"], "Ixx", p["Ixx"], "Izz", p["Izz"])
from flightloop.sim import sixdof as S
z = np.zeros(len(df)); mach = df["tas"].to_numpy() / 340.0
out = S.SG.coeffs if hasattr(S.SG, "coeffs") else S.coeffs
CX, CY, CZ, Cl, Cm_sg, Cn, CLs, CDs, sig = S.coeffs(p, a, z, z, qh.to_numpy(), z, np.degrees(df["de"].to_numpy()), z, z, mach)
print("SG Cm vs my model: polyfit", np.polyfit(Cm_mod, Cm_sg, 1), "sig max", sig.max())
print("SG Cm vs measured: polyfit", np.polyfit(Cm_sg, Cm_meas, 1))
for win in (5, 7, 11):
    qd2 = deriv(df["q"].to_numpy(), 0.02, win, 3) if win > 5 else np.gradient(df["q"].to_numpy(), 0.02)
    Cmm = AC["Iyy"] * qd2 / (qd * AC["S"] * AC["c"])
    print("win", win, "slope", np.polyfit(Cm_sg, Cmm, 1)[0])
# spectrum of model Cm
f = np.fft.rfftfreq(len(Cm_sg), 0.02); P = np.abs(np.fft.rfft(Cm_sg - Cm_sg.mean())) ** 2
print("Cm power fraction above 2 Hz", P[f > 2].sum() / P.sum(), "above 5 Hz", P[f > 5].sum() / P.sum())
from flightloop.sysid.equation_error import _lp, d4
qf = _lp(df["q"].to_numpy(), 50, 5.0); qdf = d4(qf, 0.02)
Cmf_meas = AC["Iyy"] * qdf / (qd * AC["S"] * AC["c"])
Cmf_sg = _lp(Cm_sg, 50, 5.0)
print("filtered: meas vs SG slope", np.polyfit(Cmf_sg, Cmf_meas, 1))
# exact qdot from the model, lp filtered, vs d4 of lp q
qdot_true = qd * AC["S"] * AC["c"] * Cm_sg / AC["Iyy"]
print("qdot: d4(lp q) vs lp(qdot model) slope", np.polyfit(_lp(qdot_true, 50, 5.0), qdf, 1))
print("unfiltered qdot: gradient(q) vs model", np.polyfit(qdot_true, np.gradient(df['q'].to_numpy(), 0.02), 1))
from flightloop.sysid import output_error as OE
mid = rec["seg"][:, 0]; segs = rec["segments"]
exc = [(int(np.where(mid == i)[0][0]), int(np.where(mid == i)[0][-1]) + 1) for i, s in enumerate(segs) if s[1] in ("3211", "211", "doublet")]
t1 = time.time()
th0 = np.array([r["CLa"], r["Cma"], r["Cmq"], r["Cmde"]])
o = OE.estimate(df, AC, a, exc, r["CLde"], th0)
print("OE s", round(time.time() - t1, 1), {k: round(o[k], 4) for k in ("CLa", "Cma", "Cmq", "Cmde")}, "crb", {k: round(o["crb_" + k], 4) for k in ("CLa", "Cma", "Cmq", "Cmde")}, o["success"], o["resid_rms_alpha_rad"])
