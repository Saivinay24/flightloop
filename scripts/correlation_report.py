"""Correlation report: StratoGlide's 25 pre-flight test-card predictions vs a FlightLoop flight-log identification.

No StratoGlide aircraft has flown. The "measured" column is a SYNTHETIC flight: blind seed 1 of the loop demo
(StratoGlide design model + an injected, sealed discrepancy, flown in FlightLoop's simulator with turbulence and sensor
noise). This shows the report FlightLoop produces when real logs arrive; it is not a flight result.
Columns: StratoGlide prediction and band (inputs/stratoglide/predictions.json), FlightLoop recomputation from the
design model (same formula as the measurement, to separate formula differences from aircraft differences), synthetic
measurement, status (inside / outside band / n/a with reason), and the next-iteration delta for outside-band items.
"""
import json, math, sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from flightloop.config import put_number
from flightloop.logs import schema
from flightloop.sysid.common import coeffs_from_accel

SCRIPT = "scripts/correlation_report.py"
G = 9.80665
pred = json.loads((ROOT / "inputs/stratoglide/predictions.json").read_text())
sg = json.loads((ROOT / "inputs/stratoglide/numbers.json").read_text())["numbers"]
num = lambda k: (sg[k]["value"] if isinstance(sg[k], dict) else sg[k])
blind = json.loads((ROOT / "results/loop/blind_summary.json").read_text())
rid = blind["run_id"]; seed = blind["seeds"][0]["seed"]
sd = ROOT / "runs" / rid / f"seed_{seed}"
eng = json.loads((sd / "engine" / "engine_result.json").read_text())
p = json.loads((ROOT / "inputs/stratoglide/derivs.json").read_text())["p"]
I = eng["identified"]
DEG = math.pi / 180
W = p["mass"] * G; S = p["S"]; c = p["c"]
V_sl = num("p4.v_tas_sl"); CL_loiter = num("p3.cl_loiter"); x_cg = num("p4.x_cg_m")


def perf(CD0, k):
    base = math.sqrt(2 * W / (1.225 * S))
    v_ms = base * (k / (3 * CD0)) ** 0.25; cl_ms = math.sqrt(3 * CD0 / k)
    return {"V_bestLD": base * (k / CD0) ** 0.25, "LD": 1 / (2 * math.sqrt(k * CD0)), "sink": v_ms * (4 * CD0) / cl_ms}


def stab(CLa, Cma, Cmq, Cm0, Cmde, CL0):
    qS = 0.5 * 1.225 * V_sl ** 2 * S
    Za = -qS * CLa / (p["mass"] * V_sl); Ma = qS * c * Cma / p["Iyy"]; Mq = qS * c * c * Cmq / (2 * V_sl * p["Iyy"])
    ev = np.linalg.eigvals(np.array([[Za, 1.0], [Ma, Mq]]))
    wn = float(abs(ev[0])); zeta = float(-ev[0].real / wn)
    a = (CL_loiter - CL0) / CLa
    de = -(Cm0 + Cma * a) / Cmde
    return {"wn": wn, "zeta": zeta, "trim_de_deg": math.degrees(de), "sm_pct": 100 * (-Cma / CLa),
            "dde_dcl": math.degrees(-Cma / (Cmde * CLa)), "x_np": x_cg + (-Cma / CLa) * c}


des = {"CD0": p["CD0"], "k": p["k_ind"], "CLa": p["CLa"], "Cma": p["Cma"], "Cmq": p["Cmq"], "Cm0": p["Cm0"], "Cmde": p["Cmde"] / DEG, "CL0": p["CL0"]}
pd_, pm = perf(des["CD0"], des["k"]), perf(I["CD0"], I["k"])
sd_, sm = stab(des["CLa"], des["Cma"], des["Cmq"], des["Cm0"], des["Cmde"], des["CL0"]), stab(I["CLa"], I["Cma"], I["Cmq"], I["Cm0"], I["Cmde"], I["CL0"])

# model-free trim gradient from the synthetic flight's steady-glide segments: mean de vs mean CL
df, meta = schema.read(sd / "flight1" / "flight1.parquet")
CL, CD, qd = coeffs_from_accel(df, {"S": S, "c": c, "mass": p["mass"]}, df["alpha"].to_numpy())
mid = df["maneuver_id"].to_numpy(); xs, ys = [], []
for i, sgm in enumerate(meta["segments"]):
    if sgm[1] == "glide":
        m = mid == i; xs.append(CL[m].mean()); ys.append(math.degrees(df["de"][m].mean()))
slope = float(np.polyfit(xs, ys, 1)[0])

MAP = {
    "V_bestLD EAS (SL)": ("V_bestLD", "perf"), "L/D max (SL)": ("LD", "perf"), "min sink rate (SL)": ("sink", "perf"),
    "trim elevator at loiter CL (SL)": ("trim_de_deg", "stab"), "static margin": ("sm_pct", "stab"),
    "neutral point x (from wing root LE)": ("x_np", "stab"), "short-period natural frequency (SL)": ("wn", "stab"),
    "short-period damping ratio (SL)": ("zeta", "stab"),
}
NA = {"stall EAS": "n/a (no stall test in the synthetic schedule; the identification envelope stops at alpha 9.6 deg)",
      "electrical power at SL cruise (loiter CL)": "n/a (blind flights are unpowered)", "climb rate at full power (SL)": "n/a (blind flights are unpowered)",
      "RPM at SL cruise": "n/a (blind flights are unpowered)", "Cl_p": "n/a (no roll excitation flown)",
      "altitude lost in pull-out": "n/a (pull-out cards FT-5/FT-6 not flown)", "peak load factor": "n/a (pull-out cards FT-5/FT-6 not flown)",
      "peak EAS": "n/a (pull-out cards FT-5/FT-6 not flown)", "time to recover": "n/a (pull-out cards FT-5/FT-6 not flown)"}
rows = []
for card in pred:
    q = card["quantity"]
    r = {"card": card["card"], "quantity": q, "unit": card["unit"], "predicted": card["value"], "band": [card["band_lo"], card["band_hi"]],
         "stratoglide_key": card["numbers_key"]}
    if q in MAP:
        k, grp = MAP[q]
        r["flightloop_design"] = (pd_ if grp == "perf" else sd_)[k]; r["measured"] = (pm if grp == "perf" else sm)[k]
    elif q == "trim gradient d(delta_e)/dCL":
        r["flightloop_design"] = sd_["dde_dcl"]; r["measured"] = slope; r["method"] = "regression of steady-glide mean elevator on mean CL (model-free)"
    elif q == "Cm_alpha":
        r["flightloop_design"] = des["Cma"]; r["measured"] = I["Cma"]
    elif q == "Cm_q":
        r["flightloop_design"] = des["Cmq"]; r["measured"] = I["Cmq"]
    elif q == "Cm_delta_e":
        r["flightloop_design"] = des["Cmde"] * DEG; r["measured"] = I["Cmde"] * DEG
    else:
        r["status"] = NA.get(q, "n/a (not measurable from the synthetic logs)"); rows.append(r); continue
    lo, hi = sorted(r["band"])
    r["status"] = "inside band" if lo <= r["measured"] <= hi else "outside band"
    rows.append(r)
deltas = {x["rule"]: x for x in eng["rules"]}
for r in rows:
    if r["status"] != "outside band":
        continue
    if r["quantity"] in ("V_bestLD EAS (SL)", "L/D max (SL)", "min sink rate (SL)"):
        r["next_iteration"] = deltas["R1 CD0"]["discrepancy"] + "; " + deltas["R1 CD0"]["formula"]
    elif r["quantity"] in ("static margin", "neutral point x (from wing root LE)", "Cm_alpha", "short-period natural frequency (SL)", "short-period damping ratio (SL)"):
        d = deltas["R2 static margin"]["delta"]
        r["next_iteration"] = f"move CG forward {1000 * d['cg_forward_m']:.1f} mm ({100 * d['cg_forward_mac']:.2f}% MAC) or raise V_h by {d.get('tail_volume_increase', float('nan')):.3f}"
    elif r["quantity"] in ("trim elevator at loiter CL (SL)", "trim gradient d(delta_e)/dCL", "Cm_delta_e"):
        d = deltas["R4 trim elevator offset"]["delta"]
        r["next_iteration"] = f"tail incidence change {d['tail_incidence_change_deg']:+.2f} deg (first-order, tau_e assumed)"
n_meas = sum(r["status"] in ("inside band", "outside band") for r in rows)
n_in = sum(r["status"] == "inside band" for r in rows)
out = {"label": "SYNTHETIC flight (blind loop seed); not a flight result", "loop_run_id": rid, "seed": seed, "n_cards": len(rows),
       "n_measured": n_meas, "n_inside": n_in, "rows": rows}
(ROOT / "results/loop/correlation_report.json").write_text(json.dumps(out, indent=1, default=float))
put_number("loop.correlation.n_cards", len(rows), "-", SCRIPT, "synthetic", rid)
put_number("loop.correlation.n_measured", n_meas, "-", SCRIPT, "synthetic", rid)
put_number("loop.correlation.n_inside_band", n_in, "-", SCRIPT, "synthetic", rid)
for r in rows:
    print(f"{r['card']:5s} {r['quantity'][:40]:40s} pred {r['predicted']:>9.4g}  design {r.get('flightloop_design', float('nan')):>9.4g}  meas {r.get('measured', float('nan')):>9.4g}  {r['status'][:30]}")


def _f(v):
    return "n/a" if v is None or (isinstance(v, float) and math.isnan(v)) else (f"{v:.4g}" if isinstance(v, (int, float)) else str(v))


L = ["# Correlation report (generated by scripts/correlation_report.py)", "",
     f"Measured column: SYNTHETIC blind-loop flight, run {rid}, seed {seed}. Not a flight result.", "",
     "## Table", "", "| Card | Quantity | Unit | StratoGlide prediction | Band | FlightLoop design recompute | Synthetic measured | Status |", "|---|---|---|---|---|---|---|---|"]
for r in rows:
    L.append(f"| {r['card']} | {r['quantity']} | {r.get('unit', '')} | {_f(r.get('predicted'))} | {_f(r['band'][0])} to {_f(r['band'][1])} | "
             f"{_f(r.get('flightloop_design'))} | {_f(r.get('measured'))} | {r.get('status', '')} |")
(ROOT / "results/loop/CORRELATION.md").write_text("\n".join(L) + "\n")
