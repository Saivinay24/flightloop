"""P4 post-processing: wedge shock angle + ramp pressure ratio, diamond Cl/Cd, 3-grid convergence, vs theory.py.

Wedge beta: on each vertical node column 0.3 <= x <= 1.3 (block-2 columns are exactly vertical), the outermost y where
p crosses the mid pressure (p_inf + p_ramp)/2; beta = atan of the least-squares slope through the corner (0, 0).
Ramp p2/p1: mean wall pressure for 0.2 <= x <= 1.2 divided by the freestream 101325 Pa.
Convergence state (HONEST_LOG 11, 13): "converged-residual" if rms[Rho] drops >= 6 orders; "converged-cauchy" if SU2
stopped before the iteration cap on its CD Cauchy criterion (1e-7 over 500) AND the CD max relative deviation over the
last 100 iterations is < 1e-6 (independent check from history.csv); otherwise "not-converged". CL is not used in the
relative check because it is 0 by symmetry at alpha 0.
"""
import glob, json, sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from flightloop.cfd.vtu import read_vtu
from flightloop.compressible import theory as TH
from flightloop.config import put_number

SCRIPT = "scripts/post_compressible.py"
man = json.loads((ROOT / "results/queues/compressible.json").read_text())
PINF = 101325.0
th = TH.table()


def hist(d):
    h = pd.read_csv(d / "history.csv"); h.columns = [c.strip().strip('"') for c in h.columns]
    r = h["rms[Rho]"].to_numpy()
    tail = h["CD"].iloc[-100:]
    dev = float((tail.max() - tail.min()) / abs(tail.iloc[-1]))
    drop = float(r[:20].max() - r[-1])
    cap = 6000  # ITER in every compressible case (su2_euler.build default; unchanged across both stopping rules)
    # what SU2 actually enforced is read from its own log (case.cfg of the first runs was rewritten later, see su2_euler.run)
    log = (d / "log.txt").read_text(errors="ignore")
    cauchy = any("Cauchy[CD]" in l and "Yes" in l for l in log.splitlines())
    if drop >= 6:
        state = "converged-residual"
    elif cauchy and len(h) < cap and dev < 1e-6:
        state = "converged-cauchy"
    else:
        state = "not-converged"
    return h, {"iters": int(len(h)), "iter_cap": cap, "rms_drop_orders": drop, "cd_rel_dev_last100": dev, "state": state}


def wedge_post(d):
    s = read_vtu(d / "surface_flow.vtu"); P = s["Points"]
    m = (P[:, 0] >= 0.2) & (P[:, 0] <= 1.2) & (P[:, 1] > 1e-9)
    p2 = float(s["Pressure"][m].mean())
    v = read_vtu(d / "flow.vtu"); X = v["Points"][:, :2]; p = v["Pressure"]
    pmid = 0.5 * (PINF + p2)
    xs, ys = [], []
    xr = np.round(X[:, 0], 9)
    for xc in np.unique(xr[(xr >= 0.3) & (xr <= 1.3)]):
        k = np.where(xr == xc)[0]; k = k[np.argsort(X[k, 1])]
        pp = p[k]; yy = X[k, 1]
        idx = np.where((pp[:-1] - pmid) * (pp[1:] - pmid) <= 0)[0]
        if len(idx) == 0:
            continue
        i = idx[-1]
        ys.append(yy[i] + (pmid - pp[i]) * (yy[i + 1] - yy[i]) / (pp[i + 1] - pp[i])); xs.append(xc)
    xs, ys = np.array(xs), np.array(ys)
    slope = float((xs @ ys) / (xs @ xs))
    return {"p2_p1": p2 / PINF, "beta_deg": float(np.degrees(np.arctan(slope))), "n_columns": int(len(xs))}


def gci3(f1, f2, f3, r=2.0):
    """f1 fine. Returns order, extrapolated, GCI_fine %, monotone flag."""
    e21, e32 = f2 - f1, f3 - f2
    mono = e21 * e32 > 0
    if not mono or e21 == 0:
        return {"monotone": bool(mono), "p": None, "f_ext": None, "gci_fine_pct": None}
    pp = np.log(abs(e32 / e21)) / np.log(r)
    fext = f1 + (f1 - f2) / (r ** pp - 1)
    return {"monotone": True, "p": float(pp), "f_ext": float(fext), "gci_fine_pct": float(1.25 * abs((f1 - f2) / f1) / (r ** pp - 1) * 100)}


out = {"wedge": {}, "diamond": {}, "runs": {}}
LV = {0: "coarse", 1: "medium", 2: "fine"}
for k, rid in man.items():
    case, lv, M, a = k.split("|"); lv = int(lv); M = float(M); a = float(a)
    d = ROOT / "runs" / rid
    if not (d / "exit_code").exists() or (d / "exit_code").read_text().strip() != "0":
        out["runs"][rid] = {"state": "missing"}; continue
    h, st = hist(d)
    st.update({"case": case, "level": lv, "mach": M, "aoa": a, "cells": None})
    out["runs"][rid] = st
    if case == "wedge":
        r = wedge_post(d); r.update(st); r["run_id"] = rid
        out["wedge"].setdefault(f"M{M:g}", {})[LV[lv]] = r
    else:
        r = {"cl": float(h["CL"].iloc[-1]), "cd": float(h["CD"].iloc[-1]), "run_id": rid, **st}
        out["diamond"].setdefault(f"M{M:g}_a{a:g}", {})[LV[lv]] = r

for key, e in out["wedge"].items():
    t = th["wedge"][key]
    e["theory"] = t
    if all(l in e for l in LV.values()):
        e["gci_p2p1"] = gci3(e["fine"]["p2_p1"], e["medium"]["p2_p1"], e["coarse"]["p2_p1"])
        e["fine_err"] = {"beta_deg": e["fine"]["beta_deg"] - t["beta_deg"], "p2_p1_pct": (e["fine"]["p2_p1"] / t["p2_p1"] - 1) * 100}
for key, e in out["diamond"].items():
    t = th["diamond"][key]; e["theory"] = t
    if all(l in e for l in LV.values()):
        e["gci_cd"] = gci3(e["fine"]["cd"], e["medium"]["cd"], e["coarse"]["cd"])
        e["monotone_cd"] = e["gci_cd"]["monotone"]
        if t["cl"] != 0:
            e["gci_cl"] = gci3(e["fine"]["cl"], e["medium"]["cl"], e["coarse"]["cl"])
            e["monotone_cl"] = e["gci_cl"]["monotone"]
            e["fine_err_cl_pct"] = (e["fine"]["cl"] / t["cl"] - 1) * 100
        else:
            e["monotone_cl"] = True  # alpha 0: Cl = 0 by symmetry, |Cl| reported
        e["fine_err_cd_pct"] = (e["fine"]["cd"] / t["cd"] - 1) * 100
        e["fine_err_cd_vs_ackeret_pct"] = (e["fine"]["cd"] / t["cd_ackeret"] - 1) * 100

dst = ROOT / "results/compressible"; dst.mkdir(parents=True, exist_ok=True)
(dst / "summary.json").write_text(json.dumps(out, indent=1))

# figures
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
for M, c in (("2", "C0"), ("3", "C1")):
    e = out["wedge"].get(f"M{M}", {})
    lv = [l for l in LV.values() if l in e]
    axs[0].plot(range(len(lv)), [e[l]["p2_p1"] for l in lv], "o-", c=c, label=f"CFD M{M}")
    axs[0].axhline(th["wedge"][f"M{M}"]["p2_p1"], c=c, ls="--", label=f"theory M{M}")
axs[0].set_xticks(range(3)); axs[0].set_xticklabels(list(LV.values())); axs[0].set_ylabel("p2/p1"); axs[0].legend(fontsize=7); axs[0].set_title("10 deg wedge")
for M, c in (("2", "C0"), ("3", "C1")):
    al = [0, 2, 4]
    ef = [out["diamond"].get(f"M{M}_a{a}", {}).get("fine", {}).get("cd", np.nan) for a in al]
    axs[1].plot(al, ef, "o", c=c, label=f"CFD fine M{M}")
    axs[1].plot(al, [th["diamond"][f"M{M}_a{a}"]["cd"] for a in al], "-", c=c, label=f"shock-expansion M{M}")
    axs[1].plot(al, [th["diamond"][f"M{M}_a{a}"]["cd_ackeret"] for a in al], ":", c=c, label=f"Ackeret M{M}")
axs[1].set_xlabel("alpha (deg)"); axs[1].set_ylabel("Cd"); axs[1].legend(fontsize=7); axs[1].set_title("5% diamond")
fig.tight_layout(); fig.savefig(ROOT / "results/figures/compressible.png", dpi=130)

src = "textbook verification (Euler CFD vs own theory.py)"
for M in ("2", "3"):
    e = out["wedge"].get(f"M{M}", {})
    if "fine_err" in e:
        put_number(f"comp.wedge.m{M}.beta_err_deg", e["fine_err"]["beta_deg"], "deg", SCRIPT, src, e["fine"]["run_id"])
        put_number(f"comp.wedge.m{M}.p2p1_err_pct", e["fine_err"]["p2_p1_pct"], "%", SCRIPT, src, e["fine"]["run_id"])
e = out["diamond"].get("M2_a2", {})
if "fine_err_cd_pct" in e:
    put_number("comp.diamond.m2a2.cl_err_pct", e["fine_err_cl_pct"], "%", SCRIPT, src, e["fine"]["run_id"])
    put_number("comp.diamond.m2a2.cd_err_pct", e["fine_err_cd_pct"], "%", SCRIPT, src, e["fine"]["run_id"])
print(json.dumps({k: {kk: (vv.get("fine_err") or {kx: vv.get(kx) for kx in ("fine_err_cl_pct", "fine_err_cd_pct", "monotone_cd")}) for kk, vv in v.items()} for k, v in out.items() if k != "runs"}, indent=1, default=str))

# markdown table for the report (scripts/build_docs.py includes it)
L = ["# Compressible verification (generated by scripts/post_compressible.py)", "", "## Table", "",
     "| Case | Quantity | Coarse | Medium | Fine | Theory | Fine error | States (c/m/f) |", "|---|---|---|---|---|---|---|---|"]
for k, w in out["wedge"].items():
    st_ = "/".join(w[l]["state"] for l in ("coarse", "medium", "fine"))
    L.append(f"| wedge {k}, 10 deg | shock angle (deg) | {w['coarse']['beta_deg']:.3f} | {w['medium']['beta_deg']:.3f} | {w['fine']['beta_deg']:.3f} | {w['theory']['beta_deg']:.3f} | {w['fine_err']['beta_deg']:+.3f} deg | {st_} |")
    L.append(f"| wedge {k}, 10 deg | p2/p1 | {w['coarse']['p2_p1']:.4f} | {w['medium']['p2_p1']:.4f} | {w['fine']['p2_p1']:.4f} | {w['theory']['p2_p1']:.4f} | {w['fine_err']['p2_p1_pct']:+.3f}% | |")
for k, d in out["diamond"].items():
    st_ = "/".join(d[l]["state"] for l in ("coarse", "medium", "fine"))
    if d["theory"]["cl"]:
        L.append(f"| diamond {k} | Cl | {d['coarse']['cl']:.5f} | {d['medium']['cl']:.5f} | {d['fine']['cl']:.5f} | {d['theory']['cl']:.5f} | {d['fine_err_cl_pct']:+.2f}% | |")
    L.append(f"| diamond {k} | Cd | {d['coarse']['cd']:.6f} | {d['medium']['cd']:.6f} | {d['fine']['cd']:.6f} | {d['theory']['cd']:.6f} | {d['fine_err_cd_pct']:+.2f}% | {st_} |")
(dst / "TABLE.md").write_text("\n".join(L) + "\n")
