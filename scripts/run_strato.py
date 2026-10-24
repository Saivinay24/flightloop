"""P7 Module B: stratospheric-Re disagreement table for SD7037 (StratoGlide's airfoil) at Re 60k, 80k, 100k.

Tools: SU2 SST-LM and SA (Tu 0.06%, medium mesh, alpha 2 and 4), XFOIL Ncrit 9 and 11, NeuralFoil xxlarge,
corrected XFOIL (Module C final model, refit deterministically from results/mf/points.csv), UIUC LSAT experiment
(SD7037 (B), Re 60.5k and 101.4k; public experimental). Metrics: mean |Cd error| over CL 0.5 to 0.8 vs experiment,
Cl_max, upper-surface transition / bubble location. Re 80k has no experiment: errors are n/a (no experiment).
"""
import json, sys, pathlib, importlib.util
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import pandas as pd
from flightloop.config import put_number
from flightloop.mf import uiuc as U, xfoil as X
from flightloop.geometry import read_selig, with_le_point

SCRIPT = "scripts/run_strato.py"
spec = importlib.util.spec_from_file_location("mf", ROOT / "scripts/run_mf.py"); MF = importlib.util.module_from_spec(spec); spec.loader.exec_module(MF)
OUT = ROOT / "results/strato"; OUT.mkdir(parents=True, exist_ok=True)
CLS = np.array([0.5, 0.6, 0.7, 0.8])
RES = [60000, 80000, 100000]
EXP_FOIL = "SD7037 (B)"

C = U.coords(); xy = C[U.norm(EXP_FOIL)]["xy"]
P = U.polars(); P = P[P.airfoil == EXP_FOIL]
exp = {}
for re_ in (60000, 100000):
    g = P[(P.re - re_).abs() < 3000].sort_values("alpha")
    exp[re_] = g
mf_sum = json.loads((ROOT / "results/mf/summary.json").read_text())
pts = pd.read_csv(ROOT / "results/mf/points.csv")
dev = pts[~pts.base.str.replace(" ", "").isin(MF.TEST)].reset_index(drop=True)
Xd, yd = MF.feats(dev)
make = MF.candidates()[mf_sum["test_airfoils"]["model"]]
model = make()
if "gp" in mf_sum["test_airfoils"]["model"] and len(Xd) > 900:
    idx = np.random.default_rng(0).choice(len(Xd), 900, replace=False); model.fit(Xd[idx], yd[idx])
else:
    model.fit(Xd, yd)
gm = MF.geom(xy)


def cd_at(rows, cols, cl):
    r = MF.at_cl(rows, cl, cols)
    return r


def exp_cd(re_, cl):
    g = exp.get(re_)
    if g is None:
        return None
    a = g[g.cl < 0.9 * g.cl.max()].sort_values("cl")
    if cl < a.cl.min() or cl > a.cl.max():
        return None
    return float(np.interp(cl, a.cl, a.cd))


# Low-order tools run at the experiment's own Re where one exists (60.5k, 101.4k): at Re 60000 exactly, XFOIL Ncrit 9
# converged to a non-physical branch (alpha 4: Cd 0.0102 with friction part 0.0021, below the one-sided laminar
# flat-plate value 0.0054), while at 60500 it gives Cd 0.0250 (HONEST_LOG 15). SU2 rows stay at the run Re.
RE_TOOL = {re_: (float(round(exp[re_].re.median(), -2)) if re_ in exp else float(re_)) for re_ in RES}


def physical(pol):
    """Drop XFOIL points whose friction drag (cd - cdp) is below the one-sided laminar flat-plate value 1.328/sqrt(Re)."""
    lim = 1.328 / np.sqrt(pol["re"])
    keep = [r for r in pol["rows"] if r[2] - r[3] >= lim]
    return dict(pol, rows=keep, n_nonphysical=len(pol["rows"]) - len(keep))


tools = {}
for re_ in RES:
    x9, x11 = physical(X.polar(xy, "sd7037b", RE_TOOL[re_], 9)), physical(X.polar(xy, "sd7037b", RE_TOOL[re_], 11))
    nfp = MF.nf_polar(xy, RE_TOOL[re_])
    T = {"xfoil_n9": {}, "xfoil_n11": {}, "neuralfoil": {}, "xfoil_corrected": {}}
    for cl in CLS:
        a9 = cd_at(x9["rows"], x9["cols"], cl); a11 = cd_at(x11["rows"], x11["cols"], cl); an = cd_at(nfp, ["alpha", "cl", "cd", "conf"], cl)
        T["xfoil_n9"][cl] = a9["cd"] if a9 else None
        T["xfoil_n11"][cl] = a11["cd"] if a11 else None
        T["neuralfoil"][cl] = an["cd"] if an else None
        if a9:
            row = pd.DataFrame([{"re": RE_TOOL[re_], "alpha_xf": a9["alpha"], "cl": cl, "cd_xf9": a9["cd"], "cm_xf": a9["cm"], "xtr_top": a9["xtr_top"],
                                 "xtr_bot": a9["xtr_bot"], "cd_exp": 1.0, **gm}])
            T["xfoil_corrected"][cl] = float(a9["cd"] * np.exp(model.predict(MF.feats(row)[0])[0]))
        else:
            T["xfoil_corrected"][cl] = None
    clmax = {"xfoil_n9": max(r[1] for r in x9["rows"]) if x9["rows"] else None, "xfoil_n11": max(r[1] for r in x11["rows"]) if x11["rows"] else None, "neuralfoil": max(r[1] for r in nfp)}
    xtr4 = {k: (MF.at_cl(pp["rows"], 0.7, pp["cols"]) or {}).get("xtr_top") for k, pp in (("xfoil_n9", x9), ("xfoil_n11", x11))}
    tools[re_] = {"re_tool": RE_TOOL[re_], "xfoil_nonphysical_points": {"n9": x9["n_nonphysical"], "n11": x11["n_nonphysical"]}, "cd_at_cl": T, "clmax": clmax, "xtr_top_at_cl0.7": xtr4, "nf_conf_median": float(np.median([r[3] for r in nfp]))}

# SU2 runs (registry built by analyze_e387.py includes sd7037 runs)
reg = json.loads((ROOT / "results/cfd/registry.json").read_text())["runs"]
su2 = {}
for e in reg:
    if e["airfoil"] != "sd7037" or e["level"] != "medium":
        continue
    su2.setdefault(e["re"], {}).setdefault(e["model"], {})[e["alpha"]] = e


def cfd_cd_at_cl(re_, m):
    d = su2.get(re_, {}).get(m, {})
    pts_ = sorted((e["cl"], e["cd"], e["state"], e["lsb"]) for e in d.values() if e["state"] in ("converged", "unsteady"))
    return pts_


rows = []
for re_ in RES:
    for tool in ("SU2 SST-LM", "SU2 SA", "XFOIL Ncrit 9", "XFOIL Ncrit 11", "NeuralFoil xxlarge", "XFOIL corrected (Module C)", "UIUC experiment"):
        r = {"re": re_, "tool": tool}
        if tool.startswith("SU2"):
            m = "SST_LM" if "LM" in tool else "SA"
            pp = cfd_cd_at_cl(re_, m)
            if not pp:
                r.update(cd_err="n/a (run not in compute budget)" if re_ == 80000 else "n/a (run not finished)", clmax="n/a (alpha sweep not run)", bubble="n/a")
            else:
                errs = []
                for cl_, cd_, st, lsb in pp:
                    e_ = exp_cd(re_, cl_)
                    if e_ is not None:
                        errs.append((100 * (cd_ - e_) / e_, st))
                r["cd_err"] = (f"{np.mean([abs(x) for x, _ in errs]):.1f}% (at CL " + ", ".join(f"{c:.2f}" for c, *_ in pp) + "; " + ", ".join(sorted({s for _, s in errs})) + ")") if errs else ("n/a (no experiment)" if re_ == 80000 else "n/a (CFD CL outside experimental range)")
                r["clmax"] = "n/a (alpha sweep not run)"
                l = [x[3] for x in pp if x[3] and x[3].get("bubble")]
                r["bubble"] = (f"x_s {l[0]['x_s']:.2f}c, x_r " + (f"{l[0]['x_r']:.2f}c" if l[0].get('x_r') is not None else "none") + f" at CL {pp[0][0]:.2f}") if l else "no bubble in CFD"
        elif tool == "UIUC experiment":
            g = exp.get(re_)
            r.update(cd_err="reference" if g is not None else "n/a (no experiment at this Re)", clmax=f"{g.cl.max():.3f}" if g is not None else "n/a (no experiment)",
                     bubble="n/a (no flow visualisation in LSAT vol. 3 for SD7037)")
        else:
            key = {"XFOIL Ncrit 9": "xfoil_n9", "XFOIL Ncrit 11": "xfoil_n11", "NeuralFoil xxlarge": "neuralfoil", "XFOIL corrected (Module C)": "xfoil_corrected"}[tool]
            vals = tools[re_]["cd_at_cl"][key]
            errs = [100 * (vals[cl] - exp_cd(re_, cl)) / exp_cd(re_, cl) for cl in CLS if vals[cl] is not None and exp_cd(re_, cl) is not None]
            r["cd_err"] = f"{np.mean(np.abs(errs)):.1f}%" if errs else "n/a (no experiment)" if re_ == 80000 else ("n/a (XFOIL did not converge)" if all(vals[c] is None for c in CLS) else "n/a (no overlap)")
            r["cd_err_signed_mean"] = float(np.mean(errs)) if errs else None
            cm_ = tools[re_]["clmax"].get(key, "corr")
            r["clmax"] = "n/a (correction is for Cd only)" if cm_ == "corr" else ("n/a (XFOIL halted: 4 non-converged points from alpha 0)" if cm_ is None else f"{cm_:.3f}")
            r["re_tool"] = RE_TOOL[re_]
            xt = tools[re_]["xtr_top_at_cl0.7"].get(key)
            r["bubble"] = f"x_tr upper {xt:.2f}c at CL 0.7" if xt is not None else "n/a (no transition output)"
        rows.append(r)

# Mack: Ncrit = -8.43 - 2.4 ln(Tu)
mack = {n: float(np.exp(-(n + 8.43) / 2.4) * 100) for n in (9, 11)}
expl = [
    f"XFOIL Ncrit 9 vs 11: Mack's relation Ncrit = -8.43 - 2.4 ln(Tu) maps Ncrit 9 to Tu {mack[9]:.3f}% and Ncrit 11 to Tu {mack[11]:.3f}%; the spread between the two rows is the sensitivity to the (unknown) disturbance level, which is the argument for quieter stratospheric air, not a measurement of it.",
    "XFOIL's e^N envelope with viscous-inviscid coupling is calibrated on natural Tollmien-Schlichting transition in quiet tunnels, which is the mechanism that sets the bubble at this Re; its Cd error here is the size of the residual that Module C learns.",
    "SU2 SST-LM uses Langtry-Menter correlations calibrated mostly on bypass transition; the result depends on the Tu that survives to the leading edge (freestream decay is logged per run).",
    "SU2 SA is fully turbulent: no laminar run and no bubble, so it over-predicts friction drag by construction at Re 60k to 100k.",
    "NeuralFoil is trained on XFOIL output, so it inherits XFOIL's bias; its authors flag Re near 80k as a hard regime (analysis_confidence median reported per Re).",
    "Corrected XFOIL applies the Module C model trained on every non-test family; SD7037 never entered training, so its row is a genuine out-of-sample test.",
    "UIUC LSAT vol. 3 has drag polars at Re 60.5k and 101.4k for SD7037 (B); there is no Re 80k experiment, so the 80k rows report tool spread only.",
]
L = ["# DISAGREEMENT: SD7037 at stratospheric Reynolds numbers (generated by scripts/run_strato.py)", "",
     "Airfoil: SD7037 (B), as-tested UIUC coordinates (StratoGlide's selected airfoil family). Data: UIUC LSAT (public experimental), SU2 (CFD), XFOIL/NeuralFoil (low-order). Cd error = mean |Cd_tool - Cd_exp| / Cd_exp over CL 0.5, 0.6, 0.7, 0.8 (SU2: at the CL of each run).", "",
     "Low-order tools run at the experiment's Re (60.5k, 101.4k); SU2 runs at 60k and 100k. XFOIL points with friction drag below the one-sided laminar flat-plate value 1.328/sqrt(Re) are dropped as non-physical: " + ", ".join(f"Re {re_ // 1000}k Ncrit 9: {tools[re_]['xfoil_nonphysical_points']['n9']}, Ncrit 11: {tools[re_]['xfoil_nonphysical_points']['n11']}" for re_ in RES) + ".", "",
     "| Re | Tool | Cd error vs experiment | Cl_max | Bubble / transition |", "|---|---|---|---|---|"]
for r in rows:
    rr = r.get("re_tool", r["re"])
    L.append(f"| Re {r['re'] // 1000}k" + (f" (run at {rr / 1000:.1f}k)" if rr != r["re"] else "") + f" | {r['tool']} | {r['cd_err']} | {r['clmax']} | {r['bubble']} |")
L += ["", "## Why the tools disagree (each item is tested against the table, not asserted)", ""] + [f"{i + 1}. {e}" for i, e in enumerate(expl)]
L += ["", f"NeuralFoil analysis_confidence (median over the alpha sweep): " + ", ".join(f"Re {re_ // 1000}k {tools[re_]['nf_conf_median']:.2f}" for re_ in RES) + "."]
(OUT / "DISAGREEMENT.md").write_text("\n".join(L) + "\n")
(OUT / "disagreement.json").write_text(json.dumps({"rows": rows, "tools": {str(k): v for k, v in tools.items()}, "mack_tu_pct": mack}, indent=1, default=lambda o: float(o) if hasattr(o, "__float__") else str(o)))
for re_ in (60000, 100000):
    for r in rows:
        if r["re"] == re_ and r["tool"] == "XFOIL Ncrit 9" and r.get("cd_err_signed_mean") is not None:
            put_number(f"strato.sd7037.re{re_ // 1000}k.xfoil_n9_cd_err_pct", r["cd_err_signed_mean"], "%", SCRIPT, "public experimental (UIUC) vs XFOIL", "n/a")
        if r["re"] == re_ and r["tool"] == "XFOIL corrected (Module C)" and r.get("cd_err_signed_mean") is not None:
            put_number(f"strato.sd7037.re{re_ // 1000}k.corrected_cd_err_pct", r["cd_err_signed_mean"], "%", SCRIPT, "public experimental (UIUC) vs corrected XFOIL", "n/a")
print("\n".join(L))
