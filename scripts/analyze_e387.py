"""P2 analysis: run registry, E387 vs TM-4062 (Cl, Cd, Cm, upper Cp, LSB vs oil flow, drag split), Celik GCI,
Trust Map (JSON + figure), LIMITS.md (generated, every number from this script), numbers.json keys.

Definitions follow PREREG.md fixed analysis choices 2 to 7. Unsteady runs (status rule 6) are reported with the
mean and half-amplitude of Cl/Cd over their last 2000 iterations, appear in the Trust Map, and are excluded from
error statistics (and therefore cannot pass a validation bar).
"""
import glob, json, sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from flightloop.config import cfg, put_number
from flightloop.cfd import postproc as PP
from flightloop.cfd.gci import celik

SCRIPT = "scripts/analyze_e387.py"
C = cfg()
TU = C["validation"]["e387"]["headline_tu_percent"]
OUT = ROOT / "results/cfd"; OUT.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "results/figures"
SRC = "public experimental (NASA TM-4062) vs CFD"
TM = ROOT / "data/tm4062/canonical"


def polar(re):
    p = pd.read_csv(TM / f"polar_re{int(re / 1000)}k.csv").groupby("alpha_deg", as_index=False).mean().sort_values("alpha_deg")
    return p


def exp_at(re, a):
    p = polar(re)
    return {k: float(np.interp(a, p["alpha_deg"], p[k])) for k in ("cl", "cd", "cm")}


def alpha_clmax(re):
    p = polar(re)
    return float(p["alpha_deg"].iloc[int(np.argmax(p["cl"].to_numpy()))])


# ------------------------------------------------------------------ registry
reg = []
for d in sorted(glob.glob(str(ROOT / "runs/cfd-*"))):
    d = pathlib.Path(d)
    if not (d / "case.json").exists():
        continue
    rj = d / "result.json"
    res = json.loads(rj.read_text()) if rj.exists() else None
    if res is None or res["status"]["state"] == "incomplete":
        try:
            res = PP.process(d)
        except Exception as e:
            res = {"case": json.loads((d / "case.json").read_text()), "status": {"state": "failed", "reason": str(e)[:120]}}
    c, st = res["case"], res["status"]
    unsteady = st["state"] == "unsteady"
    e = {"run_id": d.name, "airfoil": c["airfoil"], "re": c["re"], "alpha": c["alpha"], "model": c["model"],
         "level": c["level"], "tu_pct": c["tu_pct"], "tag_tu": c.get("tag_tu", False), "state": st["state"],
         "criterion": st.get("criterion"), "iters": st.get("iters"),
         "cl": st.get("cl_mean_last2000") if unsteady else st.get("cl"), "cd": st.get("cd_mean_last2000") if unsteady else st.get("cd"),
         "cl_amp": st.get("cl_amp_last2000"), "cd_amp": st.get("cd_amp_last2000"),
         "cm": res.get("cp_integral", {}).get("cm_c4_noseup"), "yplus_max": res.get("yplus_max"),
         "mesh_sha256_16": c.get("mesh_sha256_16"), "lsb": res.get("lsb_upper"), "breakdown": res.get("breakdown"),
         "reported": st["state"] in ("converged", "unsteady")}
    e["in_error_stats"] = bool(e["state"] == "converged" and e["level"] == "medium" and not e["tag_tu"] and abs(e["tu_pct"] - TU) < 1e-9)
    reg.append(e)
(OUT / "registry.json").write_text(json.dumps({"runs": reg}, indent=1, default=float))


def find(af, re, a, model, level="medium", tu=TU):
    for e in reg:
        if (e["airfoil"], e["re"], e["alpha"], e["model"], e["level"]) == (af, re, a, model, level) and abs(e["tu_pct"] - tu) < 1e-9 and not e["tag_tu"]:
            return e
    return None


def cp_rms(e, col):
    z = np.load(ROOT / "runs" / e["run_id"] / "surface.npz")
    x, cp = z["upper_x"], z["upper_cp"]; o = np.argsort(x)
    ex = pd.read_csv(TM / "cp_upper_re200k.csv")
    cfd = np.interp(ex["x_c"], x[o], cp[o])
    # the x/c = 0 orifice sits on the leading edge itself: compare with the CFD wall node of minimum x (either side),
    # not with an upper-surface interpolation clamped to the first upper node (that read Cp -0.11 vs the LE value)
    xa, cpa = np.r_[z["upper_x"], z["lower_x"]], np.r_[z["upper_cp"], z["lower_cp"]]
    cfd = np.where(ex["x_c"].to_numpy() == 0.0, cpa[np.argmin(xa)], cfd)
    return float(np.sqrt(np.mean((cfd - ex[col]) ** 2))), ex["x_c"].tolist(), ex[col].tolist(), cfd.tolist()


# ------------------------------------------------------------------ validation
oil = pd.read_csv(TM / "oilflow.csv")
val = {"definitions": "PREREG fixed choices 2-5", "tu_pct": TU}
CPCOL = {0: "cp_a0.01", 4: "cp_a3.99", 8: "cp_a8.00"}
for re in (200000, 100000):
    for a in (0, 4, 8):
        ex = exp_at(re, a)
        for m in ("SA", "SST_LM", "SST", "SA_BC"):
            e = find("e387", re, a, m)
            if e is None or not e["reported"]:
                continue
            v = {"run_id": e["run_id"], "state": e["state"], "in_error_stats": e["in_error_stats"], "exp": ex,
                 "cl": e["cl"], "cd": e["cd"], "cm": e["cm"], "cl_amp": e["cl_amp"], "cd_amp": e["cd_amp"],
                 "cl_err_pct": 100 * (e["cl"] - ex["cl"]) / ex["cl"] if abs(ex["cl"]) > 0.05 else None,
                 "cd_err_pct": 100 * (e["cd"] - ex["cd"]) / ex["cd"], "cm_err_abs": (e["cm"] - ex["cm"]) if e["cm"] is not None else None,
                 "drag_split": e["breakdown"]["CD"] if e["breakdown"] else None}
            if re == 200000:
                v["cp_rms_upper"], v["cp_x"], v["cp_exp"], v["cp_cfd"] = cp_rms(e, CPCOL[a])
            row = oil[(oil.re == re) & (np.isclose(oil.alpha_deg, a, atol=0.15))]
            lsb = e["lsb"] or {}
            if len(row):
                r0 = row.iloc[0]
                xs_e = float(r0.x_ls) if pd.notna(r0.x_ls) else None
                xr_e = float(r0.x_tr) if pd.notna(r0.x_tr) else None
                v["lsb"] = {"exp_x_s": xs_e, "exp_x_r": xr_e, "cfd_x_s": lsb.get("x_s"), "cfd_x_r": lsb.get("x_r"),
                            "xs_err": (lsb["x_s"] - xs_e) if (lsb.get("x_s") is not None and xs_e is not None) else None,
                            "xr_err": (lsb["x_r"] - xr_e) if (lsb.get("x_r") is not None and xr_e is not None) else None,
                            "exp_note": str(r0.note) if pd.notna(r0.note) else ""}
            val.setdefault(f"re{re // 1000}k_a{a}", {})[m] = v
val["re200k_a4"] = val.get("re200k_a4", {})
val["lsb_re200k"] = {str(a): (val.get(f"re200k_a{a}", {}).get("SST_LM", {}).get("lsb") or {"xs_err": None, "xr_err": None}) for a in (0, 4, 8)}
for a in ("0", "4"):
    e = val["lsb_re200k"][a]
    e["in_error_stats"] = val.get(f"re200k_a{a}", {}).get("SST_LM", {}).get("in_error_stats", False)
(OUT / "e387_validation.json").write_text(json.dumps(val, indent=1, default=float))

# ------------------------------------------------------------------ GCI
mesh = json.loads((ROOT / "results/mesh_e387.json").read_text())
N = {lv: mesh[lv]["cells"] for lv in ("coarse", "medium", "fine")}
gci = {"cells": N}
for m in ("SA", "SST_LM"):
    f = {lv: find("e387", 200000, 4, m, lv) for lv in ("fine", "medium", "coarse")}
    if all(x is not None and x["reported"] for x in f.values()):
        states = {lv: f[lv]["state"] for lv in f}
        gci[m] = {"states": states, "run_ids": {lv: f[lv]["run_id"] for lv in f}}
        for q in ("cl", "cd"):
            gci[m][q] = celik(N["fine"], N["medium"], N["coarse"], f["fine"][q], f["medium"][q], f["coarse"][q])
            if any(s != "converged" for s in states.values()):
                gci[m][q]["note"] = "at least one level unsteady: values are last-2000-iteration means"
                gci[m][q]["type"] = gci[m][q]["type"] + " (unsteady inputs)" if gci[m][q]["type"] == "monotone" else gci[m][q]["type"]
(OUT / "gci_e387.json").write_text(json.dumps(gci, indent=1, default=float))

# ------------------------------------------------------------------ Trust Map
cells = []
for e in reg:
    if e["airfoil"] != "e387":
        # other airfoils (SD7037) have no TM-4062 reference: listed with their convergence state; scored vs UIUC in
        # results/strato/DISAGREEMENT.md
        cells.append({"run_id": e["run_id"], "airfoil": e["airfoil"], "re": e["re"], "alpha": e["alpha"], "model": e["model"],
                      "level": e["level"], "tu_pct": e["tu_pct"],
                      "status": e["state"] if e["state"] != "converged" else "converged (scored in DISAGREEMENT.md)"})
        continue
    st = e["state"]
    label = st
    errs = {}
    if st in ("converged", "unsteady"):
        try:
            ex = exp_at(e["re"], e["alpha"])
            errs = {"cl_err_pct": 100 * (e["cl"] - ex["cl"]) / ex["cl"] if abs(ex["cl"]) > 0.05 else None,
                    "cd_err_pct": 100 * (e["cd"] - ex["cd"]) / ex["cd"]}
        except FileNotFoundError:
            pass
        if st == "unsteady":
            label = "unsteady"
        elif e["alpha"] > alpha_clmax(e["re"]) - 2:
            label = "near-stall"
        elif errs and (errs["cl_err_pct"] is None or abs(errs["cl_err_pct"]) < 5) and abs(errs["cd_err_pct"]) < 20:
            label = "validated"
        else:
            label = "out-of-bar"
    cells.append({"run_id": e["run_id"], "airfoil": e["airfoil"], "re": e["re"], "alpha": e["alpha"], "model": e["model"], "level": e["level"],
                  "tu_pct": e["tu_pct"], "status": label, **errs})
tm = {"labels": ["validated", "out-of-bar", "unsteady", "near-stall", "incomplete", "failed"],
      "bars": "validated = converged and |Cl err| < 5% and |Cd err| < 20% (V06 bars) vs TM-4062", "cells": cells}
(OUT / "trust_map_e387.json").write_text(json.dumps(tm, indent=1, default=float))
COL = {"validated": "#3fb37f", "out-of-bar": "#e0a030", "unsteady": "#c05050", "near-stall": "#7a6ad8", "incomplete": "#888", "failed": "#444"}
med = [c for c in cells if c["airfoil"] == "e387" and c["level"] == "medium" and abs(c["tu_pct"] - TU) < 1e-9]
models = sorted({c["model"] for c in med}); res_ = sorted({c["re"] for c in med})
fig, ax = plt.subplots(figsize=(7, 1.2 + 0.5 * len(models) * len(res_)))
yl = []
for i, (re, m) in enumerate([(r, m) for r in res_ for m in models]):
    yl.append(f"{m} Re {re // 1000}k")
    for c in med:
        if c["re"] == re and c["model"] == m:
            ax.scatter(c["alpha"], i, s=260, marker="s", c=COL.get(c["status"], "#888"))
            if c.get("cd_err_pct") is not None:
                ax.text(c["alpha"], i, f"{c['cd_err_pct']:+.0f}", ha="center", va="center", fontsize=6, color="white")
ax.set_yticks(range(len(yl))); ax.set_yticklabels(yl, fontsize=7); ax.set_xlabel("alpha (deg)")
for k, v in COL.items():
    ax.scatter([], [], c=v, marker="s", label=k)
ax.legend(fontsize=6, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.35))
ax.set_title("E387 Trust Map (number = Cd error % vs NASA TM-4062; public data)", fontsize=8)
fig.savefig(FIG / "trust_map_e387.png", dpi=140, bbox_inches="tight")

# ------------------------------------------------------------------ figures: Cp, polar
for a in (0, 4, 8):
    k = f"re200k_a{a}"
    if k not in val:
        continue
    fig, ax = plt.subplots(figsize=(6, 4))
    first = True
    for m, v in val[k].items():
        if not isinstance(v, dict) or "cp_x" not in v:
            continue
        z = np.load(ROOT / "runs" / v["run_id"] / "surface.npz")
        o = np.argsort(z["upper_x"])
        ax.plot(z["upper_x"][o], z["upper_cp"][o], label=f"{m} ({v['state']})")
        if first:
            ax.plot(v["cp_x"], v["cp_exp"], "ko", ms=3, label="TM-4062 App. D (public, uncorrected)")
            first = False
    ax.invert_yaxis(); ax.set_xlabel("x/c"); ax.set_ylabel("Cp (upper)"); ax.legend(fontsize=7)
    ax.set_title(f"E387 Re 200k alpha {a} deg, Tu {TU}%")
    fig.savefig(FIG / f"cp_e387_re200k_a{a}.png", dpi=130, bbox_inches="tight")

# ------------------------------------------------------------------ numbers
lm4 = val["re200k_a4"].get("SST_LM"); sa4 = val["re200k_a4"].get("SA")
if lm4:
    rid = lm4["run_id"]
    put_number("e387.lm.re200k.a4.cl_err_pct", lm4["cl_err_pct"], "%", SCRIPT, SRC, rid, f"state {lm4['state']}")
    put_number("e387.lm.re200k.a4.cd_err_pct", lm4["cd_err_pct"], "%", SCRIPT, SRC, rid, f"state {lm4['state']}")
    put_number("e387.lm.re200k.a4.cp_rms_upper", lm4["cp_rms_upper"], "-", SCRIPT, SRC, rid, f"state {lm4['state']}")
    put_number("e387.lm.re200k.a4.state", lm4["state"], "-", SCRIPT, SRC, rid)
if sa4:
    put_number("e387.sa.re200k.a4.cd_err_pct", sa4["cd_err_pct"], "%", SCRIPT, SRC, sa4["run_id"], f"state {sa4['state']}")
    put_number("e387.sa.re200k.a4.cl_err_pct", sa4["cl_err_pct"], "%", SCRIPT, SRC, sa4["run_id"], f"state {sa4['state']}")
xr = [val["lsb_re200k"][a]["xr_err"] for a in ("0", "4")]
xs = [val["lsb_re200k"][a]["xs_err"] for a in ("0", "4")]
if all(x is not None for x in xr):
    put_number("e387.lm.re200k.xr_err_c", max(abs(x) for x in xr), "c", SCRIPT, SRC, "multiple", "max |x_r err| over alpha 0, 4")
else:
    put_number("e387.lm.re200k.xr_err_c", "n/a (no reattachment found in at least one SST-LM run)", "c", SCRIPT, SRC, "multiple")
l4 = val["lsb_re200k"]["4"]
for nm in ("xs", "xr"):
    v4 = l4[f"{nm}_err"]
    put_number(f"e387.lm.re200k.a4.{nm}_err_c", v4 if v4 is not None else "n/a (no bubble found in CFD)", "c", SCRIPT, SRC,
               lm4["run_id"] if lm4 else "n/a", f"SST-LM alpha 4 vs TM-4062 oil flow; state {lm4['state'] if lm4 else 'missing'}")
if all(x is not None for x in xs):
    put_number("e387.lm.re200k.xs_err_c", max(abs(x) for x in xs), "c", SCRIPT, SRC, "multiple", "max |x_s err| over alpha 0, 4")
g = [gci[m]["cd"]["gci_fine_pct"] for m in ("SA", "SST_LM") if m in gci and gci[m]["cd"].get("gci_fine_pct") is not None]
if g:
    put_number("e387.gci.cd_fine_pct", max(g), "%", SCRIPT, "CFD numerical uncertainty", "multiple", "largest of SA, SST-LM with a GCI")
for m, k in (("SA", "sa"), ("SST_LM", "lm")):
    for a in (0, 4, 8):
        e = val[f"re200k_a{a}"].get(m)
        if e:
            put_number(f"e387.{k}.re200k.a{a}.cd_err_pct", e["cd_err_pct"], "%", SCRIPT, SRC, e["run_id"], f"state {e['state']}")
            put_number(f"e387.{k}.re200k.a{a}.cl_err_pct", e["cl_err_pct"], "%", SCRIPT, SRC, e["run_id"], f"state {e['state']}")
for m in ("SA", "SST_LM"):
    if m in gci:
        for q in ("cl", "cd"):
            v = gci[m][q]
            put_number(f"e387.gci.{m.lower()}.{q}_fine_pct", v["gci_fine_pct"] if v.get("gci_fine_pct") is not None else f"n/a ({v['type']})",
                       "%", SCRIPT, "CFD numerical uncertainty", "multiple")
rep = [e for e in reg if e["reported"] and e["airfoil"] == "e387"]
if rep:
    put_number("e387.yplus_max", max(e["yplus_max"] for e in rep), "-", SCRIPT, "CFD", "multiple", f"over {len(rep)} reported E387 runs")
validated_re = sorted({e["re"] for e in reg if e["in_error_stats"] and e["airfoil"] == "e387"})
put_number("e387.re_validated", " and ".join(f"{r // 1000}k" for r in validated_re) if validated_re else "n/a (no converged run)", "-", SCRIPT, SRC, "multiple")
put_number("e387.n_runs_reported", len(rep), "-", SCRIPT, "CFD", "multiple")
put_number("e387.n_runs_unsteady", sum(e["state"] == "unsteady" for e in rep), "-", SCRIPT, "CFD", "multiple")

# ------------------------------------------------------------------ LIMITS.md (generated)
L = ["# LIMITS: where the E387 CFD fails (generated by scripts/analyze_e387.py; every number is from results/cfd/*.json)", "",
     f"Reference: NASA TM-4062 (public experimental), Re 200k unless stated, headline Tu {TU}%. Error = 100 (CFD - exp) / exp.", "",
     "## Quantified failure regions", ""]
n = 0
for a in (0, 4, 8):
    v = val.get(f"re200k_a{a}", {}).get("SA")
    if v:
        n += 1
        L.append(f"{n}. SA (fully turbulent) at alpha {a} deg: Cd error {v['cd_err_pct']:+.1f}% (friction share {100 * v['drag_split']['friction'] / v['drag_split']['total']:.0f}% of Cd), because SA has no laminar run and no bubble.")
for a in (0, 4, 8):
    v = val.get(f"re200k_a{a}", {}).get("SST_LM")
    if v:
        n += 1
        amp = f", Cd oscillation half-amplitude {100 * v['cd_amp'] / v['cd']:.1f}% of mean over the last 2000 iterations" if v["state"] == "unsteady" else ""
        lsb = v.get("lsb") or {}
        lsbt = ""
        if lsb:
            lsbt = (f"; bubble x_s {lsb['cfd_x_s']:.3f}c vs oil flow {lsb['exp_x_s']}" if lsb.get("cfd_x_s") is not None and lsb.get("exp_x_s") is not None else "; no separation found in CFD") + \
                   (f", x_r {lsb['cfd_x_r']:.3f}c vs {lsb['exp_x_r']}" if lsb.get("cfd_x_r") is not None and lsb.get("exp_x_r") is not None else "")
        L.append(f"{n}. SST-LM at alpha {a} deg: state {v['state']}{amp}; Cl error {v['cl_err_pct']:+.1f}%, Cd error {v['cd_err_pct']:+.1f}%{lsbt}.")
for m in ("SA", "SST_LM"):
    if m in gci:
        n += 1
        L.append(f"{n}. Grid dependence {m} alpha 4 deg: Cl {gci[m]['cl']['type']}" + (f" GCI {gci[m]['cl']['gci_fine_pct']:.2f}%" if gci[m]['cl'].get('gci_fine_pct') is not None else "") +
                 f", Cd {gci[m]['cd']['type']}" + (f" GCI {gci[m]['cd']['gci_fine_pct']:.2f}% (p {gci[m]['cd']['p']:.2f})" if gci[m]['cd'].get('gci_fine_pct') is not None else "") + ".")
L += ["", "## Not covered", "", "- No URANS, no 3D, no wind-tunnel wall or Pt corrections beyond what TM-4062 Table B1 already applies (Appendix D Cp is uncorrected).",
      "- Lower-surface Cp was not transcribed (GAPS.md).", ""]
(ROOT / "LIMITS.md").write_text("\n".join(L))
print(json.dumps({"registry": len(reg), "states": {s: sum(e["state"] == s for e in reg) for s in {e["state"] for e in reg}},
                  "re200k_a4": {m: {k: v.get(k) for k in ("state", "cl_err_pct", "cd_err_pct", "cp_rms_upper")} for m, v in val["re200k_a4"].items()},
                  "lsb": val["lsb_re200k"], "gci": {m: {q: {k: gci[m][q].get(k) for k in ("type", "p", "gci_fine_pct")} for q in ("cl", "cd")} for m in ("SA", "SST_LM") if m in gci}}, indent=1, default=str))

# markdown tables for the report (scripts/build_docs.py includes them)
def _f(v, s):
    return "n/a" if v is None else format(v, s)


T = ["# E387 Re 200k validation (generated by scripts/analyze_e387.py)", "", "## Validation", "",
     "| alpha | Model | State | Cl | Cl exp | Cl err | Cd | Cd exp | Cd err | Cd half-amplitude | RMS dCp upper |", "|---|---|---|---|---|---|---|---|---|---|---|"]
for a in (0, 4, 8):
    for m, e in val[f"re200k_a{a}"].items():
        T.append(f"| {a} | {m} | {e['state']} | {e['cl']:.4f} | {e['exp']['cl']:.4f} | {e['cl_err_pct']:+.1f}% | {e['cd']:.5f} | {e['exp']['cd']:.5f} | "
                 f"{e['cd_err_pct']:+.1f}% | {_f(e['cd_amp'], '.5f') if e['cd_amp'] else 'steady'} | {e['cp_rms_upper']:.3f} |")
T += ["", "## Bubble", "", "| alpha | Oil flow x_s | Oil flow x_r | SST-LM x_s | SST-LM x_r | Error x_s | Error x_r | Note |", "|---|---|---|---|---|---|---|---|"]
for a, e in val["lsb_re200k"].items():
    T.append(f"| {a} | {_f(e['exp_x_s'], '.2f')} | {_f(e['exp_x_r'], '.2f')} | {_f(e['cfd_x_s'], '.3f')} | {_f(e['cfd_x_r'], '.3f')} | "
             f"{_f(e['xs_err'], '+.3f')} | {_f(e['xr_err'], '+.3f')} | {e.get('exp_note') or ('run not in error statistics' if e.get('in_error_stats') is False else '')} |")
T += ["", "## GCI", "", "| Model | Quantity | States (fine/medium/coarse) | Fine / medium / coarse | Type | p | GCI fine |", "|---|---|---|---|---|---|---|"]
for m in ("SA", "SST_LM"):
    if m not in gci:
        continue
    for q in ("cl", "cd"):
        e = gci[m][q]
        T.append(f"| {m} | {q.upper()} | {'/'.join(gci[m]['states'][l] for l in ('fine', 'medium', 'coarse'))} | {' / '.join(f'{x:.5f}' for x in e['f'])} | "
                 f"{e['type']} | {_f(e['p'], '.2f')} | {_f(e['gci_fine_pct'], '.2f') + '%' if e['gci_fine_pct'] is not None else 'n/a (' + e['type'] + ')'} |")
(OUT / "VALIDATION.md").write_text("\n".join(T) + "\n")
