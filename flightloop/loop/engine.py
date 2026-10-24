"""Next-iteration engine (Module F). Input: a directory of canonical flight logs + the design model. Output: discrepancy
table, first-order design deltas, an updated model, and next-flight predictions.

The engine reads only: the log directory it is given (Parquet + JSON meta) and the design model file. It never reads
truth files; scripts/run_loop_blind.py records every file it opens.
"""
import json, math, pathlib
import numpy as np
import pandas as pd

from ..logs import schema
from ..sysid import glide, equation_error as EE
from ..atmos import isa
from . import deltas as D

DEG = math.pi / 180
# Identification envelope: the design model is linear only for alpha <= 9 deg (flightloop/sim/truth.py, PREREG 12).
# Samples (equation error) and glide windows (force method) outside it are excluded from the regressions. Added after
# the first blind run (HONEST_LOG 14): an injected aft CG pushed 10 to 12% of flight-1 samples to 16 to 17 deg.
ALPHA_ENV_DEG = 9.0


def design_model(path):
    p = json.loads(pathlib.Path(path).read_text())["p"]
    p = dict(p); p["blend_w"] = 0.5 * DEG
    return p


def as_ident(p):
    return {"CD0": p["CD0"], "k": p["k_ind"], "CL0": p["CL0"], "CLa": p["CLa"], "Cm0": p["Cm0"] + p["CL0"] * p.get("dcg_mac", 0.0),
            "Cma": p["Cma"] + p["CLa"] * p.get("dcg_mac", 0.0), "Cmq": p["Cmq"], "Cmde": p["Cmde"] / DEG}


def identify(log_dir, ac):
    rows = []
    for f in sorted(pathlib.Path(log_dir).glob("*.parquet")):
        df, meta = schema.read(f)
        alpha = df["alpha"].to_numpy()
        segs = meta["segments"]; mid = df["maneuver_id"].to_numpy()
        wins, dropped = [], 0
        for i, sg in enumerate(segs):
            if sg[1] == "glide":
                ii = np.where(mid == i)[0]
                e = np.linspace(ii[0], ii[-1] + 1, 4).astype(int)
                for j in range(3):
                    a0, a1 = int(e[j]), int(e[j + 1])
                    if np.degrees(alpha[a0:a1]).max() <= ALPHA_ENV_DEG:
                        wins.append((a0, a1))
                    else:
                        dropped += 1
        g = glide.estimate(df, ac, alpha, "force", wins)
        inside = np.degrees(alpha) <= ALPHA_ENV_DEG
        ee = EE.estimate(df, ac, alpha, mask=inside)
        rows.append({"CD0": g["CD0"], "k": g["k"], "se_CD0": g["se_CD0"], "CL0": ee["CL0"], "CLa": ee["CLa"], "Cm0": ee["Cm0"],
                     "Cma": ee["Cma"], "Cmq": ee["Cmq"], "Cmde": ee["Cmde"], "CLde": ee["CLde"], "se_Cma": ee["se_Cma"], "se_CLa": ee["se_CLa"],
                     "file": f.name, "glide_windows_used": len(wins), "glide_windows_outside_envelope": dropped,
                     "ee_samples_outside_envelope_frac": float(1 - inside.mean())})
    return {k: float(np.mean([r[k] for r in rows])) for k in rows[0] if k != "file"} | {"n_flights": len(rows), "alpha_envelope_deg": ALPHA_ENV_DEG}


def updated_model(p_design, ident):
    p = dict(p_design)
    p.update({"CD0": ident["CD0"], "k_ind": ident["k"], "CL0": ident["CL0"], "CLa": ident["CLa"], "Cm0": ident["Cm0"],
              "Cma": ident["Cma"], "Cmq": ident["Cmq"], "Cmde": ident["Cmde"] * DEG, "dcg_mac": 0.0})
    return p


def predict(p, V, rho):
    """Steady glide at TAS V: (sink m/s, trim elevator deg) from the model's trim solution."""
    from ..sim.sixdof import trim_glide
    a, de, gam, th = trim_glide(p, V, rho)
    return float(-V * math.sin(gam)), float(de)


def run(log_dir, design_path, out_dir):
    out_dir = pathlib.Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    p_d = design_model(design_path)
    ac = {k: p_d[k] for k in ("S", "c", "b", "mass", "Iyy")}
    ident = identify(log_dir, ac)
    des = as_ident(p_d)
    disc = {k: {"design": des[k], "flight": ident[k], "diff_pct": 100 * (ident[k] / des[k] - 1)} for k in des}
    tail = {"eta_t": 0.9, "a_t": 0.8 * p_d["CLa"], "deda": 2 * p_d["CLa"] / (math.pi * p_d["b"] ** 2 / p_d["S"])}
    rl = D.rules(des, ident, ac, rho=1.1, tail=tail)
    p_u = updated_model(p_d, ident)
    res = {"identified": ident, "discrepancy": disc, "rules": rl, "tail_assumptions": tail,
           "cd0_delta_pct": 100 * (ident["CD0"] / des["CD0"] - 1),
           "sm_delta_mac": (-ident["Cma"] / ident["CLa"]) - (-des["Cma"] / des["CLa"]), "first_order": True}
    (out_dir / "engine_result.json").write_text(json.dumps(res, indent=1, default=float))
    (out_dir / "model_updated.json").write_text(json.dumps({"p": p_u}, indent=1, default=float))
    return res, p_d, p_u
