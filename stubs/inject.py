"""Blind discrepancy injection for the Module F demo (MASTER_PLAN P5).

draw(seed) samples a discrepancy set, applies it to the design model and writes SEALED_truth.json. The loop engine
never opens this file (scripts/run_loop_blind.py logs every file the engine opens via a Python audit hook and
verify.py V18 greps that log). Ranges from the plan: CD0 +10 to +40%, CG aft 2 to 5% MAC, Cmde -15 to +5%,
prop efficiency -5 to -15 points (recorded; the blind flights are unpowered, so the prop rule is not scored).
"""
import json, pathlib
import numpy as np


def draw(seed):
    rng = np.random.default_rng(10_000 + seed)
    return {"cd0_pct": float(rng.uniform(10, 40)), "cg_aft_mac": float(rng.uniform(0.02, 0.05)),
            "cmde_pct": float(rng.uniform(-15, 5)), "prop_eta_pts": float(rng.uniform(-15, -5))}


def apply(p_design, d):
    """StratoGlide's coefficient model adds CL * dcg_mac to Cm, so an aft CG shift of d (fraction of MAC) is dcg_mac = +d."""
    p = dict(p_design)
    p["CD0"] = p_design["CD0"] * (1 + d["cd0_pct"] / 100)
    p["dcg_mac"] = p_design.get("dcg_mac", 0.0) + d["cg_aft_mac"]
    p["Cmde"] = p_design["Cmde"] * (1 + d["cmde_pct"] / 100)
    return p


def seal(path, seed, d, p_true):
    pathlib.Path(path).write_text(json.dumps({"seed": seed, "discrepancy": d, "truth_params": p_true}, indent=1, default=float))
