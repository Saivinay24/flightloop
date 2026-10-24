"""Truth aircraft for synthetic flights: StratoGlide's design model (inputs/stratoglide/derivs.json).

derivs.json is exported once from StratoGlide's aero_tables.nominal() (their AVL derivatives + sizing polar) by
scripts/export_stratoglide.py. For identification the stall blend width is narrowed to 0.5 deg (StratoGlide uses
2 deg) so the model is linear inside the identification envelope (alpha <= 9 deg) and the 'true' derivatives are
well defined; this is logged in stubs/TRUTH.md.
"""
import json, math
from ..config import ROOT

DEG = math.pi / 180


def load(blend_w_deg=0.5):
    p = json.loads((ROOT / "inputs/stratoglide/derivs.json").read_text())["p"]
    p = dict(p)
    p["blend_w"] = blend_w_deg * DEG
    return p


def truth_params(p):
    """Identification targets in the units the estimators report (per rad for control derivatives)."""
    return {"CD0": p["CD0"], "k": p["k_ind"], "CL0": p["CL0"], "CLa": p["CLa"], "CLq": p["CLq"],
            "CLde": p["CLde"] / DEG, "Cm0": p["Cm0"], "Cma": p["Cma"], "Cmq": p["Cmq"], "Cmde": p["Cmde"] / DEG}
