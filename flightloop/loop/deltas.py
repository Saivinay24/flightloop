"""Rule table: flight-vs-model discrepancies -> first-order design deltas (MASTER_PLAN P5). Each rule cites its formula.

All deltas are first-order (linear, small-perturbation) and are labelled so in the outputs.
"""
import math

G = 9.80665


def glide_speeds(W, rho, S, CD0, k):
    """V_best_glide = sqrt(2W/(rho S)) (k/CD0)^(1/4); V_min_sink = sqrt(2W/(rho S)) (k/(3 CD0))^(1/4); (L/D)max = 1/(2 sqrt(k CD0))."""
    base = math.sqrt(2 * W / (rho * S))
    return {"V_best_glide": base * (k / CD0) ** 0.25, "V_min_sink": base * (k / (3 * CD0)) ** 0.25, "LD_max": 1 / (2 * math.sqrt(k * CD0))}


def rules(design, ident, ac, rho=1.225, tail=None, tau_e=0.5):
    """design/ident: dicts with CD0, k, CLa, Cma, Cm0, Cmde (per rad). ac: S, c, mass. Returns list of rule outputs."""
    W = ac["mass"] * G
    out = []
    # R1 CD0
    d_cd0 = 100 * (ident["CD0"] / design["CD0"] - 1)
    gs_d, gs_i = glide_speeds(W, rho, ac["S"], design["CD0"], design["k"]), glide_speeds(W, rho, ac["S"], ident["CD0"], ident["k"])
    out.append({"rule": "R1 CD0", "discrepancy": f"CD0 {d_cd0:+.1f}% vs design", "value_pct": d_cd0,
                "delta": {"CD0_new": ident["CD0"], "V_best_glide_ms": gs_i["V_best_glide"], "V_min_sink_ms": gs_i["V_min_sink"],
                          "LD_max": gs_i["LD_max"], "LD_max_design": gs_d["LD_max"]},
                "formula": "V_bg = sqrt(2W/(rho S)) (k/CD0)^(1/4); V_ms = sqrt(2W/(rho S)) (k/(3 CD0))^(1/4); (L/D)max = 1/(2 sqrt(k CD0))",
                "attribution_note": "excrescence/interference and profile drag at the flight Re; LSB penalty if flight Re < design Re"})
    # R2 static margin
    sm_d, sm_i = -design["Cma"] / design["CLa"], -ident["Cma"] / ident["CLa"]
    dsm = sm_i - sm_d
    r2 = {"rule": "R2 static margin", "discrepancy": f"SM {100 * dsm:+.2f}% MAC vs design (SM = -Cma/CLa)", "value_mac": dsm,
          "delta": {"cg_forward_m": max(-dsm, 0.0) * ac["c"], "cg_forward_mac": max(-dsm, 0.0)},
          "formula": "move CG forward by -dSM x MAC (first-order: Xnp fixed)",
          "attribution_note": "Cma scales with assumed Iyy; a 10% Iyy error moves identified SM by 10% of SM"}
    if tail:
        # dXnp/dVh ~ eta_t (a_t/a) (1 - de/da) (Etkin and Reid 1996, eq. 2.3,13), first-order
        sens = tail["eta_t"] * tail["a_t"] / design["CLa"] * (1 - tail["deda"])
        r2["delta"]["tail_volume_increase"] = max(-dsm, 0.0) / sens
        r2["delta"]["dXnp_dVh"] = sens
        r2["formula"] += "; or raise V_h by -dSM / (dXnp/dV_h), dXnp/dV_h = eta_t (a_t/a)(1 - d eps/d alpha)"
    out.append(r2)
    # R3 induced drag factor
    d_k = 100 * (ident["k"] / design["k"] - 1)
    out.append({"rule": "R3 k (span efficiency)", "discrepancy": f"k {d_k:+.1f}% vs design", "value_pct": d_k,
                "delta": {"note": "span-load / twist review; trim-drag estimate at the measured CG"}, "formula": "k = 1/(pi e AR)"})
    # R4 trim offset (Cm0)
    dcm0 = ident["Cm0"] - design["Cm0"]
    de_eq = -dcm0 / ident["Cmde"]
    out.append({"rule": "R4 trim elevator offset", "discrepancy": f"Cm0 {dcm0:+.4f}", "value": dcm0,
                "delta": {"elevator_equiv_deg": math.degrees(de_eq), "tail_incidence_change_deg": math.degrees(de_eq) * tau_e},
                "formula": f"d i_t = tau_e x (-dCm0 / Cmde), tau_e = {tau_e} (elevator effectiveness, assumed)"})
    out.append({"rule": "R5 prop re-match", "discrepancy": "n/a (blind demo flights are unpowered)", "delta": {},
                "formula": "T = D + W sin(gamma) + m dV/dt, eta = T V / (volt x curr)"})
    return out
