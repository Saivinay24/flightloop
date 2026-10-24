"""Flight-test schedule and closed-loop flight generation (vectorised over N seeds).

Per airspeed: settle 12 s, steady glide 15 s, 3-2-1-1 (dt = T_sp/3), 2-1-1 (dt = T_sp/2), doublet (dt = T_sp/2),
each followed by 4 s free response. T_sp = StratoGlide AVL short-period period at sea level.
Pitch-attitude-hold glide controller (de = de_trim + K_th (theta - theta_ref) + K_q q) and wings-level roll hold;
excitation inputs are added to the controller output (closed loop, as a real test pilot/autopilot would fly).
"""
from __future__ import annotations

import json
import numpy as np

from . import sixdof as S
from .dryden import Dryden
from ..config import ROOT
from ..atmos import isa

DEG = np.pi / 180


def t_sp():
    n = json.loads((ROOT / "inputs/stratoglide/numbers.json").read_text())["numbers"]
    v = n["p4.short_period_period_s_sl"]
    return float(v["value"] if isinstance(v, dict) else v)


def schedule(speeds, amp_deg=1.5):
    T = t_sp()
    segs = []
    for V in speeds:
        segs += [(V, "settle", 12.0, None), (V, "glide", 15.0, None)]
        for name, pattern, dt in (("3211", [3, -2, 1, -1], T / 3), ("211", [2, -1, -1], T / 2), ("doublet", [1, -1], T / 2)):
            seq = []
            for v in pattern:
                seq += [np.sign(v)] * abs(v)
            segs.append((V, name, len(seq) * dt + 4.0, (np.array(seq), dt, amp_deg)))
    return segs


def fly(p, n, seed, speeds, alt0=700.0, rate=200, log_rate=50, turb=True, vertical_gusts=True, thrust_fn=None,
        steady_wind=(0.0, 0.0, 0.0), level="light"):
    rng = np.random.default_rng(seed)
    from ..config import cfg as _cfg
    _ctl = _cfg()["sysid"].get("controller", {"k_theta": 0.3, "k_q": 0.05})
    k_th, k_q = _ctl["k_theta"], _ctl["k_q"]
    dt = 1.0 / rate
    segs = schedule(speeds)
    rho0 = float(isa([alt0])["rho"][0])
    X, de0 = S.initial_state(p, speeds[0], alt0, rho0, n)
    wind = S.GustWind(n, steady_wind)
    turbm = Dryden(n, dt, rng, level=level, vertical=vertical_gusts, on=turb)
    rec = {k: [] for k in ("t", "ax", "ay", "az", "p", "q", "r", "phi", "theta", "psi", "vn", "ve", "vd", "h", "tas",
                           "alpha", "beta", "de", "da", "dr", "throttle", "wn", "we", "gust_w", "seg", "qdot_true")}
    every = rate // log_rate
    t = 0.0
    k = 0
    seg_id = 0
    for (V, name, dur, exc) in segs:
        nsteps = int(round(dur / dt))
        for j in range(nsteps):
            R, z, va, Vair, alpha, beta, at, wn = S.SG.aero_state(X, wind)
            rho = at["rho"]
            a_tr, de_tr, gam_tr, th_tr = S.trim_glide(p, V, float(np.mean(rho)))
            phi, th, psi = S.q2euler(X[:, 6:10])
            de = de_tr + k_th * (th - th_tr) / DEG + k_q * X[:, 11] / DEG
            if exc is not None:
                seq, pdt, amp = exc
                idx = int((j * dt) // pdt)
                if idx < len(seq):
                    de = de + amp * seq[idx]
            de = np.clip(de, -15, 15)
            da = np.clip(2.0 * phi / DEG + 0.3 * X[:, 9] / DEG, -15, 15)
            dr = np.zeros(n)
            thr = thrust_fn(t, X, Vair) if thrust_fn else None
            ug, wg = turbm.step(np.maximum(Vair, 5.0), z)
            # gusts along the (northward) flight path and vertical (NED down positive => updraft negative)
            wind.gust = np.stack([ug, np.zeros(n), wg], 1)
            Xn, info = S.rk4(X, p, (de, da, dr), wind, dt, thr)
            if k % every == 0:
                F = info["F"]
                m = p["mass"]
                fx = F[:, 0] + (thr if thr is not None else 0.0)
                v_ned = np.einsum("nij,nj->ni", R, X[:, 3:6])
                rec["t"].append(np.full(n, t)); rec["ax"].append(fx / m); rec["ay"].append(F[:, 1] / m); rec["az"].append(F[:, 2] / m)
                rec["p"].append(X[:, 10].copy()); rec["q"].append(X[:, 11].copy()); rec["r"].append(X[:, 12].copy())
                rec["phi"].append(phi); rec["theta"].append(th); rec["psi"].append(psi)
                rec["vn"].append(v_ned[:, 0]); rec["ve"].append(v_ned[:, 1]); rec["vd"].append(v_ned[:, 2])
                rec["h"].append(z.copy()); rec["tas"].append(Vair.copy()); rec["alpha"].append(alpha.copy()); rec["beta"].append(beta.copy())
                rec["de"].append(de * DEG); rec["da"].append(da * DEG); rec["dr"].append(dr)
                rec["throttle"].append(np.zeros(n) if thr is None else np.asarray(thr) * 0 + 1.0)
                rec["wn"].append(wind.steady[:, 0].copy()); rec["we"].append(wind.steady[:, 1].copy())
                rec["gust_w"].append(wg.copy()); rec["seg"].append(np.full(n, seg_id)); rec["qdot_true"].append(info["qdot"].copy())
            X = Xn
            t += dt
            k += 1
        seg_id += 1
    out = {kk: np.array(v) for kk, v in rec.items()}
    out["segments"] = [(V, name, dur) for (V, name, dur, _) in segs]
    return out
