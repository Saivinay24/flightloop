"""Post-process one SU2 run directory: convergence status, coefficients, drag split, surface Cp/Cf/y+, LSB.

Status rules (PREREG item 6): converged if rms[P] dropped >= 6 orders from its first value, or Cl and Cd max
relative deviation over the last 1000 iterations < 1e-5; else `unsteady` (mean and amplitude of Cl over the last
2000 iterations recorded) once the iteration cap is reached; `running`/`failed` otherwise.
Cm is integrated from surface Cp with the TM-4062 definition (about c/4, nose-up positive) as an independent check
of SU2's CMz sign convention.
"""
from __future__ import annotations

import json
import pathlib
import re

import numpy as np
import pandas as pd

from .vtu import read_vtu
from ..geometry import Airfoil, read_selig, with_le_point
from ..config import ROOT


def read_history(d: pathlib.Path) -> pd.DataFrame:
    parts = sorted(d.glob("history_part*.csv"), key=lambda p: int(re.findall(r"\d+", p.stem)[-1])) + \
        ([d / "history.csv"] if (d / "history.csv").exists() else [])
    frames = []
    offset = 0
    for p in parts:
        df = pd.read_csv(p)
        df.columns = [c.strip().strip('"') for c in df.columns]
        if df.empty:
            continue
        df["iter_global"] = df["Inner_Iter"] - df["Inner_Iter"].iloc[0] + offset
        offset = int(df["iter_global"].iloc[-1]) + 1
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def breakdown(d):
    out = {}
    p = d / "forces_breakdown.dat"
    if not p.exists():
        return out
    for line in p.read_text().splitlines():
        m = re.match(r"Total (CL|CD|CMz):\s+(-?[\d.eE+-]+) \| Pressure \(\s*-?\d+%\):\s+(-?[\d.eE+-]+) \| Friction \(\s*-?\d+%\):\s+(-?[\d.eE+-]+)", line)
        if m:
            out[m.group(1)] = {"total": float(m.group(2)), "pressure": float(m.group(3)), "friction": float(m.group(4))}
    return out


def status_from_history(h: pd.DataFrame, iter_cap: int):
    if h.empty:
        return {"state": "failed", "reason": "no history"}
    rp = h["rms[P]"].to_numpy()
    cl, cd = h["CL"].to_numpy(), h["CD"].to_numpy()
    n = len(h)
    drop = float(rp[-1] - rp[: min(20, n)].max())
    st = {"iters": int(n), "rms_p_drop_orders": -drop}
    if n >= 1000:
        w = slice(n - 1000, n)
        st["cl_maxreldev_1000"] = float(np.max(np.abs(cl[w] - cl[-1])) / max(abs(cl[-1]), 1e-9))
        st["cd_maxreldev_1000"] = float(np.max(np.abs(cd[w] - cd[-1])) / max(abs(cd[-1]), 1e-9))
    res_ok = -drop >= 6.0
    cau_ok = n >= 1000 and st["cl_maxreldev_1000"] < 1e-5 and st["cd_maxreldev_1000"] < 1e-5
    if res_ok or cau_ok:
        st["state"] = "converged"
        st["criterion"] = "residual" if res_ok else "cauchy"
    elif n >= iter_cap:
        st["state"] = "unsteady"
        w = slice(max(0, n - 2000), n)
        st["cl_mean_last2000"] = float(cl[w].mean())
        st["cl_amp_last2000"] = float(0.5 * (cl[w].max() - cl[w].min()))
        st["cd_mean_last2000"] = float(cd[w].mean())
        st["cd_amp_last2000"] = float(0.5 * (cd[w].max() - cd[w].min()))
    else:
        st["state"] = "incomplete"
    st["cl"], st["cd"], st["cmz_su2"] = float(cl[-1]), float(cd[-1]), float(h["CMz"].iloc[-1])
    return st


_AF = {}


def airfoil(name):
    if name not in _AF:
        _AF[name] = Airfoil(with_le_point(read_selig(ROOT / f"data/airfoils/{name}.dat")))
    return _AF[name]


def surface(d: pathlib.Path, af_name: str):
    s = read_vtu(d / "surface_flow.vtu")
    P = s["Points"][:, :2].astype(float)
    cp = s["Pressure_Coefficient"].astype(float)
    cf = s["Skin_Friction_Coefficient"][:, :2].astype(float)
    yp = s["Y_Plus"].astype(float)
    af = airfoil(af_name)
    # classify upper/lower by nearest spline parameter
    tt = np.linspace(0, 1, 8001)
    sp = af.point(tt)
    idx = np.argmin(((P[:, None, :] - sp[None, :, :]) ** 2).sum(-1), axis=1)
    t = tt[idx]
    res = {}
    for side, mask in (("upper", t < af.t_le), ("lower", t >= af.t_le)):
        o = np.argsort(P[mask, 0])
        x, y = P[mask][o, 0], P[mask][o, 1]
        tan = np.gradient(np.stack([x, y], 1), axis=0)
        tan /= np.linalg.norm(tan, axis=1, keepdims=True)   # points aft (x increasing)
        cft = np.sum(cf[mask][o] * tan, 1)
        res[side] = {"x": x, "y": y, "cp": cp[mask][o], "cf_t": cft, "yplus": yp[mask][o]}
    # Cm about c/4 from the closed Cp integral (TM-4062 symbols list definition), nose-up positive
    up, lo = res["upper"], res["lower"]
    X = np.r_[up["x"][::-1], lo["x"]]
    Y = np.r_[up["y"][::-1], lo["y"]]
    C = np.r_[up["cp"][::-1], lo["cp"]]
    xm, ym, cm_ = 0.5 * (X[1:] + X[:-1]), 0.5 * (Y[1:] + Y[:-1]), 0.5 * (C[1:] + C[:-1])
    dx, dy = np.diff(X), np.diff(Y)
    # contour TE -> upper -> LE -> lower -> TE is counter-clockwise; outward normal (dy, -dx) per unit length
    nx, ny = dy, -dx
    fx, fy = -np.sum(cm_ * nx), -np.sum(cm_ * ny)          # force on body, per unit q c
    mz_ccw = np.sum((xm - 0.25) * (-cm_ * ny) - ym * (-cm_ * nx))
    mz = -mz_ccw   # flow along +x: counter-clockwise rotation is nose-down, so nose-up Cm = -Mz_ccw
    res["cp_integral"] = {"cn": float(fy), "cc": float(fx), "cm_c4_noseup": float(mz)}
    return res


def lsb(up, x_min=None):
    """Separation/reattachment from wall-tangential Cf on the upper surface aft of the suction peak."""
    x, cf, cp = up["x"], up["cf_t"], up["cp"]
    i_peak = int(np.argmin(cp))
    lo = max(i_peak, int(np.searchsorted(x, x_min))) if x_min is not None else i_peak
    sel = np.arange(lo, len(x))
    sel = sel[x[sel] < 0.995]
    xs = xr = None
    for k in range(len(sel) - 1):
        i, j = sel[k], sel[k + 1]
        if xs is None and cf[i] > 0 and cf[j] <= 0:
            xs = float(x[i] + (x[j] - x[i]) * cf[i] / (cf[i] - cf[j]))
        elif xs is not None and cf[i] <= 0 and cf[j] > 0:
            xr = float(x[i] + (x[j] - x[i]) * (-cf[i]) / (cf[j] - cf[i]))
            break
    return {"x_s": xs, "x_r": xr, "bubble": xs is not None}


def process(d: pathlib.Path, write=True):
    d = pathlib.Path(d)
    case = json.loads((d / "case.json").read_text())
    h = read_history(d)
    st = status_from_history(h, case["iters"])
    out = {"run_id": case["run_id"], "case": case, "status": st}
    if (d / "surface_flow.vtu").exists() and st["state"] in ("converged", "unsteady", "incomplete"):
        s = surface(d, case["airfoil"])
        out["yplus_max"] = float(max(s["upper"]["yplus"].max(), s["lower"]["yplus"].max()))
        out["lsb_upper"] = lsb(s["upper"])
        out["cp_integral"] = s["cp_integral"]
        out["breakdown"] = breakdown(d)
        np.savez(d / "surface.npz", **{f"{k}_{f}": v for k in ("upper", "lower") for f, v in s[k].items()})
        pd.DataFrame({"x": np.r_[s["upper"]["x"], s["lower"]["x"]], "y": np.r_[s["upper"]["y"], s["lower"]["y"]],
                      "side": ["upper"] * len(s["upper"]["x"]) + ["lower"] * len(s["lower"]["x"]),
                      "cp": np.r_[s["upper"]["cp"], s["lower"]["cp"]], "cf_t": np.r_[s["upper"]["cf_t"], s["lower"]["cf_t"]],
                      "yplus": np.r_[s["upper"]["yplus"], s["lower"]["yplus"]]}).to_csv(d / "surface.csv", index=False)
        if not h.empty:
            h.to_csv(d / "history_all.csv", index=False)
    if write:
        (d / "result.json").write_text(json.dumps(out, indent=1, default=float))
        (d / "status.json").write_text(json.dumps(st, indent=1))
    return out
