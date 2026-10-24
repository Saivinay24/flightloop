"""Structured C-mesh for a sharp-TE airfoil, written straight to SU2 format.

Own marching generator (not gmsh transfinite) so that the three levels are a *systematic* refinement:
streamwise nodes come from equidistributing one fixed target-spacing function, and wall-normal nodes from one
fixed exponential map d(eta) = R (exp(b eta) - 1) / (exp(b) - 1) whose b is set so the MEDIUM level has the
configured first-cell height. Every level therefore samples the same mapping at N nodes (ratio ~sqrt 2 per
direction). Layers are marched along smoothed normals with tangential-only Laplacian smoothing that grows with
distance from the wall, which keeps cells orthogonal at the wall and prevents crossing far out.
Markers: `airfoil`, `farfield` (outer layer + both outlet columns). The wake cut is a merged node line.
"""
from __future__ import annotations

import json
import numpy as np
from scipy.optimize import brentq

from ..geometry import Airfoil


def _equidistribute(s_fine, f_fine, n):
    """n nodes on s_fine so that each interval holds the same integral of 1/f."""
    w = np.r_[0, np.cumsum(0.5 * (1 / f_fine[1:] + 1 / f_fine[:-1]) * np.diff(s_fine))]
    return np.interp(np.linspace(0, w[-1], n), w, s_fine)


def wrap_nodes(af: Airfoil, Ns: int):
    """Ns nodes from TE (lower, t=1) round the LE to TE (upper, t=0). Target spacing is finer on the upper
    surface, clustered at LE and TE, and nearly uniform over 0.2c-0.8c upper (where the bubble lives)."""
    t = np.linspace(0, 1, 40001)
    xy = af.point(t)
    s = np.r_[0, np.cumsum(np.hypot(*np.diff(xy, axis=0).T))]
    s_le = np.interp(af.t_le, t, s)
    d_le = np.abs(s - s_le)                      # arc distance from LE
    d_te = np.minimum(s, s[-1] - s)              # arc distance from TE
    upper = t < af.t_le
    mid = np.where(upper, 1.0, 1.9)              # lower surface coarser
    f = np.minimum(np.minimum(0.0006 + 0.06 * d_le, 0.0012 + 0.05 * d_te), 0.0045 * mid)
    s_nodes = _equidistribute(s, f, Ns)
    t_nodes = np.interp(s_nodes, s, t)[::-1]     # t=1 (TE lower) -> 0 (TE upper)
    return af.point(t_nodes), t_nodes


def wake_nodes(ds_te: float, Nw: int, W: float):
    x = np.linspace(1.0, 1.0 + W, 20001)
    f = ds_te + 0.12 * (x - 1.0)
    return _equidistribute(x, f, Nw)


def normal_map_beta(R, h0, n_cells):
    g = lambda b: R * (np.exp(b / n_cells) - 1) / (np.exp(b) - 1) - h0
    return brentq(g, 1e-3, 60)


def build(af: Airfoil, Ns, Nw, Nn, R=100.0, W=20.0, beta=None, h0=2e-5, sigma_k=1.0, sigma_far=0.15, blend_d=0.2):
    wrap, t_nodes = wrap_nodes(af, Ns)
    ds_te = float(np.hypot(*(wrap[-1] - wrap[-2])))
    xw = wake_nodes(ds_te, Nw, W)
    wake_up = np.stack([xw[1:], np.zeros(Nw - 1)], 1)
    wake_lo = wake_up[::-1]
    P = np.vstack([wake_lo, wrap, wake_up])      # C-curve, M points
    M = len(P)
    if beta is None:
        beta = normal_map_beta(R, h0, Nn - 1)
    eta = np.linspace(0, 1, Nn)
    d = R * (np.exp(beta * eta) - 1) / (np.exp(beta) - 1)
    X = np.empty((Nn, M, 2))
    X[0] = P

    far_wake = np.r_[P[:Nw - 1, 0] > 1.0 + 3.0, np.zeros(Ns, bool), P[M - Nw + 1:, 0] > 1.0 + 3.0]
    lo_wake = np.r_[np.ones(Nw - 1, bool), np.zeros(M - Nw + 1, bool)]

    def gauss_normals(curve, sig):
        tan = np.gradient(curve, axis=0)
        tan /= np.linalg.norm(tan, axis=1, keepdims=True)
        n = np.stack([-tan[:, 1], tan[:, 0]], 1)
        if sig > 0:
            s_arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(curve, axis=0), axis=1))]
            Wg = np.exp(-((s_arc[:, None] - s_arc[None, :]) / sig) ** 2)
            n = Wg @ n
            n /= np.linalg.norm(n, axis=1, keepdims=True)
        return n

    # straight-ray directions for the outer region: smoothed wall normals with the angle forced monotone along
    # the C-curve, so neighbouring rays are parallel or diverging and can never cross.
    m0 = gauss_normals(P, sigma_far)
    th = np.unwrap(np.arctan2(m0[:, 1], m0[:, 0]))
    th = 0.5 * (np.minimum.accumulate(th) + np.maximum.accumulate(th[::-1])[::-1])
    m_far = np.stack([np.cos(th), np.sin(th)], 1)

    cur = P.copy()
    for k in range(1, Nn):
        dk = 0.5 * (d[k] + d[k - 1])
        n_loc = gauss_normals(cur, min(sigma_k * dk, 0.05))
        w = np.clip(dk / blend_d, 0, 1)
        w = w * w * (3 - 2 * w)
        n = (1 - w) * n_loc + w * m_far
        n /= np.linalg.norm(n, axis=1, keepdims=True)
        cur = cur + (d[k] - d[k - 1]) * n
        X[k] = cur
    return X, {"Nw": Nw, "Ns": Ns, "Nn": Nn, "beta": beta, "h0": float(d[1]), "R": R, "W": W, "ds_te": ds_te}


def quality(X):
    a, b, c, dd = X[:-1, :-1], X[:-1, 1:], X[1:, 1:], X[1:, :-1]
    def cross(u, v):
        return u[..., 0] * v[..., 1] - u[..., 1] * v[..., 0]
    area = 0.5 * (cross(c - a, dd - b))
    edges = [b - a, c - b, dd - c, a - dd]
    L = np.stack([np.linalg.norm(e, axis=-1) for e in edges])
    ar = L.max(0) / L.min(0)
    angs = []
    for i in range(4):
        e1, e2 = -edges[i - 1], edges[i]
        cosv = np.sum(e1 * e2, -1) / (np.linalg.norm(e1, axis=-1) * np.linalg.norm(e2, axis=-1))
        angs.append(np.degrees(np.arccos(np.clip(cosv, -1, 1))))
    angs = np.stack(angs)
    skew = np.maximum((angs.max(0) - 90) / 90, (90 - angs.min(0)) / 90)
    return {"cells": int(area.size), "min_area": float(area.min()), "n_negative": int((area <= 0).sum()),
            "max_aspect": float(ar.max()), "max_skew": float(skew.max()),
            "max_skew_near_wall": float(skew[:20].max())}


def write_su2(X, Nw, Ns, path):
    Nn, M, _ = X.shape
    idx = -np.ones((Nn, M), int)
    pts = []
    for k in range(Nn):
        for i in range(M):
            if k == 0 and i > M - 1 - i and (M - 1 - i) < Nw:   # merged wake cut / TE (upper copy -> lower id)
                idx[k, i] = idx[0, M - 1 - i]
                continue
            idx[k, i] = len(pts)
            pts.append(X[k, i])
    with open(path, "w") as f:
        f.write("NDIME= 2\n")
        f.write(f"NELEM= {(Nn - 1) * (M - 1)}\n")
        e = 0
        for k in range(Nn - 1):
            for i in range(M - 1):
                f.write(f"9 {idx[k, i]} {idx[k, i + 1]} {idx[k + 1, i + 1]} {idx[k + 1, i]} {e}\n")
                e += 1
        f.write(f"NPOIN= {len(pts)}\n")
        for j, p in enumerate(pts):
            f.write(f"{p[0]:.12e} {p[1]:.12e} {j}\n")
        f.write("NMARK= 2\n")
        a0, a1 = Nw - 1, Nw - 1 + Ns - 1
        f.write("MARKER_TAG= airfoil\n")
        f.write(f"MARKER_ELEMS= {a1 - a0}\n")
        for i in range(a0, a1):
            f.write(f"3 {idx[0, i]} {idx[0, i + 1]}\n")
        far = [(idx[Nn - 1, i], idx[Nn - 1, i + 1]) for i in range(M - 1)]
        far += [(idx[k + 1, 0], idx[k, 0]) for k in range(Nn - 1)]
        far += [(idx[k, M - 1], idx[k + 1, M - 1]) for k in range(Nn - 1)]
        f.write("MARKER_TAG= farfield\n")
        f.write(f"MARKER_ELEMS= {len(far)}\n")
        for a, b in far:
            f.write(f"3 {a} {b}\n")
    return len(pts)
