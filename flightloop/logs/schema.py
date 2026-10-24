"""Canonical flight-log schema: one Parquet per flight, 50 Hz, SI, body FRD; metadata JSON alongside."""
from __future__ import annotations

import json
import pathlib
import numpy as np
import pandas as pd

COLUMNS = ["t", "ax", "ay", "az", "p", "q", "r", "phi", "theta", "psi", "vn", "ve", "vd", "h_baro", "tas", "wn", "we",
           "de", "da", "dr", "throttle", "volt", "curr", "alpha", "beta"]
REQUIRED = ["t", "ax", "ay", "az", "p", "q", "r", "phi", "theta", "psi", "vn", "ve", "vd", "h_baro", "tas", "de", "throttle"]
SOURCES = {"synthetic", "alfa", "px4_public", "ardupilot_sitl", "ardupilot_public"}


def write(df: pd.DataFrame, path, meta: dict):
    path = pathlib.Path(path)
    assert meta.get("source") in SOURCES, meta.get("source")
    for c in COLUMNS:
        if c not in df.columns:
            df[c] = np.nan
    df = df[COLUMNS + [c for c in df.columns if c not in COLUMNS]]
    df.to_parquet(path, index=False)
    path.with_suffix(".json").write_text(json.dumps(meta, indent=1, default=str))
    return path


def read(path):
    path = pathlib.Path(path)
    return pd.read_parquet(path), json.loads(path.with_suffix(".json").read_text())


def validate(df: pd.DataFrame, rate_hz=50.0):
    """Schema checks used by verify.py V14."""
    issues = []
    for c in REQUIRED:
        if c not in df.columns:
            issues.append(f"missing {c}")
        elif df[c].isna().all():
            issues.append(f"all-NaN {c}")
    t = df["t"].to_numpy()
    if not np.all(np.diff(t) > 0):
        issues.append("t not monotone")
    dt = np.median(np.diff(t))
    if abs(dt - 1 / rate_hz) > 1e-6:
        issues.append(f"rate {1/dt:.3f} Hz != {rate_hz}")
    return issues


def resample(df: pd.DataFrame, rate_hz=50.0):
    """Linear resample of every numeric column onto a uniform grid from t[0]."""
    t = df["t"].to_numpy(float)
    order = np.argsort(t); t = t[order]
    keep = np.r_[True, np.diff(t) > 0]
    t = t[keep]
    grid = np.arange(t[0], t[-1], 1.0 / rate_hz)
    out = {"t": grid - grid[0]}
    for c in df.columns:
        if c == "t":
            continue
        v = pd.to_numeric(df[c], errors="coerce").to_numpy(float)[order][keep]
        good = ~np.isnan(v)
        out[c] = np.interp(grid, t[good], v[good]) if good.sum() > 1 else np.full(len(grid), np.nan)
    return pd.DataFrame(out)
