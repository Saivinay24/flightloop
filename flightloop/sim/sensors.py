"""Sensor models at 50 Hz: white noise per config, pitot scale error (per flight), baro lag, alpha vane noise."""
import numpy as np


def apply(rec, cfgs, rng, n):
    """rec: dict of true signals, arrays (T, n). Returns noisy dict."""
    nz = cfgs["sensor_noise"]
    out = {}
    T = rec["t"].shape[0]
    g = lambda s: rng.standard_normal((T, n)) * s
    for k in ("ax", "ay", "az"):
        out[k] = rec[k] + g(nz["accel_ms2"])
    for k in ("p", "q", "r"):
        out[k] = rec[k] + g(nz["gyro_rads"])
    scale = 1 + rng.standard_normal(n) * nz["pitot_scale"]
    out["tas"] = rec["tas"] * scale + g(nz["pitot_ms"])
    # first-order baro lag (tau 0.2 s) + noise
    h = rec["h"].copy(); tau = 0.2; dt = rec["t"][1, 0] - rec["t"][0, 0]
    hl = np.empty_like(h); hl[0] = h[0]
    a = dt / (tau + dt)
    for i in range(1, T):
        hl[i] = hl[i - 1] + a * (h[i] - hl[i - 1])
    out["h_baro"] = hl + g(nz["baro_m"])
    out["alpha_vane"] = rec["alpha"] + g(nz.get("alpha_vane_rad", 0.0087))
    for k in ("phi", "theta", "psi"):
        out[k] = rec[k] + g(0.002)      # EKF attitude
    for k in ("vn", "ve", "vd"):
        out[k] = rec[k] + g(0.05)       # EKF velocity
    for k in ("de", "da", "dr", "throttle", "t", "wn", "we"):
        out[k] = rec[k]
    out["pitot_scale_true"] = scale
    return out
