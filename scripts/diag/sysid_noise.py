"""Diagnostic: which sensor-noise term biases equation-error CLa? One turbulent flight, noise terms toggled."""
import sys, pathlib, copy
ROOT = pathlib.Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT))
import numpy as np, pandas as pd
from flightloop.config import cfg
from flightloop.sim import truth as TR, maneuvers as M, sensors as SEN
from flightloop.sysid import equation_error as EE
p = TR.load(); T = TR.truth_params(p); AC = {k: p[k] for k in ("S", "c", "b", "mass", "Iyy")}
SC = cfg()["sysid"]
rec = M.fly(p, 1, 1, [14.5, 17.0], alt0=900)
base = SC["sensor_noise"]
for name in ("none", "all", "accel_ms2", "gyro_rads", "baro_m", "pitot_ms", "pitot_scale", "alpha_vane_rad"):
    sc = copy.deepcopy(SC)
    nz = {k: 0.0 for k in base}
    if name == "all":
        nz = dict(base)
    elif name != "none":
        nz[name] = base[name]
    sc["sensor_noise"] = nz
    sen = SEN.apply(rec, sc, np.random.default_rng(3), 1)
    df = pd.DataFrame({k: sen[k][:, 0] for k in ("t", "ax", "az", "q", "tas", "de", "h_baro")})
    r = EE.estimate(df, AC, sen["alpha_vane"][:, 0])
    print(f"{name:15s} CLa {r['CLa']:.3f}  Cma {r['Cma']:.4f}  Cmde {r['Cmde']:.3f}  (truth {T['CLa']:.3f} {T['Cma']:.4f} {T['Cmde']:.3f})")
a = rec["alpha"][:, 0]
X = np.c_[np.ones_like(a), rec["q"][:, 0] * AC["c"] / (2 * rec["tas"][:, 0]), rec["de"][:, 0]]
res = a - X @ np.linalg.lstsq(X, a, rcond=None)[0]
print("alpha std", a.std(), "alpha residual std after [1,qhat,de]", res.std(), "vane noise", base["alpha_vane_rad"])
print("predicted attenuation 1/(1+s_n^2/s_r^2) =", 1 / (1 + base["alpha_vane_rad"] ** 2 / res.var()))
