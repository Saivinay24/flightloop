"""Diagnostic: measured q_dot (filtered, differentiated logs) vs the simulator's true q_dot at log instants."""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT))
import numpy as np
from flightloop.sim import truth as TR, maneuvers as M
from flightloop.sysid.equation_error import _lp, d4
p = TR.load()
rec = M.fly(p, 1, 1, [17.0], alt0=900, turb=False)
q = rec["q"][:, 0]; qt = rec["qdot_true"][:, 0]
for hz in (3.0, 5.0, 8.0):
    a = d4(_lp(q, 50, hz), 0.02); b = _lp(qt, 50, hz)
    print(hz, "Hz: slope d4(lp q) vs lp(true qdot)", np.polyfit(b[20:-20], a[20:-20], 1)[0], "corr", np.corrcoef(a, b)[0, 1])
print("unfiltered: slope gradient(q) vs true", np.polyfit(qt, np.gradient(q, 0.02), 1)[0])
rec = M.fly(p, 1, 1, [17.0], alt0=900, turb=False, log_rate=200)
q = rec["q"][:, 0]; qt = rec["qdot_true"][:, 0]
fwd = (q[1:] - q[:-1]) / 0.005
print("200 Hz: slope fwd diff vs true k1", np.polyfit(qt[:-1], fwd, 1)[0], "vs trapezoid", np.polyfit(0.5 * (qt[1:] + qt[:-1]), fwd, 1)[0])
i = np.argmax(np.abs(fwd - 0.5 * (qt[1:] + qt[:-1])))
print("worst idx", i, fwd[i - 2:i + 3], qt[i - 2:i + 3], rec["de"][i - 2:i + 3, 0])
