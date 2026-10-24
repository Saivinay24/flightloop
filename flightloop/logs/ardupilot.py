"""ArduPilot DataFlash .bin -> canonical schema (pymavlink DFReader). Used on the ALFA dataflash logs (real ArduPlane
flights, public, CC0). IMU (body FRD specific force, rates), ATT (deg), NKF1 (EKF NED velocity), BARO Alt, ARSP
Airspeed (EAS, ArduPilot convention; TAS = EAS x sqrt(1.225/rho(BARO Press))), AETR Elev/Thr (normalized
-4500..4500 / 0..100), AOA (ArduPlane's own wind-triangle estimate; kept as alpha_est, never as a vane)."""
import numpy as np
import pandas as pd
from pymavlink import DFReader
from .alfa import rho_from_p

WANT = {"IMU": ["AccX", "AccY", "AccZ", "GyrX", "GyrY", "GyrZ"], "ATT": ["Roll", "Pitch", "Yaw"], "NKF1": ["VN", "VE", "VD"],
        "BARO": ["Alt", "Press"], "ARSP": ["Airspeed"], "AETR": ["Elev", "Thr"], "AOA": ["AOA"], "BAT": ["Volt", "Curr"], "CURR": ["Volt", "Curr"]}


def parse(path, rate=50.0):
    r = DFReader.DFReader_binary(str(path), zero_time_base=True)
    rows = {k: [] for k in WANT}
    while True:
        m = r.recv_msg()
        if m is None:
            break
        t = m.get_type()
        if t in WANT:
            d = m.to_dict()
            if t == "IMU" and d.get("I", 0) not in (0,):
                continue
            rows[t].append([d["TimeUS"] * 1e-6] + [d.get(c, np.nan) for c in WANT[t]])
    D = {k: pd.DataFrame(v, columns=["t"] + WANT[k]) for k, v in rows.items() if v}
    t0 = max(D[k].t.iloc[0] for k in ("IMU", "ATT", "NKF1", "BARO", "ARSP", "AETR")); t1 = min(D[k].t.iloc[-1] for k in ("IMU", "ATT", "NKF1", "BARO", "ARSP", "AETR"))
    grid = np.arange(t0, t1, 1 / rate)
    I = lambda k, c: np.interp(grid, D[k].t.to_numpy(), D[k][c].to_numpy(float)) if k in D else np.full(len(grid), np.nan)
    rho = rho_from_p(I("BARO", "Press"))
    deg = np.pi / 180
    bat = "BAT" if "BAT" in D else ("CURR" if "CURR" in D else None)
    df = pd.DataFrame({"t": grid - grid[0], "ax": I("IMU", "AccX"), "ay": I("IMU", "AccY"), "az": I("IMU", "AccZ"),
                       "p": I("IMU", "GyrX"), "q": I("IMU", "GyrY"), "r": I("IMU", "GyrZ"),
                       "phi": I("ATT", "Roll") * deg, "theta": I("ATT", "Pitch") * deg, "psi": I("ATT", "Yaw") * deg,
                       "vn": I("NKF1", "VN"), "ve": I("NKF1", "VE"), "vd": I("NKF1", "VD"), "h_baro": I("BARO", "Alt"),
                       "tas": I("ARSP", "Airspeed") * np.sqrt(1.225 / rho), "de": I("AETR", "Elev") / 4500.0, "throttle": I("AETR", "Thr") / 100.0,
                       "volt": I(bat, "Volt") if bat else np.nan, "curr": I(bat, "Curr") if bat else np.nan, "alpha_est": I("AOA", "AOA") * deg})
    return df, {"format": "bin", "de_units": "per unit command (AETR Elev / 4500)", "alpha_note": "alpha_est = ArduPlane AOA estimate, not a vane",
                "duration_s": float(t1 - t0)}
