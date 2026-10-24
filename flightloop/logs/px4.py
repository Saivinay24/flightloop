"""PX4 ULog -> canonical schema (pyulog). Public PX4 Flight Review logs; parser verification only.

Elevator: the servo whose control-allocation pitch torque CA_SV_CSi_TRQ_P has the largest magnitude, as a normalized
command in [-1, 1] (no surface calibration: 'de' is labelled 'per unit command'). Throttle: actuator_motors control[0].
TAS: airspeed_validated.true_airspeed_m_s (else airspeed.true_airspeed_m_s). Attitude from the vehicle_attitude
quaternion (FRD body, NED). h_baro from vehicle_air_data.baro_alt_meter.
"""
import numpy as np
import pandas as pd
from pyulog import ULog
from . import schema

TOPICS = ["sensor_combined", "vehicle_attitude", "vehicle_local_position", "vehicle_air_data", "airspeed_validated", "airspeed",
          "actuator_servos", "actuator_motors", "battery_status"]


def _df(u, name):
    for d in u.data_list:
        if d.name == name and d.multi_id == 0:
            return pd.DataFrame(d.data).assign(t=lambda x: x.timestamp * 1e-6)
    return None


def q2e(q0, q1, q2, q3):
    phi = np.arctan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1 ** 2 + q2 ** 2))
    th = np.arcsin(np.clip(2 * (q0 * q2 - q3 * q1), -1, 1))
    psi = np.arctan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2 ** 2 + q3 ** 2))
    return phi, th, psi


def parse(path, rate=50.0):
    u = ULog(str(path), TOPICS)
    sc = _df(u, "sensor_combined")
    t0, t1 = sc.t.iloc[0], sc.t.iloc[-1]
    grid = np.arange(t0, t1, 1 / rate)
    out = {"t": grid}
    I = lambda df, col: np.interp(grid, df.t.to_numpy(), df[col].to_numpy(float)) if df is not None and col in df else np.full(len(grid), np.nan)
    for i, c in enumerate(("ax", "ay", "az")):
        out[c] = I(sc, f"accelerometer_m_s2[{i}]")
    for i, c in enumerate(("p", "q", "r")):
        out[c] = I(sc, f"gyro_rad[{i}]")
    att = _df(u, "vehicle_attitude")
    qq = [I(att, f"q[{i}]") for i in range(4)]
    out["phi"], out["theta"], out["psi"] = q2e(*qq)
    lp = _df(u, "vehicle_local_position")
    out["vn"], out["ve"], out["vd"] = I(lp, "vx"), I(lp, "vy"), I(lp, "vz")
    out["h_baro"] = I(_df(u, "vehicle_air_data"), "baro_alt_meter")
    asv = _df(u, "airspeed_validated")
    if asv is not None:
        out["tas"] = I(asv, "true_airspeed_m_s")
    else:
        out["tas"] = I(_df(u, "airspeed"), "true_airspeed_m_s")
    pr = u.initial_parameters
    n = int(pr.get("CA_SV_CS_COUNT", 0))
    elev = max(range(n), key=lambda i: abs(pr.get(f"CA_SV_CS{i}_TRQ_P", 0.0))) if n else None
    sv = _df(u, "actuator_servos")
    out["de"] = I(sv, f"control[{elev}]") if elev is not None else np.full(len(grid), np.nan)
    out["throttle"] = I(_df(u, "actuator_motors"), "control[0]")
    bs = _df(u, "battery_status")
    out["volt"], out["curr"] = I(bs, "voltage_v"), I(bs, "current_a")
    df = pd.DataFrame(out)
    df["t"] = df["t"] - df["t"].iloc[0]
    meta = {"format": "ulog", "elevator_servo_index": elev, "de_units": "per unit command (normalized servo)",
            "sys_name": u.msg_info_dict.get("sys_name"), "ver_sw_release": str(u.msg_info_dict.get("ver_sw_release")),
            "airframe": pr.get("SYS_AUTOSTART"), "duration_s": float(t1 - t0)}
    return df, meta
