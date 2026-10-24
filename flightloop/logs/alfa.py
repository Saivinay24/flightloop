"""CMU ALFA processed CSVs (public experimental, CC0) -> canonical schema at 50 Hz.

ROS/mavros conventions are converted: imu/data is FLU body and ENU-referenced orientation; local_position velocity is
ENU. Body FRD: (x, -y, -z); NED velocity: (vy, vx, -vz). Airspeed: mavros nav_info/airspeed 'measured' (ArduPlane
airspeed, EAS); TAS = EAS sqrt(1.225 / rho) with rho from the static pressure (imu/atm_pressure) and an ISA
temperature at that pressure (approximation, labelled). Elevator: rc-out channel 2 PWM (ArduPlane default RC2 =
elevator), mapped to (PWM - 1500) / 500 'per unit command'. Throttle: vfr_hud throttle (0 to 1).
"""
import glob, pathlib
import numpy as np
import pandas as pd
from ..atmos import isa

G = 9.80665


def _read(d, topic):
    f = glob.glob(str(pathlib.Path(d) / f"*-{topic}.csv"))
    if not f:
        return None
    df = pd.read_csv(f[0])
    df["t"] = df["%time"] * 1e-9
    return df


def rho_from_p(p):
    """ISA: invert pressure for altitude, density at that altitude (temperature approximation)."""
    hs = np.linspace(-500, 5000, 2000); a = isa(hs)
    return np.interp(-p, -a["p"], a["rho"])


def load(d, rate=50.0):
    imu = _read(d, "mavros-imu-data"); vel = _read(d, "mavros-local_position-velocity"); alt = _read(d, "mavros-global_position-rel_alt")
    asp = _read(d, "mavros-nav_info-airspeed"); rc = _read(d, "mavros-rc-out"); hud = _read(d, "mavros-vfr_hud")
    pr = _read(d, "mavros-imu-atm_pressure"); fail = _read(d, "failure_status-engines")
    t0 = max(x.t.iloc[0] for x in (imu, vel, alt, asp, rc, hud, pr)); t1 = min(x.t.iloc[-1] for x in (imu, vel, alt, asp, rc, hud, pr))
    grid = np.arange(t0, t1, 1 / rate)
    I = lambda df, c: np.interp(grid, df.t.to_numpy(), df[c].to_numpy(float))
    qx, qy, qz, qw = (I(imu, f"field.orientation.{k}") for k in "xyzw")
    # ENU/FLU quaternion -> roll, pitch, yaw in NED/FRD (roll same, pitch same sign, yaw = 90 deg - yaw_enu)
    roll = np.arctan2(2 * (qw * qx + qy * qz), 1 - 2 * (qx ** 2 + qy ** 2))
    pitch = np.arcsin(np.clip(2 * (qw * qy - qz * qx), -1, 1))
    yaw_enu = np.arctan2(2 * (qw * qz + qx * qy), 1 - 2 * (qy ** 2 + qz ** 2))
    p_st = I(pr, "field.fluid_pressure"); rho = rho_from_p(p_st)
    eas = I(asp, "field.measured")
    out = pd.DataFrame({"t": grid - grid[0],
                        "ax": I(imu, "field.linear_acceleration.x"), "ay": -I(imu, "field.linear_acceleration.y"), "az": -I(imu, "field.linear_acceleration.z"),
                        "p": I(imu, "field.angular_velocity.x"), "q": -I(imu, "field.angular_velocity.y"), "r": -I(imu, "field.angular_velocity.z"),
                        "phi": roll, "theta": pitch, "psi": np.pi / 2 - yaw_enu,
                        "vn": I(vel, "field.twist.linear.y"), "ve": I(vel, "field.twist.linear.x"), "vd": -I(vel, "field.twist.linear.z"),
                        "h_baro": I(alt, "field.data"), "tas": eas * np.sqrt(1.225 / rho), "eas": eas, "rho": rho,
                        "de": (I(rc, "field.channels1") - 1500) / 500.0, "throttle": I(hud, "field.throttle")})
    t_fail = float(fail.t.iloc[0] - grid[0]) if fail is not None and len(fail) else None
    meta = {"source": "alfa", "sequence": pathlib.Path(d).name, "t_engine_failure_s": t_fail,
            "imu_rate_hz": float(1 / np.median(np.diff(imu.t))), "airspeed_rate_hz": float(1 / np.median(np.diff(asp.t))),
            "derived_fields": {"tas": "EAS x sqrt(1.225/rho), rho from static pressure + ISA temperature (approx.)",
                               "de": "rc-out ch2 PWM -> (PWM-1500)/500, per unit command"},
            "license": "CC0 (CMU ALFA, doi 10.1184/R1/12707963)"}
    return out, meta
