"""P8: ALFA (public experimental, CC0) glide drag polar from the unpowered segments after engine failure.

Method: energy method (no alpha vane on ALFA). Window = 3 s pieces of the post-failure segment (first 2 s skipped as
transient) with throttle < 0.05 and |bank| < 20 deg. Per window: D = -W dE/dt / V (dE/dt = least-squares slope of
E = h + V^2/(2g) over the window), CL = W n_z / (q S) with n_z = -a_z/g from the IMU. Fit CD = CD0 + k CL^2 pooled
over flights; 95% CI by bootstrap over flights (flight = block); plus a systematic interval from the mass (+/-15%)
and area (+/-5%) uncertainty. Approximations labelled: still air in the vertical (vertical gusts bias the energy
method, shown on synthetic data), stopped/windmilling propeller drag is inside CD0 ("airframe + stopped prop").
"""
import glob, json, sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from flightloop.logs import alfa, schema
from flightloop.config import put_number, run_id as make_rid

SCRIPT = "scripts/run_alfa.py"; SRC = "public experimental (CMU ALFA)"
G = 9.80665
AC = {"S": 0.692, "b": 1.98, "mass": 5.4}
AR = AC["b"] ** 2 / AC["S"]


def windows(df, t_fail, rate=50.0):
    out = []
    if t_fail is None:
        return out
    t = df.t.to_numpy()
    start = t_fail + 2.0
    while start + 3.0 <= t[-1]:
        m = (t >= start) & (t < start + 3.0)
        if df.throttle[m].max() < 0.05 and np.degrees(np.abs(df.phi[m])).max() < 20 and df.tas[m].min() > 8:
            out.append(m)
        start += 3.0
    return out


def window_points(df, wins, mass, S):
    W = mass * G; pts = []
    for m in wins:
        t = df.t[m].to_numpy(); V = df.tas[m].to_numpy(); h = df.h_baro[m].to_numpy()
        E = h + V * V / (2 * G)
        dEdt = np.polyfit(t, E, 1)[0]
        Vm = V.mean(); q = 0.5 * df.rho[m].mean() * Vm * Vm
        CD = -W * dEdt / Vm / (q * S)
        CL = W * np.mean(-df.az[m].to_numpy() / G) / (q * S)
        pts.append((CL, CD, Vm))
    return pts


def fit(pts):
    P = np.array(pts); X = np.c_[np.ones(len(P)), P[:, 0] ** 2]
    b, *_ = np.linalg.lstsq(X, P[:, 1], rcond=None)
    r = P[:, 1] - X @ b
    return b, float(np.sqrt(np.mean(r * r)))


def main():
    rid = make_rid("alfa", "glide")
    seqs = sorted(glob.glob(str(ROOT / "data/alfa/processed/*engine_failure*")))
    pq = ROOT / "runs" / rid / "logs"; pq.mkdir(parents=True, exist_ok=True)
    flights, per = [], {}
    for d in seqs:
        df, meta = alfa.load(d)
        schema.write(df.copy(), pq / f"{meta['sequence']}.parquet", meta)
        wins = windows(df, meta["t_engine_failure_s"])
        pts = window_points(df, wins, AC["mass"], AC["S"])
        flights.append({"sequence": meta["sequence"], "t_engine_failure_s": meta["t_engine_failure_s"], "n_windows": len(pts),
                        "imu_rate_hz": meta["imu_rate_hz"], "used": len(pts) >= 1})
        if pts:
            per[meta["sequence"]] = pts
    used = [k for k in per]
    allp = [p for k in used for p in per[k]]
    b, rms = fit(allp)
    rng = np.random.default_rng(0); boots = []
    for _ in range(2000):
        ks = rng.choice(used, len(used), replace=True)
        pp = [p for k in ks for p in per[k]]
        if len(pp) >= 3 and np.ptp([x[0] for x in pp]) > 0.05:
            boots.append(fit(pp)[0])
    B = np.array(boots)
    ci = {"CD0": np.percentile(B[:, 0], [2.5, 97.5]).tolist(), "k": np.percentile(B[:, 1], [2.5, 97.5]).tolist()}
    # systematic: mass +/-15%, S +/-5% (CL and CD both scale with m/S for the energy method: CD ~ m/S, CL ~ m/S)
    sysr = {}
    for fm, fs in ((0.85, 1.0), (1.15, 1.0), (1.0, 0.95), (1.0, 1.05)):
        f = fm / fs   # CL and CD both scale with m/S: CD0' = f CD0, k' = k / f
        sysr[f"m{fm}_S{fs}"] = {"CD0": b[0] * f, "k": b[1] / f}
    e_flight = 1 / (np.pi * AR * b[1]) if b[1] > 0 else None
    avl = {"CLa_per_rad": (0.68670 - 0.36866) / np.radians(4.0), "e_inviscid": 0.9795, "k_inviscid": 1 / (np.pi * 0.9795 * AR)}
    summ = {"run_id": rid, "aircraft": AC, "n_flights_available": len(seqs), "n_flights_used": len(used), "n_segments_used": len(allp),
            "flights": flights, "method": "energy method on post-engine-failure glide windows (3 s)",
            "estimates": [{"name": "CD0 (airframe + stopped prop)", "value": float(b[0]), "ci95": ci["CD0"], "resid_rms": rms,
                           "systematic_range": [min(v["CD0"] for v in sysr.values()), max(v["CD0"] for v in sysr.values())]},
                          {"name": "k", "value": float(b[1]), "ci95": ci["k"], "resid_rms": rms,
                           "systematic_range": [min(v["k"] for v in sysr.values()), max(v["k"] for v in sysr.values())]},
                          {"name": "Oswald e (from k, AR 5.67)", "value": float(e_flight) if e_flight else None,
                           "ci95": sorted([1 / (np.pi * AR * x) for x in ci["k"]]) if min(ci["k"]) > 0 else None, "resid_rms": rms}],
            "avl_vs_flight": [{"quantity": "k (induced drag factor)", "avl": avl["k_inviscid"], "flight": float(b[1]), "flight_ci95": ci["k"]},
                              {"quantity": "span efficiency e", "avl": avl["e_inviscid"], "flight": float(e_flight) if e_flight else None},
                              {"quantity": "CD0", "avl": "n/a (AVL is inviscid; CDvis 0.025 is an input, not a result)", "flight": float(b[0]), "flight_ci95": ci["CD0"]},
                              {"quantity": "CLa (per rad)", "avl": avl["CLa_per_rad"], "flight": "n/a (no alpha measurement on ALFA)"}],
            "approximations": ["mass 5.4 kg +/-15% and S 0.692 m^2 +/-5% are retailer/review figures, not ALFA measurements",
                               "TAS from EAS with rho from static pressure and ISA temperature", "vertical air motion assumed zero (energy method)",
                               "stopped or windmilling propeller drag is included in CD0", "AVL geometry approximate (data/alfa/AIRCRAFT.md)",
                               f"mavros imu/data in the processed CSVs is logged at a median {np.median([f['imu_rate_hz'] for f in flights]):.1f} Hz, so CL (from az) is a mean over few IMU samples per 3 s window"],
            "points": {k: v for k, v in per.items()}}
    (ROOT / "results/sysid").mkdir(exist_ok=True)
    (ROOT / "results/sysid/alfa_summary.json").write_text(json.dumps(summ, indent=1, default=float))
    put_number("sysid.alfa.n_flights_used", len(used), "-", SCRIPT, SRC, rid, f"of {len(seqs)} engine-failure sequences")
    put_number("sysid.alfa.n_segments_used", len(allp), "-", SCRIPT, SRC, rid)
    put_number("sysid.alfa.cd0", float(b[0]), "-", SCRIPT, SRC, rid, "airframe + stopped prop, energy method")
    put_number("sysid.alfa.k", float(b[1]), "-", SCRIPT, SRC, rid)
    put_number("sysid.alfa.cd0_ci_lo", float(ci["CD0"][0]), "-", SCRIPT, SRC, rid, "95% CI, bootstrap over flights")
    put_number("sysid.alfa.cd0_ci_hi", float(ci["CD0"][1]), "-", SCRIPT, SRC, rid, "95% CI, bootstrap over flights")
    put_number("sysid.alfa.k_ci_lo", float(ci["k"][0]), "-", SCRIPT, SRC, rid, "95% CI, bootstrap over flights")
    put_number("sysid.alfa.k_ci_hi", float(ci["k"][1]), "-", SCRIPT, SRC, rid, "95% CI, bootstrap over flights")
    put_number("sysid.alfa.e_flight", float(e_flight) if e_flight else "n/a (k <= 0)", "-", SCRIPT, SRC, rid, "1/(pi AR k), includes all lift-dependent drag")
    put_number("sysid.alfa.e_avl", avl["e_inviscid"], "-", SCRIPT, "AVL (approximate geometry)", rid)
    put_number("sysid.alfa.k_avl", avl["k_inviscid"], "-", SCRIPT, "AVL (approximate geometry)", rid)
    print(json.dumps({k: v for k, v in summ.items() if k not in ("points", "flights")}, indent=1, default=float))


if __name__ == "__main__":
    main()
