"""P5 blind loop-closure demo on >= 5 seeds (synthetic). Per seed:
1. stubs/inject.py draws a discrepancy, applies it to the design model, writes runs/<id>/seed_k/SEALED_truth.json.
2. Flight 1 (truth model, turbulence, sensor noise) -> canonical Parquet logs.
3. Engine (flightloop/loop/engine.py) sees only those logs + the design model; every file it opens is logged by a
   Python audit hook to runs/<id>/engine_file_access.log.
4. Flight 2 (truth model, new seed) = the next flight. Sink rate and trim elevator per speed are measured from its
   steady-glide segments and compared with predictions from the design model (before) and the updated model (after).
Score (PREREG choice 10): reduction = 1 - |err_after| / |err_before| (mean |err| over speeds), median over seeds.
"""
import json, sys, pathlib, time
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import pandas as pd
from flightloop.config import cfg, run_id as make_rid, put_number
from flightloop.sim import truth as TR, maneuvers as M, sensors as SEN
from flightloop.logs import schema
from flightloop.loop import engine as E
from flightloop.atmos import isa
sys.path.insert(0, str(ROOT / "stubs"))
import inject  # noqa: E402

SC = cfg()["sysid"]
SCRIPT = "scripts/run_loop_blind.py"
DESIGN = ROOT / "inputs/stratoglide/derivs.json"
SEEDS = [1, 2, 3, 4, 5, 6]
LOGGING = {"on": False, "fh": None}


def audit(event, args):
    if LOGGING["on"] and event == "open" and args and isinstance(args[0], (str, bytes, pathlib.Path)):
        LOGGING["fh"].write(f"{args[0]}\n")


def write_flight(p, seed, d, tag):
    rec = M.fly(p, 1, seed, SC["speeds_ms"], alt0=SC["start_alt_m"], rate=SC["sim_rate_hz"], log_rate=SC["log_rate_hz"])
    sen = SEN.apply(rec, SC, np.random.default_rng(seed + 7), 1)
    d.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame({k: sen[k][:, 0] for k in ("t", "ax", "ay", "az", "p", "q", "r", "phi", "theta", "psi", "vn", "ve", "vd",
                                                "h_baro", "tas", "wn", "we", "de", "da", "dr", "throttle")})
    df["alpha"] = sen["alpha_vane"][:, 0]; df["beta"] = rec["beta"][:, 0]; df["maneuver_id"] = rec["seg"][:, 0]
    schema.write(df, d / f"{tag}.parquet", {"source": "synthetic", "aircraft": "StratoGlide design model + injected discrepancy (blind)",
                                             "seed": seed, "segments": rec["segments"], "derived_fields": {"alpha": "vane (simulated)"}})
    return df, rec["segments"]


def measured_glide(df, segs):
    mid = df["maneuver_id"].to_numpy(); out = []
    for i, sg in enumerate(segs):
        if sg[1] != "glide":
            continue
        m = mid == i
        h = df["h_baro"].to_numpy()[m]
        a, th, V = df["alpha"].to_numpy()[m], df["theta"].to_numpy()[m], df["tas"].to_numpy()[m]
        out.append({"V": float(V.mean()), "rho": float(isa([h.mean()])["rho"][0]),
                    "sink": float(df["vd"][m].mean()), "de_deg": float(np.degrees(df["de"][m].mean())),
                    # secondary, post hoc (HONEST_LOG 14): air-relative sink V sin(alpha - theta) from measured TAS, vane alpha
                    # and pitch attitude (wings level), which removes the air-mass vertical motion that inertial vd contains
                    "sink_air": float(np.mean(V * np.sin(a - th))),
                    "inside_envelope": bool(np.degrees(a).max() <= E.ALPHA_ENV_DEG)})
    return out


def main():
    rid = make_rid("loop", "blind6")
    base = ROOT / "runs" / rid; base.mkdir(parents=True, exist_ok=True)
    LOGGING["fh"] = open(base / "engine_file_access.log", "w", buffering=1)
    sys.addaudithook(audit)
    p_design = TR.load()
    per = []
    t0 = time.time()
    for s in SEEDS:
        sd = base / f"seed_{s}"
        d = inject.draw(s)
        p_true = inject.apply(p_design, d)
        sd.mkdir(parents=True, exist_ok=True)
        inject.seal(sd / "SEALED_truth.json", s, d, p_true)
        write_flight(p_true, 100 + s, sd / "flight1", "flight1")
        LOGGING["on"] = True
        res, p_d, p_u = E.run(sd / "flight1", DESIGN, sd / "engine")
        LOGGING["on"] = False
        df2, segs2 = write_flight(p_true, 500 + s, sd / "flight2", "flight2")
        meas = measured_glide(df2, segs2)
        eb_s, ea_s, eb_d, ea_d = [], [], [], []
        sec = {"air_all": ([], []), "air_env": ([], []), "trim_env": ([], [])}
        for mm in meas:
            sb, db = E.predict(p_d, mm["V"], mm["rho"]); sa, da = E.predict(p_u, mm["V"], mm["rho"])
            mm.update(pred_before={"sink": sb, "de_deg": db}, pred_after={"sink": sa, "de_deg": da})
            eb_s.append(abs(sb - mm["sink"])); ea_s.append(abs(sa - mm["sink"])); eb_d.append(abs(db - mm["de_deg"])); ea_d.append(abs(da - mm["de_deg"]))
            sec["air_all"][0].append(abs(sb - mm["sink_air"])); sec["air_all"][1].append(abs(sa - mm["sink_air"]))
            if mm["inside_envelope"]:
                sec["air_env"][0].append(abs(sb - mm["sink_air"])); sec["air_env"][1].append(abs(sa - mm["sink_air"]))
                sec["trim_env"][0].append(abs(db - mm["de_deg"])); sec["trim_env"][1].append(abs(da - mm["de_deg"]))
        sm_true = -d["cg_aft_mac"]
        r = {"seed": s, "cd0_delta_true_pct": d["cd0_pct"], "cd0_delta_rec_pct": res["cd0_delta_pct"],
             "cd0_delta_err_pp": res["cd0_delta_pct"] - d["cd0_pct"], "sm_delta_true_mac": sm_true, "sm_delta_rec_mac": res["sm_delta_mac"],
             "sm_delta_err_mac_pct": 100 * (res["sm_delta_mac"] - sm_true),
             "sink_err_before": float(np.mean(eb_s)), "sink_err_after": float(np.mean(ea_s)),
             "trim_err_before_deg": float(np.mean(eb_d)), "trim_err_after_deg": float(np.mean(ea_d)), "per_speed": meas}
        r["sink_err_reduction_pct"] = 100 * (1 - r["sink_err_after"] / r["sink_err_before"])
        r["trim_err_reduction_pct"] = 100 * (1 - r["trim_err_after_deg"] / r["trim_err_before_deg"])
        r["secondary_post_hoc"] = {k: {"n_glides": len(b), "err_before": float(np.mean(b)), "err_after": float(np.mean(a_)),
                                       "reduction_pct": float(100 * (1 - np.mean(a_) / np.mean(b)))} for k, (b, a_) in sec.items() if b}
        per.append(r)
        print(json.dumps({k: v for k, v in r.items() if k != "per_speed"}, default=float), flush=True)
    LOGGING["fh"].close()
    acc = (base / "engine_file_access.log").read_text()
    summ = {"run_id": rid, "n_seeds": len(per), "access_log": str((base / "engine_file_access.log").relative_to(ROOT)),
            "sealed_opened_by_engine": "SEALED_truth" in acc,
            "max_cd0_delta_err_pp": float(max(abs(r["cd0_delta_err_pp"]) for r in per)),
            "max_sm_delta_err_mac_pct": float(max(abs(r["sm_delta_err_mac_pct"]) for r in per)),
            "median_sink_err_reduction_pct": float(np.median([r["sink_err_reduction_pct"] for r in per])),
            "median_trim_err_reduction_pct": float(np.median([r["trim_err_reduction_pct"] for r in per])),
            "wall_s": time.time() - t0, "seeds": per}
    summ["secondary_post_hoc_median_reduction_pct"] = {k: float(np.median([r["secondary_post_hoc"][k]["reduction_pct"] for r in per if k in r["secondary_post_hoc"]]))
                                                      for k in ("air_all", "air_env", "trim_env")}
    summ["secondary_note"] = ("post hoc (HONEST_LOG 14): air_all = sink measured as V sin(alpha - theta); air_env and trim_env = same, "
                              "restricted to flight-2 glides with alpha <= 9 deg (identification envelope). V18 is scored on the primary metric.")
    out = ROOT / "results/loop"; out.mkdir(parents=True, exist_ok=True)
    (out / "blind_summary.json").write_text(json.dumps(summ, indent=1, default=float))
    src = "synthetic (blind injected discrepancies)"
    put_number("loop.blind.n_seeds", len(per), "-", SCRIPT, src, rid)
    put_number("loop.blind.max_cd0_delta_err_pp", summ["max_cd0_delta_err_pp"], "pp", SCRIPT, src, rid)
    put_number("loop.blind.max_sm_delta_err_mac_pct", summ["max_sm_delta_err_mac_pct"], "% MAC", SCRIPT, src, rid)
    put_number("loop.blind.sink_err_reduction_pct", summ["median_sink_err_reduction_pct"], "%", SCRIPT, src, rid, "median over seeds")
    put_number("loop.blind.trim_err_reduction_pct", summ["median_trim_err_reduction_pct"], "%", SCRIPT, src, rid, "median over seeds")
    for k, v in summ["secondary_post_hoc_median_reduction_pct"].items():
        put_number(f"loop.blind.secondary.{k}_reduction_pct", v, "%", SCRIPT, src, rid, "post hoc secondary metric, median over seeds (HONEST_LOG 14)")
    put_number("loop.blind.error_reduction_pct", min(summ["median_sink_err_reduction_pct"], summ["median_trim_err_reduction_pct"]), "%", SCRIPT, src, rid,
               "smaller of the sink-rate and trim-elevator median reductions")
    print(json.dumps({k: v for k, v in summ.items() if k != "seeds"}, indent=1, default=float))


if __name__ == "__main__":
    main()
