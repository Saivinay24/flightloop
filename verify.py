"""FlightLoop verifier: 20 numbered checks (10 integrity [I], 10 result [R] per PREREG.md) pre-registered in PREREG.md.

Output: one line per check, then `INTEGRITY a/10, RESULTS b/10`. Exit code nonzero if any [I] check fails.
A failed [R] check is acceptable only when HONEST_LOG.md records it with the real number (checked here).
Writes results/verify_report.json. Reads only pipeline artifacts; never recomputes physics.
"""
from __future__ import annotations

import glob
import json
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
R = ROOT / "results"
EM = "—"
rows = []


def J(p):
    p = R / p if not str(p).startswith("/") else pathlib.Path(p)
    return json.loads(p.read_text()) if p.exists() else None


def num(key):
    d = J("numbers.json") or {}
    return d.get(key, {}).get("value")


def add(cid, cls, ok, msg):
    rows.append({"id": cid, "class": cls, "pass": bool(ok), "msg": msg})


def honest_logged(cid):
    t = (ROOT / "HONEST_LOG.md").read_text() if (ROOT / "HONEST_LOG.md").exists() else ""
    return bool(re.search(rf"\b{cid}\b.*\d", t))


def guard(cid, cls, fn):
    try:
        fn()
    except Exception as e:  # missing artifact = check not satisfied
        add(cid, cls, False, f"not available ({type(e).__name__}: {str(e)[:90]})")


# ---------------------------------------------------------------- V01 ISA
def v01():
    from flightloop.atmos import isa
    r11, r20, n18 = isa([11000])["rho"][0], isa([20000])["rho"][0], isa([18000])["nu"][0]
    e = [abs(r11 / 0.3639 - 1), abs(r20 / 0.0889 - 1), abs(n18 / 1.18e-4 - 1)]
    ok = e[0] < 0.005 and e[1] < 0.005 and e[2] < 0.01
    add("V01", "I", ok, f"rho11 {r11:.4f} ({e[0]*100:+.2f}%), rho20 {r20:.4f} ({e[1]*100:+.2f}%), nu18 {n18:.3e} ({e[2]*100:+.2f}%)")


# ---------------------------------------------------------------- V02 provenance
def v02():
    src = (ROOT / "data/SOURCES.md").read_text()
    paths = [m.group(1).strip() for m in re.finditer(r"^\|\s*([^|]+?)\s*\|", src, re.M)]
    pats = [p for p in paths if p not in ("path", "---")]
    missing = []
    for f in sorted((ROOT / "data").rglob("*")):
        if f.is_dir() or f.name in ("SOURCES.md", ".DS_Store") or "/work/" in str(f) or "/pages/" in str(f):
            continue
        rel = str(f.relative_to(ROOT / "data"))
        if not any(rel == p or (p.endswith("*") and rel.startswith(p[:-1])) for p in pats):
            missing.append(rel)
    diff = subprocess.run([sys.executable, str(ROOT / "scripts/diff_transcriptions.py")], capture_output=True, text=True)
    m = re.search(r"unresolved after page-image resolution (\d+)", diff.stdout)
    unresolved = int(m.group(1)) if m else -1
    syn_bad = 0; syn_n = 0
    from flightloop.logs.schema import SOURCES
    for j in glob.glob(str(ROOT / "runs/*/**/*.json"), recursive=True):
        jp = pathlib.Path(j)
        if jp.parent.name in ("engine",) or jp.name in ("case.json", "status.json", "result.json", "summary.json", "config_snapshot.json", "SEALED_truth.json"):
            continue
        meta = json.loads(jp.read_text())
        if not isinstance(meta, dict) or "source" not in meta and "segments" not in meta and "format" not in meta:
            continue
        syn_n += 1
        rel = str(jp.relative_to(ROOT / "runs"))
        must_syn = rel.startswith(("sysid-", "loop-"))
        if meta.get("source") not in SOURCES or (must_syn and meta.get("source") != "synthetic"):
            syn_bad += 1
    ok = not missing and unresolved == 0 and syn_bad == 0
    add("V02", "I", ok, f"data files without SOURCES row: {len(missing)}{' ' + str(missing[:3]) if missing else ''}; TM-4062 unresolved diffs {unresolved}; canonical logs {syn_n} (simulator outputs must say synthetic), mislabelled {syn_bad}")


# ---------------------------------------------------------------- V03 mesh
def v03():
    reg = J("cfd/registry.json")
    runs = [r for r in reg["runs"] if r["reported"]]
    ymax = max(r["yplus_max"] for r in runs)
    nosha = [r["run_id"] for r in runs if not r.get("mesh_sha256_16")]
    ratios = []
    for af in ("e387", "sd7037"):
        m = J(f"mesh_{af}.json")
        if not m:
            continue
        ratios += [(m["medium"]["cells"] / m["coarse"]["cells"]) ** 0.5, (m["fine"]["cells"] / m["medium"]["cells"]) ** 0.5]
    ok = ymax < 1 and not nosha and all(1.3 <= r <= 1.6 for r in ratios)
    add("V03", "I", ok, f"max y+ {ymax:.3f} over {len(runs)} reported runs; per-direction ratios {', '.join(f'{r:.3f}' for r in ratios)}; runs without mesh SHA {len(nosha)}")


# ---------------------------------------------------------------- V04 convergence
def v04():
    reg = J("cfd/registry.json")
    stats = [r for r in reg["runs"] if r["in_error_stats"]]
    bad = [r["run_id"] for r in stats if r["state"] != "converged"]
    uns = [r for r in reg["runs"] if r["state"] == "unsteady"]
    tm = J("cfd/trust_map_e387.json")
    uns_in_tm = all(any(c["run_id"] == r["run_id"] for c in tm["cells"]) for r in uns) if tm else False
    ok = not bad and uns_in_tm
    add("V04", "I", ok, f"{len(stats)} runs in error statistics, non-converged among them {len(bad)}; unsteady runs {len(uns)} (all in Trust Map: {uns_in_tm})")


# ---------------------------------------------------------------- V05 GCI
def v05():
    g = J("cfd/gci_e387.json")
    parts, ok = [], True
    for m in ("SA", "SST_LM"):
        e = g[m]
        if e["cl"]["type"] != "monotone" or e["cd"]["type"] != "monotone":
            ok = False
            parts.append(f"{m}: Cl {e['cl']['type']}, Cd {e['cd']['type']}")
            continue
        okm = e["cl"]["gci_fine_pct"] < 2 and e["cd"]["gci_fine_pct"] < 5 and 0.5 <= e["cd"]["p"] <= 3 and 0.5 <= e["cl"]["p"] <= 3
        ok &= okm
        parts.append(f"{m}: GCI Cl {e['cl']['gci_fine_pct']:.2f}% (p {e['cl']['p']:.2f}), Cd {e['cd']['gci_fine_pct']:.2f}% (p {e['cd']['p']:.2f})")
    add("V05", "R", ok, "; ".join(parts) + " (bars Cl < 2%, Cd < 5%, p in [0.5, 3])")


# ---------------------------------------------------------------- V06-V08 validation
def v06():
    v = J("cfd/e387_validation.json")["re200k_a4"]["SST_LM"]
    ok = v["in_error_stats"] and abs(v["cl_err_pct"]) < 5 and abs(v["cd_err_pct"]) < 20 and v["cp_rms_upper"] < 0.05
    add("V06", "R", ok, f"SST-LM Re 200k a4 (state {v['state']}, in error stats {v['in_error_stats']}): Cl {v['cl_err_pct']:+.1f}%, Cd {v['cd_err_pct']:+.1f}%, RMS dCp upper {v['cp_rms_upper']:.3f} (bars 5%, 20%, 0.05)")


def v07():
    v = J("cfd/e387_validation.json")["re200k_a4"]["SA"]
    ok = v["cd_err_pct"] > 25
    add("V07", "R", ok, f"SA Cd {v['cd_err_pct']:+.1f}% vs TM-4062 (bar > +25%)")


def v08():
    v = J("cfd/e387_validation.json")["lsb_re200k"]
    parts, ok = [], True
    for a in ("0", "4"):
        e = v[a]
        if e["xs_err"] is None or e["xr_err"] is None:
            ok = False; parts.append(f"a{a}: no bubble/reattachment in CFD"); continue
        ok &= bool(e.get("in_error_stats")) and abs(e["xs_err"]) < 0.05 and abs(e["xr_err"]) < 0.08
        parts.append(f"a{a} ({'in' if e.get('in_error_stats') else 'NOT in'} error stats): x_s {e['xs_err']:+.3f}c, x_r {e['xr_err']:+.3f}c")
    add("V08", "R", ok, "; ".join(parts) + " (bars 0.05c, 0.08c)")


# ---------------------------------------------------------------- V09 limits
def v09():
    tm = J("cfd/trust_map_e387.json")
    reg = J("cfd/registry.json")
    covered = {c["run_id"] for c in tm["cells"]}
    allruns = {r["run_id"] for r in reg["runs"] if r["airfoil"] == "e387"}
    lim = (ROOT / "LIMITS.md").read_text()
    regions = [l for l in lim.splitlines() if re.match(r"^\d+\.\s", l) and re.search(r"\d+(\.\d+)?\s*(%|c\b|deg)", l)]
    ok = allruns <= covered and len(regions) >= 3
    add("V09", "I", ok, f"Trust Map covers {len(covered & allruns)}/{len(allruns)} E387 runs; LIMITS.md quantified failure regions {len(regions)}")


# ---------------------------------------------------------------- V10 strato
def v10():
    t = (ROOT / "results/strato/DISAGREEMENT.md").read_text()
    table = [l for l in t.splitlines() if l.startswith("|") and not set(l) <= set("|- ")]
    body = table[1:]
    empty = [l for l in body if re.search(r"\|\s*\|", l)]
    res = set(re.findall(r"Re (\d+)k", t))
    ok = len(body) >= 7 and not empty and {"60", "100"} <= res
    add("V10", "I", ok, f"DISAGREEMENT.md rows {len(body)}, empty cells {len(empty)}, Re covered {sorted(res)}")


# ---------------------------------------------------------------- V11 multi-fidelity
def v11():
    s = J("mf/summary.json")
    corr, raw, nf = s["cd_mape"]["corrected"], s["cd_mape"]["xfoil_n9"], s["cd_mape"]["neuralfoil"]
    leak = s["leakage_assert_passed"]
    ok = leak and corr <= 0.75 * raw and corr < nf
    add("V11", "R", ok and leak, f"Cd MAPE corrected {corr:.1f}% vs raw XFOIL {raw:.1f}% (ratio {corr/raw:.2f}, bar <= 0.75) vs NeuralFoil {nf:.1f}%; leakage assert [I] {'PASS' if leak else 'FAIL'}")
    return leak


# ---------------------------------------------------------------- V12-V13 compressible
def v12():
    c = J("compressible/summary.json")
    parts, ok = [], True
    for M, b0, p0 in (("2", 39.31, 1.707), ("3", 27.38, 2.055)):
        e = c["wedge"][f"M{M}"]["fine"]
        db, dp = e["beta_deg"] - b0, (e["p2_p1"] / p0 - 1) * 100
        ok &= abs(db) < 0.5 and abs(dp) < 1
        parts.append(f"M{M}: beta {e['beta_deg']:.2f} ({db:+.2f} deg), p2/p1 {e['p2_p1']:.4f} ({dp:+.2f}%)")
    add("V12", "R", ok, "; ".join(parts) + " (bars 0.5 deg, 1%)")


def v13():
    c = J("compressible/summary.json")["diamond"]["M2_a2"]
    dcl, dcd = (c["fine"]["cl"] / 0.0810 - 1) * 100, (c["fine"]["cd"] / 0.00863 - 1) * 100
    ok = abs(dcl) < 2 and abs(dcd) < 2 and c["monotone_cl"] and c["monotone_cd"]
    add("V13", "R", ok, f"M2 a2 fine: Cl {c['fine']['cl']:.4f} ({dcl:+.2f}%), Cd {c['fine']['cd']:.5f} ({dcd:+.2f}%), monotone Cl {c['monotone_cl']}, Cd {c['monotone_cd']}")


# ---------------------------------------------------------------- V14 parsers
def v14():
    s = J("sysid/synthetic_summary.json")
    p = J("logs/parse_report.json") or {"files": []}
    px4 = [f for f in p["files"] if f["format"] == "ulog" and f["schema_ok"]]
    ap = [f for f in p["files"] if f["format"] == "bin" and f["schema_ok"]]
    ok = s["parquet_roundtrip_exact"] and len(px4) >= 1 and len(ap) >= 1
    add("V14", "I", ok, f"synthetic Parquet round-trip exact {s['parquet_roundtrip_exact']}; PX4 ULogs to full schema {len(px4)}; ArduPilot .bin to full schema {len(ap)}")


# ---------------------------------------------------------------- V15-V16 synthetic sys-ID
def v15():
    s = J("sysid/synthetic_summary.json")["main"]
    bars = {"CD0": 5, "k": 10, "CLa": 5, "Cma": 10, "Cmde": 10, "Cmq": 25}
    ok = all(s[k]["median_abs_err_pct"] < b for k, b in bars.items())
    add("V15", "R", ok, ", ".join(f"{k} {s[k]['median_abs_err_pct']:.1f}% (<{b})" for k, b in bars.items()) + " median |err|, 20 seeds")


def v16():
    c = J("sysid/synthetic_summary.json")["coverage_pct"]
    ok = all(80 <= c[k] <= 100 for k in ("CD0", "CLa", "Cma"))
    add("V16", "R", ok, ", ".join(f"{k} {c[k]:.0f}%" for k in ("CD0", "CLa", "Cma")) + " of seeds inside 95% CI (bar 80-100%)")


# ---------------------------------------------------------------- V17 real data
def v17():
    a = J("sysid/alfa_summary.json")
    est = a["estimates"]
    no_ci = [e["name"] for e in est if e.get("ci95") is None or e.get("resid_rms") is None]
    ok = a["n_flights_used"] >= 1 and bool(a["flights"]) and not no_ci and a.get("avl_vs_flight") and a.get("approximations")
    add("V17", "I", ok, f"ALFA flights used {a['n_flights_used']}, segments {a['n_segments_used']}; estimates without CI/residual {len(no_ci)}; AVL-vs-flight table {'present' if a.get('avl_vs_flight') else 'missing'}")


# ---------------------------------------------------------------- V18 loop closure
def v18():
    b = J("loop/blind_summary.json")
    acc = (ROOT / b["access_log"]).read_text()
    sealed_opened = "SEALED_truth" in acc
    ok = (b["n_seeds"] >= 5 and b["max_cd0_delta_err_pp"] < 5 and b["max_sm_delta_err_mac_pct"] < 1
          and b["median_sink_err_reduction_pct"] >= 70 and b["median_trim_err_reduction_pct"] >= 70 and not sealed_opened)
    add("V18", "R", ok, f"{b['n_seeds']} blind seeds: CD0 delta err max {b['max_cd0_delta_err_pp']:.2f} pp, SM delta err max {b['max_sm_delta_err_mac_pct']:.2f}% MAC, "
        f"sink err cut {b['median_sink_err_reduction_pct']:.0f}%, trim err cut {b['median_trim_err_reduction_pct']:.0f}%; sealed truth opened by engine [I]: {sealed_opened}")
    return not sealed_opened


# ---------------------------------------------------------------- V19 outward surfaces
def v19():
    need = ["report/FlightLoop_Report.pdf", "report/FlightLoop_Report.md", "site/index.html", "site/data/bundle.js", "README.md", "RESUME_CONTRACT.md"]
    missing = [n for n in need if not (ROOT / n).exists()]
    scan = ["README.md", "RESUME_CONTRACT.md", "report/FlightLoop_Report.md", "LIMITS.md", "results/strato/DISAGREEMENT.md"] + \
        [str(p.relative_to(ROOT)) for p in (ROOT / "site").rglob("*") if p.is_file() and p.suffix in (".html", ".js", ".css", ".json")]
    em = [f for f in scan if (ROOT / f).exists() and EM in (ROOT / f).read_text(errors="ignore")]
    chk = subprocess.run([sys.executable, str(ROOT / "scripts/check_numbers.py")], capture_output=True, text=True)
    m = re.search(r"UNTRACEABLE (\d+)", chk.stdout)
    untr = int(m.group(1)) if m else -1
    ok = not missing and not em and untr == 0
    add("V19", "I", ok, f"missing outward files {missing}; files with U+2014 {em}; numbers in resume/report not equal to numbers.json: {untr}")


# ---------------------------------------------------------------- V20 hygiene
def v20():
    secret_re = r"(AKIA[0-9A-Z]{16}|sk-[A-Za-z0-9]{20,}|ghp_[A-Za-z0-9]{30,}|-----BEGIN (RSA|OPENSSH) PRIVATE KEY|password\s*=\s*['\"][^'\"]+)"
    files = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True).stdout.split()
    skip_bin = (".pdf", ".png", ".parquet", ".npz", ".npy", ".dat", ".su2", ".vtu", ".ulg", ".bin", ".zip", ".tgz")
    secrets, weapons = [], []
    wre = re.compile(r"warhead|seeker|guidance law|homing|ramjet|intake", re.I)
    excl = {"MASTER_PLAN.md", "PREREG.md", "verify.py"}
    for f in files:
        p = ROOT / f
        if f.endswith(skip_bin) or not p.exists():
            continue
        t = p.read_text(errors="ignore")
        if re.search(secret_re, t):
            secrets.append(f)
        if f not in excl and not f.startswith("data/tm4062/tm4062.txt") and wre.search(t):
            weapons.append(f)
    log = subprocess.run(["git", "log", "--format=%an%n%ae%n%B"], cwd=ROOT, capture_output=True, text=True).stdout
    attrib = bool(re.search(r"co-authored-by|claude|anthropic", log, re.I))
    # PREREG committed before any results commit
    first_prereg = subprocess.run(["git", "log", "--diff-filter=A", "--format=%at", "--", "PREREG.md"], cwd=ROOT, capture_output=True, text=True).stdout.split()
    first_res = subprocess.run(["git", "log", "--format=%at", "--reverse", "--", "results/numbers.json"], cwd=ROOT, capture_output=True, text=True).stdout.split()
    order_ok = bool(first_prereg) and (not first_res or int(first_prereg[-1]) <= int(first_res[0]))
    # every outward data label present
    ok = not secrets and not attrib and not weapons and order_ok
    add("V20", "I", ok, f"secrets {len(secrets)}; attribution in git log {attrib}; weapon-term hits {weapons}; PREREG committed before first results commit {order_ok}")


def main():
    leak_ok = True; sealed_ok = True
    for cid, cls, fn in [("V01", "I", v01), ("V02", "I", v02), ("V03", "I", v03), ("V04", "I", v04), ("V05", "R", v05),
                         ("V06", "R", v06), ("V07", "R", v07), ("V08", "R", v08), ("V09", "I", v09), ("V10", "I", v10),
                         ("V11", "R", v11), ("V12", "R", v12), ("V13", "R", v13), ("V14", "I", v14), ("V15", "R", v15),
                         ("V16", "R", v16), ("V17", "I", v17), ("V18", "R", v18), ("V19", "I", v19), ("V20", "I", v20)]:
        guard(cid, cls, fn)
    for r in rows:
        tag = "PASS" if r["pass"] else ("FAIL" if r["class"] == "I" or not honest_logged(r["id"]) else "FAIL (logged)")
        r["tag"] = tag
        print(f"{r['id']} [{r['class']}] {tag:13s} {r['msg']}")
    ni = sum(r["pass"] for r in rows if r["class"] == "I")
    nr = sum(r["pass"] for r in rows if r["class"] == "R")
    unlogged = [r["id"] for r in rows if r["class"] == "R" and not r["pass"] and not honest_logged(r["id"])]
    # denominators are counted from the rows (PREREG.md classes give 10 [I] and 10 [R]; the plan's "12/8" split was a miscount)
    nI, nR = sum(r["class"] == "I" for r in rows), sum(r["class"] == "R" for r in rows)
    print(f"INTEGRITY {ni}/{nI}, RESULTS {nr}/{nR}" + (f"  (R fails not in HONEST_LOG: {unlogged})" if unlogged else ""))
    R.mkdir(exist_ok=True)
    (R / "verify_report.json").write_text(json.dumps({"rows": rows, "integrity": ni, "results": nr, "n_integrity": nI, "n_results": nR, "unlogged_r_fails": unlogged}, indent=1))
    sys.exit(0 if ni == 12 and not unlogged else 1)


if __name__ == "__main__":
    main()
