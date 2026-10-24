"""Build site/data/bundle.js from pipeline outputs. The site renders only what is in this bundle; no number is typed
into site/index.html or site/assets/app.js. Run after every pipeline stage that changes results:

    .venv/bin/python scripts/build_site.py

The bundle is a plain script (window.FL = {...}) so the page works from file:// and under any subpath.
"""
import csv, hashlib, json, pathlib, re, subprocess, sys, datetime
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

R = ROOT / "results"


def J(rel, default=None):
    p = ROOT / rel
    if not p.exists():
        return default
    return json.loads(p.read_text())


def git(*a):
    return subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout.strip()


def rows_csv(rel, cols, filt=None):
    p = ROOT / rel
    if not p.exists():
        return []
    out = []
    with open(p) as fh:
        for r in csv.DictReader(fh):
            if filt and not filt(r):
                continue
            o = {}
            for c in cols:
                v = r.get(c)
                try:
                    o[c] = float(v)
                except (TypeError, ValueError):
                    o[c] = v
            out.append(o)
    return out


def section(md, title):
    """Return the markdown list items of one '## title' section of a markdown file."""
    if not md:
        return []
    m = re.search(rf"^## {re.escape(title)}\s*$(.*?)(?=^## |\Z)", md, re.S | re.M)
    if not m:
        return []
    items, cur = [], None
    for ln in m.group(1).splitlines():
        if ln.startswith("- ") or re.match(r"^\d+\. ", ln):
            if cur:
                items.append(cur)
            cur = re.sub(r"^(- |\d+\. )", "", ln).strip()
        elif cur is not None and ln.strip():
            cur += " " + ln.strip()
    if cur:
        items.append(cur)
    return items


def main():
    numbers = J("results/numbers.json", {})
    nh = hashlib.sha256(json.dumps({k: v["value"] for k, v in numbers.items()}, sort_keys=True).encode()).hexdigest()[:16]
    val = J("results/cfd/e387_validation.json", {})
    reg = (J("results/cfd/registry.json", {}) or {}).get("runs", [])
    keep = ("run_id", "airfoil", "re", "alpha", "model", "level", "tu_pct", "tag_tu", "state", "criterion", "iters", "cl", "cd",
            "cl_amp", "cd_amp", "yplus_max", "lsb", "in_error_stats", "reported")
    blind = J("results/loop/blind_summary.json")
    deltas = None
    if blind:
        s1 = ROOT / "runs" / blind["run_id"] / f"seed_{blind['seeds'][0]['seed']}" / "engine" / "engine_result.json"
        if s1.exists():
            e = json.loads(s1.read_text())
            deltas = {"seed": blind["seeds"][0]["seed"], "run_id": blind["run_id"], "rules": e["rules"], "discrepancy": e["discrepancy"],
                      "tail_assumptions": e.get("tail_assumptions")}
    honest = (ROOT / "HONEST_LOG.md").read_text() if (ROOT / "HONEST_LOG.md").exists() else ""
    limits = (ROOT / "LIMITS.md").read_text() if (ROOT / "LIMITS.md").exists() else ""
    tags = [t for t in git("tag", "--sort=creatordate").splitlines() if t.startswith("flightloop")]
    comp = J("results/compressible/summary.json", {})
    if comp:
        comp = {k: v for k, v in comp.items() if k != "runs"}
    mf = J("results/mf/summary.json", {})
    if mf:
        mf = {k: mf[k] for k in ("run_id", "meta", "cd_mape", "by_band", "test_airfoils", "nf_conf_median", "leakage_assert_passed") if k in mf}
        mf["n_families"] = len((J("results/mf/summary.json") or {}).get("families", []))
    alfa = J("results/sysid/alfa_summary.json")
    if alfa:
        alfa = {k: alfa[k] for k in ("run_id", "aircraft", "n_flights_available", "n_flights_used", "n_segments_used", "method", "estimates",
                                     "avl_vs_flight", "approximations", "points")}
    synth = J("results/sysid/synthetic_summary.json")
    strato = J("results/strato/disagreement.json")
    corr = J("results/loop/correlation_report.json")
    ver = J("results/verify_report.json", {})

    bundle = {
        "meta": {"generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC"), "git_sha": git("rev-parse", "--short", "HEAD"),
                 "tags": tags, "numbers_hash": nh, "n_numbers": len(numbers)},
        "numbers": numbers,
        "verify": ver,
        "cfd": {"validation": {k: v for k, v in val.items() if k != "definitions"}, "definitions": val.get("definitions"),
                "gci": J("results/cfd/gci_e387.json"), "trust": J("results/cfd/trust_map_e387.json"),
                "registry": [{k: e.get(k) for k in keep} for e in reg],
                "limits": section(limits, "Quantified failure regions"), "not_covered": section(limits, "Not covered")},
        "strato": strato,
        "mf": mf,
        "mf_test": rows_csv("results/mf/test_predictions.csv", ["airfoil", "base", "re", "cl", "cd_exp", "cd_xf9", "cd_corr", "cd_nf"]),
        "mf_lofo": rows_csv("results/mf/lofo_predictions.csv", ["base", "family", "re", "cl", "cd_exp", "cd_xf9", "cd_corr", "cd_nf"]),
        "comp": comp,
        "synth": synth,
        "alfa": alfa,
        "loop": blind,
        "deltas": deltas,
        "corr": corr,
        "honest": {"failed": section(honest, "Failed bars (real numbers)"),
                   "deviations": section(honest, "P0 / P1 deviations and fallbacks") + section(honest, "P3 deviations (all made before the first complete synthetic Monte Carlo)")},
        "gaps": [ln[2:].strip() for ln in ((ROOT / "GAPS.md").read_text().splitlines() if (ROOT / "GAPS.md").exists() else []) if ln.startswith("- ")],
    }
    out = ROOT / "site/data"; out.mkdir(parents=True, exist_ok=True)
    txt = "/* generated by scripts/build_site.py; do not edit */\nwindow.FL = " + json.dumps(bundle, separators=(",", ":"), default=float) + ";\n"
    (out / "bundle.js").write_text(txt)
    print(f"site/data/bundle.js {len(txt) / 1024:.0f} KB, numbers {len(numbers)}, hash {nh}")


if __name__ == "__main__":
    main()
