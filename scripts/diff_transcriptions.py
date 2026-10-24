"""Compare pass A (machine text layer / second read) with pass B (manual) for every TM-4062 table used.

Every disagreement must have a row in data/tm4062/RESOLUTIONS.csv (value read from the page image) and the
canonical file data/tm4062/canonical/*.csv must equal the resolved value. Exit 0 only when the unresolved
difference count is zero. Writes data/tm4062/diff_report.json.
"""
import csv, json, pathlib, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = ROOT / "data/tm4062"


def load(p):
    with open(p) as f:
        return list(csv.DictReader(f))


def resolutions():
    p = D / "RESOLUTIONS.csv"
    if not p.exists():
        return {}
    return {(r["file"], int(r["row"]), r["column"]): r["resolved_value"] for r in load(p)}


def same(a, b):
    if a in ("", None) and b in ("", None):
        return True
    try:
        return abs(float(a) - float(b)) < 5e-5
    except (TypeError, ValueError):
        return a == b


def compare(fname, align):
    A, B = load(D / "pass_a" / fname), load(D / "pass_b" / fname)
    cols = list(B[0].keys())
    diffs = []
    for i, rb in enumerate(B):
        if align == "index":
            ra = A[i] if i < len(A) else None
        else:  # nearest alpha among unused A rows, within 0.015 deg
            cand = [r for r in A if same(r.get("alpha_deg"), rb["alpha_deg"])]
            ra = cand[0] if cand else None
            if ra is not None:
                A.remove(ra)
        for c in cols:
            va = None if ra is None else ra.get(c)
            if not same(va, rb[c]):
                diffs.append({"file": fname, "row": i, "column": c, "pass_a": va, "pass_b": rb[c]})
    return diffs, len(B) * len(cols)


def main():
    res = resolutions()
    files = [("polar_re200k.csv", "alpha"), ("polar_re100k.csv", "index"), ("cp_upper_re200k.csv", "index"),
             ("oilflow.csv", "index")]
    report, unresolved, cells = [], [], 0
    for f, al in files:
        d, n = compare(f, al)
        cells += n
        canon = load(D / "canonical" / f) if (D / "canonical" / f).exists() else None
        for x in d:
            key = (x["file"], x["row"], x["column"])
            x["resolved"] = res.get(key)
            ok = x["resolved"] is not None and canon is not None and same(canon[x["row"]][x["column"]], x["resolved"])
            if not ok:
                unresolved.append(x)
        report += d
        print(f"{f}: {len(d)} raw disagreements out of {n} cells")
    out = {"cells_compared": cells, "raw_disagreements": len(report), "unresolved": len(unresolved),
           "disagreements": report}
    (D / "diff_report.json").write_text(json.dumps(out, indent=1))
    print(f"TOTAL cells {cells}, raw disagreements {len(report)}, unresolved after page-image resolution {len(unresolved)}")
    for x in unresolved[:80]:
        print("  UNRESOLVED", x["file"], "row", x["row"], x["column"], "A=", x["pass_a"], "B=", x["pass_b"])
    sys.exit(0 if not unresolved else 1)


if __name__ == "__main__":
    main()
