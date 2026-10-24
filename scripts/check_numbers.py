"""V19 traceability check. Prints 'UNTRACEABLE <n>'; n must be 0.

Counted as untraceable:
1. an outward document (README.md, RESUME_CONTRACT.md, report/FlightLoop_Report.md) that differs from a fresh render of
   its template against the current results/numbers.json (a hand edit, or numbers that moved after the last build);
2. a placeholder whose key is missing from numbers.json;
3. a result-like number typed into a template outside placeholders: a number followed by %, pp, 'deg' or
   'core-h' (except the '95% CI' confidence level), on a line that does not state a pre-registered bar or a plan setting (marked with 'bar', 'Bar' or 'setting');
4. site/data/bundle.js whose numbers differ from numbers.json.
"""
import json, pathlib, re, sys, importlib.util
ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("bd", ROOT / "scripts/build_docs.py"); BD = importlib.util.module_from_spec(spec); spec.loader.exec_module(BD)

bad = []
numbers = json.loads((ROOT / "results/numbers.json").read_text())
try:
    fresh = BD.build(write=False)
except KeyError as e:
    bad.append(f"missing key: {e}"); fresh = {}
for dst, txt in fresh.items():
    p = ROOT / dst
    if not p.exists() or p.read_text() != txt:
        bad.append(f"{dst} differs from a fresh render of its template")
LIT = re.compile(r"(?<![\w.{:])[-+]?\d+(?:\.\d+)?\s?(?:%|pp\b|deg\b|core-h\b)")
for src in BD.TEMPLATES:
    for i, ln in enumerate((ROOT / src).read_text().splitlines(), 1):
        stripped = BD.PH.sub("", ln)
        if re.search(r"\b(bar|Bar|bars|setting)\b", ln):
            continue
        for m in LIT.finditer(stripped):
            if stripped[m.end():m.end() + 3] == " CI":
                continue  # confidence level of an interval (a method setting, not a result)
            bad.append(f"{src}:{i} literal '{m.group(0)}'")
b = ROOT / "site/data/bundle.js"
if b.exists():
    t = b.read_text(); FL = json.loads(t[t.index("window.FL = ") + 12:].rstrip().rstrip(";"))
    if {k: v["value"] for k, v in FL["numbers"].items()} != {k: v["value"] for k, v in numbers.items()}:
        bad.append("site/data/bundle.js numbers differ from results/numbers.json (rebuild the site)")
else:
    bad.append("site/data/bundle.js missing")
for x in bad:
    print(" -", x)
print(f"UNTRACEABLE {len(bad)}")
sys.exit(1 if bad else 0)
