"""Render outward documents from templates/ with numbers pulled from results/numbers.json, then the report PDF.

Placeholders in templates:
    {{key}}            value of numbers.json[key], default format (3 significant figures for floats)
    {{key:+.1f}}       Python format spec applied to the value (strings are inserted unchanged)
    {{file:PATH}}      the full text of a generated file (e.g. results/strato/DISAGREEMENT.md table)
    {{file:PATH#Section}}  only the '## Section' block of that file
    {{verify}}         the verify table (id, class, PASS/FAIL, message) from results/verify_report.json
No result number is typed into a template; scripts/check_numbers.py re-renders and fails on any drift and on bare
result-like numbers outside placeholders.

    .venv/bin/python scripts/build_docs.py
"""
import json, pathlib, re, subprocess, sys
ROOT = pathlib.Path(__file__).resolve().parents[1]
TEMPLATES = {"templates/README.md.tmpl": "README.md", "templates/RESUME_CONTRACT.md.tmpl": "RESUME_CONTRACT.md",
             "templates/FlightLoop_Report.md.tmpl": "report/FlightLoop_Report.md"}
PH = re.compile(r"\{\{([^{}]+)\}\}")
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"


def fmt_value(v, spec=None):
    if isinstance(v, bool) or isinstance(v, str) or v is None:
        return str(v)
    if spec:
        return format(v, spec)
    if isinstance(v, int):
        return str(v)
    return f"{v:.3g}"


def file_block(arg):
    path, _, sec = arg.partition("#")
    t = (ROOT / path).read_text()
    if not sec:
        return t.strip()
    m = re.search(rf"^## {re.escape(sec)}\s*$(.*?)(?=^## |\Z)", t, re.S | re.M)
    if not m:
        raise KeyError(f"section {sec} not in {path}")
    return m.group(1).strip()


def verify_table():
    v = json.loads((ROOT / "results/verify_report.json").read_text())
    L = ["| Check | Class | Result | Measured |", "|---|---|---|---|"]
    for r in v["rows"]:
        L.append(f"| {r['id']} | {r['class']} | {r['tag']} | {r['msg'].replace('|', '/')} |")
    L.append("")
    L.append(f"INTEGRITY {v['integrity']}/{v['n_integrity']}, RESULTS {v['results']}/{v['n_results']}.")
    return "\n".join(L)


def render(text, numbers):
    def sub(m):
        tok = m.group(1).strip()
        if tok == "verify":
            return verify_table()
        if tok.startswith("file:"):
            return file_block(tok[5:])
        key, _, spec = tok.partition(":")
        if key not in numbers:
            raise KeyError(f"numbers.json has no key {key}")
        return fmt_value(numbers[key]["value"], spec or None)
    return PH.sub(sub, text)


def build(write=True):
    numbers = json.loads((ROOT / "results/numbers.json").read_text())
    out = {}
    for src, dst in TEMPLATES.items():
        out[dst] = render((ROOT / src).read_text(), numbers)
        if write:
            (ROOT / dst).parent.mkdir(parents=True, exist_ok=True)
            (ROOT / dst).write_text(out[dst])
    return out


def pdf():
    import markdown
    md = (ROOT / "report/FlightLoop_Report.md").read_text()
    body = markdown.markdown(md, extensions=["tables", "fenced_code"])
    css = """body{font:10.5pt/1.45 -apple-system,Helvetica,Arial,sans-serif;color:#0b1324;margin:0}
h1{font-size:20pt;margin:0 0 4pt;color:#0b1324;border-bottom:3pt solid #FED176;padding-bottom:4pt}
h2{font-size:13.5pt;margin:16pt 0 4pt;color:#111c33;border-bottom:1pt solid #c9d3e6;padding-bottom:2pt}
h3{font-size:11pt;margin:10pt 0 3pt}
table{border-collapse:collapse;width:100%;margin:6pt 0;font-size:8.6pt;page-break-inside:auto}
th,td{border:0.6pt solid #c9d3e6;padding:2.5pt 4pt;text-align:left;vertical-align:top}
th{background:#eef2f9}
code{font:8.6pt Menlo,monospace;background:#f3f5fa;padding:0 2pt}
pre{background:#f3f5fa;padding:6pt;font-size:8.4pt;white-space:pre-wrap}
@page{size:A4;margin:16mm 14mm}"""
    html = f"<!doctype html><html><head><meta charset='utf-8'><title>FlightLoop Report</title><style>{css}</style></head><body>{body}</body></html>"
    hp = ROOT / "report/FlightLoop_Report.html"
    hp.write_text(html)
    subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--no-pdf-header-footer", f"--print-to-pdf={ROOT / 'report/FlightLoop_Report.pdf'}",
                    hp.as_uri()], capture_output=True, timeout=120)
    return (ROOT / "report/FlightLoop_Report.pdf").exists()


if __name__ == "__main__":
    build()
    print("rendered", ", ".join(TEMPLATES.values()), "| pdf", pdf())
