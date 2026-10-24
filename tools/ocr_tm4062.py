"""Pass A: independent machine transcription of the TM-4062 tables used by FlightLoop.

Polars (Table B1): the PDF's embedded OCR text layer (Pdf.Capture engine, independent of pass B).
Cp listings (Appendix D, rotated pages): tesseract 5 at 300 dpi on the rotated page.
Both are parsed digit-only with the known fixed-decimal format of the NASA listing (alpha x.xx, cl x.xxx,
cd .xxxx, cm .xxxx, Cp x.xxxx, x/c .xxx) so OCR'd decimal points do not matter.
Writes data/tm4062/pass_a/*.csv; scripts/diff_transcriptions.py compares with pass_b.
"""
import pathlib, re, subprocess

ROOT = pathlib.Path(__file__).resolve().parents[1]
PDF = ROOT / "data/tm4062/tm4062.pdf"
WORK = ROOT / "data/tm4062/work"; WORK.mkdir(parents=True, exist_ok=True)
OUT = ROOT / "data/tm4062/pass_a"; OUT.mkdir(exist_ok=True)
MAP = str.maketrans({"O": "0", "o": "0", "D": "0", "Q": "0", "I": "1", "l": "1", "i": "1", "q": "9", "g": "9",
                     "S": "5", "s": "5", "B": "8", "Z": "2", "z": "2", "?": "7", "~": "-", "=": "-", "—": "-"})


def fixed(field, dec):
    f = field.translate(MAP)
    if "*" in f or "#" in f or "W" in f:
        return None
    neg = "-" in f
    d = re.sub(r"\D", "", f)
    if not d:
        raise ValueError(field)
    v = int(d) / 10 ** dec
    return -v if neg else v


def textlayer(page):
    return subprocess.run(["pdftotext", "-f", str(page), "-l", str(page), "-layout", str(PDF), "-"],
                          capture_output=True, text=True).stdout


def polar_layout(page, fname):
    """Layout text layer, one table row per line (works when the scan is not skewed)."""
    txt = subprocess.run(["pdftotext", "-f", str(page), "-l", str(page), "-layout", str(PDF), "-"],
                         capture_output=True, text=True).stdout.split("HYSTER")[0]
    out = []
    for line in txt.splitlines():
        fields = [f for f in re.split(r"\s{3,}", line.strip()) if f]
        if len(fields) != 4:
            continue
        vals = []
        for fld, d in zip(fields, [2, 3, 4, 4]):
            try:
                v = fixed(fld, d)
                vals.append("" if v is None else f"{v:.{d}f}")
            except ValueError:
                vals.append("?")
        if vals[0] != "?":
            out.append(vals)
    with open(OUT / fname, "w") as f:
        f.write("alpha_deg,cl,cd,cm\n")
        for v in out:
            f.write(",".join(v) + "\n")
    print(fname, len(out), "rows (layout mode)")


def polar(page, fname):
    """Raw text layer lists each table column top to bottom; split on the column headers."""
    raw = subprocess.run(["pdftotext", "-f", str(page), "-l", str(page), "-raw", str(PDF), "-"],
                         capture_output=True, text=True).stdout.split("HYSTER")[0]
    lines = [l.strip() for l in raw.splitlines() if l.strip()]
    cols, cur = [], None
    for l in lines:
        u = l.upper().replace(" ", "")
        if u.startswith("ALPHA") or u in ("CL", "CT", "C0", "CD", "CO", "CM", "CH"):
            cur = []; cols.append(cur); continue
        if cur is None or re.search(r"PSI|PTINF|RUN|TABLE|^M[-=]|^[•R]*.?[I1]?[0O]{2}[,.][0O]{3}", u):
            continue
        if len(l) <= 9:
            cur.append(l)
    cols = [c for c in cols if c][:4]
    decs = [2, 3, 4, 4]
    n = max(len(c) for c in cols)
    with open(OUT / fname, "w") as f:
        f.write("alpha_deg,cl,cd,cm\n")
        for i in range(n):
            vals = []
            for c, d in zip(cols, decs):
                try:
                    v = fixed(c[i], d) if i < len(c) else "?"
                except ValueError:
                    v = "?"
                vals.append("" if v is None else (v if v == "?" else f"{v:.{d}f}"))
            f.write(",".join(vals) + "\n")
    print(fname, [len(c) for c in cols], "entries per column")


def tess_rotated(page):
    stem = WORK / f"ocr{page}"
    subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page), "-r", "300", "-gray", "-singlefile", "-png", str(PDF), str(stem)], check=True)
    subprocess.run(["sips", "-r", "90", f"{stem}.png"], capture_output=True)
    return subprocess.run(["tesseract", f"{stem}.png", "-", "--psm", "6"], capture_output=True, text=True).stdout


def cp(pages_cols, fname):
    cols, xs = {}, None
    for page, (col_idx, label) in pages_cols:
        txt = tess_rotated(page)
        (WORK / f"tess{page}.txt").write_text(txt)
        x_here, v_here = [], []
        for line in txt.splitlines():
            toks = line.split()
            if len(toks) < 9:
                continue
            try:
                vals = [fixed(t, 4) for t in toks[1:9]]
                x = fixed(toks[0], 3)
            except ValueError:
                continue
            x_here.append(x); v_here.append(vals[col_idx])
        xs = x_here
        cols[label] = v_here
    n = min(len(v) for v in cols.values())
    with open(OUT / fname, "w") as f:
        f.write("x_c," + ",".join(f"cp_a{k}" for k in cols) + "\n")
        for i in range(n):
            f.write(f"{xs[i]:.3f}," + ",".join(f"{cols[k][i]:.4f}" for k in cols) + "\n")
    print(fname, n, "rows")


def oilflow():
    # Table III is typeset (not a listing); the text layer is column-ordered. Keep raw for the diff script.
    txt = subprocess.run(["pdftotext", "-f", "15", "-l", "15", "-raw", str(PDF), "-"], capture_output=True, text=True).stdout
    (OUT / "oilflow_textlayer.txt").write_text(txt)
    print("oilflow text layer saved")


if __name__ == "__main__":
    polar_layout(96, "polar_re200k.csv")
    polar(97, "polar_re100k.csv")
    # Cp: tesseract on the rotated dot-matrix listing was unusable (see HONEST_LOG); pass A for Cp is a second
    # manual read of 300 dpi column crops made by tools/crop_columns.py, kept in pass_a/cp_upper_re200k.csv.
    oilflow()
