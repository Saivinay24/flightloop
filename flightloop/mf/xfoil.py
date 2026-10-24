"""Headless XFOIL 6.99 batch (binary from the StratoGlide repo, used read-only). One polar per (airfoil, Re, Ncrit).

Panels: PPAR N 200. Sweep alpha 0 -> 14 then a fresh BL init 0 -> -4, step 0.5 (= ASEQ -4 14 0.5 split in two
branches so each starts from an attached solution). If XFOIL hangs it is killed after `timeout` s (15 s; a full 37-point sweep takes about 2.5 s); converged points
already written are kept and the branch restarts one step past the hang (the hang alpha is logged as
non-converged). Cached as JSON under results/mf/xfoil/ keyed by coordinate SHA, Re, Ncrit.
"""
import hashlib, json, pathlib, subprocess, tempfile
import numpy as np
from ..config import ROOT, cfg

XFOIL = pathlib.Path(cfg()["stratoglide_root"]).expanduser() / "tools/bin/xfoil"
CACHE = ROOT / "results/mf/xfoil"


def _parse(path):
    rows, go = [], False
    if not path.exists():
        return rows
    for ln in path.read_text().splitlines():
        if ln.strip().startswith("------"):
            go = True; continue
        if go:
            t = ln.split()
            if len(t) >= 7:
                rows.append([float(x) for x in t[:7]])
    return rows


def _branch(wd, re_, ncrit, a0, a1, step, timeout):
    rows, bad = [], []
    sgn = 1.0 if a1 >= a0 else -1.0
    start, tries = a0, 0
    while sgn * (a1 - start) >= -1e-9 and tries < 30:
        tries += 1
        pol = wd / f"p{'u' if sgn > 0 else 'l'}{tries}.txt"
        if pol.exists():
            pol.unlink()
        cmd = ["PLOP", "G F", "", "LOAD foil.dat", "PPAR", "N 200", "", "", "OPER", f"VISC {re_:.0f}", "VPAR", f"N {ncrit}", "",
               "ITER 150", "PACC", pol.name, "", f"ASEQ {start:.2f} {a1:.2f} {sgn * step:.2f}", "PACC", "", "QUIT"]
        to = False
        try:
            subprocess.run([str(XFOIL)], input="\n".join(cmd) + "\n", text=True, capture_output=True, timeout=timeout, cwd=wd)
        except subprocess.TimeoutExpired:
            to = True
        got = _parse(pol); rows += got
        if not to:
            break
        last = got[-1][0] if got else start - sgn * step
        bad.append(round(last + sgn * step, 2)); start = last + 2 * sgn * step
    return rows, bad


def polar(xy, name, re_, ncrit, amin=-4.0, amax=14.0, step=0.5, timeout=15.0):
    xy = np.asarray(xy)
    txt = name + "\n" + "\n".join(f"{x:.6f} {y:.6f}" for x, y in xy) + "\n"
    key = hashlib.sha256(json.dumps([txt, re_, ncrit, amin, amax, step, "v2"]).encode()).hexdigest()[:20]
    CACHE.mkdir(parents=True, exist_ok=True)
    cp = CACHE / f"{key}.json"
    if cp.exists():
        return json.loads(cp.read_text())
    with tempfile.TemporaryDirectory(dir=ROOT / "results/mf") as td:
        wd = pathlib.Path(td)
        (wd / "foil.dat").write_text(txt)
        up, b1 = _branch(wd, re_, ncrit, 0.0, amax, step, timeout)
        lo, b2 = _branch(wd, re_, ncrit, 0.0, amin, step, timeout)
    rows = sorted({round(r[0], 3): r for r in up + lo}.values())
    out = {"name": name, "re": re_, "ncrit": ncrit, "cols": ["alpha", "cl", "cd", "cdp", "cm", "xtr_top", "xtr_bot"],
           "rows": rows, "nonconverged_alpha": b1 + b2, "n_requested": int(round((amax - amin) / step)) + 1}
    cp.write_text(json.dumps(out))
    return out
