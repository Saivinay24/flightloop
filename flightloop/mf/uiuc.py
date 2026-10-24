"""UIUC LSAT volumes 1 and 3 (public experimental): drag polars (DRAGxx.TXT) and as-tested coordinates (COORDxx.TXT).

Only clean configurations are used: the polar file's Comment must start with "clean" and mention no flap, trip,
bay or vortex generator. Coordinates: the as-tested .dap whose title matches the polar's airfoil name.
"""
import re, pathlib
import numpy as np
import pandas as pd
from ..config import ROOT

D = ROOT / "data/uiuc_lsat"
BAD = re.compile(r"flap|trip|bay|vortex|gurney|tape|covering", re.I)


def _blocks(path):
    txt = pathlib.Path(path).read_text(errors="ignore")
    parts = re.split(r"^:{6,}\s*\n(.+?)\n:{6,}\s*\n", txt, flags=re.M)
    return [(parts[i].strip(), parts[i + 1]) for i in range(1, len(parts) - 1, 2)]


def norm(name):
    return re.sub(r"[^a-z0-9()]", "", name.lower().replace("gottingen", "goe"))


def polars():
    rows = []
    for f in ("DRAG01.TXT", "DRAG03.TXT"):
        vol = f[4:6]
        for fname, body in _blocks(D / f):
            m = re.search(r"Airfoil:\s*(.+)", body); c = re.search(r"Comment:\s*(.*)", body)
            if not m:
                continue
            name, comment = m.group(1).strip(), (c.group(1).strip() if c else "")
            if not comment.lower().startswith("clean") or BAD.search(comment):
                continue
            lines = body.splitlines()
            i = 0
            while i < len(lines):
                if lines[i].startswith("Average Reynolds"):
                    re_ = float(lines[i + 1].split()[0]); n = int(lines[i + 3].split()[0])
                    j = i + 4
                    while j < len(lines) and not re.match(r"^\s*-?\d", lines[j]):
                        j += 1
                    for k in range(n):
                        t = lines[j + k].split()
                        rows.append({"file": fname, "vol": vol, "airfoil": name, "re": re_, "alpha": float(t[0]), "cl": float(t[1]), "cd": float(t[2])})
                    i = j + n
                else:
                    i += 1
    return pd.DataFrame(rows)


def coords():
    out = {}
    for f in ("COORD01.TXT", "COORD03.TXT"):
        for fname, body in _blocks(D / f):
            if not fname.lower().endswith(".dap"):
                continue
            lines = [l for l in body.splitlines() if l.strip()]
            title = lines[0].strip()
            pts = []
            for l in lines[1:]:
                t = l.split()
                try:
                    pts.append((float(t[0]), float(t[1])))
                except (ValueError, IndexError):
                    continue
            P = np.array(pts)
            if len(P) < 20:
                continue
            if P[0, 0] > 1.5:  # Lednicer header (counts)
                nu, nl = int(P[0, 0]), int(P[0, 1]); up, lo = P[1:1 + nu], P[1 + nu:1 + nu + nl]
                P = np.r_[up[::-1], lo[1:]]
            out[norm(title)] = {"file": fname, "title": title, "xy": P}
    return out


FAMILY = [("SD", r"^SD\d"), ("SA", r"^SA\d"), ("SG", r"^SG\d"), ("S", r"^S\d"), ("E", r"^E\d"), ("FX", r"^FX"),
          ("NACA", r"^NACA"), ("RG", r"^RG\d"), ("MH", r"^MH\d"), ("GM", r"^GM\d")]


def family(name):
    for fam, pat in FAMILY:
        if re.match(pat, name.upper()):
            return fam
    return "other:" + re.sub(r"\s*\(.*\)", "", name).strip().upper()


def base_name(name):
    return re.sub(r"\s*\(.*\)", "", name).strip().upper()
