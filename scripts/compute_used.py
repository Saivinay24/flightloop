"""Actual compute used, from the run folders (SU2 per-iteration wall time summed from history.csv x MPI ranks) and the
recorded wall times of the XFOIL batch, sys-ID Monte Carlo and blind loop. Writes compute.* numbers and
results/compute_used.json; COMPUTE.md quotes them through the report build."""
import csv, glob, json, pathlib, re, sys
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from flightloop.config import put_number

SCRIPT = "scripts/compute_used.py"


def su2_wall(d):
    h = d / "history.csv"
    if not h.exists():
        return 0.0, 0
    with open(h) as fh:
        r = csv.reader(fh); head = [c.strip().strip('"') for c in next(r)]
        i = head.index("Time(sec)")
        tot, n = 0.0, 0
        for row in r:
            try:
                tot += float(row[i]); n += 1
            except (ValueError, IndexError):
                pass
    return tot, n


def euler_wall(d):
    """Euler histories have no time column: wall = exit_code mtime - mesh.su2 symlink creation (file-time estimate)."""
    e, m = d / "exit_code", d / "mesh.su2"
    if not e.exists() or not m.is_symlink():
        return 0.0
    return max(e.stat().st_mtime - m.lstat().st_mtime, 0.0)


out = {"groups": {}}
for g, pat, nr in (("cfd_incompressible", "runs/cfd-*", 2), ("cfd_euler", "runs/comp-*", 1)):
    # incompressible runs: mpirun -np 2 (flightloop/cfd/su2_runner.py); Euler runs: serial SU2_CFD (su2_euler.run)
    w = ch = 0.0; n = 0
    for d in sorted(glob.glob(str(ROOT / pat))):
        d = pathlib.Path(d); t = su2_wall(d)[0] if nr == 2 else euler_wall(d)
        w += t; ch += t * nr / 3600; n += 1
    out["groups"][g] = {"runs": n, "wall_h_single_job": w / 3600, "core_h": ch}
mf = json.loads((ROOT / "results/mf/summary.json").read_text())["meta"]
out["groups"]["xfoil_batch"] = {"wall_h": mf["xfoil_wall_s"] / 3600, "core_h": 3 * mf["xfoil_wall_s"] / 3600, "note": "Pool(3)"}
lp = json.loads((ROOT / "results/loop/blind_summary.json").read_text())
out["groups"]["blind_loop"] = {"wall_h": lp["wall_s"] / 3600, "core_h": lp["wall_s"] / 3600}
tot = sum(v["core_h"] for v in out["groups"].values())
out["total_core_h"] = tot
out["cloud_gpu_h"] = 0
(ROOT / "results/compute_used.json").write_text(json.dumps(out, indent=1))
src = "measured from run folders"
put_number("compute.total_core_h", tot, "core-h", SCRIPT, src, "multiple", "SU2 + XFOIL batch + blind loop; excludes sys-ID Monte Carlo and post-processing (minutes)")
put_number("compute.su2_core_h", out["groups"]["cfd_incompressible"]["core_h"] + out["groups"]["cfd_euler"]["core_h"], "core-h", SCRIPT, src, "multiple")
put_number("compute.su2_runs", out["groups"]["cfd_incompressible"]["runs"] + out["groups"]["cfd_euler"]["runs"], "-", SCRIPT, src, "multiple", "including superseded runs")
put_number("compute.cloud_gpu_h", 0, "h", SCRIPT, src, "n/a", "no cloud VM used")
print(json.dumps(out, indent=1))
