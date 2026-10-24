"""P4: wedge (M 2, 3) and diamond (M 2, 3; alpha 0, 2, 4) on 3 nested grids; serial SU2, 2 cases in parallel."""
import json, sys, pathlib
from concurrent.futures import ThreadPoolExecutor
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from flightloop.config import cfg, run_id as make_rid
from flightloop.compressible import su2_euler as E

C = cfg()["compressible"]
man_p = ROOT / "results/queues/compressible.json"; man_p.parent.mkdir(parents=True, exist_ok=True)
man = json.loads(man_p.read_text()) if man_p.exists() else {}
jobs = []
for M in C["wedge"]["mach"]:
    for lv in (0, 1, 2):
        jobs.append(("wedge", lv, M, 0.0))
for M in C["diamond"]["mach"]:
    for a in C["diamond"]["alpha_deg"]:
        for lv in (0, 1, 2):
            jobs.append(("diamond", lv, M, float(a)))
for j in jobs:
    k = "|".join(map(str, j))
    if k not in man:
        man[k] = make_rid("comp", f"{j[0]}-M{j[2]:g}-a{j[3]:g}-L{j[1]}")
        import time; time.sleep(1.05)
man_p.write_text(json.dumps(man, indent=1))
with ThreadPoolExecutor(2) as ex:
    list(ex.map(lambda j: E.run(j[0], j[1], j[2], j[3], man["|".join(map(str, j))]), jobs))
print("done", len(jobs))
