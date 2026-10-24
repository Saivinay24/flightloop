"""Resumable SU2 queue. Usage: cfd_queue.py <queue_name>   (cases defined in QUEUES below)

Manifest results/queues/<name>.json maps each case key to its run_id (created once, reused on resume).
At most compute.max_parallel_solver_procs/mpi_np concurrent jobs; a job is (re)launched from its restart file
until postproc says converged/unsteady (iteration cap). Skips finished runs. Checks free memory before launching.
"""
import json, pathlib, subprocess, sys, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from flightloop.config import ROOT, cfg, run_id as make_run_id
from flightloop.cfd import su2_runner as R, postproc as PP

C = cfg()
TU = C["validation"]["e387"]["headline_tu_percent"]
CAP = {"SA": 12000, "SST": 15000, "SA_BC": 15000, "SST_LM": 30000}


def case(af, re, a, model, level, tu=TU, tag_tu=False):
    return {"airfoil": af, "re": re, "alpha": a, "model": model, "level": level, "tu_pct": tu, "iters": CAP[model],
            "tag_tu": tag_tu}


def key(c):
    return f"{c['airfoil']}|{c['re']}|{c['alpha']}|{c['model']}|{c['level']}|{c['tu_pct']}"


QUEUES = {
    "cutline": [case("e387", 200000, 4, m, "medium") for m in ("SST_LM", "SA")]
    + [case("e387", 200000, a, m, "medium") for a in (0, 8) for m in ("SST_LM", "SA")]
    + [case("e387", 200000, 4, m, lv) for m in ("SST_LM", "SA") for lv in ("coarse", "fine")],
    "fullA": [case("e387", 200000, a, m, "medium") for m in ("SST", "SA_BC") for a in (0, 4, 8)]
    + [case("e387", 100000, a, m, "medium") for m in ("SST_LM", "SA") for a in (0, 4, 8)]
    + [case("e387", 200000, a, "SST_LM", "medium", tu, True) for a in (0, 4) for tu in (0.05, 0.1, 0.2)]
    + [case("e387", 200000, a, m, "medium") for m in ("SST_LM", "SA") for a in (-2, 2, 6, 10)]
    + [case("e387", 100000, a, m, "medium") for m in ("SST", "SA_BC") for a in (0, 4, 8)],
    # compute-bounded (COMPUTE.md): Module B needs Re 60k and 100k at the cruise Cl band; alpha 2 and 4 only.
    "strato": [case("sd7037", re, a, m, "medium") for re in (100000, 60000) for a in (4, 2) for m in ("SA", "SST_LM")],
    # Module A extensions in priority order (turbulence-model comparison, Re 100k, Tu sensitivity), alpha 4 first
    # plan-order cut (GAPS.md): the two remaining fully turbulent models at the headline point, run instead of SD7037 Re 60k
    "fullA_min": [case("e387", 200000, 4, m, "medium") for m in ("SST", "SA_BC")],
    "fullA_core": [case("e387", 200000, 4, m, "medium") for m in ("SST", "SA_BC")]
    + [case("e387", 100000, 4, m, "medium") for m in ("SA", "SST_LM")]
    + [case("e387", 200000, 4, "SST_LM", "medium", tu, True) for tu in (0.1, 0.2)]
    + [case("e387", 200000, a, m, "medium") for m in ("SST", "SA_BC") for a in (0, 8)],
}


def mem_free_pct():
    out = subprocess.run(["memory_pressure"], capture_output=True, text=True).stdout
    for line in out.splitlines():
        if "free percentage" in line:
            return float(line.split(":")[1].strip().rstrip("%"))
    return 100.0


def main(name):
    qdir = ROOT / "results/queues"; qdir.mkdir(parents=True, exist_ok=True)
    mpath = qdir / f"{name}.json"
    man = json.loads(mpath.read_text()) if mpath.exists() else {}
    cases = QUEUES[name]
    for c in cases:
        if key(c) not in man:
            nm = f"{c['airfoil']}-re{int(c['re']/1000)}k-a{c['alpha']:g}-{c['model']}-{c['level']}" + (f"-tu{c['tu_pct']:g}" if c["tag_tu"] else "")
            man[key(c)] = make_run_id("cfd", nm)
            time.sleep(1.1)
    mpath.write_text(json.dumps(man, indent=1))
    np_ = C["compute"]["mpi_np"]
    slots = C["compute"]["max_parallel_solver_procs"] // np_ if C["compute"]["max_parallel_solver_procs"] >= 2 * np_ else 2
    slots = 2
    running = {}
    log = open(qdir / f"{name}.log", "a")
    def say(*a):
        print(time.strftime("%H:%M:%S"), *a, file=log, flush=True)
    pending = list(cases)
    while pending or running:
        for k, (c, d) in list(running.items()):
            if (d / "exit_code").exists():
                r = PP.process(d)
                st = r["status"]["state"]
                say("finished chunk", k, st, r["status"].get("iters"), r["status"].get("cl"), r["status"].get("cd"))
                del running[k]
                if st in ("incomplete",) and r["status"]["iters"] < c["iters"]:
                    pending.insert(0, c)
                elif st == "failed":
                    say("FAILED", k)
        while pending and len(running) < slots:
            if mem_free_pct() < 25:
                say("memory low, waiting"); break
            c = pending.pop(0)
            k = key(c)
            d = ROOT / "runs" / man[k]
            if (d / "result.json").exists():
                st = json.loads((d / "result.json").read_text())["status"]
                if st["state"] in ("converged", "unsteady"):
                    continue
                done = st.get("iters", 0)
            else:
                done = 0
            remaining = c["iters"] - done
            if remaining <= 0:
                continue
            cc = dict(c, iters=remaining)
            d = R.prepare(cc, rid=man[k])
            cj = json.loads((d / "case.json").read_text()); cj["iters"] = c["iters"]; (d / "case.json").write_text(json.dumps(cj, indent=1))
            ec = d / "exit_code"
            if ec.exists():
                ec.unlink()
            pid = R.launch(d, np_=np_, background=True)
            say("launched", k, man[k], "pid", pid, "remaining", remaining)
            running[k] = (c, d)
        time.sleep(30)
    say("queue done")


if __name__ == "__main__":
    main(sys.argv[1])
