"""Wait for already-launched SU2 run folders to finish and post-process them exactly as cfd_queue.py does
(used when a queue manager is stopped while its last single-chunk jobs are still running)."""
import sys, time, pathlib, importlib.util
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location("q", ROOT / "scripts/cfd_queue.py"); Q = importlib.util.module_from_spec(spec); spec.loader.exec_module(Q)
dirs = [ROOT / "runs" / a for a in sys.argv[1:]]
left = list(dirs)
while left:
    for d in list(left):
        if (d / "exit_code").exists():
            r = Q.PP.process(d)
            print(time.strftime("%H:%M:%S"), "finished", d.name, r["status"]["state"], r["status"].get("iters"), flush=True)
            left.remove(d)
    time.sleep(30)
