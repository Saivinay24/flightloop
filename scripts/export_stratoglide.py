"""Copy StratoGlide's outputs into inputs/stratoglide/ (read-only use of the sister repo) with SHA-256 provenance."""
import json, shutil, sys, pathlib, hashlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from flightloop.config import cfg
SG = pathlib.Path(cfg()["stratoglide_root"]).expanduser()
sys.path.insert(0, str(SG / "src"))
from stratoglide.release import aero_tables  # noqa: E402

out = ROOT / "inputs/stratoglide"; out.mkdir(parents=True, exist_ok=True)
files = {"airfoil.dat": "outputs/geometry/E387.dat", "aircraft.avl": "outputs/geometry/stratoglide.avl",
         "aircraft.mass": "outputs/geometry/stratoglide.mass", "numbers.json": "results/numbers.json",
         "predictions.json": "outputs/test_cards/predictions.json", "tail_airfoil.dat": "outputs/geometry/HT12.dat"}
prov = {}
for dst, src in files.items():
    shutil.copy2(SG / src, out / dst)
    prov[dst] = {"source": f"~/Desktop/stratoglide/{src}", "sha256_16": hashlib.sha256((SG / src).read_bytes()).hexdigest()[:16]}
p = aero_tables.nominal()
p = {k: (float(v) if not isinstance(v, str) else v) for k, v in p.items()}
import subprocess
tag = subprocess.run(["git", "-C", str(SG), "describe", "--tags", "--always"], capture_output=True, text=True).stdout.strip()
(out / "derivs.json").write_text(json.dumps({"source": "stratoglide.release.aero_tables.nominal() (AVL 20 km case derivatives, P3 sizing polar)",
                                              "stratoglide_git": tag, "p": p}, indent=1))
prov["derivs.json"] = {"source": "aero_tables.nominal()", "stratoglide_git": tag}
(out / "PROVENANCE.json").write_text(json.dumps(prov, indent=1))
print(json.dumps(prov, indent=1))
