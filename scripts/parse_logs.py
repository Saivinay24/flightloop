"""P8: parse public real-format logs to canonical Parquet and report schema validity (verify.py V14).
PX4: data/px4/*.ulg (public Flight Review logs). ArduPilot: data/alfa/dataflash/**/*.bin (CMU ALFA, real ArduPlane)."""
import glob, json, sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from flightloop.logs import schema, px4, ardupilot
from flightloop.config import put_number

out = ROOT / "results/logs"; out.mkdir(parents=True, exist_ok=True)
pq = ROOT / "runs/logs_real"; pq.mkdir(parents=True, exist_ok=True)
files = []
for f in sorted(glob.glob(str(ROOT / "data/px4/*.ulg"))):
    df, m = px4.parse(f); iss = schema.validate(df)
    schema.write(df, pq / (pathlib.Path(f).stem + ".parquet"), {"source": "px4_public", **m, "file": pathlib.Path(f).name})
    files.append({"file": str(pathlib.Path(f).relative_to(ROOT)), "format": "ulog", "rows": len(df), "issues": iss, "schema_ok": not iss,
                  "monotone_t": bool(np.all(np.diff(df.t) > 0)), "tas_mean": float(df.tas.mean()) if df.tas.notna().any() else None})
bins = sorted(glob.glob(str(ROOT / "data/alfa/dataflash/**/*.bin"), recursive=True))[:3]
for f in bins:
    df, m = ardupilot.parse(f); iss = schema.validate(df)
    schema.write(df, pq / (pathlib.Path(f).stem.replace(" ", "_") + ".parquet"), {"source": "ardupilot_public", **m, "file": pathlib.Path(f).name})
    files.append({"file": str(pathlib.Path(f).relative_to(ROOT)), "format": "bin", "rows": len(df), "issues": iss, "schema_ok": not iss,
                  "monotone_t": bool(np.all(np.diff(df.t) > 0)), "tas_mean": float(df.tas.mean())})
(out / "parse_report.json").write_text(json.dumps({"files": files}, indent=1))
put_number("logs.px4_ok", sum(f["format"] == "ulog" and f["schema_ok"] for f in files), "-", "scripts/parse_logs.py", "public PX4 Flight Review logs", "n/a")
put_number("logs.ardupilot_ok", sum(f["format"] == "bin" and f["schema_ok"] for f in files), "-", "scripts/parse_logs.py", "public ArduPlane logs (ALFA)", "n/a")
for f in files:
    print(f["format"], f["file"][-40:], f["rows"], f["issues"])
