"""Config loader, run-id helper and numbers.json writer (single source of truth with provenance)."""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import pathlib
import subprocess

import yaml

ROOT = pathlib.Path(__file__).resolve().parents[1]


def cfg() -> dict:
    with open(ROOT / "config/flightloop.yaml") as f:
        return yaml.safe_load(f)


def git_sha() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip() or "nogit"
    except Exception:
        return "nogit"


def run_id(module: str, case: str) -> str:
    ts = _dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    return f"{module}-{case}-{ts}-{git_sha()}"


def sha_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


NUMBERS = ROOT / "results/numbers.json"


def put_number(key: str, value, unit: str, script: str, data_source: str, run_id: str = "n/a", note: str = ""):
    """Write one number with provenance. Called only by pipeline scripts, never by hand."""
    NUMBERS.parent.mkdir(parents=True, exist_ok=True)
    d = json.loads(NUMBERS.read_text()) if NUMBERS.exists() else {}
    if isinstance(value, float):
        value = float(f"{value:.6g}")
    d[key] = {"value": value, "unit": unit, "script": script, "run_id": run_id, "git_sha": git_sha(),
              "data_source": data_source, **({"note": note} if note else {})}
    NUMBERS.write_text(json.dumps(dict(sorted(d.items())), indent=1))


def get_number(key: str):
    d = json.loads(NUMBERS.read_text())
    return d[key]["value"]
