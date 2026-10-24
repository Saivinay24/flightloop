"""SU2 case writer + resumable launcher for 2D airfoil RANS (INC_RANS, non-dimensional: c = U = rho = 1, mu = 1/Re).

Base template: SU2 Tutorials incompressible_flow/Inc_Turbulent_NACA0012/turb_naca0012.cfg (v8). Only the
keys in OVERRIDES / model blocks are changed; the diff is documented in cfd/templates/DIFFS.md.
Every run lives in runs/<run_id>/ with cfg snapshot, mesh SHA, log, history.csv, surface_flow.csv, status.json.
Resume: if restart_flow.dat exists the case restarts from it (RESTART_SOL= YES) and keeps the iteration count.
"""
from __future__ import annotations

import json
import math
import os
import pathlib
import re as _re
import shutil
import subprocess

from ..config import ROOT, cfg, git_sha, run_id as make_run_id, sha_file

TEMPLATE = ROOT / "cfd/templates/naca0012_inc/turb_naca0012.cfg"
SU2_BIN = pathlib.Path.home() / "opt/su2/bin"

MODELS = {
    "SA": {"KIND_TURB_MODEL": "SA", "SA_OPTIONS": "NONE", "KIND_TRANS_MODEL": "NONE"},
    "SA_BC": {"KIND_TURB_MODEL": "SA", "SA_OPTIONS": "BCM", "KIND_TRANS_MODEL": "NONE"},
    "SST": {"KIND_TURB_MODEL": "SST", "SST_OPTIONS": "V2003m, SUSTAINING", "KIND_TRANS_MODEL": "NONE"},
    "SST_LM": {"KIND_TURB_MODEL": "SST", "SST_OPTIONS": "V2003m, SUSTAINING", "KIND_TRANS_MODEL": "LM", "LM_OPTIONS": "NONE"},
}


def build_cfg(mesh_file: str, alpha_deg: float, re: float, model: str, tu_pct: float, iters: int, restart: bool) -> str:
    a = math.radians(alpha_deg)
    ov = {
        "INC_DENSITY_INIT": "1.0",
        "INC_VELOCITY_INIT": f"( {math.cos(a):.12f}, {math.sin(a):.12f}, 0.0 )",
        "MU_CONSTANT": f"{1.0 / re:.10e}",
        "RESTART_SOL": "YES" if restart else "NO",
        "MESH_FILENAME": mesh_file,
        "SOLUTION_FILENAME": "restart_flow",
        "ITER": str(iters),
        "CFL_NUMBER": "10.0",
        "CFL_ADAPT": "YES",
        "CFL_ADAPT_PARAM": "( 0.5, 1.05, 5.0, 200.0, 0.001 )",
        "NUM_METHOD_GRAD": "WEIGHTED_LEAST_SQUARES",
        "LINEAR_SOLVER_ERROR": "1E-6",
        "LINEAR_SOLVER_ITER": "20",
        "MUSCL_TURB": "NO",
        "CONV_FIELD": "( DRAG, LIFT )",
        "CONV_RESIDUAL_MINVAL": "-14",
        "CONV_STARTITER": "2000",
        "CONV_CAUCHY_ELEMS": "1000",
        "CONV_CAUCHY_EPS": "1E-7",
        "FREESTREAM_TURBULENCEINTENSITY": f"{tu_pct / 100.0:.6f}",
        "FREESTREAM_TURB2LAMVISCRATIO": "1.0" if model.startswith("SA") else "10.0",
        "OUTPUT_FILES": "( RESTART, SURFACE_CSV, SURFACE_PARAVIEW, PARAVIEW )",
        "VOLUME_OUTPUT": "( COORDINATES, SOLUTION, PRIMITIVE )",
        "HISTORY_OUTPUT": "( ITER, RMS_RES, AERO_COEFF, CFL_NUMBER, WALL_TIME )",
        "SCREEN_OUTPUT": "( INNER_ITER, WALL_TIME, RMS_PRESSURE, LIFT, DRAG, MOMENT_Z )",
        "OUTPUT_WRT_FREQ": "500",
        "SCREEN_WRT_FREQ_INNER": "50",
        "WRT_RESTART_OVERWRITE": "YES",
        "WRT_FORCES_BREAKDOWN": "YES",
    }
    ov.update(MODELS[model])
    lines, seen = [], set()
    for line in TEMPLATE.read_text().splitlines():
        m = _re.match(r"^\s*([A-Z0-9_]+)\s*=", line)
        if m and m.group(1) in ov:
            k = m.group(1)
            if k in seen:
                continue
            lines.append(f"{k}= {ov[k]}")
            seen.add(k)
        elif m and m.group(1) == "SCREEN_OUTPUT":
            continue
        else:
            lines.append(line)
    for k, v in ov.items():
        if k not in seen:
            lines.append(f"{k}= {v}")
    return "\n".join(lines) + "\n"


def prepare(case: dict, rid: str | None = None) -> pathlib.Path:
    """case: {airfoil, level, alpha, re, model, tu_pct, iters}"""
    name = f"{case['airfoil']}-re{int(case['re'] / 1000)}k-a{case['alpha']:g}-{case['model']}-{case['level']}" + \
           (f"-tu{case['tu_pct']:g}" if case.get("tag_tu") else "")
    rid = rid or make_run_id("cfd", name)
    d = ROOT / "runs" / rid
    d.mkdir(parents=True, exist_ok=True)
    mesh = ROOT / f"cfd/meshes/{case['airfoil']}/{case['airfoil']}_{case['level']}.su2"
    link = d / "mesh.su2"
    if not link.exists():
        os.symlink(mesh, link)
    restart = (d / "restart_flow.dat").exists()
    (d / "case.cfg").write_text(build_cfg("mesh.su2", case["alpha"], case["re"], case["model"], case["tu_pct"], case["iters"], restart))
    meta = {**case, "run_id": rid, "mesh": str(mesh.relative_to(ROOT)), "mesh_sha256_16": sha_file(mesh),
            "git_sha": git_sha(), "su2_version": "v8.6.0"}
    (d / "case.json").write_text(json.dumps(meta, indent=1))
    st = d / "status.json"
    if not st.exists():
        st.write_text(json.dumps({"state": "prepared"}))
    return d


def launch(d: pathlib.Path, np_: int = 2, background: bool = True):
    env = dict(os.environ, PATH=f"{SU2_BIN}:{os.environ['PATH']}")
    if (d / "history.csv").exists() and (d / "restart_flow.dat").exists():
        n = 0
        while (d / f"history_part{n}.csv").exists():
            n += 1
        shutil.move(d / "history.csv", d / f"history_part{n}.csv")
    cmd = f"cd '{d}' && mpirun -np {np_} {SU2_BIN}/SU2_CFD case.cfg > log.txt 2>&1; echo $? > exit_code"
    (d / "status.json").write_text(json.dumps({"state": "running"}))
    if background:
        p = subprocess.Popen(["bash", "-c", cmd], start_new_session=True)
        return p.pid
    subprocess.run(["bash", "-c", cmd])
    return None
