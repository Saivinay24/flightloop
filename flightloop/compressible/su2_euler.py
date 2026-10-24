"""SU2 EULER cases for the wedge and diamond verification. Base template: SU2 Tutorials
compressible_flow/Inviscid_Wedge/inv_wedge_HLLC.cfg (v8); changed keys listed in cfd/templates/DIFFS.md.

Stopping rule: Cauchy on CD (relative change < 1e-7 over 500 iterations) or the 6000-iteration cap. The first 7 runs
used RMS_DENSITY < -11; on the finer grids the Venkatakrishnan limiter holds the residual in a limit cycle near -4.7
while CL/CD are stationary to 1e-9, so those runs burned the full cap (HONEST_LOG 11). post_compressible.py reports
both the residual drop and the force stationarity of every run."""
import json, math, os, re, subprocess, pathlib
from ..config import ROOT, git_sha, sha_file

TEMPLATE = ROOT / "cfd/templates/wedge/inv_wedge_HLLC.cfg"
SU2 = pathlib.Path.home() / "opt/su2/bin"


def build(mesh, mach, aoa, wall, iters=6000):
    ov = {"MACH_NUMBER": f"{mach}", "AOA": f"{aoa}", "FREESTREAM_PRESSURE": "101325.0", "FREESTREAM_TEMPERATURE": "288.15",
          "MARKER_EULER": f"( {wall} )", "MARKER_FAR": "( farfield )", "MARKER_PLOTTING": f"( {wall} )",
          "MARKER_MONITORING": f"( {wall} )", "MESH_FILENAME": mesh, "MESH_FORMAT": "SU2", "ITER": str(iters),
          "CONV_NUM_METHOD_FLOW": "ROE", "MUSCL_FLOW": "YES", "SLOPE_LIMITER_FLOW": "VENKATAKRISHNAN", "VENKAT_LIMITER_COEFF": "0.05",
          "MGLEVEL": "0", "CFL_NUMBER": "5.0", "CFL_ADAPT": "YES", "CFL_ADAPT_PARAM": "( 0.5, 1.2, 1.0, 50.0 )",
          "CONV_FIELD": "DRAG", "CONV_CAUCHY_ELEMS": "500", "CONV_CAUCHY_EPS": "1E-7", "CONV_STARTITER": "500", "REF_DIMENSIONALIZATION": "DIMENSIONAL",
          "OUTPUT_FILES": "( RESTART, PARAVIEW, SURFACE_PARAVIEW )", "VOLUME_OUTPUT": "( COORDINATES, SOLUTION, PRIMITIVE )",
          "HISTORY_OUTPUT": "( ITER, RMS_RES, AERO_COEFF )", "SCREEN_OUTPUT": "( INNER_ITER, RMS_DENSITY, LIFT, DRAG )",
          "SCREEN_WRT_FREQ_INNER": "100", "OUTPUT_WRT_FREQ": "2000", "SOLUTION_FILENAME": "restart_flow"}
    skip = {"MARKER_SUPERSONIC_INLET", "MARKER_OUTLET", "DV_KIND", "DV_PARAM", "DV_VALUE"}
    lines, seen = [], set()
    for line in TEMPLATE.read_text().splitlines():
        m = re.match(r"^\s*([A-Z0-9_]+)\s*=", line)
        if m and m.group(1) in skip:
            continue
        if m and m.group(1) in ov:
            lines.append(f"{m.group(1)}= {ov[m.group(1)]}"); seen.add(m.group(1))
        else:
            lines.append(line)
    lines += [f"{k}= {v}" for k, v in ov.items() if k not in seen]
    return "\n".join(lines) + "\n"


def run(case, level, mach, aoa, rid):
    d = ROOT / "runs" / rid; d.mkdir(parents=True, exist_ok=True)
    mesh = ROOT / f"cfd/meshes/supersonic/{case}_L{level}.su2"
    if not (d / "mesh.su2").exists():
        os.symlink(mesh, d / "mesh.su2")
    if (d / "exit_code").exists() and (d / "exit_code").read_text().strip() == "0":
        return d  # finished runs are never rewritten (provenance: case.cfg must match what ran)
    wall = "wall" if case == "wedge" else "airfoil"
    (d / "case.cfg").write_text(build("mesh.su2", mach, aoa, wall))
    (d / "case.json").write_text(json.dumps({"case": case, "level": level, "mach": mach, "aoa": aoa, "run_id": rid,
                                             "mesh_sha256_16": sha_file(mesh), "git_sha": git_sha()}, indent=1))
    if (d / "exit_code").exists() and (d / "exit_code").read_text().strip() == "0":
        return d
    subprocess.run(["bash", "-c", f"cd '{d}' && {SU2}/SU2_CFD case.cfg > log.txt 2>&1; echo $? > exit_code"])
    return d
