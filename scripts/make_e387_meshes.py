"""Build coarse/medium/fine E387 C-meshes, log quality + SHA to results/mesh_e387.json, plot the medium mesh."""
import json, sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from flightloop.config import cfg, ROOT, sha_file
from flightloop.geometry import Airfoil, read_selig, with_le_point
from flightloop.mesh import cmesh

airfoil = sys.argv[1] if len(sys.argv) > 1 else "e387"
C = cfg()["mesh"]
af = Airfoil(with_le_point(read_selig(ROOT / f"data/airfoils/{airfoil}.dat")))
out = ROOT / f"cfd/meshes/{airfoil}"; out.mkdir(parents=True, exist_ok=True)
lv = C["levels"]
beta = cmesh.normal_map_beta(C["farfield_radius_c"], C["first_cell_height_c"], lv["medium"]["Nn"] - 1)
log = {}
for name, L in lv.items():
    X, meta = cmesh.build(af, L["Ns"], L["Nw"], L["Nn"], R=C["farfield_radius_c"], W=C["wake_length_c"], beta=beta)
    q = cmesh.quality(X)
    path = out / f"{airfoil}_{name}.su2"
    npts = cmesh.write_su2(X, L["Nw"], L["Ns"], path)
    wrap = X[0, L["Nw"] - 1: L["Nw"] - 1 + L["Ns"]]
    up = wrap[np.argmin(wrap[:, 0]):]
    sel = up[(up[:, 0] > 0.2) & (up[:, 0] < 0.8)]
    ds_mid = float(np.mean(np.hypot(*np.diff(sel, axis=0).T)))
    rad = float(np.min(np.hypot(X[-1, :, 0] - 0.5, X[-1, :, 1])))
    log[name] = {**meta, **q, "points": npts, "sha256_16": sha_file(path), "upper_mid_spacing_c": ds_mid,
                 "farfield_min_dist_from_midchord_c": rad}
    print(name, json.dumps(log[name]))
    np.save(out / f"{airfoil}_{name}_X.npy", X[:3])   # first layers (for y+ / wall geometry checks)
    if name == "medium":
        fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
        for ax, lim in zip(axs, [((-1, 2.5), (-1.2, 1.2)), ((-0.02, 0.06), (-0.03, 0.03)), ((0.9, 1.08), (-0.04, 0.04))]):
            for k in range(0, X.shape[0], 1):
                ax.plot(X[k, :, 0], X[k, :, 1], "k-", lw=0.2)
            for i in range(0, X.shape[1], 1):
                ax.plot(X[:, i, 0], X[:, i, 1], "k-", lw=0.2)
            ax.set_xlim(*lim[0]); ax.set_ylim(*lim[1]); ax.set_aspect("equal")
        fig.suptitle(f"{airfoil.upper()} medium C-mesh (CFD grid, {q['cells']} quads)")
        fig.savefig(ROOT / f"results/figures/mesh_{airfoil}_medium.png", dpi=130, bbox_inches="tight")
if airfoil == "e387":
    # cross-check against TM-4062 Table I design ordinates (text layer, PDF p. 14, first 16 rows each side)
    up = [(0.00043,0.00233),(0.00518,0.00932),(0.01423,0.01727),(0.02748,0.02562),(0.04493,0.03408),(0.06643,0.04238),(0.09185,0.05033),(0.12093,0.05775),(0.15345,0.06448),(0.18907,0.07037),(0.22742,0.07528),(0.26813,0.07908),(0.31078,0.08157),(0.35505,0.08247),(0.40077,0.08173),(0.44767,0.07937)]
    lo = [(0.00092,-0.00287),(0.00717,-0.00682),(0.01890,-0.01017),(0.03597,-0.01265),(0.05827,-0.01425),(0.08568,-0.01500),(0.11800,-0.01502),(0.15490,-0.01442),(0.19598,-0.01328),(0.24083,-0.01177),(0.28892,-0.00998),(0.33968,-0.00803),(0.39252,-0.00605),(0.44678,-0.00410),(0.50182,-0.00228),(0.55693,-0.00065)]
    dz = [abs(af.sy(af.upper_t(x)) - z) for x, z in up] + [abs(af.sy(af.lower_t(x)) - z) for x, z in lo]
    log["geometry_check"] = {"max_abs_dz_vs_tm4062_design_c": float(max(dz)), "n_points": len(dz),
                             "source": "UIUC e387.dat + LE (0,0) from TM-4062 Table I"}
    print("geometry max |dz| vs TM-4062 design:", max(dz))
(ROOT / f"results/mesh_{airfoil}.json").write_text(json.dumps(log, indent=1))
