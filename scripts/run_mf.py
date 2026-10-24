"""P6 Module C: multi-fidelity correction of XFOIL drag with UIUC LSAT (public experimental) data.

Points: every clean UIUC drag-polar point with 50k <= Re <= 320k and Cl < 0.9 Clmax_exp(airfoil, Re), where XFOIL
(Ncrit 9) and NeuralFoil both give a pre-stall polar covering that Cl. Cd is compared at the experimental Cl
(UIUC drag polars are (Cl, Cd) pairs). Target: log(Cd_exp / Cd_XFOIL,N9). Models: ridge, GP (Matern + white),
HistGradientBoosting; model and hyperparameters chosen per outer fold by inner GroupKFold over training families.
Outer loop: leave-one-family-out over all non-test families. E387 and SD7037 (all model letters) are held out of
every training set and evaluated once at the end. Baselines on identical points: raw XFOIL N9, raw XFOIL N11,
constant mean log-bias (fit on the training families), NeuralFoil xxlarge (Ncrit 9).
"""
import json, sys, pathlib, time
from multiprocessing import Pool
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import pandas as pd
from flightloop.config import put_number, run_id as make_rid
from flightloop.mf import uiuc as U, xfoil as X

SCRIPT = "scripts/run_mf.py"
SRC = "public experimental (UIUC LSAT vol. 1 and 3) vs XFOIL/NeuralFoil"
TEST = {"E387", "SD7037"}
OUT = ROOT / "results/mf"; OUT.mkdir(parents=True, exist_ok=True)


def geom(xy):
    i = int(np.argmin(xy[:, 0])); up, lo = xy[: i + 1][::-1], xy[i:]
    x = np.linspace(0.005, 0.995, 199)
    yu = np.interp(x, up[:, 0], up[:, 1]); yl = np.interp(x, lo[:, 0], lo[:, 1])
    t, c = yu - yl, 0.5 * (yu + yl)
    return {"tc": float(t.max()), "xtc": float(x[t.argmax()]), "camber": float(c.max()), "xcamber": float(x[c.argmax()])}


def job(a):
    key, xy, re_, nc = a
    return (key, re_, nc), X.polar(xy, key, re_, nc)


def at_cl(rows, cl, cols):
    """Interpolate polar columns at Cl on the attached branch (from min alpha up to the alpha of max Cl)."""
    if len(rows) < 6:
        return None
    R = np.array(rows); i = int(np.argmax(R[:, 1])); B = R[: i + 1]
    if not (np.all(np.diff(B[:, 1]) > 0)):
        keep = np.r_[True, np.diff(np.maximum.accumulate(B[:, 1])) > 0]; B = B[keep]
    if len(B) < 4 or cl < B[0, 1] or cl > B[-1, 1]:
        return None
    return {c: float(np.interp(cl, B[:, 1], B[:, j])) for j, c in enumerate(cols)}


def nf_polar(xy, re_):
    import neuralfoil as nf
    al = np.arange(-4, 14.01, 0.25)
    r = nf.get_aero_from_coordinates(coordinates=xy, alpha=al, Re=re_, n_crit=9.0, model_size="xxlarge")
    return [[a, float(cl), float(cd), float(cf)] for a, cl, cd, cf in zip(al, r["CL"], r["CD"], r["analysis_confidence"])]


def build_points():
    P = U.polars(); C = U.coords()
    P = P[(P.re >= 50000) & (P.re <= 320000)].copy()
    P["has_xy"] = P.airfoil.map(lambda n: U.norm(n) in C)
    missing = sorted(P.loc[~P.has_xy, "airfoil"].unique().tolist())
    P = P[P.has_xy].copy()
    P["re_r"] = (P.re / 100).round() * 100
    groups = P.groupby(["airfoil", "re_r"])
    jobs = []
    for (af, re_), g in groups:
        xy = C[U.norm(af)]["xy"]
        for nc in (9, 11):
            jobs.append((U.norm(af), xy, float(re_), nc))
    t0 = time.time()
    with Pool(3) as pool:
        res = dict(pool.map(job, jobs))
    xf_time = time.time() - t0
    rows, excl = [], {"no_xfoil_n9": 0, "no_xfoil_n11": 0, "no_nf": 0, "post_stall": 0}
    for (af, re_), g in groups:
        k = U.norm(af); xy = C[k]["xy"]; gm = geom(xy)
        p9, p11 = res[(k, float(re_), 9)], res[(k, float(re_), 11)]
        nfp = nf_polar(xy, float(re_))
        clmax = g.cl.max()
        for _, r in g.iterrows():
            if r.cl >= 0.9 * clmax:
                excl["post_stall"] += 1; continue
            a9 = at_cl(p9["rows"], r.cl, p9["cols"]); a11 = at_cl(p11["rows"], r.cl, p11["cols"])
            an = at_cl(nfp, r.cl, ["alpha", "cl", "cd", "conf"])
            if a9 is None:
                excl["no_xfoil_n9"] += 1; continue
            if a11 is None:
                excl["no_xfoil_n11"] += 1; continue
            if an is None:
                excl["no_nf"] += 1; continue
            rows.append({"airfoil": af, "base": U.base_name(af), "family": U.family(af), "file": r.file, "re": r.re, "alpha_exp": r.alpha,
                         "cl": r.cl, "cd_exp": r.cd, "cd_xf9": a9["cd"], "cd_xf11": a11["cd"], "cd_nf": an["cd"], "nf_conf": an["conf"],
                         "alpha_xf": a9["alpha"], "cm_xf": a9["cm"], "xtr_top": a9["xtr_top"], "xtr_bot": a9["xtr_bot"], **gm})
    nonconv = {f"{k[0]}|{k[1]:.0f}|N{k[2]}": len(v["nonconverged_alpha"]) for k, v in res.items()}
    n_req = sum(v["n_requested"] for v in res.values()); n_got = sum(len(v["rows"]) for v in res.values())
    return pd.DataFrame(rows), {"excluded_points": excl, "airfoils_without_coords": missing, "xfoil_polars": len(res),
                                "xfoil_points_requested": n_req, "xfoil_points_converged": n_got, "xfoil_wall_s": xf_time,
                                "xfoil_hangs": int(sum(nonconv.values()))}


FEATS = ["logre", "alpha_xf", "cl", "logcd_xf", "cm_xf", "xtr_top", "xtr_bot", "tc", "xtc", "camber", "xcamber"]


def feats(df):
    d = df.copy(); d["logre"] = np.log10(d.re); d["logcd_xf"] = np.log(d.cd_xf9)
    return d[FEATS].to_numpy(), np.log(d.cd_exp / d.cd_xf9).to_numpy()


def candidates():
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import Ridge
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel
    from sklearn.ensemble import HistGradientBoostingRegressor
    c = {}
    for a in (0.1, 1.0, 10.0):
        c[f"ridge_a{a}"] = lambda a=a: make_pipeline(StandardScaler(), Ridge(alpha=a))
    for nu in (1.5, 2.5):
        c[f"gp_nu{nu}"] = lambda nu=nu: make_pipeline(StandardScaler(), GaussianProcessRegressor(
            ConstantKernel(0.1) * Matern(length_scale=np.ones(len(FEATS)), nu=nu) + WhiteKernel(0.01), normalize_y=True, random_state=0))
    for lr, d in ((0.05, 3), (0.1, None)):
        c[f"hgb_lr{lr}_d{d}"] = lambda lr=lr, d=d: HistGradientBoostingRegressor(learning_rate=lr, max_depth=d, max_iter=300, random_state=0)
    return c


def mape(cd_pred, cd_exp):
    return float(np.mean(np.abs(cd_pred - cd_exp) / cd_exp) * 100)


def fit_predict(make, Xtr, ytr, Xte, cap=900, seed=0):
    m = make()
    if len(Xtr) > cap and "gp" in str(type(m[-1] if hasattr(m, "__getitem__") else m)).lower():
        idx = np.random.default_rng(seed).choice(len(Xtr), cap, replace=False); Xtr, ytr = Xtr[idx], ytr[idx]
    m.fit(Xtr, ytr)
    return m.predict(Xte)


def select(dtr, cands):
    from sklearn.model_selection import GroupKFold
    Xa, ya = feats(dtr); g = dtr.family.to_numpy()
    k = min(5, len(set(g)))
    score = {}
    for name, make in cands.items():
        errs = []
        for tr, va in GroupKFold(k).split(Xa, ya, g):
            yp = fit_predict(make, Xa[tr], ya[tr], Xa[va])
            errs.append(mape(dtr.cd_xf9.to_numpy()[va] * np.exp(yp), dtr.cd_exp.to_numpy()[va]))
        score[name] = float(np.mean(errs))
    return min(score, key=score.get), score


def main():
    rid = make_rid("mf", "lofo")
    pts, meta = build_points()
    pts.to_csv(OUT / "points.csv", index=False)
    is_test = pts.base.str.replace(" ", "").isin(TEST)
    dev, test = pts[~is_test].reset_index(drop=True), pts[is_test].reset_index(drop=True)
    cands = candidates()
    fams = sorted(dev.family.unique())
    preds = np.full(len(dev), np.nan); bias_pred = np.full(len(dev), np.nan); chosen = {}
    leak_ok = True
    for f in fams:
        tr, te = dev[dev.family != f], dev[dev.family == f]
        leak_ok &= set(tr.family) .isdisjoint(set(te.family)) and not tr.base.str.replace(" ", "").isin(TEST).any()
        assert set(tr.family).isdisjoint(set(te.family)), f"family leak in fold {f}"
        assert not tr.base.str.replace(" ", "").isin(TEST).any(), "test airfoil in training"
        best, sc = select(tr, cands)
        Xtr, ytr = feats(tr); Xte, _ = feats(te)
        preds[te.index] = dev.cd_xf9.to_numpy()[te.index] * np.exp(fit_predict(cands[best], Xtr, ytr, Xte))
        bias_pred[te.index] = dev.cd_xf9.to_numpy()[te.index] * np.exp(np.mean(ytr))
        chosen[f] = {"model": best, "inner_mape": sc[best], "n_test": int(len(te))}
        print(f, best, round(sc[best], 1), len(te), flush=True)
    dev["cd_corr"] = preds; dev["cd_bias"] = bias_pred
    # final test airfoils: model selected and trained on all dev families
    best, sc = select(dev, cands)
    Xd, yd = feats(dev); Xt, _ = feats(test)
    test["cd_corr"] = test.cd_xf9 * np.exp(fit_predict(cands[best], Xd, yd, Xt)); test["cd_bias"] = test.cd_xf9 * np.exp(np.mean(yd))
    def block(d):
        return {"n_points": int(len(d)), "n_airfoils": int(d.airfoil.nunique()), "corrected": mape(d.cd_corr, d.cd_exp),
                "xfoil_n9": mape(d.cd_xf9, d.cd_exp), "xfoil_n11": mape(d.cd_xf11, d.cd_exp), "const_bias": mape(d.cd_bias, d.cd_exp),
                "neuralfoil": mape(d.cd_nf, d.cd_exp)}
    bands = {"60k": (50000, 80000), "100k": (80000, 160000), "200k-300k": (160000, 320001)}
    summ = {"run_id": rid, "meta": meta, "leakage_assert_passed": bool(leak_ok), "families": fams, "folds": chosen,
            "cd_mape": block(dev), "by_band": {b: block(dev[(dev.re >= lo) & (dev.re < hi)]) for b, (lo, hi) in bands.items()},
            "test_airfoils": {"model": best, **block(test), "by_airfoil_re": {}}, "nf_conf_median": float(dev.nf_conf.median())}
    for (b, re_), g in test.groupby(["base", (test.re / 1000).round()]):
        summ["test_airfoils"]["by_airfoil_re"][f"{b} Re {re_:.0f}k"] = {"n": int(len(g)), "xfoil_n9_cd_err_pct": float(np.mean((g.cd_xf9 - g.cd_exp) / g.cd_exp) * 100),
                                                                         "corrected_cd_err_pct": float(np.mean((g.cd_corr - g.cd_exp) / g.cd_exp) * 100),
                                                                         "neuralfoil_cd_err_pct": float(np.mean((g.cd_nf - g.cd_exp) / g.cd_exp) * 100)}
    dev.to_csv(OUT / "lofo_predictions.csv", index=False); test.to_csv(OUT / "test_predictions.csv", index=False)
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1, default=float))
    m = summ["cd_mape"]
    put_number("mf.cd_mape_corrected_pct", m["corrected"], "%", SCRIPT, SRC, rid, "leave-one-family-out, pooled")
    put_number("mf.cd_mape_xfoil_pct", m["xfoil_n9"], "%", SCRIPT, SRC, rid)
    put_number("mf.cd_mape_neuralfoil_pct", m["neuralfoil"], "%", SCRIPT, SRC, rid)
    put_number("mf.cd_mape_reduction_pct", 100 * (1 - m["corrected"] / m["xfoil_n9"]), "%", SCRIPT, SRC, rid, "1 - corrected/raw XFOIL N9")
    put_number("mf.n_points", m["n_points"], "-", SCRIPT, SRC, rid); put_number("mf.n_families", len(fams), "-", SCRIPT, SRC, rid)
    t = summ["test_airfoils"]
    put_number("mf.test.cd_mape_corrected_pct", t["corrected"], "%", SCRIPT, SRC, rid, "E387 + SD7037, untouched")
    put_number("mf.test.cd_mape_xfoil_pct", t["xfoil_n9"], "%", SCRIPT, SRC, rid, "E387 + SD7037, untouched")
    for k, v in t["by_airfoil_re"].items():
        if k.startswith("E387 Re 10"):
            put_number("mf.e387.re100k.cd_err_xfoil_pct", v["xfoil_n9_cd_err_pct"], "%", SCRIPT, SRC, rid, k)
            put_number("mf.e387.re100k.cd_err_corrected_pct", v["corrected_cd_err_pct"], "%", SCRIPT, SRC, rid, k)
    print(json.dumps({k: v for k, v in summ.items() if k not in ("folds",)}, indent=1, default=float)[:5000])


if __name__ == "__main__":
    main()
