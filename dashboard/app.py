"""FlightLoop dashboard: `streamlit run dashboard/app.py`.

Reads the same generated bundle as the static site (site/data/bundle.js, built by scripts/build_site.py from
results/*.json and results/numbers.json), so the dashboard, the site and the report cannot disagree.
Every tile names its data class: public experimental, CFD, low-order, synthetic.
"""
import json, pathlib
import pandas as pd
import streamlit as st

ROOT = pathlib.Path(__file__).resolve().parents[1]
txt = (ROOT / "site/data/bundle.js").read_text()
FL = json.loads(txt[txt.index("window.FL = ") + len("window.FL = "):].rstrip().rstrip(";"))
NUM = FL["numbers"]
VER = {r["id"]: r for r in FL["verify"].get("rows", [])}


def n(key, d=2):
    e = NUM.get(key)
    if e is None:
        return "n/a (not produced)"
    v = e["value"]
    return f"{v:.{d}f}" if isinstance(v, float) else str(v)


def check(vid):
    r = VER.get(vid)
    return "not run" if r is None else f"{vid} {'PASS' if r['pass'] else 'FAIL'}"


st.set_page_config(page_title="FlightLoop", layout="wide")
st.title("FlightLoop")
st.caption(f"Bundle {FL['meta']['generated_utc']}, commit {FL['meta']['git_sha']}, numbers hash {FL['meta']['numbers_hash']}. "
           "Failed pre-registered bars are shown as failed with the real number.")
tabs = st.tabs(["E387 Validation", "Trust Map", "Stratospheric Re", "Multi-fidelity", "Compressible", "Sys-ID", "Correlation Report", "Checks"])

with tabs[0]:
    st.caption("CFD (SU2) vs NASA TM-4062 (public experimental)")
    V = FL["cfd"]["validation"]
    rows = []
    for a in (0, 4, 8):
        for m, e in V.get(f"re200k_a{a}", {}).items():
            rows.append({"alpha": a, "model": m, "state": e["state"], "Cl": e["cl"], "Cl exp": e["exp"]["cl"], "Cl err %": e["cl_err_pct"],
                         "Cd": e["cd"], "Cd exp": e["exp"]["cd"], "Cd err %": e["cd_err_pct"], "RMS dCp upper": e["cp_rms_upper"]})
    st.dataframe(pd.DataFrame(rows), width="stretch")
    a = st.radio("alpha (deg)", [0, 4, 8], index=1, horizontal=True)
    cp = {}
    for m, e in V.get(f"re200k_a{a}", {}).items():
        cp["x/c"] = e["cp_x"]; cp["TM-4062"] = e["cp_exp"]; cp[m] = e["cp_cfd"]
    if cp:
        df = pd.DataFrame(cp).set_index("x/c")
        st.line_chart(-df, x_label="x/c", y_label="-Cp (upper)")
    st.subheader("GCI (alpha 4)")
    g = FL["cfd"]["gci"]
    st.dataframe(pd.DataFrame([{"model": m, "qty": q, "type": g[m][q]["type"], "p": g[m][q]["p"], "GCI fine %": g[m][q]["gci_fine_pct"],
                                "fine/medium/coarse": g[m][q]["f"]} for m in ("SA", "SST_LM") if m in g for q in ("cl", "cd")]), width="stretch")
    st.write(f"{check('V05')}, {check('V06')}, {check('V07')}, {check('V08')}")

with tabs[1]:
    st.caption("CFD runs classified against TM-4062 bars; unsteady runs listed, excluded from error statistics")
    st.dataframe(pd.DataFrame(FL["cfd"]["trust"]["cells"]), width="stretch")
    for x in FL["cfd"]["limits"]:
        st.markdown(f"- {x}")

with tabs[2]:
    st.caption("SD7037: CFD, low-order (XFOIL, NeuralFoil, corrected XFOIL) vs UIUC LSAT (public experimental)")
    if FL.get("strato"):
        st.dataframe(pd.DataFrame(FL["strato"]["rows"]), width="stretch")
    else:
        st.write("n/a (scripts/run_strato.py not run)")

with tabs[3]:
    st.caption("Low-order (XFOIL) corrected with a model trained on UIUC (public experimental); leave-one-family-out")
    c = st.columns(4)
    c[0].metric("Corrected Cd MAPE %", n("mf.cd_mape_corrected_pct")); c[1].metric("Raw XFOIL Cd MAPE %", n("mf.cd_mape_xfoil_pct"))
    c[2].metric("NeuralFoil Cd MAPE %", n("mf.cd_mape_neuralfoil_pct")); c[3].metric("Check", check("V11"))
    t = pd.DataFrame(FL["mf_test"])
    if len(t):
        st.scatter_chart(t.melt(id_vars=["cd_exp"], value_vars=["cd_xf9", "cd_corr", "cd_nf"], var_name="tool", value_name="cd_pred"),
                         x="cd_exp", y="cd_pred", color="tool")
    st.dataframe(pd.DataFrame(FL["mf"]["by_band"]).T, width="stretch")

with tabs[4]:
    st.caption("CFD (SU2 Euler) vs oblique-shock and shock-expansion theory; textbook verification only")
    rows = []
    for k, w in FL["comp"].get("wedge", {}).items():
        rows.append({"case": f"wedge {k}", "beta fine": w["fine"]["beta_deg"], "beta theory": w["theory"]["beta_deg"], "p2/p1 fine": w["fine"]["p2_p1"], "p2/p1 theory": w["theory"]["p2_p1"]})
    st.dataframe(pd.DataFrame(rows), width="stretch")
    st.dataframe(pd.DataFrame([{"case": k, "Cl fine": d["fine"]["cl"], "Cl theory": d["theory"]["cl"], "Cd fine": d["fine"]["cd"], "Cd theory": d["theory"]["cd"],
                                "Cd err %": d.get("fine_err_cd_pct")} for k, d in FL["comp"].get("diamond", {}).items()]), width="stretch")
    st.write(f"{check('V12')}, {check('V13')}")

with tabs[5]:
    st.caption("Synthetic 20-seed Monte Carlo (labelled synthetic) and CMU ALFA engine-failure glides (public experimental)")
    S = FL["synth"]
    st.dataframe(pd.DataFrame({m: {k: v["median_abs_err_pct"] for k, v in S[m].items()} for m in ("main", "output_error", "alt18km") if m in S}), width="stretch")
    st.write(f"95% CI coverage: {S['coverage_pct']} ({check('V16')}); {check('V15')}")
    A = FL.get("alfa")
    if A:
        st.dataframe(pd.DataFrame(A["estimates"]), width="stretch")
        st.dataframe(pd.DataFrame(A["avl_vs_flight"]).astype(str), width="stretch")
        for x in A["approximations"]:
            st.markdown(f"- {x}")

with tabs[6]:
    st.caption("StratoGlide test-card predictions vs a SYNTHETIC blind-loop flight (not a flight result)")
    if FL.get("corr"):
        st.dataframe(pd.DataFrame(FL["corr"]["rows"]).astype(str), width="stretch")
    if FL.get("loop"):
        st.subheader("Blind loop (synthetic)")
        st.dataframe(pd.DataFrame([{k: v for k, v in s.items() if k not in ("per_speed", "secondary_post_hoc")} for s in FL["loop"]["seeds"]]), width="stretch")
        st.write(check("V18"))
    if FL.get("deltas"):
        st.dataframe(pd.DataFrame(FL["deltas"]["rules"]).astype(str), width="stretch")

with tabs[7]:
    st.dataframe(pd.DataFrame(FL["verify"].get("rows", [])), width="stretch")
    for x in FL["honest"]["failed"]:
        st.markdown(f"- {x}")
