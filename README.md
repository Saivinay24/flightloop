# FlightLoop

Low-Reynolds-number CFD validated against wind-tunnel data, flight-log system identification, and a blind next-iteration design engine for small UAVs. It takes StratoGlide (a balloon-launched glider design, used read-only at tag `stratoglide-v1`) as its design input and closes the loop from design prediction to flight data to design change.

Every number below is rendered from `results/numbers.json` by `scripts/build_docs.py`; `scripts/check_numbers.py` fails if a document drifts from its template or carries a hand-typed result. Failed pre-registered bars are reported as failed, with the real number, in `HONEST_LOG.md`.

**Data classes.** Public experimental: NASA TM-4062 (E387, Langley LTPT), UIUC LSAT airfoil polars, CMU ALFA flight logs (CC0), PX4 Flight Review logs. Synthetic: every flight produced by this repo's simulator, labelled `source: synthetic` in its metadata. No aircraft was built or flown for this project and no wind tunnel was used; the experimental data is published by the institutions named above.

## Results at a glance

| Module | Result | Check |
|---|---|---|
| A. E387 Re 200k, SU2 SA (fully turbulent) | Cd +25.1% vs TM-4062 at alpha 4, Cl -0.13% | V07 |
| A. E387 Re 200k, SU2 SST-LM (transition) | Cd +4.4%, Cl +0.34%, bubble x_s +0.013c, x_r +0.021c, state unsteady | V06, V08 |
| A. Grid convergence (SA, 3 grids) | GCI Cl 0.32%, Cd 1.00% | V05 |
| B. SD7037 stratospheric Re | XFOIL Ncrit 9 Cd -8.8% vs corrected -1.1% at Re 100k (signed mean over CL 0.5 to 0.8) | V10 |
| C. Multi-fidelity Cd correction | MAPE 9.54% vs raw XFOIL 12.59% vs NeuralFoil 10.71%, 1356 points, 18 families, leave-one-family-out | V11 |
| D. Supersonic wedge / diamond | shock angle +0.004 / +0.012 deg (M2 / M3), diamond Cd -1.85% | V12, V13 |
| E. Sys-ID, synthetic (20 seeds) | median error CLa 1.96%, Cma 2.00%, Cmq 3.42%, CD0 1.45% | V15, V16 |
| E. Sys-ID, CMU ALFA (8 flights) | CD0 0.0212, k 0.159 (AVL 0.0574) | V17 |
| F. Blind loop (6 sealed seeds) | CD0 delta error max 2.89 pp, SM delta error max 0.17% MAC, next-flight trim error cut 51.1%, sink 3.9% | V18 |

Full verifier output: see `report/FlightLoop_Report.md`, section 1.

## Run it

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.lock   # SU2 v8.6.0 and XFOIL: see ENV.md
.venv/bin/python verify.py                    # 20 numbered checks, writes results/verify_report.json
.venv/bin/python scripts/build_site.py        # site/data/bundle.js from results/
.venv/bin/python scripts/build_docs.py        # README, RESUME_CONTRACT, report .md and .pdf
.venv/bin/streamlit run dashboard/app.py      # dashboard on the same bundle
open site/index.html                          # static page, works from file:// or any subpath
```

Pipeline order (each script writes its own results and numbers): `scripts/cfd_queue.py cutline` then `scripts/analyze_e387.py` (A), `scripts/run_compressible.py` then `scripts/post_compressible.py` (D), `scripts/run_sysid_synth.py` (E synthetic), `scripts/run_loop_blind.py` and `scripts/correlation_report.py` (F), `scripts/run_mf.py` (C), `scripts/cfd_queue.py strato` then `scripts/run_strato.py` (B), `scripts/parse_logs.py` and `scripts/run_alfa.py` (E real), `scripts/compute_used.py`.

## Layout

```
flightloop/   cfd (SU2 runner, GCI, VTU), mesh, compressible, sim (6-DOF, Dryden, sensors), sysid, loop, mf, logs
scripts/      one script per pipeline stage; scripts/diag holds the root-cause diagnostics cited in HONEST_LOG.md
results/      numbers.json (every number, with script, run id, git sha, data source), per-module JSON, figures
runs/<id>/    one folder per run: case, config, mesh hash, history, status
site/         static page; dashboard/ Streamlit app; report/ Markdown and PDF
PREREG.md     bars committed before any result; HONEST_LOG.md, LIMITS.md, GAPS.md, COMPUTE.md, ENV.md
```

Compute: 43.8 core-hours on a MacBook Air M4 (42 SU2 runs including superseded ones), 0 cloud hours.
