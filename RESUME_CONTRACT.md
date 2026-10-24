# RESUME_CONTRACT

What may be claimed about FlightLoop, word for word, and the evidence behind each number. Every number is rendered from `results/numbers.json` (`scripts/build_docs.py`) and checked by `scripts/check_numbers.py` (V19). If a number changes, the bullet changes with it; nobody edits this file by hand.

## Bullets

1. Built a flight-test system-identification pipeline and a blind design-iteration engine for small UAVs: on 20 synthetic flights with turbulence and sensor noise it recovered lift and pitch derivatives within 2.0% (CLa), 2.0% (Cma) and 3.4% (Cmq) median error, and on 6 sealed blind cases it recovered injected drag growth within 2.9 percentage points and static-margin shift within 0.17% MAC; applied to 8 public CMU ALFA engine-failure flights.
2. Validated SU2 low-Reynolds CFD against NASA TM-4062 wind-tunnel data with a 3-grid Celik GCI (1.0% on drag) and a quantified Trust Map, and trained a leave-one-family-out drag correction for XFOIL that cut drag error from 12.6% to 9.5% MAPE across 18 airfoil families (1356 points; NeuralFoil 10.7%).

## Evidence and limits (say these if asked)

| Claim | Key in numbers.json | Check | Limit to state |
|---|---|---|---|
| synthetic sys-ID errors | sysid.synth.*_median_err_pct | V15 PASS | synthetic flights from this repo's simulator; 95% CI coverage failed its bar (CD0 25%, CLa 60%), see HONEST_LOG |
| blind loop recovery | loop.blind.max_cd0_delta_err_pp, loop.blind.max_sm_delta_err_mac_pct | V18 (recovery parts pass) | the next-flight prediction part of V18 failed: sink error cut 3.9%, trim 51.1% |
| ALFA flights | sysid.alfa.n_flights_used | V17 | public data flown by CMU, not by the author; mass and area are retailer figures; CD0 95% CI spans zero |
| GCI on drag | e387.gci.sa.cd_fine_pct | V05 | SA model; the SST-LM grids did not reach a steady state, so V05 is a FAIL overall |
| drag correction | mf.cd_mape_* | V11 | reduction is 24.3%, just short of its pre-registered bar; beats NeuralFoil on identical points |

Never claim: that an aircraft was built or flown for this project, that a wind tunnel was used for this project, or that SST-LM matched TM-4062 within its bars (its runs were unsteady).
