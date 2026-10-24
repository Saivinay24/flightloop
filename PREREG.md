# PREREG: pre-registered success bars for FlightLoop

Committed before any result run. verify.py checks that the commit adding this file predates every commit that touches `runs/` or `results/numbers.json`. Bars are never edited after results exist. A failed [R] bar is reported with its real number in HONEST_LOG.md.

## Bars (copied verbatim from MASTER_PLAN.md section 6)

| # | Class | Check | Bar |
|---|---|---|---|
| V01 | I | ISA module | rho(11 km) = 0.3639 +/-0.5%; rho(20 km) = 0.0889 +/-0.5%; nu(18 km) = 1.18e-4 +/-1% |
| V02 | I | Provenance | every file in `data/` has a SOURCES.md row; TM-4062 double transcription diff = 0; every synthetic file has `source: synthetic` |
| V03 | I | Mesh | max y+ < 1 on every reported E387/SD7037 run; per-direction refinement ratio between levels in [1.3, 1.6]; mesh SHA logged per run |
| V04 | I | Convergence | every run in the error statistics meets the residual or Cauchy criterion; others are status `unsteady` and appear only in the Trust Map |
| V05 | R | GCI, E387 Re 200k alpha 4, SA and SST-LM | GCI_fine(Cl) < 2%, GCI_fine(Cd) < 5%, observed order p in [0.5, 3]; oscillatory reported as such |
| V06 | R | E387 Re 200k alpha 4, SST-LM vs TM-4062 | |Cl error| < 5%; |Cd error| < 20%; upper-surface RMS delta Cp < 0.05 |
| V07 | R | Fully turbulent penalty, same point | SA Cd over-predicts TM-4062 by > 25% (pre-registered direction and size; if it fails, the physics story changes and HONEST_LOG says how) |
| V08 | R | LSB, Re 200k, alpha 0 and 4, SST-LM vs oil flow | |x_s error| < 0.05 c and |x_r error| < 0.08 c |
| V09 | I | Limits | Trust Map covers every (Re, alpha, model) run with a status; LIMITS.md has >= 3 failure regions each with a number |
| V10 | I | Stratospheric case | DISAGREEMENT.md table complete for all tools at Re 60k and 100k, every cell a number or `n/a (reason)` |
| V11 | R | Multi-fidelity, leave-one-family-out | corrected Cd MAPE <= 0.75 x raw-XFOIL MAPE AND < NeuralFoil MAPE on identical points; code asserts no family in both train and test of any fold (that assertion is [I]) |
| V12 | R | Wedge, M 2 and 3, theta 10 deg | |beta - theory| < 0.5 deg (39.31, 27.38); p2/p1 within 1% (1.707, 2.055) |
| V13 | R | Diamond, M 2, alpha 2 deg | Cl within 2% of 0.0810 and Cd within 2% of 0.00863 on the fine grid; monotone convergence over 3 grids |
| V14 | I | Parsers | synthetic Parquet round-trip exact; >= 1 PX4 public ULog and >= 1 ArduPilot .bin parse to the full required schema; timestamps monotone; resample to 50 Hz |
| V15 | R | Synthetic sys-ID, 20 seeds, turbulence + noise | median |error|: CD0 < 5%, k < 10%, CLa < 5%, Cma < 10%, Cmde < 10%, Cmq(+Cm_alpha_dot) < 25% |
| V16 | R | CI calibration | 95% CI covers truth in 80-100% of seeds for CD0, CLa, Cma |
| V17 | I | Real-data sys-ID | ALFA results list flights and segments used, every estimate has a CI and residual RMS, AVL-vs-flight table present, approximations labelled |
| V18 | R | Loop closure, >= 5 blind seeds | CD0 delta recovered within 5 percentage points; SM delta within 1% MAC; next-flight error in sink rate and trim elevator reduced >= 70%; engine never opened `SEALED_truth.json` ([I] sub-check) |
| V19 | I | Outward surfaces | report PDF, dashboard and static HTML exist; zero U+2014 em dashes in README, report, dashboard strings, resume contract; every number in the resume contract and report equals `results/numbers.json` |
| V20 | I | Hygiene | secrets grep clean; no Co-Authored-By/Claude/Anthropic in `git log`; data labelled public/synthetic; grep for `warhead|seeker|guidance law|homing|ramjet|intake` returns nothing in code, data, report or dashboard (MASTER_PLAN.md, PREREG.md and verify.py are excluded because they state the negative scope) |

## Fixed analysis choices (decided before any result)

1. **Headline LM inlet turbulence:** Tu = 0.06% (TM-4062 measured LTPT turbulence at Re 100k and 200k, Pt 15 psi). Freestream eddy viscosity ratio 1 (SA_BC) / 10 (SST family). The Tu sweep {0.05, 0.1, 0.2}% is reported as a sensitivity and never replaces the headline.
2. **Experimental reference:** TM-4062 Table B1 (corrected), runs 9/10/13 for Re 200k, runs 15/16 for Re 100k. Where duplicate alpha rows exist the value is linearly interpolated in alpha from the full table (duplicates averaged first). Cp from Appendix D (uncorrected, as published), upper surface, alpha columns 0.01, 3.99, 8.00 deg for Re 200k. Oil flow from Table III.
3. **Alpha:** CFD runs at the nominal alpha {0, 4, 8}; no offset is applied. Experimental Cl/Cd are interpolated to the CFD alpha.
4. **Error definitions:** Cl error % = 100 (Cl_cfd - Cl_exp) / Cl_exp; Cd error % likewise (signed; bars use the absolute value except V07, which requires SA over-prediction > +25%). RMS delta Cp over the 29 upper-surface orifices with CFD Cp linearly interpolated in x/c.
5. **LSB:** separation = first upstream-to-downstream zero crossing of the wall-tangential skin friction to negative on the upper surface aft of the suction peak; reattachment = next return to positive. Compared with Table III LS and TR.
6. **Convergence:** RMS pressure residual drop >= 6 orders, or Cauchy relative change of Cl and Cd < 1e-6 over 1000 iterations (implemented as max relative deviation over the last 1000 iterations < 1e-5 when the solver does not report Cauchy). Otherwise status `unsteady` and excluded from error statistics.
7. **GCI:** Celik et al. 2008 (J. Fluids Eng. 130, 078001), Fs = 1.25, r from cell counts sqrt(N1/N2) in 2D.
8. **V11 MAPE:** Cd MAPE over identical (airfoil, Re, Cl) points, held-out family per fold; reduction = 1 - corrected/raw(Ncrit 9).
9. **V15/V16:** truth = parameters of the simulation model (StratoGlide-derived design model, section 4.3 contract); 20 seeds, Dryden light turbulence on, sensor noise per config.
10. **V18:** blind seeds drawn by `stubs/inject.py`; error reduction = 1 - |err_after| / |err_before| for sink rate and trim elevator, median over seeds; pass requires both >= 70%.

## Addendum A (written before any sys-ID or loop result; no bar changed)

11. **V15/V16 estimators:** CD0 and k from `sysid/glide.py` *force* method (specific force rotated by vane alpha, window means, OLS) on quasi-steady windows; the energy method is reported alongside as a sensitivity (it is biased by vertical air motion by construction). CLa, Cma, Cmde, Cmq from `sysid/equation_error.py` (vane alpha). Output-error results are reported, not used for the bar. CIs for V16: CD0 from the window-regression SE, CLa and Cma from the Morelli-Klein corrected SE.
12. **Truth for V15/V16:** StratoGlide design model `inputs/stratoglide/derivs.json` with the stall blend narrowed to 0.5 deg (linear inside the alpha <= 9.6 deg flown envelope). Prandtl-Glauert factor at the flown Mach (about 0.04 to 0.06) changes CLa by < 0.2% and is ignored.
13. **Simulation conditions:** sea level start 900 m, speeds 13.5 to 21 m/s, MIL-F-8785C light low-altitude intensities (u and w first-order), sensor noise per config incl. 0.5 deg vane noise, pitot scale error N(1, 0.01) per flight.
