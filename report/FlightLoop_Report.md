# FlightLoop Report

Sai Vinay Bhoomireddy. Low-Re CFD validation, flight-log system identification and a closed design loop for small UAVs, with StratoGlide (tag `stratoglide-v1`, read-only) as the design input.

Every number in this report is rendered from `results/numbers.json` or from a generated results file by `scripts/build_docs.py`; `scripts/check_numbers.py` (check V19) fails on any drift. Data classes are stated per section: public experimental, CFD, low-order, synthetic. No aircraft was built or flown and no wind tunnel was used for this project.

## 1. Verifier status

Command: `.venv/bin/python verify.py`. [I] = integrity check, [R] = result bar pre-registered in `PREREG.md` before any result was committed.

| Check | Class | Result | Measured |
|---|---|---|---|
| V01 | I | PASS | rho11 0.3648 (+0.25%), rho20 0.0889 (+0.01%), nu18 1.169e-04 (+0.96%) |
| V02 | I | PASS | data files without SOURCES row: 0; TM-4062 unresolved diffs 0; canonical logs 210 (simulator outputs must say synthetic), mislabelled 0 |
| V03 | I | PASS | max y+ 0.653 over 14 reported runs; per-direction ratios 1.419, 1.417, 1.419, 1.417; runs without mesh SHA 0 |
| V04 | I | PASS | 6 runs in error statistics, non-converged among them 0; unsteady runs 5 (all in Trust Map: True) |
| V05 | R | FAIL (logged) | SA: GCI Cl 0.32% (p 1.97), Cd 1.00% (p 1.97); SST_LM: Cl monotone (unsteady inputs), Cd oscillatory (bars Cl < 2%, Cd < 5%, p in [0.5, 3]) |
| V06 | R | FAIL (logged) | SST-LM Re 200k a4 (state unsteady, in error stats False): Cl +0.3%, Cd +4.4%, RMS dCp upper 0.078 (bars 5%, 20%, 0.05) |
| V07 | R | PASS | SA Cd +25.1% vs TM-4062 (bar > +25%) |
| V08 | R | FAIL (logged) | a0: no bubble/reattachment in CFD; a4 (NOT in error stats): x_s +0.013c, x_r +0.021c (bars 0.05c, 0.08c) |
| V09 | I | PASS | Trust Map covers 12/12 E387 runs; LIMITS.md quantified failure regions 8 |
| V10 | I | PASS | DISAGREEMENT.md rows 21, empty cells 0, Re covered ['100', '60', '80'] |
| V11 | R | FAIL (logged) | Cd MAPE corrected 9.5% vs raw XFOIL 12.6% (ratio 0.76, bar <= 0.75) vs NeuralFoil 10.7%; leakage assert [I] PASS |
| V12 | R | PASS | M2: beta 39.32 (+0.01 deg), p2/p1 1.7065 (-0.03%); M3: beta 27.40 (+0.02 deg), p2/p1 2.0544 (-0.03%) (bars 0.5 deg, 1%) |
| V13 | R | PASS | M2 a2 fine: Cl 0.0803 (-0.85%), Cd 0.00847 (-1.87%), monotone Cl True, Cd True |
| V14 | I | PASS | synthetic Parquet round-trip exact True; PX4 ULogs to full schema 1; ArduPilot .bin to full schema 3 |
| V15 | R | PASS | CD0 1.5% (<5), k 1.9% (<10), CLa 2.0% (<5), Cma 2.0% (<10), Cmde 2.1% (<10), Cmq 3.4% (<25) median /err/, 20 seeds |
| V16 | R | FAIL (logged) | CD0 25%, CLa 60%, Cma 60% of seeds inside 95% CI (bar 80-100%) |
| V17 | I | PASS | ALFA flights used 8, segments 21; estimates without CI/residual 0; AVL-vs-flight table present |
| V18 | R | FAIL (logged) | 6 blind seeds: CD0 delta err max 2.89 pp, SM delta err max 0.17% MAC, sink err cut 4%, trim err cut 51%; sealed truth opened by engine [I]: False |
| V19 | I | PASS | missing outward files []; files with U+2014 []; numbers in resume/report not equal to numbers.json: 0 |
| V20 | I | PASS | secrets 0; attribution in git log False; weapon-term hits []; PREREG committed before first results commit True |

INTEGRITY 10/10, RESULTS 4/10.

## 2. Module A: E387 at Re 200k vs NASA TM-4062 (CFD vs public experimental)

SU2 v8.6.0 incompressible RANS, SA (fully turbulent) and SST with Langtry-Menter transition (SST-LM, Tu as measured in the LTPT), own C-mesh with first-cell y+ max 0.65 over 11 reported runs, three systematically refined grids at alpha 4. A run counts in error statistics only when it meets the residual or Cauchy criterion; 4 runs ended unsteady at the iteration cap and are reported with their oscillation half-amplitude.

| alpha | Model | State | Cl | Cl exp | Cl err | Cd | Cd exp | Cd err | Cd half-amplitude | RMS dCp upper |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | SA | converged | 0.3615 | 0.3510 | +3.0% | 0.01380 | 0.01055 | +30.8% | steady | 0.051 |
| 0 | SST_LM | unsteady | 0.3697 | 0.3510 | +5.3% | 0.01090 | 0.01055 | +3.3% | 0.00078 | 0.041 |
| 4 | SA | converged | 0.7851 | 0.7860 | -0.1% | 0.01665 | 0.01330 | +25.1% | steady | 0.096 |
| 4 | SST_LM | unsteady | 0.7887 | 0.7860 | +0.3% | 0.01389 | 0.01330 | +4.4% | 0.00047 | 0.078 |
| 4 | SST | converged | 0.7721 | 0.7860 | -1.8% | 0.01677 | 0.01330 | +26.0% | steady | 0.091 |
| 8 | SA | converged | 1.1679 | 1.1794 | -1.0% | 0.02403 | 0.01744 | +37.8% | steady | 0.212 |
| 8 | SST_LM | unsteady | 1.1579 | 1.1794 | -1.8% | 0.02302 | 0.01744 | +32.0% | 0.00224 | 0.247 |

**Laminar separation bubble (SST-LM vs TM-4062 oil flow, chord fractions).**

| alpha | Oil flow x_s | Oil flow x_r | SST-LM x_s | SST-LM x_r | Error x_s | Error x_r | Note |
|---|---|---|---|---|---|---|---|
| 0 | 0.48 | 0.74 | n/a | n/a | n/a | n/a | run not in error statistics |
| 4 | 0.40 | 0.62 | 0.413 | 0.641 | +0.013 | +0.021 | run not in error statistics |
| 8 | n/a | n/a | 0.004 | 0.010 | n/a | n/a | NT at 0.32 |

**Grid convergence (Celik et al. 2008).**

| Model | Quantity | States (fine/medium/coarse) | Fine / medium / coarse | Type | p | GCI fine |
|---|---|---|---|---|---|---|
| SA | CL | converged/converged/converged | 0.78701 / 0.78505 / 0.78114 | monotone | 1.97 | 0.32% |
| SA | CD | converged/converged/converged | 0.01651 / 0.01665 / 0.01691 | monotone | 1.97 | 1.00% |
| SST_LM | CL | converged/unsteady/unsteady | 0.79185 / 0.78873 / 0.78636 | monotone (unsteady inputs) | 0.81 | 1.51% |
| SST_LM | CD | converged/unsteady/unsteady | 0.01404 / 0.01389 / 0.01437 | oscillatory | n/a | n/a (oscillatory) |

**Reading.** The fully turbulent SA model over-predicts drag at every alpha because it has no laminar run and no bubble; this is the pre-registered physics story (V07). SST-LM reproduces the alpha 4 bubble location within +0.013c (separation) and +0.021c (reattachment) and the time-mean Cd within +4.4%, but on the coarse and medium grids it never reached a steady state, so by the pre-registered rule those runs are outside error statistics and V05, V06 and V08 fail. A URANS treatment of the bubble is deferred (GAPS.md).

**Quantified failure regions (LIMITS.md).**

1. SA (fully turbulent) at alpha 0 deg: Cd error +30.8% (friction share 80% of Cd), because SA has no laminar run and no bubble.
2. SA (fully turbulent) at alpha 4 deg: Cd error +25.1% (friction share 66% of Cd), because SA has no laminar run and no bubble.
3. SA (fully turbulent) at alpha 8 deg: Cd error +37.8% (friction share 38% of Cd), because SA has no laminar run and no bubble.
4. SST-LM at alpha 0 deg: state unsteady, Cd oscillation half-amplitude 7.2% of mean over the last 2000 iterations; Cl error +5.3%, Cd error +3.3%; no separation found in CFD.
5. SST-LM at alpha 4 deg: state unsteady, Cd oscillation half-amplitude 3.4% of mean over the last 2000 iterations; Cl error +0.3%, Cd error +4.4%; bubble x_s 0.413c vs oil flow 0.4, x_r 0.641c vs 0.62.
6. SST-LM at alpha 8 deg: state unsteady, Cd oscillation half-amplitude 9.7% of mean over the last 2000 iterations; Cl error -1.8%, Cd error +32.0%; no separation found in CFD.
7. Grid dependence SA alpha 4 deg: Cl monotone GCI 0.32%, Cd monotone GCI 1.00% (p 1.97).
8. Grid dependence SST_LM alpha 4 deg: Cl monotone (unsteady inputs) GCI 1.51%, Cd oscillatory.

## 3. Module B: SD7037 at stratospheric Reynolds numbers (CFD, low-order vs public experimental)

# DISAGREEMENT: SD7037 at stratospheric Reynolds numbers (generated by scripts/run_strato.py)

Airfoil: SD7037 (B), as-tested UIUC coordinates (StratoGlide's selected airfoil family). Data: UIUC LSAT (public experimental), SU2 (CFD), XFOIL/NeuralFoil (low-order). Cd error = mean |Cd_tool - Cd_exp| / Cd_exp over CL 0.5, 0.6, 0.7, 0.8 (SU2: at the CL of each run).

Low-order tools run at the experiment's Re (60.5k, 101.4k); SU2 runs at 60k and 100k. XFOIL points with friction drag below the one-sided laminar flat-plate value 1.328/sqrt(Re) are dropped as non-physical: Re 60k Ncrit 9: 0, Ncrit 11: 0, Re 80k Ncrit 9: 0, Ncrit 11: 0, Re 100k Ncrit 9: 0, Ncrit 11: 0.

| Re | Tool | Cd error vs experiment | Cl_max | Bubble / transition |
|---|---|---|---|---|
| Re 60k | SU2 SST-LM | n/a (run not finished) | n/a (alpha sweep not run) | n/a |
| Re 60k | SU2 SA | n/a (run not finished) | n/a (alpha sweep not run) | n/a |
| Re 60k (run at 60.5k) | XFOIL Ncrit 9 | 10.8% | 1.255 | x_tr upper 0.72c at CL 0.7 |
| Re 60k (run at 60.5k) | XFOIL Ncrit 11 | 40.0% | 1.211 | x_tr upper 0.70c at CL 0.7 |
| Re 60k (run at 60.5k) | NeuralFoil xxlarge | 7.7% | 1.314 | n/a (no transition output) |
| Re 60k (run at 60.5k) | XFOIL corrected (Module C) | 4.2% | n/a (correction is for Cd only) | n/a (no transition output) |
| Re 60k | UIUC experiment | reference | 1.221 | n/a (no flow visualisation in LSAT vol. 3 for SD7037) |
| Re 80k | SU2 SST-LM | n/a (run not in compute budget) | n/a (alpha sweep not run) | n/a |
| Re 80k | SU2 SA | n/a (run not in compute budget) | n/a (alpha sweep not run) | n/a |
| Re 80k | XFOIL Ncrit 9 | n/a (no experiment) | 1.282 | x_tr upper 0.74c at CL 0.7 |
| Re 80k | XFOIL Ncrit 11 | n/a (no experiment) | 1.270 | x_tr upper 0.75c at CL 0.7 |
| Re 80k | NeuralFoil xxlarge | n/a (no experiment) | 1.335 | n/a (no transition output) |
| Re 80k | XFOIL corrected (Module C) | n/a (no experiment) | n/a (correction is for Cd only) | n/a (no transition output) |
| Re 80k | UIUC experiment | n/a (no experiment at this Re) | n/a (no experiment) | n/a (no flow visualisation in LSAT vol. 3 for SD7037) |
| Re 100k | SU2 SST-LM | 39.2% (at CL 0.68; unsteady) | n/a (alpha sweep not run) | x_s 0.28c, x_r none at CL 0.68 |
| Re 100k | SU2 SA | 7.6% (at CL 0.53, 0.74; converged) | n/a (alpha sweep not run) | no bubble in CFD |
| Re 100k (run at 101.4k) | XFOIL Ncrit 9 | 8.8% | 1.285 | x_tr upper 0.75c at CL 0.7 |
| Re 100k (run at 101.4k) | XFOIL Ncrit 11 | n/a (XFOIL did not converge) | n/a (XFOIL halted: 4 non-converged points from alpha 0) | n/a (no transition output) |
| Re 100k (run at 101.4k) | NeuralFoil xxlarge | 11.0% | 1.342 | n/a (no transition output) |
| Re 100k (run at 101.4k) | XFOIL corrected (Module C) | 1.1% | n/a (correction is for Cd only) | n/a (no transition output) |
| Re 100k | UIUC experiment | reference | 1.245 | n/a (no flow visualisation in LSAT vol. 3 for SD7037) |

## Why the tools disagree (each item is tested against the table, not asserted)

1. XFOIL Ncrit 9 vs 11: Mack's relation Ncrit = -8.43 - 2.4 ln(Tu) maps Ncrit 9 to Tu 0.070% and Ncrit 11 to Tu 0.030%; the spread between the two rows is the sensitivity to the (unknown) disturbance level, which is the argument for quieter stratospheric air, not a measurement of it.
2. XFOIL's e^N envelope with viscous-inviscid coupling is calibrated on natural Tollmien-Schlichting transition in quiet tunnels, which is the mechanism that sets the bubble at this Re; its Cd error here is the size of the residual that Module C learns.
3. SU2 SST-LM uses Langtry-Menter correlations calibrated mostly on bypass transition; the result depends on the Tu that survives to the leading edge (freestream decay is logged per run).
4. SU2 SA is fully turbulent: no laminar run and no bubble, so it over-predicts friction drag by construction at Re 60k to 100k.
5. NeuralFoil is trained on XFOIL output, so it inherits XFOIL's bias; its authors flag Re near 80k as a hard regime (analysis_confidence median reported per Re).
6. Corrected XFOIL applies the Module C model trained on every non-test family; SD7037 never entered training, so its row is a genuine out-of-sample test.
7. UIUC LSAT vol. 3 has drag polars at Re 60.5k and 101.4k for SD7037 (B); there is no Re 80k experiment, so the 80k rows report tool spread only.

NeuralFoil analysis_confidence (median over the alpha sweep): Re 60k 0.98, Re 80k 0.98, Re 100k 0.98.

## 4. Module C: multi-fidelity drag correction (low-order vs public experimental)

A correction to XFOIL (Ncrit 9) drag, log(Cd_exp / Cd_xfoil) regressed on XFOIL outputs and section geometry, chosen by inner GroupKFold and scored leave-one-family-out over 18 UIUC airfoil families (1356 identical airfoil, Re, Cl points). E387 and SD7037 families never enter training.

| Tool | Cd MAPE |
|---|---|
| corrected XFOIL | 9.54% |
| raw XFOIL Ncrit 9 | 12.59% |
| NeuralFoil 0.3.3 xxlarge | 10.71% |
| held-out test airfoils: corrected / raw XFOIL | 9.76% / 11.09% |
| E387 Re 100k signed Cd error: raw XFOIL / corrected | -10.4% / -6.8% |

The reduction against raw XFOIL is 24.3%; the pre-registered bar asks for a ratio of at most 0.75, so V11 fails by a small margin while the corrected model beats NeuralFoil on the same points. StratoGlide's XFOIL validation failed its own 15% bar at Re 100k on E387 with Ncrit 5 (its pre-registered rule); FlightLoop runs Ncrit 9 and 11 at matched UIUC points, so the two numbers are not the same metric.

## 5. Module D: supersonic wedge and diamond (CFD vs textbook theory)

Verification of the Euler solver setup only. SU2 Euler (Roe, Venkatakrishnan limiter), three nested grids, oblique-shock and shock-expansion theory from `flightloop/compressible/theory.py`.

| Case | Quantity | Coarse | Medium | Fine | Theory | Fine error | States (c/m/f) |
|---|---|---|---|---|---|---|---|
| wedge M2, 10 deg | shock angle (deg) | 39.324 | 39.321 | 39.318 | 39.314 | +0.004 deg | converged-cauchy/converged-residual/converged-cauchy |
| wedge M2, 10 deg | p2/p1 | 1.7064 | 1.7065 | 1.7065 | 1.7066 | -0.002% | |
| wedge M3, 10 deg | shock angle (deg) | 27.425 | 27.406 | 27.395 | 27.383 | +0.012 deg | converged-residual/converged-residual/converged-residual |
| wedge M3, 10 deg | p2/p1 | 2.0540 | 2.0542 | 2.0544 | 2.0545 | -0.005% | |
| diamond M2_a0 | Cd | 0.005278 | 0.005520 | 0.005640 | 0.005780 | -2.41% | converged-residual/converged-residual/converged-residual |
| diamond M2_a2 | Cl | 0.07819 | 0.07959 | 0.08032 | 0.08100 | -0.84% | |
| diamond M2_a2 | Cd | 0.008043 | 0.008325 | 0.008469 | 0.008629 | -1.85% | converged-residual/converged-cauchy/converged-cauchy |
| diamond M2_a4 | Cl | 0.15665 | 0.15940 | 0.16089 | 0.16225 | -0.84% | |
| diamond M2_a4 | Cd | 0.016307 | 0.016744 | 0.016974 | 0.017207 | -1.35% | converged-cauchy/converged-cauchy/converged-cauchy |
| diamond M3_a0 | Cd | 0.003400 | 0.003466 | 0.003500 | 0.003546 | -1.29% | converged-cauchy/converged-residual/converged-cauchy |
| diamond M3_a2 | Cl | 0.04709 | 0.04846 | 0.04915 | 0.04986 | -1.44% | |
| diamond M3_a2 | Cd | 0.005084 | 0.005188 | 0.005242 | 0.005308 | -1.24% | converged-residual/converged-residual/converged-residual |
| diamond M3_a4 | Cl | 0.09449 | 0.09722 | 0.09860 | 0.10006 | -1.45% | |
| diamond M3_a4 | Cd | 0.010114 | 0.010360 | 0.010486 | 0.010626 | -1.33% | converged-residual/converged-residual/converged-cauchy |

## 6. Module E: system identification

**Synthetic (labelled synthetic).** 20 seeds of the StratoGlide design model in FlightLoop's 6-DOF with MIL-F-8785C light turbulence and sensor noise; glide force method for CD0 and k, equation error for the longitudinal derivatives, output error reported alongside.

| Parameter | Median abs. error, equation error / force method | Median abs. error, output error |
|---|---|---|
| CD0 | 1.45% | n/a (glide method only) |
| k | 1.88% | n/a (glide method only) |
| CLa | 1.96% | 3.16% |
| Cma | 2.00% | 4.87% |
| Cmde | 2.07% | 9.15% |
| Cmq | 3.42% | 17.08% |

95% CI coverage over seeds: CD0 25%, CLa 60%, Cma 60% (V16 fails; the per-flight pitot scale error is not in the regression standard errors). Findings: the energy-method CD0 is biased by vertical gusts (16.4% with gusts, 3.0% calm); alpha estimated without a vane biases CLa by -98%; at 18 km the glide CD0 median error is 0.50%.

**CMU ALFA (public experimental).** Energy-method glide polar on post-engine-failure windows of 8 flights (21 windows), bootstrap over flights.

| Quantity | Flight | 95% CI | AVL (approximate geometry) |
|---|---|---|---|
| CD0 (airframe + stopped prop) | 0.0212 | -0.0222 to 0.0478 | n/a (inviscid) |
| k | 0.159 | 0.112 to 0.213 | 0.0574 |
| span efficiency e | 0.353 | from k | 0.9795 |

The flight k includes every lift-dependent drag term (profile drag growth, trim, prop wake), not only induced drag, and mass and wing area are retailer figures; the CD0 interval spans zero. Parsers: 1 public PX4 ULog and 3 ArduPilot .bin files reach the full canonical schema.

## 7. Module F: next-iteration engine and correlation report (synthetic)

`stubs/inject.py` draws a discrepancy (CD0 growth, aft CG, elevator effectiveness), applies it to the design model and seals it; the engine sees only the flight-1 logs and the design model (its file access is audited), identifies, applies first-order design rules and predicts flight 2.

| Metric over 6 blind seeds | Value |
|---|---|
| CD0 delta error, max | 2.89 pp |
| static margin delta error, max | 0.17% MAC |
| next-flight sink error reduction, median | 3.9% |
| next-flight trim error reduction, median | 51.1% |
| secondary (post hoc): air-relative sink, all glides / inside envelope | 55.5% / 80.3% |
| secondary (post hoc): trim inside envelope | 97.0% |

The recovery parts of V18 pass; the next-flight prediction part fails on the pre-registered metric (HONEST_LOG 14 and the V18 entry give the cause: inertial sink carries the air-mass vertical velocity, and the two slowest glides fly beyond the linear envelope).

**Correlation against StratoGlide's test cards.** 25 cards, 12 measurable from the synthetic logs, 5 inside StratoGlide's band.

| Card | Quantity | Unit | StratoGlide prediction | Band | FlightLoop design recompute | Synthetic measured | Status |
|---|---|---|---|---|---|---|---|
| FT-1 | stall EAS | m/s | 11 | 10.45 to 11.55 | n/a | n/a | n/a (no stall test in the synthetic schedule; the identification envelope stops at alpha 9.6 deg) |
| FT-1 | V_bestLD EAS (SL) | m/s | 14.13 | 12.72 to 15.54 | 12.04 | 11.31 | outside band |
| FT-1 | L/D max (SL) | - | 11.9 | 10.71 to 13.09 | 9.875 | 8.607 | outside band |
| FT-1 | min sink rate (SL) | m/s | 1.08 | 0.9722 to 1.188 | 1.07 | 1.153 | inside band |
| FT-1 | trim elevator at loiter CL (SL) | deg | -1.705 | -2.132 to -1.279 | -1.617 | -0.5763 | outside band |
| FT-1 | static margin | %MAC | 11.49 | 8.491 to 14.49 | 11.49 | 7.49 | outside band |
| FT-2 | electrical power at SL cruise (loiter CL) | W | 90.43 | 81.38 to 99.47 | n/a | n/a | n/a (blind flights are unpowered) |
| FT-2 | climb rate at full power (SL) | m/s | 11.03 | 9.931 to 12.14 | n/a | n/a | n/a (blind flights are unpowered) |
| FT-2 | RPM at SL cruise | rpm | 1983 | 1785 to 2181 | n/a | n/a | n/a (blind flights are unpowered) |
| FT-3 | trim gradient d(delta_e)/dCL | deg/CL | -5.065 | -6.331 to -3.799 | -5.065 | -3.744 | outside band |
| FT-3 | neutral point x (from wing root LE) | m | 0.1584 | 0.1505 to 0.1664 | 0.1584 | 0.1474 | outside band |
| FT-4 | short-period natural frequency (SL) | rad/s | 7.926 | 6.341 to 9.511 | 8.394 | 7.305 | inside band |
| FT-4 | short-period damping ratio (SL) | - | 0.6763 | 0.5411 to 0.8116 | 0.6618 | 0.7244 | inside band |
| FT-4 | Cm_alpha | 1/rad | -0.5462 | -0.6828 to -0.4097 | -0.5462 | -0.3493 | outside band |
| FT-4 | Cm_q | 1/rad | -8.921 | -11.15 to -6.691 | -8.921 | -8.299 | inside band |
| FT-4 | Cm_delta_e | 1/deg | -0.02269 | -0.02836 to -0.01702 | -0.02269 | -0.01962 | inside band |
| FT-4 | Cl_p | 1/rad | -0.4384 | -0.548 to -0.3288 | n/a | n/a | n/a (no roll excitation flown) |
| FT-5 | altitude lost in pull-out | m | 40.45 | 27.47 to 69.28 | n/a | n/a | n/a (pull-out cards FT-5/FT-6 not flown) |
| FT-5 | peak load factor | g | 2.003 | 1.777 to 2.302 | n/a | n/a | n/a (pull-out cards FT-5/FT-6 not flown) |
| FT-5 | peak EAS | m/s | 21.39 | 18.49 to 23.86 | n/a | n/a | n/a (pull-out cards FT-5/FT-6 not flown) |
| FT-5 | time to recover | s | 11.27 | 7.307 to 19.66 | n/a | n/a | n/a (pull-out cards FT-5/FT-6 not flown) |
| FT-6 | altitude lost in pull-out | m | 90.17 | 59.72 to 170.8 | n/a | n/a | n/a (pull-out cards FT-5/FT-6 not flown) |
| FT-6 | peak load factor | g | 2.265 | 2.101 to 2.521 | n/a | n/a | n/a (pull-out cards FT-5/FT-6 not flown) |
| FT-6 | peak EAS | m/s | 21.03 | 18.14 to 23.83 | n/a | n/a | n/a (pull-out cards FT-5/FT-6 not flown) |
| FT-6 | time to recover | s | 15.81 | 9.957 to 30.62 | n/a | n/a | n/a (pull-out cards FT-5/FT-6 not flown) |

## 8. Failed bars and deviations

- **V05 FAIL (GCI, SST-LM).** E387 Re 200k alpha 4: SA passes the bar on its own (Cl GCI 0.32%, Cd GCI 1.00%, observed p 1.97 for both, monotone). SST-LM coarse and medium ended `unsteady` at the 30000-iteration cap (Cd half-amplitude over the last 2000 iterations 0.00047 on medium) while fine converged (Cl 0.79185, Cd 0.014036 at 16328 iterations); with unsteady inputs Cl is monotone with p 0.81 (GCI 1.51%; p is inside [0.5, 3] but the inputs are not converged) and Cd is oscillatory, so no GCI is defined. Not fixed by changing iteration caps, Tu or model constants; URANS on the medium grid is deferred (GAPS.md).
- **V06 FAIL (SST-LM vs TM-4062 at alpha 4).** Medium grid state `unsteady`, so the run is outside error statistics by the pre-registered rule. Mean over the last 2000 iterations: Cl +0.3%, Cd +4.4% (within the 5% and 20% bars), RMS dCp upper 0.078 (bar 0.05). Largest residuals at x/c 0.60 (+0.072, bubble plateau end) and 0.70 (-0.042, reattachment recovery).
- **V07 PASS, barely.** SA (fully turbulent) Cd +25.1% vs TM-4062 0.0133 at alpha 4 (bar > +25%); alpha 0 +30.8%, alpha 8 +37.8%.
- **V08 FAIL (bubble).** Alpha 4: SST-LM separation x_s 0.413c vs oil flow 0.40c (+0.013c), reattachment x_r 0.641c vs 0.62c (+0.021c), inside the 0.05c/0.08c bars, but the run is `unsteady` so it does not count. Alpha 0: no separation found in the CFD (oil flow x_s 0.48c, x_r 0.74c). Alpha 8: TM-4062 reports natural transition at 0.32c (no bubble number to score); CFD shows a short leading-edge separation at 0.004c to 0.010c.

- **V16 FAIL (CI calibration).** Run `sysid-synth20-20261118-120733-030dcbe`, 20 seeds: the pre-registered 95% CIs cover truth in CD0 25%, CLa 60%, Cma 60% of seeds (bar 80 to 100%). Cause: the window-regression SE (CD0) and the Morelli-Klein SE (CLa, Cma) contain only random error; every flight also carries a pitot scale error N(1, 0.01) that scales all coefficients by about 2x the scale error, which is not in either SE. The bar is not changed. Secondary, post hoc: adding the 2 x 1% calibration term in quadrature gives CD0 100%, CLa 95%, Cma 85% (`coverage_with_calibration_pct`); the block-bootstrap CI gives CLa 80%, Cma 75%, Cmde 75%. V15 passes on the same run (median |error| CD0 1.45%, k 1.88%, CLa 1.96%, Cma 2.00%, Cmde 2.07%, Cmq 3.42%).
- **No-vane alpha is unusable in turbulence (finding, not a bar).** Alpha from EKF attitude and GPS velocity misses gust-induced alpha; equation-error CLa median bias -98%, Cma +97% on the same flights. The energy-method glide polar is biased by vertical gusts: CD0 16.4% and k 74.5% median error with gusts, 3.0% and 9.8% with vertical gusts off.
- **V11 FAIL by 0.007.** Run `mf-lofo-20261118-121208-030dcbe`, leave-one-family-out over 18 families, 1356 identical (airfoil, Re, Cl) points on 38 airfoils: corrected Cd MAPE 9.54% vs raw XFOIL Ncrit 9 12.59% (ratio 0.757, bar <= 0.75) and NeuralFoil 0.3.3 xxlarge 10.71% (beats NeuralFoil, second clause passes). Constant-bias baseline 12.51%. Held-out test airfoils (E387, SD7037 and family members, model gp_nu2.5 chosen by inner GroupKFold): corrected 9.76% vs XFOIL 11.09% vs NeuralFoil 14.56%. On E387 at Re 100k this XFOIL setup is -10.4% in Cd (corrected -6.8%), not the 38.9% CD MAPE StratoGlide reported: StratoGlide ran XFOIL at Ncrit 5 (its pre-registered rule) over its own CL sweep (stratoglide HONEST_LOG, P2 validation), FlightLoop runs Ncrit 9 and 11 at matched UIUC points (signed error at points, not MAPE over a sweep). The bar is not changed and the model family is not re-searched.

- **V18 FAIL (next-flight error reduction).** Run `loop-blind6-20261118-182042-0402430`, 6 blind seeds, engine never opened SEALED_truth.json. Recovery passes: max |CD0 delta error| 2.89 pp (bar 5), max |SM delta error| 0.17% MAC (bar 1). Prediction bars fail: median sink-rate error reduction 3.9%, trim-elevator 51.1% (bar 70% each). Cause: flight-2 inertial sink is a 15 s mean of vd that carries the air-mass vertical velocity, and the 13.5 and 14.5 m/s glides fly beyond the 9 deg linear envelope where no linear model predicts trim. Secondary, post hoc (HONEST_LOG 14): air-relative sink over all glides 55.5%, air-relative sink inside the envelope 80.3%, trim inside the envelope 97.0%. The bar and the primary metric are not changed.

All deviations from the plan, each with its reason and timing relative to results, are numbered in `HONEST_LOG.md`.

## 9. Deferred

# GAPS (deferred, with reason)

- 3D CFD of a full aircraft, LES/DNS, URANS bubble-bursting, optimizer: out of scope per plan section 0.
- Tripped, flapped and Gurney-flap UIUC configurations: excluded from phase 1 (plan P1.2).
- TM-4062 lower-surface Cp: not transcribed (V06 uses the upper surface only).
- OpenFOAM kOmegaSSTLM cross-check of E387: not run (plan order puts it after site/; compute went to the cut-line SST-LM grids, which ran to the 30,000-iteration cap).
- SD7037 SU2 at Re 60k (SA and SST-LM, alpha 2 and 4): deferred, last in plan order; the low-order tools and UIUC data at Re 60.5k are in DISAGREEMENT.md.
- SD7037 SU2 at Re 80k and alpha sweeps for Cl_max: not in the compute budget (n/a rows in DISAGREEMENT.md).
- URANS for the E387 SST-LM bubble: the coarse and medium grids end unsteady; a time-accurate run is needed to score V05, V06, V08 on them.
- E387 Tu sweep (0.1, 0.2%), Re 100k E387 CFD, SST and SA-BC at alpha 0 and 8: queued as `fullA_core` in scripts/cfd_queue.py, not run. SST and SA-BC at alpha 4 run as `fullA_min` if they finish before the final build.
- Cl correction in Module C (Cd only), lift-curve correction deferred.
- ArduPilot SITL .bin generation: replaced by 3 public ArduPlane .bin files from the CMU ALFA dataflash archive (V14).
- PX4 public logs: 1 of 3 downloaded ULogs carries airspeed; the other two parse but have no TAS (parse_report.json).
- V16 CI calibration: the per-flight pitot scale error is not in the regression SE; a calibration-aware CI is reported as secondary only.
- Wrap-up stops: E387 SA-BC alpha 4 (stopped at about 8,200 of 15,000 iterations) and SD7037 Re 100k SST-LM alpha 2 (about 20,400 of 30,000) were stopped by hand at the final wrap-up; both are `incomplete`, listed in the registry, and excluded from every statistic. E387 SST (fully turbulent) alpha 4 converged and is in the registry.

## 10. Compute and reproduction

43.8 core-hours measured from the run folders (35.8 in SU2 over 42 runs including superseded ones), 0 cloud hours; MacBook Air M4, 16 GB. Environment pins: `requirements.lock`, `ENV.md`.

```bash
.venv/bin/python verify.py
.venv/bin/python scripts/build_site.py && .venv/bin/python scripts/build_docs.py
.venv/bin/python scripts/check_numbers.py
```
