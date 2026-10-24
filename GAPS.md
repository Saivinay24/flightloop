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
