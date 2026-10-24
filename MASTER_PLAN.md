# FlightLoop MASTER PLAN: a multi-fidelity aero correlation engine (CFD <-> XFOIL/AVL <-> wind tunnel <-> flight logs -> next-iteration design deltas)

> Paste the /goal block at the bottom into a fresh Claude Code session opened in `~/Desktop/flightloop/`. The session that receives it does the work. The planning session that wrote this file is NOT executing it.
> Target role: Aerodynamics Engineer, Kalam Labs (`~/Desktop/kalam-labs-application/JD.txt`). The JD line this project answers word for word: "Run CFD (ANSYS Fluent, OpenFOAM, SU2 or similar) and lower-order tools (XFLR5, AVL), and correlate the results with wind-tunnel and flight data ... use flight logs to improve the next iteration." Requirement: "Hands-on CFD, with a working understanding of meshing, turbulence models and the limits of the results."
> Written 2026-11-18. Theory numbers in this file (oblique shock, shock-expansion, y+, Mack relation) were computed with a script in this planning session, not quoted from memory. Experimental numbers are NOT given here on purpose: the agent extracts them from the primary sources and records provenance.

---

## 0. GOAL (completion condition)

FlightLoop is done when `python3 verify.py` prints **INTEGRITY 12/12** and every RESULT check is either PASS or a logged honest FAIL (actual number in `HONEST_LOG.md`, `results/numbers.json` updated, resume contract rewritten to the honest number), demonstrated with actual command output in the transcript, not described from memory; and the following exist: SU2 E387 validation vs NASA TM-4062 with 3-grid GCI and laminar separation bubble comparison (A), the stratospheric-Re disagreement study (B), the multi-fidelity correction model with leave-one-airfoil-out validation (C), the supersonic wedge and diamond verification (D), the flight-log sys-ID pipeline on synthetic and real public logs (E), the next-iteration engine demonstrated on blind synthetic cases (F), a dashboard plus report PDF (G), and milestone tag `flightloop-v1`. Stop when that is true, or after you have surfaced a blocker you cannot resolve without the user, or after 60 agent turns, in which case `BLOCKERS.md` exists explaining precisely what is missing and why.

**Do NOT push to any remote, create a GitHub repo, post anywhere, or deploy the dashboard in this run.** Publishing is a later phase that Vinay triggers himself (clean single-commit push, secrets grep first).
**Do NOT claim ownership of any wind tunnel, aircraft or flight hardware.** All experimental data is public (NASA, UIUC, CMU) or synthetic, and every file, plot and table says which.
**Do NOT fabricate or "fill in" any number.** A missing number is written as `n/a (reason)`.
**Do NOT write any weapon, missile, intake, seeker, guidance, targeting or warhead content.** Module D is a textbook wedge and diamond airfoil verification and nothing more.
**Do NOT build 3D CFD of a full aircraft, LES/DNS, URANS bubble-bursting studies, or an optimizer.** These are listed in GAPS.md as deferred.

---

## 1. Context you need (read in this order)

1. This file. It is the source of truth. `PROGRESS.md` in the repo is the source of truth for state.
2. `~/Desktop/kalam-labs-application/JD.txt` (the role).
3. `~/Desktop/kalam-labs-application/research_notes/Kalam Labs aerodynamics role/aero_domain_problems.md` (ISA/Re tables, validation data sources, SU2/OpenFOAM transition model notes, sys-ID literature, NeuralFoil caveats). Treat items marked "[domain knowledge, unverified here]" as hypotheses to check.
4. `~/Desktop/vinay-clone/profile/01-thinking-style.md` and `03-technical-dna.md` (how Vinay builds: plan -> phases -> verifier, run IDs, checkpoint/resume, commit milestones).
5. If it exists: `~/Desktop/stratoglide/` (sister project). Read its `PROGRESS.md` and `exports/`. If it does not exist or exports are missing, use the stubs in section 4.3. Never block on StratoGlide.

Environment facts:
- Machine: MacBook Air M4, 16 GB, macOS (arm64). 4 performance + 6 efficiency cores. Use `mpirun -np 4` for SU2 on the M4.
- Optional: Azure Linux VM (Ubuntu 22.04) with many cores, via `az` CLI on Vinay's education account. Deallocate when idle (`az vm deallocate -g <rg> -n <vm>`). Never leave a VM running at end of run.
- Python 3.11 in a `uv` venv at `~/Desktop/flightloop/.venv`.
- Git identity: Sai Vinay Bhoomireddy (personal). **No Co-Authored-By, no Claude/Anthropic attribution in any commit, ever.**

---

## 2. Hardcoded config

```yaml
# config/flightloop.yaml  (all code reads this; nothing below is hardcoded in modules)
project_root: ~/Desktop/flightloop
run_id_format: "{module}-{case}-{YYYYmmdd-HHMMSS}-{shortsha}"   # every run gets isolated runs/<run_id>/
validation:
  e387:
    re: [100000, 200000]          # 60000 is stretch, see A.7
    alpha_deg: [-2, 0, 2, 4, 6, 8, 10, 12]
    alpha_gci_deg: 4
    models: [SA, SST, SA_BC, SST_LM]
    tu_percent_sweep: [0.05, 0.1, 0.2]   # LM inlet Tu sensitivity
  strato:
    airfoil: stratoglide_or_sd7037
    re: [60000, 100000]           # UIUC data exists at these; 80000 is CFD/XFOIL/NeuralFoil only
mesh:
  first_cell_height_c: 2.0e-5     # y+ 0.06-0.19 at Re 60k-200k (computed, flat-plate Cf), margin for 3x local Cf at suction peak
  growth_normal: 1.10
  farfield_radius_c: 100
  wake_length_c: 20
  levels:                          # C-mesh nodes: Ns airfoil wrap, Nw wake per side, Nn normal
    coarse: {Ns: 200, Nw: 50,  Nn: 80}    # ~24k quads
    medium: {Ns: 283, Nw: 71,  Nn: 113}   # ~48k quads
    fine:   {Ns: 400, Nw: 100, Nn: 160}   # ~96k quads; refinement ratio r = sqrt(2) per direction
xfoil: {panels: 200, iter: 200, ncrit: [9, 11], aseq: [-4, 14, 0.5]}
neuralfoil: {model_size: xxlarge}
sysid:
  seeds: 20
  log_rate_hz: 50
  sensor_noise: {accel_ms2: 0.05, gyro_rads: 0.005, baro_m: 0.3, pitot_ms: 0.3, pitot_scale: 0.01}
  turbulence: dryden_milf8785c_light
compressible:
  wedge:   {mach: [2.0, 3.0], theta_deg: 10}
  diamond: {mach: [2.0, 3.0], t_over_c: 0.05, alpha_deg: [0, 2, 4]}
```

Rationale and fallbacks:
- **Incompressible RANS (`INC_RANS`) is the default for E387/SD7037** (TM-4062 ran at M 0.03 to 0.13). If P0 shows LM is not supported or unstable in SU2's incompressible solver, fall back to compressible `RANS` at M = 0.1 with low-Mach correction and log the switch in HONEST_LOG.md.
- **SA_BC** is SU2's `KIND_TRANS_MODEL= BC` (Bas-Cakmakcioglu algebraic intermittency on SA). **SST_LM** is `KIND_TURB_MODEL= SST` + `KIND_TRANS_MODEL= LM` (Langtry-Menter gamma-ReTheta). Exact option names are taken from the pinned SU2 version's `config_template.cfg`, not from this file.
- **Base every SU2 cfg on an in-repo SU2 TestCase** (incompressible turbulent NACA0012, transition T3A, inviscid wedge), diff only what is needed, and save the diff in `cfd/templates/DIFFS.md`.

---

## 3. Architecture

```
FlightLoop
├── A  CFD validation ........ gmsh C-mesh -> SU2 (SA | SST | SA-BC | SST-LM) -> E387 vs NASA TM-4062
│                               -> GCI (Celik 2008) -> LSB from Cf_x sign change -> Trust Map + LIMITS.md
│                               -> OpenFOAM kOmegaSSTLM cross-check (2 points)
├── B  Stratospheric Re ...... SD7037 (or StratoGlide foil) Re 60k/100k: SU2-LM vs XFOIL(Ncrit 9, 11) vs NeuralFoil vs UIUC
├── C  Multi-fidelity ........ UIUC LSAT polars x XFOIL -> learn delta(Cl, Cd) -> leave-one-airfoil-out -> vs raw XFOIL, NeuralFoil
├── D  Compressible verif .... SU2 Euler wedge + diamond, M 2 and 3 -> theta-beta-M and shock-expansion theory -> 3-grid convergence
├── E  Flight-log sys-ID ..... parsers (pymavlink .bin, pyulog .ulg, ALFA CSV) -> canonical Parquet
│                               -> glide/energy drag polar, equation-error, output-error, CIs
│                               -> synthetic (6-DOF, known truth) and real (CMU ALFA, PX4 public logs)
│                               -> correlate vs AVL (+ XFLR5 manual, optional)
├── F  Next-iteration engine . discrepancy table -> attribution -> design deltas (CG, tail volume, airfoil, prop)
│                               -> updated model -> blind synthetic re-fly -> error reduction
└── G  UI face ............... Streamlit dashboard (+ static HTML export) + report PDF, numbers from results/numbers.json only
```

### 3.1 Repo layout

```
~/Desktop/flightloop/
├── MASTER_PLAN.md (copy of this file)  PROGRESS.md  HONEST_LOG.md  GAPS.md  ENV.md  COMPUTE.md
├── PREREG.md            # success bars, committed BEFORE any result run (verify.py checks commit order)
├── LIMITS.md            # where each tool is wrong, quantified
├── config/flightloop.yaml
├── data/                # raw public data, read-only after download; SOURCES.md lists URL, page/figure, method, date
│   ├── tm4062/  uiuc_lsat/  airfoils/  alfa/  px4_logs/  ardupilot_logs/
├── inputs/stratoglide/  # symlink to StratoGlide exports, or stubs/ (section 4.3)
├── flightloop/          # python package
│   ├── atmos.py  geometry.py
│   ├── mesh/gmsh_cmesh.py  mesh/gmsh_wedge.py  mesh/gmsh_diamond.py
│   ├── cfd/su2_runner.py  cfd/postproc.py  cfd/gci.py  cfd/lsb.py  cfd/openfoam_case/
│   ├── lowfi/xfoil_runner.py  lowfi/neuralfoil_runner.py  lowfi/avl_runner.py
│   ├── mf/dataset.py  mf/models.py  mf/evaluate.py
│   ├── compressible/theory.py  compressible/postproc.py
│   ├── logs/schema.py  logs/ardupilot.py  logs/px4.py  logs/alfa.py
│   ├── sim/sixdof.py  sim/dryden.py  sim/sensors.py  sim/maneuvers.py
│   ├── sysid/glide.py  sysid/equation_error.py  sysid/output_error.py  sysid/ci.py
│   └── loop/correlate.py  loop/deltas.py  loop/report.py
├── runs/<run_id>/       # config snapshot, mesh hash, solver log, history.csv, surface.csv, status.json
├── results/numbers.json # SINGLE SOURCE OF TRUTH for every number in report, dashboard, resume
├── results/figures/
├── dashboard/app.py  dashboard/export_static.py
├── report/FlightLoop_Report.md -> FlightLoop_Report.pdf
├── scripts/run_all.sh   # one command, resumable (skips runs whose status.json says converged)
└── verify.py
```

### 3.2 Canonical flight-log schema (`flightloop/logs/schema.py`)

One Parquet per flight, resampled to 50 Hz, SI units, body axes FRD:

| column | unit | notes |
|---|---|---|
| `t` | s | monotone, from log start |
| `ax, ay, az` | m/s^2 | specific force, body |
| `p, q, r` | rad/s | body rates |
| `phi, theta, psi` | rad | EKF attitude |
| `vn, ve, vd` | m/s | EKF inertial velocity |
| `h_baro` | m | |
| `tas` | m/s | pitot, converted with ISA density from baro/temp |
| `wn, we` | m/s | EKF wind estimate if logged, else NaN |
| `de, da, dr` | rad or normalized | normalized [-1, 1] when PWM-to-angle calibration is unknown; flagged in metadata |
| `throttle` | 0-1 | |
| `volt, curr` | V, A | optional, for prop efficiency |
| `alpha, beta` | rad | vane if present, else derived (section E.4) and flagged |

Metadata JSON next to each Parquet: `source` (`synthetic` / `alfa` / `px4_public` / `ardupilot_sitl`), aircraft id, license, original file hash, derived-field flags.

---

## 4. Install and environment (P0 does this; record exact versions in ENV.md)

### 4.1 Commands

```bash
# Python
cd ~/Desktop/flightloop && uv venv --python 3.11 && source .venv/bin/activate
uv pip install numpy scipy pandas pyarrow matplotlib plotly scikit-learn streamlit \
  gmsh neuralfoil aerosandbox pymavlink pyulog pyyaml weasyprint markdown
uv pip freeze > requirements.lock

# Toolchain
brew install open-mpi gcc ninja        # gcc gives gfortran for XFOIL/AVL

# SU2 native build on M4 (pin the newest v8.x release tag that is at least 2 weeks old; record it)
git clone --depth 1 --branch <v8.x.y> https://github.com/su2code/SU2.git ~/opt/SU2-src
cd ~/opt/SU2-src
./meson.py build -Dwith-mpi=enabled -Denable-pywrapper=false --buildtype=release --prefix=$HOME/opt/su2
./ninja -C build install
export SU2_RUN=$HOME/opt/su2/bin PATH=$HOME/opt/su2/bin:$PATH   # add to scripts/env.sh
SU2_CFD -h | head -5
```

| Tool | Primary (M4) | Fallback 1 | Fallback 2 |
|---|---|---|---|
| SU2 v8.x | native meson build above | conda-forge `su2` if an osx-arm64 build exists | Azure VM: `sudo apt install build-essential openmpi-bin libopenmpi-dev python3` then the same meson commands |
| gmsh | `pip install gmsh` (Python API, writes `.su2` and `.msh` directly) | Homebrew `gmsh` | VM pip |
| XFOIL 6.99 | build from source with `gfortran`, plotting disabled or XQuartz installed | `docker run --platform linux/arm64 ubuntu:22.04` + `apt install xfoil` | VM `sudo apt install xfoil` |
| AVL 3.40 | build from source with `gfortran` (plot lib needs XQuartz; only text output is used) | Docker/VM Linux build | AeroSandbox `VortexLatticeMethod`, labelled "AeroSandbox VLM, not AVL" |
| OpenFOAM (v2406 or newer .com release) | Docker: check `docker manifest inspect opencfd/openfoam-default:<ver>` lists arm64 | VM: `curl -s https://dl.openfoam.com/add-debian-repo.sh \| sudo bash && sudo apt-get install openfoam2406-default` | skip OpenFOAM cross-check, log in HONEST_LOG |
| XFLR5 | GUI only, macOS build; manual by Vinay | none | omit from claims |
| NeuralFoil | `pip install neuralfoil` | `neuralfoil-torch` | none needed |

**P0 smoke tests (must pass before any module):** SU2 runs the incompressible NACA0012 TestCase to convergence; SU2 runs the T3A LM TestCase and the agent plots Cf vs Re_x against the ERCOFTAC data shipped with the tutorial (transition onset visibly in the right place); `xfoil` returns a viscous polar for NACA 2412 at Re 1e6; `avl` returns stability derivatives for its bundled sample `.avl`; gmsh writes a `.su2` that SU2 reads.

### 4.2 Compute budget

Per-run wall time model: `wall = cells x iters x c_ci x model_factor / (cores x eff)`. Assumed `c_ci` = 3 to 6 microseconds per cell-iteration-core (SU2 2D implicit RANS), `eff` = 0.8, `model_factor` 1.0 (SA, SST), 1.3 (LM). **These are estimates. P0 measures c_ci on the medium E387 mesh with 200 iterations and overwrites this table in COMPUTE.md before the full matrix launches.**

| Case | Cells | Iters to converge [est] | M4, 4 P-cores [est] | Notes |
|---|---|---|---|---|
| E387 SA/SST, medium | 48k | 3k-5k | 4-8 min | |
| E387 SST-LM, medium | 48k | 10k-20k | 15-30 min | bubble cases converge slowly; some never go steady |
| E387 SST-LM, fine | 96k | 15k-25k | 45-80 min | GCI only |
| E387 SST-LM, coarse | 24k | 8k-15k | 5-10 min | GCI only |
| Wedge / diamond Euler | 20k-80k | 1k-3k | 1-3 min | |
| OpenFOAM simpleFoam kOmegaSSTLM | 48k | 5k-15k | 20-40 min | |
| XFOIL polar | n/a | n/a | 1-3 s | 800 polars for module C under 1 h |
| 6-DOF synthetic flight, 300 s at 200 Hz | n/a | n/a | under 5 s | numpy RK4 |

| Matrix | Runs | M4 sequential [est] | Azure 32+ physical cores [est] |
|---|---|---|---|
| **Cut-line A**: Re 200k, alpha {0, 4, 8}, {SA, SST-LM}, medium + GCI at alpha 4 (coarse, fine) | 6 + 4 | ~3 h | ~30 min |
| Full A: 2 Re x 8 alpha x 4 models medium + GCI + Tu sweep | ~76 | ~16 h | ~2 h |
| A.7 Re 60k LM, 8 alpha | 8 | ~3.5 h | ~25 min |
| B: SD7037, 2 Re x 6 alpha x {SA, LM} | 24 | ~6 h | ~45 min |
| D: wedge 2 M x 3 grids + diamond 2 M x 3 alpha x 3 grids | 24 | ~1 h | ~10 min |
| OpenFOAM cross-check | 2-3 | ~1.5 h | ~20 min |
| **Total** | ~140 | **~30 h** | **~4 h** |

Decision rule (frugal by default): run the cut-line on the M4 overnight. If the measured `c_ci` puts full A + B above 12 h on the M4, spin up the Azure VM for the matrix only (a compute-optimized 64 vCPU size; check the current price in the portal and the education credit balance first), launch with `nohup`, rsync results back, deallocate. RAM is not a constraint: a 100k-cell 2D case needs well under 1 GB. **Per-unit time math goes in COMPUTE.md before every big launch** so kill/continue is a cheap decision.

### 4.3 StratoGlide interface contract (stubs when missing)

| File | From StratoGlide | Stub if missing |
|---|---|---|
| `inputs/stratoglide/airfoil.dat` | chosen airfoil, Selig format | `data/airfoils/sd7037.dat` (UIUC) |
| `inputs/stratoglide/aircraft.avl`, `aircraft.mass` | AVL model | `stubs/generic_glider.avl`: b = 2.0 m, S = 0.40 m^2 (AR 10), c = 0.20 m, SD7037 wing, conventional tail with V_h = 0.45, V_v = 0.03; mass 2.0 kg, Iyy 0.15 kg m^2 |
| `inputs/stratoglide/derivs.json` | CLa, Cma, Cmq, Cmde, CD0, k, Xnp, mass, inertia, S, c, b | generated by running AVL on the stub + XFOIL profile drag + Raymer component build-up for CD0 |
| `inputs/stratoglide/sixdof.py` with `simulate(controls, wind, seed) -> DataFrame` | their 6-DOF | `flightloop/sim/sixdof.py` (own) |

Stub truth for synthetic sys-ID (documented in `stubs/TRUTH.md`): CD0 = 0.025, e = 0.85, k = 1/(pi e AR) = 0.0374; longitudinal derivatives taken from the AVL run of the stub. FlightLoop must pass every check with stubs only.

---

## 5. Phases

Every phase: done-state, negative scope, steps, verifier checks it unlocks, milestone. Work one phase end to end before the next, except long solver runs, which go to `nohup` while another workstream advances.

### P0 Environment, benchmark, scaffold (2-3 h)

**Done-state:** all P0 smoke tests pass with output pasted in the transcript; ENV.md lists exact versions and build flags; COMPUTE.md has measured `c_ci`; repo scaffolded; PREREG.md written and committed (bars from section 6, copied verbatim); `verify.py` skeleton runs and reports which checks are not yet implemented.
**Negative scope:** no physics runs beyond smoke tests and the benchmark.
**Milestone:** `git tag flightloop-m0-env`.

### P1 Data acquisition and provenance (2-3 h, partly parallel with P2)

**Done-state:** every dataset downloaded or transcribed, with `data/SOURCES.md` rows (URL, document page/figure/table, method, date, license note).
**Negative scope:** no modelling.

Steps:
1. **NASA TM-4062** (NTRS 19890001471): download the PDF. Extract for E387 at Re 100k and 200k (and 60k): Cl, Cd, Cm vs alpha tables; Cp vs x/c at alpha 0, 4, 8; oil-flow laminar separation, transition and turbulent reattachment locations vs alpha. Tables: OCR (`tesseract`) plus a second independent transcription pass; `scripts/diff_transcriptions.py` must report zero differences, any disagreement resolved by reading the page image. Figures: `tools/digitize.py` with four-point axis calibration; record digitizing uncertainty (target +/-0.01 c). If oil-flow locations are not in TM-4062 itself, use the companion paper (NTRS 19900059314) and say so; last resort is Cp-plateau-inferred bubble location, labelled "inferred from Cp, not oil flow".
2. **UIUC LSAT** Vols 1-5 and SoarTech 8: download all ASCII polar files and the as-built model coordinates. Parse into one table `(airfoil, volume, Re, alpha, Cl, Cd, Cm, source_file)`. Exclude tripped, flapped and Gurney-flap configurations in phase 1 (park them in GAPS.md). Map each airfoil to a family id (for example all SD70xx variants of one shape share an id) for leakage-free splits.
3. **CMU ALFA** (doi 10.1184/R1/12707963): check the archive size first; download the processed data, then raw only if needed. Read the dataset README for the exact topic list; do not assume field names. Record aircraft facts (Carbon-Z T-28 airframe, ArduPlane on Pixhawk per the paper) with sources in `data/alfa/AIRCRAFT.md`; any geometry not found is "approximate, from manufacturer spec" with an uncertainty.
4. **PX4 public ULogs**: from the PX4 Flight Review public database, pick 3 to 5 fixed-wing logs with airspeed sensor data. Record log ids and URLs.
5. **ArduPilot .bin**: generate 2 logs with ArduPilot SITL (ArduPlane) on the Linux VM or Docker. Label as `ardupilot_sitl` (real format, simulated physics). If SITL is unavailable, find a public ArduPlane .bin with a clear license; else log the gap.

### P2 Module A: E387 CFD validation (critical path)

**Done-state:** GCI table, polar comparison, Cp comparison, LSB comparison, Trust Map and LIMITS.md for E387 at Re 200k (cut-line) and Re 100k (full), with SA, SST, SA-BC, SST-LM; OpenFOAM cross-check at 2 points.
**Negative scope:** no URANS, no 3D, no tuning of model constants, no adjusting alpha to match experiment (alpha offsets are reported, not applied).

Steps:
1. **Geometry:** E387 design coordinates from TM-4062 (cross-check against UIUC `e387.dat`), cubic spline by arc length, sharp trailing edge, unit chord.
2. **Mesh** (`mesh/gmsh_cmesh.py`): structured C-mesh with transfinite blocks and `Recombine`, farfield radius 100c, wake 20c, first cell 2e-5 c, normal growth 1.10. Streamwise clustering at LE and TE, and nearly uniform spacing (about 0.003c on medium) over 0.2c to 0.8c on the upper surface where the bubble lives, because LM needs streamwise resolution through the bubble. Three levels per config, refinement ratio sqrt(2) per direction. Export `.su2` with markers `airfoil`, `farfield`. Log cell count, max skewness, max aspect ratio, mesh SHA.
3. **y+ check:** after each run, read `Y_Plus` from the surface output; max over the airfoil must be < 1 (expected 0.1 to 0.4).
4. **Runs** (`cfd/su2_runner.py`): non-dimensional (chord 1, U 1, rho 1, mu = 1/Re). Implicit Euler pseudo-time, adaptive CFL, FGMRES + ILU, second-order MUSCL for flow. Convergence: RMS pressure residual drop >= 6 orders, or Cauchy criterion on Cl and Cd (relative change < 1e-6 over 1000 iterations). If neither is met after the iteration cap, record mean and amplitude of Cl over the last 2000 iterations and mark status `unsteady`. LM inlet: set freestream Tu and viscosity ratio; check the Tu that actually reaches x/c = -0.5 in the solution (freestream decay is a known LM pitfall) and log it; run the Tu sweep {0.05, 0.1, 0.2}% at alpha 0 and 4, Re 200k.
5. **GCI** (`cfd/gci.py`, Celik et al. 2008, J. Fluids Eng. 130, 078001): 3 grids, r = sqrt(N_fine/N_medium) in 2D, observed order p from the three solutions (fixed-point iteration if r differs between pairs), Richardson-extrapolated value, GCI_fine = 1.25 |e21| / (r^p - 1). Report apparent convergence type (monotone / oscillatory / divergent). Oscillatory -> report the band, no GCI claim.
6. **Comparisons** (`cfd/postproc.py`): Cl, Cd, Cm vs alpha against TM-4062 with experimental uncertainty bands; Cp at alpha 0, 4, 8 (RMS delta Cp on the upper surface); drag split into pressure and friction for every model (this is where fully turbulent models go wrong).
7. **LSB** (`cfd/lsb.py`): on the upper surface downstream of the suction peak, separation x_s = first zero crossing of Cf_x to negative, reattachment x_r = return to positive; cross-check with the Cp plateau. Fully turbulent models are expected to show no bubble; report that as a result, not a bug. Re 60k stretch: run LM; expect long bubbles and possible non-steady behaviour; whatever happens is reported.
8. **Trust Map** (`results/figures/trust_map_e387.png` + JSON): grid over (Re, alpha, model) with status in {validated (inside bar), out-of-bar, unsteady, near-stall (alpha > alpha_Clmax_exp - 2 deg)}. LIMITS.md names at least 3 quantified failure regions with the numbers (for example "SA over-predicts Cd by X% at Re 100k because it has no laminar run; LM bubble reattachment off by Y c at alpha 8 deg; no steady solution beyond alpha Z deg").
9. **OpenFOAM cross-check** (`cfd/openfoam_case/`): same medium mesh extruded one cell in z by gmsh (front/back `empty`), `gmshToFoam`, `simpleFoam`, `RASModel kOmegaSSTLM`; walls k = 0, omega `omegaWallFunction` (acts low-Re at y+ < 1), nut `nutLowReWallFunction`, gammaInt and ReThetat `zeroGradient`; inlet Tu matched to SU2. Re 200k, alpha 0 and 4. Done when Cl and Cd vs SU2-LM are reported and any difference above 5% (Cl) / 15% (Cd) is explained or logged.

**Unlocks:** V03-V09. **Milestone:** `flightloop-m1-e387-validated` (cut-line version tagged `flightloop-m1a` if only Re 200k exists).

### P3 Module E-synthetic: 6-DOF truth, sensors, sys-ID (runs while P2 solves)

**Done-state:** 20-seed Monte Carlo of synthetic flights; glide/energy, equation-error and output-error estimators recover truth within the bars; CI coverage measured.
**Negative scope:** no real data in this phase; no Kalman smoother beyond what is needed; no lateral-directional ID (park in GAPS.md unless time remains).

Steps:
1. `sim/sixdof.py`: quaternion rigid-body 6-DOF, RK4 at 200 Hz, linear aero model from `derivs.json` plus quadratic drag polar; ISA atmosphere from `atmos.py`; thrust model with prop efficiency and electric power (for the prop rule in F). Validate first: trimmed glide holds within 0.1 m/s and 0.1 deg for 60 s with turbulence off; energy is conserved with aero and thrust off.
2. `sim/dryden.py`: MIL-F-8785C Dryden turbulence (light) plus a steady wind. `sim/sensors.py`: noise per config, 50 Hz output, baro lag, pitot scale error.
3. `sim/maneuvers.py`: steady glides at 6 airspeeds (throttle 0), 3-2-1-1 and 2-1-1 elevator inputs with the pulse width set near the short-period period from the AVL eigenvalues, doublets. Write logs to the canonical schema with `source: synthetic`. Also one 18 km altitude case (rho 0.1207) to show the same pipeline works at stratospheric density.
4. `sysid/glide.py`: energy method. Specific energy E = h + V^2/(2g) using TAS; with zero thrust, D = -W (dE/dt) / V, CD = D / (q S), CL = W n_z / (q S). Fit CD = CD0 + k CL^2 by weighted least squares over quasi-steady windows (|dV/dt| and |q| below thresholds). Vertical air mass motion is the main bias: show it on synthetic by turning on vertical gusts.
5. `sysid/equation_error.py`: Cm = Iyy q_dot / (q_bar S c) with q_dot from a Savitzky-Golay or smoothing-spline derivative of gyro q (never raw finite differences); regress on [1, alpha, q c/(2V), de]. Same for CL on [1, alpha, q c/(2V), de]. OLS standard errors with the Morelli-Klein correction for colored residuals.
6. `sysid/output_error.py`: short-period model (alpha, q states), parameters by `scipy.optimize.least_squares` on output residuals, Cramer-Rao bounds from the Fisher matrix scaled by residual variance. Report the parameter correlation matrix and flag |rho| > 0.9 (Cm_q vs Cm_alpha_dot are not separately identifiable; estimate their sum).
7. `sysid/ci.py`: block bootstrap over maneuvers as a second CI. Coverage test: fraction of 20 seeds whose 95% CI contains truth.
8. Two alpha cases: vane (truth + noise) and no vane (alpha from attitude, inertial velocity and wind estimate). Report the bias the no-vane path introduces.

**Unlocks:** V14 (synthetic part), V15, V16. **Milestone:** `flightloop-m2-sysid-synthetic`.

**Minimum viable product = P0 + P1 (TM-4062 only) + cut-line P2 + P3 + verify.py. Everything below is above the cut-line.**

### P4 Module D: compressible verification (2-3 h, cheap and high signal)

**Done-state:** wedge shock angle and pressure ratio vs theory, diamond Cl/Cd vs shock-expansion theory, 3-grid convergence for each.
**Negative scope:** academic geometry only (10 deg wedge, 5% diamond). No bodies of revolution, fins, inlets, or vehicle shapes. No viscous supersonic runs in phase 1.

Theory targets (computed in planning, gamma = 1.4; `compressible/theory.py` must reproduce them, which is itself a check):

| Case | Theory |
|---|---|
| Wedge M 2.0, theta 10 deg | beta = 39.31 deg, p2/p1 = 1.707, M2 = 1.641 |
| Wedge M 3.0, theta 10 deg | beta = 27.38 deg, p2/p1 = 2.055, M2 = 2.505 |
| Diamond t/c 0.05, M 2.0, alpha 0 / 2 / 4 deg | shock-expansion Cl = 0 / 0.0810 / 0.1622; Cd = 0.00578 / 0.00863 / 0.01721 (Ackeret: Cl 0.0806 at 2 deg, Cd 0.00577 at 0) |
| Diamond t/c 0.05, M 3.0, alpha 0 / 2 / 4 deg | Cl = 0 / 0.0499 / 0.1001; Cd = 0.00355 / 0.00531 / 0.01063 |

Steps: gmsh meshes (wedge: structured, shock-aligned refinement optional; diamond: O-mesh or blocked, farfield 20c); SU2 `EULER`, Roe or HLLC with MUSCL and Venkatakrishnan limiter, supersonic inlet/outlet or farfield; 3 grids. Shock angle: fit a line through the locations where p crosses (p1 + p2)/2 along rays from the wedge apex. Note in LIMITS.md that shock-expansion theory is exact for surface pressures here (attached shocks, no wave reflections on the surface), so the residual is pure numerical error, while Ackeret is first-order.

**Unlocks:** V12, V13.

### P5 Module F: next-iteration engine (3 h)

**Done-state:** `loop/` turns a flight-vs-model discrepancy table into design deltas, updates the model, and re-flies; demonstrated on at least 5 blind synthetic cases.
**Negative scope:** no optimizer, no new airframe design from scratch; deltas are first-order and say so.

Rule table (`loop/deltas.py`, each rule cites its formula in the report):

| Discrepancy (flight vs model) | Attribution checks | Design delta |
|---|---|---|
| CD0 higher by X% | excrescence/interference build-up factor; profile drag at the flight-identified Re (corrected polars from C); LSB penalty if flight Re is below design Re | update CD0 build-up; recompute V_best_glide = sqrt(2W/(rho S)) (k/CD0)^(1/4), V_min_sink with (k/(3 CD0))^(1/4), (L/D)max = 1/(2 sqrt(k CD0)); if the LSB penalty dominates, recommend the airfoil from the candidate set with lowest corrected Cd at the flight CL and Re |
| Static margin lower by d (SM = -Cma/CLa) | Iyy uncertainty propagated (Cma scales with assumed Iyy); CG measurement | move CG forward by d x MAC, or raise V_h by d / (dXnp/dV_h) with the sensitivity from AVL finite differences |
| k higher (e lower) | trim drag; tail download | span-load / twist note, trim-drag estimate at measured CG |
| Trim elevator offset | Cm0 / incidence | tail incidence change = -Cm0_error / (Cmde-equivalent incidence sensitivity) |
| Powered climb short at matched drag | prop efficiency from T = D + W sin(gamma) + m dV/dt and eta = T V / (volt x curr) | prop re-match note (pitch/diameter direction) at the measured advance ratio |

Blind demo: `stubs/inject.py` draws discrepancies from a seed (CD0 +10 to +40%, CG aft 2 to 5% MAC, Cmde -15 to +5%, prop eta -5 to -15 points) into the truth model and writes them to `runs/<id>/SEALED_truth.json`, which the engine never reads (verify.py greps the engine's file access log). Engine sees only the synthetic logs and the design model, outputs deltas, updates the model, re-simulates the next flight. Score: recovered CD0 delta, SM delta, and reduction of next-flight prediction error in sink rate and trim elevator.

**Unlocks:** V18. **Milestone:** `flightloop-m3-loop-closed`.

### P6 Module C: multi-fidelity correction (3-4 h)

**Done-state:** leave-one-airfoil-out (by family) evaluation of a delta model vs raw XFOIL, constant-bias correction, and NeuralFoil, on the same points, against the bars pre-registered in PREREG.md.
**Negative scope:** no deep nets beyond a small MLP; no tuning on held-out families; E387 and SD7037 are held out from all training as final test airfoils.

Steps:
1. XFOIL batch on as-built coordinates, panels 200, Ncrit 9 (and 11 as feature ablation), `ASEQ -4 14 0.5`, BL re-init on failures. Log non-converged points; exclude them with counts reported.
2. Alignment: UIUC lift data is vs alpha and drag polars are (Cl, Cd) pairs measured separately. Align Cl corrections on alpha and Cd corrections on Cl (interpolate XFOIL at the experimental Cl). Restrict to pre-stall points (Cl < 0.9 Clmax_exp).
3. Features: log10 Re, alpha, XFOIL Cl, Cd, Cm, upper/lower transition x_tr (very informative for the bubble), max t/c and location, max camber and location. Targets: Cl_exp - Cl_xfoil and log(Cd_exp / Cd_xfoil).
4. Models: ridge (baseline), GP (Matern + white noise), HistGradientBoosting. Hyperparameters by nested GroupKFold on families. Final: leave-one-family-out over all non-test families, plus E387 and SD7037 as untouched test.
5. Baselines on identical points: raw XFOIL (Ncrit 9), raw XFOIL (Ncrit 11), constant mean-bias correction, NeuralFoil xxlarge (record its `analysis_confidence`).
6. Report Cd MAPE and Cl MAE per Re band (60k / 100k / 200k-300k). The 60k band is where everything is worst; say so with the number.

**Unlocks:** V11.

### P7 Module B: stratospheric-Re case (2 h)

**Done-state:** `results/strato/DISAGREEMENT.md` with a table: SD7037 (or StratoGlide foil) at Re 60k and 100k (plus 80k without experiment): SU2-LM, SU2-SA, XFOIL Ncrit 9 and 11, NeuralFoil, corrected XFOIL (from C), UIUC experiment; Cd error at CL 0.5 to 0.8, Clmax, bubble location (CFD and XFOIL x_tr), and a numbered explanation of each disagreement.
**Negative scope:** no new experiments, no claims about stratospheric turbulence levels beyond the Ncrit sensitivity.

Explanations to test, not assert: XFOIL's e^N envelope with viscous-inviscid coupling is calibrated for natural TS transition in quiet tunnels and handles bubbles well at this Re; LM correlations are calibrated largely on bypass transition and depend on the Tu that survives to the leading edge; NeuralFoil inherits XFOIL's bias and its authors flag a hard regime near Re 80k; Mack's relation Ncrit = -8.43 - 2.4 ln(Tu) maps Ncrit 9 to Tu 0.07% and Ncrit 11 to about Tu 0.03%, the argument for quieter stratospheric air.

**Unlocks:** V10.

### P8 Module E-real: ALFA and public logs (3 h)

**Done-state:** parsers for ArduPilot `.bin` (pymavlink `DFReader`), PX4 `.ulg` (pyulog), and ALFA processed CSVs produce canonical Parquet; glide drag polar from ALFA unpowered segments (engine-failure flights) and equation-error longitudinal derivatives from flights with usable elevator excitation; AVL (approximate T-28 geometry) vs flight table with CIs.
**Negative scope:** no fault-detection ML (ALFA's main use) in this phase; no claims of absolute accuracy (no truth exists); powered segments are not used for drag because ALFA has no thrust measurement.

Rules: stopped or windmilling propeller drag is inside the glide CD0 and is labelled "airframe + stopped prop"; elevator in normalized PWM units unless a calibration is sourced, so Cmde is "per unit command"; Iyy from a radius-of-gyration estimate with stated uncertainty, propagated into Cma and SM. Every number has a CI and a residual RMS, or it is not reported. Number of flights and segments used goes into `numbers.json` as the real count.

**Unlocks:** V14 (real part), V17.

### P9 Module G: dashboard and report (3 h)

**Done-state:** `streamlit run dashboard/app.py` shows tabs: E387 Validation (polars, Cp, Cf with bubble markers, GCI table), Trust Map, Stratospheric Re, Multi-fidelity, Compressible, Sys-ID (fits, residuals, CIs, coverage), Correlation Report (discrepancies -> deltas -> re-fly). Every tile labels its data as public-experimental, CFD, low-order or synthetic. `dashboard/export_static.py` writes a self-contained HTML. `report/FlightLoop_Report.pdf` in Vinay's document style (tables, numbered checks, exact commands), every number pulled from `results/numbers.json` by a template, zero em dashes.
**Negative scope:** no deploy, no hosting. Design bar: clean HUD/glass look, SVG icons, no emoji, plots readable in dark and light.

Flow-field plots: SU2 volume output (Paraview VTU) read with `meshio` or SU2's CSV; pressure contours and streamlines near the trailing half to show the bubble at alpha 4, Re 200k.

**Milestone:** `flightloop-v1` after `verify.py` is green.

---

## 6. verify.py (20 numbered checks, pre-registered)

Two classes. **[I] Integrity** must pass for any milestone; failing one means the system is wrong, fix the system. **[R] Result** bars are pre-registered in PREREG.md; a FAIL is allowed only if logged in HONEST_LOG.md with the actual number and the resume contract updated. Never edit a bar after results exist (verify.py checks that the commit adding PREREG.md predates every `runs/` result commit).

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

Output format: one line per check `V07 [R] PASS  SA Cd +61.3% vs TM-4062 (bar > 25%)`, then `INTEGRITY a/12, RESULTS b/8`, exit code nonzero if any [I] fails. Writes `results/verify_report.json`.

---

## 7. Timeline (2-3 days, autonomous agent) and cut-line

| When | Workstream 1 (solver queue) | Workstream 2 (python) | Gate |
|---|---|---|---|
| Day 1, h0-3 | P0 env, smoke tests, c_ci benchmark | scaffold, config, PREREG.md, verify.py skeleton | `m0-env` |
| Day 1, h3-6 | P1 TM-4062 extraction (critical path) | P1 UIUC, ALFA download (background) | SOURCES.md |
| Day 1, h6-8 | P2 meshes, y+ check, launch cut-line runs with nohup | P3 6-DOF + sensors | |
| Day 1, night | cut-line A (~3 h M4) | P3 estimators | |
| Day 2, h0-4 | GCI, LSB, Trust Map on cut-line; decide Azure from COMPUTE.md | P3 Monte Carlo, coverage | **`m1a` + `m2-sysid-synthetic` = MVP** |
| Day 2, h4-8 | full A + B matrix (VM or M4 overnight) | P4 compressible, P5 loop engine | `m3-loop-closed` |
| Day 2, night | A.7 Re 60k, OpenFOAM | P6 XFOIL batch | |
| Day 3, h0-4 | collect matrix, finish A, P7 B | P6 multi-fidelity eval | `m1-e387-validated` |
| Day 3, h4-9 | P8 E-real | P9 dashboard + report, final verify | `flightloop-v1` |

**CUT-LINE (minimum viable, never skip):** P0 + P1 (TM-4062) + P2 cut-line (Re 200k, alpha {0, 4, 8}, SA + SST-LM, GCI at alpha 4) + P3 (E-synthetic) + verify.py with V01-V09, V14 (synthetic), V15, V16, V19, V20 implemented.
**Above the cut-line, in this priority order:** D (cheap, shows compressible competence) > F (closes the loop, needs only synthetic) > C > B > E-real > G full > OpenFOAM cross-check > Re 60k > SST and SA-BC full polars > XFLR5 (manual, by Vinay).

---

## 8. Rules and guardrails

- **No guess work.** Root-cause every failed run (read the solver log, check the mesh, check the cfg diff) before changing anything. One change at a time, logged.
- **Anti-gaming.** Never hardcode answers or tune on test sets; never tune turbulence model constants or Tu to hit a bar (the Tu sweep is reported as a sensitivity, and the headline uses the pre-registered Tu); never shift alpha to match experiment; never weaken, delete or game a verifier check (fix the system or extend the data instead); null results go into HONEST_LOG.md immediately; smaller honest numbers beat bigger soft ones.
- **Every assertion carries provenance.** Every number in `numbers.json` has `{value, unit, script, run_id, git_sha, data_source}`.
- **Label the data.** Public-experimental, CFD, low-order, synthetic. Plot titles and table captions say which. No sentence anywhere implies Vinay ran a tunnel or flew these aircraft.
- **Run discipline.** Unique run IDs, isolated `runs/<run_id>/`, config snapshot per run, resumable `scripts/run_all.sh` that skips converged runs, nohup for anything over 10 minutes, per-unit time math in COMPUTE.md before big launches, preserve every artifact.
- **Degrade and continue.** A blocked tool (OpenFOAM image, AVL build, ALFA size) gets a HONEST_LOG entry and the fallback from section 4.1; do not stall the run.
- **Commit at every milestone** with conventional messages that state the mechanism, for example `feat: E387 SST-LM Re 200k validated, GCI 1.8% on Cd, bubble x_r off 0.03c`. No Co-Authored-By. Back up before anything destructive.
- **Keep PROGRESS.md current** (phase, last good commit, runs in flight with PIDs, next action) so a cold session can resume. Re-read PROGRESS.md and this plan after any context compaction; trust git and verify.py over summaries.
- **No em dashes** in any outward file (report, README, dashboard, resume contract). Commit messages may use them.
- **No weapon content.** See section 0.

---

## 9. Interview defence (15 hardest questions)

| # | Question | Concise correct answer | Repo evidence |
|---|---|---|---|
| 1 | Why does transition matter so much at Re 100k-200k? | A large share of the chord is laminar. Flat-plate Cf at Re_x = 1e5 is 0.0042 laminar (Blasius) vs 0.0074 turbulent, about 1.8x. A fully turbulent model adds that friction everywhere and misses the laminar separation bubble, which sets pressure drag, the Cp plateau and, when it bursts, stall. | V07, drag split plot in `results/figures/e387_drag_split.png` |
| 2 | Why SA, SST and LM, and what does each get wrong here? | SA: one equation, robust, calibrated for attached fully turbulent aero flows. SST: k-omega near the wall, k-epsilon outside, shear-stress limiter, better in adverse pressure gradients. Both are fully turbulent. LM adds intermittency and transition-onset Re_theta transport on SST, so it can produce natural and separation-induced transition, but it relies on empirical correlations and on the Tu reaching the leading edge. SA-BC is algebraic and built for attached transition, weak on bubbles. | Trust Map, LIMITS.md |
| 3 | Why y+ < 1? | SA, SST (omega wall condition) and LM are integrated to the wall here. Wall functions assume an equilibrium turbulent log layer, which does not exist in laminar, transitional or separated regions. LM's formulation assumes y+ around 1. I used a first cell of 2e-5 c, which gives y+ 0.1-0.4, and checked max y+ per run. | V03, `runs/*/status.json` |
| 4 | What is GCI and how did you compute it? | Grid Convergence Index (Roache, Celik 2008 procedure). Three systematically refined grids, r = sqrt(2) per direction, observed order p from the three solutions, Richardson extrapolation, GCI_fine = 1.25 |e21| / (r^p - 1). It is a numerical-error band, not a validation error. Oscillatory convergence means no GCI claim. | V05, `results/gci_e387.json` |
| 5 | Verification vs validation in your project? | Verification: are the equations solved right (GCI, wedge vs theta-beta-M, diamond vs shock-expansion, which is exact for surface pressures with attached shocks). Validation: are they the right equations (E387 vs TM-4062 tunnel data, sys-ID vs flight). | V05, V12, V13 vs V06, V08 |
| 6 | Why does XFOIL sometimes beat RANS at low Re? | XFOIL couples a panel method with integral boundary layer equations and an e^N envelope transition model calibrated for natural Tollmien-Schlichting transition in quiet tunnels, exactly this regime. RANS-LM uses correlations largely calibrated on bypass transition and is sensitive to the freestream Tu decay. XFOIL fails at massive separation and post-stall, and it is 2D only. | DISAGREEMENT.md, Trust Map |
| 7 | What Ncrit and Tu, and how did you justify them? | Mack's relation Ncrit = -8.43 - 2.4 ln(Tu): Ncrit 9 is about Tu 0.07%, a typical tunnel; quiet stratospheric air argues for higher Ncrit (11 is about Tu 0.03%). I show Ncrit 9 and 11 and an LM Tu sweep as sensitivities instead of picking the value that matches. | B table, Tu sweep in `results/e387_tu_sweep.json` |
| 8 | What are the limits of your CFD results? | Steady RANS cannot capture bubble unsteadiness, bursting or stall hysteresis; near Clmax the steady solution is not unique or does not converge; 2D CFD vs a tunnel with sidewall boundary layers (TM-4062 applies tunnel corrections, which carry their own uncertainty); LM results move with Tu. I mapped where each model is inside the bar. | LIMITS.md, V09 |
| 9 | How did you make the tunnel comparison fair? | Corrected TM-4062 data, its stated uncertainty bands, same alpha (no offset applied, offset reported), drag compared as total profile drag (the tunnel uses a wake survey), Cp compared point by point, convergence criteria fixed before comparing. | V06, `cfd/postproc.py` |
| 10 | How do you know the wedge result is right and not just pretty? | Theory at M 2, 10 deg: beta 39.31 deg, p2/p1 1.707. Shock angle fitted from the mid-pressure contour along rays, not eyeballed; three grids show convergence toward theory. Shock smearing over 2-3 cells sets the error floor. Beyond theta_max (about 23 deg at M 2) the shock detaches and the attached-shock theory does not apply. | V12, `results/figures/wedge_convergence.png` |
| 11 | What is identifiability and where did it bite? | Whether the data can separate parameters. In the short period, alpha_dot and q are nearly collinear, so Cm_q and Cm_alpha_dot cannot be separated; I estimate their sum. Doublets with poor frequency content inflate the Cm_alpha / Cm_q correlation; 2-1-1 inputs with pulse width near the short-period period fix that. I report the correlation matrix and the regressor condition number. | V15, `results/sysid/correlation_matrix.png` |
| 12 | Equation-error vs output-error? | Equation-error: linear regression on measured states, fast, but biased when regressors are noisy and it needs pitch acceleration, which I take from a smoothed derivative. Output-error: simulate the model and fit outputs, handles measurement noise, not process noise (turbulence), for which filter-error is the next step. | `sysid/equation_error.py`, `sysid/output_error.py`, V15 |
| 13 | How do you get alpha without a vane? | Air-relative velocity = inertial velocity from the EKF minus the wind estimate, rotated into body axes; alpha = atan(w/u). Its error is dominated by the wind estimate; I measured the bias on synthetic with and without a vane. | P3 step 8, `results/sysid/no_vane_bias.json` |
| 14 | How do you get a drag polar from flight logs? | Energy method in unpowered glide: D = -W (dE/dt) / V with E = h + V^2/2g using TAS, CL from load factor, fit CD = CD0 + k CL^2. Biases: vertical air mass motion, and on ALFA the stopped prop is inside CD0, so it is labelled airframe plus stopped prop. Powered segments need a thrust model, which ALFA does not have, so I did not use them. | V15, V17 |
| 15 | How does a flight-log discrepancy become a design change? | Example: flight Cma/CLa gives SM 3% MAC below design, so move CG forward 3% MAC or raise tail volume by 0.03 / (dXnp/dV_h) from AVL. CD0 30% high: update the build-up and recompute best-glide and min-sink speeds; if the low-Re bubble penalty dominates, swap the airfoil using corrected polars. I proved it on blind synthetic cases where the injected truth was sealed from the engine. | V18, `loop/deltas.py`, Correlation Report tab |

Backup probes to prepare from the repo: why the stratospheric Re drops while TAS rises (Re scales with sqrt(sigma)/mu at fixed W/S and CL, research notes Q1); why short-period damping drops at altitude (relative density parameter 2m/(rho S c) grows as 1/rho); what you would do with Kalam's own tunnel and flight data on day one (swap the data loaders, keep the verifier, rerun the Trust Map).

---

## 10. Resume contract

Stack header (one line):

```
FlightLoop - SU2, OpenFOAM, gmsh, XFOIL, AVL, NeuralFoil, pymavlink/pyulog, scikit-learn, Streamlit
```

Bullets (numbers are TARGETS in [~X]; final text uses the value in `results/numbers.json`, and a clause is dropped if its source check failed and the honest number is not worth stating):

1. Built FlightLoop, a multi-fidelity aero correlation engine. Ran SU2 RANS (SA, SST, Langtry-Menter transition) on gmsh y+<1 meshes for the Eppler 387 at Re 100k to 200k; matched NASA TM-4062 tunnel lift within [~4%] and bubble reattachment within [~0.05c], GCI [~<3%].
2. Closed the flight-log loop: parsed ArduPilot and PX4 logs, recovered synthetic CD0 within [~5%], fit drag polars on [~10] real CMU ALFA flights, and auto-generated CG and tail deltas cutting simulated next-flight error [~70%]; corrected XFOIL drag error [~35%] on held-out airfoils.

| Number | Produced by | Artifact / key in `results/numbers.json` | Check |
|---|---|---|---|
| Re 100k to 200k | config + completed runs | `e387.re_validated` | V04 |
| lift within [~4%] | `cfd/postproc.py` | `e387.lm.re200k.a4.cl_err_pct` | V06 |
| reattachment within [~0.05c] | `cfd/lsb.py` | `e387.lm.re200k.xr_err_c` (max over alpha 0, 4) | V08 |
| GCI [~<3%] | `cfd/gci.py` | `e387.gci.cd_fine_pct` (largest of SA, LM) | V05 |
| y+<1 | `cfd/postproc.py` | `e387.yplus_max` | V03 |
| CD0 within [~5%] | `sysid/glide.py` Monte Carlo | `sysid.synth.cd0_median_err_pct` | V15 |
| [~10] ALFA flights | `logs/alfa.py` + `sysid/glide.py` | `sysid.alfa.n_flights_used` | V17 |
| next-flight error cut [~70%] | `loop/correlate.py` blind demo | `loop.blind.error_reduction_pct` (median over seeds) | V18 |
| XFOIL drag error cut [~35%] | `mf/evaluate.py` | `mf.cd_mape_reduction_pct` = 1 - corrected/raw | V11 |

Words: bullet 1 = 41, bullet 2 = 42. Both checked for em dashes. "ArduPilot and PX4 logs" is only true if V14 passes on a real-format .bin and a public ULog; otherwise rewrite the clause.

---

## 11. When done (final message shape)

1. `verify.py` output pasted in full (INTEGRITY a/12, RESULTS b/8).
2. Milestone tags created, with commit SHAs.
3. Table of every headline number: value, bar, PASS/FAIL, source key.
4. HONEST_LOG.md summary: every failed bar and every fallback taken.
5. Final resume bullets with real numbers substituted, word counts, em-dash scan result.
6. What is deferred (GAPS.md) and the compute actually used (M4 hours, VM hours, VM deallocated: yes/no).

---

# THE /goal COMMAND (paste this)

```text
/goal Build FlightLoop in ~/Desktop/flightloop: a multi-fidelity aero correlation engine (SU2/OpenFOAM CFD <-> XFOIL/AVL <-> public tunnel data <-> flight logs -> design deltas). DONE when `python3 verify.py` prints INTEGRITY 12/12 and every RESULT check is PASS or a logged honest FAIL (number in HONEST_LOG.md, numbers.json and resume contract updated), shown with real command output in the transcript, and tag flightloop-v1 exists. Or stop after 60 turns with BLOCKERS.md saying exactly what is missing and why.

DIRECTIVE: follow the full spec in ~/Desktop/kalam-labs-application/master_plans/FlightLoop_MASTER_PLAN.md. Copy it to MASTER_PLAN.md in the repo and read it FIRST, then the JD and research_notes/.../aero_domain_problems.md. Re-read MASTER_PLAN.md and PROGRESS.md after any context compaction; trust git and verify.py over summaries.

PHASES (done-state + negative scope in the plan): P0 env+benchmark+PREREG.md -> P1 data (TM-4062 double transcription, UIUC LSAT, ALFA, PX4 ULogs, ArduPilot SITL .bin) -> P2 A: E387 SU2 SA/SST/SA-BC/SST-LM, gmsh C-mesh y+<1, 3-grid GCI, LSB vs oil flow, Trust Map, LIMITS.md, OpenFOAM kOmegaSSTLM cross-check -> P3 E-synthetic: own 6-DOF, Dryden, sensor noise, glide/equation-error/output-error, 20 seeds, CI coverage -> P4 D: Euler wedge+diamond vs theta-beta-M and shock-expansion -> P5 F: next-iteration engine on blind sealed synthetic cases -> P6 C: XFOIL->UIUC delta model, leave-one-family-out vs raw XFOIL and NeuralFoil -> P7 B: SD7037 Re 60k/100k disagreement -> P8 E-real: ALFA + public logs -> P9 G: Streamlit + static HTML + report PDF.

CUT-LINE (never skip): P0 + P1(TM-4062) + P2 at Re 200k, alpha 0/4/8, SA + SST-LM, GCI at alpha 4 + P3 + verify.py. Then D > F > C > B > E-real > G > OpenFOAM > Re 60k.

HARD RULES:
- Pre-registered bars live in PREREG.md, committed before any result. Never edit a bar, tune model constants/Tu/alpha to hit one, tune on test airfoils, or weaken/delete a verifier check. Fix the system or extend the data. Null results go in HONEST_LOG.md immediately. Smaller honest numbers beat bigger soft ones.
- All data is public (NASA TM-4062, UIUC LSAT, CMU ALFA, PX4 Flight Review) or synthetic, labelled on every file, plot and table. Never imply Vinay owns a tunnel or flew these aircraft. Never invent a number: write n/a (reason).
- Module D is a textbook wedge (theta 10 deg; M2 beta 39.31, M3 27.38) and 5% diamond only. NO missile, intake, seeker, guidance, targeting or warhead content.
- Every number lives in results/numbers.json with script, run_id, git_sha, source. Report, dashboard and resume read only from it.
- Unique run IDs, isolated runs/<run_id>/, resumable scripts/run_all.sh, nohup for runs over 10 min, per-unit time math in COMPUTE.md before big launches. Measure c_ci in P0; use the Azure VM only if the M4 matrix exceeds 12 h, and deallocate it after.
- Work autonomously, never stop to ask. On a blocker, log it in HONEST_LOG.md, take the fallback from the plan's install table, and continue.
- Keep PROGRESS.md current (phase, last good commit, PIDs, next action). Commit at each milestone (flightloop-m0-env, m1a, m2-sysid-synthetic, m3-loop-closed, m1-e387-validated, v1) with mechanism in the message. NO Co-Authored-By, no Claude/Anthropic attribution, ever. NO push, NO remote, NO deploy.
- Zero em dashes in the report, README, dashboard and resume contract.

WHEN DONE: paste the full verify.py output, tags with SHAs, the headline-number table (value, bar, PASS/FAIL, numbers.json key), the HONEST_LOG summary, both resume bullets with real numbers, word counts and the em-dash scan, GAPS.md deferrals, and the compute used (M4 h, VM h, VM deallocated yes/no).
```
