# COMPUTE

Machine: MacBook Air M4, 16 GB, 10 cores (4 P + 6 E). SU2 v8.6.0 native build (meson, MPI via Homebrew open-mpi).

## Measured c_ci (P0)

| case | cells | procs | s / iteration | c_ci (us per cell-iter-core) |
|---|---|---|---|---|
| NACA0012 INC SA tutorial mesh (897x257) | 230,144 | 2 | 0.48 | 4.2 |
| E387 medium, INC SST-LM, Re 200k, alpha 4 | 47,264 | 2 | 0.18 | 7.6 |

## Per-unit time math before the cut-line launch

wall = cells x iters x c_ci / procs. Two jobs at 2 MPI ranks each (memory-bounded: each 2D job is about 200 MB).

| run | cells | iter cap | est wall at cap |
|---|---|---|---|
| SST-LM medium | 47k | 30,000 | 90 min |
| SST-LM coarse | 23k | 30,000 | 45 min |
| SST-LM fine | 95k | 30,000 | 180 min |
| SA medium | 47k | 12,000 | 20 min (c_ci 4.2) |
| SA coarse / fine | 23k / 95k | 12,000 | 10 / 40 min |

Cut-line (10 runs) worst case at caps: about 8.5 h of single-job time, about 4.3 h wall on 2 parallel jobs; runs that converge early cost less. Below the 12 h Azure threshold: M4 only.
