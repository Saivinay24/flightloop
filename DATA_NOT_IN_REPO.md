# Large files kept out of this repository

These inputs and raw solver dumps are over 10 MB each and are not committed. Every derived result built from them is in `results/`.

| File | What | Source |
|---|---|---|
| `data/alfa/raw/dataflash.zip`, `data/alfa/raw/processed.zip` | CMU ALFA fixed-wing flight logs | AirLab ALFA dataset (Keipour et al., IJRR 2021), see `data/SOURCES.md` |
| `data/px4/*.ulg`, `data/px4/dbinfo.json` | Public PX4 flight logs | PX4 Flight Review public log database, see `data/SOURCES.md` |
| `runs/*/history.csv`, `vol_solution.dat`, `restart.csv` | Raw SU2 convergence histories and volume solutions | Regenerate by running the CFD stage (see README) |
