# ENV: exact environment FlightLoop ran in

Snapshot taken 2026-11-18 on the machine that produced every result in `results/`.

| Component | Version | Where | How installed |
|---|---|---|---|
| macOS | 26.3.1 | MacBook Air M4, 10 cores, 16 GB | |
| SU2 | v8.6.0 "Harrier" | `~/opt/su2/bin` | native meson build with MPI (`scripts/env.sh` puts it on PATH) |
| Open MPI | 5.0.11 | Homebrew | `brew install open-mpi` |
| gmsh | 4.15.2 | Python wheel in `.venv` | `uv pip install gmsh` |
| Python | 3.11.14 | `.venv` | `uv venv --python 3.11` |
| uv | 0.11.7 | `~/.local/bin/uv` | |
| XFOIL | 6.99 | `~/Desktop/stratoglide/tools/bin/xfoil` | built by the StratoGlide repo, used read-only |
| AVL | 3.40 | `~/Desktop/stratoglide/tools/bin/avl` | built by the StratoGlide repo, used read-only |
| Python packages | see `requirements.lock` (38 pinned) | `.venv` | `uv pip install -r requirements.lock` |

Key pins: numpy 2.4.6, scipy 1.17.1, pandas 3.1.1, scikit-learn 1.9.1, neuralfoil 0.3.3, aerosandbox 4.2.10, pyulog 1.3.0, pymavlink 2.4.50, pyarrow 26.0.1, matplotlib 3.11.2, meshio 5.3.5.

## One-command checks

```bash
source scripts/env.sh          # SU2 on PATH
SU2_CFD --help | head -1       # SU2 v8.6.0 "Harrier"
.venv/bin/python verify.py     # 20 numbered checks
```

## Notes

- meshio 5.3.5 fails on SU2 v8 appended-binary VTU files ("buffer size"); `flightloop/cfd/vtu.py` is a small raw reader used instead.
- No cloud VM was used (COMPUTE.md has the per-unit time math that kept the matrix on the M4).
- Streamlit 1.67.0 added for the dashboard (`dashboard/app.py`). A first install of streamlit 1.39.0 silently downgraded pandas 3.1.1 to 2.3.3 (and pillow, packaging, rich, markdown); no pipeline script ran in that window, the pins were restored from requirements.lock and streamlit 1.67.0 installs without touching them (checked with `uv pip freeze` against the lock).
