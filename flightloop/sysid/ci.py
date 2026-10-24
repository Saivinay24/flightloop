"""Block bootstrap over maneuvers (segments) as a second confidence interval, and coverage helpers."""
import numpy as np


def block_bootstrap(fn, blocks, n_boot=200, seed=0):
    """fn(list_of_blocks) -> dict of estimates. Resample blocks with replacement."""
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(blocks), len(blocks))
        try:
            draws.append(fn([blocks[i] for i in idx]))
        except Exception:
            continue
    keys = draws[0].keys()
    return {k: (float(np.percentile([d[k] for d in draws], 2.5)), float(np.percentile([d[k] for d in draws], 97.5))) for k in keys}


def covers(est, se, truth, z=1.96):
    return bool(abs(est - truth) <= z * se)
