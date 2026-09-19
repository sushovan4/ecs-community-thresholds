"""EXPLORATORY: the TITAN-anchored test on the glades data.

Run AFTER the primary result (p = 0.30) was seen, so this is labeled
exploratory wherever it is reported; its pre-registered (untainted) home is
the Butte Hill protocol (paper sec. 5). The test: the pooled statistic D(s)
at the single split whose boundary midpoint lies nearest TITAN's decliner
change point (15.1 ug/L, regenerated on the identical matrix), against the
same permutation null. Permuting plot positions leaves the sorted dose values
-- hence the window positions -- unchanged, so the anchored split index is
constant across permutations: a genuine single-look test.
"""
import importlib
import time

import numpy as np

gl = importlib.import_module("09_glades")

TITAN_CP = 15.1
NPERM = 199

C, zc = gl.curves(gl.env, "chi")
K = len(zc)
mids = [(zc[s - 1] + zc[s]) / 2 for s in range(2, K - 1)]
s0 = int(np.argmin(np.abs(np.array(mids) - TITAN_CP))) + 2
print(f"windows at {np.round(zc, 1)}")
print(f"anchored split s0={s0}: boundary midpoint "
      f"{(zc[s0-1] + zc[s0]) / 2:.1f} ug/L (TITAN decliner cp {TITAN_CP})")


def d_at(Cw, s):
    a, b = Cw[:s], Cw[s:]
    pooled = np.sqrt((a.var(0, ddof=1) + b.var(0, ddof=1)) / 2) + 1e-12
    return float(np.abs((a.mean(0) - b.mean(0)) / pooled).mean())


obs = d_at(C, s0)
t0 = time.time()
hits = 0
for b in range(NPERM):
    r = np.random.default_rng(51000 + b)
    Cb, _ = gl.curves(r.permutation(gl.env), "chi")
    hits += d_at(Cb, s0) >= obs
p = (1 + hits) / (1 + NPERM)
print(f"anchored D = {obs:.2f}   p = {p:.3f} ({NPERM} perms)   "
      f"[{time.time() - t0:.0f}s]")
print("EXPLORATORY -- see docstring; confirmatory version lives in the "
      "Butte protocol.")
