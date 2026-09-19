"""Repeat-survey power at publication strength (supersedes the 15-realization
runs in 11_design_extras.py, whose Monte Carlo SE was ~0.12 per cell).

Same community model and permutation schemes as 11_design_extras.py:
  averaged  years averaged within plots, P = 80 columns
  stacked   plot-years as columns, windows of W plots carry W*T columns,
            the permutation moves whole plots (years stay together)
50 realizations per cell, 199 permutations, T in {1, 3, 5, 10},
parallel over cells.  Output: results/16_multiyear_full.txt
"""
import importlib
import os
import time
from multiprocessing import Pool

import numpy as np

sim = importlib.import_module("06_gradient_estimator")
ext = importlib.import_module("11_design_extras")

P, NREAL, NPERM = 80, 50, 199
TS = (1, 3, 5, 10)
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "16_multiyear_full.txt")


def cell(args):
    mode, T, scen, i = args
    est = sim.ESTIMATORS[sim.CHOSEN]
    if mode == "averaged":
        W, stride = sim.design(P)
        X, z = ext.community_years(P, T, scen,
                                   np.random.default_rng(73600 + 100 * T + i))
        r = np.random.default_rng(74600 + 100 * T + i)
        _, obs = est(*sim.window_curves(X, z, W, stride))
        hits = sum(est(*sim.window_curves(X, r.permutation(z), W, stride))[1]
                   >= obs for _ in range(NPERM))
    else:
        W = max(6, int(round(P / 4))) * T
        stride = max(1, W // 3)
        X, zc = ext.community_years_stacked(
            P, T, scen, np.random.default_rng(83600 + 100 * T + i))
        r = np.random.default_rng(84600 + 100 * T + i)
        _, obs = est(*sim.window_curves(X, zc, W, stride))
        hits = 0
        for _ in range(NPERM):
            zp = r.permutation(zc.reshape(P, T)).reshape(-1)
            hits += est(*sim.window_curves(X, zp, W, stride))[1] >= obs
    return mode, T, scen, (1 + hits) / (1 + NPERM)


if __name__ == "__main__":
    t0 = time.time()
    jobs = [(m, T, s, i) for m in ("averaged", "stacked") for T in TS
            for s in ("THRESH", "NULL") for i in range(NREAL)]
    with Pool(os.cpu_count()) as pool:
        res = pool.map(cell, jobs, chunksize=1)
    with open(OUT, "w") as f:
        f.write(f"repeat-survey power at P={P}, tau={ext.TAU}, {NPERM} perms, "
                f"{NREAL} realizations per cell, alpha=0.05\n")
        f.write("(Monte Carlo SE of a rate p is sqrt(p(1-p)/50) <= 0.071)\n\n")
        f.write(f"{'mode':>9} {'T':>3} {'power':>7} {'FPR':>6}\n")
        for m in ("averaged", "stacked"):
            for T in TS:
                pw = np.mean([p <= 0.05 for mm, tt, s, p in res
                              if (mm, tt, s) == (m, T, "THRESH")])
                fp = np.mean([p <= 0.05 for mm, tt, s, p in res
                              if (mm, tt, s) == (m, T, "NULL")])
                f.write(f"{m:>9} {T:3d} {pw:7.2f} {fp:6.2f}\n")
        f.write(f"\n[total {time.time() - t0:.0f}s]\n")
    print(open(OUT).read())
