"""The frozen estimator on real ecological data: the Everglades TP gradient.

Dataset: TITAN2's `glades` data -- 126 marsh sites x 164 macroinvertebrate
taxa along a surface-water total-phosphorus gradient (2.5-169.5 ug/L), the
demonstration dataset of Baker & King (2010).  Chosen for exactly that
reason, BEFORE any outcome was seen: it is the canonical dataset of the
incumbent method, it carries a real anthropogenic enrichment gradient, and a
head-to-head on it is the comparison an MEE reviewer will trust most.

Honesty header.  The estimator (window rule, edge-count grid, pooled split
statistic, permutation and bootstrap schemes) was frozen on simulation on
2026-09-01 -- commit b0ad6cc, before this dataset was touched -- and is
applied here unchanged.  The application itself is not pre-registered (the
dataset choice and the preprocessing below were fixed on 2026-09-01, before
outcomes, but by the same person who then ran them); the Butte Hill protocol
in the paper remains the pre-registered test.

Preprocessing, fixed a priori:
  - taxa occurring in >= 5 sites (TITAN's minSplt convention; all 164 pass)
  - X = log1p(abundance), taxa x sites; correlation distance then removes
    per-taxon location and scale
  - gradient z = TP.ugL; sites ordered by z; ties broken by file order

Outputs: z-hat (the TP value at which co-occurrence structure changes
fastest), permutation p (NPERM label shuffles of site positions), bootstrap
percentile CI, for both the chi summary and the distribution baseline r(m);
TITAN's own change points on the identical matrix come from
09_glades_titan.R.
"""
import os
import sys
import time

import numpy as np
import pandas as pd

from ec import corr_distance, chi_by_edge_count, cohen_d  # noqa: F401
import importlib

sim = importlib.import_module("06_gradient_estimator")

NPERM = int(os.environ.get("NPERM", 199))
NBOOT = int(os.environ.get("NBOOT", 199))
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

taxa = pd.read_csv(os.path.join(ROOT, "data", "glades_taxa.csv"))
env = pd.read_csv(os.path.join(ROOT, "data", "glades_env.csv"))["TP.ugL"].values
keep = (taxa.values > 0).sum(0) >= 5
X = np.log1p(taxa.values[:, keep].T.astype(float))     # taxa x sites
S, P = X.shape
W = max(6, int(round(P / 4)))
STRIDE = max(1, int(round(W / 3)))
M_GRID = np.unique(np.linspace(int(np.ceil(2 * S / 3)),
                               (S * (S - 1) // 2) // 2, 12).astype(int))
print(f"{S} taxa x {P} sites, TP {env.min():.1f}-{env.max():.1f} ug/L")
print(f"windows of {W} sites, stride {STRIDE}, m in "
      f"[{M_GRID[0]}, {M_GRID[-1]}]\n")


def curves(z, summary):
    o = np.argsort(z, kind="stable")
    starts = range(0, P - W + 1, STRIDE)
    C, zc = [], []
    for s in starts:
        idx = o[s:s + W]
        D = corr_distance(X[:, idx])
        if summary == "chi":
            C.append(chi_by_edge_count(D, M_GRID, max_dim=3))
        else:                                            # r(m) quantiles
            from scipy.spatial.distance import squareform
            d = np.sort(squareform(D, checks=False))
            C.append(d[np.minimum(M_GRID, len(d) - 1)])
        zc.append(z[idx].mean())
    return np.array(C), np.array(zc)


def analyse(summary):
    t0 = time.time()
    C, zc = curves(env, summary)
    zhat, obs = sim.est_split_mean(C, zc)
    hits = 0
    for b in range(NPERM):
        r = np.random.default_rng(11000 + b)
        _, s = sim.est_split_mean(*curves(r.permutation(env), summary))
        hits += s >= obs
    pval = (1 + hits) / (1 + NPERM)
    o = np.argsort(env, kind="stable")
    bhat = []
    for b in range(NBOOT):
        r = np.random.default_rng(21000 + b)
        starts = range(0, P - W + 1, STRIDE)
        Cb, zcb = [], []
        for s in starts:
            idx = r.choice(o[s:s + W], W, replace=True)
            D = corr_distance(X[:, idx])
            if summary == "chi":
                Cb.append(chi_by_edge_count(D, M_GRID, max_dim=3))
            else:
                from scipy.spatial.distance import squareform
                d = np.sort(squareform(D, checks=False))
                Cb.append(d[np.minimum(M_GRID, len(d) - 1)])
            zcb.append(env[idx].mean())
        bhat.append(sim.est_split_mean(np.array(Cb), np.array(zcb))[0])
    lo, hi = np.percentile(bhat, [2.5, 97.5])
    print(f"{summary:>4}: z-hat = {zhat:.1f} ug/L TP   stat {obs:.2f}   "
          f"p = {pval:.3f} ({NPERM} perms)   "
          f"95% CI [{lo:.1f}, {hi:.1f}] ({NBOOT} boots)   "
          f"[{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    which = sys.argv[1:] or ["chi", "rq"]
    for summary in which:
        analyse(summary)
