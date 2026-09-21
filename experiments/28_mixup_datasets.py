"""Second pass: the intersection-ECP estimator on every dataset and series.

Registered in amendment 5 before any of these runs.  The first-generation
results (19_multi_*.txt, 17_everglades_lt_primary.txt, 22_temporal_*.txt)
stand unchanged; this reports the same data through the better instrument.

  python3 28_mixup_datasets.py D1 D5 ELT T2 ...
Output: results/28_mixup_<name>.txt
"""
import os

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import importlib
import sys
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
multi = importlib.import_module("19_multi")
temporal = importlib.import_module("22_temporal")
lt = importlib.import_module("17_everglades_lt")
mix = importlib.import_module("26_mixup_split")

NPERM = 199
OUT = os.path.join(ROOT, "results")
TAG = "_reg" if os.environ.get("MIXUP_FRACS") else ""


def located(X, z):
    """the split (as a gradient value) minimizing the statistic, and that value."""
    o = np.argsort(z, kind="stable")
    best, bz = np.inf, np.nan
    for f in mix.FRACS:
        cut = max(3, int(round(f * len(z))))
        T = mix.stat_at(mix.corr_profiles(X[:, o[:cut]]),
                        mix.corr_profiles(X[:, o[cut:]]))
        if T < best:
            best, bz = T, float((z[o[cut - 1]] + z[o[cut]]) / 2)
    return bz, best


def boot_ci(X, z, groups, nboot=199, seed=11):
    """percentile interval for the located split, resampling units with
    replacement (all of a unit's columns travel together)."""
    rng = np.random.default_rng(seed)
    units = np.unique(groups)
    zu = np.array([z[groups == u][0] for u in units])
    idx = {u: np.flatnonzero(groups == u) for u in units}
    out = []
    for _ in range(nboot):
        pick = rng.choice(units, len(units), replace=True)
        cols = np.concatenate([idx[u] for u in pick])
        zb = np.concatenate([np.full(len(idx[u]), zu[u]) for u in pick])
        out.append(located(X[:, cols], zb)[0])
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def run(name):
    t0 = time.time()
    if name == "ELT":
        dens, keep, z_psu, psu = lt.build()
        X = np.log1p(dens.values.T)
        z = z_psu.loc[psu].values
        groups = pd.factorize(psu)[0]
        # the PSU key is REGION_SITE_PLOT, so the region is its first field;
        # the within-region test is the one that matters here too.
        reg = np.array([u.split("_")[0] for u in pd.unique(psu)])
        wide = dens
    elif name.startswith(("D",)):
        w, zc, reg_s, zname = multi.BUILD[name]()
        X, z, groups, units, zunit, keep, ncand, wide = multi.prepare(w, zc)
        reg = None if reg_s is None else reg_s.reindex(pd.unique(units)).fillna("NA").values
    else:
        w = temporal.SERIES[name]()
        X, z, groups, keep, ncand, wide = temporal.prepare(w)
        reg = None
    zero = float((wide.values == 0).mean())
    zhat, T = located(X, z)
    rng = np.random.default_rng(3)
    units = np.unique(groups)
    zu = np.array([z[groups == u][0] for u in units])
    hits = 0
    for _ in range(NPERM):
        zp = rng.permutation(zu)[groups]
        hits += located(X, zp)[1] <= T
    p = (1 + hits) / (1 + NPERM)
    lo, hi = boot_ci(X, z, groups)
    L = [f"{name}: {len(units)} units, {X.shape[1]} columns, {X.shape[0]} taxa, "
         f"sparsity {zero:.2f}",
         f"mixup split: z-hat = {zhat:.4g}  T = {T:.4f}  p = {p:.3f} "
         f"({NPERM} unit permutations)  95% CI [{lo:.4g}, {hi:.4g}]",
         f"split grid: {len(mix.FRACS)} quantiles from {mix.FRACS[0]} to {mix.FRACS[-1]}"]
    if reg is not None:
        rng2 = np.random.default_rng(4)
        hits = 0
        for _ in range(NPERM):
            zp_u = zu.copy()
            for b in np.unique(reg):
                i = np.flatnonzero(reg == b)
                zp_u[i] = rng2.permutation(zu[i])
            hits += located(X, zp_u[groups])[1] <= T
        L.append(f"within-region ({len(set(reg))} regions): "
                 f"p = {(1 + hits) / (1 + NPERM):.3f}")
    L.append(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(OUT, f"28_mixup_{name}{TAG}.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L), flush=True)


if __name__ == "__main__":
    for n in sys.argv[1:]:
        run(n)
