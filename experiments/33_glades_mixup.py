"""Amendment 7: the interaction test on the TITAN2 `glades` benchmark.

The registered curve statistic does not reject on this matrix (p = 0.30,
results/09_glades.txt).  This runs the amendment-5 estimator on the same
126 sites with the same preprocessing, on both split grids.

  python3 33_glades_mixup.py
Output: results/33_glades_mixup.txt
"""
import os

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import importlib
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
mix = importlib.import_module("26_mixup_split")
ds = importlib.import_module("28_mixup_datasets")

NPERM = 199
SMAX = 100                      # the cap every other dataset carries


def build():
    taxa = pd.read_csv(os.path.join(ROOT, "data", "glades_taxa.csv"))
    z = pd.read_csv(os.path.join(ROOT, "data", "glades_env.csv"))["TP.ugL"].values
    occ = (taxa.values > 0).sum(0)
    ok = np.flatnonzero(occ >= 5)
    keep = ok[np.argsort(-occ[ok], kind="stable")][:SMAX]
    X = np.log1p(taxa.values[:, keep].T.astype(float))   # taxa x sites
    return X, np.asarray(z, dtype=float), taxa.columns[keep]


def run_grid(X, z, fracs, tag):
    mix.FRACS = tuple(fracs)
    t0 = time.time()
    zhat, T = ds.located(X, z)
    rng = np.random.default_rng(3)
    hits = sum(ds.located(X, rng.permutation(z))[1] <= T for _ in range(NPERM))
    p = (1 + hits) / (1 + NPERM)
    groups = np.arange(len(z))                 # one column per site
    lo, hi = ds.boot_ci(X, z, groups)
    return (f"{tag}: z-hat = {zhat:.4g} ug/L  T = {T:.4f}  p = {p:.3f}  "
            f"95% CI [{lo:.4g}, {hi:.4g}]  "
            f"({len(fracs)} quantiles {fracs[0]:g}-{fracs[-1]:g})  "
            f"[{time.time() - t0:.0f}s]")


if __name__ == "__main__":
    X, z, names = build()
    zero = float((X == 0).mean())
    L = [f"TITAN2 glades benchmark: {X.shape[1]} sites, {X.shape[0]} taxa "
         f"(capped at {SMAX}), sparsity {zero:.2f}, TP "
         f"{z.min():.1f}-{z.max():.1f} ug/L",
         "registered curve statistic on this matrix: z-hat = 12.5, p = 0.30",
         run_grid(X, z, tuple(np.round(np.arange(0.15, 0.851, 0.05), 2)), "widened grid"),
         run_grid(X, z, (0.3, 0.4, 0.5, 0.6, 0.7), "registered grid")]
    out = os.path.join(ROOT, "results", "33_glades_mixup.txt")
    open(out, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
