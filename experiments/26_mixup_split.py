"""A second-generation split test built on the Intersection ECP (mixup).

Why this might beat the paper's estimator.  The registered statistic compares
two GROUPS OF WINDOW CURVES: its effective sample size is the number of
windows, about ten, however many plots the survey has.  The intersection
Euler characteristic profile instead measures how much two clouds of points
overlap, and here the points are TAXA: each taxon is represented by its
correlation profile on one side of a candidate split, so the sample size
becomes the number of taxa (tens to a hundred).  The same taxa appear on
both sides, so the test can be paired and exact.

Construction, for a candidate split of the gradient-ordered units:
  left/right correlation matrices  ->  taxon i is the point C_left[i] on one
  side and C_right[i] on the other  ->  PCA to d0 dimensions on the pooled
  points (label-blind, so the test stays exact)  ->  two clouds of S points
  ->  Dchi(r) = chi(U(A;r) cap U(B;r)) and T = int|Dchi| dr / r_max.
Structure that reorganizes across the split SEPARATES the clouds, so small T
is evidence, and the null permutes the units along the gradient and rebuilds
everything, as elsewhere in this paper.  (A per-taxon sign-flip null was
tried first and is INVALID: the observed clouds are internally coherent,
all-left against all-right, while a flipped cloud is a mixture, so the
observed statistic is smaller than its flipped replicates whether or not
anything changed -- it rejected on null gradients at the same rate as on
threshold gradients.)  The split is chosen by minimizing T, and every null
replicate repeats that minimization.

  python3 26_mixup_split.py compare     power against the registered statistic
Output: results/26_mixup_split.txt
"""
import os

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import importlib
import multiprocessing as mp
import sys
import time

import numpy as np

sys.path.insert(0, os.environ.get("MIXUP_PATH",
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "vendor")))
from intersection_ecp.core import (default_r_grid, intersection_profile,  # noqa: E402
                                   profile_stat)
from intersection_ecp.report import project  # noqa: E402

sim = importlib.import_module("06_gradient_estimator")

D0, NFLIP, NREAL = 4, 199, 40
PS = (24, 40, 80, 160)
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "26_mixup_split.txt")


def corr_profiles(X):
    """S x S correlation matrix; row i is taxon i's profile."""
    Xc = X - X.mean(1, keepdims=True)
    Xc /= (np.linalg.norm(Xc, axis=1, keepdims=True) + 1e-12)
    return np.clip(Xc @ Xc.T, -1.0, 1.0)


def split_points(X, z, frac):
    """taxon profiles left and right of the quantile split at `frac`."""
    o = np.argsort(z, kind="stable")
    cut = max(3, int(round(frac * len(z))))
    return corr_profiles(X[:, o[:cut]]), corr_profiles(X[:, o[cut:]])


def stat_at(A, B):
    """T = normalized integral of the intersection profile (small = separated)."""
    P = project(np.vstack([A, B]), d0=D0)
    a, b = P[:len(A)], P[len(A):]
    g = default_r_grid([a, b], num=60)
    return profile_stat(intersection_profile([a, b], g), g, "intnorm")


FRACS = (0.3, 0.4, 0.5, 0.6, 0.7)


def scan(X, z):
    """min over candidate splits of T, rebuilding profiles from (X, z)."""
    return min(stat_at(*split_points(X, z, f)) for f in FRACS)


def job(args):
    P, scen, i = args
    X, z = sim.community(P, scen, np.random.default_rng(71000 + 31 * P + i))
    obs = scan(X, z)
    rng = np.random.default_rng(72000 + 31 * P + i)
    hits = sum(scan(X, rng.permutation(z)) <= obs for _ in range(NFLIP))
    return P, scen, (1 + hits) / (1 + NFLIP)


if __name__ == "__main__":
    t0 = time.time()
    jobs = [(P, s, i) for P in PS for s in ("THRESH", "NULL")
            for i in range(NREAL)]
    with mp.get_context("spawn").Pool(int(os.environ.get("NPROC", 6))) as pool:
        res = pool.map(job, jobs, chunksize=1)
    reg = {24: (0.03, 0.10), 40: (0.03, 0.07), 80: (0.25, 0.07), 160: (0.57, 0.05)}
    L = [f"Intersection-ECP split test vs the registered statistic; "
         f"{NREAL} realizations, {NFLIP} unit permutations, d0={D0}, alpha 0.05",
         "(registered figures from results/23_statistic_upgrade.txt)", "",
         f"{'P':>5} {'mixup power':>12} {'mixup FPR':>10} "
         f"{'registered power':>17} {'registered FPR':>15}"]
    for P in PS:
        th = [p for p_, s, p in res if p_ == P and s == "THRESH"]
        nu = [p for p_, s, p in res if p_ == P and s == "NULL"]
        L.append(f"{P:5d} {np.mean([p <= .05 for p in th]):12.2f} "
                 f"{np.mean([p <= .05 for p in nu]):10.2f} "
                 f"{reg[P][0]:17.2f} {reg[P][1]:15.2f}")
    L += ["", f"[total {time.time() - t0:.0f}s]"]
    open(OUT, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
