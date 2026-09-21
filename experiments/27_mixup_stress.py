"""Does the intersection-ECP split test survive what broke the registered one?

Two questions, both at P = 160 on the merge gradient:
  data type  Gaussian abundances, Poisson-lognormal counts at mean occupancy
             0.8 / 0.4 / 0.15, and presence/absence -- the axis on which the
             registered statistic fell from 0.60 to 0.07 (18_limits.py)
  spatial    false-positive rate when community structure varies by region
             and the gradient's region means are arbitrary, for plot-level
             and within-region permutation -- where the registered statistic
             inflated to 0.15
Output: results/27_mixup_stress.txt
"""
import os

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import importlib
import multiprocessing as mp
import time

import numpy as np

lim = importlib.import_module("18_limits")
mix = importlib.import_module("26_mixup_split")

P, NPERM, NREAL, NSPAT = 160, 99, 40, 60
OBS = ("gauss", "count:0.8", "count:0.4", "count:0.15", "binary:0.4")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "27_mixup_stress.txt")


def scan_perm(X, z, rng, nperm, blocks=None):
    obs = mix.scan(X, z)
    hits = 0
    for _ in range(nperm):
        if blocks is None:
            zp = rng.permutation(z)
        else:
            zp = z.copy()
            for b in np.unique(blocks):
                i = np.flatnonzero(blocks == b)
                zp[i] = rng.permutation(z[i])
        hits += mix.scan(X, zp) <= obs
    return (1 + hits) / (1 + nperm)


def job(args):
    kind, level, i = args
    rng = np.random.default_rng(81000 + i)
    if kind == "spatial":
        K, S = 8, 60
        region = np.repeat(np.arange(K), P // K)
        a_r = rng.uniform(0, 0.9, size=K)
        z = rng.uniform(size=K)[region] + rng.normal(0, 0.05, size=P)
        X = lim.latent(P, S, a_r[region], z, rng)
        return kind, level, scan_perm(X, z, rng, NPERM), \
            scan_perm(X, z, np.random.default_rng(82000 + i), NPERM, region)
    z = np.sort(rng.uniform(size=P))
    a = np.where(z > lim.BASE["zstar"], lim.BASE["a"], 0.0)
    X = lim.observe(lim.latent(P, 60, a, z, rng), level, rng)
    return kind, level, scan_perm(X, z, rng, NPERM), None


if __name__ == "__main__":
    t0 = time.time()
    jobs = [("obs", o, i) for o in OBS for i in range(NREAL)]
    jobs += [("spatial", "8 regions", i) for i in range(NSPAT)]
    with mp.get_context("spawn").Pool(int(os.environ.get("NPROC", 16))) as pool:
        res = pool.map(job, jobs, chunksize=1)
    old = {"gauss": 0.60, "count:0.8": 0.37, "count:0.4": 0.13,
           "count:0.15": 0.07, "binary:0.4": 0.13}
    L = [f"Intersection-ECP split test under the conditions that broke the "
         f"registered statistic (P={P}, {NREAL} realizations, {NPERM} perms)", "",
         f"{'data type':>12} {'mixup power':>12} {'registered':>11}"]
    for o in OBS:
        v = [r[2] for r in res if r[0] == "obs" and r[1] == o]
        L.append(f"{o:>12} {np.mean([p <= .05 for p in v]):12.2f} {old[o]:11.2f}")
    sp = [r for r in res if r[0] == "spatial"]
    L += ["", f"spatial confounding ({NSPAT} null realizations, 8 regions):",
          f"  false positives, plot-level permutation:    "
          f"{np.mean([r[2] <= .05 for r in sp]):.2f}   (registered 0.15)",
          f"  false positives, within-region permutation: "
          f"{np.mean([r[3] <= .05 for r in sp]):.2f}   (registered 0.07)",
          "", f"[total {time.time() - t0:.0f}s]"]
    open(OUT, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
