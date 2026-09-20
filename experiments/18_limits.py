"""Where does the method break?  A stress map around the paper's P = 160
simulation (merge scenario of 06_gradient_estimator.py: six modules, rho 0.55,
pairs merge above z*; every taxon's marginal distribution constant along z).

Axes varied one at a time from the base (P=160, S=60, a=0.7, z*=0.55,
Gaussian abundances):
  effect     mixing weight a in {0.3, 0.5, 0.7, 0.9}
  taxa       S in {12, 24, 60, 120}
  realism    Gaussian; Poisson-lognormal counts at mean occupancy 0.8, 0.4,
             0.15 (log1p); presence/absence at occupancy 0.4
  edge       true threshold z* in {0.2, 0.35, 0.55, 0.8}
  spatial    FALSE-POSITIVE rate when structure varies by region but not with
             the gradient within region, and the gradient's region means are
             arbitrary: plot-level permutation vs within-region permutation
Each cell: 30 realizations (spatial: 60), 99 permutations, alpha = 0.05.
Output: results/18_limits.txt
"""
import os

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import importlib
import multiprocessing as mp
import time
import zlib

import numpy as np

from ec import corr_distance, chi_by_edge_count

sim = importlib.import_module("06_gradient_estimator")
NPERM, NREAL, NSPATIAL = 99, 30, 60
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "18_limits.txt")
BASE = dict(P=160, S=60, a=0.7, zstar=0.55, obs="gauss")


def m_grid(S):
    return np.unique(np.linspace(int(np.ceil(2 * S / 3)),
                                 (S * (S - 1) // 2) // 2, 12).astype(int))


def latent(P, S, a_of_z, z, rng, nmod=6, rho=0.55):
    mod = np.arange(S) % nmod
    u = rng.normal(size=(nmod, P))
    v = rng.normal(size=(nmod // 2, P))
    f = np.sqrt(1 - a_of_z) * u[mod] + np.sqrt(a_of_z) * v[mod // 2]
    return np.sqrt(rho) * f + np.sqrt(1 - rho) * rng.normal(size=(S, P))


def mu_for_occupancy(occ, rng):
    g = rng.normal(size=200000)
    lo, hi = -12.0, 6.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if 1 - np.mean(np.exp(-np.exp(mid + g))) < occ:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def observe(L, obs, rng):
    if obs == "gauss":
        return L
    kind, occ = obs.split(":")
    counts = rng.poisson(np.exp(mu_for_occupancy(float(occ), rng) + L))
    return (counts > 0).astype(float) if kind == "binary" else np.log1p(counts)


def curves(X, z, mg):
    P = X.shape[1]
    W, stride = max(6, int(round(P / 4))), None
    stride = max(1, int(round(W / 3)))
    o = np.argsort(z, kind="stable")
    starts = range(0, P - W + 1, stride)
    C = np.array([chi_by_edge_count(corr_distance(X[:, o[s:s + W]]), mg)
                  for s in starts])
    zc = np.array([z[o[s:s + W]].mean() for s in starts])
    return C, zc


def perm_p(X, z, mg, rng, blocks=None):
    zhat, obs = sim.est_split_mean(*curves(X, z, mg))
    hits = 0
    for _ in range(NPERM):
        if blocks is None:
            zp = rng.permutation(z)
        else:
            zp = z.copy()
            for b in np.unique(blocks):
                i = np.flatnonzero(blocks == b)
                zp[i] = rng.permutation(z[i])
        hits += sim.est_split_mean(*curves(X, zp, mg))[1] >= obs
    return zhat, (1 + hits) / (1 + NPERM)


def job(args):
    axis, level, i = args
    rng = np.random.default_rng(zlib.crc32(f"{axis}|{level}|{i}".encode()))
    if axis == "spatial":
        P, S, K = 160, 60, 8
        region = np.repeat(np.arange(K), P // K)
        a_r = rng.uniform(0, 0.9, size=K)          # structure varies by region
        z = rng.uniform(size=K)[region] + rng.normal(0, 0.05, size=P)
        X = latent(P, S, a_r[region], z, rng)
        mg = m_grid(S)
        _, p_std = perm_p(X, z, mg, rng)
        _, p_blk = perm_p(X, z, mg, rng, blocks=region)
        return axis, level, i, p_std, p_blk
    cfg = dict(BASE, **{axis: level}) if axis != "none" else dict(BASE)
    P, S = cfg["P"], cfg["S"]
    z = np.sort(rng.uniform(size=P))
    a = np.where(z > cfg["zstar"], cfg["a"], 0.0)
    X = observe(latent(P, S, a, z, rng), cfg["obs"], rng)
    zhat, p = perm_p(X, z, m_grid(S), rng)
    return axis, level, i, zhat, p


GRID = [("a", v) for v in (0.3, 0.5, 0.7, 0.9)] + \
       [("S", v) for v in (12, 24, 60, 120)] + \
       [("obs", v) for v in ("gauss", "count:0.8", "count:0.4", "count:0.15",
                             "binary:0.4")] + \
       [("zstar", v) for v in (0.2, 0.35, 0.55, 0.8)]


if __name__ == "__main__":
    t0 = time.time()
    jobs = [(ax, lv, i) for ax, lv in GRID for i in range(NREAL)]
    jobs += [("spatial", "8 regions", i) for i in range(NSPATIAL)]
    jobs.sort(key=lambda j: -(j[1] == 120))          # longest first
    with mp.get_context("spawn").Pool(int(os.environ.get("NPROC", 4))) as pool:
        res = pool.map(job, jobs, chunksize=1)
    L = [f"stress map: base P=160, S=60, a=0.7, z*=0.55, Gaussian; "
         f"{NREAL} realizations/cell ({NSPATIAL} spatial), {NPERM} perms, alpha 0.05",
         "(Monte Carlo SE of a rejection rate r is sqrt(r(1-r)/30) <= 0.09)", ""]
    L.append(f"{'axis':>6} {'level':>11} {'power':>6} {'RMSE(z-hat)':>12}")
    for ax, lv in GRID:
        rr = [r for r in res if r[0] == ax and r[1] == lv]
        zh = np.array([r[3] for r in rr])
        zs = lv if ax == "zstar" else BASE["zstar"]
        L.append(f"{ax:>6} {str(lv):>11} {np.mean([r[4] <= .05 for r in rr]):6.2f} "
                 f"{np.sqrt(np.mean((zh - zs) ** 2)):12.3f}")
    sp = [r for r in res if r[0] == "spatial"]
    L += ["", f"spatial confounding ({NSPATIAL} null realizations, 8 regions x 20 plots):",
          f"  false-positive rate, plot-level permutation:    {np.mean([r[3] <= .05 for r in sp]):.2f}",
          f"  false-positive rate, within-region permutation: {np.mean([r[4] <= .05 for r in sp]):.2f}",
          "", f"[total {time.time() - t0:.0f}s]"]
    open(OUT, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
