"""Can a covariance-aware pooling of the connectance grid beat the
pre-registered statistic?

The pooled split statistic of the paper (Section 2.4) averages |Cohen d|
across the 12 edge-count grid points.  Those points are nearly collinear --
chi(m) is a cumulative sum -- so a flat average discards information that a
covariance-aware combination keeps.  This is the transferable idea from the
ecp-stats line, where the dependence across the sweep is treated explicitly
rather than as a worst case.

Candidates (all use the same windows, splits, and permutation scheme):
  mean     |Cohen d| averaged over the grid            (pre-registered)
  maha     shrinkage Mahalanobis distance between the left and right window
           means, covariance pooled within groups, Ledoit-Wolf-style shrinkage
           toward its diagonal
  maha_d   the same on the per-grid-point standardized curves
  cusum    maximum over splits of the whitened cumulative deviation
Reported: rejection at alpha=0.05 on THRESH, false positives on NULL, and
RMSE of z-hat, at P in {24, 40, 80, 160}.  P=24 and P=40 are the regime of
the temporal series, where the pre-registered statistic detected nothing.
Output: results/23_statistic_upgrade.txt
"""
import os

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import importlib
import multiprocessing as mp
import time

import numpy as np

sim = importlib.import_module("06_gradient_estimator")

NPERM, NREAL = 99, 40
PS = (24, 40, 80, 160)
KEYS = ("mean", "maha", "maha_d", "cusum")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "23_statistic_upgrade.txt")


def _shrink_cov(A, B):
    """pooled within-group covariance with shrinkage toward its diagonal."""
    n1, n2 = len(A), len(B)
    Ac, Bc = A - A.mean(0), B - B.mean(0)
    S = (Ac.T @ Ac + Bc.T @ Bc) / max(n1 + n2 - 2, 1)
    d = np.diag(np.diag(S))
    # shrinkage weight from the effective sample size (few windows -> more shrink)
    lam = min(1.0, S.shape[0] / max(n1 + n2 - 2, 1))
    return (1 - lam) * S + lam * d + 1e-9 * np.eye(S.shape[0])


def stat_mean(A, B):
    p = np.sqrt((A.var(0, ddof=1) + B.var(0, ddof=1)) / 2) + 1e-12
    return float(np.abs((A.mean(0) - B.mean(0)) / p).mean())


def stat_maha(A, B):
    S = _shrink_cov(A, B)
    d = A.mean(0) - B.mean(0)
    try:
        v = np.linalg.solve(S, d)
    except np.linalg.LinAlgError:
        return stat_mean(A, B)
    return float(np.sqrt(max(d @ v, 0.0)))


def stat_maha_d(A, B):
    s = np.concatenate([A, B]).std(0, ddof=1) + 1e-12
    return stat_maha(A / s, B / s)


def split_stat(C, zc, kind):
    f = {"mean": stat_mean, "maha": stat_maha, "maha_d": stat_maha_d}[kind]
    best, bz = -np.inf, np.nan
    for s in range(2, len(zc) - 1):
        v = f(C[:s], C[s:])
        if v > best:
            best, bz = v, (zc[s - 1] + zc[s]) / 2
    return bz, best


def cusum_stat(C, zc):
    """whitened cumulative-sum statistic over the window sequence."""
    s = C.std(0, ddof=1) + 1e-12
    Z = (C - C.mean(0)) / s
    K = len(zc)
    best, bz = -np.inf, np.nan
    tot = Z.sum(0)
    for k in range(2, K - 1):
        part = Z[:k].sum(0)
        v = np.abs(part - k * tot / K).sum() / np.sqrt(k * (K - k) / K)
        if v > best:
            best, bz = v, (zc[k - 1] + zc[k]) / 2
    return bz, best


def evaluate(C, zc, kind):
    return cusum_stat(C, zc) if kind == "cusum" else split_stat(C, zc, kind)


def job(args):
    P, scen, i = args
    W, stride = sim.design(P)
    X, z = sim.community(P, scen, np.random.default_rng(51000 + 97 * P + i))
    C, zc = sim.window_curves(X, z, W, stride)
    obs = {k: evaluate(C, zc, k) for k in KEYS}
    hits = dict.fromkeys(KEYS, 0)
    rng = np.random.default_rng(52000 + 97 * P + i)
    for _ in range(NPERM):
        Cb, zb = sim.window_curves(X, rng.permutation(z), W, stride)
        for k in KEYS:
            hits[k] += evaluate(Cb, zb, k)[1] >= obs[k][1]
    return P, scen, {k: (obs[k][0], (1 + hits[k]) / (1 + NPERM)) for k in KEYS}


if __name__ == "__main__":
    t0 = time.time()
    jobs = [(P, s, i) for P in PS for s in ("THRESH", "NULL")
            for i in range(NREAL)]
    with mp.get_context("spawn").Pool(int(os.environ.get("NPROC", 4))) as pool:
        res = pool.map(job, jobs, chunksize=1)
    L = [f"pooling statistics compared, {NREAL} realizations, {NPERM} perms, "
         f"alpha 0.05 (MC SE of a rate <= {0.5 / np.sqrt(NREAL):.2f})", "",
         f"{'P':>5} {'statistic':>9} {'power':>7} {'FPR':>6} {'RMSE':>7}"]
    for P in PS:
        for k in KEYS:
            th = [v[k] for p_, s, v in res if p_ == P and s == "THRESH"]
            nu = [v[k] for p_, s, v in res if p_ == P and s == "NULL"]
            zh = np.array([t[0] for t in th])
            L.append(f"{P:5d} {k:>9} "
                     f"{np.mean([t[1] <= .05 for t in th]):7.2f} "
                     f"{np.mean([t[1] <= .05 for t in nu]):6.2f} "
                     f"{np.sqrt(np.mean((zh - sim.ZSTAR) ** 2)):7.3f}")
        L.append("")
    L.append(f"[total {time.time() - t0:.0f}s]")
    open(OUT, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
