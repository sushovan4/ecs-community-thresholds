"""Does the intersection-ECP test LOCALIZE, or only detect?

The registered curve statistic recovers a known threshold with RMSE about
0.10 of gradient span at P = 160.  The interaction test detects far more
often, but on real data its located split moved when the split grid changed
and its bootstrap interval spanned most of the gradient.  Here both are run
on the same simulated merge gradients and scored on recovery of z* = 0.55.
Output: results/29_mixup_localization.txt
"""
import os

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import importlib
import multiprocessing as mp
import time

import numpy as np

sim = importlib.import_module("06_gradient_estimator")
mix = importlib.import_module("26_mixup_split")

NREAL, PS = 60, (80, 160)
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "29_mixup_localization.txt")


def mixup_located(X, z):
    o = np.argsort(z, kind="stable")
    best, bz = np.inf, np.nan
    for f in mix.FRACS:
        cut = max(3, int(round(f * len(z))))
        T = mix.stat_at(mix.corr_profiles(X[:, o[:cut]]),
                        mix.corr_profiles(X[:, o[cut:]]))
        if T < best:
            best, bz = T, float((z[o[cut - 1]] + z[o[cut]]) / 2)
    return bz


def job(args):
    P, i = args
    X, z = sim.community(P, "THRESH", np.random.default_rng(91000 + 7 * P + i))
    W, stride = sim.design(P)
    z_curve = sim.ESTIMATORS[sim.CHOSEN](*sim.window_curves(X, z, W, stride))[0]
    return P, mixup_located(X, z), z_curve


if __name__ == "__main__":
    t0 = time.time()
    with mp.get_context("spawn").Pool(int(os.environ.get("NPROC", 8))) as pool:
        res = pool.map(job, [(P, i) for P in PS for i in range(NREAL)], chunksize=1)
    L = [f"recovery of z* = {sim.ZSTAR} on merge gradients, {NREAL} realizations",
         "", f"{'P':>5} {'estimator':>12} {'bias':>8} {'RMSE':>8} {'sd':>8}"]
    for P in PS:
        for k, idx in (("mixup", 1), ("curve", 2)):
            v = np.array([r[idx] for r in res if r[0] == P])
            L.append(f"{P:5d} {k:>12} {v.mean() - sim.ZSTAR:+8.3f} "
                     f"{np.sqrt(np.mean((v - sim.ZSTAR) ** 2)):8.3f} {v.std(ddof=1):8.3f}")
    L += ["", f"[total {time.time() - t0:.0f}s]"]
    open(OUT, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
