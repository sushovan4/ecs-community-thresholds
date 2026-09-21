"""Is the interaction test calibrated at the unit counts of the surveys?

Every spatial dataset (500-1900 units) returned p = 0.005, the floor.  The
test's false-positive rate has only been measured at P = 160.  Here it is
measured on NULL gradients (no threshold, no association) at P = 400, 800
and 1600, on Gaussian abundances and on sparse counts, with the same
15-quantile scan and 199-permutation null used on the data.
Output: results/30_mixup_largeP_null.txt
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

NPERM, NREAL = 199, 40
CELLS = [(P, o) for P in (400, 800, 1600) for o in ("gauss", "count:0.4")]
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "30_mixup_largeP_null.txt")


def job(args):
    P, obs, i = args
    rng = np.random.default_rng(95000 + 13 * P + i)
    z = np.sort(rng.uniform(size=P))
    a = np.full(P, lim.BASE["a"] / 2)          # constant: no change along z
    X = lim.observe(lim.latent(P, 60, a, z, rng), obs, rng)
    T = mix.scan(X, z)
    hits = sum(mix.scan(X, rng.permutation(z)) <= T for _ in range(NPERM))
    return P, obs, (1 + hits) / (1 + NPERM)


if __name__ == "__main__":
    t0 = time.time()
    jobs = [(P, o, i) for P, o in CELLS for i in range(NREAL)]
    with mp.get_context("spawn").Pool(int(os.environ.get("NPROC", 16))) as pool:
        res = pool.map(job, jobs, chunksize=1)
    L = [f"false-positive rate of the interaction test on NULL gradients, "
         f"{NREAL} realizations, {NPERM} permutations, alpha 0.05",
         "(no threshold and no association with z; nominal rate is 0.05)", "",
         f"{'P':>6} {'data':>12} {'FPR':>6} {'median p':>10}"]
    for P, o in CELLS:
        v = np.array([r[2] for r in res if r[0] == P and r[1] == o])
        L.append(f"{P:6d} {o:>12} {np.mean(v <= .05):6.2f} {np.median(v):10.3f}")
    L += ["", f"[total {time.time() - t0:.0f}s]"]
    open(OUT, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
