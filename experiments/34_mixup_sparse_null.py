"""Is the interaction test calibrated on VERY sparse matrices at survey scale?

30_mixup_largeP_null.py checked P = 400-1600 at 40% occupancy.  Several of
the detecting field matrices are far sparser than that -- D6 is 92% zeros,
D4 is 90% -- and a test that over-rejects on near-empty correlation
matrices would manufacture exactly those detections.  This measures the
false-positive rate on NULL gradients at 10% and 20% occupancy, at the unit
counts those datasets have.

A second cell varies what the field data also varies and the earlier
simulations did not: the number of COLUMNS per unit.  D6 carries 29
plot-years per unit, so its correlations rest on far more data than its
zero fraction suggests, while D4 carries about one.  `cols` replicates each
unit's latent draw into that many observed columns.
Output: results/34_mixup_sparse_null.txt
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
#        units, occupancy, columns per unit
CELLS = [(700, "count:0.1", 1),
         (700, "count:0.2", 1),
         (200, "count:0.1", 29),
         (200, "count:0.1", 1),
         (1100, "count:0.3", 1)]
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "34_mixup_sparse_null.txt")


def job(args):
    P, obs, cols, i = args
    rng = np.random.default_rng(96000 + 17 * P + 3 * cols + i)
    z = np.sort(rng.uniform(size=P))
    a = np.full(P, lim.BASE["a"] / 2)            # constant: no change along z
    L = lim.latent(P, 60, a, z, rng)             # taxa x units
    if cols > 1:                                 # repeat draws within a unit
        L = np.repeat(L, cols, axis=1)
        zc = np.repeat(z, cols)
    else:
        zc = z
    X = lim.observe(L, obs, rng)
    T = mix.scan(X, zc)
    # the null permutes UNITS, so a unit's columns travel together
    o = np.arange(P)
    hits = 0
    for _ in range(NPERM):
        pz = np.repeat(rng.permutation(z), cols) if cols > 1 else rng.permutation(z)
        hits += mix.scan(X, pz) <= T
    return P, obs, cols, (1 + hits) / (1 + NPERM)


if __name__ == "__main__":
    t0 = time.time()
    jobs = [(P, o, c, i) for P, o, c in CELLS for i in range(NREAL)]
    with mp.get_context("spawn").Pool(int(os.environ.get("NPROC", 16))) as pool:
        res = pool.map(job, jobs, chunksize=1)
    L = [f"false-positive rate of the interaction test on NULL gradients at "
         f"field sparsity, {NREAL} realizations, {NPERM} permutations, "
         f"alpha 0.05",
         "(no threshold, no association with z; nominal rate is 0.05)", "",
         f"{'units':>6} {'occupancy':>12} {'cols/unit':>10} {'zeros':>7} "
         f"{'FPR':>6} {'median p':>10}"]
    for P, o, c in CELLS:
        v = np.array([r[3] for r in res if r[0] == P and r[1] == o and r[2] == c])
        L.append(f"{P:6d} {o:>12} {c:10d} {'':>7} "
                 f"{np.mean(v <= .05):6.2f} {np.median(v):10.3f}")
    L += ["", f"[total {time.time() - t0:.0f}s]"]
    open(OUT, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
