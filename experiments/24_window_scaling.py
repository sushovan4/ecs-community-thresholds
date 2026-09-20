"""Does the window rule W = P/4 break down on large surveys?

The rule fixes the NUMBER of windows at about ten regardless of P, so a
survey with 1000 units gets windows spanning a quarter of the gradient each.
Extra units then buy precision within a window but never resolution along
the gradient, and a threshold is averaged away.  This is untested above
P = 160, yet four of the spatial datasets have 500-1900 units.

Here W is varied independently of P on the merge gradient of
06_gradient_estimator.py:
  P in {160, 400, 800};  W in {20, 40, 80, P/4}  (stride W/3 throughout)
reporting rejection at alpha=0.05 (THRESH), false positives (NULL), and RMSE
of z-hat.  A capped window should gain power and localization at large P if
the diagnosis is right; if it does not, the P/4 rule stands.
Output: results/24_window_scaling.txt
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
CELLS = [(P, W) for P in (160, 400, 800)
         for W in sorted({20, 40, 80, max(6, int(round(P / 4)))})]
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "24_window_scaling.txt")


def job(args):
    P, W, scen, i = args
    stride = max(1, int(round(W / 3)))
    est = sim.ESTIMATORS[sim.CHOSEN]
    X, z = sim.community(P, scen, np.random.default_rng(61000 + 13 * P + i))
    zhat, obs = est(*sim.window_curves(X, z, W, stride))
    rng = np.random.default_rng(62000 + 13 * P + i)
    hits = 0
    for _ in range(NPERM):
        hits += est(*sim.window_curves(X, rng.permutation(z), W, stride))[1] >= obs
    return P, W, scen, zhat, (1 + hits) / (1 + NPERM)


if __name__ == "__main__":
    t0 = time.time()
    jobs = [(P, W, s, i) for P, W in CELLS for s in ("THRESH", "NULL")
            for i in range(NREAL)]
    with mp.get_context("spawn").Pool(int(os.environ.get("NPROC", 4))) as pool:
        res = pool.map(job, jobs, chunksize=1)
    L = [f"window width vs survey size, {NREAL} realizations, {NPERM} perms, "
         f"alpha 0.05; merge gradient, true z* = {sim.ZSTAR}",
         "(MC SE of a rate <= %.2f)" % (0.5 / np.sqrt(NREAL)), "",
         f"{'P':>5} {'W':>5} {'windows':>8} {'power':>7} {'FPR':>6} {'RMSE':>7}  rule"]
    for P, W in CELLS:
        th = [(zh, p) for p_, w_, s, zh, p in res
              if (p_, w_, s) == (P, W, "THRESH")]
        nu = [p for p_, w_, s, zh, p in res if (p_, w_, s) == (P, W, "NULL")]
        nwin = len(range(0, P - W + 1, max(1, int(round(W / 3)))))
        zh = np.array([t[0] for t in th])
        tag = "<- P/4 rule" if W == max(6, int(round(P / 4))) else ""
        L.append(f"{P:5d} {W:5d} {nwin:8d} "
                 f"{np.mean([t[1] <= .05 for t in th]):7.2f} "
                 f"{np.mean([p <= .05 for p in nu]):6.2f} "
                 f"{np.sqrt(np.mean((zh - sim.ZSTAR) ** 2)):7.3f}  {tag}")
    L.append("")
    L.append(f"[total {time.time() - t0:.0f}s]")
    open(OUT, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
