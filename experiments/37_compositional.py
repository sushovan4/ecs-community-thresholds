"""Does closure manufacture detections?

Correlations among the parts of a constant-sum composition are negatively
biased by the constraint itself (Pearson 1897), which is why the microbiome
literature works in log-ratio coordinates.  The objection applies to any
method that reads a correlation matrix, including this one.

The argument that it does NOT by itself produce false detections is that the
permutation null absorbs it: permuting gradient positions leaves each unit's
composition untouched, so a closure bias that is the SAME everywhere along
the gradient appears identically in the observed statistic and in every
replicate.  What that argument does not cover is closure bias that VARIES
along the gradient -- if total abundance or richness changes with z, so does
the strength of the induced bias, and the two sides of a split differ for a
reason that has nothing to do with co-occurrence.

Cells, all with INDEPENDENT taxa (no co-occurrence structure anywhere, so
every rejection is a false positive):
  open        counts, no closure                       -- the control
  closed      counts divided by their unit total       -- homogeneous closure
  closed+tot  closure, with total abundance rising 8x along z
  closed+ric  closure, with richness rising along z
  clr         the closed+tot data in centred log-ratio coordinates

Output: results/37_compositional.txt
"""
import os

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import importlib
import multiprocessing as mp
import time

import numpy as np

mix = importlib.import_module("26_mixup_split")

NPERM, NREAL, P, S = 199, 40, 160, 60
CELLS = ("open", "closed", "closed+tot", "closed+ric", "clr")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "37_compositional.txt")


def clr(Y):
    """centred log ratio, zeros replaced multiplicatively at half the
    smallest positive value in the unit."""
    Y = np.asarray(Y, float).copy()
    for j in range(Y.shape[1]):
        col = Y[:, j]
        pos = col[col > 0]
        if len(pos) == 0:
            continue
        col[col == 0] = pos.min() / 2
        Y[:, j] = col / col.sum()
    L = np.log(Y)
    return L - L.mean(0, keepdims=True)


def make(cell, rng):
    """taxa x units of INDEPENDENT taxa, observed as the cell prescribes."""
    z = np.sort(rng.uniform(size=P))
    lam = np.exp(rng.normal(0.0, 0.8, size=(S, P)))      # independent taxa
    if cell == "closed+tot":
        lam = lam * (1.0 + 7.0 * z)[None, :]             # total rises 8x
    counts = rng.poisson(lam).astype(float)
    if cell == "closed+ric":
        # the number of taxa that can occur at all rises along the gradient
        for j in range(P):
            k = int(round(S * (0.35 + 0.6 * z[j])))
            off = rng.permutation(S)[k:]
            counts[off, j] = 0.0
    if cell == "open":
        return np.log1p(counts), z
    tot = counts.sum(0, keepdims=True)
    tot[tot == 0] = 1.0
    if cell == "clr":
        return clr(counts), z
    return np.log1p(counts / tot * 100.0), z             # percent-cover-like


def job(args):
    cell, i = args
    rng = np.random.default_rng(97000 + 101 * CELLS.index(cell) + i)
    X, z = make(cell, rng)
    T = mix.scan(X, z)
    hits = sum(mix.scan(X, rng.permutation(z)) <= T for _ in range(NPERM))
    return cell, (1 + hits) / (1 + NPERM)


if __name__ == "__main__":
    t0 = time.time()
    jobs = [(c, i) for c in CELLS for i in range(NREAL)]
    with mp.get_context("spawn").Pool(int(os.environ.get("NPROC", 10))) as pool:
        res = pool.map(job, jobs, chunksize=1)
    L = [f"Does closure manufacture detections?  {NREAL} realizations, "
         f"P = {P} units, S = {S} INDEPENDENT taxa, {NPERM} permutations, "
         f"alpha 0.05",
         "(no co-occurrence structure in any cell, so every rejection is a "
         "false positive; nominal rate 0.05)", "",
         f"{'cell':>12} {'FPR':>6} {'median p':>10}"]
    for c in CELLS:
        v = np.array([r[1] for r in res if r[0] == c])
        L.append(f"{c:>12} {np.mean(v <= .05):6.2f} {np.median(v):10.3f}")
    L += ["", f"[total {time.time() - t0:.0f}s]"]
    open(OUT, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
