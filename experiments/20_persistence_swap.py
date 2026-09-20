"""Persistence summaries through the same estimator (extends Table 2).

"Why not persistence?"  The Euler curve is the alternating sum of Betti
curves; Betti curves were already compared (08_summary_swap.py). Here the
stronger persistence representations are added, all on the SAME edge-count
(rank) filtration so connectance stays controlled:
  pl1_H0   first persistence landscape of H0, evaluated on the edge-count grid
  pl1_H1   first persistence landscape of H1
  totpers  total persistence curve, sum over H0 and H1 bars of the length
           of each bar's part below m
  chi      Euler characteristic, the paper's summary (reference)
Scenarios and statistic as Table 2 (08_summary_swap.py: P=160, REWIRE,
STRENGTH, TRI_NULL; merge THRESH/NULL; floored pooled split).  Statistics
over 100 realizations; permutation rejection over 30 realizations with 99
permutations.  Output: results/20_persistence_swap.txt
"""
import os

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import importlib
import multiprocessing as mp
import time

import gudhi as gd
import numpy as np
from scipy.spatial.distance import squareform

from ec import corr_distance

swap = importlib.import_module("08_summary_swap")
sim = importlib.import_module("06_gradient_estimator")
P, NSTAT, NREJ, NPERM = 160, 100, 30, 99
KEYS = ("chi", "pl1_H0", "pl1_H1", "totpers")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "20_persistence_swap.txt")


def rank_matrix(D):
    """Replace distances by their rank among all pairs (1 = closest), so the
    Rips filtration value of every simplex is an edge count."""
    S = D.shape[0]
    v = squareform(D, checks=False)
    r = np.empty_like(v)
    r[np.argsort(v, kind="stable")] = np.arange(1, len(v) + 1)
    return squareform(r)


def summaries(X):
    mg = sim.M_GRID
    D = corr_distance(X)
    R = rank_matrix(D)
    top = float(mg[-1])
    st = gd.RipsComplex(distance_matrix=R, max_edge_length=top
                        ).create_simplex_tree(max_dimension=3)
    f = np.fromiter((v for _, v in st.get_filtration()), float)
    sg = np.fromiter(((-1.0) ** (len(s) - 1) for s, _ in st.get_filtration()),
                     float)
    o = np.argsort(f, kind="stable")
    cum = np.cumsum(sg[o])
    i = np.searchsorted(f[o], mg, side="right") - 1
    chi = np.where(i >= 0, cum[np.clip(i, 0, None)], 0.0)
    st2 = gd.RipsComplex(distance_matrix=R, max_edge_length=top
                         ).create_simplex_tree(max_dimension=2)
    st2.compute_persistence(homology_coeff_field=2)
    out = {"chi": chi}
    tp = np.zeros(len(mg))
    for dim in (0, 1):
        iv = st2.persistence_intervals_in_dimension(dim)
        if len(iv) == 0:
            out[f"pl1_H{dim}"] = np.zeros(len(mg))
            continue
        b = iv[:, 0]
        d = np.where(np.isinf(iv[:, 1]), top + 1, iv[:, 1])
        tri = np.minimum(mg[None, :] - b[:, None], d[:, None] - mg[None, :])
        out[f"pl1_H{dim}"] = np.maximum(tri, 0).max(0)
        tp += (np.minimum(d[:, None], mg[None, :])
               - np.minimum(b[:, None], mg[None, :])).sum(0)
    out["totpers"] = tp
    return out


def curves(X, z):
    W, stride = sim.design(P)
    o = np.argsort(z, kind="stable")
    starts = range(0, len(z) - W + 1, stride)
    per = [summaries(X[:, o[s:s + W]]) for s in starts]
    zc = np.array([z[o[s:s + W]].mean() for s in starts])
    return {k: np.array([p[k] for p in per]) for k in KEYS}, zc


def stat_job(args):
    sc, i = args
    X, z = swap.community(P, sc, np.random.default_rng(8600 + i))
    C, zc = curves(X, z)
    return sc, {k: swap.est_split_floored(C[k], zc) for k in KEYS}


def rej_job(args):
    sc, i = args
    X, z = swap.community(P, sc, np.random.default_rng(9600 + i))
    C, zc = curves(X, z)
    obs = {k: swap.est_split_floored(C[k], zc)[1] for k in KEYS}
    hits = dict.fromkeys(KEYS, 0)
    r = np.random.default_rng(9900 + i)
    for _ in range(NPERM):
        Cb, zb = curves(X, r.permutation(z))
        for k in KEYS:
            hits[k] += swap.est_split_floored(Cb[k], zb)[1] >= obs[k]
    return sc, {k: (1 + hits[k]) / (1 + NPERM) for k in KEYS}


if __name__ == "__main__":
    t0 = time.time()
    ctx = mp.get_context("spawn")
    scen_s = ("REWIRE", "STRENGTH", "TRI_NULL", "THRESH", "NULL")
    scen_r = ("REWIRE", "STRENGTH", "TRI_NULL", "THRESH")
    with ctx.Pool(int(os.environ.get("NPROC", 4))) as pool:
        st = pool.map(stat_job, [(s, i) for s in scen_s for i in range(NSTAT)],
                      chunksize=1)
        rj = pool.map(rej_job, [(s, i) for s in scen_r for i in range(NREJ)],
                      chunksize=1)
    L = [f"persistence summaries vs chi, P={P}, edge-count (rank) filtration; "
         f"stats over {NSTAT} realizations, rejection over {NREJ} with {NPERM} perms", "",
         "mean pooled split statistic"]
    L.append(f"{'summary':>8}" + "".join(f"{s:>10}" for s in scen_s)
             + f"{'RMSE(THRESH)':>14}")
    for k in KEYS:
        row = f"{k:>8}"
        for s in scen_s:
            row += f"{np.mean([v[k][1] for sc, v in st if sc == s]):10.2f}"
        zh = np.array([v[k][0] for sc, v in st if sc == "THRESH"])
        row += f"{np.sqrt(np.mean((zh - sim.ZSTAR) ** 2)):14.3f}"
        L.append(row)
    L += ["", "permutation rejection rate at alpha = 0.05"]
    L.append(f"{'summary':>8}" + "".join(f"{s:>10}" for s in scen_r))
    for k in KEYS:
        L.append(f"{k:>8}" + "".join(
            f"{np.mean([v[k] <= .05 for sc, v in rj if sc == s]):10.2f}"
            for s in scen_r))
    L += ["", f"[total {time.time() - t0:.0f}s]"]
    open(OUT, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
