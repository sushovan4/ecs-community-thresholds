"""The estimator is summary-agnostic; when is topology NECESSARY?

Two reviewer questions this experiment answers with numbers:

  (a) "Why not just report the correlation distribution?" -- now asked in the
      LOCALISATION setting, where 02_synthetic_matched_density.py answered it
      only for two-group separation.
  (b) "Why the Euler characteristic and not Betti curves / modularity / any
      other summary?" -- the machinery of 06 consumes any curve-valued
      summary, so run them all through the SAME pre-registered estimator.

Scenario REWIRE is the experiment-02 move ported to a gradient.  Modules
share cross-module latent factors along a fixed 2-regular "module graph":
below z* the graph is two triangles (0-1-2, 3-4-5); above z* it is one
hexagon (0-1-2-3-4-5-0).  Every module keeps exactly two cross-links of
identical weight in both configurations, so the degree structure and hence
the expected pairwise-correlation DISTRIBUTION are identical on both sides of
the threshold; only the arrangement -- two components against one ring --
changes.  REWIRE uses rho=0.75, a=0.9 (cross-module corr 0.34, above the
20-plot sampling-noise floor of 0.22); THRESH/NULL pass through to 06 with
its parameters unchanged.  Windows that straddle z* mix the two configurations and dilute the
rewired links, so a little transition signal leaks into distribution
summaries even here; the smoke runs showed that leak clearly, and the design
answers it with a fourth scenario completing a 2x2:

  STRENGTH  the module graph stays two triangles everywhere, but the
            cross-link weight halves above z* -- a pure strength change.
            Distribution summaries see it directly; at matched edge counts
            the edge RANKING is preserved, so the flag complex is unchanged
            and shape summaries should be blind.

Scenario THRESH (the 06 merge, where distribution and shape move together) is
kept as the reference in which everything detects.  TRI_NULL is triangles
everywhere at full strength.

Summaries, all evaluated at the same matched edge counts and all fed to the
pre-registered split_mean estimator of 06:
  chi      Euler characteristic of the flag complex (the paper's summary)
  betti0   number of components alive at r(m)
  betti1   number of independent loops alive at r(m)
  modul    greedy modularity of the m-edge graph (network-ecology incumbent)
  rq       r(m), the distance quantile reaching m edges  -- distribution only
  meand    mean pairwise distance in the window          -- distribution only

Subcommands: recovery (stats + z-hat recovery, ~10 min) · rejection
(permutation rejection rates on REWIRE and NULL, background, ~1.5 h).
"""
import importlib
import sys
import time

import gudhi as gd
import networkx as nx
import numpy as np
from scipy.spatial.distance import squareform

from ec import corr_distance

sim = importlib.import_module("06_gradient_estimator")


def est_split_floored(C, zc, floor=0.05):
    """split_mean with columns standardised by their between-window deviation
    and the pooled within-group deviation floored at 5% of it.  Needed for the
    cross-summary comparison because the integer-valued Betti curves hit zero
    within-group variance at ties, where a bare Cohen d degenerates; for chi
    the floor never binds (checked against est_split_mean in the results log).
    """
    sc = C.std(0, ddof=1) + 1e-12
    Y = (C - C.mean(0)) / sc
    best, bz = -np.inf, np.nan
    for s in range(2, len(zc) - 1):
        a, b = Y[:s], Y[s:]
        pooled = np.sqrt((a.var(0, ddof=1) + b.var(0, ddof=1)) / 2)
        d = np.abs(a.mean(0) - b.mean(0)) / np.maximum(pooled, floor)
        if d.mean() > best:
            best, bz = float(d.mean()), (zc[s - 1] + zc[s]) / 2
    return bz, best

S, NMOD, PER = 60, 6, 10
ZSTAR = sim.ZSTAR
# REWIRE needs its cross-module links above the sampling-noise floor of a
# 20-plot window (corr sd ~ 1/sqrt(20) = 0.22): with rho*a/2 = 0.34 the
# module graph is recoverable; at 06's parameters (0.19) it is not, which the
# population-limit check in the results log quantifies.
RHO, A_HI = 0.75, 0.9
M_GRID = sim.M_GRID
TRI = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)]
HEX = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]


def community(P, scenario, rng):
    """REWIRE / TRI_NULL by module-graph rewiring; THRESH / NULL from 06."""
    if scenario in ("THRESH", "NULL"):
        return sim.community(P, scenario, rng)
    z = np.sort(rng.uniform(size=P))
    mod = np.repeat(np.arange(NMOD), PER)
    u = rng.normal(size=(NMOD, P))
    g = {e: rng.normal(size=P) for e in sorted(set(TRI) | set(HEX))}
    f = np.empty((NMOD, P))
    for p in range(P):
        edges = HEX if (scenario == "REWIRE" and z[p] > ZSTAR) else TRI
        a = A_HI / 2 if (scenario == "STRENGTH" and z[p] > ZSTAR) else A_HI
        for i in range(NMOD):
            mine = [g[e][p] for e in edges if i in e]
            f[i, p] = (np.sqrt(1 - a) * u[i, p]
                       + np.sqrt(a / 2) * sum(mine))
    X = np.sqrt(RHO) * f[mod] + np.sqrt(1 - RHO) * rng.normal(size=(S, P))
    return X, z


def betti_curve(st, dim, r_grid):
    iv = st.persistence_intervals_in_dimension(dim)
    if len(iv) == 0:
        return np.zeros(len(r_grid))
    b, d = iv[:, 0], iv[:, 1]
    return np.array([((b <= r) & (r < d)).sum() for r in r_grid], float)


def window_summaries(X, z, W, stride):
    """All six summary curves per sliding window."""
    o = np.argsort(z)
    starts = range(0, len(z) - W + 1, stride)
    out = {k: [] for k in ("chi", "betti0", "betti1", "modul", "rq", "meand")}
    for s in starts:
        D = corr_distance(X[:, o[s:s + W]])
        d = np.sort(squareform(D, checks=False))
        r_at_m = d[np.minimum(np.asarray(M_GRID), len(d) - 1)]
        st = gd.RipsComplex(distance_matrix=D, max_edge_length=float(r_at_m[-1])
                            ).create_simplex_tree(max_dimension=3)
        filt = np.fromiter((v for _, v in st.get_filtration()), float)
        sign = np.fromiter(((-1.0) ** (len(x) - 1)
                            for x, _ in st.get_filtration()), float)
        oo = np.argsort(filt, kind="stable")
        cum = np.cumsum(sign[oo])
        i = np.searchsorted(filt[oo], r_at_m, side="right") - 1
        out["chi"].append(np.where(i >= 0, cum[np.clip(i, 0, None)], 0.0))
        st.compute_persistence(homology_coeff_field=2)
        out["betti0"].append(betti_curve(st, 0, r_at_m))
        out["betti1"].append(betti_curve(st, 1, r_at_m))
        iu = np.triu_indices(S, 1)
        order = np.argsort(D[iu])
        mods = []
        for m in M_GRID:
            G = nx.Graph()
            G.add_nodes_from(range(S))
            G.add_edges_from(zip(iu[0][order[:m]], iu[1][order[:m]]))
            comm = nx.community.greedy_modularity_communities(G)
            mods.append(nx.community.modularity(G, comm))
        out["modul"].append(np.array(mods))
        out["rq"].append(r_at_m)
        out["meand"].append(np.array([D[iu].mean()]))
    zc = np.array([z[o[s:s + W]].mean() for s in starts])
    return {k: np.array(v) for k, v in out.items()}, zc


def run_recovery(nreal=30, P=80):
    W, stride = sim.design(P)
    print(f"S={S} rho={RHO} a_hi={A_HI} z*={ZSTAR}  P={P}, windows {W} "
          f"stride {stride}, {nreal} realisations")
    print("estimator: pre-registered split_mean of 06, per summary\n")
    scens = ("REWIRE", "STRENGTH", "TRI_NULL", "THRESH", "NULL")
    keys = ("chi", "betti0", "betti1", "modul", "rq", "meand")
    stats = {sc: {k: [] for k in keys} for sc in scens}
    zhat = {sc: {k: [] for k in keys} for sc in scens}
    for sc in scens:
        for i in range(nreal):
            C, zc = window_summaries(
                *community(P, sc, np.random.default_rng(8600 + i)), W, stride)
            for k in keys:
                z_, s_ = est_split_floored(C[k], zc)
                stats[sc][k].append(s_)
                zhat[sc][k].append(z_)
    hdr = f"{'summary':>8}" + "".join(f"{sc:>10}" for sc in scens) \
        + f"{'RMSE(REWIRE)':>14}{'RMSE(THRESH)':>14}"
    print(hdr)
    for k in keys:
        row = f"{k:>8}"
        for sc in scens:
            row += f"{np.mean(stats[sc][k]):10.2f}"  # floored stat
        for sc in ("REWIRE", "THRESH"):
            z_ = np.array(zhat[sc][k])
            row += f"{np.sqrt(((z_ - ZSTAR) ** 2).mean()):14.3f}"
        print(row)
    print("\ncolumns are the mean split_mean statistic per scenario;"
          "\nREWIRE holds the correlation distribution fixed across z* by"
          "\nconstruction, so distribution summaries (rq, meand) should sit at"
          "\ntheir TRI_NULL level there while shape summaries separate.")


def _one_rejection(args):
    sc, i, P, W, stride = args
    keys = ("chi", "betti0", "betti1", "modul", "rq", "meand")
    X, z = community(P, sc, np.random.default_rng(9600 + i))
    C, zc = window_summaries(X, z, W, stride)
    obs = {k: est_split_floored(C[k], zc)[1] for k in keys}
    hits = {k: 0 for k in keys}
    r = np.random.default_rng(9900 + i)
    for _ in range(sim.NPERM):
        Cb, zcb = window_summaries(X, r.permutation(z), W, stride)
        for k in keys:
            hits[k] += est_split_floored(Cb[k], zcb)[1] >= obs[k]
    return sc, {k: (1 + hits[k]) / (1 + sim.NPERM) for k in keys}


def run_rejection(nreal=15, P=80):
    from multiprocessing import Pool
    W, stride = sim.design(P)
    keys = ("chi", "betti0", "betti1", "modul", "rq", "meand")
    scens = ("REWIRE", "STRENGTH", "TRI_NULL")
    print(f"permutation rejection at alpha=0.05, {sim.NPERM} shuffles, "
          f"{nreal} realisations, P={P}\n", flush=True)
    jobs = [(sc, i, P, W, stride) for sc in scens for i in range(nreal)]
    with Pool(10) as pool:
        results = pool.map(_one_rejection, jobs)
    print(f"{'summary':>8}" + "".join(f"{sc:>10}" for sc in scens))
    for k in keys:
        row = f"{k:>8}"
        for sc in scens:
            ps = np.array([r[1][k] for r in results if r[0] == sc])
            row += f"{np.mean(ps <= 0.05):10.2f}"
        print(row)


if __name__ == "__main__":
    t0 = time.time()
    what = sys.argv[1] if len(sys.argv) > 1 else "recovery"
    if what in ("recovery", "all"):
        run_recovery()
    if what in ("rejection", "all"):
        run_rejection()
    print(f"[total {time.time() - t0:.0f}s]")
