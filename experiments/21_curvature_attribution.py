"""Which taxa carry a detected structural change?  Curvature attribution.

For a flag complex K, chi(K) = sum_v kappa(v) with the local (Knill) curvature
    kappa(v) = sum_{k>=0} (-1)^k f_k(v) / (k + 1),
f_k(v) = number of k-simplices containing taxon v (discrete Gauss-Bonnet;
Knill 2011; exact here for the 3-skeleton used throughout).  So a change in
the Euler characteristic between the two sides of the estimated split
decomposes exactly into per-taxon contributions.

Applied to the Everglades long-term panel detection (17_everglades_lt.py):
same data, windows, grid, and split.  Per taxon: the grid-averaged Cohen's d
of kappa between the windows left and right of the detected split, with a
permutation p-value from the same PSU-level permutation (199), and
Benjamini-Hochberg q-values across taxa.  TITAN's pure & reliable indicator
taxa on the same PSU-mean matrix are listed for comparison
(17_everglades_lt_titan.R).  Output: results/21_curvature_attribution.txt
"""
import importlib
import os
import sys

import gudhi as gd
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "ecsurf", "inst", "python"))
import ecsurf as E  # noqa: E402

lt = importlib.import_module("17_everglades_lt")
NPERM = 199


def kappa(D, mg, maxdim=3):
    """per-taxon curvature at each edge-count grid point: S x len(mg)."""
    d = np.sort(D[np.triu_indices(len(D), 1)])
    r = d[np.minimum(mg, len(d) - 1)]
    st = gd.RipsComplex(distance_matrix=D, max_edge_length=float(r[-1])
                        ).create_simplex_tree(max_dimension=maxdim)
    delta = np.zeros((len(D), len(mg)))
    for simp, f in st.get_filtration():
        j = np.searchsorted(r, f, side="left")
        if j < len(mg):
            k = len(simp) - 1
            delta[simp, j] += (-1.0) ** k / (k + 1)
    return np.cumsum(delta, axis=1)


def window_kappa(X, z, groups, W, stride, mg):
    units = np.unique(groups)
    zu = np.array([z[groups == u][0] for u in units])
    K, zc = [], []
    for w in E._windows(zu, W, stride):
        cols = np.isin(groups, units[w])
        K.append(kappa(E._corr_distance(X[:, cols]), mg))
        zc.append(float(zu[w].mean()))
    return np.array(K), np.array(zc)          # windows x S x grid


def split_index(C, zc):
    best, bs = -np.inf, None
    for s in range(2, len(zc) - 1):
        a, b = C[:s], C[s:]
        pooled = np.sqrt((a.var(0, ddof=1) + b.var(0, ddof=1)) / 2) + 1e-12
        dd = np.abs((a.mean(0) - b.mean(0)) / pooled).mean()
        if dd > best:
            best, bs = dd, s
    return bs


def taxon_d(K, s):
    """grid-averaged change in mean curvature, right minus left.  These sum
    over taxa exactly to the change in chi (Gauss-Bonnet), so each is a share
    of the detected change; a standardized version is undefined when a
    taxon's curvature is constant on one side (two windows there)."""
    return (K[s:].mean(0) - K[:s].mean(0)).mean(1)


if __name__ == "__main__":
    dens, keep, z_psu, psu = lt.build()
    X = np.log1p(dens.values.T)
    z = z_psu.loc[psu].values
    groups = pd.factorize(psu)[0]
    P = len(np.unique(groups))
    W = max(6, int(round(P / 4)))
    stride = max(1, int(round(W / 3)))
    mg = E._default_m_grid(X.shape[0])
    K, zc = window_kappa(X, z, groups, W, stride, mg)
    chi = K.sum(1)                                           # windows x grid
    C_chk, _ = E._curves(X, z, groups, W, stride, mg)
    assert np.allclose(chi, C_chk), "Gauss-Bonnet identity failed"
    s = split_index(chi, zc)
    zhat = (zc[s - 1] + zc[s]) / 2
    d_obs = taxon_d(K, s)
    units = np.unique(groups)
    zu = np.array([z[groups == u][0] for u in units])
    rng = np.random.default_rng(3)
    hits = np.zeros(len(keep))
    for _ in range(NPERM):
        zp = rng.permutation(zu)[groups]
        Kb, zb = window_kappa(X, zp, groups, W, stride, mg)
        sb = split_index(Kb.sum(1), zb)
        hits += np.abs(taxon_d(Kb, sb)) >= np.abs(d_obs)
    p = (1 + hits) / (1 + NPERM)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    q = np.minimum(q, 1)
    tf = os.path.join(ROOT, "results", "21_titan_indicators.txt")
    titan = pd.read_csv(tf) if os.path.exists(tf) else None
    L = ["Curvature attribution of the Everglades long-term detection",
         f"Gauss-Bonnet check: sum of per-taxon curvature equals chi in all "
         f"{chi.size} window-grid cells (exact)",
         f"split at z-hat = {zhat:.1f} d since last dry (windows {s} | {len(zc) - s})",
         f"per-taxon change in curvature (right minus left, grid-averaged); "
         f"sum over taxa = change in chi = {d_obs.sum():+.2f}",
         f"share = taxon change / chi change; PSU permutation p ({NPERM}), BH q", "",
         f"{'taxon':>18} {'change':>9} {'share':>7} {'p':>6} {'q':>6}  TITAN"]
    for i in np.argsort(-np.abs(d_obs)):
        t = keep[i]
        tag = ""
        if titan is not None and t in set(titan.taxon):
            r = titan[titan.taxon == t].iloc[0]
            tag = f"{r.group} indicator, cp {r.cp:.0f} d"
        L.append(f"{t:>18} {d_obs[i]:+9.2f} {d_obs[i] / d_obs.sum():+7.2f} "
                 f"{p[i]:6.3f} {q[i]:6.3f}  {tag}")
    out = os.path.join(ROOT, "results", "21_curvature_attribution.txt")
    open(out, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
