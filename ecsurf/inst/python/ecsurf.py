"""ecsurf: connectance-indexed Euler characteristic surfaces along gradients.

Library-grade core of the pipeline in ../experiments, with the estimator
exactly as pre-registered in the paper (sec. 2.4): windows of P/4 plots
sliding by W/3, chi of the flag complex at 12 edge counts from 2S/3 to half
of all pairs, pooled-|Cohen d| split statistic, plot-position permutation
test, within-window bootstrap CI.  Optional `groups` implements the stacked
repeat-survey design: columns sharing a group id (a site's years) carry one
gradient value and are permuted and windowed together.

API:  estimate_threshold(X, z, nperm=999, nboot=999, groups=None, seed=0)
      X: taxa x samples array (rows = taxa); z: gradient value per sample.
"""
import numpy as np
import gudhi as gd
from scipy.spatial.distance import squareform

MAXDIM = 3


def _corr_distance(X):
    Xc = X - X.mean(1, keepdims=True)
    Xc /= (np.linalg.norm(Xc, axis=1, keepdims=True) + 1e-12)
    return 1.0 - np.clip(Xc @ Xc.T, -1.0, 1.0)


def _chi_by_edge_count(D, m_grid):
    d = np.sort(squareform(D, checks=False))
    r_at_m = d[np.minimum(np.asarray(m_grid), len(d) - 1)]
    st = gd.RipsComplex(distance_matrix=D, max_edge_length=float(r_at_m[-1])
                        ).create_simplex_tree(max_dimension=MAXDIM)
    f = np.fromiter((v for _, v in st.get_filtration()), float)
    s = np.fromiter(((-1.0) ** (len(x) - 1) for x, _ in st.get_filtration()),
                    float)
    o = np.argsort(f, kind="stable")
    cum = np.cumsum(s[o])
    i = np.searchsorted(f[o], r_at_m, side="right") - 1
    return np.where(i >= 0, cum[np.clip(i, 0, None)], 0.0)


def _default_m_grid(S):
    return np.unique(np.linspace(int(np.ceil(2 * S / 3)),
                                 (S * (S - 1) // 2) // 2, 12).astype(int))


def _windows(z_units, W, stride):
    order = np.argsort(z_units, kind="stable")
    return [order[s:s + W] for s in range(0, len(order) - W + 1, stride)]


def _curves(X, z, groups, W, stride, m_grid):
    units = np.unique(groups)
    zu = np.array([z[groups == u][0] for u in units])
    C, zc = [], []
    for w in _windows(zu, W, stride):
        cols = np.isin(groups, units[w])
        C.append(_chi_by_edge_count(_corr_distance(X[:, cols]), m_grid))
        zc.append(float(zu[w].mean()))
    return np.array(C), np.array(zc)


def _split_stat(C, zc):
    best, bz = -np.inf, np.nan
    for s in range(2, len(zc) - 1):
        a, b = C[:s], C[s:]
        pooled = np.sqrt((a.var(0, ddof=1) + b.var(0, ddof=1)) / 2) + 1e-12
        d = float(np.abs((a.mean(0) - b.mean(0)) / pooled).mean())
        if d > best:
            best, bz = d, (zc[s - 1] + zc[s]) / 2
    return bz, best


def estimate_threshold(X, z, nperm=999, nboot=999, groups=None, seed=0):
    X = np.asarray(X, float)
    z = np.asarray(z, float)
    if groups is None:
        groups = np.arange(X.shape[1])
    groups = np.asarray(groups)
    units = np.unique(groups)
    P = len(units)
    W = max(6, int(round(P / 4)))
    stride = max(1, int(round(W / 3)))
    m_grid = _default_m_grid(X.shape[0])

    zhat, obs = _split_stat(*_curves(X, z, groups, W, stride, m_grid))

    rng = np.random.default_rng(seed)
    zu = np.array([z[groups == u][0] for u in units])
    hits = 0
    for _ in range(nperm):
        zp_units = rng.permutation(zu)
        zp = np.array([zp_units[np.searchsorted(units, g)] for g in groups])
        hits += _split_stat(*_curves(X, zp, groups, W, stride, m_grid))[1] >= obs
    p = (1 + hits) / (1 + nperm)

    if nboot == 0:
        return {"zhat": float(zhat), "stat": float(obs), "p": float(p),
                "ci": [float("nan")] * 2, "windows": int(P), "W": int(W)}
    boot = []
    for _ in range(nboot):
        C, zc = [], []
        for w in _windows(zu, W, stride):
            pick = rng.choice(units[w], W, replace=True)
            cols = np.concatenate([np.flatnonzero(groups == u) for u in pick])
            C.append(_chi_by_edge_count(_corr_distance(X[:, cols]), m_grid))
            zc.append(float(zu[w].mean()))
        boot.append(_split_stat(np.array(C), np.array(zc))[0])
    lo, hi = np.percentile(boot, [2.5, 97.5])

    return {"zhat": float(zhat), "stat": float(obs), "p": float(p),
            "ci": [float(lo), float(hi)], "windows": int(P), "W": int(W)}
