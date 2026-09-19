"""Shared Euler-characteristic machinery for the ECS-on-omics probes.

The substrate is the FLAG (Vietoris-Rips) complex of a thresholded similarity
graph, not the graph itself: a graph is 1-dimensional, so chi = V - E is
exactly the edge count and carries no topological information.

Sweeps are indexed by EDGE COUNT rather than radius, which holds edge density
fixed by construction -- so any chi difference between groups is structural
rather than a difference in the correlation distribution.
"""
import numpy as np
import gudhi as gd
from scipy.spatial.distance import squareform

MAXDIM = 3


def corr_distance(X):
    """Correlation distance between rows of X.  Bounded in [0, 2], so it is
    comparable across samples with different absolute expression scales."""
    Xc = X - X.mean(1, keepdims=True)
    Xc /= (np.linalg.norm(Xc, axis=1, keepdims=True) + 1e-12)
    return 1.0 - np.clip(Xc @ Xc.T, -1.0, 1.0)


def chi_by_edge_count(D, m_grid, max_dim=MAXDIM):
    """chi of the flag complex at each prescribed number of edges."""
    d = np.sort(squareform(D, checks=False))
    r_at_m = d[np.minimum(np.asarray(m_grid), len(d) - 1)]
    st = gd.RipsComplex(distance_matrix=D,
                        max_edge_length=float(r_at_m[-1])
                        ).create_simplex_tree(max_dimension=max_dim)
    filt = np.fromiter((v for _, v in st.get_filtration()), float)
    sign = np.fromiter(((-1.0) ** (len(s) - 1) for s, _ in st.get_filtration()), float)
    o = np.argsort(filt, kind="stable")
    cum = np.cumsum(sign[o])
    i = np.searchsorted(filt[o], r_at_m, side="right") - 1
    return np.where(i >= 0, cum[np.clip(i, 0, None)], 0.0)


def chi_by_radius(D, r_grid, max_dim=MAXDIM):
    """chi of the flag complex on a fixed radius grid (density NOT controlled)."""
    st = gd.RipsComplex(distance_matrix=D, max_edge_length=float(r_grid[-1])
                        ).create_simplex_tree(max_dimension=max_dim)
    filt = np.fromiter((v for _, v in st.get_filtration()), float)
    sign = np.fromiter(((-1.0) ** (len(s) - 1) for s, _ in st.get_filtration()), float)
    o = np.argsort(filt, kind="stable")
    cum = np.cumsum(sign[o])
    i = np.searchsorted(filt[o], r_grid, side="right") - 1
    return np.where(i >= 0, cum[np.clip(i, 0, None)], 0.0)


def edge_counts(D, r_grid):
    d = np.sort(squareform(D, checks=False))
    return np.searchsorted(d, r_grid, side="right").astype(float)


def cohen_d(a, b):
    """Standardized mean difference per column, between two replicate stacks."""
    pooled = np.sqrt((a.var(0, ddof=1) + b.var(0, ddof=1)) / 2) + 1e-12
    return (a.mean(0) - b.mean(0)) / pooled
