"""Scale-free rerun.  Correlation distance (bounded, comparable across
accessions).  Decisive question: does chi carry signal BEYOND edge density?
Compare observed chi of the flag complex against the Erdos-Renyi prediction at
MATCHED edge count -- the density-only null."""
import numpy as np, gudhi as gd
from math import comb
from scipy.spatial.distance import squareform, pdist

G, S, NMOD, MAXDIM, NACC = 120, 6, 8, 3, 15

def accession(tolerant, rng):
    mod  = rng.integers(0, NMOD, size=G)
    prof = rng.normal(size=(NMOD, S))
    X    = prof[mod].copy()
    decay = np.linspace(0.0, 0.0 if tolerant else 1.6, S)
    return X + rng.normal(scale=0.35 + decay, size=(G, S))

def corr_dist(X):
    Xc = X - X.mean(1, keepdims=True)
    Xc /= (np.linalg.norm(Xc, axis=1, keepdims=True) + 1e-12)
    C = np.clip(Xc @ Xc.T, -1, 1)
    return 1.0 - C                      # in [0, 2], scale-free

def chi_flag(D, r_grid, max_dim):
    st = gd.RipsComplex(distance_matrix=D, max_edge_length=float(r_grid[-1])
        ).create_simplex_tree(max_dimension=max_dim)
    f = np.fromiter((v for _, v in st.get_filtration()), float)
    s = np.fromiter(((-1.0)**(len(x)-1) for x, _ in st.get_filtration()), float)
    o = np.argsort(f, kind="stable"); f, s = f[o], s[o]
    cum = np.cumsum(s)
    i = np.searchsorted(f, r_grid, side="right") - 1
    return np.where(i >= 0, cum[np.clip(i, 0, None)], 0.0)

def n_edges(D, r_grid):
    d = np.sort(squareform(D, checks=False))
    return np.searchsorted(d, r_grid, side="right").astype(float)

def chi_er(n, m_grid, max_dim):
    """E[chi] for Erdos-Renyi flag complex at matched edge count."""
    out = np.empty(len(m_grid))
    for i, m in enumerate(m_grid):
        p = min(1.0, 2.0 * m / (n * (n - 1)))
        out[i] = sum((-1.0)**k * comb(n, k+1) * p**comb(k+1, 2)
                     for k in range(max_dim + 1))
    return out

r_grid = np.linspace(0.05, 0.55, 16)

feat = {}
for label, tol in (("tolerant", True), ("sensitive", False)):
    obs, dens, dev = [], [], []
    for a in range(NACC):
        D = corr_dist(accession(tol, np.random.default_rng(7000 + a + (0 if tol else 900))))
        c = chi_flag(D, r_grid, MAXDIM); m = n_edges(D, r_grid)
        obs.append(c); dens.append(m); dev.append(c - chi_er(G, m, MAXDIM))
    feat[label] = tuple(np.array(z) for z in (obs, dens, dev))

def sep(a, b):
    """Standardized mean difference, per radius."""
    s = np.sqrt((a.var(0, ddof=1) + b.var(0, ddof=1)) / 2) + 1e-12
    return (a.mean(0) - b.mean(0)) / s

d_obs  = sep(feat['tolerant'][0], feat['sensitive'][0])
d_dens = sep(feat['tolerant'][1], feat['sensitive'][1])
d_dev  = sep(feat['tolerant'][2], feat['sensitive'][2])

print(f"genes={G} conditions={S} modules={NMOD} max_dim={MAXDIM} "
      f"accessions={NACC}/group\ncorrelation distance, common radius grid\n")
print(f"{'r':>6} {'|d| chi_obs':>12} {'|d| density':>12} {'|d| chi-ER dev':>15}")
for i in range(0, len(r_grid), 2):
    print(f"{r_grid[i]:6.3f} {abs(d_obs[i]):12.2f} {abs(d_dens[i]):12.2f} "
          f"{abs(d_dev[i]):15.2f}")
print(f"\nmax |Cohen d|  chi_obs={np.abs(d_obs).max():.2f}   "
      f"density={np.abs(d_dens).max():.2f}   "
      f"chi-minus-ER={np.abs(d_dev).max():.2f}")
