"""Decisive test.  Sweep by EDGE COUNT, not radius -- density is then identical
across accessions by construction, so any chi difference is pure structure.

Two geometries with deliberately similar local statistics:
  A  isolated modules            (many components, no loop)
  B  modules arranged in a cycle (same module count/size/spread, one big loop)
If chi separates A from B at matched density, topology carries information the
correlation distribution cannot.  If it does not, the concept is in trouble."""
import numpy as np, gudhi as gd
from scipy.spatial.distance import squareform

G, NMOD, MAXDIM, NACC = 120, 8, 3, 15
SPREAD, RADIUS = 0.55, 3.0

def cloud(loop, rng):
    """NMOD gaussian blobs, same size and spread.  loop=False: blobs at random
    positions.  loop=True: blob centres equally spaced on a circle."""
    per = G // NMOD
    if loop:
        th = np.linspace(0, 2*np.pi, NMOD, endpoint=False)
        ctr = np.c_[RADIUS*np.cos(th), RADIUS*np.sin(th), np.zeros(NMOD)]
    else:
        ctr = rng.normal(scale=RADIUS, size=(NMOD, 3))
    return np.vstack([c + rng.normal(scale=SPREAD, size=(per, 3)) for c in ctr])

def chi_at_edge_counts(X, m_grid, max_dim):
    """chi of the flag complex, indexed by number of edges present."""
    from scipy.spatial.distance import pdist
    d = np.sort(pdist(X))
    r_at_m = d[np.minimum(m_grid, len(d)-1)]           # radius giving m edges
    st = gd.RipsComplex(points=X, max_edge_length=float(r_at_m[-1])
        ).create_simplex_tree(max_dimension=max_dim)
    f = np.fromiter((v for _, v in st.get_filtration()), float)
    s = np.fromiter(((-1.0)**(len(x)-1) for x, _ in st.get_filtration()), float)
    o = np.argsort(f, kind="stable"); f, s = f[o], s[o]
    cum = np.cumsum(s)
    i = np.searchsorted(f, r_at_m, side="right") - 1
    return np.where(i >= 0, cum[np.clip(i, 0, None)], 0.0)

m_grid = np.unique(np.linspace(60, 2200, 18).astype(int))

out = {}
for label, loop in (("isolated", False), ("loop", True)):
    out[label] = np.array([
        chi_at_edge_counts(cloud(loop, np.random.default_rng(4200 + a + (0 if loop else 700))),
                           m_grid, MAXDIM)
        for a in range(NACC)])

a, b = out['isolated'], out['loop']
pooled = np.sqrt((a.var(0, ddof=1) + b.var(0, ddof=1)) / 2) + 1e-12
d = (a.mean(0) - b.mean(0)) / pooled

print(f"genes={G} modules={NMOD} max_dim={MAXDIM} accessions={NACC}/group")
print("density held FIXED by construction (sweep indexed by edge count)\n")
print(f"{'edges':>7} {'chi isolated':>13} {'chi loop':>10} {'Cohen d':>9}")
for i in range(0, len(m_grid), 2):
    print(f"{m_grid[i]:7d} {a[:,i].mean():13.1f} {b[:,i].mean():10.1f} {d[i]:9.2f}")
print(f"\nmax |Cohen d| at matched density = {np.abs(d).max():.2f} "
      f"(at {m_grid[np.abs(d).argmax()]} edges)")
