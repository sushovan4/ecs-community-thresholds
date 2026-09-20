"""Export real curves from the completed probes as JSON, for the briefing page.
Nothing here is synthesised for presentation: every number is a run output
(the switchgrass d_by_tod values are transcribed from the 04 run log)."""
import os, json, glob
import numpy as np, pandas as pd
from math import comb
from ec import corr_distance, chi_by_edge_count, chi_by_radius, edge_counts, cohen_d
from scipy.spatial.distance import squareform, pdist

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA, OUT = os.path.join(ROOT, "data"), os.path.join(ROOT, "results")
rngs = np.random.default_rng
out = {}

# ---- A. synthetic modules: radius sweep, where density wins ----------------
G, S, NMOD, NACC = 120, 6, 8, 15
def accession(tol, rng):
    mod, prof = rng.integers(0, NMOD, G), rng.normal(size=(NMOD, S))
    X = prof[mod].copy()
    return X + rng.normal(scale=0.35 + np.linspace(0, 0 if tol else 1.6, S), size=(G, S))

r_grid = np.linspace(0.05, 0.55, 16)
def chi_er(n, m, md=3):
    p = min(1.0, 2.0*m/(n*(n-1)))
    return sum((-1.0)**k * comb(n, k+1) * p**comb(k+1, 2) for k in range(md+1))

rad = {}
for lab, tol in (("tolerant", True), ("sensitive", False)):
    ch, de, dv = [], [], []
    for a in range(NACC):
        D = corr_distance(accession(tol, rngs(7000+a+(0 if tol else 900))))
        c = chi_by_radius(D, r_grid); m = edge_counts(D, r_grid)
        ch.append(c); de.append(m); dv.append(c - np.array([chi_er(G, mi) for mi in m]))
    rad[lab] = tuple(np.array(z) for z in (ch, de, dv))
out["radius_sweep"] = {
    "r": r_grid.round(4).tolist(),
    "d_chi":     np.abs(cohen_d(rad['tolerant'][0], rad['sensitive'][0])).round(3).tolist(),
    "d_density": np.abs(cohen_d(rad['tolerant'][1], rad['sensitive'][1])).round(3).tolist(),
    "d_dev":     np.abs(cohen_d(rad['tolerant'][2], rad['sensitive'][2])).round(3).tolist()}

# ---- B. synthetic geometry: matched density, where topology wins -----------
GG, NM, SPREAD, RADIUS = 120, 8, 0.55, 3.0
def cloud(loop, rng):
    per = GG//NM
    if loop:
        th = np.linspace(0, 2*np.pi, NM, endpoint=False)
        ctr = np.c_[RADIUS*np.cos(th), RADIUS*np.sin(th), np.zeros(NM)]
    else:
        ctr = rng.normal(scale=RADIUS, size=(NM, 3))
    return np.vstack([c + rng.normal(scale=SPREAD, size=(per, 3)) for c in ctr])

m_grid = np.unique(np.linspace(60, 2200, 18).astype(int))
md = {}
for lab, loop in (("isolated", False), ("loop", True)):
    md[lab] = np.array([chi_by_edge_count(
        squareform(pdist(cloud(loop, rngs(4200+a+(0 if loop else 700))))), m_grid)
        for a in range(NACC)])
out["matched_density"] = {
    "m": m_grid.tolist(),
    "isolated": md['isolated'].mean(0).round(1).tolist(),
    "loop":     md['loop'].mean(0).round(1).tolist(),
    "d":        np.abs(cohen_d(md['isolated'], md['loop'])).round(3).tolist()}

# ---- C. switchgrass: the chi curves behind the diel result -----------------
txt = open(os.path.join(DATA, "gse57887_samples.txt")).read()
import re
meta = []
for b in txt.split("^SAMPLE = ")[1:]:
    g = lambda p: (re.search(p, b).group(1).strip() if re.search(p, b) else None)
    meta.append(dict(gsm=b.split("\n",1)[0].strip(), treatment=g(r"treatment:\s*(.+)"),
                     day=g(r"sampling day:\s*(.+)"), tod=g(r"sampling time:\s*(.+)")))
meta = pd.DataFrame(meta)
counts = pd.read_pickle(os.path.join(DATA, "gse57887_counts.pkl.gz"))
meta = meta[meta.gsm.isin(counts.columns)]
X = np.log2((counts/counts.sum(0)*1e6).values + 1.0)
X = X[np.argsort(X.var(1))[-150:]]
ix = {g: i for i, g in enumerate(counts.columns)}
d14 = meta[meta.day == "14 d"]
TOD = ["5AM", "10AM", "12PM", "2PM"]
mg = np.unique(np.linspace(200, 2600, 14).astype(int))
sur = {}
for tr in ("control", "drought"):
    rows = []
    for t in TOD:
        idx = [ix[g] for g in d14[(d14.treatment == tr) & (d14.tod == t)].gsm]
        r = rngs(11)
        reps = np.array([chi_by_edge_count(corr_distance(X[:, r.choice(idx, 6, replace=True)]), mg)
                         for _ in range(12)])
        rows.append(reps.mean(0))
    sur[tr] = np.array(rows)
out["switchgrass"] = {"m": mg.tolist(), "tod": TOD,
                      "control": sur['control'].round(1).tolist(),
                      "drought": sur['drought'].round(1).tolist(),
                      # transcribed from results/04_switchgrass_199perm.txt
                      "d_by_tod": [0.51, 2.26, 0.30, 0.27]}

json.dump(out, open(os.path.join(OUT, "figures.json"), "w"))
print("wrote results/figures.json")
for k, v in out.items():
    print(f"  {k}: {list(v.keys())}")
