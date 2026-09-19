"""Real data.  Does gene co-expression structure differ TOPOLOGICALLY between
early- and late-flowering Arabidopsis accessions, at matched edge density?

Data: 1001 Genomes expression (GEO GSE80744) x flowering time at 16C.
Genes are points; coordinates are expression across a group's accessions.
Control: a label-permutation null, because bootstrap replicates alone will
manufacture a large Cohen d out of nothing.
"""
import os
import time
import numpy as np
import pandas as pd
from ec import corr_distance, chi_by_edge_count, cohen_d

NGENE = int(os.environ.get("NGENE", 150))
NACC  = int(os.environ.get("NACC", 90))
BOOT  = int(os.environ.get("BOOT", 12))
NPERM = int(os.environ.get("NPERM", 199))
DATA  = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data") + os.sep

expr = pd.read_csv(DATA + "ath_tx.tsv.gz", sep="\t", index_col=0)
expr.columns = [c.lstrip("X") for c in expr.columns]

ft = pd.read_csv(DATA + "FT_16C.txt", sep="\t", index_col=0)["FT16_mean"]
ft.index = ft.index.astype(str)
ft = ft[~ft.index.duplicated()]                     # the phenotype file repeats ids

common = [a for a in expr.columns if a in ft.index]
expr, ft = expr[common], ft.loc[common]
assert expr.shape[1] == len(ft), (expr.shape, len(ft))

X = np.log2(expr.values.astype(np.float64) + 1.0)   # genes x accessions
X = X[np.argsort(X.var(1))[-NGENE:]]                # label-blind gene selection
print(f"{X.shape[0]} genes x {X.shape[1]} accessions matched to FT16")
print(f"FT16 {ft.min():.0f}-{ft.max():.0f} d, median {ft.median():.0f}\n")

order = np.argsort(ft.values)
third = len(order) // 3
early, late = order[:third], order[-third:]

m_grid = np.unique(np.linspace(200, 2600, 14).astype(int))


def replicates(pool, rng):
    return np.array([chi_by_edge_count(corr_distance(X[:, rng.choice(pool, NACC, replace=True)]),
                                       m_grid) for _ in range(BOOT)])


t0 = time.time()
rng = np.random.default_rng(0)
obs = np.abs(cohen_d(replicates(early, rng), replicates(late, rng)))
print(f"observed  max |Cohen d| = {obs.max():.2f} at {m_grid[obs.argmax()]} edges"
      f"   [{time.time()-t0:.0f}s]")

null = []
for p in range(NPERM):
    r = np.random.default_rng(1000 + p)
    perm = r.permutation(np.concatenate([early, late]))
    null.append(np.abs(cohen_d(replicates(perm[:third], r),
                               replicates(perm[third:], r))).max())
null = np.array(null)
pval = (1 + (null >= obs.max()).sum()) / (1 + NPERM)

print(f"permuted  max |Cohen d| = {null.mean():.2f} mean, "
      f"{null.min():.2f}-{null.max():.2f} range over {NPERM} label shuffles")
print(f"\nempirical p = {pval:.3f}    [total {time.time()-t0:.0f}s]")
