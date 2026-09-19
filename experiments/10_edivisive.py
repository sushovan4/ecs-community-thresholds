"""Export raw simulated communities for the e-divisive baseline.

E-divisive (Matteson & James 2014) is the statistical incumbent for
multivariate change-point detection: energy-distance tests on the raw
sequence of plot abundance vectors, no summary at all.  A reviewer will ask
why the windowed topological pipeline beats just running that; this
experiment answers with detection rates on the same simulated gradients.
Ten realisations per cell; REWIRE and TRI_NULL also at P=160, the size the
shape summaries need.  Run 10_edivisive.R afterwards.
"""
import importlib
import os

import numpy as np
import pandas as pd

swap = importlib.import_module("08_summary_swap")

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "data", "ediv_sim")
os.makedirs(OUT, exist_ok=True)
CELLS = [("THRESH", 80), ("STRENGTH", 80), ("REWIRE", 80), ("TRI_NULL", 80),
         ("REWIRE", 160), ("TRI_NULL", 160)]
for scen, P in CELLS:
    for i in range(10):
        X, z = swap.community(P, scen, np.random.default_rng(600 + i))
        o = np.argsort(z)
        pd.DataFrame(X[:, o].T).to_csv(
            os.path.join(OUT, f"{scen}_{P}_{i}.csv"), index=False)
        pd.Series(z[o]).to_csv(
            os.path.join(OUT, f"{scen}_{P}_{i}_z.csv"), index=False)
print(f"wrote {10 * len(CELLS)} sequences to {OUT}")
