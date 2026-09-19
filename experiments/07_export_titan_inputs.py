"""Export simulated gradient communities for the TITAN2 comparison.

Five realisations each of two scenarios from 06_gradient_estimator.py at the
Butte-scale design (P = 80):

  THRESH  correlation structure reorganises at z* = 0.55; every taxon's
          marginal abundance distribution is constant along the gradient
  ABUND   the converse: module-0 taxa shift mean abundance at z*, the
          correlation structure never changes

TITAN scores taxon-level abundance change, so it should find the ABUND
threshold and nothing in THRESH; the ECS estimator should do the reverse
(see the power table in results/06_estimator_slow.txt).  Abundances are
exponentiated (lognormal) so TITAN sees nonnegative values; the monotone
transform does not affect which scenario carries taxon-level signal.

Outputs data/titan_sim/{scenario}_{i}_{taxa,env}.csv -- regenerable,
gitignored with the rest of data/.
"""
import importlib
import os

import numpy as np
import pandas as pd

sim = importlib.import_module("06_gradient_estimator")

P, NREAL = 80, 5
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "data", "titan_sim")
os.makedirs(OUT, exist_ok=True)

for scen in ("THRESH", "ABUND"):
    for i in range(NREAL):
        X, z = sim.community(P, scen, np.random.default_rng(600 + i))
        taxa = pd.DataFrame(np.exp(X.T),
                            columns=[f"t{j:02d}m{j // (sim.S // sim.NMOD)}"
                                     for j in range(sim.S)])
        taxa.to_csv(os.path.join(OUT, f"{scen}_{i}_taxa.csv"), index=False)
        pd.Series(z, name="z").to_csv(os.path.join(OUT, f"{scen}_{i}_env.csv"),
                                      index=False)
print(f"wrote {2 * NREAL} realisations to {OUT} (true z* = {sim.ZSTAR})")
