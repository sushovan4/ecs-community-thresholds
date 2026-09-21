"""Amendment 9: the coastal marsh-forest ecotone, against a published TITAN answer.

Payne et al. (2026, Ecosphere; doi:10.5061/dryad.5tb2rbpcm) report TITAN
community change points at 0.8 and 7.6 PSU soil salinity.  Those numbers
were produced independently, by users of the incumbent method, before we
saw the data -- the only such comparison in this programme.

Primary: the two-stage procedure on log(1+cover), exactly as every other
dataset.  Secondary (registered in amendment 9): the identical pipeline on
centred-log-ratio cover, because percent cover is the one genuinely closed
substrate here.

  python3 38_ecotone.py
Output: results/38_ecotone.txt
"""
import os

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import importlib
import sys
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "ecsurf", "inst", "python"))
import ecsurf as E  # noqa: E402

mix = importlib.import_module("26_mixup_split")
ds = importlib.import_module("28_mixup_datasets")
comp = importlib.import_module("37_compositional")

CSV = os.path.join(ROOT, "data", "ecotone", "Ecotone_PlantSpp-EnvData.csv")
OUT = os.path.join(ROOT, "results", "38_ecotone.txt")
NPERM = 199
MINUNITS, SMAX = 5, 100
TITAN_CP = (0.8, 7.6)             # published, Payne et al. 2026


def build():
    d = pd.read_csv(CSV, encoding="utf-8-sig")
    c = list(d.columns)
    taxa = c[c.index("Poaceae_sp"):c.index("PHAU") + 1]      # README fields 6-79
    d = d[np.isfinite(d.Salinity.values)].copy()
    w = d[taxa].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    occ = (w > 0).sum()
    ok = occ[occ >= MINUNITS]
    keep = sorted(ok.index, key=lambda t: (-ok[t], str(t)))[:SMAX]
    return (w[keep].values.T.astype(float), d.Salinity.values.astype(float),
            d.Site.astype(str).values, keep, len(ok))


def run(X, z, reg, tag, fracs):
    mix.FRACS = tuple(fracs)
    t0 = time.time()
    groups = np.arange(len(z))
    zhat, T = ds.located(X, z)
    rng = np.random.default_rng(3)
    hits = sum(ds.located(X, rng.permutation(z))[1] <= T for _ in range(NPERM))
    p = (1 + hits) / (1 + NPERM)
    lo, hi = ds.boot_ci(X, z, groups)
    rng2 = np.random.default_rng(4)
    hr = 0
    for _ in range(NPERM):
        zp = z.copy()
        for b in np.unique(reg):
            i = np.flatnonzero(reg == b)
            zp[i] = rng2.permutation(z[i])
        hr += ds.located(X, zp)[1] <= T
    preg = (1 + hr) / (1 + NPERM)
    inside = TITAN_CP[0] <= zhat <= TITAN_CP[1]
    return (f"  {tag}: z-hat = {zhat:.3g} PSU  T = {T:.4f}  p = {p:.3f}  "
            f"within-site p = {preg:.3f}  95% CI [{lo:.3g}, {hi:.3g}]  "
            f"({len(fracs)} quantiles)  "
            f"[{'inside' if inside else 'outside'} TITAN's 0.8-7.6 PSU]  "
            f"[{time.time() - t0:.0f}s]")


def localize(X, z):
    """stage two: the curve statistic's split, on the same matrix."""
    P = len(z)
    W = max(6, int(round(P / 4)))
    r = E.estimate_threshold(X, z, nperm=NPERM, nboot=NPERM,
                             groups=np.arange(P), seed=0)
    return (f"  curve statistic (localizer): z-hat = {r['zhat']:.3g} PSU  "
            f"p = {r['p']:.3f}  95% CI [{r['ci'][0]:.3g}, {r['ci'][1]:.3g}]  "
            f"(windows of {W})")


if __name__ == "__main__":
    X, z, reg, keep, ncand = build()
    # X holds raw percent cover; the two analyses transform it differently.
    zero = float((X == 0).mean())
    L = [f"Coastal marsh-forest ecotone (Payne et al. 2026), amendment 9",
         f"{len(z)} sampling points at 3 sites, {X.shape[0]} taxa of {ncand} "
         f"passing the 5-point rule; salinity {z.min():.2f}-{z.max():.2f} PSU; "
         f"sparsity {zero:.2f}",
         f"published TITAN change points: {TITAN_CP[0]} and {TITAN_CP[1]} PSU",
         "",
         "PRIMARY -- log(1+cover), the pipeline used on every other dataset:"]
    Xl = np.log1p(X)
    L.append(run(Xl, z, reg, "interaction, widened grid",
                 np.round(np.arange(0.15, 0.851, 0.05), 2)))
    L.append(run(Xl, z, reg, "interaction, registered grid",
                 (0.3, 0.4, 0.5, 0.6, 0.7)))
    L.append(localize(Xl, z))
    L += ["", "SECONDARY -- centred log ratio, registered in amendment 9 "
          "because percent cover is closed:"]
    Xc = comp.clr(X)
    L.append(run(Xc, z, reg, "interaction, widened grid",
                 np.round(np.arange(0.15, 0.851, 0.05), 2)))
    L.append(localize(Xc, z))
    open(OUT, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
