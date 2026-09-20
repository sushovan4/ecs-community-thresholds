"""Everglades long-term insect panel: the pre-registered powered test.

Plan: preregistration/everglades-long-term-panel.md (rules 1-8, amendment 1).
Data: Pintar & Dorn, Dryad doi:10.5061/dryad.wstqjq32d, file
Everglades-insect-community-data.csv, CERP panel (148 PSUs, 2005-2024).

  python3 17_everglades_lt.py primary     z-hat, CI, primary p  (+ TITAN inputs)
  python3 17_everglades_lt.py anchored    secondary test at TITAN's sum(z-) cp
                                          (run after 17_everglades_lt_titan.R)
"""
import os
import sys
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "ecsurf", "inst", "python"))
import ecsurf  # noqa: E402

DATA = os.path.join(ROOT, "data", "ever_lt", "Everglades-insect-community-data.csv")
OUT = os.path.join(ROOT, "results")
NPERM = NBOOT = 199


def build():
    d = pd.read_csv(DATA)
    c = d[d.PROJECT == "CERP"].copy()
    cols = list(d.columns)
    taxa = cols[cols.index("COLEOPTERA"):cols.index("EPITHECA") + 1]
    c["PSU"] = (c.REGION.astype(str) + "_" + c.SITE.astype(str) + "_"
                + c.PLOT.astype(str))
    # rule 2 / amendment 3: PSU-year density = summed abundance / summed throws
    agg = c.groupby(["PSU", "YEAR"]).agg({**{t: "sum" for t in taxa},
                                          "THROW": "sum"})
    dens = agg[taxa].div(agg["THROW"], axis=0)
    # rule 3: taxa present in >= 5 PSUs
    present = (dens.groupby(level="PSU").sum() > 0).sum()
    keep = [t for t in taxa if present[t] >= 5]
    dens = dens[keep]
    # rule 4a / amendment 4: z = PSU mean of DSLDD over its survey years
    z_psu = c.groupby("PSU").DSLDD.mean()
    psu = dens.index.get_level_values("PSU")
    return dens, keep, z_psu, psu


def primary():
    dens, keep, z_psu, psu = build()
    X = np.log1p(dens.values.T)                     # taxa x PSU-years
    z = z_psu.loc[psu].values
    groups = pd.factorize(psu)[0]
    t0 = time.time()
    res = ecsurf.estimate_threshold(X, z, nperm=NPERM, nboot=NBOOT,
                                    groups=groups, seed=0)
    lines = [
        "Everglades long-term insect panel (CERP), pre-registered primary test",
        f"{len(keep)} taxa x {X.shape[1]} PSU-year columns, "
        f"{len(z_psu)} PSUs; z = PSU mean days since last dry "
        f"({z_psu.min():.0f}-{z_psu.max():.0f} d)",
        f"windows of {res['W']} PSUs; {NPERM} permutations, {NBOOT} bootstraps",
        f"z-hat = {res['zhat']:.1f} d   stat {res['stat']:.2f}   "
        f"p = {res['p']:.3f}   95% CI [{res['ci'][0]:.1f}, {res['ci'][1]:.1f}]"
        f"   [{time.time() - t0:.0f}s]",
    ]
    open(os.path.join(OUT, "17_everglades_lt_primary.txt"), "w").write(
        "\n".join(lines) + "\n")
    print("\n".join(lines))
    # inputs for the secondary (TITAN) test: PSU-mean densities, not logged
    tdir = os.path.join(ROOT, "data", "ever_lt")
    dens.groupby(level="PSU").mean().loc[z_psu.index].to_csv(
        os.path.join(tdir, "titan_taxa.csv"), index=False)
    pd.Series(z_psu.values, name="z").to_csv(
        os.path.join(tdir, "titan_env.csv"), index=False)


def anchored():
    cp = float(open(os.path.join(OUT, "17_everglades_lt_titan_cp.txt")).read())
    dens, keep, z_psu, psu = build()
    X = np.log1p(dens.values.T)
    z = z_psu.loc[psu].values
    groups = pd.factorize(psu)[0]
    units = np.unique(groups)
    P = len(units)
    W = max(6, int(round(P / 4)))
    stride = max(1, int(round(W / 3)))
    m_grid = ecsurf._default_m_grid(X.shape[0])
    C, zc = ecsurf._curves(X, z, groups, W, stride, m_grid)
    mids = np.array([(zc[s - 1] + zc[s]) / 2 for s in range(2, len(zc) - 1)])
    s0 = int(np.argmin(np.abs(mids - cp))) + 2

    def d_at(Cw, s):
        a, b = Cw[:s], Cw[s:]
        pooled = np.sqrt((a.var(0, ddof=1) + b.var(0, ddof=1)) / 2) + 1e-12
        return float(np.abs((a.mean(0) - b.mean(0)) / pooled).mean())

    obs = d_at(C, s0)
    zu = np.array([z[groups == u][0] for u in units])
    rng = np.random.default_rng(1)
    hits = 0
    for _ in range(NPERM):
        zpu = rng.permutation(zu)
        zp = zpu[groups]
        Cb, _ = ecsurf._curves(X, zp, groups, W, stride, m_grid)
        hits += d_at(Cb, s0) >= obs
    p = (1 + hits) / (1 + NPERM)
    line = (f"anchored split at boundary {mids[s0 - 2]:.1f} d "
            f"(TITAN sum(z-) cp {cp:.1f} d): D = {obs:.2f}   "
            f"p = {p:.3f} ({NPERM} perms)")
    open(os.path.join(OUT, "17_everglades_lt_primary.txt"), "a").write(
        line + "\n")
    print(line)


if __name__ == "__main__":
    {"primary": primary, "anchored": anchored}[sys.argv[1]]()
