"""Temporal precedence: does structure reorganize before composition?

Plan: preregistration/temporal-precedence.md.  The estimator of paper
Section 2.4 is applied unchanged with the gradient set to YEAR: units are
survey years, columns are the plot-years within them, groups = year.

  python3 22_temporal.py list                 series inventory (structure only)
  python3 22_temporal.py run <SERIES> ...     primary test + TITAN export
  python3 22_temporal.py oos <SERIES> ...     secondary, out-of-sample y*
Outputs results/22_temporal_<SERIES>.txt and, for TITAN, data/temporal/.
"""
import os

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "ecsurf", "inst", "python"))
import ecsurf as E  # noqa: E402

OUT = os.path.join(ROOT, "results")
TDIR = os.path.join(ROOT, "data", "temporal")
NPERM = NBOOT = 199
SMAX, MINUNITS = 100, 5
MINYEARS, MINCOLS, MINTAXA = 10, 40, 15   # T-series need 40+ occasions (amendment 2C)
REGIONS = ("PHD", "SRS", "TSL", "WCA")


def everglades(project, region):
    f = os.path.join(ROOT, "data", "ever_lt",
                     "Everglades-insect-community-data.csv")
    d = pd.read_csv(f, low_memory=False)
    cols = list(d.columns)
    taxa = cols[cols.index("COLEOPTERA"):cols.index("EPITHECA") + 1]
    d = d[(d.PROJECT == project) & (d.REGION == region)].copy()
    if d.empty:
        return None
    d["PLOT_KEY"] = (d.REGION.astype(str) + "_" + d.SITE.astype(str) + "_"
                     + d.PLOT.astype(str))
    d["COL"] = d.PLOT_KEY + "|" + d.YEAR.astype(str)
    agg = d.groupby(["YEAR", "COL"]).agg({**{t: "sum" for t in taxa},
                                          "THROW": "sum"})
    dens = agg[taxa].div(agg["THROW"], axis=0)
    return dens


def everglades_occasion(region):
    """Amendment 2C: the unit is the survey OCCASION (year x period), not the
    year; MWD samples five periods a year.  z = decimal date."""
    f = os.path.join(ROOT, "data", "ever_lt",
                     "Everglades-insect-community-data.csv")
    d = pd.read_csv(f, low_memory=False)
    cols = list(d.columns)
    taxa = cols[cols.index("COLEOPTERA"):cols.index("EPITHECA") + 1]
    d = d[(d.PROJECT == "MWD") & (d.REGION == region)].copy()
    d = d.dropna(subset=["PERIOD"])
    if d.empty:
        return None
    d["OCC"] = d.YEAR + (d.PERIOD.astype(float) - 1) / 5.0
    d["PLOT_KEY"] = (d.REGION.astype(str) + "_" + d.SITE.astype(str) + "_"
                     + d.PLOT.astype(str))
    d["COL"] = d.PLOT_KEY + "|" + d.OCC.astype(str)
    agg = d.groupby(["OCC", "COL"]).agg({**{t: "sum" for t in taxa},
                                         "THROW": "sum"})
    return agg[taxa].div(agg["THROW"], axis=0)


def cedar(field):
    f = os.path.join(ROOT, "data", "multi", "cdr", "e001.csv")
    d = pd.read_csv(f, low_memory=False).rename(
        columns={"Biomass(g/m2)": "Biomass"})
    bad = d.Species.astype(str).str.contains(
        "litter|miscellaneous|moss|lichen|fung|unsorted|total", case=False)
    d = d[~bad & (d.Biomass > 0) & (d.Field == field)].copy()
    if d.empty:
        return None
    d["COL"] = d.Field + "_" + d.Plot.astype(str) + "|" + d.Year.astype(str)
    w = d.pivot_table(index=["Year", "COL"], columns="Species",
                      values="Biomass", aggfunc="sum", fill_value=0.0)
    return w


def niwot():
    f = os.path.join(ROOT, "data", "multi", "nwt", "plants.csv")
    p = pd.read_csv(f, low_memory=False)
    p = p[~p.USDA_code.astype(str).str.startswith("2")].copy()
    p["one"] = 1.0
    p["COL"] = p["plot"].astype(str) + "|" + p["year"].astype(str)
    return p.pivot_table(index=["year", "COL"], columns="USDA_code",
                         values="one", aggfunc="sum", fill_value=0.0)


SERIES = {}
for i, r in enumerate(REGIONS):
    SERIES[f"E{i + 1}"] = (lambda r=r: everglades("MWD", r))
    SERIES[f"E{i + 5}"] = (lambda r=r: everglades("CERP", r))
for i, fl in enumerate("ABCD"):
    SERIES[f"C{i + 1}"] = (lambda fl=fl: cedar(fl))
SERIES["N1"] = niwot
for i, r in enumerate(REGIONS):
    SERIES[f"T{i + 1}"] = (lambda r=r: everglades_occasion(r))
LABEL = {**{f"E{i + 1}": f"Everglades MWD {r}" for i, r in enumerate(REGIONS)},
         **{f"E{i + 5}": f"Everglades CERP {r}" for i, r in enumerate(REGIONS)},
         **{f"C{i + 1}": f"Cedar Creek field {fl}" for i, fl in enumerate("ABCD")},
         **{f"T{i + 1}": f"Everglades MWD {r}, by occasion"
            for i, r in enumerate(REGIONS)},
         "N1": "Niwot Saddle grid"}


def prepare(w):
    """taxa filters and cap; returns X (taxa x columns), z (year per column)."""
    years = w.index.get_level_values(0).values.astype(float)
    occ = (w.groupby(level=0).sum() > 0).sum()
    ok = occ[occ >= MINUNITS]
    keep = sorted(ok.index, key=lambda t: (-ok[t], str(t)))[:SMAX]
    w = w[keep]
    X = np.log1p(w.values.T.astype(float))
    groups = pd.factorize(years)[0]
    return X, years, groups, keep, len(ok), w


def inventory():
    rows = []
    for s, build in SERIES.items():
        w = build()
        if w is None or len(w) == 0:
            rows.append((s, LABEL[s], 0, 0, 0, "no rows"))
            continue
        X, z, g, keep, ncand, _ = prepare(w)
        P, C, S = len(np.unique(z)), X.shape[1], X.shape[0]
        why = ("" if P >= MINYEARS and C >= MINCOLS and S >= MINTAXA
               else "excluded: " + ", ".join(
                   x for x, ok in (("<10 years", P >= MINYEARS),
                                   ("<40 columns", C >= MINCOLS),
                                   ("<15 taxa", S >= MINTAXA)) if not ok))
        rows.append((s, LABEL[s], P, C, S, why))
    L = [f"{'id':>4} {'series':>26} {'years':>6} {'cols':>6} {'taxa':>5}  note"]
    for r in rows:
        L.append(f"{r[0]:>4} {r[1]:>26} {r[2]:>6} {r[3]:>6} {r[4]:>5}  {r[5]}")
    open(os.path.join(OUT, "22_temporal_inventory.txt"), "w").write(
        "\n".join(L) + "\n")
    print("\n".join(L))


def run(s):
    t0 = time.time()
    w = SERIES[s]()
    X, z, g, keep, ncand, wide = prepare(w)
    P, C, S = len(np.unique(z)), X.shape[1], X.shape[0]
    minu = 40 if s.startswith("T") else MINYEARS
    L = [f"{s} ({LABEL[s]}): {P} survey occasions {int(z.min())}-{int(z.max())}, "
         f"{C} plot-year columns, {S} taxa of {ncand}"]
    if P < minu or C < MINCOLS or S < MINTAXA:
        L.append("EXCLUDED (inclusion rule)")
    else:
        r = E.estimate_threshold(X, z, nperm=NPERM, nboot=NBOOT, groups=g,
                                 seed=0)
        L.append(f"structural: z-hat = {r['zhat']:.1f}  stat {r['stat']:.3f}  "
                 f"p = {r['p']:.3f}  95% CI [{r['ci'][0]:.1f}, {r['ci'][1]:.1f}]"
                 f"  (windows of {r['W']} years)")
        os.makedirs(TDIR, exist_ok=True)
        ym = wide.groupby(level=0).mean()
        ym.to_csv(f"{TDIR}/{s}_taxa.csv", index=False)
        pd.Series(ym.index.values, name="z").to_csv(f"{TDIR}/{s}_env.csv",
                                                    index=False)
    L.append(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(OUT, f"22_temporal_{s}.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L), flush=True)


def oos(s):
    """secondary: earliest cut-off year at which the test already rejects."""
    t0 = time.time()
    w = SERIES[s]()
    X, z, g, keep, ncand, _ = prepare(w)
    yrs = np.unique(z)
    lines = []
    ystar = None
    for k in range(8, len(yrs) + 1):
        sub = z <= yrs[k - 1]
        Xs, zs = X[:, sub], z[sub]
        gs = pd.factorize(zs)[0]
        r = E.estimate_threshold(Xs, zs, nperm=NPERM, nboot=0, groups=gs,
                                 seed=1)
        lines.append(f"  years <= {int(yrs[k - 1])}: z-hat {r['zhat']:.1f}, "
                     f"p = {r['p']:.3f}")
        if r["p"] <= 0.05 and ystar is None:
            ystar = yrs[k - 1]
            break
    head = (f"{s} ({LABEL[s]}) out-of-sample: y* = "
            f"{int(ystar) if ystar else 'not reached'}")
    out = "\n".join([head] + lines + [f"[{time.time() - t0:.0f}s]"])
    open(os.path.join(OUT, f"22_temporal_oos_{s}.txt"), "w").write(out + "\n")
    print(out, flush=True)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "list":
        inventory()
    else:
        for s in sys.argv[2:]:
            (run if cmd == "run" else oos)(s)
