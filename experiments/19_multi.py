"""The frozen estimator across public gradient datasets (D1-D7).

Plan: preregistration/multi-dataset-program.md (rules 1-9, amendment 1).
Usage:  python3 19_multi.py D1 [D2 ...]      -> results/19_multi_<D>.txt
Each run: primary (199 perms, 199 boots); 1999-perm precision check when the
primary p is in [0.01, 0.10]; within-region permutation test (199) where a
region field is registered; unit-mean densities exported for TITAN.
"""
import io
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

M = os.path.join(ROOT, "data", "multi")
OUT = os.path.join(ROOT, "results")
NPERM = NBOOT = 199
SMAX, MINUNITS = 100, 5


def _wide(long, unit, col, taxon, abund):
    """long table -> (unit-col) x taxa abundance, summed."""
    w = long.pivot_table(index=[unit, col], columns=taxon, values=abund,
                         aggfunc="sum", fill_value=0.0)
    w.index.names = ["unit", "col"]
    return w


def _emap(f):
    lines = open(f, encoding="latin-1").read().splitlines()
    h = next(i for i, l in enumerate(lines) if l.startswith('"'))
    d = pd.read_csv(io.StringIO("\n".join(lines[h:])), na_values=["."],
                    low_memory=False, skipinitialspace=True)
    d.columns = [c.strip() for c in d.columns]
    return d


def d1(appalachia=False):
    b = pd.read_csv(f"{M}/nrsa/bent.csv", low_memory=False)
    c = pd.read_csv(f"{M}/nrsa/chem.csv", low_memory=False)
    b = b[(b.IS_DISTINCT300 == 1) & (b.TOTAL300 > 0)]
    if appalachia:
        b = b[b.AG_ECO9.isin(["NAP", "SAP"])]
    w = _wide(b, "UNIQUE_ID", "UID", "TARGET_TAXON", "TOTAL300")
    zc = c.groupby("UID").COND_RESULT.mean()
    reg_field = "STATE" if appalachia else "AG_ECO9"
    reg = b.drop_duplicates("UNIQUE_ID").set_index("UNIQUE_ID")[reg_field]
    return w, zc, reg, "specific conductance (uS/cm)"


def d3():
    b = _emap(f"{M}/emap/bencnt.txt")
    c = _emap(f"{M}/emap/chmval.txt")
    b = b[(b.DISTINCT == "Y") & (b.ABUND > 0)].copy()
    b["COL"] = b.STRM_ID + "|" + b.YEAR.astype(int).astype(str) + "|" + \
        b.VISIT_NO.astype(int).astype(str)
    c = c.dropna(subset=["STRM_ID", "YEAR", "VISIT_NO", "COND"]).copy()
    c["COL"] = c.STRM_ID + "|" + c.YEAR.astype(int).astype(str) + "|" + \
        c.VISIT_NO.astype(int).astype(str)
    w = _wide(b, "STRM_ID", "COL", "TAXANAME", "ABUND")
    zc = c.groupby("COL").COND.mean()
    units = w.index.get_level_values("unit").unique()
    reg = pd.Series([u[:2] for u in units], index=units)
    return w, zc, reg, "specific conductance (uS/cm)"


def d4():
    b = pd.read_csv(f"{M}/ncca/bent.csv", low_memory=False)
    s = pd.read_csv(f"{M}/ncca/sedchem.csv", low_memory=False)
    b = b[(b.IS_DISTINCT.astype(str).isin(["1", "1.0", "Y"])) & (b.TOTAL > 0)].copy()
    b["UID"] = pd.to_numeric(b.UID).astype("int64").astype(str)
    w = _wide(b, "SITE_ID", "UID", "TARGET_TAXON", "TOTAL")
    cu = s[s.ANALYTE == "CU"].copy()
    cu["UID"] = pd.to_numeric(cu.UID).astype("int64").astype(str)
    zc = pd.to_numeric(cu.RESULT, errors="coerce").groupby(cu.UID).mean()
    reg = b.drop_duplicates("SITE_ID").set_index("SITE_ID").NCCA_REG
    return w, zc, reg, "sediment copper (ug/dry g)"


def d5():
    b = pd.read_csv(f"{M}/nla/bent.csv", low_memory=False)
    c = pd.read_csv(f"{M}/nla/chem.csv", low_memory=False)
    s = pd.read_csv(f"{M}/nla/site.csv", low_memory=False)
    b = b[(b.IS_DISTINCT == 1) & (b.TOTAL > 0)]
    w = _wide(b, "SITE_ID", "UID", "TARGET_TAXON", "TOTAL")
    p = c[c.ANALYTE == "PTL"]
    zc = pd.to_numeric(p.RESULT, errors="coerce").groupby(p.UID).mean()
    reg = s.drop_duplicates("SITE_ID").set_index("SITE_ID").AG_ECO9
    return w, zc, reg, "total phosphorus"


def d6():
    d = pd.read_csv(f"{M}/cdr/e001.csv", low_memory=False)
    bad = d.Species.astype(str).str.contains(
        "litter|miscellaneous|moss|lichen|fung|unsorted|total", case=False)
    d = d.rename(columns={"Biomass(g/m2)": "Biomass"})
    d = d[~bad & (d.Biomass > 0)].copy()
    d["UNIT"] = d.Field + "_" + d.Plot.astype(str)
    d["COL"] = d.UNIT + "|" + d.Year.astype(str)
    w = _wide(d, "UNIT", "COL", "Species", "Biomass")
    zc = d.groupby("COL").NAdd.mean()
    reg = d.drop_duplicates("UNIT").set_index("UNIT").Field
    return w, zc, reg, "N addition (g N m-2 yr-1)"


def d7():
    p = pd.read_csv(f"{M}/nwt/plants.csv", low_memory=False)
    s = pd.read_csv(f"{M}/nwt/snow.csv", low_memory=False)
    p = p[~p.USDA_code.astype(str).str.startswith("2")].copy()
    p["one"] = 1.0
    p["COL"] = p["plot"].astype(str) + "|" + p["year"].astype(str)
    w = _wide(p, "plot", "COL", "USDA_code", "one")
    zp = pd.to_numeric(s.mean_depth, errors="coerce").groupby(s.point_ID).mean()
    cols = w.index.get_level_values("col")
    units = w.index.get_level_values("unit")
    zc = pd.Series(zp.reindex(units).values, index=cols)
    return w, zc, None, "mean snow depth (cm)"


def d8():
    """CRMS Louisiana marsh vegetation; amendments 3-4."""
    v = pd.read_csv(f"{M}/../crms/CRMS_Marsh_Vegetation.csv", low_memory=False,
                    encoding="latin-1",
                    usecols=["Station ID", "Collection Date (mm/dd/yyyy)",
                             "Scientific Name As Currently Recognized",
                             "% Cover"])
    v.columns = [c.strip() for c in v.columns]
    v = v.rename(columns={"Station ID": "unit",
                          "Collection Date (mm/dd/yyyy)": "date",
                          "Scientific Name As Currently Recognized": "taxon",
                          "% Cover": "cover"})
    v = v.dropna(subset=["unit", "taxon", "cover"])
    v = v[pd.to_numeric(v.cover, errors="coerce") > 0].copy()
    v["cover"] = pd.to_numeric(v.cover, errors="coerce")
    v["year"] = pd.to_datetime(v.date, errors="coerce").dt.year
    v = v.dropna(subset=["year"])
    v["COL"] = v.unit.astype(str) + "|" + v.year.astype(int).astype(str)
    w = _wide(v, "unit", "COL", "taxon", "cover")
    e = pd.read_csv(f"{M}/../crms/elev.csv", low_memory=False, encoding="latin-1")
    e.columns = [c.strip() for c in e.columns]
    z_unit = pd.to_numeric(e["Marsh_Elevation_NAVD88_FT(GEOID12A)"],
                           errors="coerce").groupby(e.StationID).mean()
    cols = w.index.get_level_values("col")
    units = w.index.get_level_values("unit")
    zc = pd.Series(z_unit.reindex(units).values, index=cols)
    reg = e.drop_duplicates("StationID").set_index("StationID").Basin
    return w, zc, reg, "station elevation (ft NAVD88)"


BUILD = {"D1": d1, "D2": lambda: d1(appalachia=True), "D3": d3, "D4": d4,
         "D5": d5, "D6": d6, "D7": d7, "D8": d8}


def prepare(w, zc):
    """Apply rules 3-4 and the amendment's taxon cap; unit-mean z."""
    cols = w.index.get_level_values("col")
    zcol = pd.Series(zc.reindex(cols).values, index=w.index)
    w = w[zcol.notna().values]
    zcol = zcol[zcol.notna().values]
    units = w.index.get_level_values("unit")
    zunit = zcol.groupby(level="unit").mean()
    occ = (w.groupby(level="unit").sum() > 0).sum()
    ok = occ[occ >= MINUNITS]
    ncand = len(ok)
    keep = sorted(ok.index, key=lambda t: (-ok[t], str(t)))[:SMAX]
    w = w[keep]
    X = np.log1p(w.values.T.astype(float))
    z = zunit.reindex(units).values
    groups = pd.factorize(units)[0]
    return X, z, groups, units, zunit, keep, ncand, w


def region_test(X, z, groups, region_of_group, obs, seed=5):
    units = np.unique(groups)
    P = len(units)
    W = max(6, int(round(P / 4)))
    stride = max(1, int(round(W / 3)))
    mg = E._default_m_grid(X.shape[0])
    zu = np.array([z[groups == u][0] for u in units])
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(NPERM):
        zp = zu.copy()
        for r in np.unique(region_of_group):
            i = np.flatnonzero(region_of_group == r)
            zp[i] = rng.permutation(zu[i])
        hits += E._split_stat(*E._curves(X, zp[groups], groups, W, stride,
                                         mg))[1] >= obs
    return (1 + hits) / (1 + NPERM)


def run(D):
    t0 = time.time()
    w, zc, reg, zname = BUILD[D]()
    X, z, groups, units, zunit, keep, ncand, wide = prepare(w, zc)
    P, S, C = len(zunit), X.shape[0], X.shape[1]
    zero = float((wide.values == 0).mean())
    L = [f"{D}: {P} units, {C} columns, {S} taxa kept of {ncand} passing the "
         f"5-unit rule; z = {zname}, unit range {zunit.min():.3g}-{zunit.max():.3g}",
         f"matrix sparsity: {zero:.2f} of entries are zero"]
    if P < 40 or S < 15:
        L.append("EXCLUDED under rule 1 (fewer than 40 units or 15 taxa)")
    else:
        r = E.estimate_threshold(X, z, nperm=NPERM, nboot=NBOOT,
                                 groups=groups, seed=0)
        L.append(f"primary: z-hat = {r['zhat']:.4g}  stat {r['stat']:.3f}  "
                 f"p = {r['p']:.3f}  95% CI [{r['ci'][0]:.4g}, {r['ci'][1]:.4g}]"
                 f"  (windows of {r['W']} units)")
        if 0.01 <= r["p"] <= 0.10:
            q = E.estimate_threshold(X, z, nperm=1999, nboot=0, groups=groups,
                                     seed=7)
            L.append(f"precision check (1999 perms, same test): p = {q['p']:.4f}")
        if reg is not None:
            ug = pd.unique(units)
            rg = reg.reindex(ug).fillna("NA").values
            pr = region_test(X, z, groups, rg, r["stat"])
            L.append(f"within-region permutation ({len(set(rg))} regions): "
                     f"p = {pr:.3f}")
        tdir = os.path.join(M, "titan")
        os.makedirs(tdir, exist_ok=True)
        wide.groupby(level="unit").mean().loc[zunit.index].to_csv(
            f"{tdir}/{D}_taxa.csv", index=False)
        pd.Series(zunit.values, name="z").to_csv(f"{tdir}/{D}_env.csv",
                                                 index=False)
    L.append(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(OUT, f"19_multi_{D}.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L), flush=True)


if __name__ == "__main__":
    for D in sys.argv[1:]:
        run(D)
