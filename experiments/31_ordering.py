"""The pre-registered cross-dataset ordering test (amendment 6).

Detection comes from the interaction test (28_mixup_*.txt), localization from
the curve statistic (19_multi_*.txt, 17_everglades_lt_primary.txt) because
29_mixup_localization.txt showed it to be the more accurate recoverer, and
the taxon-level change point from TITAN (25_titan_*.txt).

  R = (z_TITAN - z_struct) / IQR(z over units),   R > 0: structure first.

The confirmatory test is a two-sided Wilcoxon signed-rank test of R against
zero across INDEPENDENT detecting datasets (D1-D8 and the Everglades
long-term panel; the Everglades temporal series are excluded as not
independent of it).  It runs only if at least four such datasets detect.
Output: results/31_ordering.txt
"""
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")
sys.path.insert(0, os.path.join(ROOT, "experiments"))


def grab(path, pattern, group=1):
    if not os.path.exists(path):
        return None
    m = re.search(pattern, open(path).read())
    return float(m.group(group)) if m else None


def iqr_of(name):
    import importlib
    if name == "ELT":
        lt = importlib.import_module("17_everglades_lt")
        _, _, z_psu, _ = lt.build()
        return float(np.subtract(*np.percentile(z_psu.values, [75, 25])))
    multi = importlib.import_module("19_multi")
    w, zc, reg, zn = multi.BUILD[name]()
    X, z, g, units, zunit, keep, nc, wide = multi.prepare(w, zc)
    return float(np.subtract(*np.percentile(zunit.values, [75, 25])))


NAMES = ["D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "ELT"]

rows = []
for n in NAMES:
    mix_f = f"{R}/28_mixup_{n}.txt"
    p_mix = grab(mix_f, r"p = ([0-9.]+)")
    p_reg = grab(mix_f, r"within-region \([^)]*\): p = ([0-9.]+)")
    if n == "ELT":
        z_struct = grab(f"{R}/17_everglades_lt_primary.txt", r"z-hat = ([0-9.]+)")
        z_titan = grab(f"{R}/17_everglades_lt_titan_cp.txt", r"([0-9.]+)")
    else:
        z_struct = grab(f"{R}/19_multi_{n}.txt", r"z-hat = ([-0-9.e+]+)")
        t = f"{R}/25_titan_{n}.txt"
        z_titan = grab(t, r"sumz-\s+([-0-9.e+]+)")
    detect = (p_mix is not None and p_mix <= 0.05
              and (p_reg is None or p_reg <= 0.05))
    rows.append(dict(dataset=n, p_mix=p_mix, p_region=p_reg, detect=detect,
                     z_struct=z_struct, z_titan=z_titan))

df = pd.DataFrame(rows)
df["iqr"] = [iqr_of(n) if (r.z_struct is not None and r.z_titan is not None)
             else np.nan for n, r in zip(df.dataset, df.itertuples())]
df["Rstat"] = (df.z_titan - df.z_struct) / df.iqr

L = ["Cross-dataset ordering test (amendment 6)", "",
     f"{'dataset':>8} {'p(mixup)':>9} {'p(region)':>10} {'detect':>7} "
     f"{'z_struct':>10} {'z_TITAN':>10} {'R':>8}"]
for r in df.itertuples():
    L.append(f"{r.dataset:>8} {str(r.p_mix):>9} {str(r.p_region):>10} "
             f"{str(r.detect):>7} {r.z_struct if r.z_struct is not None else float('nan'):10.4g} "
             f"{r.z_titan if r.z_titan is not None else float('nan'):10.4g} "
             f"{r.Rstat:8.3f}")

ok = df[df.detect & df.Rstat.notna()]
L += ["", f"detecting datasets with both change points: {len(ok)}"]
if len(ok) >= 4:
    stat, p = wilcoxon(ok.Rstat.values)
    pos = int((ok.Rstat > 0).sum())
    L += [f"Wilcoxon signed-rank of R against zero: W = {stat:.1f}, p = {p:.4f}",
          f"sign test: {pos} of {len(ok)} have R > 0 (structure at a lower "
          f"gradient value)",
          f"median R = {np.median(ok.Rstat):.3f} "
          f"(in units of the gradient's interquartile range)"]
else:
    L += ["fewer than four independent detecting datasets with both change "
          "points: the confirmatory test is not run, and the R values above "
          "stand as description only."]
open(f"{R}/31_ordering.txt", "w").write("\n".join(L) + "\n")
print("\n".join(L))
