"""Amendment 8: do the two channels run on different taxa?

For each detecting dataset, decompose the detected structural change into
per-taxon curvature contributions (discrete Gauss-Bonnet, exactly as in
21_curvature_attribution.py), then ask what share of that change is carried
by TITAN's pure and reliable indicator taxa on the same matrix:

    C = ( sum_{i in I} |dkappa_i| / sum_i |dkappa_i| ) / ( |I| / S ).

C = 1 is the proportional share; C < 1 means the structural change is
carried by taxa TITAN does not flag.  The p-value draws subsets of size |I|
uniformly from the S taxa (9999 draws), which is exact under the null that
the indicator set is unrelated to the curvature shares.

  python3 36_attribution_multi.py D1 D2 ...    per-dataset attribution
  python3 36_attribution_multi.py test         the confirmatory test
Outputs: results/36_attribution_<D>.txt, results/36_attribution.txt
"""
import os

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import importlib
import re
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, wilcoxon

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")
sys.path.insert(0, os.path.join(ROOT, "ecsurf", "inst", "python"))
import ecsurf as E  # noqa: E402

att = importlib.import_module("21_curvature_attribution")
multi = importlib.import_module("19_multi")

NDRAW = 9999
NAMES = [f"D{i}" for i in range(1, 9)]


def rname(s):
    """R's make.names, enough of it for taxon labels."""
    s = re.sub(r"[^A-Za-z0-9._]", ".", str(s))
    return s if re.match(r"^[A-Za-z.]", s) else "X" + s


def detects(D):
    """amendment-5 rule: the interaction test rejects, and the within-region
    test too where a region field exists."""
    f = os.path.join(R, f"28_mixup_{D}.txt")
    if not os.path.exists(f):
        return None
    t = open(f).read()
    p = re.search(r"p = ([0-9.]+)", t)
    pr = re.search(r"within-region \([^)]*\): p = ([0-9.]+)", t)
    if not p:
        return None
    return float(p.group(1)) <= 0.05 and (pr is None or float(pr.group(1)) <= 0.05)


def indicators(D):
    """TITAN's pure and reliable taxa from the seeded run (amendment 8)."""
    f = os.path.join(R, f"35_titan_{D}_sppmax.csv")
    if not os.path.exists(f):
        return None
    sp = pd.read_csv(f)
    col = "filter" if "filter" in sp.columns else sp.columns[-2]
    return set(sp.loc[sp[col] > 0, "taxon"].astype(str))


def attribute(D):
    """|dkappa| per taxon at the curve statistic's split."""
    t0 = time.time()
    w, zc, reg_s, zname = multi.BUILD[D]()
    X, z, groups, units, zunit, keep, ncand, wide = multi.prepare(w, zc)
    P = len(np.unique(groups))
    W = max(6, int(round(P / 4)))
    stride = max(1, int(round(W / 3)))
    mg = E._default_m_grid(X.shape[0])
    K, zcw = att.window_kappa(X, z, groups, W, stride, mg)
    chi = K.sum(1)
    C_chk, _ = E._curves(X, z, groups, W, stride, mg)
    assert np.allclose(chi, C_chk), "Gauss-Bonnet identity failed"
    s = att.split_index(chi, zcw)
    d = att.taxon_d(K, s)
    out = pd.DataFrame({"taxon": [str(k) for k in keep], "dkappa": d})
    out.to_csv(os.path.join(R, f"36_attribution_{D}.csv"), index=False)
    L = [f"{D}: {P} units, {len(keep)} taxa; split at "
         f"z = {(zcw[s - 1] + zcw[s]) / 2:.4g} ({zname})",
         f"change in chi across the split: {d.sum():+.4g} "
         f"(sum of per-taxon curvature changes, exact)",
         f"[{time.time() - t0:.0f}s]"]
    open(os.path.join(R, f"36_attribution_{D}.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L), flush=True)


def ratio(D, rng):
    """C and its permutation p, or None if an input is missing."""
    f = os.path.join(R, f"36_attribution_{D}.csv")
    ind = indicators(D)
    if not os.path.exists(f) or ind is None:
        return None
    a = pd.read_csv(f)
    w = np.abs(a.dkappa.values)
    S = len(w)
    by = {rname(t): i for i, t in enumerate(a.taxon)}
    idx = np.array(sorted({by[rname(t)] for t in ind if rname(t) in by}))
    if len(idx) == 0 or len(idx) == S or w.sum() == 0:
        return dict(D=D, S=S, nI=len(idx), C=np.nan, p=np.nan, top5=np.nan,
                    p5=np.nan)
    share = w[idx].sum() / w.sum()
    C = share / (len(idx) / S)
    draws = np.array([w[rng.choice(S, len(idx), replace=False)].sum() / w.sum()
                      for _ in range(NDRAW)])
    p = (1 + min((draws <= share).sum(), (draws >= share).sum())) / (1 + NDRAW)
    top5 = set(np.argsort(-w)[:5].tolist())
    k = len(top5 & set(idx.tolist()))
    p5 = float(hypergeom.cdf(k, S, len(idx), 5))
    return dict(D=D, S=S, nI=len(idx), C=C, p=min(2 * p, 1.0), top5=k, p5=p5)


def test():
    rng = np.random.default_rng(8)
    rows = [r for r in (ratio(D, rng) for D in NAMES) if r is not None]
    det = {D: detects(D) for D in NAMES}
    L = ["Do the two channels run on different taxa? (amendment 8)",
         "",
         "C = share of |dkappa| carried by TITAN's indicators, divided by",
         "their share of the taxon list.  C < 1: structure is carried by",
         "taxa TITAN does not flag.", "",
         f"{'dataset':>8} {'detect':>7} {'S':>5} {'|I|':>5} {'C':>7} "
         f"{'p':>7} {'top5 in I':>10} {'p(hyp)':>8}"]
    ok = []
    for r in rows:
        d = det.get(r["D"])
        L.append(f"{r['D']:>8} {str(d):>7} {r['S']:>5} {r['nI']:>5} "
                 f"{r['C']:7.3f} {r['p']:7.4f} {r['top5']:10d} {r['p5']:8.4f}")
        if d and np.isfinite(r["C"]) and r["C"] > 0:
            ok.append(r["C"])
    L += ["", f"detecting datasets with both inputs: {len(ok)}"]
    if len(ok) >= 4:
        lc = np.log(np.array(ok))
        stat, p = wilcoxon(lc)
        L += [f"Wilcoxon signed-rank of log C against zero: W = {stat:.1f}, "
              f"p = {p:.4f}",
              f"sign test: {int((np.array(ok) < 1).sum())} of {len(ok)} have "
              f"C < 1",
              f"median C = {np.median(ok):.3f}"]
    else:
        L += ["fewer than four detecting datasets with both inputs: the "
              "confirmatory test is not run, and the C values above stand as "
              "description only."]
    open(os.path.join(R, "36_attribution.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    args = sys.argv[1:]
    if args == ["test"]:
        test()
    else:
        for D in args:
            attribute(D)
