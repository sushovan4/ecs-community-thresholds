"""Fix the threshold estimator BEFORE any Butte Hill data is seen.

The draft paper (sec 2.4) leaves the change statistic open among three
candidates: a discrete derivative along the gradient, a two-sample separation
across a candidate split, and a fitted breakpoint.  This experiment chooses
among them on synthetic gradients with a KNOWN threshold, so that the choice
can be pre-registered rather than made on the data it will be used on.

Community model
---------------
S taxa in NMOD modules; module memberships never change.  Modules are paired,
and above the threshold z* each pair begins sharing a latent factor with
mixing weight A_HI -- the community reorganises from NMOD blocks into NMOD/2
super-blocks.  Every taxon's marginal abundance distribution is IDENTICAL at
every gradient position by construction (unit variance throughout): only the
correlation structure changes.  This is exactly the regime the method claims
-- reorganisation without taxon-level abundance shifts -- and the regime
taxon-aggregation methods (TITAN) are by construction blind to; see
07_titan_comparison.R for the direct check.

Scenarios:  THRESH  mixing 0 below z*, A_HI above (sharp threshold at z*)
            SMOOTH  mixing rises linearly in z from 0 to A_HI (no threshold)
            NULL    mixing constant A_HI/2 (no change along z at all)
            ABUND   no correlation change; module-0 taxa shift their MEAN
                    abundance above z* -- the converse blind spot: TITAN's
                    signal, invisible to a correlation-based summary

Outcome (results/06_recovery*.txt): the two-sample split with the mean over
the connectance grid of |Cohen d| ("split_mean") is unbiased with the best
RMSE at both design sizes, and is the pre-registered choice.  Localisation
precision is set by the window width, not the statistic.

Windows of W consecutive plots (by gradient order) slide with stride; each
window gives chi at matched edge counts (ec.py), so connectance is fixed and
any change along z is structural.

Subcommands:
  recovery   estimator comparison: bias/RMSE of z-hat on THRESH, scatter on
             SMOOTH, statistic magnitudes on NULL             (~1 min)
  bootstrap  CI width of z-hat, THRESH vs SMOOTH -- sharpness diagnostic
  power      permutation power on THRESH and FPR on NULL at three design
             sizes, chosen estimator                          (~1 h)
"""
import os
import sys
import time
import numpy as np
from ec import corr_distance, chi_by_edge_count, cohen_d

S, NMOD = 60, 6                       # taxa, modules (paired: 0-1, 2-3, 4-5)
RHO, A_HI, ZSTAR = 0.55, 0.7, 0.55    # factor loading^2, merge weight, truth
NPERM = int(os.environ.get("NPERM", 99))
NBOOT = int(os.environ.get("NBOOT", 50))
M_GRID = np.unique(np.linspace(40, 900, 12).astype(int))


def community(P, scenario, rng):
    """Abundance matrix S x P and gradient positions z.  Unit marginal
    variance for every taxon at every plot; only correlations move."""
    z = np.sort(rng.uniform(size=P))
    a = {"THRESH": np.where(z > ZSTAR, A_HI, 0.0),
         "SMOOTH": z * A_HI,
         "NULL":   np.full(P, A_HI / 2),
         "ABUND":  np.zeros(P)}[scenario]
    mod = np.repeat(np.arange(NMOD), S // NMOD)
    u = rng.normal(size=(NMOD, P))            # module factors
    v = rng.normal(size=(NMOD // 2, P))       # shared pair factors
    f = np.sqrt(1 - a) * u[mod] + np.sqrt(a) * v[mod // 2]
    X = np.sqrt(RHO) * f + np.sqrt(1 - RHO) * rng.normal(size=(S, P))
    if scenario == "ABUND":
        X[mod == 0] += np.where(z > ZSTAR, 1.0, 0.0)
    return X, z


def window_curves(X, z, W, stride):
    """chi(m) per sliding window of W plots in gradient order."""
    o = np.argsort(z)
    starts = range(0, len(z) - W + 1, stride)
    C = np.array([chi_by_edge_count(corr_distance(X[:, o[s:s + W]]), M_GRID)
                  for s in starts])
    zc = np.array([z[o[s:s + W]].mean() for s in starts])
    return C, zc


# ---- the three candidate estimators ---------------------------------------
# Each returns (z-hat, change statistic).  The statistic doubles as the
# permutation-test statistic, so estimator and test are one object.

def est_deriv(C, zc):
    sc = C.std(0, ddof=1) + 1e-12
    g = np.abs(np.diff(C, axis=0) / sc).max(1)
    i = int(g.argmax())
    return (zc[i] + zc[i + 1]) / 2, float(g.max())


def est_split(C, zc):
    best, bz = -np.inf, np.nan
    for s in range(2, len(zc) - 1):
        d = np.abs(cohen_d(C[:s], C[s:])).max()
        if d > best:
            best, bz = d, (zc[s - 1] + zc[s]) / 2
    return bz, float(best)


def est_break(C, zc):
    """Continuous two-segment fit, shared break across m; statistic is the
    pooled fraction of variance the break explains beyond a single line."""
    sc = C.std(0, ddof=1) + 1e-12
    Y = (C - C.mean(0)) / sc
    A0 = np.c_[np.ones_like(zc), zc]
    sse0 = (np.linalg.lstsq(A0, Y, rcond=None)[1]).sum()
    best, bz = np.inf, np.nan
    for zb in zc[1:-1]:
        A = np.c_[A0, np.maximum(zc - zb, 0.0)]
        r = Y - A @ np.linalg.lstsq(A, Y, rcond=None)[0]
        sse = float((r ** 2).sum())
        if sse < best:
            best, bz = sse, zb
    return bz, float((sse0 - best) / (sse0 + 1e-12))


def est_split_mean(C, zc):
    """Two-sample split, |Cohen d| POOLED over the connectance grid rather
    than maximised -- pooling stabilises the location."""
    best, bz = -np.inf, np.nan
    for s in range(2, len(zc) - 1):
        d = np.abs(cohen_d(C[:s], C[s:])).mean()
        if d > best:
            best, bz = d, (zc[s - 1] + zc[s]) / 2
    return bz, float(best)


ESTIMATORS = {"deriv": est_deriv, "split": est_split, "break": est_break,
              "split_mean": est_split_mean}
CHOSEN = "split_mean"                 # chosen by the recovery run; see results


def design(P):
    W = max(6, int(round(P / 4)))
    return W, max(1, W // 3)


def perm_p(X, z, W, stride, est, rng):
    _, obs = est(*window_curves(X, z, W, stride))
    hits = 0
    for _ in range(NPERM):
        _, s = est(*window_curves(X, rng.permutation(z), W, stride))
        hits += s >= obs
    return (1 + hits) / (1 + NPERM)


def run_recovery(nreal=40, P=40):
    W, stride = design(P)
    print(f"S={S} modules={NMOD} rho={RHO} a_hi={A_HI} z*={ZSTAR}")
    print(f"P={P} plots, windows of {W} stride {stride}, "
          f"m in [{M_GRID[0]}, {M_GRID[-1]}], {nreal} realisations\n")
    for scen in ("THRESH", "SMOOTH", "NULL"):
        zh = {k: [] for k in ESTIMATORS}
        st = {k: [] for k in ESTIMATORS}
        for i in range(nreal):
            C, zc = window_curves(*community(P, scen, np.random.default_rng(600 + i)),
                                  W, stride)
            for k, est in ESTIMATORS.items():
                z_, s_ = est(C, zc)
                zh[k].append(z_)
                st[k].append(s_)
        print(f"-- {scen}" + (f"  (true z* = {ZSTAR})" if scen == "THRESH" else ""))
        for k in ESTIMATORS:
            z_, s_ = np.array(zh[k]), np.array(st[k])
            line = f"   {k:>10}  stat {s_.mean():6.2f}"
            if scen == "THRESH":
                line += (f"   bias {z_.mean() - ZSTAR:+.3f}"
                         f"   RMSE {np.sqrt(((z_ - ZSTAR) ** 2).mean()):.3f}"
                         f"   z-hat sd {z_.std(ddof=1):.3f}")
            else:
                line += f"   z-hat sd {z_.std(ddof=1):.3f}"
            print(line)
        print()


def run_bootstrap(nreal=15, P=80):
    W, stride = design(P)
    est = ESTIMATORS[CHOSEN]
    print(f"bootstrap CI width of z-hat ({CHOSEN}), {NBOOT} resamples, "
          f"{nreal} realisations, P={P}\n")
    for scen in ("THRESH", "SMOOTH"):
        widths = []
        for i in range(nreal):
            X, z = community(P, scen, np.random.default_rng(1600 + i))
            o = np.argsort(z)
            rng = np.random.default_rng(2600 + i)
            bhat = []
            for _ in range(NBOOT):
                starts = range(0, P - W + 1, stride)
                C = np.array([chi_by_edge_count(
                    corr_distance(X[:, rng.choice(o[s:s + W], W, replace=True)]),
                    M_GRID) for s in starts])
                zc = np.array([z[o[s:s + W]].mean() for s in starts])
                bhat.append(est(C, zc)[0])
            lo, hi = np.percentile(bhat, [2.5, 97.5])
            widths.append(hi - lo)
        w = np.array(widths)
        print(f"   {scen:>6}  CI width mean {w.mean():.3f}  "
              f"range {w.min():.3f}-{w.max():.3f}")
    print()


def run_power(nreal=20):
    est = ESTIMATORS[CHOSEN]
    print(f"permutation test ({CHOSEN} statistic), {NPERM} shuffles of plot "
          f"positions, {nreal} realisations per cell, alpha = 0.05\n")
    print(f"{'plots':>6} {'window':>7} {'scenario':>9} {'rejection rate':>15}")
    for P in (24, 40, 80):
        W, stride = design(P)
        scens = ("THRESH", "NULL", "ABUND") if P == 80 else ("THRESH", "NULL")
        for scen in scens:
            ps = [perm_p(*community(P, scen, np.random.default_rng(3600 + P + i)),
                         W, stride, est, np.random.default_rng(4600 + P + i))
                  for i in range(nreal)]
            print(f"{P:6d} {W:7d} {scen:>9} {np.mean(np.array(ps) <= 0.05):15.2f}",
                  flush=True)


if __name__ == "__main__":
    t0 = time.time()
    what = sys.argv[1] if len(sys.argv) > 1 else "recovery"
    if what in ("recovery", "all"):
        run_recovery(nreal=40, P=40)
        run_recovery(nreal=30, P=80)
    if what in ("bootstrap", "all", "slow"):
        run_bootstrap()
    if what in ("power", "all", "slow"):
        run_power()
    print(f"[total {time.time() - t0:.0f}s]")
