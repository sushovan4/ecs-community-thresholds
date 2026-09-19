"""Two design questions the review will ask, answered on the 06 gradients.

coverage   Does the 95% bootstrap percentile interval on z-hat actually
           cover the true threshold ~95% of the time?  (The width was
           measured in 06; coverage never was.)

multiyear  The Butte Hill network is ~80 plots but ~10 annual surveys.
           If each plot contributes T repeat surveys sharing the plot's
           latent community state (a fraction TAU of factor variance is
           year-specific and averages out), plot-level averaging shrinks
           the noise that made P=80 power weak (0.20 locally).  Power vs
           T in {1, 3, 5, 10} at P=80, THRESH scenario -- the number that
           decides whether the Butte Hill protocol is a fair test.
"""
import importlib
import os
import sys
import time

import numpy as np

sim = importlib.import_module("06_gradient_estimator")

TAU = 0.3
NPERM = sim.NPERM


def community_years(P, T, scenario, rng):
    """The 06 merge model with T annual surveys per plot, averaged."""
    z = np.sort(rng.uniform(size=P))
    a = {"THRESH": np.where(z > sim.ZSTAR, sim.A_HI, 0.0),
         "NULL": np.full(P, sim.A_HI / 2)}[scenario]
    NMOD = 6
    mod = np.repeat(np.arange(NMOD), sim.S // NMOD)
    u_p = rng.normal(size=(NMOD, P))
    v_p = rng.normal(size=(NMOD // 2, P))
    f_plot = np.sqrt(1 - a) * u_p[mod] + np.sqrt(a) * v_p[mod // 2]
    years = []
    for _ in range(T):
        u_y = rng.normal(size=(NMOD, P))
        v_y = rng.normal(size=(NMOD // 2, P))
        f_year = np.sqrt(1 - a) * u_y[mod] + np.sqrt(a) * v_y[mod // 2]
        f = np.sqrt(1 - TAU) * f_plot + np.sqrt(TAU) * f_year
        years.append(np.sqrt(sim.RHO) * f
                     + np.sqrt(1 - sim.RHO) * rng.normal(size=(sim.S, P)))
    return np.mean(years, axis=0), z


def run_coverage(nreal=20, P=80, nboot=99):
    W, stride = sim.design(P)
    est = sim.ESTIMATORS[sim.CHOSEN]
    hits, widths = 0, []
    for i in range(nreal):
        X, z = sim.community(P, "THRESH", np.random.default_rng(31600 + i))
        o = np.argsort(z)
        rng = np.random.default_rng(32600 + i)
        bhat = []
        for _ in range(nboot):
            starts = range(0, P - W + 1, stride)
            C = np.array([sim.chi_by_edge_count(
                sim.corr_distance(X[:, rng.choice(o[s:s + W], W, replace=True)]),
                sim.M_GRID) for s in starts])
            zc = np.array([z[o[s:s + W]].mean() for s in starts])
            bhat.append(est(C, zc)[0])
        lo, hi = np.percentile(bhat, [2.5, 97.5])
        hits += lo <= sim.ZSTAR <= hi
        widths.append(hi - lo)
    print(f"coverage of the 95% CI at P={P}, THRESH, {nreal} realisations, "
          f"{nboot} boots: {hits / nreal:.2f}   "
          f"(width mean {np.mean(widths):.3f})", flush=True)


def community_years_stacked(P, T, scenario, rng):
    """As community_years, but plot-years become separate columns (stacked),
    so within-window correlation is estimated over W*T samples rather than W.
    Windows must move whole plots, so columns are ordered plot-major."""
    z = np.sort(rng.uniform(size=P))
    a = {"THRESH": np.where(z > sim.ZSTAR, sim.A_HI, 0.0),
         "NULL": np.full(P, sim.A_HI / 2)}[scenario]
    NMOD = 6
    mod = np.repeat(np.arange(NMOD), sim.S // NMOD)
    u_p = rng.normal(size=(NMOD, P))
    v_p = rng.normal(size=(NMOD // 2, P))
    f_plot = np.sqrt(1 - a) * u_p[mod] + np.sqrt(a) * v_p[mod // 2]
    cols, zc = [], []
    for pth in range(P):
        for _ in range(T):
            u_y = rng.normal(size=NMOD)
            v_y = rng.normal(size=NMOD // 2)
            f_y = (np.sqrt(1 - a[pth]) * u_y[mod]
                   + np.sqrt(a[pth]) * v_y[mod // 2])
            f = np.sqrt(1 - TAU) * f_plot[:, pth] + np.sqrt(TAU) * f_y
            cols.append(np.sqrt(sim.RHO) * f
                        + np.sqrt(1 - sim.RHO) * rng.normal(size=sim.S))
            zc.append(z[pth])
    return np.column_stack(cols), np.array(zc)


def run_multiyear_stacked(nreal=15, P=80):
    """Permutation must move whole plots (years within a plot are dependent),
    so shuffle plot z-positions and propagate to that plot's T columns."""
    est = sim.ESTIMATORS[sim.CHOSEN]
    print(f"\nSTACKED multi-year power at P={P} (tau={TAU}, {NPERM} perms, "
          f"{nreal} realisations, alpha=0.05)\n")
    print(f"{'years':>6} {'power (THRESH)':>15} {'FPR (NULL)':>12}")
    for T in (3, 5, 10):
        W = max(6, int(round(P / 4))) * T          # W plots -> W*T columns
        stride = max(1, W // 3)
        rej = {}
        for scen in ("THRESH", "NULL"):
            ps = []
            for i in range(nreal):
                X, zc = community_years_stacked(
                    P, T, scen, np.random.default_rng(43600 + 100 * T + i))
                r = np.random.default_rng(44600 + 100 * T + i)
                _, obs = est(*sim.window_curves(X, zc, W, stride))
                hits = 0
                for _ in range(NPERM):
                    zp = r.permutation(zc.reshape(P, T)) .reshape(-1)
                    _, sb = est(*sim.window_curves(X, zp, W, stride))
                    hits += sb >= obs
                ps.append((1 + hits) / (1 + NPERM))
            rej[scen] = np.mean(np.array(ps) <= 0.05)
        print(f"{T:6d} {rej['THRESH']:15.2f} {rej['NULL']:12.2f}", flush=True)


def run_multiyear(nreal=15, P=80):
    W, stride = sim.design(P)
    est = sim.ESTIMATORS[sim.CHOSEN]
    print(f"multi-year power at P={P} (tau={TAU}, {NPERM} perms, "
          f"{nreal} realisations, alpha=0.05)\n")
    print(f"{'years':>6} {'power (THRESH)':>15} {'FPR (NULL)':>12}")
    for T in (1, 3, 5, 10):
        rej = {}
        for scen in ("THRESH", "NULL"):
            ps = []
            for i in range(nreal):
                X, z = community_years(P, T, scen,
                                       np.random.default_rng(33600 + 100 * T + i))
                ps.append(sim.perm_p(X, z, W, stride, est,
                                     np.random.default_rng(34600 + 100 * T + i)))
            rej[scen] = np.mean(np.array(ps) <= 0.05)
        print(f"{T:6d} {rej['THRESH']:15.2f} {rej['NULL']:12.2f}", flush=True)


if __name__ == "__main__":
    t0 = time.time()
    what = sys.argv[1] if len(sys.argv) > 1 else "coverage"
    if what in ("coverage", "all"):
        run_coverage()
    if what in ("multiyear", "all"):
        run_multiyear()
    if what in ("stacked",):
        run_multiyear_stacked()
    print(f"[total {time.time() - t0:.0f}s]")
