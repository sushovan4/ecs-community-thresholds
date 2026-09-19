"""Publication-strength simulation runs for paper sections 3.3-3.5.

Run from the experiments/ directory (the sbatch does):
    python3 ../hpc/driver_full.py
Honours NPERM / NBOOT from the environment (set before import, the modules
read them at import time) and SLURM_CPUS_PER_TASK for the pool size.
Writes results/<name>_pegasus.txt, mirroring the local reduced-strength logs.
"""
import importlib
import os
import sys
import time
from contextlib import redirect_stdout
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "experiments"))
sim = importlib.import_module("06_gradient_estimator")
swap = importlib.import_module("08_summary_swap")

NCPU = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 8))
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results")
os.makedirs(OUT, exist_ok=True)


def power_cell(args):
    P, scen, i = args
    W, stride = sim.design(P)
    est = sim.ESTIMATORS[sim.CHOSEN]
    p = sim.perm_p(*sim.community(P, scen, np.random.default_rng(3600 + P + i)),
                   W, stride, est, np.random.default_rng(4600 + P + i))
    return P, scen, p


def boot_cell(args):
    scen, i, P = args
    W, stride = sim.design(P)
    est = sim.ESTIMATORS[sim.CHOSEN]
    X, z = sim.community(P, scen, np.random.default_rng(1600 + i))
    o = np.argsort(z)
    rng = np.random.default_rng(2600 + i)
    bhat = []
    for _ in range(sim.NBOOT):
        starts = range(0, P - W + 1, stride)
        C = np.array([sim.chi_by_edge_count(
            sim.corr_distance(X[:, rng.choice(o[s:s + W], W, replace=True)]),
            sim.M_GRID) for s in starts])
        zc = np.array([z[o[s:s + W]].mean() for s in starts])
        bhat.append(est(C, zc)[0])
    lo, hi = np.percentile(bhat, [2.5, 97.5])
    return scen, P, hi - lo


def main():
    t0 = time.time()
    print(f"driver: NCPU={NCPU} NPERM={sim.NPERM} NBOOT={sim.NBOOT}", flush=True)

    with open(os.path.join(OUT, "06_estimator_recovery_pegasus.txt"), "w") as f, \
            redirect_stdout(f):
        for P in (40, 80, 160):
            sim.run_recovery(nreal=100, P=P)
    print(f"recovery done [{time.time() - t0:.0f}s]", flush=True)

    nreal = 50
    with Pool(NCPU) as pool:
        jobs = [(scen, i, P) for P in (80, 160)
                for scen in ("THRESH", "SMOOTH") for i in range(nreal)]
        boots = pool.map(boot_cell, jobs)
    with open(os.path.join(OUT, "06_bootstrap_pegasus.txt"), "w") as f, \
            redirect_stdout(f):
        print(f"bootstrap CI width of z-hat ({sim.CHOSEN}), {sim.NBOOT} "
              f"resamples, {nreal} realisations per cell\n")
        for P in (80, 160):
            for scen in ("THRESH", "SMOOTH"):
                w = np.array([b[2] for b in boots if b[:2] == (scen, P)])
                print(f"   P={P:3d}  {scen:>6}  CI width mean {w.mean():.3f}  "
                      f"range {w.min():.3f}-{w.max():.3f}")
    print(f"bootstrap done [{time.time() - t0:.0f}s]", flush=True)

    with Pool(NCPU) as pool:
        jobs = [(P, scen, i) for P in (24, 40, 80, 160)
                for scen in (("THRESH", "NULL", "ABUND") if P >= 80
                             else ("THRESH", "NULL"))
                for i in range(nreal)]
        cells = pool.map(power_cell, jobs)
    with open(os.path.join(OUT, "06_power_pegasus.txt"), "w") as f, \
            redirect_stdout(f):
        print(f"permutation test ({sim.CHOSEN}), {sim.NPERM} shuffles, "
              f"{nreal} realisations per cell, alpha = 0.05\n")
        print(f"{'plots':>6} {'window':>7} {'scenario':>9} {'rejection rate':>15}")
        for P in (24, 40, 80, 160):
            W, _ = sim.design(P)
            scens = ("THRESH", "NULL", "ABUND") if P >= 80 else ("THRESH", "NULL")
            for scen in scens:
                ps = np.array([c[2] for c in cells if c[:2] == (P, scen)])
                print(f"{P:6d} {W:7d} {scen:>9} {np.mean(ps <= 0.05):15.2f}")
    print(f"power done [{time.time() - t0:.0f}s]", flush=True)

    with open(os.path.join(OUT, "08_summary_recovery_pegasus.txt"), "w") as f, \
            redirect_stdout(f):
        swap.run_recovery(nreal=100, P=80)
        swap.run_recovery(nreal=100, P=160)
    print(f"swap recovery done [{time.time() - t0:.0f}s]", flush=True)

    keys = ("chi", "betti0", "betti1", "modul", "rq", "meand")
    scens = ("REWIRE", "STRENGTH", "TRI_NULL")
    with open(os.path.join(OUT, "08_summary_rejection_pegasus.txt"), "w") as f, \
            redirect_stdout(f):
        print(f"permutation rejection at alpha=0.05, {sim.NPERM} shuffles, "
              f"{nreal} realisations per cell\n")
        for P in (80, 160):
            W, stride = sim.design(P)
            with Pool(NCPU) as pool:
                jobs = [(sc, i, P, W, stride)
                        for sc in scens for i in range(nreal)]
                res = pool.map(swap._one_rejection, jobs)
            print(f"-- P={P}, windows of {W}")
            print(f"{'summary':>8}" + "".join(f"{sc:>10}" for sc in scens))
            for k in keys:
                row = f"{k:>8}"
                for sc in scens:
                    ps = np.array([r[1][k] for r in res if r[0] == sc])
                    row += f"{np.mean(ps <= 0.05):10.2f}"
                print(row)
            f.flush()
    print(f"all done [{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    sys.exit(main())
