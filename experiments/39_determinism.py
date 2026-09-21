"""Is the interaction statistic reproducible?  Same data, same seeds.

D6 returned T = 10.0932 on one machine and 9.3814 on another, with the
located split moving from 17 to 0.  Nothing in the pipeline draws random
numbers without a fixed seed, so either a linear-algebra backend or the
Alpha complex is the source.  This probe reports, for a dataset:

  * T at every candidate split (no permutations), to 10 decimals;
  * the same after forcing GUDHI's exact predicates;
  * the gap between the best and second-best split, which says whether the
    located split is even identifiable.

Run it on two machines and diff the output.

  python3 39_determinism.py D6 D1 ELT
Output: results/39_determinism_<name>_<host>.txt
"""
import os

for _v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import importlib
import platform
import socket
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
mix = importlib.import_module("26_mixup_split")
ds = importlib.import_module("28_mixup_datasets")
multi = importlib.import_module("19_multi")
lt = importlib.import_module("17_everglades_lt")

import gudhi  # noqa: E402
from intersection_ecp import core as _core  # noqa: E402


def load(name):
    if name == "ELT":
        dens, keep, z_psu, psu = lt.build()
        return np.log1p(dens.values.T), z_psu.loc[psu].values
    w, zc, reg_s, zname = multi.BUILD[name]()
    X, z, groups, units, zunit, keep, ncand, wide = multi.prepare(w, zc)
    return X, z


def profile(X, z):
    o = np.argsort(z, kind="stable")
    out = []
    for f in mix.FRACS:
        cut = max(3, int(round(f * len(z))))
        T = mix.stat_at(mix.corr_profiles(X[:, o[:cut]]),
                        mix.corr_profiles(X[:, o[cut:]]))
        out.append((float(f), float(T),
                    float((z[o[cut - 1]] + z[o[cut]]) / 2)))
    return out


def report(name):
    X, z = load(name)
    L = [f"{name} on {socket.gethostname()} ({platform.platform()})",
         f"numpy {np.__version__}, gudhi {gudhi.__version__}, "
         f"python {sys.version.split()[0]}",
         f"X: {X.shape}, checksum {np.sum(X):.12e}", ""]
    for tag, prec in (("fast", "fast"), ("safe", "safe"), ("exact", "exact")):
        _core._PRECISION = prec
        rows = profile(X, z)
        Ts = np.array([r[1] for r in rows])
        k = int(np.argmin(Ts))
        srt = np.sort(Ts)
        L.append(f"-- precision={tag}")
        for f, T, zc in rows:
            L.append(f"   frac {f:.2f}  T {T:.10f}  z {zc:.6g}"
                     + ("   <- min" if abs(T - Ts[k]) < 1e-12 else ""))
        L.append(f"   best {srt[0]:.10f}  second {srt[1]:.10f}  "
                 f"gap {srt[1] - srt[0]:.10f}  "
                 f"({100 * (srt[1] - srt[0]) / srt[0]:.3f}% of best)")
        L.append("")
    out = os.path.join(ROOT, "results",
                       f"39_determinism_{name}_{socket.gethostname().split('.')[0]}.txt")
    open(out, "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    for n in sys.argv[1:]:
        report(n)
