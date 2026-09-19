"""Switchgrass drought (GEO GSE57887, Panicum virgatum Alamo AP13, 119 samples).

The first test with a genuine ORDERED ENVIRONMENTAL AXIS.  Design at day 14:
treatment (control / drought / recovery) x time of day (5AM, 10AM, 12PM, 2PM).

Meyer et al. report that drought-responsive expression reverses sign across the
diel cycle -- up pre-dawn, down mid-day.  So the question is not just "does
drought differ from control" but "does the difference REORGANIZE along the
axis" -- i.e. is the surface informative where a single curve would not be.

Genes are points; coordinates are the replicate samples within one cell.
Sweep indexed by edge count, so density is fixed.  Treatment-label permutation
within each time point is the null.
"""
import os, re, glob, gzip, time
import numpy as np
import pandas as pd
from ec import corr_distance, chi_by_edge_count, cohen_d

NGENE = int(os.environ.get("NGENE", 150))
NREP  = int(os.environ.get("NREP", 6))
BOOT  = int(os.environ.get("BOOT", 15))
NPERM = int(os.environ.get("NPERM", 199))
ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA  = os.path.join(ROOT, "data")

# ---- metadata --------------------------------------------------------------
txt = open(os.path.join(DATA, "gse57887_samples.txt")).read()
blocks = txt.split("^SAMPLE = ")[1:]
meta = []
for b in blocks:
    gsm = b.split("\n", 1)[0].strip()
    def g(pat):
        m = re.search(pat, b)
        return m.group(1).strip() if m else None
    meta.append(dict(gsm=gsm,
                     treatment=g(r"treatment:\s*(.+)"),
                     day=g(r"sampling day:\s*(.+)"),
                     tod=g(r"sampling time:\s*(.+)")))
meta = pd.DataFrame(meta)

# ---- counts ----------------------------------------------------------------
cache = os.path.join(DATA, "gse57887_counts.pkl.gz")
if os.path.exists(cache):
    counts = pd.read_pickle(cache)
else:
    cols = {}
    for f in sorted(glob.glob(os.path.join(DATA, "gse57887", "GSM*.txt.gz"))):
        gsm = os.path.basename(f).split("_")[0]
        cols[gsm] = pd.read_csv(f, sep="\t", index_col=0, header=0).iloc[:, 0]
    counts = pd.DataFrame(cols).fillna(0.0)
    counts.to_pickle(cache)
meta = meta[meta.gsm.isin(counts.columns)].reset_index(drop=True)

cpm = counts / counts.sum(0) * 1e6
X = np.log2(cpm.values + 1.0)
keep = np.argsort(X.var(1))[-NGENE:]                 # label-blind selection
X = X[keep]
print(f"{X.shape[0]} genes x {X.shape[1]} samples   (from {counts.shape[0]} isogroups)")

d14 = meta[(meta.day == "14 d")]
TOD = ["5AM", "10AM", "12PM", "2PM"]
print("\ncell sizes at day 14:")
print(pd.crosstab(d14.treatment, d14.tod).reindex(columns=TOD).to_string(), "\n")

gsm_ix = {g: i for i, g in enumerate(counts.columns)}
m_grid = np.unique(np.linspace(200, 2600, 14).astype(int))


def cell(treat, tod):
    return [gsm_ix[g] for g in d14[(d14.treatment == treat) & (d14.tod == tod)].gsm]


def replicates(idx, rng):
    return np.array([chi_by_edge_count(
        corr_distance(X[:, rng.choice(idx, NREP, replace=True)]), m_grid)
        for _ in range(BOOT)])


def stat(assign, rng):
    """max |Cohen d| between the two arms, per time of day."""
    out = []
    for t in TOD:
        a, b = assign[t]
        if len(a) < 3 or len(b) < 3:
            out.append(np.nan); continue
        out.append(np.abs(cohen_d(replicates(a, rng), replicates(b, rng))).max())
    return np.array(out)


obs_assign = {t: (cell("control", t), cell("drought", t)) for t in TOD}
t0 = time.time()
obs = stat(obs_assign, np.random.default_rng(0))
print("observed max |Cohen d|, control vs drought, by time of day:")
for t, v in zip(TOD, obs):
    print(f"   {t:>5s}  {v:5.2f}")
print(f"   range across the axis = {np.nanmax(obs)-np.nanmin(obs):.2f}   [{time.time()-t0:.0f}s]\n")

# null: shuffle treatment labels WITHIN each time point, preserving cell sizes
null_range, null_max = [], []
for p in range(NPERM):
    r = np.random.default_rng(500 + p)
    asg = {}
    for t in TOD:
        a, b = obs_assign[t]
        pool = r.permutation(a + b)
        asg[t] = (list(pool[:len(a)]), list(pool[len(a):]))
    v = stat(asg, r)
    null_range.append(np.nanmax(v) - np.nanmin(v)); null_max.append(np.nanmax(v))
null_range, null_max = np.array(null_range), np.array(null_max)
if NPERM == 0:
    raise SystemExit("NPERM=0: observed statistic only, no null computed.")

p_rng = (1 + (null_range >= (np.nanmax(obs)-np.nanmin(obs))).sum()) / (1 + NPERM)
p_max = (1 + (null_max >= np.nanmax(obs)).sum()) / (1 + NPERM)
print(f"null over {NPERM} treatment-label shuffles within time point:")
print(f"   max|d|  mean {null_max.mean():.2f}, range {null_max.min():.2f}-{null_max.max():.2f}"
      f"   -> p = {p_max:.3f}   (is drought separable at all?)")
print(f"   spread  mean {null_range.mean():.2f}, range {null_range.min():.2f}-{null_range.max():.2f}"
      f"   -> p = {p_rng:.3f}   (does it REORGANIZE along the axis?)")
print(f"\n[total {time.time()-t0:.0f}s]")
