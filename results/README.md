# Results log

## 03_arabidopsis — 2026-08-23

```
150 genes x 690 accessions matched to FT16
observed  max |Cohen d| = 1.61 at 569 edges
permuted  max |Cohen d| = 0.85 mean, 0.20-2.63 range over 199 label shuffles
empirical p = 0.040
```

**Read this as a hint, not a result.**

- p = 0.040 means 7 of 199 label shuffles beat the observed value. The null is
  heavy-tailed — mean 0.85 but reaching 2.63, well above the observed 1.61 — so
  the statistic is unstable and a different seed could move this materially.
- Nothing was pre-registered. `NGENE=150`, `NACC=90`, `BOOT=12`, `MAXDIM=3` and
  the edge grid were all chosen by hand. Other settings could easily move p.
- **The dataset cannot properly test the claim.** GSE80744 is unstressed rosette
  tissue at a single condition, and flowering time is a developmental trait. The
  method is about response across a stress gradient, so there is no second axis
  here and no surface — only a single chi curve. What ran is a degenerate slice
  of the actual construction.

**What this does establish:** the pipeline runs end to end on real 1001 Genomes
data at a workable cost, the matched-density design behaves as intended, and the
permutation control is in place. That is the scaffold; the test still needs a
dataset with a real environmental axis.

**Next:** switchgrass `pvdiv` — ten common gardens on a latitudinal gradient, so
ECS(edge count, site) is an actual surface rather than one slice.

---

## 04_switchgrass_drought — 2026-08-23

```
observed max |Cohen d|, control vs drought:  5AM 0.51 · 10AM 2.26 · 12PM 0.30 · 2PM 0.27
                                             spread across the axis = 1.99

null, 199 treatment-label shuffles within time point (cell sizes preserved):
   max|d|   mean 1.61, range 0.59–4.26   ->  p = 0.135
   spread   mean 1.16, range 0.08–3.61   ->  p = 0.115
```

**Negative. The 10AM peak did not survive.**

The null routinely reaches and exceeds the observed value, so the peak is
consistent with cell-size noise — the confound the design was built to test,
and it fired. Both questions come back negative: drought is not separably
detected at all (p = 0.135), and there is no evidence it reorganizes along the
axis (p = 0.115).

**Read it as "not detected", not "shown absent".** The test is badly
underpowered: cells of 6–7 samples, bootstrap replicates within them, and a
max-statistic over four time points. A null ranging 0.59–4.26 is the signature
of a very noisy statistic, not of a well-estimated zero.

**What it costs.** The localization claim — identified in
[verdict.md](../docs/verdict.md) as the one surviving distinctive
contribution — now has **no real-data support**. It stands on synthetic
geometry alone (Cohen d = 4.71 at matched density).

**What it does not cost.** Neither public dataset was ever a fair test.
GSE80744 has no environmental axis at all; GSE57887's axis is diel, not a
stress dose, with six samples a cell. The construction has not been tested on
what it was designed for: a real dose gradient with adequate replication.
**Butte Hill is now the first real test of the claim, not a confirmation of it.**

---

## 06_gradient_estimator — 2026-09-01

The estimator-selection run that fixes §2.4 of the paper, on simulated
gradients with a known threshold (z* = 0.55): community of 60 taxa in 6
modules whose correlation structure merges pairwise above z*, every taxon's
marginal abundance distribution held constant along the gradient by
construction.

```
P=80 (windows of 20):                RMSE     bias    stat THRESH/SMOOTH/NULL
  derivative along z                 0.227   -0.011      3.15 / 3.09 / 3.12
  split, max over m                  0.177   -0.003      5.02 / 3.65 / 2.91
  two-segment fit                    0.205   -0.101      0.30 / 0.32 / 0.28
  split, mean over m   << CHOSEN >>  0.152   -0.003      2.23 / 1.74 / 1.36
```

**Pre-registered choice: the pooled split (`split_mean`).** Unbiased at both
design sizes, best RMSE at the Butte-scale design, and the only statistic whose
signal ordering THRESH > SMOOTH > NULL holds at both sizes. The max-over-grid
statistic — what experiments 03 and 04 used — is confirmed unstable at small n.
At P=40 **no** statistic separates threshold from null; detection at that size
fails, which is the switchgrass lesson made quantitative. Localization RMSE
tracks half the window width at both sizes: resolution is a survey-design
property, not an estimator property.

Full log: `06_estimator_recovery.txt`. Power/FPR/sharpness:
`06_estimator_slow.txt`. TITAN complementarity: `07_titan.txt`.

---

## 07_titan_comparison — 2026-09-01

TITAN2 on five realizations each of the two simulated gradient types at P=80
(true threshold 0.55; 250 permutations, 250 bootstraps, purity/reliability
0.95). See `07_titan.txt`.

```
reorganization gradients:  1, 5, 9, 4, 1 of 60 taxa flagged; per-run median
                           change points 0.36, 0.85, 0.26, 0.75, 0.83;
                           sum(z) intervals span most of the gradient
abundance-shift gradients: 11, 11, 13, 12, 11 of 60 flagged, centered on the
                           shifted module; median indicator change points
                           0.52, 0.49, 0.52, 0.55, 0.53
```

**TITAN finds no coherent threshold when the community reorganizes without
taxon-level abundance change — scattered flags at roughly the false-positive
rate its cutoffs admit — and recovers a pure abundance-shift threshold
cleanly in every run.** The ECS side of the two-way table comes from the
ABUND row of the power run in `06_estimator_slow.txt`.

---

## 06 power/FPR/sharpness (local, reduced strength) — 2026-09-01

```
power (THRESH) / FPR (NULL) at alpha=0.05, 99 perms, 20 realizations:
   P=24: 0.15 / 0.00      P=40: 0.05 / 0.05      P=80: 0.20 / 0.00
   ABUND converse at P=80: 0.00  (the ECS side of the TITAN 2x2)
sharpness: bootstrap CI width THRESH 0.525 vs SMOOTH 0.527 -> diagnostic FAILS
```

Error control sound, selectivity clean, absolute power modest at survey
scale — single-survey designs under ~100 plots are underpowered. The
sharpness negative is reported in the paper as such. Publication-strength
refresh + P=160: Pegasus job 73622978.

---

## 11_design_extras + publication-strength power — 2026-09-01

```
CI coverage (P=80, THRESH): 20/20 — conservative, mean width 0.49
multi-year power at P=80:  averaging T=1/3/5/10 -> 0.33/0.33/0.27/0.27
                           stacking  T=3/5/10   -> 0.47/0.73/0.53
power at 499 perms, 50 reals: 0.00/0.08/0.18/0.66 at P=24/40/80/160
                              FPR at nominal; ABUND 0.00 (P=80), 0.04 (P=160)
```

**Averaging repeat surveys buys nothing; stacking them roughly doubles power
by T=5** (binding noise is correlation sampling across plots — averaging
cleans the wrong term). **P=160 single-survey power is 0.66.** Protocol
amended (paper §5): stack plot-years, permute whole plots; 5+ survey years
puts a Butte-scale design in the fair-test regime.

---

## 09_glades — 2026-09-01

The frozen estimator on TITAN's own benchmark (126 Everglades sites x 164
macroinvertebrate taxa, TP gradient 2.5-169.5 ug/L):

```
chi: z-hat = 12.5 ug/L TP   p = 0.300 (199 perms)   95% CI [12.0, 45.3]
 rq: z-hat = 16.9 ug/L TP   p = 0.365               95% CI [11.9, 44.0]
TITAN (regenerated, 09_glades_titan.txt): decliners cp 15.1 [13.2, 19.9],
increasers 32.5 [21.7, 34.1], 87/164 pure & reliable taxa
```

**Point estimate lands inside TITAN's decliner interval and beside the
10 ug/L legal criterion; the pre-registered test does not reject, and the
paper reports no detection.** Consistent with the power table: 126
single-survey sites sit between the 0.18 (P=80) and 0.66 (P=160) power rows.
The run bounds the signal rather than refuting it, and demonstrates the
protocol discipline Butte Hill gets.

---

## Publication-strength summary swap (driver, 499 perms) — 2026-09-02

```
P=160 rejection: STRENGTH -> rq 0.34, chi 0.12 (blindness by design, confirmed)
                 REWIRE   -> modul 0.18, chi 0.08; statistics separate, no
                             summary clears alpha reliably (hardest signal)
P=80: all summaries at false-alarm level on every structured scenario
modularity edges chi on merge (stat 4.28 vs 3.13; RMSE 0.074 vs 0.102)
```

Both directions of the distribution/shape asymmetry confirmed at
publication strength; arrangement-only reorganization exceeds dependable
detection even at 160 single-survey plots; the transition-window leak is
visible (meand 0.14 on REWIRE). Reported in paper Table 3 with the
chi-vs-modularity comparison stated plainly.

---

## 15_glades_anchored (EXPLORATORY) — 2026-09-02

TITAN-anchored single-look test on glades, run AFTER the primary p=0.30 was
seen (labeled exploratory; the pre-registered version lives in the Butte
protocol): D = 0.99 at the boundary nearest TITAN's decliner cp (window
boundary 16.9 vs cp 15.1), p = 0.375. The glades non-detection is not the
free-split multiplicity tax — the covariance channel is genuinely
unresolvable at 126 single-survey sites.
