# Pre-registration: Everglades long-term insect panel — 2026-09-02

Frozen BEFORE downloading or inspecting any data file. Dataset: Pintar &
Dorn, Dryad doi:10.5061/dryad.wstqjq32d (Ecological Entomology,
doi:10.1111/een.70122). Committed prior to first download; the commit
timestamp is the registration.

**Purpose.** The first real-data test in the POWERED regime: the spatially
expansive panel (~146 sites, annual sampling 2005–2024) sits near the P=160
single-survey power point (0.66), with ~20 years available to stack.

**Frozen choices, in order of application:**

1. **Panel**: the spatially expansive annual panel (~146 sites) is primary.
   The 24-site × 5/year panel is a robustness set, analyzed only after the
   primary result is recorded.
2. **Sampling unit → column**: mean count per throw-trap within each
   site-year (density), then log(1+x). One column per site-year.
3. **Taxa**: all insect taxa provided, keeping those present in ≥5 sites.
4. **Gradient z (site-level scalar)**: the site's long-term hydrologic
   position, taken as the first available of, in this fixed preference
   order: (a) hydroperiod (days inundated / days since dry, site mean over
   years); (b) 180-day mean water depth, site mean over years. Whichever
   exists first is z; no alternative axes will be tried.
5. **Stacking**: per the paper's §5 protocol — site-years are columns;
   windows are W = P/4 consecutive SITES on the z axis (all their years'
   columns travel together); stride W/3; unbalanced year counts allowed.
6. **Estimator**: exactly paper §2.4 (split_mean pooled statistic, 12-point
   edge-count grid from 2S/3 to half of all pairs).
7. **Tests**: primary = free-split permutation test, site z-positions
   permuted with years attached, NPERM = 199. Secondary = TITAN-anchored
   single-look test at TITAN's decliner change point computed on the same
   matrix (site-mean abundances vs z, TITAN2 defaults as in 09_glades_titan).
   NBOOT = 199 for the CI on z-hat.
8. **Reporting**: z-hat, CI, primary p, secondary p, TITAN change points —
   published regardless of outcome. A null at this power is a bounded null
   and will be reported as such.

Ambiguities discovered in the files are resolved by the nearest rule above;
any decision not covered here will be recorded in this file BEFORE the
estimator is run, in a follow-up commit.

---

## Amendment 1 — 2026-09-19, before any test was run

Recorded after downloading the data and reading only its README and
structural summaries (row, site, year, and taxon counts; missing values).
No gradient association, window curve, or test statistic had been computed.
Each item resolves an ambiguity by the nearest rule above.

1. **Site unit (rule 1).** In the CERP program the `SITE` field has 34 values,
   each a cluster of numbered sampling points (`PLOT` = primary sampling
   unit, PSU). The "~146-site" annual panel of the source paper is the set of
   PSUs: 148 PSUs (key REGION_SITE_PLOT), 2005-2024, median 19 survey years.
   The PSU is the site unit; P = 148.
2. **File and taxa (rule 3).** `Everglades-insect-community-data.csv`, whose 26
   taxon columns (COLEOPTERA ... EPITHECA) are non-overlapping. Excluded:
   `Richness` and `Abund` (summaries) and the other file, whose aggregate
   columns (HETEROPTERA, ODONATA, INSECTS, ...) would double-count. The
   presence rule counts PSUs.
3. **Column (rule 2).** One column per PSU-year: summed abundance over that
   PSU-year's rows divided by summed `THROW`, then log(1 + x).
4. **Gradient (rule 4a).** `DSLDD` (days since the site was last dry) is the
   hydroperiod variable; z = PSU mean of DSLDD over its survey years. No
   values are missing, so rule 4b is not reached.
5. **Windows and resampling (rules 5-7).** W = round(148/4) = 37 PSUs, stride
   12; a window carries every PSU-year of its PSUs. The permutation permutes
   PSU z-values with years attached; the bootstrap resamples PSUs within a
   window with all their years. Implementation: `ecsurf.estimate_threshold`
   with `groups` = PSU, the released code.
6. **Secondary test (rule 7).** TITAN2 with the settings of
   `09_glades_titan.R` (minSplt 5, 250 permutations, 250 bootstraps, purity
   and reliability 0.95) on PSU-mean densities (not log-transformed, as for
   the Everglades benchmark) against z; the anchored split is the window
   boundary nearest TITAN's sum(z-) change point.
7. **Robustness panel (rule 1).** The 24-site MWD panel is deferred until the
   primary result is recorded.
