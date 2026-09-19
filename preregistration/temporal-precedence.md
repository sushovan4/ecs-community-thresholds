# Pre-registration: does structural change precede compositional change?

Committed 2026-09-19, before any within-site time series was analyzed. The
commit timestamp is the registration. Companion to
`multi-dataset-program.md`, whose estimator and reporting rules apply
unchanged; only the gradient changes, from space to time.

**Question.** Applying the frozen estimator along TIME rather than along an
environmental gradient, does the co-occurrence structure of a community
reorganize *earlier* than its taxon-level composition does? A positive
answer would make the structural summary a leading indicator, which is the
claim the early-warning literature makes for correlation structure
(Scheffer et al. 2009) and which no taxon-aggregative method can test.

## Data

Long-term datasets already in hand, each split into independent series by a
grouping that is fixed here before analysis:

| Series set | Source | Grouping | Approx. years |
|---|---|---|---|
| E1-E4 | Everglades insect panel (Dryad doi:10.5061/dryad.wstqjq32d) | the four primary regions (PHD, SRS, TSL, WCA), MWD programme | 1996-2025 |
| E5-E8 | same, CERP programme | the same four regions where present | 2005-2024 |
| C1-C4 | Cedar Creek e001 (knb-lter-cdr.14) | field (A, B, C, D) | 1982-2021 |
| N1 | Niwot Saddle grid (knb-lter-nwt.93) | whole grid | 1989-2023 |

**Inclusion.** A series is analyzed if it has at least 10 distinct survey
years, at least 40 plot-years, and at least 15 taxa after the filters of the
multi-dataset registration (which apply unchanged, including the 100-taxon
cap). Series failing any of these are listed with the reason and not
analyzed.

## Analysis

1. **Units and columns.** The unit is the survey YEAR; columns are the
   plot-years within it (stacked, never averaged). z = the year.
2. **Estimator.** Exactly Section 2.4 of the paper via the released `ecsurf`
   code, with `groups` = year: windows of W = round(P/4) years (min 6),
   stride round(W/3), 12-point edge-count grid, pooled split statistic.
3. **Structural change point.** z-hat_struct = the estimated year, with a
   199-permutation test that permutes year labels with all their plots
   attached, and a 199-resample bootstrap interval. A series is
   *structurally detected* if p <= 0.05.
4. **Compositional change point.** TITAN2 (settings of
   `09_glades_titan.R`) with the environmental variable set to year, on
   year-mean densities of the same taxa. z-hat_marg = the sum(z-) community
   change point; the sum(z+) change point is also recorded.
5. **Precedence statistic.** For each series, D = z-hat_marg - z-hat_struct,
   in years. D > 0 means structure moved first.
6. **Confirmatory test (the one test of the hypothesis).** Across the
   structurally detected series, a two-sided Wilcoxon signed-rank test of
   D against zero, with the sign test reported alongside. Series that are
   not structurally detected are reported but excluded from this test, since
   their z-hat is not interpretable; their D values are reported separately.
7. **Secondary, out-of-sample form.** For each detected series, refit the
   estimator using only years up to y, for every y from the 8th survey year
   onward, and record y*, the earliest y at which the test rejects. Record
   whether z-hat_marg computed on the FULL series falls after y*. The
   proportion of series with z-hat_marg > y* is reported with a binomial
   confidence interval. This is a forecasting statement and is reported
   whatever it shows.
8. **Sensitivity of the comparison.** Because TITAN and the estimator have
   different resolutions, D is also reported in units of the structural
   window width W.

## Reporting

Every series in the table is reported: included or excluded with reason;
P (years), columns, taxa; z-hat_struct with p and interval; z-hat_marg
(both directions); D; and for detected series, y* and the out-of-sample
comparison. The confirmatory test is reported once, whatever its outcome.
No series is reanalyzed with different choices; ambiguities are resolved by
the nearest rule in the multi-dataset registration and recorded in a dated
amendment before the estimator runs.

**Known confounders, declared now.** Calendar time carries everything that
changed at a site, including sampling effort, taxonomy revisions, invasions
(the Everglades swamp-eel invasion is documented in that dataset), fire and
fencing regimes at Cedar Creek, and climate. A detected temporal change
point is not attributed to any cause here, and precedence between two
statistics computed on the same series is a statement about the statistics,
not about mechanism.
