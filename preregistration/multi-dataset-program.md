# Pre-registration: the frozen estimator across public gradient datasets

Committed 2026-09-19, before any dataset named here was downloaded. The
commit timestamp is the registration. Companion to
`everglades-long-term-panel.md`, whose rules this generalizes.

**Purpose.** Establish how far the method reaches on real data: apply the
estimator of the paper's Section 2.4, unchanged, to every public dataset in
the candidate list below that meets the structural requirements, and report
every one of them.

## Rules (in order of application)

1. **Inclusion.** A dataset is analyzed if it is publicly downloadable and,
   after the taxon filter, has at least 15 taxa, at least 40 sampling units,
   and a per-unit measured ordered gradient. A dataset that fails any of
   these is listed in the report with the reason. No other exclusions.
2. **Gradient.** The gradient variable z for each dataset is fixed in the
   candidate list below from its metadata alone, before download, preferring
   the variable the data providers name as the primary stressor or design
   factor. z is the unit's mean over surveys. No alternative gradient
   variables are tested.
3. **Sampling unit and columns.** The unit is the finest spatial unit that
   carries its own gradient measurement (site, plot, PSU). One column per
   unit-survey (for single-survey data, per unit): abundance per unit effort
   where effort is recorded (else raw abundance or cover), then
   log(1 + x). Repeat surveys of a unit are stacked, never averaged.
4. **Taxa.** The finest non-overlapping taxonomic resolution provided;
   summary or aggregate columns (totals, richness, higher-rank sums of listed
   taxa) are excluded; taxa must be present in at least 5 units.
5. **Estimator.** Exactly Section 2.4 via the released `ecsurf` code: windows
   of W = round(P/4) units (min 6), stride round(W/3), 12-point edge-count
   grid from 2S/3 to half of all pairs, pooled split statistic.
6. **Tests.** Primary: free-split permutation test, unit z-values permuted
   with all their surveys attached, 199 permutations; 95% bootstrap interval
   from 199 unit resamples within windows. If the primary p lies in
   [0.01, 0.10], a 1999-permutation precision check of the same test is run
   and reported alongside; the 199-permutation p remains the registered
   verdict.
7. **Spatial structure.** Where the data carry a region, basin, ecoregion,
   stratum, or site-cluster field, a within-region permutation test is a
   registered secondary test, because plot-level permutation is
   anti-conservative when the gradient is spatially structured.
8. **Comparison.** TITAN2 (settings of `09_glades_titan.R`) on unit-mean
   densities versus z; its change points are reported beside z-hat.
9. **Reporting.** For every candidate: included or excluded (with reason);
   P, S, surveys, z range; z-hat, CI, primary p, precision-check p where
   triggered, within-region p where applicable, TITAN change points. Results
   are reported regardless of outcome. A dataset's analysis is not repeated
   with different choices.

Ambiguities met in a specific dataset are resolved by the nearest rule and
recorded in a dated per-dataset amendment committed before its estimator
runs.

## Candidate list

(appended before any download; see amendments below)
