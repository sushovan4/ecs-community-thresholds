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

## Candidate list — fixed 2026-09-19, before any of these files was downloaded

Gradients were chosen from landing pages and metadata only. "Region" names
the field used for the rule-7 within-region test.

| # | Dataset | Units | Gradient z | Region field |
|---|---|---|---|---|
| D1 | EPA NRSA 2018-19 benthic macroinvertebrate counts + water chemistry (national) | site (UNIQUE_ID); within-cycle revisits stacked | specific conductance (uS/cm) | EPA 9 aggregate ecoregions |
| D2 | D1 restricted to the Appalachian aggregate ecoregions (northern + southern Appalachians), the coal-mining conductivity question | as D1 | specific conductance | none (single region pair; rule 7 via state if present) |
| D3 | EPA EMAP Mid-Atlantic Highlands streams 1993-96, benthic counts + chemistry (acid mine drainage region) | site; revisits stacked | specific conductance | ecoregion or state field if present |
| D4 | EPA NCCA 2015 estuarine benthic counts + sediment chemistry | site | sediment copper (dry weight); zinc if copper is not reported | NCCA region/province |
| D5 | EPA NLA 2017 lake benthic counts + water chemistry | lake site; revisits stacked | total phosphorus | EPA 9 aggregate ecoregions |
| D6 | Cedar Creek e001 aboveground biomass by species (knb-lter-cdr.14) | plot; years stacked | N addition rate (g N m-2 yr-1) | field (A-D) |
| D7 | Niwot Ridge Saddle grid plant composition (knb-lter-nwt.93) + Saddle snow depth (knb-lter-nwt.31) | 1 m2 plot; years stacked | mean snow depth at the plot's grid point | none |

Considered and not included, with reason: England EA BIOSYS (chemistry in a
separate archive, joining rules not fixable from metadata); Konza PVC02
(treatment replicated at watershed level only); Park Grass (registration
required, factorial non-monotone treatments); Chesapeake benthos (fixed
stations likely < 40); SCCWRP Bight, Swedish MVM, CABIN, ECCC oil sands
(access or chemistry linkage unconfirmed); WV DEP, MBSS, Clark Fork,
Arkansas River, Sudbury, Harjavalta (not publicly archived as community
matrices); Blackbird Mine, NF Clear Creek, Pomeranz NZ (too few sites).

Order of execution follows the table. D6 and D7 depend on a public DataONE
mirror of EDI, which since July 2026 requires login for direct access; if the
mirror refuses, they are reported as excluded (inaccessible) under rule 1.

---

## Amendment 1 — 2026-09-19, after download, before any estimator was run

Recorded from file headers and structural counts only (units, visits,
years, taxon counts, field names, analyte names). No gradient association,
curve, or test statistic had been computed for any dataset.

**General (all datasets).**
- **Taxon cap for computability.** The flag-complex computation counts
  cliques in a graph with up to half of all taxon pairs as edges; above
  roughly 150 taxa this is infeasible at these unit counts. Where more than
  100 taxa pass rule 4, the 100 taxa present in the most units are kept
  (ties broken alphabetically). The choice uses occurrence only, never the
  gradient. This is itself recorded as a limit of the method.
- Units lacking a gradient measurement are dropped; columns (visits) lacking
  community data are dropped.

**D1 NRSA 2018-19.** Unit UNIQUE_ID; one column per visit (UID); abundance
TOTAL300 (fixed 300-count subsample) for rows with IS_DISTINCT300 = 1; taxon
TARGET_TAXON. z = COND_RESULT per visit, averaged per unit. Both sample types
(BERW, BETB) kept; each site has one. Region: AG_ECO9.
**D2** as D1 with AG_ECO9 in {NAP, SAP}; region for rule 7: STATE.
**D3 EMAP MAHA.** Unit STRM_ID; one column per (STRM_ID, YEAR, VISIT_NO),
abundance ABUND summed over that visit's samples (pool and riffle), rows with
DISTINCT = Y; taxon TAXANAME. z = COND per visit, averaged per unit. Region
for rule 7: the two-letter state prefix of STRM_ID (no ecoregion field).
**D4 NCCA 2015 estuarine.** Unit SITE_ID; one column per UID; abundance TOTAL
for IS_DISTINCT = 1; taxon TARGET_TAXON. z = CU (ug/dry g) per UID, averaged
per unit. Region: NCCA_REG.
**D5 NLA 2017.** Unit SITE_ID; one column per UID; abundance TOTAL for
IS_DISTINCT = 1; taxon TARGET_TAXON. z = PTL (total phosphorus) per UID,
averaged per unit. Region: AG_ECO9 from the site file.
**D6 Cedar Creek e001.** Unit (Field, Plot); one column per (unit, Year);
abundance Biomass (g/m2); taxon Species, excluding non-taxon categories
(names containing litter, miscellaneous, moss, lichen, fungi, unsorted, or
total, case-insensitive). z = NAdd. Region: Field. Fence removal (2004) and
burning regimes are not modeled and are reported as confounders.
**D7 Niwot Saddle.** Unit plot (1-88); one column per (plot, year);
abundance = point-intercept hit count per USDA_code; non-plant codes
(beginning with "2", e.g. litter, rock, bare ground) excluded. z = mean of
mean_depth at the matching snow stake (point_ID = plot) over all dates.
No region field.
