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

---

## Amendment 2 — 2026-09-19, cross-dataset tests, before any was computed

Three tests across datasets are registered here. At this point the spatial
programme has produced exactly one result (D5, no detection); no TITAN change
point, attribution, or ordering statistic has been computed for any dataset.

**A. Ordering (does structure change earlier on the gradient than
composition?).** For every dataset where the primary test detects
(p <= 0.05), let z_struct be z-hat and z_marg be TITAN's sum(z-) community
change point on the same unit-mean matrix (rule 8 settings). The ordering
statistic is R = (z_marg - z_struct) / (gradient interquartile range of the
units), a scale-free signed quantity; R > 0 means structure changes at a
lower gradient value. The confirmatory test across detected datasets is a
two-sided Wilcoxon signed-rank test of R against zero, reported once with
the sign test beside it. Datasets that do not detect are reported with their
R values but excluded from the test, since z-hat is not interpretable there.
The Everglades long-term panel (a separate registration) is included as one
dataset. Should fewer than four datasets detect, the test is reported as
underpowered and the individual R values stand as description.

**B. Attribution (do the taxa carrying structural change differ from the
taxa whose abundance marks the threshold?).** For every detected dataset,
per-taxon curvature changes are computed at the detected split by the exact
Gauss-Bonnet decomposition (as in `21_curvature_attribution.py`). The
comparison statistic is the overlap between the five taxa with the largest
absolute curvature share and TITAN's pure-and-reliable indicator taxa,
against the overlap expected if the two sets were drawn independently
(hypergeometric). Reported per dataset and pooled by summing observed and
expected overlaps.

**C. Sub-annual temporal series.** The temporal registration's unit is the
survey year; ten series so detected nothing, consistent with 20-35 units
being below the method's operating range. Where a programme surveys several
times per year, the survey OCCASION (year x period) is a valid unit and
multiplies the unit count. Registered now, before any such series is built:
the Everglades MWD programme's five annual periods give occasion-level
series for the same four regions (E1-E4 become T1-T4, with PERIOD4 ordinal
within a water year), z = the occasion's decimal date, columns = plot
samples within the occasion, all other rules unchanged. Inclusion still
requires 40+ occasions, 40+ columns, 15+ taxa. The year-level results
already obtained are reported unchanged beside these.

---

## Amendment 3 — 2026-09-19, dataset D8, before download

Added after the operating range of the paper's Section 3.5 was measured, and
chosen to match it: quantitative percent cover rather than sparse counts,
thousands of units, and a gradient along which a sharp transition is
ecologically expected. Registered before any CRMS file was downloaded.

| # | Dataset | Units | Gradient z | Region field |
|---|---|---|---|---|
| D8 | Coastwide Reference Monitoring System (CRMS), coastal Louisiana marsh vegetation (`CRMS_Marsh_Vegetation.zip` from cims.coastal.la.gov) | vegetation station (site x station), station-years stacked | soil porewater salinity or specific conductance at the shallowest reported depth, station mean over years (`CRMS_Soil_Properties`) | CRMS hydrologic basin |

All rules of this registration and amendment 1 apply unchanged: percent cover
per species per station-year is the column value, log(1+x); taxa present in
at least 5 stations, capped at the 100 most widespread; windows of P/4 units;
199 permutations and bootstrap resamples; the within-region permutation is a
registered secondary test because marsh salinity is strongly spatially
structured; TITAN on station-mean cover for comparison.

**Fallback fixed in advance.** If porewater salinity is not reported per
vegetation station (only per site, or on a different station set that cannot
be joined by the identifiers present), the gradient becomes station elevation
(ft NAVD88) from the USGS release
`CRMS_2014_Vegetation_Station_Elevation_Data.csv`, and that substitution is
recorded in a dated amendment before the estimator runs. No other gradient
variable will be tried, and the dataset is analyzed once.

**Declared risk.** Louisiana marshes are often near-monodominant, so despite
being quantitative, the cover matrix may be zero-heavy; per the operating
range that would lower power. The analysis is run coastwide, not restricted
to a basin, and the outcome is reported either way.

---

## Amendment 4 — 2026-09-19, D8 gradient, before the estimator ran

The fallback declared in amendment 3 is triggered, and this records it from
structural checks alone (identifier overlap and column names; no gradient
association or test statistic computed).

Soil porewater salinity in `CRMS_Soil_Properties.csv` is measured at SOIL
stations (identifiers ending `-S##`) and vegetation at VEGETATION stations
(`-V##`): 3,171 soil stations and 4,619 vegetation stations with **zero**
identifier overlap, though 392 sites are shared. Salinity is therefore not a
per-unit measurement for the vegetation stations, and rule 2 is not
satisfiable with it.

Per the pre-declared fallback, z is **station elevation, feet NAVD88
(Geoid 12A)**, from `CRMS_2014_Vegetation_Station_Elevation_Data.csv` (USGS
ScienceBase item 606dde87d34eae125e9c75a8), which joins to 3,180 vegetation
stations. Resolutions of detail, fixed here: all stations carrying an
elevation are used (the file's `Exclude_from_Site_Mean` flag governs
site-mean aggregation, not station validity, and we do not aggregate to
sites); the unit is the vegetation station; columns are station-years, with
cover summed per species within a station-year; the species field is
`Scientific Name As Currently Recognized`; the region field for the
within-region test is `Basin` (9 basins).

Marsh elevation is the proximate control on flooding duration and hence on
zonation, so it is the ecologically appropriate axis; it is also, unlike
salinity, measured at the vegetation stations themselves.
