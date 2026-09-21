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

---

## Amendment 5 — 2026-09-19, a second-generation estimator and a second pass

The estimator registered in the paper's Section 2.4 compares two groups of
window curves, so its effective sample size is the number of windows (about
ten). An alternative built on the Intersection Euler Characteristic Profile
takes the TAXA as the sample: at a candidate split, each taxon is the point
given by its correlation profile on each side, the pooled points are
projected to four dimensions by PCA (label-blind), and the statistic is the
normalized integral of the intersection profile of the two clouds, small
values meaning the two structures have separated. The null permutes units
along the gradient and rebuilds everything, as in every test here.

In simulation (`26_mixup_split.txt`, `27_mixup_stress.txt`, fixed before
this amendment and before any re-analysis) it dominates the registered
statistic: power 0.88 against 0.57 at 160 units, 0.60 against 0.25 at 80;
0.90 against 0.37 at 80% occupancy and 0.47 against 0.13 at 40%; and
false positives 0.07 against 0.15 under regional confounding. Both fail
below about 15% occupancy.

**Registered now, before it is applied to any dataset.** The second-generation
estimator is applied to every dataset and series already analyzed here:
D1-D8, the Everglades long-term panel, and the thirteen temporal series
(annual and occasion-level). Settings: candidate splits at unit quantiles
0.3, 0.4, 0.5, 0.6, 0.7; d0 = 4; 199 unit permutations; within-region
permutation additionally where a region field is registered. Reported for
each: the located split on the gradient, p, and the matrix sparsity. Every
dataset is reported whether or not it detects, and the first-generation
results stand unchanged beside them. This is a second pass with a better
instrument, not a replacement of the record.

Expectation, recorded in advance: the analyzed matrices have 70%-96% zeros
(4%-30% occupancy), which the simulations place at or below the regime where
the new estimator reaches power 0.47, so most are expected to remain null.

---

## Amendment 6 — 2026-09-19, the two-stage ordering test

Registered before the second-pass results for the independent spatial
datasets exist (D1 and D5 were running when this was written; D2, D3, D4 had
not started, and no ordering statistic has been computed for any dataset).

**Motivation.** Amendment 2A registered a cross-dataset ordering test that
required at least four detecting datasets and was reported as underpowered
when only one detected. The second-generation estimator detects far more
often, so the test may become runnable. But detection and localization are
different jobs, and the interaction statistic's located split moved between
split grids on the Everglades panel, so the ordering must not rest on it
until its localization is characterized.

**The procedure, fixed now.**
1. **Detection** is decided by the intersection-ECP test of amendment 5
   (p <= 0.05 on the registered 199-permutation test; where a region field
   exists, the within-region test must also reject).
2. **Localization** uses whichever estimator the simulation in
   `29_mixup_localization.txt` shows to be the more accurate recoverer of a
   known threshold, chosen on that simulation alone and applied uniformly to
   every dataset. If the two are within 0.02 of gradient span in RMSE, the
   registered curve statistic is used.
3. **The ordering statistic** is R = (z_TITAN - z_struct) / IQR(z over units),
   with z_TITAN the sum(z-) community change point on the same unit-mean
   matrix. R > 0 means structure changes at a lower gradient value.
4. **The confirmatory test** is a two-sided Wilcoxon signed-rank test of R
   against zero across the INDEPENDENT spatial datasets that detect (D1-D8
   and the Everglades long-term panel; the Everglades temporal series are
   not independent of it and are excluded). It runs only if at least four
   such datasets detect; otherwise the R values are reported descriptively
   and no test is claimed.
5. Every dataset's R is reported whether it detects or not.

No dataset will be re-analyzed with different settings, and this test is run
exactly once.

---

## Deviation recorded — 2026-09-20: the ordering test was run early

Amendment 6 states that the ordering test runs exactly once. It was in fact
run once with incomplete inputs, and the output is kept as
`results/31_ordering_incomplete.txt`. At that point the second-pass
detection files for the Everglades long-term panel and for D6-D8 were
missing (they had been deleted when the split grid was widened and not yet
regenerated), and TITAN's change point for D1 was still computing, so those
five rows were treated as non-detecting or incomplete. The four datasets
with complete inputs gave a median R of 0.039, two of four positive, and a
Wilcoxon p of 1.00.

The complete run follows once the missing inputs exist, and both are
reported. No specification was changed between the two runs: the estimator,
the detection rule, the localization choice, the statistic and the test are
exactly those registered in amendment 6. The early look is disclosed rather
than discarded because discarding it would leave the record incomplete.

---

## Deviation recorded — 2026-09-20: the split grid was widened

Amendment 5 registered candidate splits at unit quantiles 0.3, 0.4, 0.5, 0.6,
0.7. The runs reported use fifteen quantiles from 0.15 to 0.85 in steps of
0.05. The reason: on the Everglades long-term panel the located split sat at
the edge of the narrow grid, so the grid, not the data, was choosing the
answer, and a located split pinned at a boundary is not a location estimate.
The grid was widened before the second-pass detection results for D2, D3, D4,
D6, D7 and D8 existed, and the D1, D5 and Everglades runs made on the narrow
grid were rerun; those narrow-grid outputs were not retained, which is a gap
in the record.

The change does not affect the validity of the test. Every permutation
replicate minimizes over the same grid as the observed statistic, so the
null distribution is the distribution of the same function of the data and
the p-value stays exact whatever the grid contains. What the grid changes is
the located split and, through the minimum, the power.

Because the registered setting was departed from, the whole second pass is
also run on the registered 0.30-0.70 grid and both are reported for every
dataset (`results/28_mixup_<name>_reg.txt` against
`results/28_mixup_<name>.txt`). The widened grid is the one used for the
ordering test of amendment 6, as it was when that amendment was written.

---

## Amendment 7 — 2026-09-20, the interaction test on the TITAN benchmark

Written and committed before the run. Amendment 5 listed the datasets the
second pass would cover: D1-D8, the Everglades long-term panel and the
thirteen temporal series. It did not include the `glades` dataset shipped
with TITAN2 (126 marsh sites, 164 macroinvertebrate taxa, total-phosphorus
gradient), which the paper analyzes separately as the incumbent's own
benchmark. The registered curve statistic did not reject there (p = 0.30),
and that non-detection is reported in the paper.

**Registered now.** The interaction test of amendment 5 is applied to that
matrix once, with the settings used for every other dataset: the same
preprocessing already fixed for it (taxa in at least five sites,
log(1+abundance), sites ordered by TP), the taxon cap of 100 that every
other dataset carries, both split grids (the widened 0.15-0.85 and the
registered 0.30-0.70), 199 unit permutations. There is no region field, so
no restricted permutation. Output: `results/33_glades_mixup.txt`.

**Expectation, recorded in advance.** 126 sites of sparse counts sits below
the P = 160 row of the simulation, where the interaction test reaches power
0.47 at 40% occupancy. A detection is roughly a coin flip and a null would
be uninformative about the biology.

**Commitments.** The result is reported whichever way it falls, beside the
curve statistic's p = 0.30 on the same matrix. It is an extra look at a
dataset already analyzed, taken after the interaction test was seen to
detect widely elsewhere, and it is labeled as such in the paper. It does
not enter the amendment-6 ordering test, whose dataset list was fixed
before this amendment.

---

## Amendment 8 — 2026-09-20, do the two channels run on different taxa?

Registered before the statistic has been computed for any dataset other than
the Everglades long-term panel, and before the second-pass and TITAN inputs
for the remaining datasets exist.

**Disclosure, first.** This test generalizes a result already seen. On the
Everglades panel the five taxa carrying 54% of the structural change were
none of TITAN's eight pure and reliable indicator taxa, and those indicators
contributed at most 1% each. That observation motivates this amendment, so
this is a registered generalization of a result seen once, not a blind
prediction, and the paper will say so in those words. The Everglades panel
is reported alongside but excluded from the confirmatory test.

**The question.** When both instruments fire on the same matrix, do they
implicate the same species? The ordering test of amendment 6 asks whether
the two channels change at the same PLACE. This asks whether they are
carried by the same TAXA.

**The statistic.** For each dataset where the interaction test detects
(amendment 5 rule, including the within-region test where a region field
exists), take the curvature attribution of Section 2.5 of the paper at the
split located by the curve statistic, giving each taxon i a contribution
|dkappa_i| to the change in chi. Let I be TITAN's pure and reliable
indicator taxa on the same unit-mean matrix. Define

    C = ( sum_{i in I} |dkappa_i| / sum_i |dkappa_i| ) / ( |I| / S ).

C = 1 means the indicator taxa carry exactly their proportional share of the
structural change; C < 1 means the structural change is carried
disproportionately by taxa TITAN does not flag. Reported per dataset with a
permutation p-value from drawing subsets of size |I| uniformly at random
from the S taxa (9999 draws), which is exact.

**The confirmatory test** is a two-sided Wilcoxon signed-rank test of
log C against zero across the independent detecting datasets D1-D8,
excluding the Everglades panel for the reason given above. It runs only if
at least four such datasets detect and have both a curvature attribution and
a TITAN indicator set; otherwise the C values are reported descriptively and
no test is claimed. A secondary, coarser reading is also reported: the
overlap between the five largest contributors and I, against the
hypergeometric null.

**TITAN's stochasticity.** The existing `25_titan_<D>.txt` runs recorded only
the NUMBER of indicators, not their names, and were run without a fixed
seed. A seeded rerun with identical settings (minSplt 5, 250 permutations,
250 bootstraps, purity and reliability 0.95) writes the full sppmax table
and is the input to this test. The change points from that seeded rerun are
also reported beside the originals as a reproducibility check, and the
amendment-6 ordering test continues to use the originally registered
outputs.

**Commitments.** Every detecting dataset is reported whether its C is below
1 or not, the test runs once, and a null result is reported as a null.
Failure here means the Everglades disjointness was idiosyncratic, which is
worth knowing and will be stated as such.

---

## Amendment 9 — 2026-09-20, a tenth dataset with a published incumbent answer

Written and committed **before the data were downloaded**; every choice below
is fixed from the Dryad landing page and the published abstract alone.

**The dataset.** Payne et al., *Plant thresholds and community composition of
coastal marsh-forest ecotones in the US Northeast*, Ecosphere (2026);
data at doi:10.5061/dryad.5tb2rbpcm, published 2025-12-15, file
`Ecotone_PlantSpp-EnvData.csv`. Three sites experiencing marsh upslope
migration (Waquoit Bay MA, Pine Neck NY, Egg Harbor NJ), understory percent
cover by species, with soil salinity, moisture, bulk density, organic
matter, redox, light, water depth, elevation and flooding duration per
sampling point.

**Why this one, stated plainly.** It is the only dataset in this programme
whose incumbent answer is already in print: the authors report TITAN
community change points at **0.8 and 7.6 PSU** soil salinity. Every other
TITAN comparison here was regenerated by us. This one was produced
independently, by the method's own users, before we saw the data.

**Registered choices.**
1. *Gradient* z = soil salinity (PSU, from 5:1 water:soil extracts). Chosen
   because it is the axis the published change points are on. No other
   environmental variable is used as a gradient.
2. *Unit* = sampling point. *Taxa* = understory cover species; saplings and
   mature trees are excluded, as those are separate response variables in
   the source paper, not the understory community.
3. The standing inclusion rule applies unchanged: at least 40 units and at
   least 15 taxa after the 5-unit occurrence rule, else the dataset is
   excluded and reported as excluded.
4. *Region* = site (three sites), so the within-region permutation is the
   test of record.
5. *Primary analysis* is the two-stage procedure exactly as applied to every
   other dataset: interaction test for detection (199 unit permutations,
   the widened 0.15-0.85 split grid, also reported on the registered
   0.30-0.70 grid), curve statistic for localization, taxon cap 100,
   log(1+cover).

**Compositionality, registered as a secondary analysis.** Percent cover is
closer to closed data than any other substrate here, and correlations among
closed compositions are the classic spurious-network trap. The primary run
therefore has a pre-registered companion: the identical pipeline on
centred-log-ratio-transformed cover (zeros handled by a multiplicative
replacement at half the smallest positive value), reported beside the
primary whether or not the two agree. If they disagree, that disagreement
is the result and will be reported as such.

**What counts as what.** A detection whose located split falls inside
[0.8, 7.6] PSU would be concordance with an independently published
threshold. A detection outside it, or a null, is reported as plainly.
This dataset does not enter the amendment-6 ordering test or the
amendment-8 attribution test, whose dataset lists were fixed before it.

**Authorship note.** No author of the source dataset is involved in this
analysis at the time of writing. If that changes, it will be recorded here
with a date, and the analysis above will not be altered.

### Amendment 9, implementation note — 2026-09-20, fixed from the README

Recorded after downloading the file and reading its README, before any
statistic was computed. Amendment 9 said "understory cover species"; the
file needs that made exact.

The README states that fields 6-79 are percent cover in 0.5 m^2 quadrats of
"understory plant species, plants, wrack, dead material, or bare ground".
The taxa are therefore columns `Poaceae_sp` through `PHAU` inclusive (70
columns), and the following are excluded because they are not species:
`Wrack`, `Dead`, `Bare`, the `Total` sum, the derived cover groupings
`ST_Gram`, `ST_Shrub`, `ST_Herb`, `Phrag`, `SI Phrag_Stems`, and every
`_s` (sapling) and `_t` (tree) column, the last per amendment 9 itself.
Unknown morphospecies (`UnForb1`, `UnShrub2`, ...) are kept: they are
distinct taxa in the source data and dropping them would change the
community.

The file has 125 rows and three sites (PN 43, WB 41, EH 41). Salinity runs
0.09-42.58 PSU. `Site` is the region for the restricted permutation. The
unit is the row (a sampling point); there is one quadrat per point, so
units and columns coincide, as they do for the single-survey datasets
D1-D8.

---

## Correction recorded — 2026-09-20: every interaction-test number is being recomputed

Not a change of specification. A defect in the implementation, found while
checking why one dataset gave two different answers on two machines.

The intersection profile is computed from Alpha complexes, and the vendored
code requested GUDHI's `precision="fast"`, which uses inexact geometric
predicates. The taxon clouds here are near-degenerate -- taxa with identical
sparse occurrence patterns have identical correlation profiles -- so the
underlying Delaunay construction sits on ties that it resolves differently
depending on the linear-algebra backend. On the Cedar Creek matrix at the
0.15 split, `fast` returned T = 12.3814 where both `safe` and `exact` return
8.1186, and the located split moved between machines.

`safe` (interval arithmetic with an exact fallback) agrees with `exact` to
ten decimal places and costs about 0.7 of `exact`, roughly nine times
`fast`. It is now the default, overridable by `ECP_PRECISION`.

**Every number produced by the interaction test is therefore being
recomputed**: all of `28_mixup_*`, the TITAN benchmark run of amendment 7,
the ecotone run of amendment 9, and the simulation studies
(`26`, `27`, `29`, `30`, `34`, `37`). The superseded outputs are preserved
under `results/fast_precision_superseded/` rather than deleted, and the
paper will report any verdict that changes.

No pre-registered specification was altered: the estimator, the split
grids, the permutation schemes, the seeds and the dataset list are exactly
as registered. Only the geometric predicate changed, from one that is not
reproducible to one that is.
