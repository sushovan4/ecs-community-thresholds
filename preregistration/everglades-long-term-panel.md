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
