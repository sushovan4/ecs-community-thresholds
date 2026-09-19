# Locating community thresholds with Euler characteristic surfaces

Code, simulation scripts, result logs, and paper source for

> S. Majhi, R. Pal, A. Mitra, Y. Sun. *Locating community thresholds along
> contamination gradients with Euler characteristic surfaces.* Preprint, 2026.

The method locates where a community's **co-occurrence structure** reorganizes
along an environmental gradient. Within sliding windows of plots it computes
the Euler characteristic of the flag complex of the taxon correlation graph,
**indexed by edge count** so that connectance is fixed and only arrangement
varies; a pooled two-sample split locates the change, a plot-permutation test
gives significance, and a within-window bootstrap gives an interval. It is
complementary to taxon-level threshold methods such as TITAN: each detects a
signal the other cannot.

## Layout

```
ecsurf/          R package (thin reticulate wrapper) and reference Python core
experiments/     every analysis in the paper, numbered in the order run
  ec.py                          shared Euler-characteristic machinery
  01, 02                         the connectance confound and its fix (Fig. 1)
  06                             estimator selection, power, sharpness (Table 1)
  07, 12, 10                     TITAN, PERMANOVA/PERMDISP, e-divisive baselines
  08                             six summaries through one estimator (Table 2)
  09, 15                         Everglades benchmark, primary and anchored tests
  11, 16                         repeat-survey (stacked vs averaged) power
  14                             edge-count grid sensitivity
  03, 04                         public omics stress tests
  13                             regenerates the paper figures from results/
hpc/             publication-strength driver (SLURM script included)
results/         logs of every run the paper cites
preregistration/ analysis plans committed before their data were accessed
paper/           LaTeX source (elsarticle) and figures
data/fetch.sh    downloads public inputs; nothing is redistributed
```

## Quick start

Python (needs `numpy scipy pandas gudhi networkx`):

```bash
bash data/fetch.sh                          # public inputs, incl. the Everglades benchmark
cd experiments && python3 06_gradient_estimator.py recovery
```

R (needs Python with `gudhi` visible to reticulate):

```r
remotes::install_local("ecsurf")
library(ecsurf); library(TITAN2)
data(glades.taxa); data(glades.env)
ecs_threshold(glades.taxa, glades.env$TP.ugL)   # zhat, p, ci
```

Baseline comparisons additionally need the R packages `TITAN2`, `vegan`, `ecp`.

## License

MIT (see `LICENSE`). Public datasets remain under their original terms.
