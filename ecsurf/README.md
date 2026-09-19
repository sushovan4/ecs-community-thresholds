# ecsurf

R interface to the connectance-indexed Euler characteristic surface
threshold estimator (Majhi et al.). Thin `reticulate` wrapper over the
reference Python implementation in `inst/python/ecsurf.py`.

```r
# needs Python with numpy, scipy, gudhi on the reticulate path
install.packages("remotes")
remotes::install_local("ecsurf")

library(ecsurf); library(TITAN2)
data(glades.taxa); data(glades.env)
res <- ecs_threshold(glades.taxa, glades.env$TP.ugL)
res$zhat; res$p; res$ci
```

For stacked repeat surveys (the powered design of the paper's sec. 3.4),
pass one column per site-year and `groups = site_id`.

Note: point reticulate at a Python that has gudhi installed, e.g.
`Sys.setenv(RETICULATE_PYTHON = "/path/to/python3")` before loading.
