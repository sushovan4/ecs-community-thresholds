# e-divisive on the raw ordered plot sequences from 10_edivisive.py.
# Reports, per cell: how often a change point is found (sig.lvl 0.05) and
# where (gradient position of the estimated split), against true z* = 0.55.
suppressMessages(library(ecp))
root <- "../data/ediv_sim"
cells <- list(c("THRESH","80"), c("STRENGTH","80"), c("REWIRE","80"),
              c("TRI_NULL","80"), c("REWIRE","160"), c("TRI_NULL","160"))
for (cell in cells) {
  scen <- cell[1]; P <- cell[2]
  found <- 0; locs <- c()
  for (i in 0:9) {
    X <- as.matrix(read.csv(file.path(root, sprintf("%s_%s_%d.csv", scen, P, i))))
    z <- read.csv(file.path(root, sprintf("%s_%s_%d_z.csv", scen, P, i)))[[1]]
    res <- e.divisive(X, sig.lvl = 0.05, R = 199, min.size = 10)
    cps <- res$estimates[c(-1, -length(res$estimates))]  # interior points
    if (length(cps) > 0) { found <- found + 1; locs <- c(locs, z[cps[1]]) }
  }
  cat(sprintf("%9s P=%3s: detected %2d/10", scen, P, found))
  if (length(locs) > 0)
    cat(sprintf("   first-cp gradient position: median %.2f range %.2f-%.2f",
        median(locs), min(locs), max(locs)))
  cat("   (true z* = 0.55)\n")
}
