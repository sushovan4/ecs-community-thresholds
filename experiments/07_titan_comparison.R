# TITAN2 on the simulated gradients from 07_export_titan_inputs.py.
#
# The claim under test is COMPLEMENTARITY, in both directions: TITAN detects
# taxon-level abundance change points and is by construction blind to a pure
# reorganisation of the correlation structure; the connectance-controlled ECS
# estimator is the reverse.  Anything else would mean the two methods overlap
# and one is redundant.
#
# Run after `python3 07_export_titan_inputs.py`:
#   Rscript 07_titan_comparison.R > ../results/07_titan.txt 2>&1

suppressMessages(library(TITAN2))

root <- file.path(dirname(getwd()), "data", "titan_sim")
if (!dir.exists(root)) root <- file.path("..", "data", "titan_sim")
zstar <- 0.55

for (scen in c("THRESH", "ABUND")) {
  cat(sprintf("== %s (true z* = %.2f) ==\n", scen, zstar))
  for (i in 0:4) {
    taxa <- read.csv(file.path(root, sprintf("%s_%d_taxa.csv", scen, i)))
    env  <- read.csv(file.path(root, sprintf("%s_%d_env.csv",  scen, i)))$z
    res  <- suppressWarnings(
      titan(env, taxa, minSplt = 5, numPerm = 250, boot = TRUE, nBoot = 250,
            imax = FALSE, ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95,
            ncpus = 1, messaging = FALSE))
    filt   <- res$sppmax[, "filter"]
    n_ind  <- sum(filt > 0)
    cps    <- res$sppmax[filt > 0, "zenv.cp"]
    mods   <- sub("^t[0-9]+", "", rownames(res$sppmax)[filt > 0])
    cat(sprintf("  run %d: %2d / %2d pure & reliable indicator taxa", i,
                n_ind, ncol(taxa)))
    if (n_ind > 0)
      cat(sprintf("   indicator cps: median %.2f range %.2f-%.2f  (modules: %s)",
                  median(cps), min(cps), max(cps),
                  paste(sort(unique(mods)), collapse = " ")))
    cat("\n")
    scz <- res$sumz.cp
    cat(sprintf("         sumz-  cp %.2f [%.2f, %.2f]   sumz+  cp %.2f [%.2f, %.2f]\n",
                scz["sumz-", "cp"], scz["sumz-", "0.05"], scz["sumz-", "0.95"],
                scz["sumz+", "cp"], scz["sumz+", "0.05"], scz["sumz+", "0.95"]))
  }
  cat("\n")
}
