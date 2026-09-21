# Amendment 8: TITAN with a fixed seed, writing the full sppmax table so the
# INDICATOR TAXA (not just their count) are on the record.  Settings are
# identical to 25_titan_multi.R; only the seed and the outputs differ.
# Writes results/35_titan_<D>_sppmax.csv and results/35_titan_<D>.txt.
suppressMessages(library(TITAN2))
ds <- commandArgs(trailingOnly = TRUE)
for (D in ds) {
  f <- sprintf("../data/multi/titan/%s_taxa.csv", D)
  if (!file.exists(f)) { cat(D, "no export\n"); next }
  taxa <- read.csv(f, check.names = FALSE)
  env <- read.csv(sprintf("../data/multi/titan/%s_env.csv", D))$z
  ok <- is.finite(env); taxa <- taxa[ok, ]; env <- env[ok]
  set.seed(20260920)
  t0 <- Sys.time()
  res <- try(suppressWarnings(titan(env, taxa, minSplt = 5, numPerm = 250,
             boot = TRUE, nBoot = 250, imax = FALSE, ivTot = FALSE,
             pur.cut = 0.95, rel.cut = 0.95, ncpus = 1, messaging = FALSE)),
             silent = TRUE)
  if (inherits(res, "try-error")) {
    writeLines(paste("FAILED:", res[1]), sprintf("../results/35_titan_%s.txt", D))
    next
  }
  sp <- as.data.frame(res$sppmax)
  sp$taxon <- rownames(res$sppmax)
  write.csv(sp, sprintf("../results/35_titan_%s_sppmax.csv", D), row.names = FALSE)
  k <- sum(sp[, "filter"] > 0)
  writeLines(c(sprintf("%s (seed 20260920): %d units, %d taxa; %d pure & reliable indicators",
                       D, length(env), ncol(taxa), k),
               "sumz change points (gradient units):",
               capture.output(print(round(res$sumz.cp, 4))),
               sprintf("[%.0fs]", as.numeric(Sys.time() - t0, units = "secs"))),
             sprintf("../results/35_titan_%s.txt", D))
  cat(D, "done\n")
}
