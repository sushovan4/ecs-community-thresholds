# TITAN2 on each analyzed spatial dataset (registration rule 8): the
# taxon-level community change points reported beside z-hat.  Inputs are the
# unit-mean matrices exported by 19_multi.py.  Smallest datasets first.
suppressMessages(library(TITAN2))
ds <- commandArgs(trailingOnly = TRUE)
for (D in ds) {
  f <- sprintf("../data/multi/titan/%s_taxa.csv", D)
  if (!file.exists(f)) { cat(D, "no export\n"); next }
  taxa <- read.csv(f); env <- read.csv(sprintf("../data/multi/titan/%s_env.csv", D))$z
  ok <- is.finite(env); taxa <- taxa[ok, ]; env <- env[ok]
  t0 <- Sys.time()
  res <- try(suppressWarnings(titan(env, taxa, minSplt = 5, numPerm = 250,
             boot = TRUE, nBoot = 250, imax = FALSE, ivTot = FALSE,
             pur.cut = 0.95, rel.cut = 0.95, ncpus = 1, messaging = FALSE)),
             silent = TRUE)
  out <- sprintf("../results/25_titan_%s.txt", D)
  if (inherits(res, "try-error")) { writeLines(paste("FAILED:", res[1]), out); next }
  sp <- res$sppmax; k <- sum(sp[, "filter"] > 0)
  cp <- round(res$sumz.cp, 4)
  writeLines(c(sprintf("%s: %d units, %d taxa; %d pure & reliable indicators",
                       D, length(env), ncol(taxa), k),
               "sumz change points (gradient units):",
               capture.output(print(cp)),
               sprintf("[%.0fs]", as.numeric(Sys.time() - t0, units = "secs"))), out)
  cat(readLines(out), sep = "\n")
}
