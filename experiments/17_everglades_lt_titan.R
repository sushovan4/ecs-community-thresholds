# TITAN2 on the Everglades long-term panel (pre-registration rule 7, amendment 6):
# PSU-mean densities vs z = PSU mean days since last dry, settings as in
# 09_glades_titan.R. Writes the sum(z-) change point for the anchored test.
# Run from experiments/ after `python3 17_everglades_lt.py primary`.
suppressMessages(library(TITAN2))
taxa <- read.csv("../data/ever_lt/titan_taxa.csv")
env  <- read.csv("../data/ever_lt/titan_env.csv")$z
res <- suppressWarnings(titan(env, taxa, minSplt = 5, numPerm = 250,
       boot = TRUE, nBoot = 250, imax = FALSE, ivTot = FALSE,
       pur.cut = 0.95, rel.cut = 0.95, ncpus = 1, messaging = FALSE))
filt <- res$sppmax[, "filter"]
cat(sprintf("pure & reliable indicator taxa: %d of %d\n", sum(filt > 0), ncol(taxa)))
cat("sumz change points (days since last dry):\n")
print(round(res$sumz.cp, 1))
writeLines(format(res$sumz.cp["sumz-", "cp"]), "../results/17_everglades_lt_titan_cp.txt")
