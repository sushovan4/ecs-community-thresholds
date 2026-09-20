# TITAN2 pure & reliable indicator taxa on the Everglades long-term panel
# (same inputs and settings as 17_everglades_lt_titan.R), for comparison
# with the curvature attribution in 21_curvature_attribution.py.
suppressMessages(library(TITAN2))
taxa <- read.csv("../data/ever_lt/titan_taxa.csv")
env  <- read.csv("../data/ever_lt/titan_env.csv")$z
set.seed(1)
res <- suppressWarnings(titan(env, taxa, minSplt = 5, numPerm = 250,
       boot = TRUE, nBoot = 250, imax = FALSE, ivTot = FALSE,
       pur.cut = 0.95, rel.cut = 0.95, ncpus = 1, messaging = FALSE))
sp <- res$sppmax
keep <- sp[, "filter"] > 0
out <- data.frame(taxon = rownames(sp)[keep],
                  group = ifelse(sp[keep, "filter"] == 1, "z-", "z+"),
                  cp = sp[keep, "zenv.cp"])
write.csv(out, "../results/21_titan_indicators.txt", row.names = FALSE)
print(out)
