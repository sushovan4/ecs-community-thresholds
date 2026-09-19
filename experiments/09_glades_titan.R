# TITAN2 on the glades data -- the incumbent's own result on its own dataset,
# regenerated here so the head-to-head with 09_glades.py cites one run of
# record rather than numbers transcribed from Baker & King (2010).
suppressMessages(library(TITAN2))
data(glades.taxa); data(glades.env)
res <- suppressWarnings(titan(glades.env, glades.taxa, minSplt = 5,
       numPerm = 250, boot = TRUE, nBoot = 250, imax = FALSE, ivTot = FALSE,
       pur.cut = 0.95, rel.cut = 0.95, ncpus = 1, messaging = FALSE))
filt <- res$sppmax[, "filter"]
cat(sprintf("pure & reliable indicator taxa: %d of %d\n",
    sum(filt > 0), ncol(glades.taxa)))
cat("sumz change points (ug/L TP):\n")
print(round(res$sumz.cp, 1))
