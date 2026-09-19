# PERMANOVA and PERMDISP on the simulated gradients -- the ecology incumbents.
#
# PERMANOVA (vegan::adonis2, Anderson 2001) is the field's default test for
# "do two groups of sites differ in community composition"; PERMDISP
# (betadisper) tests dispersion differences.  Both are given the ORACLE split
# at the true threshold z* = 0.55 -- their best possible shot, since neither
# localises.  Uses the same exported sequences as the e-divisive baseline.
suppressMessages(library(vegan))
root <- "../data/ediv_sim"
cells <- list(c("THRESH","80"), c("STRENGTH","80"), c("REWIRE","80"),
              c("TRI_NULL","80"), c("REWIRE","160"), c("TRI_NULL","160"))
cat("oracle-split two-group tests at z*=0.55, Euclidean, 199 permutations\n\n")
cat(sprintf("%9s %5s  %s\n", "scenario", "P", "rejections/10: PERMANOVA  PERMDISP"))
for (cell in cells) {
  scen <- cell[1]; P <- cell[2]
  ra <- 0; rd <- 0
  for (i in 0:9) {
    X <- as.matrix(read.csv(file.path(root, sprintf("%s_%s_%d.csv", scen, P, i))))
    z <- read.csv(file.path(root, sprintf("%s_%s_%d_z.csv", scen, P, i)))[[1]]
    grp <- factor(z > 0.55)
    D <- dist(X)
    pa <- adonis2(D ~ grp, permutations = 199)$`Pr(>F)`[1]
    pd <- permutest(betadisper(D, grp), permutations = 199)$tab$`Pr(>F)`[1]
    ra <- ra + (pa <= 0.05); rd <- rd + (pd <= 0.05)
  }
  cat(sprintf("%9s %5s  %27d %9d\n", scen, P, ra, rd))
}
