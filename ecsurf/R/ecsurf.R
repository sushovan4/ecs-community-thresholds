#' Locate a community-structure threshold along an environmental gradient
#'
#' Runs the pre-registered estimator (Majhi et al.): connectance-indexed
#' Euler characteristic curves over sliding windows, pooled split statistic,
#' permutation p-value, bootstrap CI.
#'
#' @param x taxa-by-samples numeric matrix (rows = taxa). A sites-by-taxa
#'   data.frame (the vegan convention) is transposed automatically.
#' @param z numeric gradient value per sample (site).
#' @param nperm,nboot permutation and bootstrap replicates (protocol: 999).
#' @param groups optional site ids per column for stacked repeat surveys:
#'   columns sharing an id carry one z and move together.
#' @param seed integer RNG seed.
#' @return list with zhat, stat, p, ci, windows, W.
#' @examples \dontrun{
#'   library(TITAN2); data(glades.taxa); data(glades.env)
#'   ecs_threshold(glades.taxa, glades.env$TP.ugL, nperm = 999, nboot = 999)
#' }
ecs_threshold <- function(x, z, nperm = 999, nboot = 999,
                          groups = NULL, seed = 0) {
  py <- reticulate::import_from_path(
    "ecsurf", system.file("python", package = "ecsurf"))
  x <- as.matrix(x)
  if (nrow(x) == length(z) && ncol(x) != length(z)) x <- t(x)
  if (!is.null(groups)) groups <- as.integer(factor(groups)) - 1L
  res <- py$estimate_threshold(x, as.numeric(z),
                               nperm = as.integer(nperm),
                               nboot = as.integer(nboot),
                               groups = groups, seed = as.integer(seed))
  res
}
