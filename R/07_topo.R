# omnikit — 0-dimensional persistent homology, dependency-free.
#
# Vietoris-Rips H0 of a point cloud is the single-linkage MST: each finite bar is
# one merge (birth 0, death = merge height). `stats::hclust(method = "single")`
# gives those heights with no extra dependency.

#' H0 persistence diagram (deaths = single-linkage merge heights), descending
#'
#' @param points numeric matrix, one row per point
#' @return numeric vector of finite bar lengths, descending
#' @export
omni_h0_diagram <- function(points) {
  h <- stats::hclust(stats::dist(points), method = "single")$height
  sort(h, decreasing = TRUE)
}

#' Shannon entropy of a persistence diagram's normalised bar lengths
#'
#' @param deaths numeric vector of bar lengths
#' @export
omni_persistent_entropy <- function(deaths) {
  d <- deaths[deaths > 0]
  if (length(d) == 0) return(0)
  p <- d / sum(d)
  -sum(p * log(p))
}

#' Bottleneck distance between two H0 diagrams (common birth 0)
#'
#' Bars are matched in descending order; each unmatched bar matches the diagonal
#' at cost death / 2. Exact for 1-D (H0) diagrams.
#'
#' @param a,b numeric vectors of H0 bar lengths
#' @export
omni_bottleneck_h0 <- function(a, b) {
  a <- sort(a, decreasing = TRUE); b <- sort(b, decreasing = TRUE)
  m <- min(length(a), length(b))
  paired <- if (m > 0) max(abs(a[seq_len(m)] - b[seq_len(m)])) else 0
  extra <- 0
  if (length(a) > m) extra <- max(extra, max(a[(m + 1):length(a)]) / 2)
  if (length(b) > m) extra <- max(extra, max(b[(m + 1):length(b)]) / 2)
  max(paired, extra)
}
