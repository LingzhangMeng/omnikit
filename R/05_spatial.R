# omnikit — 05_spatial.R  (Gap G2: spatial niche statistics)
# docs/MATH.md §4. Weights W are row-normalised kNN (default) or fixed-radius; symmetric.

#' Symmetric kNN (or fixed-radius) spatial graph; row-normalised by default.
#' @param coords Numeric N x 2 coordinate matrix.
#' @param k Number of nearest neighbours.
#' @param radius Fixed radius (overrides k when given).
#' @param row_normalize Row-normalise the weight matrix.
#' @return A symmetric (row-normalised) sparse weight matrix.
omni_knn_graph <- function(coords, k = 6L, radius = NULL, row_normalize = TRUE) {
  coords <- as.matrix(coords); n <- nrow(coords)
  if (is.null(radius)) {
    nn <- .omni_knn(coords, coords, k, exclude_self = TRUE)
    i <- rep(seq_len(n), each = k); j <- as.vector(t(nn))
    W <- Matrix::sparseMatrix(i, j, x = 1, dims = c(n, n))
  } else {
    D <- as.matrix(stats::dist(coords))
    W <- Matrix::Matrix(1 * (D <= radius), sparse = TRUE); Matrix::diag(W) <- 0
  }
  W <- Matrix::forceSymmetric(W, uplo = "U")
  if (row_normalize) {
    rs <- Matrix::rowSums(W); rs[rs == 0] <- 1
    W <- Matrix::Diagonal(x = 1 / rs) %*% W
  }
  W
}

#' Global Moran's I with a permutation null. E[I] = -1/(n-1).
#' @param x Numeric variable.
#' @param W Spatial weight matrix.
#' @param n_perm Number of permutations.
#' @param seed RNG seed.
#' @return A list with \code{I}, \code{EI}, \code{Z}, \code{p_perm}, \code{null_sd}.
omni_moran <- function(x, W, n_perm = 999L, seed = 0L) {
  x <- as.numeric(x); n <- length(x); S0 <- sum(W)
  stat <- function(v) { z <- v - mean(v); (n / S0) * as.numeric(crossprod(z, W %*% z)) / sum(z^2) }
  I <- stat(x)
  set.seed(seed); null <- replicate(n_perm, stat(sample(x)))
  list(I = I, EI = -1 / (n - 1), Z = (I - mean(null)) / stats::sd(null),
       p_perm = omni_perm_p(I, null, "two.sided"), null_sd = stats::sd(null))
}

#' Global Geary's C. E[C] = 1.
#' @param x Numeric variable.
#' @param W Spatial weight matrix.
#' @param n_perm Number of permutations.
#' @param seed RNG seed.
#' @return A list with \code{C}, \code{EC}, \code{p_perm}.
omni_geary <- function(x, W, n_perm = 999L, seed = 0L) {
  x <- as.numeric(x); n <- length(x); S0 <- sum(W)
  stat <- function(v) {
    z <- v - mean(v); rs <- as.numeric(Matrix::rowSums(W)); cs <- as.numeric(Matrix::colSums(W))
    num <- sum(v^2 * (rs + cs)) - 2 * as.numeric(crossprod(v, W %*% v))
    ((n - 1) / (2 * S0)) * num / sum(z^2)
  }
  C <- stat(x); set.seed(seed); null <- replicate(n_perm, stat(sample(x)))
  list(C = C, EC = 1, p_perm = omni_perm_p(C, null, "two.sided"))
}

#' Local Moran's I_i (LISA) with an (approximate) global permutation null per cell.
#' @param x Numeric variable.
#' @param W Spatial weight matrix.
#' @param n_perm Number of permutations.
#' @param seed RNG seed.
#' @return A data.frame with local \code{Ii}, \code{lag}, \code{p_perm}, \code{quad}.
omni_lisa <- function(x, W, n_perm = 999L, seed = 0L) {
  x <- as.numeric(x); z <- x - mean(x)
  lag <- as.numeric(W %*% z); Ii <- z * lag
  set.seed(seed)
  null <- vapply(seq_len(n_perm), function(b) { zs <- sample(z); zs * as.numeric(W %*% zs) },
                 numeric(length(z)))
  p <- vapply(seq_along(Ii), function(i) omni_perm_p(Ii[i], null[i, ], "two.sided"), numeric(1))
  data.frame(Ii = Ii, lag = lag, p_perm = p, quad = ifelse(Ii > 0, "HH/LL", "HL/LH"))
}

#' Neighbourhood enrichment: cell-type co-occurrence z-scores on the graph.
#' @param labels Character cell-type labels.
#' @param W Spatial weight matrix.
#' @param n_perm Number of permutations.
#' @param seed RNG seed.
#' @return A list with \code{observed} and \code{zscore} co-occurrence matrices.
omni_nhood_enrichment <- function(labels, W, n_perm = 999L, seed = 0L) {
  labels <- as.character(labels); lv <- sort(unique(labels)); K <- length(lv)
  Wt <- Matrix::triu(W); s <- Matrix::summary(Wt)          # i<j edges (weight > 0)
  a <- match(labels[s$i], lv); b <- match(labels[s$j], lv)
  obs <- matrix(0, K, K, dimnames = list(lv, lv))
  for (e in seq_along(a)) { obs[a[e], b[e]] <- obs[a[e], b[e]] + 1; obs[b[e], a[e]] <- obs[b[e], a[e]] + 1 }
  set.seed(seed)
  null <- array(0, c(K, K, n_perm))
  for (p in seq_len(n_perm)) {
    lb <- sample(labels); aa <- match(lb[s$i], lv); bb <- match(lb[s$j], lv)
    nm <- matrix(0, K, K); for (e in seq_along(aa)) { nm[aa[e], bb[e]] <- nm[aa[e], bb[e]] + 1; nm[bb[e], aa[e]] <- nm[bb[e], aa[e]] + 1 }
    null[, , p] <- nm
  }
  mu <- apply(null, c(1, 2), mean); sd_ <- apply(null, c(1, 2), stats::sd)
  z <- (obs - mu) / pmax(sd_, 1e-9)
  list(observed = obs, zscore = z)
}

#' Ripley's K and L for 2-D point pattern (optional mark subset), border-corrected.
#' O(n^2) via a distance matrix; if n > max_n the pattern is randomly subsampled (seed 0).
#' @param coords Numeric N x 2 coordinate matrix.
#' @param radii Radii at which to evaluate (defaulted from the bounding box).
#' @param mark Optional logical/0-1 subset of points.
#' @param max_n Subsample size guard (default 20000).
#' @return A data.frame with \code{radius}, \code{K}, \code{L}.
omni_ripley <- function(coords, radii = NULL, mark = NULL, max_n = 20000L) {
  coords <- as.matrix(coords)
  if (!is.null(mark)) coords <- coords[as.logical(mark), , drop = FALSE]
  if (nrow(coords) > max_n) { set.seed(0); coords <- coords[sample(nrow(coords), max_n), , drop = FALSE] }
  n <- nrow(coords)
  rng <- apply(coords, 2, range); A <- prod(rng[2, ] - rng[1, ])
  if (is.null(radii)) radii <- seq(0, 0.25 * sqrt(A), length.out = 30)[-1]
  D <- as.matrix(stats::dist(coords))
  inside <- outer(coords[, 1], rng[1, 1] + radii, ">=") & outer(coords[, 1], rng[2, 1] - radii, "<=") &
            outer(coords[, 2], rng[1, 2] + radii, ">=") & outer(coords[, 2], rng[2, 2] - radii, "<=")
  K <- vapply(seq_along(radii), function(k) {
    cnt <- sum((D <= radii[k]) & (D > 0))
    wt <- sum(inside[, k]); frac <- if (n > 0) wt / n else 1
    (A / (n * (n - 1))) * (cnt / max(frac, 1e-9))
  }, numeric(1))
  data.frame(radius = radii, K = K, L = sqrt(K / pi) - radii)
}

#' Spatial ligand-receptor co-expression on the graph (permutation null).
#' @param expr cells x genes expression matrix.
#' @param lr_pairs List of c(ligand, receptor) pairs (or a 2-column data.frame).
#' @param W Spatial weight matrix.
#' @param n_perm Number of permutations.
#' @param seed RNG seed.
#' @return A data.frame of ligand-receptor pairs with \code{gamma} and \code{p_perm}.
omni_spatial_lr <- function(expr, lr_pairs, W, n_perm = 999L, seed = 0L) {
  expr <- as.matrix(omni_as_matrix(expr)); genes <- colnames(expr)
  if (is.data.frame(lr_pairs)) lr_pairs <- split(as.matrix(lr_pairs), seq_len(nrow(lr_pairs)))
  S0 <- sum(W)
  gamma <- function(s) as.numeric(crossprod(s, W %*% s)) / S0
  set.seed(seed)
  res <- do.call(rbind, lapply(lr_pairs, function(p) {
    L <- p[[1]]; R <- p[[2]]
    if (!(L %in% genes) || !(R %in% genes)) return(data.frame(ligand = L, receptor = R, gamma = NA, p_perm = NA))
    s <- 1 * (expr[, L] > 0) * 1 * (expr[, R] > 0)
    g <- gamma(s); null <- replicate(n_perm, gamma(sample(s)))
    data.frame(ligand = L, receptor = R, gamma = g, p_perm = omni_perm_p(g, null, "greater"))
  }))
  res[order(-res$gamma), ]
}
