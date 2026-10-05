# omnikit — 04_qc.R  (Gap G3: QC-first multi-omics audit)
# docs/MATH.md §3. All metrics are computed from RAW/integer counts.

#' Chunked k-nearest-neighbour search: for each row of A, indices (into B) of the k nearest.
#' Brute force with a partial sort; no external dependency.
#' @noRd
.omni_knn <- function(A, B, k, exclude_self = FALSE, chunk = 2000L) {
  A <- as.matrix(A); B <- as.matrix(B)
  k <- min(k, nrow(B) - as.integer(exclude_self))
  b2 <- rowSums(B^2); out <- matrix(0L, nrow(A), k)
  for (s in seq(1L, nrow(A), by = chunk)) {
    rows <- s:min(s + chunk - 1L, nrow(A))
    Am <- A[rows, , drop = FALSE]
    D2 <- outer(rowSums(Am^2), b2, "+") - 2 * Am %*% t(B)
    for (r in seq_len(nrow(Am))) {
      v <- D2[r, ]; if (exclude_self) v[rows[r]] <- Inf
      out[rows[r], ] <- order(v)[seq_len(k)]          # note: order(partial=) is unusable on R 4.6.1
    }
  }
  out
}

#' Per-cell QC metrics from integer counts. mt/ribo/hb detected by symbol if not supplied.
#' @param X Integer-count cells x genes matrix.
#' @param mt,ribo,hb Optional logical vectors marking mitochondrial / ribosomal / haemoglobin genes (detected by symbol if NULL).
#' @return A data.frame with nCount, nFeature and \%mt/\%ribo/\%hb.
omni_qc_metrics <- function(X, mt = NULL, ribo = NULL, hb = NULL) {
  X <- as.matrix(omni_as_matrix(X))
  genes <- colnames(X); if (is.null(genes)) genes <- as.character(seq_len(ncol(X)))
  up <- toupper(genes)
  if (is.null(mt))   mt   <- grepl("^MT[-.]", up)
  if (is.null(ribo)) ribo <- grepl("^RP[SL]", up)
  if (is.null(hb))   hb   <- grepl("^HB[ABDEGMQZ][0-9]?$", up)
  ncount <- rowSums(X)
  pct <- function(m) if (any(m)) 100 * rowSums(X[, m, drop = FALSE]) / pmax(ncount, 1) else rep(0, nrow(X))
  data.frame(nCount = ncount, nFeature = rowSums(X > 0),
             pct_mt = pct(mt), pct_ribo = pct(ribo), pct_hb = pct(hb))
}

#' MAD-based outlier flag (Gaussian-consistent scale 1.4826*MAD).
#' @param x Numeric vector.
#' @param k Number of MADs for the flag.
#' @return A list with \code{flag}, \code{median}, \code{scale}, \code{lower}, \code{upper}.
omni_mad_outlier <- function(x, k = 3) {
  x <- as.numeric(x); med <- stats::median(x, na.rm = TRUE)
  s <- 1.4826 * stats::mad(x, constant = 1, na.rm = TRUE)
  list(flag = abs(x - med) > k * s, median = med, scale = s,
       lower = med - k * s, upper = med + k * s)
}

#' Fit a 2-component 1-D Gaussian mixture (EM) and return the density-crossing threshold.
#' @param x Numeric vector.
#' @param max_iter Maximum EM iterations.
#' @param tol Convergence tolerance.
#' @return A list with \code{weights}, \code{means}, \code{sds}, \code{threshold}.
omni_gmm2 <- function(x, max_iter = 200L, tol = 1e-7) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  m <- c(stats::quantile(x, 0.25), stats::quantile(x, 0.75))
  s <- c(stats::sd(x), stats::sd(x)) / 2; w <- c(0.5, 0.5)
  for (it in seq_len(max_iter)) {
    d <- vapply(seq_len(2), function(j) w[j] * stats::dnorm(x, m[j], s[j]), numeric(length(x)))
    r <- d / pmax(rowSums(d), 1e-300)
    w0 <- w
    nk <- colSums(r) + 1e-12
    m <- colSums(r * x) / nk
    s <- sqrt(colSums(r * (x - rep(m, each = length(x)))^2) / nk) + 1e-12
    w <- nk / length(x)
    if (max(abs(w - w0)) < tol) break
  }
  # crossing: solve w1*N(m1,s1) = w2*N(m2,s2)  (the point between the two component means)
  lo <- min(m); hi <- max(m)
  f <- function(t) log(w[1]) - log(w[2]) + stats::dnorm(t, m[1], s[1], log = TRUE) -
                        stats::dnorm(t, m[2], s[2], log = TRUE)
  thr <- tryCatch(stats::uniroot(f, c(lo, hi))$root, error = function(e) mean(m))
  list(weights = w, means = m, sds = s, threshold = thr)
}

#' Scrublet-style doublet detection (Wolock 2019), re-implemented (docs/MATH.md §3.2).
#' @param X Raw-count cells x genes matrix.
#' @param expected_rate Expected doublet rate (default 0.8\% per 1000 cells).
#' @param sim_doublets Number of simulated doublets.
#' @param k Number of neighbours for the kNN doublet score.
#' @param n_pcs Number of principal components.
#' @param seed RNG seed.
#' @return A list with \code{call}, \code{score}, \code{threshold}, \code{bimodal} and GMM details.
omni_doublet_scrublet <- function(X, expected_rate = NULL, sim_doublets = NULL,
                                  k = 30L, n_pcs = 30L, seed = 0L) {
  X <- as.matrix(omni_as_matrix(X)); n <- nrow(X)
  if (is.null(expected_rate)) expected_rate <- 0.008 * n / 1000       # 10x rule of thumb
  if (is.null(sim_doublets)) sim_doublets <- n
  set.seed(seed)
  i <- sample(n, sim_doublets, replace = TRUE); j <- sample(n, sim_doublets, replace = TRUE)
  sim <- X[i, , drop = FALSE] + X[j, , drop = FALSE]                   # doublet = sum of two cells
  norm <- function(M) log1p(sweep(M, 1, pmax(rowSums(M), 1), "/") * 1e4)
  allm <- rbind(norm(X), norm(sim))
  v <- apply(allm, 2L, stats::var); hvg <- order(-v)[seq_len(min(2000L, ncol(allm)))]
  Z <- scale(allm[, hvg, drop = FALSE]); Z <- Z[, is.finite(colSums(Z)), drop = FALSE]
  sv <- svd(Z, nu = min(n_pcs, ncol(Z)), nv = 0L)
  pcs <- sv$u %*% diag(sv$d[seq_len(min(n_pcs, ncol(Z)))], min(n_pcs, ncol(Z)))
  nn <- .omni_knn(pcs[seq_len(n), , drop = FALSE], pcs, k, exclude_self = TRUE)
  score <- rowMeans(nn > n)                                            # fraction simulated neighbours
  gm <- omni_gmm2(score)
  sep <- abs(gm$means[2] - gm$means[1]) > 2 * mean(gm$sds)
  thr <- stats::quantile(score, 1 - expected_rate, names = FALSE)   # Scrublet's calibrated default
  list(call = score >= thr, score = score, threshold = thr, threshold_source = "expected_rate",
       gmm_threshold = gm$threshold, bimodal = sep, expected_rate = expected_rate, gmm = gm)
}

#' Batch-mixing score (Local Inverse Simpson Index). 1 = pure, B = fully mixed.
#' @param embedding Numeric embedding (e.g. PCs).
#' @param batch Batch labels.
#' @param k Number of neighbours.
#' @return A list with \code{lisi}, \code{mean} and \code{max} (= number of batches).
omni_mixing <- function(embedding, batch, k = 30L) {
  E <- as.matrix(embedding); b <- factor(as.character(batch)); B <- nlevels(b)
  code <- as.integer(b)
  nn <- .omni_knn(E, E, k, exclude_self = TRUE)
  lisi <- apply(nn, 1L, function(ix) {
    p <- tabulate(code[ix], nbins = B) / length(ix)
    1 / sum(p^2)
  })
  list(lisi = lisi, mean = mean(lisi), max = B)
}
