# omnikit — split & Mondrian conformal prediction primitives (classifier-agnostic).

#' Nonconformity score 1 - p(true class)
#' @param proba numeric matrix (samples x classes), rows sum to 1
#' @param y integer vector of true classes (1-based)
#' @export
omni_nonconformity <- function(proba, y) {
  proba <- as.matrix(proba)
  y <- as.integer(y)
  1 - proba[cbind(seq_along(y), y)]
}

#' Finite-sample conformal quantile of nonconformity scores
#' @param scores numeric vector
#' @param alpha target miscoverage (e.g. 0.1 for 90% coverage)
#' @export
omni_conformal_quantile <- function(scores, alpha) {
  s <- sort(as.numeric(scores))
  n <- length(s)
  k <- min(max(ceiling((n + 1) * (1 - alpha)), 1), n)
  s[k]
}

#' Per-class (Mondrian) conformal thresholds
#' @param proba,samples; y true classes; alpha
#' @export
omni_mondrian_quantiles <- function(proba, y, alpha) {
  proba <- as.matrix(proba); y <- as.integer(y)
  glob <- omni_conformal_quantile(omni_nonconformity(proba, y), alpha)
  q <- sapply(seq_len(ncol(proba)), function(c) {
    s <- omni_nonconformity(proba[y == c, , drop = FALSE], y[y == c])
    if (length(s)) omni_conformal_quantile(s, alpha) else glob
  })
  q
}

#' Prediction-set membership: include class c iff 1 - p_c <= q_c
#' @param proba; q scalar or per-class vector
#' @export
omni_prediction_sets <- function(proba, q) {
  proba <- as.matrix(proba)
  if (length(q) == 1) q <- rep(q, ncol(proba))
  (1 - proba) <= sweep(matrix(q, nrow(proba), ncol(proba), byrow = TRUE), 1, 0) + 1e-12
}

#' Empirical coverage: fraction of samples whose true class is in the set
#' @export
omni_coverage <- function(sets, y) {
  sets <- as.matrix(sets); y <- as.integer(y)
  mean(sets[cbind(seq_along(y), y)])
}

#' Mean prediction-set size
#' @export
omni_avg_set_size <- function(sets) mean(rowSums(as.matrix(sets)))
