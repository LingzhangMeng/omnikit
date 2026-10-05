# omnikit — 06_integrate.R  (Gap G1: cross-omics integration)
# docs/MATH.md §5.

#' Multi-view NMF with a shared sample-factor matrix W (Liu et al. 2013 multiNMF).
#' X^(v) ~ W H^(v) : W (samples x k) shared, H^(v) (k x features) view-specific.
#' @param views List of non-negative matrices with the SAME number of rows (samples/cells).
#' @param k Number of shared factors.
#' @param n_iter Maximum iterations.
#' @param tol Convergence tolerance.
#' @param seed RNG seed.
#' @return A list with \code{W} (samples x k), view-specific \code{H}, \code{objective}, \code{k}.
omni_joint_nmf <- function(views, k = 5L, n_iter = 200L, tol = 1e-5, seed = 0L) {
  Xs <- lapply(views, function(X) { X <- as.matrix(omni_as_matrix(X)); X[X < 0] <- 0; X })
  n <- nrow(Xs[[1]]); V <- length(Xs)
  set.seed(seed)
  W <- matrix(stats::runif(n * k, 0, 1), n, k)
  Hs <- lapply(Xs, function(X) matrix(stats::runif(k * ncol(X), 0, 1), k, ncol(X)))
  obj <- numeric(n_iter)
  for (it in seq_len(n_iter)) {
    numW <- matrix(0, n, k); denW <- matrix(0, n, k)
    for (v in seq_len(V)) { numW <- numW + Xs[[v]] %*% t(Hs[[v]]); denW <- denW + W %*% (Hs[[v]] %*% t(Hs[[v]])) }
    W <- W * (numW / pmax(denW, 1e-10))
    for (v in seq_len(V)) { Hs[[v]] <- Hs[[v]] * ((t(W) %*% Xs[[v]]) / pmax(t(W) %*% W %*% Hs[[v]], 1e-10)) }
    obj[it] <- sum(vapply(seq_len(V), function(v) sum((Xs[[v]] - W %*% Hs[[v]])^2), numeric(1)))
    if (it > 1L && abs(obj[it - 1L] - obj[it]) / pmax(obj[it - 1L], 1e-10) < tol) { obj <- obj[seq_len(it)]; break }
  }
  list(W = W, H = Hs, objective = obj, k = k)
}

#' Canonical correlation analysis: canonical correlations + weights (whitening SVD).
#' @param X,Y Two numeric matrices with the same number of rows (samples).
#' @param k Number of canonical components (default: all).
#' @return A list with \code{cor} (canonical correlations) and \code{x_weights}, \code{y_weights}.
omni_cca <- function(X, Y, k = NULL) {
  X <- scale(as.matrix(X), center = TRUE, scale = FALSE)
  Y <- scale(as.matrix(Y), center = TRUE, scale = FALSE)
  n1 <- nrow(X) - 1
  Sxx <- crossprod(X) / n1; Syy <- crossprod(Y) / n1; Sxy <- crossprod(X, Y) / n1
  ex <- eigen(Sxx, symmetric = TRUE); ey <- eigen(Syy, symmetric = TRUE)
  Wh <- function(e) { w <- 1 / sqrt(pmax(e$values, 1e-12)); e$vectors %*% (w * t(e$vectors)) }
  Wx <- Wh(ex); Wy <- Wh(ey)
  sv <- svd(Wx %*% Sxy %*% Wy)
  kk <- if (is.null(k)) length(sv$d) else min(k, length(sv$d))
  list(cor = sv$d[seq_len(kk)],
       x_weights = Wx %*% sv$u[, seq_len(kk)],
       y_weights = Wy %*% sv$v[, seq_len(kk)])
}

#' RV coefficient: global similarity of two centered configurations (0..1).
#' @param A,B Two centered configurations (same number of rows).
#' @return The RV coefficient (0..1).
omni_rv <- function(A, B) {
  A <- scale(as.matrix(A), center = TRUE, scale = FALSE)
  B <- scale(as.matrix(B), center = TRUE, scale = FALSE)
  Sxy <- crossprod(A, B); Sxx <- crossprod(A); Syy <- crossprod(B)
  as.numeric(sum(Sxy * Sxy)) / sqrt(sum(Sxx * Sxx) * sum(Syy * Syy))
}

#' Procrustes rotation/scale aligning X to Y (SVD); returns rotation, scale, residual.
#' @param X,Y Two configurations (same number of rows).
#' @return A list with \code{rotation}, \code{scale}, \code{residual}.
omni_procrustes <- function(X, Y) {
  X <- scale(as.matrix(X), center = TRUE, scale = FALSE)
  Y <- scale(as.matrix(Y), center = TRUE, scale = FALSE)
  sv <- svd(crossprod(X, Y))
  R <- sv$u %*% t(sv$v)
  list(rotation = R, scale = sum(sv$d) / sum(X^2),
       residual = 1 - sum(sv$d)^2 / (sum(X^2) * sum(Y^2)))
}

#' MOFA-lite: multi-view probabilistic factor analysis (EM), k latent factors.
#' X^(v) = Z Lambda^(v)^T + eps;  E[Z|X] via precision-weighted combination.
#' @param views List of numeric matrices (same rows).
#' @param k Number of latent factors.
#' @param n_iter Maximum EM iterations.
#' @param tol Convergence tolerance.
#' @param seed RNG seed.
#' @return A list with \code{factors}, \code{loadings}, \code{sigma2}, \code{var_explained}, \code{loglik}, \code{k}.
omni_mofa_lite <- function(views, k = 5L, n_iter = 200L, tol = 1e-6, seed = 0L) {
  Xs <- lapply(views, function(X) as.matrix(omni_as_matrix(X)))
  n <- nrow(Xs[[1]]); V <- length(Xs)
  set.seed(seed)
  Lam <- lapply(Xs, function(X) matrix(stats::rnorm(k * ncol(X), 0, 0.1), k, ncol(X)))
  sig2 <- vapply(Xs, function(X) 0.5 * stats::var(as.numeric(X)), numeric(1))
  ll <- numeric(n_iter); Ez <- matrix(0, n, k)
  for (it in seq_len(n_iter)) {
    M <- diag(k)
    for (v in seq_len(V)) M <- M + tcrossprod(Lam[[v]]) / sig2[v]      # (k x p)(p x k) = k x k
    Minv <- solve(M + diag(k) * 1e-8)
    Ez <- matrix(0, n, k)
    for (v in seq_len(V)) Ez <- Ez + Xs[[v]] %*% t(Lam[[v]]) / sig2[v]
    Ez <- Ez %*% Minv
    Ezz <- n * Minv + crossprod(Ez)
    for (v in seq_len(V)) {
      Lam[[v]] <- solve(Ezz) %*% crossprod(Ez, Xs[[v]])             # (k x k)(k x p) = k x p
      sig2[v] <- max(1e-8, sum((Xs[[v]] - Ez %*% Lam[[v]])^2) / (n * ncol(Xs[[v]])))
    }
    ll[it] <- sum(vapply(seq_len(V), function(v)
      -n * ncol(Xs[[v]]) / 2 * log(sig2[v]) -
        sum((Xs[[v]] - Ez %*% Lam[[v]])^2) / (2 * sig2[v]), numeric(1)))
    if (it > 1L && abs(ll[it] - ll[it - 1L]) / pmax(abs(ll[it - 1L]), 1e-9) < tol) { ll <- ll[seq_len(it)]; break }
  }
  ve <- vapply(seq_len(V), function(v) 1 - sig2[v] / stats::var(as.numeric(Xs[[v]])), numeric(1))
  list(factors = Ez, loadings = Lam, sigma2 = sig2, var_explained = ve, loglik = ll, k = k)
}

#' scDRS-flavoured program-GWAS association: standardise an observed program score
#' against a matched control-gene-set score distribution.
#' @param program_score Numeric observed programme scores.
#' @param control_score Numeric matched-control scores (permutation / matched gene sets).
#' @return A list with \code{stat}, \code{z}, \code{p}, \code{n_control}.
omni_program_gwas <- function(program_score, control_score = NULL) {
  obs <- mean(program_score, na.rm = TRUE)
  if (is.null(control_score)) return(list(stat = obs, z = NA_real_, p = NA_real_, n_control = 0L))
  z <- (obs - mean(control_score)) / stats::sd(control_score)
  list(stat = obs, z = z, p = 2 * stats::pnorm(-abs(z)), n_control = length(control_score))
}
