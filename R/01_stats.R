# omnikit — 01_stats.R
# Shared estimators, implemented from their definitions (see docs/MATH.md §0).

#' AUC of a score against a positive mask, via the normalised Mann-Whitney U.
#' AUC = P(x_pos > x_neg) + 0.5 P(tie).
#' @param x Numeric scores.
#' @param pos Logical positive mask, same length as \code{x}.
#' @return The AUC (scalar).
omni_auc <- function(x, pos) {
  x <- as.numeric(x); pos <- as.logical(pos)
  n1 <- sum(pos); n0 <- sum(!pos)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(x, ties.method = "average")
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Mann-Whitney U p-value (normal approximation with tie correction).
#' @param x,y Numeric vectors to compare.
#' @param alternative One of "two.sided", "greater", "less".
#' @return A list with \code{statistic} (U), \code{z} and \code{p.value}.
omni_mwu <- function(x, y, alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  x <- as.numeric(x); y <- as.numeric(y)
  n1 <- length(x); n0 <- length(y)
  if (n1 == 0L || n0 == 0L) return(list(statistic = NA_real_, p.value = NA_real_))
  r <- rank(c(x, y), ties.method = "average")
  U <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  mu <- n1 * n0 / 2
  tt <- table(c(x, y))
  tie <- sum(tt^3 - tt)
  N <- n1 + n0
  sig2 <- (n1 * n0 / 12) * ((N + 1) - tie / (N * (N - 1)))
  z <- (U - mu) / sqrt(sig2)
  p <- switch(alternative,
              two.sided = 2 * stats::pnorm(-abs(z)),
              greater   = stats::pnorm(z, lower.tail = FALSE),
              less      = stats::pnorm(z))
  list(statistic = U, z = z, p.value = p)
}

#' DeLong AUC with a variance/CI (fast midrank form; Sun & Xu 2014).
#' @param scores Numeric scores.
#' @param pos Logical positive mask.
#' @param conf.level Confidence level for the interval.
#' @return A list with \code{auc}, \code{se}, \code{lower}, \code{upper}, \code{n_pos}, \code{n_neg}.
omni_delong <- function(scores, pos, conf.level = 0.95) {
  scores <- as.numeric(scores); pos <- as.logical(pos)
  x <- scores[pos]; y <- scores[!pos]
  m <- length(x); n <- length(y)
  if (m < 2L || n < 2L) return(list(auc = omni_auc(scores, pos), se = NA_real_,
                                    lower = NA_real_, upper = NA_real_, n_pos = m, n_neg = n))
  R  <- rank(c(x, y), ties.method = "average")
  Rx <- R[seq_len(m)]; Ry <- R[m + seq_len(n)]
  rxx <- rank(x, ties.method = "average"); ryy <- rank(y, ties.method = "average")
  V10 <- (Rx - rxx) / n            # placement values for positives
  V01 <- (Ry - ryy) / m            # placement values for negatives
  auc <- mean(V10)
  S10 <- sum((V10 - auc)^2) / (m - 1)
  S01 <- sum((V01 - auc)^2) / (n - 1)
  se  <- sqrt(S10 / m + S01 / n)
  z   <- stats::qnorm(1 - (1 - conf.level) / 2)
  list(auc = auc, se = se, lower = auc - z * se, upper = auc + z * se, n_pos = m, n_neg = n)
}

#' Weighted Stouffer meta-analysis of z-scores -> global Z and two-sided p.
#' @param z Numeric vector of z-scores.
#' @param w Optional numeric weights.
#' @return A list with the combined \code{Z}, two-sided \code{p} and \code{k}.
omni_stouffer <- function(z, w = NULL) {
  z <- as.numeric(z); z <- z[is.finite(z)]
  if (length(z) == 0L) return(list(Z = NA_real_, p = NA_real_, k = 0L))
  if (is.null(w)) w <- rep(1, length(z))
  w <- as.numeric(w)[seq_along(z)]
  Z <- sum(w * z) / sqrt(sum(w^2))
  list(Z = Z, p = 2 * stats::pnorm(-abs(Z)), k = length(z))
}

#' Empirical-Bayes (method-of-moments) shrinkage of estimates toward the grand mean.
#' @param beta Numeric estimates.
#' @param se Numeric standard errors.
#' @return Shrunk estimates (with attribute \code{tau2}).
omni_eb_shrink <- function(beta, se) {
  beta <- as.numeric(beta); se <- as.numeric(se)
  ok <- is.finite(beta) & is.finite(se) & se > 0
  b <- beta[ok]; s <- se[ok]
  bbar <- stats::median(b)
  tau2 <- max(0, stats::var(b) - mean(s^2))
  w <- if (tau2 > 0) tau2 / (tau2 + s^2) else 0
  out <- rep(NA_real_, length(beta))
  out[ok] <- bbar + w * (b - bbar)
  attr(out, "tau2") <- tau2; out
}

#' Fisher's method to combine p-values; optional Brown correction given the
#' correlation matrix R of the underlying test statistics.
#' @param p Numeric vector of p-values.
#' @param R Optional correlation matrix of the underlying statistics (Brown correction).
#' @return A list with \code{X2}, \code{p}, \code{df}, \code{k} (and Brown variants when \code{R} is given).
omni_fisher <- function(p, R = NULL) {
  p <- as.numeric(p); p <- p[is.finite(p) & p > 0 & p <= 1]
  k <- length(p)
  if (k == 0L) return(list(X2 = NA_real_, p = NA_real_, df = 0L, k = 0L))
  X2 <- sum(-2 * log(p))
  out <- list(X2 = X2, p = stats::pchisq(X2, df = 2 * k, lower.tail = FALSE), df = 2L * k, k = k)
  if (!is.null(R) && k > 1L) {
    R <- as.matrix(R)
    # Under H0, Z_i = -2 ln P_i ~ chi2_2 (mean 2, var 4) and Cov(Z_i, Z_j) ~ 4 R_ij (Brown 1975).
    V <- 4 * k + 8 * sum(R[lower.tri(R)])
    X2b <- (X2 - 2 * k) * sqrt(2 * k) / sqrt(V) + 2 * k
    out$X2_brown <- X2b
    out$p_brown <- stats::pchisq(X2b, df = 2 * k, lower.tail = FALSE)
  }
  out
}

#' Cohen's d (pooled SD) and Welch's t-test.
#' @param x,y Numeric vectors.
#' @return Cohen's d (pooled SD).
omni_cohen_d <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  n1 <- length(x); n0 <- length(y)
  sp <- sqrt(((n1 - 1) * stats::var(x) + (n0 - 1) * stats::var(y)) / (n1 + n0 - 2))
  (mean(x) - mean(y)) / sp
}

#' Welch's t-test (unequal variances).
#' @param x,y Numeric vectors.
#' @return The \code{htest} object from \code{stats::t.test}.
omni_welch <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  stats::t.test(x, y, var.equal = FALSE)
}

#' Spearman rho (Pearson on tie-averaged ranks).
#' @param x,y Numeric vectors.
#' @return Spearman's rho.
omni_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  stats::cor(rank(x[ok]), rank(y[ok]), method = "pearson")
}

#' Empirical permutation p-value from a null distribution and an observed statistic.
#' @param obs Observed statistic.
#' @param null Numeric null distribution.
#' @param alternative One of "greater", "two.sided".
#' @return The empirical p-value.
omni_perm_p <- function(obs, null, alternative = c("greater", "two.sided")) {
  alternative <- match.arg(alternative)
  null <- null[is.finite(null)]; B <- length(null)
  if (B == 0L) return(NA_real_)
  if (alternative == "greater") (1 + sum(null >= obs)) / (B + 1)
  else {
    m <- mean(null); (1 + sum(abs(null - m) >= abs(obs - m))) / (B + 1)
  }
}
