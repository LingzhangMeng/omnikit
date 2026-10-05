# omnikit — 03_program.R  (Gap G4: cross-species program scoring & conservation)
# docs/MATH.md §2.

#' UCell (Andreatta & Carmona 2021), re-implemented. Rank-based per-cell score in [0,1].
#' The bound is maxU = k*N with N = min(max_rank, G)  (NOT the distinct-rank sum).
#' @param X cells x genes matrix (counts or log-normalised; ranks are per cell).
#' @param gene_sets Named list of character gene vectors (or a single character vector).
#' @param max_rank Rank cap (UCell default 1500).
#' @param seed RNG seed for the random tie-break.
#' @param chunk Rows processed per block.
#' @return A cells x \code{length(gene_sets)} matrix of UCell scores in [0,1].
omni_ucell <- function(X, gene_sets, max_rank = 1500L, seed = 0L, chunk = 2000L) {
  X <- as.matrix(omni_as_matrix(X))
  G <- ncol(X); N <- min(max_rank, G)
  genes <- colnames(X); if (is.null(genes)) genes <- as.character(seq_len(G))
  if (is.character(gene_sets)) gene_sets <- list(signature = gene_sets)
  idx_list <- lapply(gene_sets, function(s) which(genes %in% s))
  k_list <- vapply(idx_list, length, integer(1))
  out <- matrix(NA_real_, nrow(X), length(gene_sets), dimnames = list(rownames(X), names(gene_sets)))
  set.seed(seed)
  n <- nrow(X)
  for (b in seq(1L, n, by = chunk)) {
    rows <- b:min(b + chunk - 1L, n)
    M <- X[rows, , drop = FALSE]
    M <- M + matrix(stats::runif(length(M), 0, 1e-4), nrow = nrow(M))  # random tie-break
    R <- t(apply(-M, 1L, rank, ties.method = "first"))                 # rank by descending expression
    R <- pmin(R, max_rank)
    for (j in seq_along(idx_list)) {
      k <- k_list[j]
      if (k == 0L) { out[rows, j] <- NA_real_; next }
      U <- rowSums(R[, idx_list[[j]], drop = FALSE])
      minU <- k * (k + 1) / 2; maxU <- k * N
      out[rows, j] <- (maxU - U) / (maxU - minU)
    }
  }
  out
}

#' Signed UCell programme score = UCell(up) - UCell(down).
#' @param X cells x genes matrix.
#' @param up,down Character vectors of up- and down-programme genes.
#' @param ... Passed to \code{omni_ucell}.
#' @return A numeric vector (UCell(up) - UCell(down)).
omni_ucell_signed <- function(X, up, down, ...) {
  s <- omni_ucell(X, list(up = up, down = down), ...)
  s[, "up"] - s[, "down"]
}

#' Control-bin (expression-matched) score, scanpy/Pagès style (docs/MATH.md §2.2).
#' @param X cells x genes matrix (log-normalised).
#' @param gene_set Character vector of signature genes.
#' @param ctrl_size Number of control genes per bin.
#' @param n_bins Number of expression bins.
#' @param seed RNG seed for control-gene sampling.
#' @param zscore Z-score the per-bin scores (otherwise return the mean).
#' @return A numeric vector of scores.
omni_score_genes <- function(X, gene_set, ctrl_size = 50L, n_bins = 25L,
                             seed = 0L, zscore = FALSE) {
  X <- as.matrix(omni_as_matrix(X))
  genes <- colnames(X); if (is.null(genes)) genes <- as.character(seq_len(ncol(X)))
  gene_mean <- colMeans(X)
  qs <- unique(stats::quantile(gene_mean, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  if (length(qs) < 2) qs <- range(gene_mean)
  bins <- cut(gene_mean, breaks = qs, include.lowest = TRUE, labels = FALSE)
  set.seed(seed)
  sig <- which(genes %in% gene_set)
  bin_scores <- list()
  for (b in sort(unique(bins))) {
    inb <- which(bins == b); sig_b <- intersect(sig, inb)
    if (length(sig_b) == 0L) next
    ctrl_pool <- setdiff(inb, sig_b); nc <- min(ctrl_size, length(ctrl_pool))
    if (nc == 0L) next
    ctrl <- sample(ctrl_pool, nc)
    bin_scores[[length(bin_scores) + 1L]] <-
      rowMeans(X[, sig_b, drop = FALSE]) - rowMeans(X[, ctrl, drop = FALSE])
  }
  if (length(bin_scores) == 0L) return(rep(NA_real_, nrow(X)))
  M <- do.call(cbind, bin_scores)
  score <- rowMeans(M)
  if (zscore) score <- score / pmax(apply(M, 1L, stats::sd), 1e-9)
  score
}

#' Diagnostic: Spearman(score, log library size). Report this for every score.
#' @param score Numeric score.
#' @param lib Numeric library sizes.
#' @return A list with \code{spearman_log_lib} and \code{n}.
omni_score_diagnostics <- function(score, lib) {
  ok <- is.finite(score) & is.finite(lib)
  list(spearman_log_lib = stats::cor(score[ok], log1p(lib[ok]), method = "spearman"),
       n = sum(ok))
}

#' Cross-species conservation call: per-gene Fisher combine + sign consistency + thresholds.
#' @param df data.frame of per-gene statistics.
#' @param t_col Column name of the human effect (t).
#' @param rho_col Column name of the mouse trend (rho).
#' @param ph_col,pm_col Column names of the human / mouse p-values.
#' @param sign_col Optional column with the sign-consistency flag.
#' @param t0,rho0 Effect-size thresholds.
#' @param alpha FDR threshold.
#' @param gene_col Column name of the gene identifier.
#' @return A data.frame ordered by combined q, with \code{conserved} flags.
omni_conservation <- function(df, t_col = "human_t", rho_col = "mouse_rho",
                              ph_col = "human_p", pm_col = "mouse_p",
                              sign_col = NULL, t0 = 2, rho0 = 0.27, alpha = 0.05,
                              gene_col = "gene") {
  t <- df[[t_col]]; rho <- df[[rho_col]]
  ph <- df[[ph_col]]; pm <- df[[pm_col]]
  X2 <- -2 * (log(ph) + log(pm))
  p_comb <- stats::pchisq(X2, df = 4, lower.tail = FALSE)
  q <- omni_bh(p_comb)
  sign_ok <- if (!is.null(sign_col)) as.logical(df[[sign_col]]) else sign(t) == sign(rho)
  out <- data.frame(gene = df[[gene_col]], human_t = t, mouse_rho = rho,
                    combined_p = p_comb, combined_q = q,
                    sign_consistent = sign_ok)
  out$conserved <- with(out, sign_consistent & abs(human_t) >= t0 &
                          abs(mouse_rho) >= rho0 & combined_q < alpha)
  out[order(out$combined_q), ]
}

#' Carrier inference: per-cell-type AUC vs rest (+ BH), and a donor-aware within-donor contrast.
#' @param scores Numeric score.
#' @param celltype Character cell-type labels.
#' @param donor Optional donor labels (enables the donor-aware contrast).
#' @param min_cells Minimum cells per type to test.
#' @return A data.frame of per-cell-type AUC, p, q and (optionally) donor contrast.
omni_carrier <- function(scores, celltype, donor = NULL, min_cells = 10L) {
  scores <- as.numeric(scores); celltype <- as.character(celltype)
  rows <- list()
  for (ct in unique(celltype)) {
    m <- celltype == ct
    if (sum(m) < min_cells || sum(m) == length(m)) next
    r <- data.frame(cell_type = ct, n_cells = sum(m), mean_score = mean(scores[m]),
                    auc_vs_rest = omni_auc(scores, m),
                    p_mwu = omni_mwu(scores[m], scores[!m])$p.value,
                    donor_diff_mean = NA_real_, donor_diff_median = NA_real_, donor_p = NA_real_)
    if (!is.null(donor)) {
      d <- as.character(donor)
      diffs <- vapply(split(seq_along(scores), d), function(ix) {
        own <- ix[celltype[ix] == ct]; oth <- ix[celltype[ix] != ct]
        if (length(own) == 0L || length(oth) == 0L) return(NA_real_)
        mean(scores[own]) - mean(scores[oth])
      }, numeric(1))
      diffs <- diffs[is.finite(diffs)]
      if (length(diffs)) {
        r$donor_diff_mean <- mean(diffs); r$donor_diff_median <- stats::median(diffs)
        r$donor_p <- if (length(diffs) >= 3L) stats::wilcox.test(diffs)$p.value else NA_real_
      }
    }
    rows[[ct]] <- r
  }
  res <- do.call(rbind, rows)
  res$q_BH <- omni_bh(res$p_mwu)
  res[order(-res$auc_vs_rest), ]
}
