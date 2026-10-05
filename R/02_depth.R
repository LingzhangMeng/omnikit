# omnikit — 02_depth.R
# Depth control and ortholog handling (docs/MATH.md §1, §0).

#' Regress a score on log-library size and return the depth-controlled residual.
#' Fits per group when `group` is supplied (recommended), else globally.
#' @param score Numeric score.
#' @param log_lib Numeric log library size.
#' @param group Optional grouping (control fits per group).
#' @return A data.frame (score, log_lib, group, fit, resid) with attributes r_before, r_after.
omni_depth_control <- function(score, log_lib, group = NULL) {
  score <- as.numeric(score); log_lib <- as.numeric(log_lib)
  n <- length(score)
  if (is.null(group)) group <- rep("all", n)
  group <- as.character(group)
  resid <- rep(NA_real_, n); fit <- rep(NA_real_, n)
  for (g in unique(group)) {
    idx <- which(group == g & is.finite(score) & is.finite(log_lib))
    if (length(idx) < 3L) next
    m <- stats::lm(score[idx] ~ log_lib[idx])
    resid[idx] <- stats::residuals(m)
    fit[idx] <- stats::fitted(m)
  }
  keep <- is.finite(score) & is.finite(log_lib)
  r_before <- if (sum(keep) > 2) stats::cor(score[keep], log_lib[keep], method = "spearman") else NA_real_
  k2 <- keep & is.finite(resid)
  r_after <- if (sum(k2) > 2) stats::cor(resid[k2], log_lib[k2], method = "spearman") else NA_real_
  out <- data.frame(score = score, log_lib = log_lib, group = group,
                    fit = fit, resid = resid)
  attr(out, "r_before") <- r_before
  attr(out, "r_after") <- r_after
  out
}

#' Depth-stratified (matched) AUC: per-quantile-bin AUC vs rest + size-weighted mean.
#' This is the correct correction when the score-depth relation is non-linear.
#' @param score Numeric score.
#' @param log_lib Numeric log library size.
#' @param pos Logical positive mask.
#' @param nq Number of depth quantile bins.
#' @return A list with \code{per_bin}, \code{auc_matched}, \code{auc_min}, \code{auc_raw}.
omni_depth_stratified_auc <- function(score, log_lib, pos, nq = 5L) {
  score <- as.numeric(score); log_lib <- as.numeric(log_lib); pos <- as.logical(pos)
  keep <- is.finite(score) & is.finite(log_lib)
  score <- score[keep]; log_lib <- log_lib[keep]; pos <- pos[keep]
  q <- stats::quantile(log_lib, probs = seq(0, 1, length.out = nq + 1L), na.rm = TRUE)
  q <- unique(q)
  bin <- cut(log_lib, breaks = q, include.lowest = TRUE, labels = FALSE)
  rows <- do.call(rbind, lapply(sort(unique(bin)), function(b) {
    m <- bin == b
    data.frame(bin = b, n = sum(m), auc = omni_auc(score[m], pos[m]))
  }))
  wmean <- stats::weighted.mean(rows$auc, rows$n, na.rm = TRUE)
  list(per_bin = rows, auc_matched = wmean, auc_min = min(rows$auc, na.rm = TRUE),
       auc_raw = omni_auc(score, pos))
}

#' Case-fold symbol map: fold -> representative original symbol (drops case collisions).
#' @param symbols Character vector of gene symbols.
#' @return A list with \code{map} (fold -> symbol) and \code{collisions}.
omni_casefold_map <- function(symbols) {
  symbols <- as.character(symbols)
  m <- list(); amb <- character(0)
  for (g in symbols) {
    k <- toupper(g)
    if (!is.null(m[[k]]) && m[[k]] != g) amb <- c(amb, k) else m[[k]] <- g
  }
  for (k in unique(amb)) m[[k]] <- NULL
  list(map = m, collisions = unique(amb))
}

#' Collapse many-to-many ortholog pairs to ONE representative mouse gene per human gene.
#' Preference: case-fold-identical symbol > highest optional `score`.
#' @param pairs data.frame(human_symbol, mouse_symbol, score[optional])
#' @return A two-column data.frame (human_symbol, mouse_symbol), one row per human gene.
omni_collapse_orthologs <- function(pairs) {
  stopifnot(all(c("human_symbol", "mouse_symbol") %in% names(pairs)))
  has_score <- "score" %in% names(pairs)
  ident <- toupper(pairs$human_symbol) == toupper(pairs$mouse_symbol)
  key <- data.frame(human_symbol = pairs$human_symbol, mouse_symbol = pairs$mouse_symbol,
                    ident = ident, score = if (has_score) pairs$score else 0,
                    stringsAsFactors = FALSE)
  key <- key[order(key$human_symbol, -key$ident, -key$score), ]
  key <- key[!duplicated(key$human_symbol), c("human_symbol", "mouse_symbol")]
  rownames(key) <- NULL
  key
}
