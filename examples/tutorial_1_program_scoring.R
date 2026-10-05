#!/usr/bin/env Rscript
# examples/tutorial_1_program_scoring.R
# Tutorial 1 — score a gene programme, diagnose/remove the depth confound, find its carriers.
# Self-contained (synthetic). Run from the package root:  Rscript examples/tutorial_1_program_scoring.R
suppressMessages(library(omnikit))
dir.create("docs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs/tutorial_outputs", recursive = TRUE, showWarnings = FALSE)
set.seed(0)

## synthetic single-cell matrix: 3000 cells x 400 genes (sparse, variable sequencing depth)
n_cells <- 3000; n_genes <- 400
genes <- sprintf("G%03d", seq_len(n_genes))
gene_rate <- rgamma(n_genes, 0.6, 1) + 0.2
depth <- sample(300:4000, n_cells, replace = TRUE)
myeloid <- runif(n_cells) < 0.30
prog_idx <- sample(n_genes, 20)
gene_rate[prog_idx] <- 0.3                                    # low expression -> detection depends on depth
X <- matrix(rpois(n_cells * n_genes, rep(gene_rate, each = n_cells) * (depth / 1000)), n_cells, n_genes)
X[matrix(runif(n_cells * n_genes) < pmax(0, 1 - depth / 4000) * 0.5, n_cells, n_genes)] <- 0L  # 10x dropout
for (j in prog_idx) X[myeloid, j] <- X[myeloid, j] + rpois(sum(myeloid), 3)   # embed programme
colnames(X) <- genes

## (1) UCell programme score
score <- omni_ucell(X, list(prog = genes[prog_idx]))[, "prog"]
lib <- rowSums(X)
cat("Spearman(score, log counts) =", round(omni_score_diagnostics(score, lib)$spearman_log_lib, 3), "\n")

## (2) depth control + the non-parametric alternative
dc <- omni_depth_control(score, log1p(lib))
ds <- omni_depth_stratified_auc(score, log1p(lib), myeloid, nq = 5)
cat("after depth control: r_after =", round(attr(dc, "r_after"), 3), "\n")
cat("myeloid AUC  raw:", round(ds$auc_raw, 3), "| depth-stratified:", round(ds$auc_matched, 3), "\n")

## (3) which cells carry the programme?
tab <- omni_carrier(dc$resid, ifelse(myeloid, "Myeloid", "Other"))
write.csv(tab, "docs/tutorial_outputs/tutorial_1_carriers.csv", row.names = FALSE)
write.csv(data.frame(
  metric = c("spearman_raw", "spearman_after_depth_control", "auc_myeloid_raw", "auc_myeloid_stratified"),
  value  = round(c(omni_score_diagnostics(score, lib)$spearman_log_lib, attr(dc, "r_after"), ds$auc_raw, ds$auc_matched), 3)),
  "docs/tutorial_outputs/tutorial_1_diagnostics.csv", row.names = FALSE)

## figure
col <- ifelse(myeloid, "#d62728", "#9aa7b3")
omni_save("docs/figures/tutorial_1", function() {
  omni_theme(); par(mfrow = c(1, 3))
  plot(log1p(lib), score, col = col, pch = 16, cex = .4, xlab = "log(1+counts)", ylab = "UCell score", main = "(a) raw score vs depth")
  plot(log1p(lib), dc$resid, col = col, pch = 16, cex = .4, xlab = "log(1+counts)", ylab = "residual", main = "(b) depth-controlled")
  a <- setNames(tab$auc_vs_rest, tab$cell_type)
  barplot(c(Myeloid = a[["Myeloid"]], Other = a[["Other"]]), col = c("#d62728", "#9aa7b3"),
          ylab = "AUC vs rest", main = "(c) carrier AUC"); abline(h = .5, lty = 3)
}, width = 13, height = 4)
cat("wrote docs/figures/tutorial_1.{pdf,jpeg} + docs/tutorial_outputs/tutorial_1_*.csv\n")
