#!/usr/bin/env Rscript
# examples/tutorial_3_qc_doublets.R
# Tutorial 3 — QC-first audit: per-cell metrics, Scrublet-style doublets, batch mixing (LISI).
# Self-contained (synthetic). Run from the package root:  Rscript examples/tutorial_3_qc_doublets.R
suppressMessages(library(omnikit))
dir.create("docs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs/tutorial_outputs", recursive = TRUE, showWarnings = FALSE)
set.seed(2)

## synthetic counts: 2 batches, mitochondrial genes, 40 injected doublets
n <- 2000; g <- 300
genes <- sprintf("G%03d", seq_len(g)); genes[1:2] <- c("MT-CO1", "MT-ND1")
X <- matrix(rpois(n * g, 2), n, g); X[, 1:2] <- X[, 1:2] + rpois(n * 2, 5)
for (i in 1:40) { ab <- sample(n, 2); X[i, ] <- X[ab[1], ] + X[ab[2], ] }
colnames(X) <- genes
batch <- rep(c("batchA", "batchB"), each = n / 2)

## (1) QC metrics
qc <- omni_qc_metrics(X)
write.csv(cbind(qc, batch = batch), "docs/tutorial_outputs/tutorial_3_qc_metrics.csv", row.names = FALSE)
cat("median %mt:", round(median(qc$pct_mt), 3), "\n")

## (2) Scrublet-style doublets  (integer counts)
d <- omni_doublet_scrublet(X, k = 30, n_pcs = 20, seed = 0)
write.csv(data.frame(score = d$score, call = d$call), "docs/tutorial_outputs/tutorial_3_doublets.csv", row.names = FALSE)
cat(sprintf("doublets called: %d/%d (%.1f%%)  threshold=%s\n", sum(d$call), n, 100 * mean(d$call),
            if (is.null(d$threshold)) "NA" else sprintf("%.3f", d$threshold)))

## (3) batch mixing (LISI) on a 10-PC embedding
Xn  <- log1p(sweep(X, 1, pmax(rowSums(X), 1), "/") * 1e4)
pcs <- svd(scale(Xn), nu = 10, nv = 0)$u
mx  <- omni_mixing(pcs, batch, k = 30)
cat(sprintf("LISI mean = %.2f (max %d)\n", mx$mean, mx$max))

## figure
omni_save("docs/figures/tutorial_3", function() {
  omni_theme(); par(mfrow = c(1, 3))
  plot(qc$nCount, qc$pct_mt, col = ifelse(d$call, "#d62728", "#9aa7b3"), pch = 16, cex = .5,
       xlab = "nCount", ylab = "% mitochondrial", main = "(a) QC (red = doublet)")
  hist(d$score, breaks = 50, col = "#6f8fb0", xlab = "doublet score", main = "(b) doublet score")
  abline(v = d$threshold, col = "red", lty = 2)
  hist(mx$lisi, breaks = 40, col = "#8fae6f", xlab = "LISI", main = sprintf("(c) mixing (mean %.2f)", mx$mean))
  abline(v = c(1, mx$max), lty = 3)
}, width = 13, height = 4)
cat("wrote docs/figures/tutorial_3.{pdf,jpeg} + docs/tutorial_outputs/tutorial_3_*.csv\n")
