#!/usr/bin/env Rscript
# examples/tutorial_7_conformal.R
# Tutorial 7 — split & Mondrian conformal prediction: finite-sample coverage and set size.
# Self-contained (synthetic). Run from the package root:  Rscript examples/tutorial_7_conformal.R
suppressMessages(library(omnikit))
dir.create("docs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs/tutorial_outputs", recursive = TRUE, showWarnings = FALSE)
set.seed(7)

## synthetic 3-class probabilistic classifier: separated class centres
K <- 3; n <- 1500
centres <- rbind(c(2, 0), c(0, 2), c(-2, -1))
make <- function(n) {
  y <- sample(seq_len(K), n, replace = TRUE)
  z <- centres[y, , drop = FALSE] + matrix(rnorm(n * 2, 0, 0.7), n, 2)
  d <- -sapply(seq_len(K), function(c) rowSums((z - matrix(centres[c, ], n, 2, byrow = TRUE))^2)) / 2
  p <- exp(d - apply(d, 1, max)); p <- p / rowSums(p)
  list(y = y, p = p)
}
cal <- make(n); te <- make(n)
alpha <- 0.10

q <- omni_conformal_quantile(omni_nonconformity(cal$p, cal$y), alpha)
sets_split <- omni_prediction_sets(te$p, q)
qm <- omni_mondrian_quantiles(cal$p, cal$y, alpha)
sets_mondrian <- omni_prediction_sets(te$p, qm)

cov_split <- omni_coverage(sets_split, te$y); cov_mond <- omni_coverage(sets_mondrian, te$y)
write.csv(data.frame(method = c("split", "Mondrian"), target = 1 - alpha,
                     coverage = round(c(cov_split, cov_mond), 3),
                     avg_set_size = round(c(omni_avg_set_size(sets_split), omni_avg_set_size(sets_mondrian)), 2)),
          "docs/tutorial_outputs/tutorial_7_coverage.csv", row.names = FALSE)
cat(sprintf("split coverage = %.3f | Mondrian coverage = %.3f (target %.2f)\n", cov_split, cov_mond, 1 - alpha))

## figure
omni_save("docs/figures/tutorial_7", function() {
  omni_theme(); par(mfrow = c(1, 3))
  boxplot(cal$p[cbind(seq_len(n), cal$y)] ~ factor(cal$y), col = c("#1f77b4", "#d62728", "#2ca02c"),
          xlab = "true class", ylab = "p(true class)", main = "(a) calibration probabilities")
  barplot(table(factor(rowSums(sets_split), levels = 0:K)), col = "#7f9fc4",
          xlab = "set size", ylab = "count", main = "(b) split-conformal set size")
  bp <- barplot(c(cov_split, cov_mond), names.arg = c("split", "Mondrian"),
                col = c("#7f9fc4", "#c0392b"), ylim = c(0, 1.05), ylab = "empirical coverage",
                main = "(c) coverage vs target")
  abline(h = 1 - alpha, lty = 3, lwd = 2)
}, width = 14, height = 4.3)
cat("wrote docs/figures/tutorial_7.{pdf,jpeg} + docs/tutorial_outputs/tutorial_7_coverage.csv\n")
