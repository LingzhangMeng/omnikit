#!/usr/bin/env Rscript
# examples/tutorial_5_integration.R
# Tutorial 5 — cross-omics integration: RV, CCA, joint NMF and MOFA-lite over shared samples.
# Self-contained (synthetic). Run from the package root:  Rscript examples/tutorial_5_integration.R
suppressMessages(library(omnikit))
dir.create("docs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs/tutorial_outputs", recursive = TRUE, showWarnings = FALSE)
set.seed(4)

## three "omics" views over the same 3000 samples, driven by 2 shared factors
n  <- 3000
z  <- matrix(rnorm(n * 2), n, 2)
v1 <- z %*% matrix(rnorm(2 * 8), 2, 8) + matrix(rnorm(n * 8, 0, .8), n, 8)   # transcriptome (8 features)
v2 <- z %*% matrix(rnorm(2 * 4), 2, 4) + matrix(rnorm(n * 4, 0, .8), n, 4)   # mouse transcriptome (4)
v3 <- z %*% matrix(rnorm(2 * 1), 2, 1) + matrix(rnorm(n, 0, 1), n, 1)        # genetics (1)

## RV coefficient + CCA
rv <- c(human_mouse = omni_rv(v1, v2), human_genetics = omni_rv(v1, v3), mouse_genetics = omni_rv(v2, v3))
cc <- omni_cca(v1, v3, k = 3)
write.csv(data.frame(pair = names(rv), RV = round(rv, 3)), "docs/tutorial_outputs/tutorial_5_rv.csv", row.names = FALSE)

## joint multi-view NMF + MOFA-lite
nmf <- omni_joint_nmf(list(abs(v1), abs(v2), abs(v3)), k = 2, n_iter = 300, seed = 0)
mf  <- omni_mofa_lite(list(v1, v2, v3), k = 2, n_iter = 200, seed = 0)
write.csv(data.frame(view = c("human", "mouse", "genetics"), var_explained = mf$var_explained),
          "docs/tutorial_outputs/tutorial_5_mofa.csv", row.names = FALSE)
cat("RV:", paste(names(rv), round(rv, 3), sep = "=", collapse = " "), "\n")
cat("CCA top:", round(cc$cor[1], 3), "| MOFA VE:", paste(round(mf$var_explained, 3), collapse = ", "), "\n")

## figure
omni_save("docs/figures/tutorial_5", function() {
  omni_theme(); par(mfrow = c(1, 3))
  barplot(rv, col = c("#4c72b0", "#dd8452", "#55a868"), ylab = "RV coefficient", main = "(a) RV similarity")
  xv <- v1 %*% cc$x_weights[, 1]; yv <- v3 %*% cc$y_weights[, 1]
  plot(xv, yv, col = "#9aa7b3", pch = 16, cex = .4, xlab = "human canonical variate 1",
       ylab = "genetics canonical variate 1", main = sprintf("(b) CCA (r = %.2f)", cc$cor[1]))
  barplot(mf$var_explained, names.arg = c("human", "mouse", "genetics"), col = "#6f8fb0",
          ylab = "variance explained", main = "(c) MOFA-lite per view")
}, width = 13.5, height = 4.2)
cat("wrote docs/figures/tutorial_5.{pdf,jpeg} + docs/tutorial_outputs/tutorial_5_*.csv\n")
