#!/usr/bin/env Rscript
# examples/tutorial_6_topo.R
# Tutorial 6 — topological conservation: 0-dim persistent homology (H0) + bottleneck distance.
# Self-contained (synthetic). Run from the package root:  Rscript examples/tutorial_6_topo.R
suppressMessages(library(omnikit))
dir.create("docs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs/tutorial_outputs", recursive = TRUE, showWarnings = FALSE)
set.seed(6)

## two related point clouds (shared latent manifold + independent noise)
n <- 200
latent <- matrix(rnorm(n * 3), n, 3)
A <- latent + matrix(rnorm(n * 3, 0, 0.2), n, 3)
B <- latent + matrix(rnorm(n * 3, 0, 0.2), n, 3)
dA <- omni_h0_diagram(A); dB <- omni_h0_diagram(B)
bd <- omni_bottleneck_h0(dA, dB)

## null: unrelated random clouds -> bottleneck distribution
null <- vapply(1:100, function(i) omni_bottleneck_h0(dA, omni_h0_diagram(matrix(rnorm(n * 3), n, 3))), numeric(1))
p <- (1 + sum(null <= bd)) / (1 + length(null))

write.csv(data.frame(n_points = n, h0_bars_A = length(dA), h0_bars_B = length(dB),
                     entropy_A = omni_persistent_entropy(dA), entropy_B = omni_persistent_entropy(dB),
                     bottleneck_AB = bd, null_mean = mean(null), p_emp = p),
          "docs/tutorial_outputs/tutorial_6_topo.csv", row.names = FALSE)
cat(sprintf("H0 bars: %d/%d | bottleneck(A,B) = %.3f | null mean = %.3f | p = %.3f\n",
            length(dA), length(dB), bd, mean(null), p))

## figure
omni_save("docs/figures/tutorial_6", function() {
  omni_theme(); par(mfrow = c(1, 3))
  plot(A[, 1], A[, 2], pch = 16, cex = .5, col = "#1f77b4", xlab = "x", ylab = "y",
       main = "(a) two related point clouds")
  points(B[, 1], B[, 2], pch = 16, cex = .5, col = "#d62728")
  plot(rep(0, length(dA)), dA, pch = 16, cex = .6, col = "#1f77b4", ylim = range(c(dA, dB)),
       xlab = "birth", ylab = "death", main = "(b) H0 persistence diagrams")
  points(rep(0, length(dB)), dB, pch = 16, cex = .6, col = "#d62728")
  hist(null, breaks = 20, col = "#bbbbbb", main = sprintf("(c) bottleneck null (p=%.3f)", p),
       xlab = "bottleneck distance")
  abline(v = bd, col = "#c0392b", lwd = 2)
}, width = 14, height = 4.3)
cat("wrote docs/figures/tutorial_6.{pdf,jpeg} + docs/tutorial_outputs/tutorial_6_topo.csv\n")
