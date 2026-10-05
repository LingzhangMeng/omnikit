#!/usr/bin/env Rscript
# examples/tutorial_2_conservation.R
# Tutorial 2 — call a cross-species conserved gene programme from per-gene statistics.
# Self-contained (synthetic). Run from the package root:  Rscript examples/tutorial_2_conservation.R
suppressMessages(library(omnikit))
dir.create("docs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs/tutorial_outputs", recursive = TRUE, showWarnings = FALSE)
set.seed(1)

## synthetic per-gene statistics over 2000 genes
n <- 2000; genes <- sprintf("G%04d", seq_len(n))
true_prog <- runif(n) < 0.15
human_t   <- rnorm(n); mouse_rho <- rnorm(n, sd = .15)
human_t[true_prog]   <- human_t[true_prog]   + rnorm(sum(true_prog), 3, 1)
mouse_rho[true_prog] <- mouse_rho[true_prog] + rnorm(sum(true_prog), .5, .1)
human_p <- 2 * pnorm(-abs(human_t))
mouse_p <- 2 * pt(-abs(mouse_rho * sqrt(50) / sqrt(1 - mouse_rho^2)), 50)
df <- data.frame(gene = genes, human_t = human_t, mouse_rho = mouse_rho, human_p = human_p, mouse_p = mouse_p)

## conservation call (Fisher combine + sign consistency + thresholds)
call <- omni_conservation(df, t_col = "human_t", rho_col = "mouse_rho", ph_col = "human_p", pm_col = "mouse_p")
write.csv(call, "docs/tutorial_outputs/tutorial_2_conservation.csv", row.names = FALSE)
prec <- mean(call$gene[call$conserved] %in% genes[true_prog])
cat(sprintf("conserved genes: %d | precision vs truth: %.2f\n", sum(call$conserved), prec))

## figure: human_t vs mouse_rho, conserved genes in red
omni_save("docs/figures/tutorial_2", function() {
  omni_theme()
  plot(call$human_t, call$mouse_rho, col = "#cccccc", pch = 16, cex = .4,
       xlab = "human effect (t)", ylab = "mouse trend (rho)",
       main = sprintf("Tutorial 2 - conserved programme (%d genes)", sum(call$conserved)))
  points(call$human_t[call$conserved], call$mouse_rho[call$conserved], col = "#d62728", pch = 16, cex = .5)
  abline(h = 0, v = 0, col = "grey70")
}, width = 5.2, height = 4.6)
cat("wrote docs/figures/tutorial_2.{pdf,jpeg} + docs/tutorial_outputs/tutorial_2_conservation.csv\n")
