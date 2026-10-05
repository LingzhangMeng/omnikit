#!/usr/bin/env Rscript
# examples/tutorial_4_spatial.R
# Tutorial 4 — spatial niche statistics: Moran's I, LISA, neighbourhood graph, Ripley's L.
# Self-contained (synthetic). Run from the package root:  Rscript examples/tutorial_4_spatial.R
suppressMessages(library(omnikit))
dir.create("docs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs/tutorial_outputs", recursive = TRUE, showWarnings = FALSE)
set.seed(3)

## synthetic 2-D tissue: 1500 cells, a smooth spatial gradient + noise
n <- 1500
coords <- matrix(runif(n * 2) * 100, n, 2)
value  <- sin(coords[, 1] / 15) * cos(coords[, 2] / 15) + rnorm(n, 0, .3)

## spatial graph + Moran / Geary / LISA
W  <- omni_knn_graph(coords, k = 6)
mo <- omni_moran(value, W, n_perm = 999, seed = 0)
ge <- omni_geary(value, W, n_perm = 999)
li <- omni_lisa(value, W, n_perm = 99)
rp <- omni_ripley(coords)
cat(sprintf("Moran I = %.3f (p %.3f) | Geary C = %.3f\n", mo$I, mo$p_perm, ge$C))
write.csv(data.frame(radius = rp$radius, K = rp$K, L = rp$L), "docs/tutorial_outputs/tutorial_4_ripley.csv", row.names = FALSE)

## neighbourhood enrichment (two synthetic cell types)
labels <- ifelse(runif(n) < .5, "A", "B")
ne <- omni_nhood_enrichment(labels, W, n_perm = 199, seed = 0)
write.csv(ne$zscore, "docs/tutorial_outputs/tutorial_4_nhood.csv")

## figure
rng <- colorRampPalette(c("#4575b4", "white", "#d73027"))(100)
omni_save("docs/figures/tutorial_4", function() {
  omni_theme(); par(mfrow = c(1, 3))
  plot(coords, col = rng[cut(value, 100)], pch = 16, cex = .5, asp = 1, main = sprintf("(a) value (Moran I=%.2f)", mo$I))
  plot(coords, col = rng[cut(li$Ii, 100)], pch = 16, cex = .5, asp = 1, main = "(b) local Moran (LISA)")
  plot(rp$radius, rp$L, type = "l", col = "#1f77b4", xlab = "radius", ylab = "L(r) - r", main = "(c) Ripley's L")
  abline(h = 0, lty = 3)
}, width = 14, height = 4.4)
cat("wrote docs/figures/tutorial_4.{pdf,jpeg} + docs/tutorial_outputs/tutorial_4_*.csv\n")
