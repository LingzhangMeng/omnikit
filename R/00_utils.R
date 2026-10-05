# omnikit — 00_utils.R
# Small helpers shared by every module. Nothing here that begins with "." is exported.

#' Save a plot as BOTH vector PDF and 300-dpi JPEG, no grid lines.
#' @param base Path without extension; the PDF and JPEG are written.
#' @param plot_fun A zero-argument function that draws the plot.
#' @param width,height Figure size in inches.
#' @param dpi JPEG resolution.
#' @return (Invisibly) the two written file paths.
omni_save <- function(base, plot_fun, width = 6, height = 4.5, dpi = 300) {
  grDevices::pdf(paste0(base, ".pdf"), width = width, height = height)
  plot_fun(); grDevices::dev.off()
  grDevices::jpeg(paste0(base, ".jpeg"), width = width, height = height,
                  units = "in", res = dpi, bg = "white")
  plot_fun(); grDevices::dev.off()
  invisible(c(pdf = paste0(base, ".pdf"), jpeg = paste0(base, ".jpeg")))
}

#' A minimal grid-free base theme for omnikit figures.
#' @param mar Plot margins (see \code{graphics::par}).
#' @param cex Character expansion.
#' @return (Invisibly) \code{NULL}; called for its side effect on the device.
omni_theme <- function(mar = c(4, 4, 2, 1), cex = 1) {
  graphics::par(mar = mar, cex = cex, las = 1, bty = "n", xaxs = "r", yaxs = "r")
}

#' Coerce to a base numeric matrix (dense) — Matrix objects are expanded.
#' @param X A matrix, data.frame, or Matrix object.
#' @return A base \code{matrix}.
omni_as_matrix <- function(X) {
  if (inherits(X, "Matrix")) return(as.matrix(X))
  if (is.data.frame(X)) return(as.matrix(X))
  X
}

#' Row / column names with a safe fallback.
#' @param X A matrix-like object.
#' @return A character vector of names (or indices if unnamed).
#' @rdname omni_names
omni_rnames <- function(X) if (!is.null(rownames(X))) rownames(X) else as.character(seq_len(nrow(X)))
#' @rdname omni_names
omni_cnames <- function(X) if (!is.null(colnames(X))) colnames(X) else as.character(seq_len(ncol(X)))

#' BH / BY FDR on a numeric p vector (NA-safe).
#' @param p Numeric vector of p-values.
#' @return Adjusted p-values.
#' @rdname omni_fdr
omni_bh <- function(p) stats::p.adjust(p, method = "BH")
#' @rdname omni_fdr
omni_by <- function(p) stats::p.adjust(p, method = "BY")
