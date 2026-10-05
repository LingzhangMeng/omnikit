# omnikit — smoke tests (run: Rscript tests/test_smoke.R)
suppressMessages(library(omnikit))
ok <- function(label, cond) cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", label))
set.seed(1)

## ---- G4: UCell bound (top-ranked ~1, bottom-ranked ~0) ----
# deterministic gradient: gene j has the (61-j)-th highest expression in every cell
X <- matrix(rep(61 - (1:60), each = 300), 300, 60, dimnames = list(NULL, paste0("g", 1:60)))
X <- X + matrix(stats::runif(300 * 60, 0, 0.1), 300, 60)      # break ties
uc <- omni_ucell(X, list(top = paste0("g", 1:5), bottom = paste0("g", 56:60)), chunk = 100)
ok("UCell top signature ~ 1",  all(uc[, "top"] > 0.99))
ok("UCell bottom signature ~ 0", all(uc[, "bottom"] < 0.1))
ok("UCell in [0,1]", all(uc >= 0 & uc <= 1))

## ---- stats: AUC / DeLong / Stouffer / EB / Fisher ----
s <- c(stats::rnorm(200, 1), stats::rnorm(200, 0))
p <- c(rep(TRUE, 200), rep(FALSE, 200))
d <- omni_delong(s, p)
ok("AUC in (0,1) and CI brackets it", d$auc > 0.5 && d$lower < d$auc && d$auc < d$upper)
ok("AUC == omni_auc", abs(d$auc - omni_auc(s, p)) < 1e-9)
st <- omni_stouffer(c(2, 2, 2)); ok("Stouffer Z of 3x2 ~ 3.46", abs(st$Z - 2 * sqrt(3)) < 1e-6)
eb <- omni_eb_shrink(c(1, 2, 3, 4), rep(1, 4)); ok("EB shrink finite", all(is.finite(eb)))
fi <- omni_fisher(c(0.01, 0.02, 0.03)); ok("Fisher p < min(p)", fi$p < 0.01)

## ---- depth control recovers a clean score ----
lib <- exp(stats::rnorm(500, 8, 1)); truth <- stats::rnorm(500)
sc <- truth + 0.4 * log(lib)
dc <- omni_depth_control(sc, log(lib))
ok("depth control reduces |rho|", abs(attr(dc, "r_after")) < abs(attr(dc, "r_before")))
ds <- omni_depth_stratified_auc(sc, log(lib), truth > 0, nq = 5)
ok("depth-stratified AUC returns bins", nrow(ds$per_bin) >= 3)

## ---- G3: QC metrics / MAD / Scrublet / mixing ----
Xc <- matrix(stats::rpois(400 * 50, 3), 400, 50, dimnames = list(NULL, paste0("g", 1:50)))
colnames(Xc)[1:3] <- c("MT-CO1", "RPS3", "HBB")
qc <- omni_qc_metrics(Xc)
ok("QC metrics columns", all(c("nCount", "nFeature", "pct_mt") %in% names(qc)))
ok("MAD flags outliers", any(omni_mad_outlier(c(rep(0, 100), 100))$flag))
Xc[1:5, ] <- Xc[1:5, ] + Xc[6:10, ]           # make 5 artificial doublets
dbl <- omni_doublet_scrublet(Xc, sim_doublets = 200, k = 15, n_pcs = 15)
ok("Scrublet returns scores in [0,1]", all(dbl$score >= 0 & dbl$score <= 1))
mx <- omni_mixing(matrix(stats::rnorm(300 * 5), 300, 5), rep(c("a", "b", "c"), 100), k = 15)
ok("LISI in [1,B]", all(mx$lisi >= 1 - 1e-9 & mx$lisi <= mx$max + 1e-9))

## ---- G2: spatial ----
coords <- cbind(x = stats::runif(400), y = stats::runif(400))
W <- omni_knn_graph(coords, k = 6)
fld <- coords[, 1] + 0.05 * stats::rnorm(400)  # smooth spatial field
mo <- omni_moran(fld, W); ok("Moran I > 0 for smooth field", mo$I > 0)
ge <- omni_geary(fld, W); ok("Geary C < 1 for smooth field", ge$C < 1)
li <- omni_lisa(fld, W, n_perm = 49); ok("LISA returns Ii", nrow(li) == 400)
lab <- ifelse(coords[, 1] < 0.5, "A", "B")
ne <- omni_nhood_enrichment(lab, W, n_perm = 49)
ok("nhood z matrix square", all(dim(ne$zscore) == c(2, 2)))
rp <- omni_ripley(coords); ok("Ripley K monotone-ish", all(is.finite(rp$K)))
E <- cbind(L = stats::rpois(400, 2), R = stats::rpois(400, 2), G = stats::rpois(400, 2))
slr <- omni_spatial_lr(E, list(c("L", "R")), W, n_perm = 49)
ok("spatial LR returns p", is.finite(slr$p_perm[1]))

## ---- G1: integration ----
n <- 200; z <- matrix(stats::rnorm(n * 2), n, 2)
v1 <- z %*% matrix(stats::rnorm(2 * 30), 2, 30) + matrix(stats::rnorm(n * 30), n, 30)
v2 <- z %*% matrix(stats::rnorm(2 * 20), 2, 20) + matrix(stats::rnorm(n * 20), n, 20)
nmf <- omni_joint_nmf(list(abs(v1), abs(v2)), k = 2, n_iter = 60)
ok("joint NMF shapes", nrow(nmf$W) == n && ncol(nmf$W) == 2 && length(nmf$H) == 2)
cc <- omni_cca(v1, v2, k = 3); ok("CCA cor descending", all(diff(cc$cor) <= 1e-9) && length(cc$cor) == 3)
ok("RV in (0,1)", omni_rv(v1, v2) > 0 && omni_rv(v1, v2) < 1)
pr <- omni_procrustes(v1, v1 + matrix(stats::rnorm(n * 30, 0, 0.01), n, 30))
ok("Procrustes residual small", pr$residual < 0.1)
mf <- omni_mofa_lite(list(v1, v2), k = 2, n_iter = 60)
ok("MOFA factors shape", all(dim(mf$factors) == c(n, 2)) && length(mf$var_explained) == 2)
pg <- omni_program_gwas(stats::rnorm(50, 1), stats::rnorm(200, 0))
ok("program-GWAS z > 0", pg$z > 0)

cat("\nomnikit smoke test complete.\n")
