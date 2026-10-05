# omnikit

**Algorithm-first toolkit for kidney-disease multi-omics** — depth-robust gene-programme scoring,
cross-species conservation, spatial-niche statistics, and cross-omics integration.

![R-CMD-check](https://github.com/LingzhangMeng/omnikit/actions/workflows/R-CMD-check.yml/badge.svg)
![License: MIT](https://img.shields.io/badge/license-MIT-blue)
![R](https://img.shields.io/badge/R-%3E%3D4.1-blue)

`omnikit` is a small, **dependency-light** R package (base R + `Matrix`) in which every estimator is
implemented **from its definition** — nothing is a thin wrapper, so the behaviour is auditable and
reproducible offline. It is the **R canonical** implementation of the `omnikit` toolkit; a Python
mirror with the same mathematics lives at [**pyomnikit**](https://github.com/LingzhangMeng/pyomnikit).

---

## Why this exists — purpose

Small, heterogeneous disease-omics cohorts (e.g. lupus-nephritis kidney biopsies) keep failing the
**same four ways**. `omnikit` packages the correct, tested solutions so they are applied
consistently instead of being re-implemented (and re-broken) in every project.

| # | Problem (failure mode) | What `omnikit` does |
|---|---|---|
| **G4** | Programme scores are confounded by **sequencing depth** (a rank-based signature score tracks how many signature genes are detected — observed **ρ = +0.72** in real 10x data). | `omni_score_diagnostics`, `omni_depth_control` (linear), `omni_depth_stratified_auc` (non-parametric). |
| **G4** | No standard **cross-species conservation** test; case-fold symbol matching drops alias/many-to-many orthologs. | `omni_conservation` (Fisher + sign consistency + thresholds) + `omni_collapse_orthologs`. |
| **G3** | **QC is skipped / doublet thresholds are unstable**; batch mixing is rarely quantified. | `omni_qc_metrics`, `omni_doublet_scrublet` (calibrated), `omni_mixing` (LISI). |
| **G2** | **Spatial coordinates are ignored** because processed objects often ship none. | `omni_knn_graph`, `omni_moran`, `omni_geary`, `omni_lisa`, `omni_nhood_enrichment`, `omni_ripley`, `omni_spatial_lr`. |
| **G1** | **Omics layers are analysed in isolation.** | `omni_joint_nmf`, `omni_cca`, `omni_rv`, `omni_procrustes`, `omni_mofa_lite`, `omni_program_gwas`. |
| — | Fragile bespoke statistics (AUC CIs, meta-analysis, FDR). | `omni_delong`, `omni_stouffer`, `omni_eb_shrink`, `omni_fisher`, `omni_bh`, `omni_by`, `omni_perm_p`, `omni_auc`, `omni_mwu`. |

**Design principles:** no hidden depth confound · name the null (permutation / matched / analytic) ·
deterministic (`seed=`) · auditable (implemented from the definitions).

## Core algorithms

Every estimator is defined in **[`docs/MATH.md`](docs/MATH.md)** and implemented from scratch. Highlights:

- **UCell** rank-based score with the correct bound `maxU = k·N` (`N = min(maxRank, G)`) — *not* the
  distinct-rank sum, which yields out-of-range scores when ranks are capped.
- **Depth control:** residual of a per-group regression of score on `log(1+counts)`; plus the
  non-parametric **depth-stratified AUC** for non-linear cases.
- **Conservation:** Fisher combine of human/mouse p-values (`χ²₄`) + BH, called with sign
  consistency and `|t| ≥ t₀`, `|ρ| ≥ ρ₀`, `q < α` (defaults `2`, `0.27`, `0.05`).
- **DeLong** AUC variance via O(n log n) midrank placement values.
- **Scrublet** doublet score (simulated doublets → PCA → kNN), calibrated at the expected rate.
- **LISI** mixing; **Moran's I / Geary's C / LISA / Ripley's L** / neighbourhood enrichment.
- **Joint NMF / CCA / RV / Procrustes / MOFA-lite / program–GWAS**.

## Installation

```r
# from GitHub (recommended)
install.packages("remotes")
remotes::install_github("LingzhangMeng/omnikit")
```

```bash
# or from a clone
git clone https://github.com/LingzhangMeng/omnikit.git
R CMD INSTALL omnikit
```

Requirements: **R ≥ 4.1**; imports `stats`, `Matrix` (base/recommended only).

```r
library(omnikit)
?omnikit              # package overview; per-function help e.g. ?omni_ucell
```

## Quickstart

```r
library(omnikit)

# -- G4: score a programme, check + remove depth confounding, find carriers -------------
X         <- as.matrix(counts)                      # cells x genes (raw counts)
genes     <- colnames(X)
score     <- omni_ucell(X, list(prog = program_genes), ...)   # rank-based UCell score
o <- omni_score_diagnostics(score[, "prog"], rowSums(X))      # Spearman(score, log counts)
resid <- omni_depth_control(score[, "prog"], log1p(rowSums(X)))$resid
omni_carrier(resid, celltype)                                 # AUC vs rest per cell type

# -- G4: cross-species conservation ----------------------------------------------------
call <- omni_conservation(df, t_col = "human_meta_t", rho_col = "mouse_rho",
                          ph_col = "human_p", pm_col = "mouse_p", sign_col = "sign_consistent")

# -- G3: QC + doublets + mixing --------------------------------------------------------
omni_qc_metrics(X, genes = genes)
omni_doublet_scrublet(X, k = 30, n_pcs = 30)
omni_mixing(pca_embedding, batch = lane, k = 30)

# -- G2: spatial -----------------------------------------------------------------------
W  <- omni_knn_graph(coords, k = 6)
omni_moran(value, W, n_perm = 999)
omni_lisa(value, W)
omni_ripley(coords)

# -- G1: cross-omics integration -------------------------------------------------------
omni_rv(view1, view2); omni_cca(view1, view3, k = 3)
omni_joint_nmf(list(abs(v1), abs(v2), abs(v3)), k = 2)
omni_mofa_lite(list(v1, v2, v3), k = 2)
```

## Modules

| module | functions |
|---|---|
| `stats` | `omni_auc`, `omni_mwu`, `omni_delong`, `omni_stouffer`, `omni_eb_shrink`, `omni_fisher`, `omni_cohen_d`, `omni_spearman`, `omni_perm_p`, `omni_bh`, `omni_by` |
| `depth` | `omni_depth_control`, `omni_depth_stratified_auc`, `omni_casefold_map`, `omni_collapse_orthologs` |
| `program` | `omni_ucell`, `omni_ucell_signed`, `omni_score_genes`, `omni_conservation`, `omni_carrier`, `omni_score_diagnostics` |
| `qc` | `omni_qc_metrics`, `omni_mad_outlier`, `omni_gmm2`, `omni_doublet_scrublet`, `omni_mixing` |
| `spatial` | `omni_knn_graph`, `omni_moran`, `omni_geary`, `omni_lisa`, `omni_nhood_enrichment`, `omni_ripley`, `omni_spatial_lr` |
| `integrate` | `omni_joint_nmf`, `omni_cca`, `omni_rv`, `omni_procrustes`, `omni_mofa_lite`, `omni_program_gwas` |
| `utils` | `omni_save` (PDF + JPEG, no gridlines), `omni_theme` |

## Validation

```bash
Rscript tests/test_smoke.R        # 27 checks
```

- Reproduces a published **760-gene** cross-species conserved programme with **Jaccard = 1.0**.
- On a raw 10x kidney lane, `omni_depth_control` removes the depth confound (**ρ = +0.72 → +0.09**)
  while preserving the myeloid carriers.
- Deterministic: every randomised estimator takes `seed=`.

## Python mirror

The Python package **pyomnikit** (`pip install pyomnikit`, `import omnikit`) mirrors this package
function-for-function (`omni_<name>` in R vs `omnikit.<name>` in Python). See
<https://github.com/LingzhangMeng/pyomnikit>.

## Citation

See `inst/CITATION` (R `citation("omnikit")`) and cite the underlying methods (UCell, Scrublet,
DeLong, Moran, Lee–Seung NMF, MOFA, scDRS) listed in [`docs/MATH.md`](docs/MATH.md).

## License

MIT © 2026 Lingzhang Meng — see [`LICENSE.md`](LICENSE.md).
