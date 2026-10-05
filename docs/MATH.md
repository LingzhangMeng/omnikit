# omnikit — mathematical specification

Every estimator below is implemented directly from its definition. Notation: a dataset has
`n` cells/samples × `G` genes; `S` is a signature (gene set); `x_i` is the value of the score in
cell `i`; `L_i = log(1 + library size_i)` is the depth covariate.

---

## 0. Shared statistics (`stats`)

### 0.1 AUC as a normalised Mann–Whitney U
For a score `x` and a positive mask `m` (|m| = n₁, ¬m = n₀):

    U  = R⁺ − n₁(n₁+1)/2 ,   R⁺ = Σ_{i∈m} rank(x_i)
    AUC = U / (n₁ n₀)

AUC = P(x_pos > x_neg) + ½P(x_pos = x_neg). Under H₀, AUC = ½.

### 0.2 DeLong variance & CI
With placement values `V₁₀(i) = (1/n₀)Σ_{j∉m} ψ(x_i, x_j)`, `V₀₁(j) = (1/n₁)Σ_{i∈m} ψ(x_i, x_j)`,
and `ψ(a,b) = 1[a>b] + ½·1[a=b]`, the DeLong covariance is

    S₁₀ = (1/(n₁−1)) Σ (V₁₀(i) − AUC)² ,   S₀₁ likewise,
    Var(AUC) = S₁₀/n₁ + S₀₁/n₀ .

CI = AUC ± z_{1−α/2}·√Var (Wald; logit-scale option in the code). The implementation uses the
O((n₁+n₀) log(n₁+n₀)) midrank formulation (Sun & Xu 2014), not the O(n²) double loop.

### 0.3 Weighted Stouffer meta
    Z = Σ_k w_k z_k / √(Σ_k w_k²) ,   p = 2·Φ(−|Z|)
If a study reports p (not z), `z = −Φ⁻¹(p/2)·sign(effect)`.

### 0.4 Empirical-Bayes (hierarchical) shrinkage
Given per-unit estimates `β_j` with standard errors `se_j`: the method-of-moments prior variance is
`τ² = max(0, Var(β) − mean(se²))`; the posterior mean shrinks toward the grand mean `β̄` with weight

    β̃_j = β̄ + (β_j − β̄) · τ² / (τ² + se_j²) .

(Limited-translation / limma-flavoured; `omni_eb_shrink`.)

### 0.5 Multiple testing
BH: sort p ascending, `q_(k) = min_{j≥k} p_(j)·m/j`. BY: divide by `Σ 1/i` (arbitrary dependence).

### 0.6 Effect sizes
Welch t, Cohen's d (pooled SD), Spearman ρ (Pearson on ranks, tie-averaged), Fisher combine
`X² = −2Σ ln p ~ χ²_{2k}` with the Brown correction (mean/variance standardisation).

---

## 1. Depth control (`stats`)

### 1.1 Global regression control
Fit `x_i = a + b·L_i + ε_i` (ordinary least squares); the depth-controlled score is the residual
`x*_i = x_i − (â + b̂ L_i)`. Diagnostics report `ρ(x, L)` before and `ρ(x*, L)` after.

### 1.2 Depth-stratified AUC (matching, not modelling)
Partition cells into `Q` quantile bins of `L`. Within bin `q` compute `AUC_q`; the matched
estimate is the size-weighted mean `Σ_q (n_q/n) AUC_q`. This is valid even when the `x–L`
relation is non-linear (e.g. UCell's rank-cap effect). Report `min_q AUC_q` as the conservative
floor.

### 1.3 Why raw UCell is depth-sensitive on sparse counts
`omni_ucell` caps ranks at `maxRank`. A signature gene that is **undetected** receives rank
`maxRank` (the cap), which contributes the *maximum* `U`; therefore a cell that detects *more*
signature genes (higher depth) gets a *lower* `U` → *higher* score. Hence `ρ(x, L) > 0` on
zero-inflated 10x data, and `1.2` is the correct correction. See §2.1.

---

## 2. Program scoring (G4, `program`)

### 2.1 UCell (rank-based, depth-robust *in principle*)
Per cell, rank the `G` genes by descending expression (ties broken at random), cap ranks at
`maxRank` (UCell default 1500): `r_i = min(rank_i, maxRank)`. For a signature `S` of size `k`:

    U = Σ_{g∈S} r_g ,   N = min(maxRank, G) ,
    minU = k(k+1)/2 ,   maxU = k·N ,        (★ the correct bound)
    score = (maxU − U)/(maxU − minU) ∈ [0,1] .

(★) `maxU = k·N`, **not** the distinct-rank sum `k(2N−k+1)/2`; the latter yields negative /
out-of-range scores once ranks are capped. A top-ranked signature → ≈1, bottom-ranked → ≈0.

For a **signed** programme, `score = UCell(up) − UCell(down)`.

### 2.2 Control-bin (expression-matched) score
scanpy/Pagès-style: place genes into `n_bins` bins of mean expression; for signature gene `g`,
draw `ctrl_size` control genes from the neighbouring bins *excluding* the signature; the score is
the mean(signature) − mean(control), averaged over bins, optionally z-scaled by the binned SD.
Robust to depth because signature and control genes share the expression distribution.

### 2.3 Conservation call
Per shared ortholog `g`: human statistic `t_g` (per-compartment, combined by Stouffer/Fisher
across compartments), mouse trend `ρ_g` (Spearman of expression vs ordered severity). Combine by
Fisher: `X²_g = −2(ln p^H_g + ln p^M_g)`, `q_g = BH(X² p)`. The **conserved** set requires
(i) sign consistency — both human compartments agree with each other *and* with `sign(ρ_g)` —
and (ii) `|t_g| ≥ t₀`, `|ρ_g| ≥ ρ₀`, `q_g < α`. Defaults `t₀=2`, `ρ₀=0.27`, `α=0.05`
(Paper_1 pre-registration).

### 2.4 Permutation null
Stratified label permutation (preserve compartment/severity design), full pipeline recompute,
`B` draws. Empirical p = `(1 + #{null ≥ obs})/(B + 1)`.

### 2.5 Carrier inference
Per cell type `c`: AUC of the score vs all other cells, Mann–Whitney p, BH q. **Donor-aware**
contrast: pseudo-bulk mean per (donor × cell type); the within-donor difference of type `c` from
the donor's other cells, tested by a paired Wilcoxon / permutation across donors (removes
pseudo-replication without collapsing opposing donors the way naive pooling does).

---

## 3. QC-first audit (G3, `qc`)

### 3.1 Metrics
`nFeature_i = #{g: X_ig>0}`, `nCount_i = Σ_g X_ig`, `%mt/%ribo/%hb` by gene-set membership.
Detected over the **raw/integer** counts.

### 3.2 Doublet score (Scrublet, Wolock 2019), re-implemented
1. Expected doublet rate `π ≈ n·r/1000` (10x loading); simulate `round(π·n/…)` doublets by
   summing expression of randomly paired observed cells: `x_sim = x_a + x_b`.
2. Pool observed + simulated; normalise (library-size), log1p, scale, PCA (via truncated SVD).
3. For each observed cell, kNN (default k=30) in PCA space; `d_i =` fraction of neighbours that
   are simulated.
4. Threshold `θ` = the minimum of a 2-component Gaussian KDE of `d` (the dip between the modes);
   cells with `d_i ≥ θ` are doublets. `omni_doublet_scrublet`.

### 3.3 Robust thresholds
`omni_mad_outlier(x, k)`: flag `|x − median| > k·1.4826·MAD` (MAD ⇒ Gaussian-consistent scale).

### 3.4 Batch mixing
Local Inverse Simpson Index (LISI): for cell `i` over its `k` neighbours,
`LISI_i = 1/Σ_b p_{ib}²` where `p_{ib}` is the neighbour batch-b proportion; `E[LISI]=1` (pure),
`=B` (fully mixed). `omni_mixing` (uses the kNN graph from §5.1 or an embedding).

---

## 4. Spatial niche statistics (G2, `spatial`)

Graph `W` = symmetric kNN (default k=6) with row-normalised weights (or fixed-radius).

### 4.1 Moran's I (global autocorrelation)
    I = (n / S₀) · (zᵀ W z) / (zᵀ z) ,   z = x − x̄ ,  S₀ = Σ_ij W_ij .
`E[I] = −1/(n−1)`. Analytic variance (randomisation) + a permutation null (shuffle `x` over the
graph). `omni_geary`: `C = ((n−1)/2S₀)·Σ W_ij(x_i−x_j)²/Σ(x_i−x̄)²`, `E[C]=1`.

### 4.2 LISA (local Moran)
    I_i = z_i · Σ_j W_ij z_j  (z standardised, `Σ z²/n = 1`);  conditional permutation p per cell.

### 4.3 Neighbourhood enrichment (cell-type co-occurrence)
Count `O_ab` = #(edges with labels a–b) on the kNN graph; null by permuting labels `B×`;
`z_ab = (O_ab − μ_ab)/σ_ab`. `omni_nhood_enrichment`.

### 4.4 Ripley's K / L (clustering of a marked point pattern)
    K(r) = (|A| / n(n−1)) · Σ_{i≠j} 1[d_ij ≤ r] · e_ij ,   L(r) = √(K(r)/π) − r ,
with a border (or Ripley) edge correction `e_ij`. 2D CSRs.

### 4.5 Spatial ligand–receptor co-expression
For each LR pair (L,R): co-expression `s_i = 1[x^L_i>0]·1[x^R_i>0]`; neighbourhood statistic
`Γ = Σ_ij W_ij s_i s_j / Σ_ij W_ij`; permutation null (shuffle positions). `omni_spatial_lr`.

---

## 5. Cross-omics integration (G1, `integrate`)

### 5.1 Multi-view (joint) NMF
Given views `X^{(v)} ≥ 0` (e.g. gene-expression, program scores, GWAS gene-z), factorise
`X^{(v)} ≈ W^{(v)} H` sharing the **common** factor matrix `H` (`k` factors). Multiplicative
updates (Lee–Seung), with `λ`-regularised couplings for partially-shared structure:

    H ← H ∘ ( Σ_v W^{(v)ᵀ} X^{(v)} ) ⊘ ( Σ_v W^{(v)ᵀ} W^{(v)} H + λ (H − H^{(v)}) )

`omni_joint_nmf` (symmetric variant `omni_joint_nmf_shared`).

### 5.2 Canonical correlation (CCA)
Given centred `X∈ℝ^{n×p}`, `Y∈ℝ^{n×q}` with covariances `S_xx, S_yy`, `S_xy`:
`T = S_xx^{−½} S_xy S_yy^{−½}`, `σ_k =` singular values of `T` = canonical correlations.
Weights = whitening vectors ▸ right singular vectors. `omni_cca`.

### 5.3 RV coefficient (global similarity of two configurations)
    RV = tr(S_xy S_yx) / √( tr(S_xx²) · tr(S_yy²) ) ,  S_xy = XᵀY (both centred).
= cosine similarity of the two configuration matrices flattened. `omni_rv`.

### 5.4 Procrustes
Optimal rotation/scale minimising `‖sXR − Y‖_F`: via SVD `XᵀY = UΣVᵀ` ⇒ `R = UVᵀ`,
`s = tr(Σ)/‖X‖_F²`. Residual `= 1 − (tr Σ)²/(‖X‖²‖Y‖²)`. `omni_procrustes`.

### 5.5 MOFA-lite (multi-view factor analysis, EM)
Model `X^{(v)}_i = Z_i Λ^{(v)ᵀ} + ε^{(v)}`, `Z_i ∼ N(0, I_k)`, `ε^{(v)} ∼ N(0, σ_v² I)`.
E-step: `E[Z|X] = ( Σ_v Λ^{(v)ᵀ}Λ^{(v)}/σ_v² + I )^{−1} Σ_v Λ^{(v)ᵀ} x^{(v)}/σ_v²`.
M-step: `Λ^{(v)} = (Σ_i x_i^{(v)} E[z_i]ᵀ)(Σ_i E[z_i z_iᵀ])^{−1}`, `σ_v²` = residual variance.
Iterate to convergence; per-view variance explained = `1 − σ_v²/var(X^{(v)})`. `omni_mofa_lite`.

### 5.6 Program–GWAS association (scDRS-flavoured)
Given a program score `s_i` per cell and a gene-level GWAS z-vector `γ_g`: control gene sets are
drawn matching the program's mean-expression/co-expression bin; the association statistic is the
standardised program score `ẑ = (mean(s_program) − mean(s_control)) / sd_control`, aggregated
across cells to a per-sample z with a permutation null. `omni_program_gwas`.

---

## References
- UCell: Andreatta & Carmona 2021, *CSBJ* 19:3320.
- Scrublet: Wolock, Lopez & Klein 2019, *Cell Syst* 8:281.
- DeLong ER, DeLong DM, Clarke-Pearson 1988, *Biometrics* 44:837 (fast midrank: Sun & Xu 2014).
- Scanpy `score_genes`: Satija/Satpathy (Pagès control-bin approach).
- NMF updates: Lee & Seung 2001, *NIPS*.
- MOFA: Argelaguet et al. 2018, *Mol Syst Biol* 14:e8124.
- scDRS: Zhang et al. 2022, *Nat Genet* 54:1572.
- Moran: Moran 1950; Geary: Geary 1954; LISA: Anselin 1995; Ripley 1976.
- Stouffer 1949; Benjamini–Hochberg 1995; Benjamini–Yekutieli 2001; James–Stein 1961.
