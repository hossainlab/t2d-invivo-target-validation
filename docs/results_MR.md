# Mendelian randomization: are the candidates causal for type 2 diabetes?

**Status (2026-09-16):** complete (decisions R20, R21). **Scripts:** `19_MR_druggability.R` (first pass, distance-pruned), `20a_tissue_instruments.R` + `20b_tissue_MR.R` (tissue instruments and LD clumping), shared estimators in `mr_functions.R`.
**Reference-paper analogue:** Fig 4 — the paper ran MR of *PAK1* expression on rheumatoid arthritis (IVW, medians, modes, Cochran's Q, MR-Egger intercept, Steiger) and combined it with the Finan et al. 2017 druggable-genome list to nominate PAK1 as a druggable causal target.

> **Headline:** with LD-clumped instruments, **no gene of the selected AMPK axis or of the p53 arrest axis is causally supported for type 2 diabetes**. Only SREBF1 and SIRT1 survive multiple testing, and both are whole-blood instruments pointing against hepatic biology, so neither supports a target claim. The earlier distance-pruned result (PRKAA1, PPARGC1A, SERPINE1, CDKN1A) did not hold up and is reported as a superseded sensitivity analysis.

## 1. Data and instruments

| Element | This study |
|---|---|
| Exposure (primary) | cis-eQTLs, eQTLGen whole blood, n = 31,684 (Võsa et al. 2021), FDR < 0.05 file, P < 5e-8 |
| Exposure (attempted) | GTEx liver (n = 208) and kidney cortex (n = 73) from the eQTL Catalogue — see §5, not usable |
| Outcome | Type 2 diabetes, Xue et al. 2018, GWAS Catalog **GCST006867** (62,892 cases / 596,424 controls, European), harmonised summary statistics |
| Genes | 19-gene candidate set + p53 axis (CDKN1A, PHLDA3, ZMAT3, BAX, MDM2, TP53) + AMPK core (PRKAA1, PRKAA2, SIRT1, STK11, CAMKK2, SREBF1, FASN) = 32 requested, **24 with usable clumped instruments** |
| Instrument pruning | **LD clumping, r² < 0.1 within 500 kb, 1000 Genomes EUR via the Ensembl REST LD API** (primary). Distance pruning at 250 kb (script 19) is kept as a sensitivity analysis |
| Harmonisation | alleles aligned to the eQTL effect allele (strand flips handled); palindromic SNPs with EAF 0.42–0.58 dropped; eQTLGen Z converted to beta/se with the outcome allele frequency (Zhu et al. 2016); F ≥ 10 |
| Estimators | IVW (multiplicative random effects), MR-Egger, weighted median, weighted mode; Wald ratio for single-instrument genes. Cross-checked against the `MendelianRandomization` package (identical point estimates) |

Genes without a genome-wide-significant cis-eQTL in blood could not be tested: **VSNL1, PNPLA3, IGFBP1, IGF1, PRKAA2, SREBF2, STK11, GAS6**.

## 2. Primary result — LD-clumped (`results/mr/mr_results_blood_ldclumped.csv`)

OR per SD of genetically predicted expression; FDR across the 24 primary tests.

| Gene | SNPs | OR (95 % CI) | P | FDR | Q P |
|---|---|---|---|---|---|
| **SREBF1** | 2 | **0.903 (0.872–0.935)** | 9.0e-9 | **2.2e-7** | 0.70 |
| **SIRT1** | 4 | **1.027 (1.011–1.044)** | 9.5e-4 | **0.011** | 0.88 |
| SERPINE1 | 3 | 0.813 (0.680–0.972) | 0.023 | 0.19 | 0.26 |
| FASN | 1 | 1.044 (0.997–1.093) | 0.068 | 0.41 | — |
| TP53 | 1 | 0.766 (0.545–1.075) | 0.12 | 0.50 | — |
| FBP1 | 3 | 0.891 (0.770–1.032) | 0.12 | 0.50 | 0.032 |
| **PPARGC1A** | 12 | 0.983 (0.959–1.007) | 0.16 | 0.55 | 0.63 |
| **PRKAA1** | 5 | 0.865 (0.682–1.098) | 0.23 | 0.62 | 1.1e-4 |
| **CDKN1A** | 9 | 1.025 (0.956–1.099) | 0.49 | 0.78 | 0.043 |
| **CCND1** | 5 | 1.033 (0.928–1.149) | 0.55 | 0.78 | 0.034 |
| CAMKK2 | 5 | 1.008 (0.919–1.106) | 0.86 | 0.86 | 3.1e-9 |

**Reading.**
- **The axis genes are null.** CDKN1A (OR 1.03, P 0.49) and CCND1 (1.03, P 0.55) give no evidence of causality for T2D, and neither does the AMPK core (PRKAA1 P 0.23, PPARGC1A P 0.16, CAMKK2 P 0.86 with extreme heterogeneity). MR therefore does **not** nominate a druggable causal target in this study, unlike the reference paper's PAK1 result.
- **SREBF1 and SIRT1 survive** but cannot be used for a target claim: they are whole-blood instruments for genes whose relevant biology is hepatocyte-specific, and their directions (higher SREBF1 protective, higher SIRT1 harmful) run against that biology. They are reported, not interpreted.
- The AMPK selection in `docs/pathway_selection.md` rests on the human bulk transcriptome, the KGML edge analysis and the planned protein-level work; **it gains no genetic support here, and the write-up must say so.**

### 2.1 Why this differs from the first pass (`mr_blood_pruning_comparison.csv`, Fig 5F)

| Gene | LD-clumped: SNPs, OR, P | Distance-pruned (script 19): SNPs, OR, P |
|---|---|---|
| SREBF1 | 2, 0.903, 9.0e-9 | 11, 0.886, 1.8e-10 |
| SIRT1 | 4, 1.027, 9.5e-4 | 10, 1.039, 2.1e-4 |
| SERPINE1 | 3, 0.813, 0.023 | 1, 0.729, 2.5e-3 |
| PPARGC1A | 12, 0.983, 0.16 | 3, 0.953, 7.3e-3 |
| PRKAA1 | 5, 0.865, 0.23 | 1, 0.766, 3.5e-3 |
| CDKN1A | 9, 1.025, 0.49 | 4, 1.107, 0.016 |

Distance pruning keeps one SNP per 250 kb window, which for PRKAA1 and SERPINE1 meant a **single lead SNP** carrying the whole estimate. Once the locus is clumped properly and the other independent signals enter, those estimates attenuate to the null. The earlier reading — "MR nominates the AMPK arm" — was an artefact of that pruning and is withdrawn (decision R21).

## 3. Sensitivity and limitations

- **Steiger**: every gene passes directionality (variance explained in expression ≫ in T2D); `mr_steiger.csv`.
- **Heterogeneity**: PRKAA1 (Q P 1.1e-4), CAMKK2 (3.1e-9), CDKN1A (0.043) and CCND1 (0.034) have heterogeneous instruments, a further reason not to read their point estimates.
- **Egger intercepts** are null for the surviving genes (SREBF1 0.45, SIRT1 0.16 in the first pass).
- **Instrument recovery**: of 374 clumped instrument rsIDs, 116 were present in the outcome GWAS, so some genes are tested on fewer SNPs than were selected.
- **Tissue**: whole blood only — see §5. This remains the main limitation.
- The exposure is *transcript abundance*, so these estimates say nothing about AMPK or p53 phosphorylation, which is where the human bulk data locate the defect.
- LD clumping used the Ensembl 1000G EUR panel through its REST API; no SNP failed to resolve in the tissue set, and unresolved lead SNPs are flagged in `instruments_blood_ldclumped.csv`.

## 4. Druggability

The paper's step here was to combine a causal estimate with the Finan et al. 2017 druggable genome. Since no candidate is causally supported, no target nomination is made from MR. AMPK remains a clinically drugged node (metformin upstream; direct activators A-769662, PF-06409577), which is the rationale for the metformin arm of the in vivo experiment — but that rationale is pharmacological and mechanistic, **not** genetic.

## 5. Tissue-specific instruments were attempted and are not available (script 20a)

The whole-blood limitation was tested rather than assumed. Tissue cis-eQTLs were pulled by remote tabix from the eQTL Catalogue (the GTEx portal's own eQTL API returns empty results for every gene, including positive controls, and the eQTL Catalogue REST API now returns HTTP 410, so the indexed FTP files are the remaining route):

| Dataset | n | cis rows for the 32 genes | Genes with cis-eQTL P < 5e-8 | Genes at relaxed P < 1e-5 |
|---|---|---|---|---|
| GTEx liver (QTD000266) | 208 | 53,299 | **1** (GAS6, P 1.1e-9) | 3 more: VSNL1, PLIN1, CTSD |
| GTEx kidney cortex (QTD000261) | 73 | 42,540 | **0** | 0 |

**No axis gene — CDKN1A, CCND1, PRKAA1, PPARGC1A, SREBF1, SIRT1 — has a usable liver or kidney instrument.** At n = 208 and n = 73 these datasets are not powered for cis-MR. The three liver genes that did clear the relaxed threshold (CTSD, PLIN1, VSNL1; GAS6 was lost at harmonisation) are all null for T2D (P 0.43–0.55, `mr_results_tissue.csv`). Raw pulls are in `results/mr/tissue_eqtl_<dataset>.csv`.

## 6. Figures (`figures/Fig5_MR/`)

| Panel | Contents |
|---|---|
| Fig5A | Forest plot, all 25 genes, distance-pruned first pass (superseded; kept for the comparison) |
| Fig5B | SNP-effect scatter with the IVW slope (SREBF1, SIRT1, PPARGC1A, CDKN1A) |
| Fig5C | Leave-one-out IVW |
| Fig5D | Funnel plot of the Wald ratios against instrument strength |
| **Fig5E** | Tissue-instrument forest (GTEx liver; kidney cortex has no instruments) |
| **Fig5F** | LD-clumped vs distance-pruned ORs — the panel that shows why the first pass over-called |
