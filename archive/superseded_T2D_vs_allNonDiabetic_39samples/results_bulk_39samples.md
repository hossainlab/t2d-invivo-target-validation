# Bulk Liver Results: T2D vs Control (GSE15653 + GSE64998)

**Status (2026-09-14):** Re-run under decisions M8–M11. Sections 1–4 are final. Sections 5–7 are pending until scripts 03 (module preservation), 04, 05 and 11 finish.
**Scripts:** `scripts/01`–`05`, `11`. All decisions are in `docs/decisions_log.md`.
**Comparison (D1):** T2D = obese type 2 diabetic (n = 16) vs Control = all non-diabetic, i.e. lean/non-obese + obese non-diabetic (n = 23).

| Decision | Setting |
|---|---|
| M8 | Raw CEL + RMA per platform, **no expression filter** |
| M9 | DEG = nominal **P < 0.05 and \|log2FC\| > 0.5** |
| M10 | Candidates = **DEGs ∩ WGCNA hub genes** (\|kME\| > 0.7, \|GS\| > 0.2); hyperglycemia gene set = annotation |
| M11 | Key modules = **all modules correlated with T2D at p < 0.1** (not dataset-driven) |

---

## 1. Data and preprocessing (Fig 1A)

| Item | Result |
|---|---|
| Arrays | GSE15653 (HG-U133A, 18) + GSE64998 (HuGene 1.1 ST, 21), raw CEL files, RMA per platform |
| Genes | 12,650 and 18,865 annotated genes → **11,972 common genes, all kept** |
| Samples | 39 (GSE15653: 9 Control / 9 T2D; GSE64998: 14 Control / 7 T2D); no outliers met the ≥2-flag rule |
| Batch | ComBat (batch = dataset, condition protected) for PCA / WGCNA; DE uses dataset as a covariate |

## 2. Differential expression (Fig 1B)

Model: limma `~condition + dataset`.

| Analysis | Up | Down |
|---|---|---|
| FDR < 0.05 (any model: primary, paper-style ComBat, + sex, per dataset) | 0 | 0 |
| **Nominal P < 0.05, \|log2FC\| > 0.5 (used downstream)** | **15** | **12** |
| Supplementary: obese T2D vs lean only (FDR < 0.05, \|log2FC\| > 0.5) | 31 | 31 |

**The 27 DEGs:**
- **Up:** CDKN1A, GOLM1, HPS5, CCND1, PLIN2, SPP1, HBB, SERPINE1, ANGPTL8, INHBE, FABP4, ID1, IL32, CHI3L1, PGGHG
- **Down:** RETREG1, CYP2C19, SLC16A4, SLC16A7, VIL1, SLITRK3, GPR88, RND2, CYP26A1, P4HA1, SDS, TSPAN8

HBB is probably blood contamination.

**Why this differs from the earlier project** (`invivo-target-validation_old`, 197 nominal DEGs): that analysis used GEO series-matrix values, while this one starts from raw CEL files. Between the two runs, logFC rank correlation is 0.68 and −log10 P rank correlation is 0.34 (decision M8).

## 3. WGCNA (Fig 1C–F)

Signed-hybrid network with bicor; soft power 5; 11,972 genes; modules merged at 0.25.

**Key modules (T2D p < 0.1, M11):**

| Module | Genes | r (T2D) | p | r (Obese) | Candidates |
|---|---|---|---|---|---|
| green | 860 | −0.44 | 0.005 | −0.61 | — |
| blue | 1,546 | +0.34 | 0.037 | +0.36 | — |
| turquoise | 2,065 | +0.31 | 0.058 | +0.51 | CDKN1A, CCND1, IL32 |
| grey60 | 98 | −0.29 | 0.077 | −0.40 | RETREG1 |
| brown | 1,308 | −0.27 | 0.093 | −0.34 | — |
| yellow | 1,001 | −0.27 | 0.094 | −0.12 | — |

- **912 hub genes** (\|kME\| > 0.7, \|GS\| > 0.2) across these modules.
- Most modules track obesity at least as strongly as T2D. Yellow is the most T2D-specific.
- Fig 1F shows the MM–GS panels for turquoise and grey60. The other key modules are in `results/supplementary_figures/bulk/`.

## 4. Candidate genes: DEGs ∩ hub genes (Fig 1G)

| Gene | T2D vs Control logFC | P | Module | kME | GS | GSE15653 logFC | GSE64998 logFC | Hyperglycemia gene set |
|---|---|---|---|---|---|---|---|---|
| **CDKN1A** | +0.62 | 1.2×10⁻⁴ | turquoise | 0.71 | 0.57 | +0.55 | +0.69 | no |
| **CCND1** | +0.61 | 0.001 | turquoise | 0.72 | 0.48 | +0.67 | +0.54 | yes |
| **IL32** | +0.69 | 0.028 | turquoise | 0.75 | 0.37 | +0.39 | +0.98 | no |
| **RETREG1** | −0.50 | 0.001 | grey60 | 0.77 | −0.51 | −0.32 | −0.68 | no |

- All four change in the same direction in both discovery datasets.
- The paper's 3-way overlap (DEGs ∩ key-module members ∩ hyperglycemia gene set, 772 genes, 618 in the data) is **2 genes**, both non-hubs. This is in line with the earlier project's own 3-way result (0–1 genes).
- **Venn counts (Fig 1G):** 27 DEGs, 912 hub genes, 4 shared (23 DEG-only, 908 hub-only).

**Mouse evidence (scRNA STZ vs Control; details in `docs/results_scRNA.md` §6):**

| Candidate | Evidence |
|---|---|
| CDKN1A | **Supported.** Up in kidney endothelium and liver cholangiocytes (FDR < 0.01) and liver LSEC (nominal); opposite direction in liver NK cells |
| CCND1 | **Supported.** Up in kidney endothelium (FDR 0.0015), liver cholangiocytes and LSEC (nominal) |
| RETREG1 | Not supported |
| IL32 | No mouse ortholog; cannot be tested in the GV1 mouse model |

## 5. Module preservation (GSE15653 → GSE64998)

All six key modules are preserved. Zsummary values: blue 13.7, turquoise 13.7, yellow 14.8, green 16.4, grey60 19.0, brown 34.4. Zsummary > 10 counts as strong preservation, so the module structure is not an artefact of one dataset.

## 6. Functional enrichment (Fig 2)

With only 4 candidate genes, over-representation analysis is skipped. The Fig 2 evidence is GSEA on the full T2D vs Control ranking.

| Collection | Significant (padj < 0.05) | Top sets (NES) |
|---|---|---|
| Hallmark | 18, all up in T2D | Cholesterol homeostasis (2.21), **p53 pathway (2.18)**, TGF-β signalling (1.98), apoptosis (1.96), apical junction (1.96), IL6–JAK–STAT3 (1.95), UV response up, hypoxia, myogenesis, **unfolded protein response (1.68)**, TNFα via NF-κB (1.65), MYC targets v2, EMT (1.50), mTORC1 (1.49) |
| Reactome | 44 | **IRE1α activates chaperones (2.44)**, RMTs methylate histone arginines, RNA Pol I / rRNA regulation, granulopoiesis TFs, scavenging by class A receptors, fat-soluble vitamin metabolism, DNA methylation, lipoprotein assembly/clearance, UPR, ECM proteoglycans |
| KEGG (MEDICUS) | 0 | — |

The p53 pathway ties the candidates together: CDKN1A is a p53 target and CCND1 is part of the G1/S control that p53 regulates. The ER stress / UPR and TGF-β / EMT themes are also present.

## 7. External validation: GSE23343 (10 T2D vs 7 NGT; raw CEL + RMA)

| Test | Result |
|---|---|
| **CCND1** | Concordant ↑, logFC 0.54, **P = 0.016**, AUC 0.81 → **replicated** |
| **CDKN1A** | Concordant ↑, logFC 0.40, P = 0.064, AUC 0.74 → same direction, borderline |
| IL32 | Concordant ↑, logFC 0.23, P = 0.47 → not replicated |
| RETREG1 | Discordant (↑ in GSE23343) → not replicated |
| 27 DEGs, direction | 19 / 27 concordant, sign test **p = 0.026** |
| 27-DEG signature score | AUC 0.69 (95% CI 0.38–0.99), Wilcoxon p = 0.23 |
| Genome-wide logFC correlation | ρ = −0.27 (cohort / biopsy differences; see decision R6) |

The genome-wide profiles of the two cohorts disagree. Even so, the top T2D DEGs agree in direction more often than chance, and CCND1 replicates independently with CDKN1A close behind.

## 8. Limitations

- **Weak bulk signal.** No gene reaches FDR < 0.05; DEGs are nominal.
- **Near-tie module selection.** Candidates depend on selecting modules at p < 0.1: turquoise has p = 0.058.
- **Obesity confounding.** Sex is partly confounded with group, and key modules correlate with obesity as strongly as with T2D.
- **Replication.** The earlier signature did not replicate in GSE23343 (cohort and biopsy differences), so human evidence is supportive only. Cross-species concordance (CDKN1A, CCND1) carries more weight.
