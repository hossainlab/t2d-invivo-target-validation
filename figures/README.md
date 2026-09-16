# Figures: Reference-Framework Panels Only

This folder holds only the panels that mirror the analytical framework of the reference paper (Xu et al., *Phytomedicine* 2026, Figs 1–3). Everything else lives in `results/`:

| Location | Contents |
|---|---|
| `results/supplementary_figures/` | QC, diagnostic and sensitivity plots (bulk QC, scRNA demultiplexing/QC, pseudobulk DE counts, validation, prioritization) |
| `results/figure_exports/` | TIFF exports of main panels |
| `archive/superseded_figures/` | Outputs from superseded designs (3-group bulk contrast, 12-gene DEG run) |

## Fig 1: Bulk DEG + WGCNA (`Fig1_bulk_DEG_WGCNA/`)

| Panel | Reference paper | This project | Script |
|---|---|---|---|
| Fig1A | PCA before/after batch correction | PCA before/after ComBat (GSE15653 + GSE64998), T2D vs Control | 01 |
| Fig1B | Volcano of DEGs | T2D (16) vs lean Control (11), P < 0.05 and \|log2FC\| > 0.5 (M9) | 02 |
| Fig1C | Soft-threshold scale independence / mean connectivity | Same | 03 |
| Fig1D | Module dendrogram | Same | 03 |
| Fig1E | Module–trait heatmap | Same (traits T2D, Dataset, HbA1c) | 03 |
| Fig1F | MM vs GS in key module | MM vs GS for the T2D-associated key modules (T2D p < 0.1; M11) that contain 19-gene candidates. Other modules are in `results/supplementary_figures/bulk/` | 03/04 |
| Fig1G | Venn: DEG ∩ module ∩ phenotype gene set | Paper 3-way overlap (M13): DEGs (P < 0.05, \|log2FC\| > 0.5) ∩ all genes of the T2D-associated key modules (T2D p < 0.1) ∩ hyperglycemia gene set. Counts only; the caption lists the up/down genes and the DEG ∩ hub count. 27 samples, T2D vs lean (D1b) | 04 |

## Fig 2: Functional enrichment (`Fig2_enrichment/`)

| Panel | Reference paper | This project | Script |
|---|---|---|---|
| Fig2A | GO enrichment | GO over-representation of the 19 candidate genes (230 terms; top: response to insulin) | 05 |
| Fig2B | KEGG bubble plot | KEGG over-representation (9 pathways; top: AMPK signalling) | 05 |
| Fig2C | Pathway interaction network | `Fig2C_emap_KEGG.pdf` (KEGG) and `Fig2C_emap_GOBP.pdf` (GO BP) | 05 |
| Fig2D | pathview pathway diagram | `Fig2D_pathview/hsa04152.T2D_vs_Control_logFC.png`: AMPK signalling coloured by T2D vs Control logFC | 05 |
| Fig2_GSEA_hallmark | — (additional) | GSEA Hallmark on the full T2D vs Control ranking (22 sets; p53 pathway, UPR) | 05 |

## Fig 3: scRNA-seq, STZ vs Control (`Fig3_scRNA/`)

Selected map: **AMPK signalling** (stage 1, script 16). Reported axis is **candidate-anchored** (decision R23): **AMPKα2 (PRKAA2) → PGC-1α (PPARGC1A)** with **CCND1** as output node — both anchors are Fig 1 candidates. Fig 3D–G show that axis (script 23); the pair is not co-localized or correlated, consistent with post-translational control. Fig 3H–J show the downstream partner **Cdkn1a (p21) – Ccnd1**, the only pair that satisfies all three of the paper's single-cell criteria (script 18). Supplementary: AMPK eyeballed-axis panels and the p53 per-cell-type score (`results/supplementary_figures/pathway_selection/`), shortlist-wide confirmation matrix (`17_edge_confirmation.pdf`), selection criteria (`16_selection_criteria.pdf`), chain audit (`results/supplementary_figures/consistency/22_chain_consistency.pdf`). The earlier 19-gene candidate panels are in `results/supplementary_figures/candidates_19gene/`.

| Panel | Reference paper | This project | Script |
|---|---|---|---|
| Fig3A | UMAP cell atlas | `Fig3A_umap_celltype_<liver/kidney>` | 08 |
| Fig3B | Marker bubble plot | `Fig3B_marker_bubble_<tissue>` | 08 |
| Fig3C | Cell proportions per sample | `Fig3C_cell_proportion_<tissue>` | 08 |
| Fig3D | Rac1/Pak1 expression by group in neutrophils (cell-level Wilcoxon) | `Fig3D_anchored_axis_expression_<tissue>`: AMPK-map axis genes (Prkaa1/2, Stk11, Camkk2, Sirt1, Ppargc1a, Irs2, Igf1, Fbp1, Ccnd1) in the key cell type; Fig 1 candidates marked * | 23 |
| Fig3E | Pathway elevation in the key cell type | `Fig3E_anchored_score_<tissue>`: UCell score of the Fig 1 candidates that lie on the AMPK map, split by human direction, per cell type | 23 |
| Fig3F | Rac1+Pak1 co-localization UMAP | `Fig3F_anchored_colocalization_<tissue>`: Ppargc1a+Ccnd1 joint density and % co-expressing cells | 23 |
| Fig3G | Rac1–Pak1 correlation in neutrophils | `Fig3G_anchored_correlation_<tissue>`: cell-level Pearson, per-mouse Pearson, and a detection-matched background null — the pair does **not** correlate | 23 |
| Fig3H | — (downstream readout) | `Fig3H_p53_partner_expression_<tissue>`: Cdkn1a, Ccnd1, Phlda3, Zmat3, Bax, Serpine1, Mdm2, Trp53 in the key cell type | 18 |
| Fig3I | — | `Fig3I_p53_partner_colocalization_<tissue>`: Cdkn1a+Ccnd1 joint density and % co-expressing cells | 18 |
| Fig3J | — | `Fig3J_p53_partner_correlation_<tissue>`: the pair that does satisfy all three criteria in kidney endothelium, with its background null | 18 |

## Fig 4: Machine learning (`Fig4_ML/`)

Not a reference-paper panel — an addition (script 21, `docs/results_ML.md`, decision R22). Classification of T2D vs lean from the 19-gene candidate panel, built so the small-n failure mode is visible: optimistic CV, honest nested CV, a permutation null, and an independent cohort. Assembled panel: `Fig4_ML.pdf`.

| Panel | Contents | Script |
|---|---|---|
| Fig4A | AUC of panel-fixed CV vs nested CV (selection re-run inside folds) vs label-permutation null, for LASSO / random forest / SVM | 21 |
| Fig4B | ROC in the independent cohort GSE23343 (10 T2D / 7 NGT), three models plus a parameter-free signature score | 21 |
| Fig4C | LASSO selection frequency over 100 folds (gene-selection stability) | 21 |
| Fig4D | Random-forest permutation importance | 21 |
| Fig4E | External predicted probabilities by true group | 21 |
| Fig4F | Top consensus genes, z-scored, in both cohorts | 21 |

## Fig 5: Mendelian randomization (`Fig5_MR/`)

Exposure: eQTLGen whole-blood cis-eQTLs. Outcome: type 2 diabetes (GCST006867, Xue et al. 2018). Primary analysis is LD-clumped (1000G EUR, Ensembl); no candidate is causally supported. Details and caveats in `docs/results_MR.md` (decisions R20, R21).

| Panel | Reference paper | This project | Script |
|---|---|---|---|
| Fig5A | Forest plot of the causal estimate for PAK1 on RA | `Fig5A_forest.pdf`: all 25 instrumentable candidate genes, OR per SD of genetically predicted expression — **superseded first pass (distance-pruned)**; SERPINE1, PRKAA1 and PPARGC1A do not survive LD clumping (see Fig5F) | 19 |
| Fig5B | Radial / SNP plot | `Fig5B_scatter_<gene>.pdf`: SNP effect on expression vs on T2D with the IVW slope | 19 |
| Fig5C | Per-SNP forest | `Fig5C_leaveoneout_<gene>.pdf`: leave-one-out IVW | 19 |
| Fig5D | Steiger directionality | `Fig5D_funnel_<gene>.pdf` (funnel); Steiger results are tabular, `results/mr/mr_steiger.csv` | 19 |
| Fig5E | — (added) | `Fig5E_tissue_forest.pdf`: MR with GTEx liver instruments (kidney cortex has none) | 20b |
| Fig5F | — (added) | `Fig5F_blood_clumping_comparison.pdf`: LD-clumped vs distance-pruned ORs; shows why the first pass over-called (decision R21) | 20b |
