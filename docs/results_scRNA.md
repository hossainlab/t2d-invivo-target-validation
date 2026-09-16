# scRNA-seq Results: GSE244475 Liver and Kidney, STZ vs Control

**Status:** Phases 6–8 and candidate prioritization complete (2026-09-14), integrated with the final bulk design: 27 samples, T2D vs lean (D1b), 19-gene 3-way candidate set (M13).
**Scripts:** `06a`, `06b`, `07`, `08` (with `docs/annotation/<tissue>_overrides.csv`), `09`, `10`, `13`, `12`.
**Decisions:** `docs/decisions_log.md` entries S1–S5, Q3, R7 (scRNA); D1b, M13, R12–R14 (integration with bulk).
**Design:** C57BL/6J male mice, 3 Control vs 4 STZ. Each 10x lane holds one mouse, with heart, kidney, liver and spleen pooled by hashtag.

---

## 1. Tissue assignment (hashtag + RNA)

The hashtag (HTO) signal is weak and uneven between cell types:

| Cell population | HTO label correct / status |
|---|---|
| Liver-restricted clusters | LSEC 96%, Kupffer 89%, hepatocyte 85% |
| Kidney endothelium | 90% |
| Kidney tubular epithelium | 64–79% HTO-Negative; singlet calls spread over all four tissues |
| Stromal cells | Mostly HTO-Negative |

A hybrid assignment was used instead (script 06b):

| Step | Rule |
|---|---|
| RNA override | Clusters dominated by kidney-tubular or hepatocyte programs are assigned by RNA (4 kidney and 3 liver clusters) |
| HTO | HTO singlets keep their tag unless the RNA kNN probability for that tissue is < 0.2 |
| kNN rescue | HTO-Negatives are rescued only in clusters whose CV accuracy is ≥ 0.90 (12 clusters: immune, LSEC, Kupffer, kidney and heart endothelium) |

kNN CV accuracy: overall 84% (93% at pmax ≥ 0.8); Spleen 92%, Liver 88%, Kidney 40%, Heart 28%.

## 2. QC (script 07)

| Tissue | Assigned | Pass QC | After doublets / ambient | After annotation exclusions |
|---|---|---|---|---|
| Liver | 20,713 | 19,432 | 16,991 | **16,861** |
| Kidney | 7,960 | 6,639 | 6,248 | **6,224** |

QC rules:
- ≥ 200 genes; genes detected in ≥ 3 cells.
- Hemoglobin < 3%.
- Mito ≤ max(median + 3 MAD, 10%), capped at 20% (liver) / 30% (kidney). This is never stricter than the reference paper's 10%.
- DecontX for ambient RNA; scDblFinder for doublets.
- Harmony by mouse; UMAP on 20 PCs, as in the reference paper.

## 3. Annotation (script 08, Fig 3A–C)

Automated marker-set scoring with a SingleR cross-check. Eight clusters were corrected manually from canonical markers:

| Tissue | Cluster | Automatic label | Final label | Evidence |
|---|---|---|---|---|
| Liver | 17 | LSEC | Exclude | Mito/Malat1-dominated, low quality |
| Liver | 12 | Macrophage | Monocyte (non-classical) | Fcgr4, Clec4a3, Gngt2 |
| Liver | 10 | Proliferating | cDC | Psmb8/9, H2 genes, no Mki67 |
| Kidney | 11 | Plasma cell | LoH | Cryab, Epcam, Apela |
| Kidney | 17 | Podocyte | Exclude | Cardiomyocyte transcripts (heart cells mis-assigned) |
| Kidney | 8 | CD_PC | Distal tubule | Atp1b1, Atp1a1, Wfdc2 |
| Kidney | 13 | Fibroblast | Urothelium | Krt7, Sprr1a, Foxq1 |
| Kidney | 15 | PT | Tubular epithelium (other) | Cdh16, Pkhd1, Bicc1 |

**Liver cell types:** B cells (dominant), T, NK, monocytes, cDC, pDC, neutrophils, Kupffer cells, other macrophages, LSEC, HSC/fibroblasts, cholangiocytes, hepatocytes (rare: 69 cells; 10x largely excludes hepatocytes).

**Kidney cell types:** endothelium (dominant), macrophages, PT, LoH, distal tubule, CD-IC, urothelium, fibroblast/pericyte, T, B, cDC, neutrophils.

**Cell numbers limit what can be tested:**
- Control_3 has only 185 kidney cells.
- Fibroblasts come mainly from STZ_4: 1,354 / 2,024 in liver and 523 / 591 in kidney.

## 4. Differential expression, STZ vs Control (script 09)

Pseudobulk per mouse × cell type (≥ 20 cells in ≥ 3 mice per group), edgeR quasi-likelihood test with a DESeq2 check; FDR < 0.05 and |log2FC| > 0.5.

| Tissue | Cell type | DEGs (up / down) | DESeq2-concordant | Without STZ_4 (same direction, P < 0.05) |
|---|---|---|---|---|
| Liver | NK cells | 62 (18 / 44) | 52 | 60 / 62 |
| Liver | Cholangiocytes | 17 (3 / 14) | 17 | 16 / 17 |
| Liver | B cells | 6 (5 / 1) | 5 | 6 / 6 |
| Liver | T cells | 2 (1 / 1) | 1 | 2 / 2 |
| Liver | HSC/fibroblast, LSEC, monocyte, neutrophil, cDC, pDC | 0 | – | – |
| Kidney | Endothelium (only eligible type) | 88 (38 / 50) | 55 | 87 / 88 |

**Top genes:**
- **Liver cholangiocytes:** up Phlda3, Cdkn1a; down Hsph1, Hspa8, Klf6, Ccn1, Atf3, Hspa1b, Dnajb1, Egr1, Jun.
- **Liver NK cells:** up Nrgn, Cd160, Cd7, Tox; down Irf8, Eomes, Nfkbia, Cdkn1a, Nr4a1/3, Gzma, Ifng.
- **Kidney endothelium:** up Phlda3, Ifi27l2a, Acy3, Ccnd1, Lgals1, Zmat3, Ltc4s; down Esm1, Hspa1b, Hspd1, Hspb1, Spry1, Efna1.

**Interpretation:**
- **p53 / cell-cycle arrest is the convergent signal** (Cdkn1a, Phlda3, Zmat3, Ccnd1). It matches the human bulk results: CDKN1A and CCND1 up in T2D, and the Hallmark p53 pathway at the top of GSEA.
- **Many down-regulated genes are immediate-early / heat-shock genes.** This is a dissociation-stress signature and may partly reflect processing differences between lanes, so those genes are down-weighted.

## 5. Differential abundance (propeller; exploratory)

Results after annotation overrides. Proportions are means of per-mouse cell-type fractions.

| Tissue | Cell type | Control → STZ | All mice (3 vs 4): P / FDR | Without STZ_4 (3 vs 3): P / FDR | Verdict |
|---|---|---|---|---|---|
| Liver | Hepatocytes | 0.14% → 0.68% | 0.003 / **0.039** | 0.030 / 0.17 | Only FDR hit; very rare population (69 cells), interpret with care |
| Liver | Other macrophages | 0.5% → 1.6% | 0.022 / 0.13 | 0.053 / 0.17 | Suggestive |
| Liver | B cells | 58.5% → 35.0% | 0.029 / 0.13 | 0.029 / 0.17 | Consistent, nominal |
| Liver | HSC/fibroblasts | 3.6% → 16.8% | 0.068 / 0.22 | not in top | Driven by STZ_4 |
| Kidney | Fibroblast/pericyte | 1.0% → 8.8% | 0.037 / 0.30 | 0.20 / 0.41 | Driven by STZ_4 |
| Kidney | LoH | 1.5% → 3.9% | 0.088 / 0.30 | 0.018 / 0.24 | Suggestive |
| Kidney | Macrophages | 7.9% → 17.3% | 0.11 / 0.30 | 0.064 / 0.41 | Suggestive |

Kidney podocytes are no longer tested: the former "podocyte" cluster was cardiomyocyte contamination and was excluded. No abundance change is robust. The apparent fibroblast expansion in both organs depends on STZ_4.

## 6. Mouse evidence for the bulk candidate genes (27-sample design, 19-gene set; decisions D1b, M13, R14)

**Bulk candidates:** the reference paper's 3-way overlap. DEGs (T2D vs lean, P < 0.05, |log2FC| > 0.5) ∩ genes in the T2D-associated WGCNA modules ∩ hyperglycemia gene set.

**Genes that cannot be tested in mouse:**
- VSNL1 and PLIN1: no 1:1 mouse ortholog detected in either tissue.
- IGFBP1: not detected in the liver data.
- PNPLA3: not detected in the kidney data.

Most hepatocyte-restricted genes (IGFBP1, PNPLA3, FBP1, KHK, IGF1) are hard to test because the 10x data contain almost no hepatocytes.

## 7. Bulk–scRNA integration (script 10; Fig 3D–E)

Tier 1 = pseudobulk FDR < 0.05, same direction as in human. Tier 1b = nominal P < 0.05, same direction. Tier 2 = hub or gene-set gene with cell-type-restricted expression, but not DE in STZ.

| Candidate (human direction) | Liver | Kidney | GSE23343 |
|---|---|---|---|
| **CCND1** (↑) | Tier 1b: cholangiocytes, logFC 1.10, FDR 0.12 | **Tier 1: endothelium, logFC 0.97, FDR 0.002** | **Replicated** (P = 0.016, AUC 0.81) |
| PPARGC1A (↓) | Tier 1b: cholangiocytes, logFC −1.03 | Tier 2 | Concordant, n.s. |
| GPX1 (↑) | Tier 1b: cholangiocytes, logFC +0.37 | not DE | Concordant, n.s. |
| IRS2 (↓) | Tier 1b: pDC, logFC −0.86 | not DE | Discordant |
| PNPLA3 (↑) | Tier 2 (not DE-testable) | no ortholog in data | **Replicated** (P = 0.017, AUC 0.74) |
| IGFBP2 (↓) | Tier 2 | Tier 2 | Concordant, P = 0.061 |
| FBP1, IGF1, ENPP1, GAS6, SERPINE1 | Tier 2 | not DE / Tier 2 | Mixed |
| SREBF2, CTSD, HMOX1, LZTFL1, KHK | not supported | not supported | Mixed |

- **Signature scores** (human 19-gene and top-DEG signatures scored in mouse cells): no cell type reaches padj < 0.05.
  - The 19-gene down-signature is nominally higher in STZ kidney PT (p = 0.005, padj 0.07).
  - The same down-signature is lower in STZ liver neutrophils (p = 0.017).
- **Conclusion:** cross-species support is gene-level, not signature-wide. **CCND1** is the only candidate with FDR-level mouse support (kidney endothelium) plus human replication.

## 8. Co-localization and correlation (script 13; Fig 3F–G)

The focus cell type is the one with the most candidates changed in the human direction (nominal P < 0.05).

| Tissue | Focus cell type (concordant candidates) | Top pair | Per-mouse pseudobulk Pearson (n = 7) | Co-expressing cells, Spearman |
|---|---|---|---|---|
| Liver | Cholangiocytes (Ccnd1, Gpx1, Ppargc1a) | Ppargc1a–Lztfl1 | r = 0.82, p = 0.025 | ρ = 0.56, p = 5×10⁻⁴ (35 cells) |
| Liver | Cholangiocytes | Ccnd1–Gpx1 | r = 0.41, p = 0.36 | ρ = 0.17, p = 0.001 (369 cells) |
| Kidney | Endothelium (Ccnd1) | Ccnd1–Lztfl1 | r = 0.73, p = 0.06 | ρ = 0.30, p = 5×10⁻⁴ (127 cells) |

Correlations based on 7 mice are underpowered; cell-level results are secondary evidence.

## 9. Candidate prioritization (script 12)

Weighted evidence score, reported as % of the maximum. Literature scores still need to be filled in by the team.

**Pools:**
- **A:** bulk 3-way candidates.
- **B:** scRNA-first genes, i.e. mouse STZ DEGs with a 1:1 human ortholog that are not in the bulk set.

| Rank | Liver | Kidney |
|---|---|---|
| 1 | **CCND1** 71% (A) | **CCND1** 82% (A) |
| 2 | **CDKN1A** 67% (B; Cdkn1a FDR 0.009 in cholangiocytes, human DEG FDR 0.016) | **CDKN1A** 67% (B; Cdkn1a FDR 0.008 in endothelium) |
| 3 | HMOX1 58% (A) | ENPEP 64% (B) |
| 4 | PNPLA3 58% (A; human replication, no mouse DE) | VIM 60% (B) |
| 5 | CD7 56% (B) | BAX 56% (B; p53 target) |

**Recommended wet-lab panel for the GV1 Ob-T2D model:**

| Candidate | Priority | Evidence | Assays |
|---|---|---|---|
| **CCND1 / cyclin D1** | Primary | Human T2D liver DEG and hub gene, replicated in GSE23343, mouse Tier 1 | Kidney and liver: qRT-PCR, Western blot, IHC/IF with CD31/Emcn (endothelium) and CK19 (cholangiocytes) |
| **CDKN1A / p21** | Primary | Strongest human DEG outside the 3-way set; FDR-level mouse evidence in both organs; p53 pathway | Same as CCND1 |
| **PNPLA3**, **PPARGC1A**, **IGFBP2** | Secondary, liver | Human evidence and liver biology (PNPLA3 replicated); hepatocyte genes poorly captured by scRNA | Whole-liver qRT-PCR |
| **HMOX1** | Optional | Hub gene in the oxidative-stress module | Tissue qRT-PCR |

Expected outcome: changed in Ob-T2D vs control in the human direction, and attenuated by GV1.

## 10. Pathway confirmation in the mouse atlas (Fig 3D–G; decision R19)

Full details are in `docs/pathway_selection.md`. Stage 1 selected **AMPK signalling (hsa04152)** from the human bulk data alone (six scored criteria, script 16). Stage 2 (script 17) tested every axis candidate of all eight shortlisted KEGG maps against the reference paper's three Fig 3 criteria, in every pseudobulk-eligible cell type, BH-corrected within a tissue.

**Screen result — no map-drawn axis is confirmable here.** 81 pairs tested; 22 pair × cell-type combinations cleared co-localization and correlation, but **none had both nodes changed in the human direction**, and the largest correlation was R = 0.22 (the paper's benchmark: R = 0.34 with both nodes changed). Gene level: 478 tests, 17 concordant, concentrated in liver cholangiocytes, liver LSEC and kidney endothelium.

Three structural reasons, all reportable: AMPK/p53/SIRT1 activity is post-translational (the human map's PRKAA2 → PPARGC1A edge is discordant: subunit up, target down); the AMPK-relevant parenchyma is missing (69 hepatocytes, 150/174 PT cells); and STZ is insulin-deficient whereas the human samples are obese and insulin-resistant, so the lipogenic arm inverts (Srebf1, Fasn, Fbp1 fall in STZ cholangiocytes).

**Reported axis (candidate-anchored, script 23): AMPKα2 → PGC-1α with CCND1 as output.** Both anchors are Fig 1 candidates. In liver cholangiocytes Ppargc1a↓ (cell BH 2.0e-8; pseudobulk P 0.016) and Ccnd1↑ (BH 0.0066; P 6.4e-4); in kidney endothelium Ccnd1↑ (BH 8.2e-17; P 1.2e-6) and the human-up candidate score rises specifically in that cell type (P 0.0018, BH 0.011). The pair is **not** co-localized or correlated (liver R 0.03, empirical P 0.60; kidney 0.17 % co-expression), consistent with AMPK being controlled post-translationally (Fig 3D–G).

**The downstream partner that does pass all three criteria: Cdkn1a (p21) – Ccnd1 (cyclin D1), in kidney endothelium** (script 18, Fig 3H–J):

| Criterion | Kidney endothelium (1,112 / 1,763 cells) | Liver cholangiocytes (164 / 472) |
|---|---|---|
| Changed, human direction, both tests | **5 / 5**: Phlda3 (pseudobulk P 1.1e-10), Ccnd1 (1.2e-6), Zmat3 (1.6e-6), Cdkn1a (1.4e-5), Bax (4.7e-4) | 4 / 6: Phlda3 (2.0e-6), Cdkn1a (6.4e-6), Ccnd1 (6.4e-4), Bax (7.0e-4) |
| p53-target UCell score | up only here (P 1.6e-4, BH 9.4e-4) — cell-type-specific | P 0.0055 (LSEC 0.010; BH 0.055) |
| Co-localized | 11.1 % of cells co-express the pair; not the atlas maximum (cDC 48 %, proliferating 24 %, urothelium 16 %, PT 12 %, all on 43–152 cells) | 18.7 % |
| Correlated | R 0.07 (P 1.2e-4); co-expressing cells ρ 0.21 (P 1.8e-4); per-mouse r 0.77 (P 0.045); **empirical P 0.009** vs detection-matched random pairs | **fails**: R −0.10, empirical P 0.99 |

Serpine1, Mdm2 and Trp53 are unchanged in mouse, as expected for a program controlled by p53 protein stabilization. The per-mouse r = 0.77 is carried by the group split (both genes rise with STZ; within-group slopes are flat or negative), so the co-expressing-cell ρ = 0.21 is the closer equivalent of the paper's Fig 3G.

The AMPK Fig 3 panels (Sirt1–Ppargc1a, R 0.03, not correlated) are kept as supplementary in `results/supplementary_figures/pathway_selection/AMPK_axis/`; the earlier 19-gene candidate panels are in `results/supplementary_figures/candidates_19gene/`.
