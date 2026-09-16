# Bulk Liver Results: T2D vs Lean Control (GSE15653 + GSE64998, 27 samples)

**Status (2026-09-14):** Sections 1–4 are final. Sections 5–7 are pending until scripts 05, 11 and 10 finish.
**Scripts:** `scripts/01`–`05`, `11`. All decisions are in `docs/decisions_log.md`.
**Earlier designs (archived):** T2D vs all non-diabetic, 39 samples → `archive/superseded_T2D_vs_allNonDiabetic_39samples/`.

| Decision | Setting |
|---|---|
| D1b | **T2D (n = 16) vs lean / non-obese Control (n = 11)**; obese non-diabetic samples excluded |
| M8 | Raw CEL + RMA per platform on the selected samples only; **no expression filter** |
| M9 | DEG = nominal **P < 0.05 and \|log2FC\| > 0.5** (FDR also reported) |
| M11 | WGCNA key modules = all modules correlated with T2D at **p < 0.1** (not dataset-driven) |
| M13 | Candidates = **paper 3-way overlap**: DEGs ∩ key-module genes ∩ hyperglycemia gene set; hub status is annotation |

In this design every T2D sample is obese and every control is lean, so obesity and diabetes effects cannot be separated. Results describe **obese T2D vs lean** liver.

---

## 1. Samples and preprocessing (Fig 1A)

| Dataset | Control | T2D |
|---|---|---|
| GSE15653 (HG-U133A) | GSM391693–GSM391697 (lean, n = 5) | GSM391702–GSM391710 (obese DM well + poorly controlled, n = 9) |
| GSE64998 (HuGene 1.1 ST) | GSM1585592–GSM1585597 (non-obese, n = 6) | GSM1585585–GSM1585591 (obese T2D, n = 7) |

- **Genes:** 12,650 (GSE15653) and 18,865 (GSE64998) annotated → **11,972 common genes, all kept**.
- **Outliers:** none removed. GSM391703 and GSM391704 were each flagged by connectivity only (Z.k < −2.5); removal needs ≥ 2 flags.
- **Batch correction:** ComBat (batch = dataset, condition protected) for PCA and WGCNA; DE uses dataset as a covariate.

## 2. Differential expression (Fig 1B)

Model: limma `~condition + dataset`.

| Analysis | Up | Down |
|---|---|---|
| **FDR < 0.05 & \|log2FC\| > 0.5 (reference-paper criteria)** | **12** | **19** |
| Paper-style ComBat → limma, FDR < 0.05 | 34 | 25 |
| + sex covariate, FDR < 0.05 | 0 | 0 |
| GSE15653 only / GSE64998 only, FDR < 0.05 | 0 / 4 | 0 / 4 |
| **Nominal P < 0.05 & \|log2FC\| > 0.5 (used downstream)** | **96** | **82** |

- 84 genes reach FDR < 0.05 regardless of fold change.
- **Robustness of the 178 nominal DEGs:**
  - 97.8% (GSE15653) and 98.3% (GSE64998) have the same direction in each dataset alone.
  - 93% stay nominal after adding sex.
  - Jaccard overlap with paper-style ComBat → limma is 0.67.
- The sex-adjusted model has no FDR genes because sex is partly confounded: all GSE15653 lean controls are female.
- **Top DEGs:** HBB↑ (blood contamination), VSNL1↓, MYOT↓, CDKN1A↑, HPS5↑, DAPK1↓, KMO↓, LZTFL1↓, GAS6↑, CA14↓, TUT7↓, GPR88↓, CYP2C19↓, SREBF2↑, HSPA5↑.
- **All 10 genes of the earlier project are DEGs here** (INHBE, SERPINE1, SLC16A4, P4HA1, PPARGC1A, FASN, PLIN1, IGFBP1, TUT7/ZCCHC6, SLCO1A2). The earlier 10-gene result therefore depended mainly on using lean controls (decision R12).

## 3. WGCNA (Fig 1C–F)

- **Input:** top 5,000 MAD genes (pre-specified for < 30 samples).
- **Network:** signed-hybrid, bicor; **soft power 7** (scale-free R² = 0.82).

**Key modules (T2D p < 0.1; all have p ≤ 0.016):**

| Module | Genes | r (T2D) | p | r (HbA1c, GSE15653) | Preservation Zsummary |
|---|---|---|---|---|---|
| red | 135 | **+0.65** | 2×10⁻⁴ | +0.50 | 24.1 |
| yellow | 170 | **−0.60** | 0.001 | −0.21 | 25.9 |
| tan | 46 | −0.60 | 0.001 | −0.37 | 13.6 |
| greenyellow | 67 | −0.52 | 0.006 | −0.55 (p = 0.04) | 11.2 |
| blue | 1,671 | −0.48 | 0.012 | −0.57 (p = 0.035) | 8.2 (moderate) |
| green | 137 | +0.47 | 0.014 | +0.53 (p = 0.05) | 20.3 |
| turquoise | 1,796 | +0.46 | 0.016 | +0.51 | 8.2 (moderate) |

- None of these modules is dataset-driven (\|r_Dataset\| ≤ 0.08).
- **Hub genes:** 1,770 (\|kME\| > 0.7, \|GS\| > 0.2).
- **Fig 1F** shows MM–GS panels for the six modules that contain candidate genes (red, yellow, tan, blue, green, turquoise).

## 4. Candidate genes: paper 3-way overlap (Fig 1G)

**Venn counts:**
- 178 DEGs
- 4,022 key-module genes
- 618 hyperglycemia genes present in the data
- **19 genes in the triple overlap** (DEG ∩ module only: 156; DEG ∩ gene set only: 1)

| Gene | Direction | logFC | P | FDR | Module | Hub | Gene-set source(s) |
|---|---|---|---|---|---|---|---|
| VSNL1 | ↓ | −0.57 | 3×10⁻⁶ | **0.016** | yellow | yes | GO response to monosaccharide / carbohydrate |
| LZTFL1 | ↓ | −0.53 | 1×10⁻⁵ | **0.034** | blue | yes | HPO insulin resistance |
| GAS6 | ↑ | +0.57 | 1×10⁻⁵ | **0.037** | turquoise | no | GO response to monosaccharide / carbohydrate |
| SREBF2 | ↑ | +0.51 | 2×10⁻⁴ | **0.047** | red | yes | CTD insulin resistance |
| SERPINE1 | ↑ | +0.94 | 5×10⁻⁴ | 0.053 | turquoise | no | KEGG AGE-RAGE |
| PPARGC1A | ↓ | −0.61 | 6×10⁻⁴ | 0.054 | yellow | no | KEGG insulin resistance |
| HMOX1 | ↑ | +0.53 | 6×10⁻⁴ | 0.054 | green | yes | CTD insulin resistance |
| CCND1 | ↑ | +0.80 | 7×10⁻⁴ | 0.058 | red | yes | KEGG AGE-RAGE |
| ENPP1 | ↓ | −0.52 | 9×10⁻⁴ | 0.060 | blue | yes | GO response to insulin; HPO hyperglycemia / insulin resistance |
| PNPLA3 | ↑ | +0.63 | 9×10⁻⁴ | 0.060 | turquoise | yes | GO response to insulin |
| GPX1 | ↑ | +0.69 | 0.001 | 0.064 | green | yes | CTD hyperglycemia; GO carbohydrate |
| IGFBP1 | ↓ | −1.66 | 0.002 | 0.067 | yellow | yes | GO response to insulin |
| FBP1 | ↑ | +0.62 | 0.003 | 0.075 | turquoise | yes | GO carbohydrate; response to insulin |
| IRS2 | ↓ | −0.58 | 0.003 | 0.075 | tan | yes | CTD, GO, HPO, KEGG insulin resistance / T2D (8 sources) |
| PLIN1 | ↑ | +0.53 | 0.006 | 0.092 | turquoise | no | HPO insulin resistance |
| CTSD | ↑ | +0.55 | 0.007 | 0.096 | turquoise | yes | GO response to insulin |
| KHK | ↑ | +0.53 | 0.011 | 0.108 | turquoise | yes | GO response to insulin; HPO hyperglycemia |
| IGF1 | ↓ | −0.75 | 0.016 | 0.122 | yellow | yes | GO carbohydrate / insulin; HPO insulin resistance |
| IGFBP2 | ↓ | −0.86 | 0.030 | 0.154 | yellow | no | CTD insulin resistance |

**Summary of the 19 genes:**
- **Direction:** all 19 change the same way in both discovery datasets.
- **Network:** 14 of 19 are WGCNA hubs; 4 reach FDR < 0.05.
- **Overlap with earlier sets:** 4 of the earlier project's 10 genes are included (SERPINE1, PPARGC1A, IGFBP1, PLIN1), and CCND1 from the interim set.
- **CDKN1A is absent:** it is a strong DEG (FDR 0.016) but is not in the hyperglycemia gene set.
- **Biology:** the set centres on insulin / IGF signalling (IRS2, IGF1, IGFBP1, IGFBP2, ENPP1, PPARGC1A), hepatic lipid metabolism (SREBF2, PNPLA3, PLIN1), gluconeogenesis and fructose metabolism (FBP1, KHK), and oxidative stress / fibrinolysis (HMOX1, GPX1, SERPINE1).
- **Caveat:** the 7 key modules contain 4,022 of the 5,000 WGCNA genes, so the module filter barely narrows the list (175 of 178 DEGs pass). The 3-way set is effectively DEG ∩ hyperglycemia gene set, and hub status is the stricter network evidence.

## 5. Functional enrichment (Fig 2)

**Over-representation analysis of the 19 candidate genes** (clusterProfiler; universe = 11,972 genes; BH adj.P < 0.05).

| Collection | Significant | Top terms (genes) |
|---|---|---|
| GO (Fig 2A) | 230 terms (BP 221, CC 5, MF 4) | Response to insulin (IRS2, IGF1, IGFBP1, ENPP1, CTSD, KHK, FBP1, PNPLA3; adj.P 6×10⁻⁷); cellular response to insulin stimulus; regulation of carbohydrate / glycogen metabolism; IGF receptor signalling (IRS2, IGF1, IGFBP1, IGFBP2); insulin receptor signalling; response to nutrient levels |
| KEGG (Fig 2B–C) | 9 pathways | **AMPK signalling** (IRS2, IGF1, PPARGC1A, CCND1, FBP1; adj.P 2×10⁻⁴); Apelin signalling; **p53 signalling** (IGF1, CCND1, SERPINE1); longevity regulating; HIF-1; FoxO; fructose and mannose metabolism (KHK, FBP1); insulin signalling |

**Fig 2D:** pathview map of **AMPK signalling (hsa04152)**, coloured by T2D vs Control logFC (`figures/Fig2_enrichment/Fig2D_pathview/hsa04152.T2D_vs_Control_logFC.png`). This is the analogue of the reference paper's pathway diagram.

**GSEA on the full T2D vs Control ranking** (moderated t):

| Collection | Significant (padj < 0.05) | Top sets (NES) |
|---|---|---|
| Hallmark | 22 (19 up, 3 down) | Up: cholesterol homeostasis 2.21, **p53 pathway 1.92**, MYC targets v2 1.80, apical junction, **unfolded protein response 1.77**, IL6–JAK–STAT3 1.75, UV response up, apoptosis, myogenesis, complement. Down: UV response down −1.62, protein secretion −1.61 |
| Reactome | 32 (29 up, 3 down) | Up: **IRE1α activates chaperones 2.60**, UPR 2.10, rRNA expression regulation, DNA methylation, histone arginine methylation, lipoprotein clearance, COPI anterograde transport. Down: snRNP assembly, nuclear envelope breakdown |
| KEGG (MEDICUS) | 0 | — |

**Biological themes carried to the wet-lab rationale:**
1. Insulin / IGF signalling and hepatic glucose metabolism: IRS2↓, IGF1↓, IGFBP1/2↓, ENPP1↓, FBP1↑, KHK↑.
2. AMPK / PGC-1α energy sensing: PPARGC1A↓.
3. p53 / cell-cycle control: CCND1↑; CDKN1A↑ as the top non-3-way DEG.
4. ER stress / UPR, IRE1α.
5. Cholesterol and lipid homeostasis: SREBF2↑, PNPLA3↑.

**Pathway and axis selection** (paper logic: enrichment → pathview → scRNA co-localization and correlation) is in `docs/pathway_selection.md` (decision R16). The recommendation is **p53 signalling with a CDKN1A–CCND1 axis**. AMPK ranks first by ORA but fails the scRNA step.

## 6. External validation: GSE23343 (10 T2D vs 7 NGT; raw CEL + RMA)

| Test | Result |
|---|---|
| **CCND1** | Concordant ↑, logFC 0.54, **P = 0.016**, AUC 0.81 → **replicated** |
| **PNPLA3** | Concordant ↑, logFC 0.46, **P = 0.017**, AUC 0.74 → **replicated** |
| IGFBP2 | Concordant ↓, logFC −0.90, P = 0.061, AUC 0.74 → borderline |
| HMOX1, GPX1, SERPINE1, PLIN1, GAS6, SREBF2, IGFBP1, PPARGC1A | Concordant, P > 0.1 |
| VSNL1, KHK, LZTFL1, CTSD, IGF1, FBP1, IRS2, ENPP1 | Discordant |
| 19 genes, direction | 11 / 19 concordant (sign test p = 0.32); 19-gene signature AUC 0.61 (95% CI 0.31–0.92) |
| 178 DEGs, direction | 76 / 178 concordant (p = 0.98); signature AUC 0.43 |
| Genome-wide logFC correlation with discovery | ρ = −0.29 |

The GSE23343 cohort (Japanese hospital cohort, percutaneous needle biopsies; BMI not given on GEO) does not reproduce the discovery profile as a whole (see R6). Two candidates replicate individually: CCND1 and PNPLA3.

## 7. Mouse STZ scRNA evidence (details in `docs/results_scRNA.md` §6–9)

| Candidate | Mouse evidence | Overall |
|---|---|---|
| **CCND1** | Kidney endothelium Tier 1 (FDR 0.002); liver cholangiocytes nominal (FDR 0.12) | **Human DEG + hub + replication + mouse** → top score in both organs (liver 71%, kidney 82%) |
| PPARGC1A, GPX1 | Liver cholangiocytes nominal, concordant | Supportive |
| IRS2 | Liver pDC nominal | Weak |
| PNPLA3, IGFBP2, FBP1, IGF1, ENPP1, GAS6, SERPINE1 | Not DE in mouse cell types (hepatocyte genes poorly captured) | Human evidence only; test by whole-liver qPCR |
| VSNL1, PLIN1 | No mouse ortholog detected | Not testable in mouse |
| *CDKN1A* (not in the 3-way set) | Cdkn1a up in liver cholangiocytes (FDR 0.009) and kidney endothelium (FDR 0.008) | Second-ranked through the scRNA route (67% in both organs) |

## 8. Limitations

- **Obesity and diabetes are fully collinear** in this design (obese T2D vs lean).
- **Sex is partly confounded:** GSE15653 lean controls are all female, and the sex-adjusted model has no FDR genes.
- **Modest cohort:** 27 samples on two array platforms, with a between-dataset logFC correlation across all genes of only 0.25.
- **Module filter:** see the caveat in Section 4; the three-way set is effectively DEG ∩ gene set.
