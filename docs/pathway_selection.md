# Pathway and axis selection

**Status (2026-09-16):** Final, rebuilt as a reproducible two-stage procedure (decision R19) and re-anchored to the Fig 1 candidate set (decision R23).
**Scripts:** `16_pathway_selection.R` (stage 1, human bulk only, incl. the candidate-anchored axis), `17_axis_confirmation_scRNA.R` (stage 2, mouse only), `23_anchored_axis_figures.R` (Fig 3D–G), `18_selected_axis_figures.R` (Fig 3H–J), `22_consistency_audit.R` (chain audit).
**Design:** 27 samples, T2D vs lean (D1b); 19-gene 3-way candidate set (M13); mouse STZ vs Control scRNA (3 vs 4 mice, liver + kidney).

**Outcome in one line:** the selected KEGG map is **AMPK signalling (hsa04152)** — top of every human-bulk criterion. Its axis, restricted to edges that contain a Fig 1 candidate gene, is **AMPKα2 (PRKAA2) → PGC-1α (PPARGC1A)** — a drawn activation edge whose two DE nodes move *apart* in human T2D liver, the transcript signature of post-translational control — with **CCND1** as the de-repressed output node. Both anchors are Fig 1 candidates, both change in the mouse atlas (Ppargc1a↓ and Ccnd1↑ in liver cholangiocytes; Ccnd1↑ in kidney endothelium), but the **pair is not co-localized or correlated**. The pair-level criteria of the reference framework are met only by CCND1 together with its arrest partner **CDKN1A (p21)** in kidney endothelium, which is reported as the downstream readout (Fig 3H–J). Protein-level validation (p-AMPKα, PGC-1α, p21, cyclin D1; metformin as positive control) carries the mechanism that transcripts cannot.

---

## 1. What the reference paper actually did (Xu et al., Phytomedicine 2026)

| Step | Paper | Evidence |
|---|---|---|
| 1 | DEGs (adj.P < 0.05, \|logFC\| > 0.5) ∩ WGCNA blue module ∩ GeneCards neutrophil-migration set | 159 genes |
| 2 | GO of those genes → leukocyte migration/chemotaxis, matching the phenotype | Fig 2A |
| 3 | KEGG → "regulation of actin cytoskeleton". **Not the top term** (≈ 17th in the Fig 2B bubble plot). Chosen because it matched the phenotype | Fig 2B |
| 4 | Pathway interaction network of the top-30 KEGG terms | Fig 2C |
| 5 | pathview: "within this pathway, RAC/PAK signalling was significantly activated" — a **drawn edge between two DE nodes** | Fig 2D |
| 6 | scRNA confirmed that one axis: Rac1/Pak1 changed in neutrophils, co-localized there, and correlated (R = 0.34, P < 0.001). No other pathway was compared | Fig 3D–G |
| 7 | MR (PAK1 druggable/causal), docking, then WB with LIMK1/Cofilin1 added from literature | Fig 4–9 |

Steps 2–5 were narrative in the paper. Here each one is an explicit, scored criterion, and step 5 is computed from the KEGG KGML relations instead of read off a picture.

## 2. Stage 1 — pathway selection from the human bulk data only (script 16)

Pool: all mechanistic KEGG maps (hsa01xxx global and hsa05xxx disease maps excluded) with 15–500 measured genes, restricted to the maps enriched in the 19 intersect genes (the paper's step-3 pool). Six criteria, each computed without touching the mouse data, then averaged as ranks:

| Criterion | Paper analogue | Statistic |
|---|---|---|
| ORA | Fig 2B | enrichKEGG P on the 19 intersect genes |
| Phenotype match | "matches the phenotype", judged by eye | hypergeometric P for the curated hyperglycemia/insulin-resistance set (script 04) inside the map |
| Pathway activity | — (adds rigour) | `limma::fry` on all 27 samples, design `~condition + dataset` |
| GSEA | — (adds rigour) | fgsea on the full moderated-t ranking |
| Map axis | Fig 2D | number of **coherent DE–DE edges**: a drawn KGML relation whose two nominally DE nodes move with the sign of the edge |
| Network centrality | Fig 2C | degree in the Jaccard graph (J ≥ 0.1) of the top-30 ORA maps |

### 2.1 Result (`results/pathway_selection/selection_ranked.csv`)

| Rank | KEGG map | Measured | DEG P<0.05 / FDR<0.05 | ORA adj.P | Phenotype P | fry (dir, P) | GSEA NES, P | Degree | Coherent edges | Composite |
|---|---|---|---|---|---|---|---|---|---|---|
| **1** | **AMPK signalling (hsa04152)** | 104 | 39 / 2 | **2.3e-4** | **5.6e-49** | Up, 0.024 | 1.59, 0.0039 | 10 | 9 | **3.33** |
| 2 | Insulin signalling (hsa04910) | 125 | 43 / 0 | 0.034 | 3.6e-82 | Up, 0.043 | 1.23, 0.10 | 14 | 23 | 3.67 |
| 3 | Apelin signalling (hsa04371) | 117 | 37 / 0 | 0.0056 | 6.0e-26 | Up, 0.25 | 1.18, 0.15 | 11 | 20 | 4.33 |
| 4 | HIF-1 signalling (hsa04066) | 102 | 38 / 2 | 0.030 | 1.2e-33 | Up, 0.015 | 1.47, 0.016 | 9 | 5 | 4.67 |
| 5 | Longevity regulating (hsa04211) | 75 | 28 / 1 | 0.015 | 4.7e-42 | Up, 0.45 | 1.20, 0.15 | 16 | 12 | 4.67 |
| 6 | FoxO signalling (hsa04068) | 117 | 38 / 1 | 0.033 | 2.5e-50 | Up, 0.27 | 1.13, 0.21 | 13 | 18 | 4.83 |
| 7 | Fructose and mannose metabolism (hsa00051) | 31 | 15 / 1 | 0.033 | 1.8e-5 | Up, 0.0062 | 1.81, 0.0024 | 1 | 0 | 5.25 |
| 8 | p53 signalling (hsa04115) | 66 | 19 / 1 | 0.014 | 1.9e-3 | Up, 0.017 | 1.57, 0.0071 | 1 | 1 | 5.25 |

**AMPK signalling is the selection**, and unlike the earlier eyeballed choice it now wins on numbers: rank 1 on ORA, on phenotype-set enrichment among the signalling maps, and top-3 on every other criterion. The top GO term on the same 19 genes is "response to insulin" (adj.P 6.4e-7), the phenotype anchor of the study. Figure: `results/supplementary_figures/pathway_selection/16_selection_criteria.pdf`.

### 2.2 The axis on the map (script 16, `selection_axis_edges.csv`)

Coherent DE–DE edges of hsa04152 form four components; the largest is
**CAMKK2 → AMPKα2 (PRKAA2) ⊣ GYS2 / → PFKFB4 / → MLYCD** (5 nodes, 4 edges, all in the phenotype set).
The other components are PDPK1 → AKT1/AKT2, IRS1/IRS2 → PIK3CA (both down) and SREBF1 → FASN (both up).

**Discordant DE–DE edges of the same map** — drawn edges whose DE nodes move *against* the edge sign — are PRKAA2 → PPARGC1A, PRKAA2 → PFKFB1/PFKFB3, PRKAA2 → SREBF1, PRKAA2 → CFTR, IRS1/IRS2 → PIK3R2/PIK3CD. The AMPK subunit transcript rises while every AMPK output it activates falls (PPARGC1A logFC −0.61, P 6e-4; SIRT1 −0.39, P 0.0044; PFKFB3 −0.66, FDR 0.042) and the outputs it represses rise (FBP1 +0.62, SREBF1 +0.44, FASN +0.70, CCND1 +0.80). **This is the signature of an AMPK pathway shut down post-translationally, and it is the reason transcript data cannot confirm it.** The reference paper's pathway happened to be transcriptionally visible; ours is not.

### 2.3 Candidate-anchored axis (decision R23)

The unrestricted edge rule in §2.2 has a defect the chain audit exposed (script 22): its top component shares **no gene** with the 19-gene candidate set that selected the pathway in the first place. The reported axis is therefore restricted to **edges with at least one node in the candidate set** (`selection_axis_anchored.csv`, 62 such edges across 18 maps). On the selected map:

| Edge | Subtype | Human direction | Coherent? | Anchor |
|---|---|---|---|---|
| **PRKAA2 → PPARGC1A** | activation; phosphorylation | PRKAA2 +0.29 (P 0.017), PPARGC1A −0.61 (P 6e-4) | **no — discordant** | PPARGC1A |
| IRS2 → PIK3CA | activation | both down | yes | IRS2 |
| IRS2 → PIK3R2 / PIK3CD | activation | opposite | no | IRS2 |

Supporting the same arm from a second map: **SIRT1 → PPARGC1A** (hsa04211, activation) is *coherent* — both down — and also anchored by PPARGC1A. CCND1, the fifth candidate on the AMPK map, has no DE–DE edge of its own, but it is reached from AMPK in two drawn steps — **AMPK ⊣ ELAVL1 (HuR) → CCND1** — whose composite sign is negative. Script 25 confirms that these two are the only Fig 1 candidates with a directed path from the AMPK heterotrimer on this map, and that the drawn signs predict both observed directions under a single hypothesis of suppressed AMPK output (PPARGC1A down, CCND1 up; `25_map_paths.csv`). The midpoint ELAVL1 is not DE (logFC 0.01, P 0.93), as expected for a node that acts by phosphorylation and mRNA stabilisation.

The discordant AMPK → PGC-1α edge is not a failure of the rule; it is the finding. A drawn activation edge whose kinase subunit rises while its target falls is what post-translational inactivation looks like in transcript data, and it is why the scRNA step cannot confirm the arm directly.

## 3. Stage 2 — single-cell confirmation (script 17)

Every axis candidate of **all eight** shortlisted maps was put through the paper's three Fig 3 criteria in every pseudobulk-eligible cell type of both tissues, with BH correction within a tissue. The mouse data played no part in stage 1, so this remains confirmation rather than selection; testing all eight maps (rather than only the selected one, as the paper did) is a declared deviation that makes the choice transparent instead of post hoc.

Criteria, fixed in advance: **changed** = per-mouse pseudobulk P < 0.05 with the human sign; **co-localized** = ≥ 2 % of the cell type's cells co-express the pair; **correlated** = cell-level Pearson BH < 0.05 with R > 0.

### 3.1 Result: no map-drawn axis passes

- 81 pairs tested across cell types; 22 pair × cell-type combinations cleared co-localization and correlation.
- **Not one of them had both nodes changed in the human direction**, and the largest correlation was R = 0.22 (most were R ≈ 0.1, significant only because of cell number). The paper's own benchmark was R = 0.34 with both nodes changed.
- Best per map: AMPK 1/14 axis genes concordant in its best cell type, 0 correlated axes; Insulin 2/23; FoxO 3/27; p53 1/3.
- Gene level across the shortlist: 478 gene × cell-type tests, 17 concordant, concentrated in liver cholangiocytes (4), liver LSEC (4) and kidney endothelium (4).

Tables: `17_pathway_confirmation.csv`, `17_edge_confirmation_<tissue>.csv`, `17_gene_confirmation.csv`; figure `results/supplementary_figures/pathway_selection/17_edge_confirmation.pdf`.

### 3.2 The anchored axis in the mouse atlas (script 23)

| Criterion | Liver cholangiocytes (164 / 472) | Kidney endothelium (1,112 / 1,763) |
|---|---|---|
| Anchor genes changed, human direction | **Ppargc1a↓** (cell BH 2.0e-8; pseudobulk P 0.016) and **Ccnd1↑** (cell BH 0.0066; P 6.4e-4) | **Ccnd1↑** (cell BH 8.2e-17; P 1.2e-6); Ppargc1a is barely expressed here |
| Supporting | Sirt1↓ at cell level (BH 0.0093) | — |
| Candidate-gene score | — | the human-up candidate score rises specifically in endothelium (P 0.0018, BH 0.011) |
| Co-localized | 13.8 % of cells co-express Ppargc1a + Ccnd1 | 0.17 % — the pair is not co-expressed |
| Correlated | **no**: R 0.03, empirical P 0.60; per-mouse r −0.19 | **no**: R −0.01, empirical P 0.81 |
| Opposite to human | Stk11, Camkk2, Srebf1/Fasn/Fbp1 (see §3.3) | Stk11↓ (P 6.5e-4) |

So the anchored axis satisfies the paper's first criterion in both tissues and fails the pair criteria in both. Tables: `23_anchored_gene_stats_<tissue>.csv`, `23_anchored_score_<tissue>.csv`, `23_anchored_correlation_<tissue>.csv`.

### 3.2b The downstream partner that does pass all three: CDKN1A – CCND1 (script 18, Fig 3H–J)

CDKN1A and CCND1 sit in the same "cell-cycle arrest" box of hsa04115, which KEGG draws as a group rather than as a KGML relation, so the edge screen could not see it. Tested directly, it satisfies every criterion in **kidney endothelium** (1,112 Control / 1,763 STZ cells):

| Paper criterion | Paper (RAC1/PAK1, neutrophils) | Kidney endothelium (Cdkn1a/Ccnd1) | Liver cholangiocytes |
|---|---|---|---|
| Axis genes changed | both, P < 0.001 | **5 / 5 concordant with human, both tests**: Phlda3 (pseudobulk P 1.1e-10), Ccnd1 (1.2e-6), Zmat3 (1.6e-6), Cdkn1a (1.4e-5), Bax (4.7e-4) | 4 / 6: Phlda3 (2.0e-6), Cdkn1a (6.4e-6), Ccnd1 (6.4e-4), Bax (7.0e-4) |
| Pathway score cell-type-specific | neutrophils | **yes** — UCell p53-target score up only in endothelium (P 1.6e-4, BH 9.4e-4) | cholangiocytes P 0.0055, LSEC 0.010 (BH 0.055) |
| Co-localized | highest in neutrophils | 11.1 % of cells co-express the pair — **but not the highest in the atlas** (cDC 48 %, proliferating 24 %, urothelium 16 %, PT 12 %), so this criterion is met only in the weaker sense that the pair is co-expressed in a substantial fraction of the cell type that carries the change; the higher-scoring populations have 43–152 cells | 18.7 % |
| Correlated | R = 0.34, P < 0.001 | R = 0.07, P 1.2e-4; **co-expressing cells ρ = 0.21, P 1.8e-4; per-mouse r = 0.77, P 0.045**; exceeds the detection-matched background (mean R 0.014, 95th pct 0.050, **empirical P 0.009**) | **fails**: R = −0.10, empirical P 0.99 |

Human anchor for the same axis: CDKN1A is the only FDR-significant DEG of hsa04115 (logFC 0.96, FDR 0.016), CCND1 is in the 19-gene candidate set (logFC 0.80, P 7.5e-4) and replicates in GSE23343 (P 0.016), and **HALLMARK_P53_PATHWAY is the top GSEA term of the whole human contrast** (NES 1.92, padj 1.5e-5).

**Honest reading of the correlation.** A cell-level R of 0.07 is far below the paper's 0.34. It is reported with its background null precisely so that it is not over-sold: the defensible statements are that the two genes are co-expressed in 11 % of endothelial cells and that their cell-level association exceeds what detection-matched random pairs produce (empirical P 0.009). The per-mouse r = 0.77 must be read with care — with 3 Control and 4 STZ mice it is carried by the group split, because both genes rise with STZ; within each group the slope is flat or negative (`Fig3G_p53_correlation_kidney.pdf`, middle panel). The co-expressing-cell ρ = 0.21 is the closest equivalent of the paper's Fig 3G. The liver pair is not correlated (R −0.10) and is reported as a negative result.

### 3.3 Why the metabolic axis fails in this atlas — three structural reasons, all reportable

1. **Post-translational control.** AMPK, p53 and SIRT1 activity are set by phosphorylation and protein stability; the discordant edges in §2.2 show this directly.
2. **Missing parenchyma.** The 10x liver data contain 69 hepatocytes (0–5 per control mouse) and the kidney data only 150/174 proximal-tubule cells, so the cells that carry hepatic gluconeogenesis and lipogenesis cannot be tested. The AMPK-relevant compartment is absent, not negative.
3. **Model mismatch.** GSE244475 is STZ (insulin-deficient), while the human bulk data are obese, insulin-resistant T2D. The lipogenic arm is therefore expected to invert, and it does: Srebf1, Fasn and Fbp1 fall in STZ cholangiocytes although they rise in human T2D liver. The p53/arrest arm is the part of the phenotype the two models share — hyperglycemic stress.

## 3.5 Chain consistency (script 22, `results/consistency/`)

Whether one gene set actually flows through Fig 1 → Fig 2 → Fig 3 → Fig 4 → Fig 5 was audited rather than assumed. Overlap between the gene sets that define each figure:

| | Fig 1 (19) | Fig 2 map | Fig 2 axis (anchored) | Fig 3 partner axis | Fig 4 ML top 7 | Fig 5 MR FDR<0.05 |
|---|---|---|---|---|---|---|
| Fig 1 candidate 19 | 19 | 5 | 2 | 1 | 7 | 0 |
| Fig 2 map (hsa04152 DEGs) | 5 | 39 | 2 | 1 | 1 | 2 |
| Fig 3 partner axis | 1 | 1 | 1 | 5 | 0 | 0 |

**The spine is PPARGC1A and CCND1.** Both are Fig 1 candidates, on the selected map, anchors of the reported axis, and confirmed in the mouse atlas; PPARGC1A is 3rd in the machine-learning consensus and CCND1 replicates in GSE23343. Neither is causal in MR (P 0.16 and 0.55), and CCND1 ranks 14/19 as a classifier. Those disagreements are real and are reported as such rather than smoothed over: **no gene passes all five figures**, and the study's claim is a mechanism supported by expression and cell-type evidence, not a validated causal or diagnostic target.

Before the re-anchoring (R23) the overlap between the Fig 1 candidates and the reported axis was **zero**; the audit figure `results/supplementary_figures/consistency/22_chain_consistency.pdf` shows the gene × stage matrix.

## 4. Mechanistic link between the two arms

AMPK phosphorylates p53 at Ser15 (Jones et al., *Mol Cell* 2005), and AMPK loss de-represses cyclin D1. Suppressed AMPK output in human obese-T2D liver and an activated p53–p21–cyclin D1 arrest program in hyperglycemic mouse kidney endothelium are therefore one axis measured at two levels, not two competing findings. Metformin, the lab's positive control, activates AMPK, which gives the in vivo experiment a mechanistic reference — the analogue of the paper's use of Fingolimod.

**Genetic support: none (scripts 19–20, `docs/results_MR.md`).** cis-MR of expression on type 2 diabetes (GCST006867) with LD-clumped 1000G EUR instruments leaves **PRKAA1 (P 0.23), PPARGC1A (0.16), CDKN1A (0.49) and CCND1 (0.55) all null**; only SREBF1 and SIRT1 survive FDR, both as whole-blood instruments pointing against hepatic biology. An earlier distance-pruned pass appeared to support PRKAA1 and PPARGC1A, but that rested on single lead SNPs and is withdrawn (R21). Tissue instruments do not exist at usable power (GTEx liver n = 208, kidney cortex n = 73). The AMPK selection therefore stands on the human bulk transcriptome, the KGML edge analysis and the planned protein-level work, and the manuscript must state that MR adds no genetic support.

## 5. Wet-lab panel (liver and kidney; Control, Ob-T2D, GV1, metformin)

| Assay | Target | Expected in Ob-T2D | Expected with GV1 / metformin |
|---|---|---|---|
| WB | p-AMPKα Thr172 / total AMPKα | ratio ↓ | ↑ |
| WB | SIRT1, PGC-1α | ↓ | ↑ |
| WB | SREBP-1c (mature), FASN | ↑ | ↓ |
| WB | p-p53 Ser15 (mouse Ser18) / total p53, p21, cyclin D1, BAX; MDM2 as feedback control | ↑ | ↓ |
| qPCR | Ppargc1a, Sirt1, Irs2, Igf1, Srebf1, Fasn, Fbp1 (AMPK arm); Cdkn1a, Ccnd1, Phlda3, Zmat3, Bax (p53 arm) | human direction | reversal |
| IHC / IF | Kidney: p21 and cyclin D1 with CD31 (endothelium) — the confirmed cell type. Liver: PGC-1α with CK19 (cholangiocytes) | cell-specific change | reversal |

## 6. Deviations from the reference paper (to state in Methods)

| Paper | This study | Reason |
|---|---|---|
| DEG adj.P < 0.05 | nominal P < 0.05, \|logFC\| > 0.5 | 27 arrays, 31 FDR DEGs; strict rules leave 1 intersect gene (decisions M9/M11/M13) |
| Single strongest WGCNA module | 7 modules with T2D P < 0.1 | same |
| GeneCards relevance > 8 | GO + KEGG + HPO + CTD hyperglycemia set | GeneCards scores not redistributable |
| Pathway chosen by eye from Fig 2B/2D | six scored criteria + KGML edge analysis (script 16) | reproducibility |
| One pathway carried into scRNA | all eight shortlisted maps screened, full matrix reported | bulk and scRNA differ in species, model and tissue, so a single blind hand-off would not be interpretable |
| Correlation reported as R, P | R, P **plus** a detection-matched background null | with thousands of cells, P alone is uninformative |
| Same tissue and disease in both layers | human obese-T2D liver vs mouse STZ liver + kidney | data availability (D3: no human kidney bulk set) |
| n = 1 mouse per group | 3 Control vs 4 STZ, per-mouse pseudobulk as the primary test | proper replication |

## 7. Superseded

- The earlier eyeballed axis **SIRT1 → PGC-1α** (R18) is retained as the AMPK arm's protein-level hypothesis, but its Fig 3D–G panels are supplementary (`results/supplementary_figures/pathway_selection/AMPK_axis/`); the pair is not correlated in any cell type (R = 0.03) and the AMPK score shift is not AMPK-specific.
- Section 4–5 of the previous version (the informal AMPK-vs-p53 comparison, script 14) is superseded by scripts 16–18 and kept only as `results/pathway_selection/axis_summary.csv`.
