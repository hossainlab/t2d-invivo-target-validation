# Potential targets, how they were selected, and why the method is defensible

**Status (2026-09-16).** Script: `24_target_scoring.R`. Tables: `results/targets/target_tiers.csv`, `target_evidence_long.csv`.
**Related:** `docs/pathway_selection.md` (pathway + axis), `docs/results_scRNA.md` (cell-type confirmation), `docs/results_ML.md` (classification), `docs/results_MR.md` (causal genetics), `results/consistency/` (chain audit, script 22).

---

## 1. Short answer

| Tier | Genes | What they are | What to do with them in the GV1 Ob-T2D experiment |
|---|---|---|---|
| **Tier 1 — primary targets** | **PPARGC1A (PGC-1α)**, **CCND1 (cyclin D1)** | Human-DE, in the Fig 1 candidate set, on the selected AMPK map, anchors of the reported axis, and confirmed in the mouse atlas in a named cell type | Full panel: qPCR + Western blot in liver and kidney; IF with a cell-type marker (PGC-1α/CK19; cyclin D1/CD31) |
| **Tier 2 — readouts, not targets** | **PFKFB3**, IRS2, PPP2R5A, and the arrest program CDKN1A (p21), ZMAT3, BAX, PHLDA3 | Human-DE and confirmed in the mouse atlas, but each missing one Tier 1 anchor (PFKFB3, CDKN1A: outside the Fig 1 candidate set; IRS2, PPP2R5A: confirmed only in a minor cell type) | qPCR + WB: PFKFB3 (glycolysis arm), p21 and BAX (arrest arm); p21/CD31 IF in kidney |
| **Tier 2b — anchored but untestable here** | IGF1, FBP1 | Fig 1 candidates on the selected map whose cell of origin (hepatocyte, proximal tubule) is missing from the scRNA atlas | Whole-tissue qPCR only; do not claim cell-type specificity |
| **Tier 3 — protein-level nodes** | PRKAA1/PRKAA2 (AMPKα), SIRT1, STK11, CAMKK2 | The mechanism itself. Activity is set by phosphorylation/NAD⁺, which transcript data cannot see | WB only: p-AMPKα Thr172 / total AMPKα, SIRT1; metformin arm as the pharmacological positive control |
| **Not targets** | VSNL1, LZTFL1, SREBF2, GAS6, HMOX1, SERPINE1, PNPLA3, … and the MR hits SREBF1, SIRT1 (blood) | Supported on one evidence axis only | Report; do not carry into validation |

**The two genes to put on the grant figure are PPARGC1A and CCND1.** They are the only ones that survive every stage that was actually able to test them.

---

## 2. Why these two

### PPARGC1A (PGC-1α) — the effector the pathway points at

| Evidence | Value |
|---|---|
| Human obese T2D vs lean liver | logFC −0.61, P 6.0e-4 (FDR 0.054) |
| WGCNA | member of a T2D-associated module; in the hyperglycemia gene set |
| Fig 1 candidate set | yes (1 of 19) |
| Selected KEGG map | yes, hsa04152, **one drawn step from the AMPK heterotrimer**: AMPK →(activation; phosphorylation) PPARGC1A; anchor of the reported axis |
| Second map, coherent | SIRT1 → PPARGC1A on hsa04211, both down |
| Mouse scRNA | **Ppargc1a↓ in liver cholangiocytes**, cell-level BH 2.0e-8, pseudobulk P 0.016 |
| Machine learning | **rank 3** of 19 in the consensus (LASSO selection frequency 0.89) |
| Independent human cohort | same direction, not significant (P 0.63) |
| cis-MR on T2D | null (OR 0.98, P 0.16) |

PGC-1α is the transcriptional co-activator through which AMPK controls mitochondrial biogenesis and gluconeogenic tone. It is down in human T2D liver, down in the mouse, ranks third as a classifier, and sits at the receiving end of the one axis edge that both starts inside our candidate set and is drawn on the selected map. Its weakness is genetic: MR gives no causal support, so it is presented as a mechanistic node, not a validated causal target.

### CCND1 (cyclin D1) — the de-repressed output and the most reproducible readout

| Evidence | Value |
|---|---|
| Human obese T2D vs lean liver | logFC +0.80, P 7.5e-4 |
| Fig 1 candidate set | yes; hub gene in the red module |
| Selected KEGG map | yes, **two drawn steps from AMPK**: AMPK ⊣(inhibition) ELAVL1/HuR →(activation) CCND1; net sign negative |
| Mouse scRNA | **Ccnd1↑ in kidney endothelium**, pseudobulk P 1.2e-6, cell-level BH 8.2e-17; also ↑ in liver cholangiocytes (P 6.4e-4) |
| Independent human cohort GSE23343 | **replicated**, P 0.016, AUC 0.81 |
| Machine learning | rank 14 of 19 — weak as a classifier |
| cis-MR on T2D | null (OR 0.96, P 0.55) |

CCND1 is the only candidate that reproduces in a second human cohort *and* in a second species *and* in a defined cell type. Its link to AMPK is drawn on the map itself — AMPK inhibits HuR (ELAVL1), HuR stabilises cyclin D1 mRNA — so loss of AMPK output predicts cyclin D1 up. The intermediate, ELAVL1, is not differentially expressed (logFC 0.01, P 0.93), which is expected for a node that acts through phosphorylation and mRNA stabilisation rather than through its own transcript level; the two-step route is therefore drawn and sign-consistent, but not DE-supported at its midpoint. Because MR is null and its classifier weight is low, it is proposed as a **pharmacodynamic readout of the axis**, not as a druggable target.

### Are the two consistent with each other? (script 25, `25_map_paths.csv`)

Yes on mechanism and direction, no at the level of cells. Both statements matter.

**Consistent — one hypothesis explains both directions.** Of the 19 Fig 1 candidates, only these two have a *directed path from the AMPK heterotrimer* drawn on hsa04152. Multiplying the edge signs along each path and assuming a single upstream cause — AMPK output suppressed in obese T2D liver — predicts each gene's direction:

| Gene | Drawn path from AMPK | Edge signs | Net | Predicted if AMPK output is lost | Observed in human | Match |
|---|---|---|---|---|---|---|
| PPARGC1A | AMPK → PPARGC1A | + | +1 | **down** | −0.61 (P 6.0e-4) | ✔ |
| CCND1 | AMPK ⊣ ELAVL1 → CCND1 | − then + | −1 | **up** | +0.80 (P 7.5e-4) | ✔ |

Two genes, opposite directions, one cause, both predicted correctly by the map before the data were consulted. That is the strongest internal consistency the study has, and it is why the axis is stated as *AMPK output suppression* rather than as a claim about AMPK transcript levels (PRKAA2 mRNA is in fact slightly **up**, +0.29, P 0.017 — the discordance that marks post-translational control).

**Not consistent — they are not the same cells.** In the mouse atlas the two anchors are confirmed in different compartments: Ppargc1a↓ in **liver cholangiocytes** (PGC-1α is essentially restricted to them in this non-parenchymal atlas), Ccnd1↑ mainly in **kidney endothelium**. Accordingly the pair is not co-expressed or correlated anywhere: 13.8 % co-expression and R 0.03 (empirical P 0.60) in cholangiocytes, 0.17 % co-expression and R −0.01 in endothelium. They are two readouts of one upstream event in two tissues, **not** a co-regulated module, and the manuscript must not draw them as an arrow between two cell-level measurements.

**Not consistent — the downstream evidence disagrees.** PPARGC1A is ML rank 3 but does not replicate in GSE23343 (P 0.63); CCND1 replicates (P 0.016) but is ML rank 14; both are null in MR (P 0.16 and 0.55). Each is strong on the axes the other is weak on, which is why both are carried and why neither is called a validated target.

### Why CDKN1A (p21) is Tier 2 and not Tier 1

CDKN1A is the strongest single gene in the whole study on several axes — the only FDR-significant human DEG in the p53 map (logFC 0.96, FDR 0.016), up in mouse cholangiocytes (P 6.4e-6) and kidney endothelium, and together with CCND1 it forms the only pair that satisfies all three of the reference framework's single-cell criteria. It is **not** in the Fig 1 candidate set for one reason: it is absent from the curated hyperglycemia gene set, although it is a DEG and a member of a T2D-associated module. Promoting it would mean changing the Fig 1 selection rule after seeing the downstream result, which is exactly the circularity this project has been avoiding (decision R17). It is therefore reported as the arrest readout downstream of the axis, with its provenance stated.

## 2.4 Why only two, when the map has 122 genes

The map is the starting pool, not the answer. Each filter below is one of the five in §3, applied to the map itself (script 25, `25_map_paths.csv`). Two DE definitions are shown because they answer different questions — the Fig 1 rule uses both a P-value and an effect-size cut, while the pathway-level statistics use the P-value alone:

| Step | DE = P < 0.05 **and** \|logFC\| > 0.5 (Fig 1 rule) | DE = P < 0.05 only |
|---|---|---|
| 1. On KEGG hsa04152 | **122** | 122 |
| 2. Measured on both array platforms | 104 | 104 |
| 3. Differentially expressed | **8** | 39 |
| 4. + reachable from the AMPK node by a *drawn directed path* | **4** | 14 |
| 5. + the drawn sign predicts the observed direction | **4** | 9 |
| 6. + confirmed in the mouse atlas | **3** | 3 |
| 7. + inside the Fig 1 candidate set | **2** | 2 |

Step 4 is what does most of the work, and it is the step most pathway papers skip: being *annotated to* AMPK signalling is not the same as being *downstream of* AMPK. Of the 122 map genes only 49 are reachable from the AMPK node at all; IRS2 and IGF1, for instance, sit on the map as insulin-branch **inputs**, not as AMPK outputs.

Under the looser P-only definition the nine genes surviving step 5 are PFKFB3, PPARGC1A, CFTR, PFKFB1, SREBF1 (1 step) and CCND1, EEF2, FASN, RPS6KB2 (2 steps) — the map's own prediction of what suppressed AMPK output should look like. Five of the nine are lipogenic or glycolytic nodes that the STZ mouse cannot test, because it is insulin-deficient rather than insulin-resistant (`docs/pathway_selection.md` §3.3).

**The last step is a judgement call, and it costs one good gene.** PFKFB3 clears steps 1–6: FDR 0.042 (one of only two FDR-significant genes on the whole map), one drawn step from AMPK, sign-consistent, and down in mouse liver LSEC (pseudobulk P 0.049) and kidney LoH (cell-level BH 0.03). It is excluded from Tier 1 only because it is absent from the curated hyperglycemia gene set — the same reason CDKN1A is excluded, and the same rule that was fixed before any of this was run. **It is therefore Tier 2 and should be measured**, with its provenance stated: it entered on pathway and cross-species evidence, not through the Fig 1 funnel.

ULK1, the other FDR-significant map gene (+0.47, FDR 0.042, hub of a T2D module), fails step 5: AMPK activates ULK1, so under AMPK suppression it should fall, and it rises. That is a genuine disagreement with the model and is reported rather than dropped.

## 2.5 Is the funnel's output better than chance? (script 26)

The disease labels were permuted 1,000 times **within dataset** (preserving the batch structure) and the entire funnel re-run on each permutation. Observed counts against that null:

| Step | Observed | Null mean | Null 95th pct | Empirical P |
|---|---|---|---|---|
| Fig 1 candidate set | **19** | 1.0 | 4 | **0.002** |
| 3. DE on the AMPK map | **8** | 0.46 | 2 | **0.002** |
| 4. + reachable from AMPK | **4** | 0.27 | 1 | **0.004** |
| 5. + sign matches the map | **4** | 0.12 | 1 | **0.001** |
| 6. + confirmed in the mouse atlas | **3** | 0.05 | 0 | **0.001** |
| 7. + Fig 1 candidate (Tier 1) | **2** | 0.04 | 0 | **0.003** |

Every stage is far outside its null: a permuted dataset produces a Tier 1 list of 2 genes about 4 times in 1,000. The empirical false-discovery rate for the statement *"two genes survive the whole funnel"* is therefore ≈ 0.003, and the 19-gene candidate set itself is not a threshold artefact (permuted sets average 1 gene).

**What the permutation does not establish.** The *fraction* of reachable DE genes whose direction the map predicts is 4/4 under the Fig 1 rule, but with only four genes the null reaches 1.0 often enough that this is not significant on its own (empirical P 0.088; under the looser P-only rule the fraction is 9/14 = 0.64). The map's directional agreement is supportive, not proven. Figure: `results/supplementary_figures/consistency/26_funnel_null.pdf`.

## 2.6 How much does this depend on the DEG threshold? (script 26)

This is the study's most consequential arbitrary choice, and it is now tabulated (`26_threshold_sensitivity.csv`):

| DEG rule | Modules | Candidates | On the AMPK map | Reachable from AMPK | Tier 1/2 genes kept |
|---|---|---|---|---|---|
| **P < 0.05 & \|logFC\| > 0.5 (used)** | 7 key | **19** | 5 | 2 | PPARGC1A, CCND1, IRS2, VSNL1, LZTFL1, GAS6, SREBF2 |
| P < 0.01 & \|logFC\| > 0.5 | 7 key | 16 | 4 | 2 | same |
| FDR < 0.10 & \|logFC\| > 0.5 | 7 key | 16 | 4 | 2 | same |
| **FDR < 0.05 & \|logFC\| > 0.5** | 7 key | **4** | **0** | **0** | VSNL1, LZTFL1, GAS6, SREBF2 only |
| FDR < 0.05, no logFC cut | 7 key | 6 | 0 | 0 | VSNL1, LZTFL1, GAS6, SREBF2 |
| Paper-strict (FDR + strongest module) | red | 1 | 0 | 0 | SREBF2 |

**Read this honestly.** The AMPK/PGC-1α/cyclin D1 result is stable from nominal P < 0.05 through FDR < 0.10, which spans the thresholds most papers use. It disappears entirely at FDR < 0.05: PPARGC1A (FDR 0.054) and CCND1 (0.058) fall just outside, and with them every AMPK-map candidate — the strict candidate set contains no gene on the map at all.

The four genes that do survive FDR < 0.05 — **VSNL1, LZTFL1, GAS6, SREBF2** — are precisely the top of the machine-learning consensus (ranks 1, 2, 5, 4). So the study contains two internally coherent readings, and the manuscript must state both:

- **the mechanistic reading** (nominal P to FDR 0.10): suppressed AMPK output, PGC-1α down, cyclin D1 up, confirmed in defined cell types — the basis of the in vivo experiment;
- **the statistically strict reading** (FDR 0.05): VSNL1, LZTFL1, GAS6, SREBF2 — robust and predictive, but with no AMPK connection, no cell-type confirmation and no mechanism.

With 16 T2D and 11 lean livers the data cannot decide between them. The protein-level experiment (p-AMPKα Thr172/AMPKα, PGC-1α, cyclin D1, p21 across Control / Ob-T2D / GV1 / metformin) is what discriminates, and that is the argument for running it.

---

## 3. The selection method

Five filters, applied in a fixed order, each answering a different question. A gene is a Tier 1 target only if it passes every filter that its data type allows.

| # | Filter | Question | Rule |
|---|---|---|---|
| 1 | **Disease association** | Does the gene differ in the disease, in the right direction? | limma on 27 human liver samples (16 obese T2D vs 11 lean), P < 0.05 and \|log2FC\| > 0.5, two-dataset batch model |
| 2 | **Network + phenotype context** | Is it part of a coordinated module and of the phenotype being studied? | WGCNA module with T2D p < 0.1, intersected with a curated hyperglycemia / insulin-resistance gene set (GO + KEGG + HPO + CTD) → **19 candidates** |
| 3 | **Mechanism** | Does it sit on a pathway the data select, and on an edge that is actually drawn? | Six-criterion scoring of KEGG maps on human bulk data only → AMPK signalling; axis = KGML edges with **≥ 1 node inside the 19-gene set** (decision R23) |
| 4 | **Cell-type confirmation** | Does the same change happen in a hyperglycemic animal, and in which cells? | Mouse STZ vs control atlas, per-mouse pseudobulk (3 vs 4 mice), same direction as human, P < 0.05, plus the reference framework's co-localization and correlation tests with a detection-matched background null |
| 5 | **Independent human evidence** | Does it hold outside the discovery cohort? | Replication in GSE23343; multivariate contribution by LASSO/RF/SVM consensus; causal test by LD-clumped cis-MR on a 660,000-sample T2D GWAS |

Tiering is then mechanical (`24_target_scoring.R`): Tier 1 = filters 1–4 plus at least one arm of filter 5; Tier 2 = confirmed in filter 4 but not a Fig 1 candidate; Tier 2b = Fig 1 candidate on the map that filter 4 cannot test; Tier 3 = required mechanistic node that transcripts cannot measure.

---

## 4. Why this is scientifically correct

**1. The funnel is closed.** Every stage reports on genes that entered from the stage before. When the audit (script 22) showed the axis shared zero genes with the candidate set, the axis rule was changed rather than the story — the reported axis must contain a Fig 1 candidate. Chains that quietly change gene set between figures can always produce a clean-looking result.

**2. Selection and confirmation use different data.** The pathway and the axis were chosen from human bulk data alone; the mouse single-cell data were only ever used to test that choice. When the first attempt reversed this (choosing the pathway that happened to pass in mouse), it was recorded as an error and undone (decision R17). This is what keeps the single-cell step a test rather than a second selection.

**3. Every "significant" claim is referred to the right null.**
- Single-cell correlations are compared with 2,000 detection-matched random gene pairs, because with 2,875 cells a Pearson R of 0.07 is "significant" while being meaningless (empirical P 0.009 is the statement that survives).
- The machine-learning AUC is reported against a 200-shuffle permutation null (0.489) and against a nested CV in which feature selection is re-run inside every training fold.
- MR uses LD clumping against 1000 Genomes EUR, Steiger directionality, Cochran's Q and MR-Egger intercepts; a first pass using distance pruning produced four "causal" genes that dissolved once LD was handled properly, and that correction is documented (R21).

**4. Multiple testing is controlled at the level decisions are made.** BH within tissue for the single-cell edge screen (81 pairs), BH across the 24 MR genes, BH across cell types for the pathway scores.

**5. The animal evidence is matched to the cell type, not to the organ.** Pseudobulk per mouse (not per cell) is the primary test, because cells within a mouse are not independent replicates — the common inflation in single-cell papers. Cell-level Wilcoxon statistics are shown alongside for comparability with the reference framework.

**6. Where the data cannot answer, the tier says so.** AMPK, SIRT1 and p53 are regulated by phosphorylation and protein stability. The human map shows this directly: the drawn activation edge PRKAA2 → PPARGC1A is *discordant* (subunit up, target down). Rather than report a transcript-level AMPK result that cannot exist, these genes are placed in Tier 3 with protein assays specified. Likewise IRS2, IGF1 and FBP1 are Tier 2b because the atlas has 69 hepatocytes and 324 proximal-tubule cells — absence of evidence, stated as such.

**7. Negative results are kept.** MR nominates no causal target; the 19-gene panel does not transport to GSE23343 (external AUC 0.41–0.67 against 0.87–0.89 internal); the anchored pair is not co-localized or correlated. These are reported in the same tables as the positive findings, and they bound the claim: **this study proposes a mechanism to test in vivo, not a validated drug target or a diagnostic signature.**

**8. Species and model mismatch is handled explicitly, not ignored.** Human data are obese, insulin-resistant T2D; the mouse atlas is STZ, insulin-deficient. The lipogenic arm therefore inverts (Srebf1, Fasn, Fbp1 fall in STZ) and this is reported rather than dropped. Only the part of the phenotype the two share — hyperglycemic stress and the arrest program — is carried forward as cross-species evidence.

---

## 5. What would move a gene up a tier

| Gene | What is missing | Experiment that would supply it |
|---|---|---|
| PPARGC1A | independent human replication; causal genetics | PGC-1α protein in the GV1 cohort; liver eQTL MR when a better-powered liver eQTL dataset exists (GTEx liver n = 208 is not enough) |
| CCND1 | mechanism (currently a readout) | cyclin D1 in endothelium by IF with CD31, with and without AMPK activation (metformin arm) |
| CDKN1A | Fig 1 membership | would require changing the phenotype gene set *a priori*, not after the fact |
| PRKAA1/2, SIRT1 | any transcript-level evidence — by design | p-AMPKα Thr172 / AMPKα ratio and SIRT1 protein across the four groups; this is the decisive experiment for the whole model |

## 6. Files

| File | Contents |
|---|---|
| `results/targets/target_tiers.csv` | every gene with all evidence columns and its tier |
| `results/targets/target_evidence_long.csv` | the same in long form, for plotting |
| `results/consistency/gene_stage_matrix.csv` | gene × stage matrix behind §4.1 |
| `scripts/24_target_scoring.R` | the tier rule, applied mechanically |
