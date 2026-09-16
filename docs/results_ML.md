# Machine learning: how well does the candidate panel classify type 2 diabetes?

**Status (2026-09-16):** complete (decision R22). **Script:** `21_ML_classifier.R`. **Figure:** `figures/Fig4_ML/Fig4_ML.pdf`.
**Not a reference-paper step** — Xu et al. ran no machine learning. This is an addition, and it is built so that the usual failure mode of small-n signature papers is visible rather than hidden.

> **Headline:** the 19-gene panel separates T2D from lean liver **within** the discovery cohort (nested-CV AUC 0.87–0.89, permutation P = 0.005) but **does not transport**: in the independent cohort GSE23343 the AUC is 0.41 (LASSO), 0.67 (random forest), 0.63 (SVM) and 0.61 (parameter-free signature score). The panel is a description of these 27 samples, not a diagnostic classifier.

## 1. Design

| Element | Detail |
|---|---|
| Discovery | 27 liver samples (16 obese T2D, 11 lean), ComBat-merged GSE15653 + GSE64998 |
| Features | the 19-gene candidate panel (M13); all 19 are measured in both cohorts |
| External | GSE23343, human liver, **10 T2D / 7 NGT**, hgu133plus2, RMA, MaxMean collapse — never used for training or selection |
| Scaling | gene-wise z-scoring **within each dataset**, since discovery is ComBat-corrected two-platform data and the external set is a third platform |
| Models | LASSO logistic regression (glmnet), random forest (ranger, 1000 trees), SVM with RBF kernel (e1071); plus a parameter-free signature score = mean z(up genes) − mean z(down genes), directions from the discovery contrast |
| Resampling | stratified 5-fold CV, 20 repeats (panel-fixed) and 10 repeats (nested) |

Four evaluations, reported side by side (Fig 4A):

1. **Panel-fixed CV** — the usual published design. Optimistic, because the 19 genes were chosen using all 27 samples.
2. **Nested CV** — the selection rule (limma DEG P < 0.05 and |logFC| > 0.5 ∩ hyperglycemia gene set) is re-run inside every training fold, so no test sample influences feature choice.
3. **Permutation null** — 200 label shuffles through the same CV, to show what this design yields by chance at n = 27.
4. **External cohort** — trained on all 27 discovery samples, tested once on GSE23343.

## 2. Results (`results/ml/ml_summary.csv`)

| Evaluation | LASSO | Random forest | SVM | Signature score |
|---|---|---|---|---|
| Panel-fixed CV (optimistic) | 0.918 ± 0.028 | 0.907 ± 0.020 | 0.906 ± 0.020 | — |
| **Nested CV (honest)** | **0.872 ± 0.052** | **0.888 ± 0.026** | **0.878 ± 0.032** | — |
| Permutation null | 0.489 ± 0.087 (95th pct 0.659) | — | — | — |
| **External GSE23343** | **0.414** | **0.671** | **0.629** | **0.614** |

- Permutation P for the observed LASSO CV AUC: **0.005**. The within-cohort separation is real, not a small-n artefact.
- The gap between panel-fixed (0.918) and nested (0.872) CV is modest, so **selection bias is not the main problem** — nested folds re-selected 15.5 genes on average, largely the same panel.
- The external drop is the finding: **0.87–0.89 internal → 0.41–0.67 external**, with confidence intervals spanning 0.5 for every model (e.g. RF 0.67, 95 % CI 0.37–0.97). LASSO is below chance.
- The parameter-free score does no better (0.61), so the failure is not overfitting of a complex model — the **panel's direction of effect does not reproduce** in the second cohort.

**Why this is consistent with the rest of the study.** GSE23343 compares T2D with normal glucose tolerance without matching obesity, whereas the discovery contrast is obese T2D vs lean. The MR analysis (`docs/results_MR.md`) independently found no causal support for these genes. Both point the same way: the panel describes a liver state in this particular contrast rather than a portable T2D signature.

## 3. Which genes carry the within-cohort signal (`feature_ranking.csv`, Fig 4C–D)

Consensus of LASSO selection frequency over 100 folds, random-forest permutation importance and SVM-RFE rank:

| Rank | Gene | LASSO freq | RF importance | SVM-RFE |
|---|---|---|---|---|
| 1 | VSNL1 | 0.96 | 0.053 | 3 |
| 2 | LZTFL1 | 0.92 | 0.017 | 1 |
| 3 | PPARGC1A | 0.89 | 0.011 | 2 |
| 4 | SREBF2 | 0.68 | 0.012 | 7 |
| 5 | GAS6 | 0.45 | 0.029 | 13 |
| 6 | HMOX1 | 0.77 | 0.009 | 8 |
| 7 | SERPINE1 | 0.66 | 0.007 | 9 |

**CCND1 and the other p53-axis genes rank near the bottom** (CCND1 LASSO frequency 0.02, FBP1 0.03, KHK 0.02, PNPLA3 0.02). The genes that drive classification are not the genes the pathway and single-cell work nominates as mechanism — worth stating plainly, because it is a common unstated tension in this kind of paper. PPARGC1A appearing third is the one point of contact with the AMPK arm.

## 4. Limitations

- n = 27 discovery / 17 external. Every estimate is wide; the external AUCs cannot exclude 0.5 in either direction.
- ComBat was fitted on all 27 discovery samples before CV, a mild optimism that the nested design does not remove.
- The hyperglycemia gene set used inside the nested folds is external prior knowledge, held fixed by design.
- The signature-score directions come from the discovery contrast, so only its external AUC is a clean test.
- No hyperparameter search beyond glmnet's internal lambda CV; with this n, tuning would add variance rather than performance.

## 5. Figure 4 panels (`figures/Fig4_ML/`)

| Panel | Contents |
|---|---|
| Fig4A | AUC across the four evaluations, three models |
| Fig4B | External ROC curves, three models plus the signature score |
| Fig4C | LASSO selection frequency (stability) |
| Fig4D | Random-forest permutation importance |
| Fig4E | External predicted probabilities by true group |
| Fig4F | Top consensus genes, z-scored, in both cohorts |
| Fig4_ML.pdf | assembled figure |
