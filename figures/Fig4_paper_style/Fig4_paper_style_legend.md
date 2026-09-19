# Figure 4, reference-paper style - legend

Drawn in the visual language of Xu et al., *Phytomedicine* 154 (2026) 158050.
The Nature/Cell-contract version of the same figure is in `figures/Fig4_ML/`.
**Ship one or the other, not both.**

## This figure has no panel-for-panel source in the reference

The reference paper runs no classifier, and its Figures 5 to 11 are chemistry and wet-lab validation that this study does not have. Figures 1, 2, 3 and 5 of this project each reproduce a reference figure panel for panel; this one reproduces only the reference's **visual language**, so that the shipped set reads as one figure set.

The panel is the 8-gene Fig 1 candidate set. Nothing about the analysis changes: every number is what script 21 wrote.

**Fig. 4. Classification of T2D from the candidate panel.** (A) Cross-validated AUC for LASSO, random forest and SVM under two schemes: panel-fixed CV, in which the gene panel selected on the full data is held constant across folds, and nested CV, in which selection is repeated inside each fold. Grey, AUC from 200 label permutations; dashed line, chance. (B) ROC curves in the independent cohort GSE23343 (17 samples); the key gives AUC. Score is a parameter-free signature score. (C) LASSO selection frequency across folds and random-forest permutation importance for each candidate gene. (D) Predicted probability, or signature score, in GSE23343 by true group; points are samples. (E) Candidate-gene expression, z-scored within cohort and clipped at +/-2, in the discovery data and in GSE23343.

## Read this with the figure

- **The two CV estimates in (A) disagree, and only the lower one is honest.** Panel-fixed CV gives LASSO a mean AUC of 0.92; nested CV, which repeats gene selection inside every fold, gives 0.87. The gap is selection bias, not model quality.
- **No model separates the groups in the independent cohort at conventional confidence.** The best external AUC is 0.64 (Score), 95% CI 0.33-0.96, and every confidence interval includes 0.5. LASSO is below chance at 0.43. This is a negative result and should be reported as one.
- **(C) ranks the genes the models lean on.** VSNL1 and PPARGC1A are selected in 100% and 97% of LASSO folds. With 8 features and 27 discovery samples, that stability says the panel is small and internally consistent, not that these are validated markers.
- With 16 vs 11 samples in discovery and 10 vs 7 in validation, every interval on this figure is wide. Read it as a negative result, not as a ranking of models.
- **(E) is z-scored within cohort**, so it shows relative pattern; absolute expression is not comparable across the two cohorts.

## Statistics to quote

| Quantity | Value |
|---|---|
| Panel-fixed CV AUC (optimistic) | LASSO 0.92, RF 0.94, SVM 0.90 |
| Nested CV AUC (honest) | LASSO 0.87, RF 0.89, SVM 0.88 |
| Permutation null AUC | 0.48 (sd 0.09), 200 permutations |
| External AUC, GSE23343 | LASSO 0.43, RF 0.61, SVM 0.61, Score 0.64 |
| Best external model | Score, AUC 0.64 (95% CI 0.33-0.96) |
| Most stable features | VSNL1 (LASSO 100%), PPARGC1A (LASSO 97%), SREBF2 (LASSO 76%) |
