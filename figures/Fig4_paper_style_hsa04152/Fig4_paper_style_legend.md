# Figure 4, reference-paper style - legend

Drawn in the visual language of Xu et al., *Phytomedicine* 154 (2026) 158050.
The Nature/Cell-contract version of the same figure is in `figures/Fig4_ML/`.
**Ship one or the other, not both.**

## This figure has no panel-for-panel source in the reference

The reference paper runs no classifier, and its Figures 5 to 11 are chemistry and wet-lab validation that this study does not have. Figures 1, 2, 3 and 5 of this project each reproduce a reference figure panel for panel; this one reproduces only the reference's **visual language**, so that the shipped set reads as one figure set.

**The panel is restricted to the selected KEGG map (hsa04152).** Of the 8 Fig 1 candidate genes, 4 are nodes of that map: CCND1, IGF1, IRS2, PPARGC1A. The other four (VSNL1, SREBF2, IGFBP1, IGFBP2) appear in no KEGG pathway at all and are excluded, so the classifier is tested on the same genes the pathway, single-cell and MR figures are about. The restriction is applied inside the nested-CV selection rule as well, not only to the fixed panel.

**Fig. 4. Classification of T2D from the 4-gene hsa04152 panel.** (A) Cross-validated AUC for LASSO, random forest and SVM under two schemes: panel-fixed CV, in which the gene panel selected on the full data is held constant across folds, and nested CV, in which selection is repeated inside each fold. Grey, AUC from 200 label permutations; dashed line, chance. (B) ROC curves in the independent cohort GSE23343 (17 samples); the key gives AUC. Score is a parameter-free signature score. (C) LASSO selection frequency across folds and random-forest permutation importance for each candidate gene. (D) Predicted probability, or signature score, in GSE23343 by true group; points are samples. (E) Candidate-gene expression, z-scored within cohort and clipped at +/-2, in the discovery data and in GSE23343.

## Read this with the figure

- **The two CV estimates in (A) disagree, and only the lower one is honest.** Panel-fixed CV gives LASSO a mean AUC of 0.83; nested CV, which repeats gene selection inside every fold, gives 0.78. The gap is selection bias, not model quality.
- **No model separates the groups in the independent cohort at conventional confidence.** The best external AUC is 0.71 (LASSO), 95% CI 0.44-0.98, and every confidence interval includes 0.5. This is a negative result and should be reported as one.
- **(C) ranks the genes the models lean on.** CCND1 and PPARGC1A are selected in 97% and 100% of LASSO folds. With 4 features and 27 discovery samples, that stability says the panel is small and internally consistent, not that these are validated markers.
- With 16 vs 11 samples in discovery and 10 vs 7 in validation, every interval on this figure is wide. Read it as a negative result, not as a ranking of models.
- **(E) is z-scored within cohort**, so it shows relative pattern; absolute expression is not comparable across the two cohorts.

## Statistics to quote

| Quantity | Value |
|---|---|
| Panel-fixed CV AUC (optimistic) | LASSO 0.83, RF 0.82, SVM 0.76 |
| Nested CV AUC (honest) | LASSO 0.78, RF 0.78, SVM 0.77 |
| Permutation null AUC | 0.48 (sd 0.08), 200 permutations |
| External AUC, GSE23343 | LASSO 0.71, RF 0.59, SVM 0.60, Score 0.67 |
| Best external model | LASSO, AUC 0.71 (95% CI 0.44-0.98) |
| Most stable features | CCND1 (LASSO 97%), PPARGC1A (LASSO 100%), IGF1 (LASSO 39%) |
| Same models on the unrestricted 8-gene panel | LASSO 0.43, RF 0.61, SVM 0.61, Score 0.64 |
