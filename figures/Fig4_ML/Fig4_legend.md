# Figure 4 legend - classification of T2D from the candidate panel

**Fig. 4 | A 19-gene candidate panel does not generalise to an independent cohort.**
**a**, Cross-validated AUC for LASSO, random forest and SVM under two schemes: panel-fixed CV, in which the gene panel selected on the full data is held constant across folds, and nested CV, in which selection is repeated inside each fold. Grey, AUC from 200 label permutations. Boxes, median and interquartile range across CV repeats.
**b**, ROC curves in the independent cohort GSE23343 (17 samples); the legend gives AUC. Score is a parameter-free signature score.
**c**, LASSO selection frequency across folds and random-forest permutation importance for each candidate gene.
**d**, Predicted probability (or signature score) in GSE23343 by true group; points are samples.
**e**, Candidate-gene expression, z-scored within cohort, in the discovery data and in GSE23343.

## Caveats that belong in the text, not the figure

- The gap between panel-fixed and nested CV in **a** is the selection bias: panel-fixed median AUC 0.93 vs nested 0.88 (LASSO). Only the nested estimate is honest, and it is close to the permutation null.
- No model separates groups in the independent cohort at conventional confidence: best AUC 0.67 (RF), 95% CI 0.37-0.97, and every CI includes 0.5.
- With 16 vs 11 samples in discovery and 10 vs 7 in validation, all of these intervals are wide; the figure should be read as a negative result, not a ranking of models.
- **e** is z-scored within cohort, so it shows relative pattern, not comparable absolute expression.

## Statistics to quote

| Quantity | Value |
|---|---|
| Permutation null AUC (median) | 0.49 |
| Panel-fixed CV AUC, median | LASSO 0.93; RF 0.91; SVM 0.91 |
| Nested CV AUC, median | LASSO 0.88; RF 0.88; SVM 0.86 |
| External AUC (95% CI) | LASSO 0.41 (0.06-0.77); RF 0.67 (0.37-0.97); SVM 0.63 (0.32-0.93); Score 0.61 (0.31-0.92) |
| Most stable LASSO feature | VSNL1 (selected in 96% of folds) |
| Top RF feature | VSNL1 |
