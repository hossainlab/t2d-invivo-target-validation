# Figure 4 — interpretation and defence

Companion to `Fig4_legend.md`. Panel: the 4 candidates on hsa04152 — CCND1, IGF1, IRS2, PPARGC1A.

---

## 1. The claim this figure supports

> The 4-gene AMPK-map panel does **not** classify type 2 diabetes in an independent human cohort at
> conventional confidence. Restricting the panel to the selected pathway improves transportability
> substantially (external LASSO AUC 0.43 → 0.71), but every confidence interval still includes 0.5.

This is a **negative result reported as one**, plus a methodological finding about pathway
restriction. Presenting it as a positive would be indefensible; presenting it as a failure of the
project would be wrong. It is the figure that tells you how far the candidate set generalises, and
the answer is "not far, on this sample size".

---

## 2. Panel by panel

**A — the two cross-validation estimates disagree, and only the lower one is honest.**
Panel-fixed CV gives LASSO 0.83, RF 0.82, SVM 0.76. Nested CV — in which gene selection is repeated
*inside every training fold*, so no test sample influences feature choice — gives 0.78, 0.78, 0.77.
The gap is **selection bias**, not model quality. The permutation null (200 label shuffles) sits at
0.48 ± 0.08, confirming the design produces chance performance when the labels carry no information.

**B — no model separates the groups in GSE23343.**
LASSO 0.71, Score 0.67, SVM 0.60, RF 0.59. The best interval is 0.44–0.98. Every interval spans 0.5.

**C — the panel is internally consistent but that is not validation.**
PPARGC1A is selected in 100 % of LASSO folds and CCND1 in 97 %. With **4 features and 27 samples**,
near-total selection frequency is close to inevitable — it says the panel is small and the genes are
mutually non-redundant, not that they are validated markers. RF permutation importance ranks CCND1
and PPARGC1A top, IGF1 essentially zero (slightly negative).

**D — the group distributions overlap in the external cohort.** Medians differ in the expected
direction for three of four models, but the distributions are not separated.

**E — the expression pattern is visible but not transportable.** Within the discovery cohort the
four genes show the expected contrast; in GSE23343 the pattern is weaker and noisier. Since both are
z-scored *within* cohort, this panel shows relative pattern only.

---

## 3. The finding that is worth reporting

**Restricting the panel to the selected pathway improved external performance.**

| Panel | External AUC (GSE23343) |
|---|---|
| 8 Fig-1 candidates (unrestricted) | LASSO **0.43**, RF 0.61, SVM 0.61, Score 0.64 |
| 4 genes on hsa04152 | LASSO **0.71**, RF 0.59, SVM 0.60, Score 0.67 |

LASSO moves from *below chance* to 0.71. The four genes dropped (VSNL1, SREBF2, IGFBP1, IGFBP2)
appear in **no KEGG pathway at all**; they were contributing variance that did not transport across
platforms. This is a small-sample observation and should be stated as suggestive rather than
established — with n = 17 in the validation cohort, a 0.28 AUC difference is well within sampling
noise. But it is consistent with the general principle that pathway-constrained signatures generalise
better than unconstrained ones, and it is the kind of result that is only visible because both
analyses were run.

---

## 4. Why the analysis was done this way

**Why nested CV is on the figure at all.** Most small-n signature papers report the panel-fixed
number and stop. Showing both, side by side, with the gap labelled as selection bias, is the single
most defensible feature of this figure. The honest estimate is the lower one, and it is the one
quoted.

**Why a permutation null.** At n = 27 with a 4-gene panel, it is reasonable to ask what AUC the
*design* produces from noise. 200 label shuffles answer that: 0.48 ± 0.08. The observed nested
estimates (≈ 0.78) are clearly above it, so the discovery-side signal is real — it simply does not
transport.

**Why the restriction is applied inside the nested-CV rule, not only to the fixed panel.** If the
pathway filter were applied only to the fixed panel, the nested estimate would still be free to
select genes outside hsa04152, and the two numbers in panel A would not be measuring the same thing.
`scripts/21_ML_classifier.R` applies the map restriction inside `select_in_fold()` as well.

**Why a parameter-free signature score is included.** "Score" is the mean z of the human-up genes
minus the mean z of the human-down genes. It has no fitted parameters, so it transports across
platforms better than a trained classifier. Its external AUC (0.67) is comparable to LASSO's, which
tells you the trained models are not extracting much beyond direction.

---

## 5. Anticipated questions

**"Why include a figure that fails?"** Because the alternative is to report the panel-fixed CV
number (0.83) and omit the external validation, which is what the small-n signature literature
routinely does. The figure's purpose is to establish the *limit* of the candidate set's evidence, and
that limit is load-bearing for the thesis: it is why the project claims a mechanism supported by
expression and cell-type evidence rather than a validated biomarker. A negative result that
constrains the claim is more useful than a positive result that cannot be reproduced.

**"AUC 0.71 with a CI of 0.44–0.98 — isn't that just noise?"** It is consistent with noise, and the
legend says every interval includes 0.5. The correct statement is "no model separates the groups at
conventional confidence", not "the panel has AUC 0.71". Quote the interval, never the point estimate
alone.

**"n = 17 in validation is too small to conclude anything."** Correct, and the conclusion drawn is
correspondingly weak: *no evidence of generalisation*, not *evidence of no generalisation*. The
distinction matters and should be made explicitly. A larger independent liver cohort would be the
obvious next step; GSE23343 is what exists.

**"The reference paper has no ML figure. Why do you?"** It is an addition, declared as one in the
legend. The reference's chain runs straight from pathway to wet lab. Inserting a classification step
tests something the reference never asks: whether the candidate set carries enough information to be
useful outside the discovery data. The answer here is no, and knowing that prevents the candidate set
from being over-sold as a signature.

**"Doesn't a failed classifier undermine Figures 1 and 2?"** No, and this is worth separating
carefully. Classification asks whether 4 genes can *assign individual patients* to groups.
Figures 1 and 2 ask whether those genes are *differentially expressed and pathway-coherent at the
group level*. A gene can be robustly differentially expressed on average and still be useless for
individual prediction when within-group variance is large — which at 16 vs 11 samples it is. The two
results are compatible, and `docs/results_ML.md` states this.

---

## 6. What this figure does not show

- It does not show that the 4 genes are *not* differentially expressed — Figure 1 shows that they are.
- It does not show that pathway restriction reliably improves generalisation; n = 17 cannot support
  that claim, only suggest it.
- It does not rank the three model families meaningfully; their intervals overlap almost entirely.
- It does not test the AMPK mechanism. Classification and mechanism are different questions.
