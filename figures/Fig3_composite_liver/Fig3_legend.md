# Figure 3 legend - single-cell atlas of liver, STZ vs Control

**Fig. 3 | Single-cell transcriptomes of STZ-diabetic and control mouse liver.**
**a**, UMAP of 16,861 cells from 7 mice (3 Control, 4 STZ) coloured by cell type; the legend gives the proportion of cells in each type.
**b**, Canonical markers used for annotation. Dot size, percentage of cells in the type with non-zero expression; colour, mean log-normalized expression.
**c**, Left, proportion of each cell type contributed by Control and STZ mice. Right, number of cells per type (log scale).
**d**,**e**, *Cdkn1a* (**d**) and *Ccnd1* (**e**) expression in cholangiocytes (636 cells), Control vs STZ. Violins are width-scaled over all cells of that type; white points are per-mouse means; the percentage under each violin is the fraction of cells with non-zero expression. Brackets give the per-mouse t-test (n = 3 vs 4 mice).
**f**, Joint kernel density of *Cdkn1a* and *Ccnd1* co-expression (Nebulosa) on the UMAP of **a**; grey, cells below 2% of the maximum density.
**g**, *Cdkn1a* against *Ccnd1* in cholangiocytes, the cell type with the strongest concordance with the bulk axis; line, linear fit with 95% confidence band.

## Statistics to quote

| Quantity | Value |
|---|---|
| Cells (total) | 16,861 |
| Cell types | 13 |
| Mice | 3 Control, 4 STZ |
| Cdkn1a, mean expression (Control / STZ) | 0.252 / 0.533 |
| Cdkn1a, cells detected (Control / STZ) | 25.6% / 28.4% |
| Panels d,e cell type | Cholangiocyte (n = 636 cells) |
| Cdkn1a, cell-level Wilcoxon | P = 0.073 |
| Cdkn1a, per-mouse t-test (n = 3 vs 4) | P = 0.002 |
| Ccnd1, mean expression (Control / STZ) | 1.224 / 1.522 |
| Ccnd1, cells detected (Control / STZ) | 72.6% / 66.7% |
| Ccnd1, cell-level Wilcoxon | P = 0.003 |
| Ccnd1, per-mouse t-test (n = 3 vs 4) | P = 0.18 |
| Panel g cell type | Cholangiocyte (n = 636 cells) |
| Panel g Pearson R | −0.10 (P = 0.010) |
| Panel g, % cells co-expressing both | 18.7% |
| Panel g background null (detection-matched random pairs) | mean R = 0.06, 95th percentile = 0.15, empirical P = 0.992 |

## Caveats that belong in the text, not the figure

- The cell-level Wilcoxon in **d**,**e** treats cells as independent replicates, which they are not; the brackets therefore report the per-mouse test. The two disagree in both directions here: Cdkn1a gives P = 0.073 across cells and P = 0.002 across mice, Ccnd1 gives P = 0.003 across cells and P = 0.18 across mice. With 3 vs 4 mice the per-mouse test is the valid one but has very little power, so neither should be read as strong evidence.
- The Cdkn1a-Ccnd1 correlation in **g** is negative (R = −0.10) and does not exceed a detection-matched background (empirical P = 0.992), i.e. the two genes are not co-expressed within cholangiocytes above chance.
- Cell numbers are uneven between mice (STZ_4 contributes most fibroblasts); proportions in **c** should not be read as differential abundance without the propeller test in `results/sc/<tissue>_propeller_DA.csv`.

## Cells per type

| Cell type | Cells | % |
|---|---|---|
| B cell | 7,386 | 43.81 |
| HSC/fibroblast | 2,024 | 12.00 |
| T cell | 2,018 | 11.97 |
| Monocyte | 920 | 5.46 |
| pDC | 748 | 4.44 |
| LSEC | 745 | 4.42 |
| cDC | 731 | 4.34 |
| Cholangiocyte | 636 | 3.77 |
| NK cell | 596 | 3.53 |
| Neutrophil | 572 | 3.39 |
| Macrophage (other) | 214 | 1.27 |
| Kupffer cell | 202 | 1.20 |
| Hepatocyte | 69 | 0.41 |
