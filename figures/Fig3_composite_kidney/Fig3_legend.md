# Figure 3 legend - single-cell atlas of kidney, STZ vs Control

**Fig. 3 | Single-cell transcriptomes of STZ-diabetic and control mouse kidney.**
**a**, UMAP of 6,224 cells from 7 mice (3 Control, 4 STZ) coloured by cell type; the legend gives the proportion of cells in each type.
**b**, Canonical markers used for annotation. Dot size, percentage of cells in the type with non-zero expression; colour, mean log-normalized expression.
**c**, Left, proportion of each cell type contributed by Control and STZ mice. Right, number of cells per type (log scale).
**d**,**e**, *Cdkn1a* (**d**) and *Ccnd1* (**e**) expression in endothelials (2,875 cells), Control vs STZ. Violins are width-scaled over all cells of that type; white points are per-mouse means; the percentage under each violin is the fraction of cells with non-zero expression. Brackets give the per-mouse t-test (n = 3 vs 4 mice).
**f**, Joint kernel density of *Cdkn1a* and *Ccnd1* co-expression (Nebulosa) on the UMAP of **a**; grey, cells below 2% of the maximum density.
**g**, *Cdkn1a* against *Ccnd1* in endothelials, the cell type with the strongest concordance with the bulk axis; line, linear fit with 95% confidence band.

## Statistics to quote

| Quantity | Value |
|---|---|
| Cells (total) | 6,224 |
| Cell types | 14 |
| Mice | 3 Control, 4 STZ |
| Cdkn1a, mean expression (Control / STZ) | 0.375 / 0.627 |
| Cdkn1a, cells detected (Control / STZ) | 21.6% / 29.8% |
| Panels d,e cell type | Endothelial (n = 2,875 cells) |
| Cdkn1a, cell-level Wilcoxon | P = 1.5 × 10⁻⁹ |
| Cdkn1a, per-mouse t-test (n = 3 vs 4) | P = 0.015 |
| Ccnd1, mean expression (Control / STZ) | 0.469 / 0.797 |
| Ccnd1, cells detected (Control / STZ) | 27.1% / 40.4% |
| Ccnd1, cell-level Wilcoxon | P < 2.2 × 10⁻¹⁶ |
| Ccnd1, per-mouse t-test (n = 3 vs 4) | P = 2.9 × 10⁻⁴ |
| Panel g cell type | Endothelial (n = 2,875 cells) |
| Panel g Pearson R | 0.07 (P = 1.2 × 10⁻⁴) |
| Panel g, % cells co-expressing both | 11.1% |
| Panel g background null (detection-matched random pairs) | mean R = 0.01, 95th percentile = 0.05, empirical P = 0.009 |

## Caveats that belong in the text, not the figure

- The cell-level Wilcoxon in **d**,**e** treats cells as independent replicates, which they are not; the brackets therefore report the per-mouse test. The two disagree in both directions here: Cdkn1a gives P = 1.5 × 10⁻⁹ across cells and P = 0.015 across mice, Ccnd1 gives P < 2.2 × 10⁻¹⁶ across cells and P = 2.9 × 10⁻⁴ across mice. With 3 vs 4 mice the per-mouse test is the valid one but has very little power, so neither should be read as strong evidence.
- The Cdkn1a-Ccnd1 correlation in **g** is negative (R = 0.07) and does not exceed a detection-matched background (empirical P = 0.009), i.e. the two genes are not co-expressed within endothelials above chance.
- Cell numbers are uneven between mice (STZ_4 contributes most fibroblasts); proportions in **c** should not be read as differential abundance without the propeller test in `results/sc/<tissue>_propeller_DA.csv`.

## Cells per type

| Cell type | Cells | % |
|---|---|---|
| Endothelial | 2,875 | 46.19 |
| Macrophage | 907 | 14.57 |
| Fibroblast/pericyte | 591 | 9.50 |
| T cell | 345 | 5.54 |
| PT | 324 | 5.21 |
| CD-IC | 296 | 4.76 |
| Distal tubule | 234 | 3.76 |
| B cell | 189 | 3.04 |
| LoH | 163 | 2.62 |
| cDC | 152 | 2.44 |
| Urothelium | 43 | 0.69 |
| Proliferating | 37 | 0.59 |
| Tubular epithelium (other) | 34 | 0.55 |
| Neutrophil | 34 | 0.55 |
