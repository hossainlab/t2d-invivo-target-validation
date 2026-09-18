# Figure 3 (axis panels) legend - kidney

**Fig. 3 (continued) | Pathway-axis evidence in STZ-diabetic mouse kidney.**
**a**, AMPK-PGC-1a axis genes in endothelial, the cell type with the most human-concordant axis genes. Points are pseudobulk log2 fold changes (STZ vs Control, edgeR quasi-likelihood); colour marks agreement with the human T2D direction and shape marks significance. An asterisk marks a Fig. 1g candidate gene.
**b**, UCell score of the candidate genes lying on the AMPK map, per cell type, as the STZ minus Control difference of per-mouse means.
**c**, Cell-level Pearson correlation of Ppargc1a–Ccnd1 in endothelial against 2,000 detection-matched random gene pairs from the same cell type. Bar, the 95th percentile of the null.
**d**,**e**,**f**, As **a**,**b**,**c** for the p53 arrest axis and the Cdkn1a–Ccnd1 pair in endothelial.

## Statistics to quote

| Quantity | Value |
|---|---|
| Key cell type (anchored axis) | Endothelial |
| Ppargc1a-Ccnd1, cells | 2,875 |
| Ppargc1a-Ccnd1, co-expressing cells | 0.2% |
| Ppargc1a-Ccnd1, cell-level R | −0.013 (P = 0.50) |
| Ppargc1a-Ccnd1, null mean / 95th pct / empirical P | 0.005 / 0.038 / 0.810 |
| Cdkn1a-Ccnd1, co-expressing cells | 11.1% |
| Cdkn1a-Ccnd1, cell-level R | 0.072 (P = 1.2 × 10⁻⁴) |
| Cdkn1a-Ccnd1, null mean / 95th pct / empirical P | 0.014 / 0.050 / 0.009 |
| Axis genes concordant with human (anchored) | 1 of 7 |
| Axis genes concordant with human (p53) | 5 of 8 |

## Caveats that belong in the text, not the figure

- Neither pair is co-expressed above chance: Ppargc1a-Ccnd1 gives empirical P = 0.810 and Cdkn1a-Ccnd1 gives P = 0.009 against detection-matched null pairs. The axis claim rests on the gene-level concordance in **a** and **d**, not on co-expression.
- Pseudobulk tests use 3 Control and 4 STZ mice; cell-type-level effects with few eligible mice are unstable, and the per-mouse cell counts are uneven (STZ_4 dominates the fibroblast compartment).
- The cross-species direction check compares mouse pseudobulk against the human bulk contrast, which is an association, not a validation of the human effect.
