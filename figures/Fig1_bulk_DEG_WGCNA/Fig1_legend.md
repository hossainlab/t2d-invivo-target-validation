# Figure 1 legend - bulk liver transcriptome: differential expression and co-expression modules

**Fig. 1 | Differentially expressed genes and T2D-associated co-expression modules in human liver.**
**a**, Principal components of the merged GSE15653 + GSE64998 liver arrays before and after ComBat batch correction (27 samples: 16 T2D, 11 lean control).
**b**, Volcano plot of T2D vs control (limma, ~condition + dataset). Dashed lines, |log2 fold change| > 0.5 and nominal P < 0.05; 96 genes up and 82 down. The 10 most significant genes in each direction are labelled.
**c**, Scale-free topology fit and mean connectivity against soft-thresholding power; power 7 (red) was the lowest reaching a fit of 0.8 (dashed).
**d**, Gene dendrogram from the signed-hybrid network with the assigned module colours beneath.
**e**, Pearson correlation between each module eigengene and T2D status, HbA1c (GSE15653 only, n = 14) and dataset of origin. Each cell gives r above P. The 7 key modules taken forward are those with P < 0.1 for T2D: red, yellow, tan, greenyellow, blue, green, turquoise.
**f**, Gene significance for T2D against module membership within each key module; line, linear fit.
**g**, Overlap of the 178 DEGs, the 4,022 genes in the key modules and the 618 hyperglycaemia-associated genes present in the expression universe; the 19 genes shared by all three (bold) are the candidate set carried into Fig. 2.

## Caveats that belong in the text, not the figure

- † **e**, the Dataset column is computed on ComBat-corrected expression. ComBat removes the batch mean by construction, so this column is bounded near zero (max |r| = 0.08) and is *not* an independent test for residual batch effect. The honest batch check is **a**.
- No gene reaches FDR < 0.05 for T2D vs control; **b** therefore uses nominal P < 0.05 and the DEG set is exploratory (decision M9).
- **e**, module-trait P values are uncorrected across 14 modules x 3 traits.
- **f**, MM and GS are both computed from the same 27 samples, so the correlation is not independent evidence of module relevance.

## Statistics to quote

| Quantity | Value |
|---|---|
| Samples | 16 T2D, 11 lean control (GSE15653 + GSE64998) |
| Genes tested | 11,972 |
| DEGs (P < 0.05, |log2FC| > 0.5) | 96 up, 82 down |
| Genes at FDR < 0.05 | 84 |
| Soft-thresholding power | 7 |
| Modules (excluding grey) | 14 |
| Key modules (T2D P < 0.1) | red, yellow, tan, greenyellow, blue, green, turquoise |
| Strongest module-trait correlation | red, r = 0.65, P = 2.2 × 10⁻⁴ |
| Three-way intersection | 19 genes |
