# Figure 1 legend - bulk liver transcriptome: differential expression and co-expression modules

**Fig. 1 | Differentially expressed genes and T2D-associated co-expression modules in human liver.**
**a**, Principal components of the merged GSE15653 + GSE64998 liver arrays before and after ComBat batch correction (27 samples: 16 T2D, 11 lean control).
**b**, Volcano plot of T2D vs control (limma, ~condition + dataset). Dashed lines, |log2 fold change| > 0.5 and nominal P < 0.05; 96 genes up and 82 down. The 10 most significant genes in each direction are labelled.
**c**, Scale-free topology fit and mean connectivity against soft-thresholding power; power 7 (red) was the lowest reaching a fit of 0.8 (dashed).
**d**, Gene dendrogram from the signed-hybrid network with the assigned module colours beneath.
**e**, Pearson correlation between each module eigengene and group (control and T2D, all 27 samples), T2D within each cohort separately, HbA1c (GSE15653 only, n = 14) and dataset of origin. Each cell gives r above P. Control is the complement of T2D, so its correlation is the exact negative and its P value identical. Key modules, in bold, are those with P < 0.1 for T2D pooled AND the same direction at P < 0.1 in both cohorts (decision M11b): red, yellow, tan, greenyellow.
**f**, Gene significance for T2D against module membership within each key module; line, linear fit.
**g**, Overlap of the 178 DEGs, the 418 genes in the WGCNA key modules and the 618 hyperglycaemia-associated genes present in the expression universe; the 8 genes shared by all three (bold) are the candidate set carried into Fig. 2.

## Caveats that belong in the text, not the figure

- † **e**, the Dataset column is computed on ComBat-corrected expression. ComBat removes the batch mean by construction, so this column is bounded near zero (max |r| = 0.08) and is *not* an independent test for residual batch effect. It is reported, not used as a filter; the per-cohort columns are what the key-module rule tests, and the honest batch check is **a**.
- Three modules with a strong pooled T2D correlation are excluded because they do not replicate: blue (r = −0.72 in GSE15653 vs +0.11 in GSE64998, i.e. the sign flips), green (0.69 vs 0.23) and turquoise (0.63 vs 0.16). blue and turquoise are also the two largest modules, so this is what cuts the key-module gene pool from 4,022 to 418.
- 84 genes reach FDR < 0.05 for T2D vs control and only 31 of those also clear |log2FC| > 0.5; carrying that set through the strict downstream rules leaves one intersecting gene (SREBF2, decision R17). **b** therefore uses nominal P < 0.05 and the DEG set is exploratory (decision M9).
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
| Key modules (T2D P < 0.1) | red, yellow, tan, greenyellow |
| Strongest module-trait correlation | red, r = 0.65, P = 2.2 × 10⁻⁴ |
| Three-way intersection | 8 genes |
