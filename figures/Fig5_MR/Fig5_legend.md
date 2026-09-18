# Figure 5 legend - Mendelian randomisation of candidate gene expression on T2D

**Fig. 5 | Cis-eQTL Mendelian randomisation does not support a causal role for most candidates.**
**a**, Inverse-variance-weighted causal estimates for the 21 instrumentable candidate genes, as odds ratio for T2D per standard deviation of genetically predicted whole-blood expression (eQTLGen exposures, GCST006867 outcome). Bars, 95% confidence interval; the number at the right of each row is the instrument count. Coloured, FDR < 0.05.
**b**, SNP-level effect on expression against effect on T2D for the four best-instrumented genes; line, the IVW slope. Error bars, standard errors.
**c**, Leave-one-out IVW estimates; the vertical line is the estimate using all SNPs.
**d**, Estimates from LD-clumped instruments against the superseded distance-pruned pass, per gene; point size, instrument count after clumping.
**c**, The 3 genes with usable GTEx liver instruments (Wald ratio, single SNP each); kidney cortex yielded none.
**d**, Steiger directionality test: variance in expression explained by the instruments (exposure) against variance explained in T2D (outcome), per gene. Exposure exceeds outcome for 25 of 25 genes, i.e. the instruments act on expression first.

## Caveats that belong in the text, not the figure

- Panel a is the LD-clumped analysis (decision R21). The first pass used distance pruning and over-called: SERPINE1, PPARGC1A, PRKAA1 were significant there and not after clumping, which is what panel d shows.
- 2 of 21 genes reach FDR < 0.05, and the strongest is SREBF1 (OR 0.90, 95% CI 0.87-0.94, P = 9.0 × 10⁻⁹).
- Exposures are whole-blood eQTLs, not liver; **e** is the liver-specific check and is limited to single-SNP Wald ratios, which cannot be tested for pleiotropy.
- Instrument counts are small for most genes, so horizontal pleiotropy tests (Q, Egger intercept) have little power; they are tabulated in `results/mr/mr_results_blood_ldclumped.csv` rather than plotted.

## Statistics to quote

| Quantity | Value |
|---|---|
| Genes with instruments | 21 |
| Genes at FDR < 0.05 (LD-clumped) | 2 |
| Strongest estimate | SREBF1, OR 0.90 (95% CI 0.87-0.94), P = 9.0 × 10⁻⁹, FDR 0.000 |
| Lost after LD clumping | SERPINE1, PPARGC1A, PRKAA1 |
| Median instrument count | 4 |
| GTEx liver genes tested | CTSD, PLIN1, VSNL1 |
