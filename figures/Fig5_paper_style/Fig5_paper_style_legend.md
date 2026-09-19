# Figure 5, reference-paper style - legend

Drawn to match Xu et al., *Phytomedicine* 154 (2026) 158050, Fig. 4, with this project's data.
The Nature/Cell-contract version of the same figure is in `figures/Fig5_MR/`.
**Ship one or the other, not both.**

## Read this before the caption

**The reference's Fig. 4 nominates a causal druggable target. This figure does not.** PAK1 reached IVW P = 0.039 there and was carried forward. Here neither Tier 1 target is causal for T2D: PPARGC1A OR 0.983, P = 0.16; CCND1 OR 1.033, P = 0.55. The 2 gene(s) that do reach FDR < 0.05 (SREBF1, SIRT1) are listed in `docs/target_selection.md` as *not* targets, because each is supported on this axis alone. Panel A therefore carries three exposures rather than one, so the null is on the figure rather than behind it.

**Fig. 5. Results of MR analysis.** (A) Forest plot illustrating the causal inference between the two Tier 1 targets (PPARGC1A, CCND1), the strongest MR hit (SREBF1) and type 2 diabetes. (B) SNP-level effect on the outcome against effect on the exposure for PPARGC1A, with one fitted line per MR method. (C) Forest plot for each SNP of the PPARGC1A analysis. (D) Results of Steiger filtering for directionality testing.

## Where this data forces a deviation from the reference

- **(A) has three exposures, not one,** for the reason given above. Rows are bold where P < 0.05.
- **(B) and (C) feature PPARGC1A,** the better-instrumented of the two Tier 1 targets (12 SNPs against 5 for CCND1). Both panels are diagnostics of pleiotropy and heterogeneity, and they are worth showing on a null estimate: the reference's versions carry the same information about a positive one.
- **CCND1 is the one exposure whose methods disagree.** IVW gives OR 1.033 (P = 0.55) while MR-Egger gives OR 0.788 (P = 0.011) with an intercept at P = 0.006, i.e. directional pleiotropy on 5 instruments. It is reported, not interpreted.
- **(D) labels the three exposures of panel A only;** all 25 instrumented genes are plotted. Exposure R2 exceeds outcome R2 for 25 of 25, so the instruments act on expression first.
- Exposures are whole-blood eQTLs (eQTLGen), the outcome is GCST006867. The primary analysis is the LD-clumped one; the earlier distance-pruned pass over-called SERPINE1, PPARGC1A and PRKAA1 and is withdrawn (decision R21).

## Statistics to quote

| Quantity | Value |
|---|---|
| Genes with usable instruments | 21 |
| Genes at FDR < 0.05 | 2 (SREBF1, SIRT1) |
| Strongest estimate | SREBF1, OR 0.903 (95% CI 0.872-0.935), P = 9.01E-09 |
| PPARGC1A (Tier 1) | 12 SNPs, IVW OR 0.983 (95% CI 0.959-1.007), P = 0.16 |
| CCND1 (Tier 1) | 5 SNPs, IVW OR 1.033 (95% CI 0.928-1.149), P = 0.55 |
| PPARGC1A heterogeneity | Cochran Q = 8.87, P = 0.63 |
| PPARGC1A pleiotropy | MR-Egger intercept P = 0.80 |
| Steiger, correct direction | 25 of 25 genes |
