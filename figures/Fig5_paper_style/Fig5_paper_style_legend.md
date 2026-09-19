# Figure 5, reference-paper style - legend

Drawn to match Xu et al., *Phytomedicine* 154 (2026) 158050, Fig. 4, with this project's data.
The Nature/Cell-contract version of the same figure is in `figures/Fig5_MR/`.
**Ship one or the other, not both.**

## What this figure is restricted to

The reference runs MR on a gene taken from its own selected KEGG map: PAK1, a node of regulation of actin cytoskeleton. The map selected here is **AMPK signalling (hsa04152)**, so the exposures on this figure are the **9 genes of that map with usable instruments**, out of 122 genes on the map. The 15 instrumented genes that are not on the map (p53-arrest genes and Fig 1 candidates without KEGG annotation) are excluded from A, B and C, and appear in D only as grey points.

**Unlike the wider gene set, the selected map does contain causal genes.** 2 genes reach FDR < 0.05: **SREBF1** (OR 0.903, 95% CI 0.872-0.935, P = 9.01E-09); **SIRT1** (OR 1.027, 95% CI 1.011-1.044, P = 9.47E-04). Both are nodes of hsa04152, so the reference's structure - a causal gene drawn from the selected pathway - is reproduced rather than only imitated.

**Fig. 5. Results of MR analysis.** (A) Forest plot illustrating the causal inference between the 9 instrumented genes of the selected KEGG map (hsa04152) and type 2 diabetes, by every applicable MR estimator. (B) SNP-level effect on the outcome against effect on the exposure for SIRT1, with one fitted line per MR method. (C) Forest plot for each SNP of the SIRT1 analysis. (D) Results of Steiger filtering for directionality testing; genes of the selected map are named.

## Where this data forces a deviation from the reference

- **(A) carries 9 exposures, not one.** The reference had a single prior druggable candidate (PAK1, from Finan et al. 2017) and tested it. There is no equivalent prior nomination here, so every instrumented gene of the map is tested and the whole screen is shown. Rows are bold where P < 0.05.
- **(B) and (C) feature SIRT1, not the strongest hit.** SREBF1 has the smaller P (9.01E-09) but only 2 instruments, which is too few for MR-Egger, the weighted median or the weighted mode, so its pleiotropy and heterogeneity diagnostics would be empty. SIRT1 carries 4 instruments and 4 estimators, and is significant in 3 of them.
- **Neither Tier 1 target is causal.** PPARGC1A OR 0.983, P = 0.16; CCND1 OR 1.033, P = 0.55. `docs/target_selection.md` nominates these two on expression and cell-type evidence; MR does not support either, and that disagreement is reported rather than smoothed over.
- **SIRT1 is causal in the direction opposite to the usual reading.** Higher genetically predicted expression associates with *higher* T2D risk (OR 1.027 per SD). It is reported as the estimate came out; this figure does not interpret the direction.
- **(D) plots the 25 genes with a Steiger test** and names the 9 on the map. That is one more than the 24 with a causal estimate: GAS6 passes Steiger but has no usable instrument after clumping. Exposure R2 exceeds outcome R2 for 25 of 25, so the instruments act on expression first.
- **FDR is the one computed over all 24 genes with a causal estimate, not recomputed over the 9 shown.** Restricting the figure to the selected map does not undo the multiple testing that was actually done, so the q values here are the conservative ones.
- **Methods disagree for PRKAA1 and CAMKK2.** PRKAA1 is null by IVW (OR 0.865, P = 0.23) but significant by the weighted median and weighted mode (OR 0.765, P = 2.19E-05), and CAMKK2 only by the weighted mode. With 5 instruments each, that pattern is what a single outlying instrument produces; it is reported, not interpreted.
- Exposures are whole-blood eQTLs (eQTLGen), the outcome is GCST006867. The primary analysis is the LD-clumped one; the earlier distance-pruned pass over-called SERPINE1, PPARGC1A and PRKAA1 and is withdrawn (decision R21).

## Statistics to quote

| Quantity | Value |
|---|---|
| Selected KEGG map | hsa04152, 122 genes, 104 measured in the bulk data |
| Map genes with usable instruments | 9 (CAMKK2, CCND1, FASN, FBP1, IRS2, PPARGC1A, PRKAA1, SIRT1, SREBF1) |
| Map genes at FDR < 0.05 | 2 (SREBF1, SIRT1) |
| Strongest on-map estimate | SREBF1, OR 0.903 (95% CI 0.872-0.935), P = 9.01E-09, 2 SNPs |
| SIRT1 (featured in B, C) | 4 SNPs, IVW OR 1.027 (95% CI 1.011-1.044), P = 9.47E-04 |
| SIRT1 heterogeneity | Cochran Q = 0.68, P = 0.88 |
| SIRT1 pleiotropy | MR-Egger intercept P = 0.62 |
| Tier 1 targets | PPARGC1A OR 0.983, P = 0.16; CCND1 OR 1.033, P = 0.55 |
| Steiger, correct direction | 25 of 25 genes |
