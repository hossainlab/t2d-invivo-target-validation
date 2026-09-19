# Figure 3, reference-paper style (liver) - legend

Drawn to match Xu et al., *Phytomedicine* 154 (2026) 158050, Fig. 3, with this project's data.
The Nature/Cell-contract version of the same figure is in `figures/Fig3_composite_liver_ppargc1a_ccnd1/`.
**Ship one or the other, not both.**

**Fig. 3. Results of scRNA-seq analysis.** (A) UMAP projection illustrating the single-cell atlas of mouse liver tissue. (B) Bubble plot showing marker genes for each cell type. (C) Bar plot displaying the proportion and number of cells from each sample source. (D-E) Analysis of Ppargc1a and Ccnd1 expression levels in cholangiocyte across experimental groups. Per-mouse t-test (n = 3 vs 4 mice); ns: P > 0.05, *: P < 0.05, **: P < 0.01, ***: P < 0.001, ****: P < 0.0001. (F) UMAP-based heatmap showing co-localization of Ppargc1a and Ccnd1 across the atlas. (G) Correlation between Ppargc1a and Ccnd1 expression levels in cholangiocyte.

## Where this data forces a deviation from the reference

- **Two groups, not three.** The reference contrasts Ctrl, CIA and Treat, so its D and E carry three brackets. This atlas is Control vs STZ, so each panel carries one.
- **The bracket is the per-mouse t-test (n = 3 vs 4 mice), not the cell-level Wilcoxon.** Cells within a mouse are not independent replicates, so a cell-level P is anticonservative. Both are tabulated below; the reference reported the cell-level test.
- **(F) and (G) do not agree, and that is the result.** The pair is co-localized in cholangiocyte (13.8% of cells, above the 2% gate), but the cell-level correlation is 0.03 and the background-matched empirical P is 0.600, so the correlation criterion of the reference framework is **not** met in this tissue.
- **(A) labels every cluster on the plot** as the reference does; the reference's atlas is bone marrow with 11 types, this one is liver with 13.
- **D, E and G are computed inside Cholangiocyte only,** the cell type the pair was selected in. The atlas-wide test is a different and weaker number (Ppargc1a over all 16861 cells gives per-mouse P = 0.63), and is not what these panels show.

## Statistics to quote

| Quantity | Value |
|---|---|
| Cells in the atlas | 16861 |
| Cell types | 13 |
| Most abundant type | B cell (43.81%) |
| Key cell type for the pair | Cholangiocyte (636 cells, 7 mice) |
| Ppargc1a in Cholangiocyte, Control vs STZ | cell-level P = 2e-09; per-mouse t-test (n = 3 vs 4 mice) P = 0.011 |
| Ccnd1 in Cholangiocyte, Control vs STZ | cell-level P = 3.3e-03; per-mouse t-test (n = 3 vs 4 mice) P = 0.18 |
| Detection rate in Cholangiocyte | Ppargc1a 31.7% / 11.4%, Ccnd1 72.6% / 66.7% (Control / STZ) |
| Cells co-expressing the pair | 13.8% |
| Pair correlation in Cholangiocyte | R = 0.03, P = 0.41, empirical P = 0.600 |
