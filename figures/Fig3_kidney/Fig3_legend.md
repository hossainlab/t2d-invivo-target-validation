# Figure 3, reference-paper style (kidney) - legend

Drawn to match Xu et al., *Phytomedicine* 154 (2026) 158050, Fig. 3, with this project's data.
Pair: **Ppargc1a and Ccnd1**, both nodes of the selected KEGG map hsa04152.

**Fig. 3. Results of scRNA-seq analysis.** (A) UMAP projection illustrating the single-cell atlas of mouse kidney tissue. (B) Bubble plot showing marker genes for each cell type. (C) Bar plot displaying the proportion and number of cells from each sample source. (D-E) Analysis of Ppargc1a and Ccnd1 expression levels in endothelial across experimental groups. Per-mouse t-test (n = 3 vs 4 mice); ns: P > 0.05, *: P < 0.05, **: P < 0.01, ***: P < 0.001, ****: P < 0.0001. (F) UMAP-based heatmap showing co-localization of Ppargc1a and Ccnd1 across the atlas. (G) Correlation between Ppargc1a and Ccnd1 expression levels in endothelial.

## Where this data forces a deviation from the reference

- **Two groups, not three.** The reference contrasts Ctrl, CIA and Treat, so its D and E carry three brackets. This atlas is Control vs STZ, so each panel carries one.
- **The bracket is the per-mouse t-test (n = 3 vs 4 mice), not the cell-level Wilcoxon.** Cells within a mouse are not independent replicates, so a cell-level P is anticonservative. Both are tabulated below; the reference reported the cell-level test.
- **(F) fails before (G) is reached.** Only 0.17% of endothelial cells co-express the pair, below the 2% co-localization gate of the reference framework, so the correlation that follows (R -0.01, empirical P 0.810) is not interpretable here. Ppargc1a is barely expressed in this cell type.
- **(A) labels every cluster on the plot** as the reference does; the reference's atlas is bone marrow with 11 types, this one is kidney with 14.
- **D, E and G are computed inside Endothelial only,** the cell type the pair was selected in. The atlas-wide test is a different and weaker number (Ppargc1a over all 6224 cells gives per-mouse P = 0.048), and is not what these panels show.

## Statistics to quote

| Quantity | Value |
|---|---|
| Cells in the atlas | 6224 |
| Cell types | 14 |
| Most abundant type | Endothelial (46.19%) |
| Key cell type for the pair | Endothelial (2875 cells, 7 mice) |
| Ppargc1a in Endothelial, Control vs STZ | cell-level P = 5.8e-01; per-mouse t-test (n = 3 vs 4 mice) P = 0.96 |
| Ccnd1 in Endothelial, Control vs STZ | cell-level P = 0e+00; per-mouse t-test (n = 3 vs 4 mice) P = 0.00029 |
| Detection rate in Endothelial | Ppargc1a 0.5% / 0.4%, Ccnd1 27.1% / 40.4% (Control / STZ) |
| Cells co-expressing the pair | 0.2% |
| Pair correlation in Endothelial | R = -0.01, P = 0.5, empirical P = 0.810 |
