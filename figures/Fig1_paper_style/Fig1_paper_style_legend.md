# Figure 1, reference-paper style - legend

Drawn to match Xu et al., *Phytomedicine* 154 (2026) 158050, Fig. 1, with this project's data.
The Nature/Cell-contract version of the same figure is in `figures/Fig1_bulk_DEG_WGCNA/`.
**Ship one or the other, not both.**

**Fig. 1. Screening of DEGs related to T2D in human liver.** (A) PCA plot before and after batch effect correction. (B) Volcano plot of DEGs. Blue indicates significantly downregulated genes, red indicates significantly upregulated genes, and grey indicates non-significant genes. (C) Scale independence and mean connectivity in the weighted gene co-expression network. At a soft threshold of 7, the network exhibits high scale independence and mean connectivity close to zero. (D) Different co-expression modules within the weighted gene co-expression network. (E) Heatmap showing the correlation of 14 modules and the unassigned grey bin with T2D case and control groups. (F) Correlation between MM and GS for all genes in the red module. (G) Venn diagram illustrating the intersection among DEGs, key-module genes, and hyperglycaemia-related genes.

## Where this data forces a deviation from the reference

- **(B) uses nominal P, not adjusted P** (decision M9). The paper reported 1,648 DEGs at adjusted P < 0.05. Here 84 genes reach FDR < 0.05 and only 31 of those also clear |log2FC| > 0.5; carrying that set through the paper's strict module and gene-set rules leaves a single intersecting gene (SREBF2, decision R17), which no downstream analysis can rest on. The panel therefore shows nominal P < 0.05 and |log2FC| > 0.5: 96 up, 82 down, labelled exploratory. The adjusted-P column is reported in `DEG_T2D_vs_Control_all.csv` throughout.
- **(C) reference line at 0.80,** not the paper's 0.90; 0.80 is the fit this project selects power at. Power 7 happens to match the paper.
- **(E) has 14 modules plus the unassigned grey bin, not seven,** and the strongest is **red** (r = 0.65) rather than the paper's blue (r = 0.83). Control is the complement of T2D, so its column is the exact negative with an identical P, exactly as in the reference. Grey is drawn last, as in the reference, and is display-only: it is the unassigned-gene bin, so script 03 excludes it from the key-module rule.
- **(F) shows the red module** and this project's hub thresholds (|kME| > 0.7, |GS| > 0.2, decision M10), not the paper's 0.7 / 0.7. Points take the module colour as in the reference; because that colour is red here, the two threshold lines are drawn black and dashed rather than the reference's red, which would be invisible against the points.
- **(G)** the third set is hyperglycaemia-related genes (618 in the expression universe), standing in for the paper's neutrophil-migration set. The three-way intersection is **8 genes**.
- Key modules follow decision **M11b** (T2D P < 0.1 pooled and replicated in both cohorts), which is stricter than the paper's single-best-module rule.

## Statistics to quote

| Quantity | Value |
|---|---|
| Samples | 16 T2D, 11 lean control (GSE15653 + GSE64998) |
| DEGs (nominal P < 0.05, \|log2FC\| > 0.5) | 96 up, 82 down |
| Genes at FDR < 0.05 | 84 |
| Genes at FDR < 0.05 and \|log2FC\| > 0.5 | 31 |
| Soft-thresholding power | 7 |
| Modules (excluding grey) | 14 |
| Key modules | red, yellow, tan, greenyellow |
| Strongest module | red, n = 135 genes, r = 0.65, P = 2.2e-04 |
| MM vs GS in the red module (panel F) | cor = 0.64, P = 3.8e-17 |
| Three-way intersection | 8 genes |
