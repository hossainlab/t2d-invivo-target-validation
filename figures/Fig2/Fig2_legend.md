# Figure 2, reference-paper style - legend

Drawn to match Xu et al., *Phytomedicine* 154 (2026) 158050, Fig. 2, with this project's data.

**Fig. 2. Functional enrichment analysis of intersecting genes.** (A) Results of GO enrichment analysis. Blue represents biological process (BP) terms, yellow represents cellular component (CC) terms, and purple represents molecular function (MF) terms. (B) Bubble chart displaying the results of KEGG analysis. (C) Pathway interaction network of the top 30 enriched KEGG terms. (D) Pathway diagram generated using the pathview package, illustrating signal activation within AMPK signaling pathway.

## Where this data forces a deviation from the reference

- **(A) has no CC column.** Over-representation of the 8 candidate genes returns BP and MF terms only; not one cellular-component term is returned, so the paper's yellow column has nothing to draw. The colour key is kept as the reference states it.
- **(B) and (C) rest on 4 genes, not 159.** Only 4 of the 8 candidates (CCND1, IGF1, IRS2, PPARGC1A) map into any KEGG pathway, so every GeneRatio on this figure has denominator 4.
- **The 21 pathways at adjusted P < 0.05 are not 21 independent findings.** They are recombinations of those same 4 genes: AMPK signaling pathway has all 4, and 18 of the 21 rest on a single pair. The cancer maps in the list are annotation artefacts of CCND1, not evidence of cancer biology. Read the panel as one signal seen through 21 annotation sets.
- **(C) shows the top 30 KEGG terms by nominal P,** as the reference does, because only 21 reach adjusted P < 0.05; significance is carried by node colour. Edges join pathways sharing at least half of their candidate genes (Jaccard >= 0.5), which is gene-overlap similarity, not a curated pathway-pathway relation.
- **(D) shows AMPK signaling pathway (hsa04152)**, the map selected in `docs/pathway_selection.md`, standing in for the paper's regulation of actin cytoskeleton. Node colour is the T2D vs control log2 fold change, so the panel shows the *transcript* state of the map, which for this pathway is the finding: the AMPK subunit transcript rises while its outputs fall.
- The candidate set was defined partly by a hyperglycaemia gene set (Fig. 1G), which shares annotation with the GO and KEGG terms tested here, so (A) and (B) are not independent of that selection step.

## Statistics to quote

| Quantity | Value |
|---|---|
| Candidate genes tested | 8 |
| Background genes | 11209 (GO), 6259 (KEGG) |
| GO terms at adjusted P < 0.05 | 205 (BP and MF) |
| Top GO term | insulin-like growth factor receptor signaling pathway, adjusted P = 6.7e-06 |
| KEGG pathways at adjusted P < 0.05 | 21 |
| Selected pathway | AMPK signaling pathway (hsa04152), GeneRatio 4/4, adjusted P = 5.8e-06 |
| Genes driving the KEGG result | CCND1, IGF1, IRS2, PPARGC1A |
