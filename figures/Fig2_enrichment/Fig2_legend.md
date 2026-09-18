# Figure 2 legend - functional enrichment of the candidate set

**Fig. 2 | Pathways enriched in the T2D candidate genes.**
**a**, GO over-representation of the 19 candidate genes from Fig. 1g, top 12 terms per ontology by P. Point size, genes from the candidate set annotated to the term; colour, BH-adjusted q.
**b**, KEGG over-representation, all 9 pathways at q < 0.05.
**c**, Enrichment map of the 20 most significant GO biological-process terms; nodes are terms, edges join terms sharing at least a quarter of their candidate genes (Jaccard >= 0.25).
**d**, GSEA of the full T2D vs control ranking against the Hallmark collection, sets at q < 0.25 (29 of 50). Bars, normalised enrichment score; points, set size.

## Caveats that belong in the text, not the figure

- **a**,**b** test 19 genes against a 6259-gene background; with a set this small, an over-representation q value is fragile and a single gene changes several terms.
- The candidate set was defined partly by a hyperglycaemia gene set (Fig. 1g), which shares annotation with the GO and KEGG terms tested here, so **a** and **b** are not independent of that selection step.
- **c** edges are gene-overlap similarity, not curated pathway relationships.

## Statistics to quote

| Quantity | Value |
|---|---|
| Candidate genes tested | 19 |
| GO terms at q < 0.05 | 230 |
| Top GO term | Response to insulin (P = 5.7 × 10⁻¹⁰) |
| KEGG pathways at q < 0.05 | 9 |
| Top KEGG pathway | AMPK signaling pathway (P = 2.0 × 10⁻⁶) |
| Hallmark sets at q < 0.25 | 29 of 50 |
| Strongest Hallmark NES | Cholesterol homeostasis, NES = 2.19 (q = 0.000) |
