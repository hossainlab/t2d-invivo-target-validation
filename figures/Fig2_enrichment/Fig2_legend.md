# Figure 2 legend - functional enrichment of the candidate set

**Fig. 2 | Pathways enriched in the T2D candidate genes.**
**a**, GO over-representation of the 8 candidate genes from Fig. 1g, top 12 terms per ontology by P. Point size, genes from the candidate set annotated to the term; colour, BH-adjusted q.
**b**, KEGG over-representation, all 21 pathways at q < 0.05.
**c**, Pathway interaction network of the 30 most enriched KEGG terms; nodes are pathways, edges join pathways sharing at least half of their candidate genes (Jaccard >= 0.5). Node colour, BH-adjusted q; size, candidate genes in the pathway.
**d**, GSEA of the full T2D vs control ranking against the Hallmark collection, sets at q < 0.25 (29 of 50). Bars, normalised enrichment score; points, set size.

## Caveats that belong in the text, not the figure

- **a**,**b** test 8 genes against a 6259-gene background; with a set this small, an over-representation q value is fragile and a single gene changes several terms.
- The candidate set was defined partly by a hyperglycaemia gene set (Fig. 1g), which shares annotation with the GO and KEGG terms tested here, so **a** and **b** are not independent of that selection step.
- **c** edges are gene-overlap similarity between pathways, not curated pathway-pathway relationships; with only 19 candidate genes most KEGG pathways share the same few genes, which is why the overlap threshold is high.

## Statistics to quote

| Quantity | Value |
|---|---|
| Candidate genes tested | 8 |
| GO terms at q < 0.05 | 205 |
| Top GO term | Insulin-like growth factor receptor signaling pathway (P = 9.6 × 10⁻⁹) |
| KEGG pathways at q < 0.05 | 21 |
| Top KEGG pathway | AMPK signaling pathway (P = 7.2 × 10⁻⁸) |
| Hallmark sets at q < 0.25 | 29 of 50 |
| Strongest Hallmark NES | Cholesterol homeostasis, NES = 2.19 (q = 0.000) |
