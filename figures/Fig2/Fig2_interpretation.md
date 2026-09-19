# Figure 2 — interpretation and defence

Companion to `Fig2_legend.md`.

---

## 1. The claim this figure supports

> The 8 candidate genes are not a scattered list. Four of them converge on a single KEGG map — **AMPK
> signalling (hsa04152)** — and the transcript state of that map in human T2D liver has a specific,
> interpretable shape: the kinase subunit rises while the outputs it activates fall.

What is *not* claimed: that 21 independent pathways are dysregulated, that the cancer maps mean
anything about cancer, or that pathway enrichment on 8 genes constitutes evidence independent of how
those 8 genes were selected.

---

## 2. Panel by panel

**A — GO enrichment, BP and MF only.**
205 terms reach adjusted P < 0.05. The top term is *insulin-like growth factor receptor signalling
pathway* (adjusted P = 6.7 × 10⁻⁶). The terms cluster around insulin/IGF signalling and carbohydrate
metabolism — i.e. the phenotype. There is **no CC column** because over-representation of 8 genes
returns no cellular-component term at all; the yellow key is retained only because the figure
reproduces the reference's colour scheme.

**B — KEGG enrichment, and the selection.**
21 pathways at adjusted P < 0.05. **AMPK signalling is rank 1** (adjusted P = 5.8 × 10⁻⁶, GeneRatio
4/4) and is boxed in red. It is the only map that contains *all four* KEGG-annotated candidates. The
next best map is FoxO signalling at 6.8 × 10⁻⁴ — roughly 117-fold weaker.

**C — the pathway network shows the redundancy, not extra evidence.**
The top 30 KEGG terms, with edges where two pathways share ≥ 50 % of their candidate genes. The dense
central cluster is the point: these maps are not independent findings, they are the same few genes
re-annotated.

**D — the pathview map is where the biology actually is.**
hsa04152 coloured by T2D-vs-control log2 fold change. The pattern to read out: **PRKAA2 (the AMPK α2
subunit transcript) is up (+0.29, P = 0.017) while PPARGC1A is down (−0.61, P = 6 × 10⁻⁴), SIRT1 is
down (−0.39, P = 0.0044) and PFKFB3 is down (−0.66, FDR 0.042); meanwhile the outputs AMPK represses
rise — FBP1 +0.62, SREBF1 +0.44, FASN +0.70, CCND1 +0.80.**

---

## 3. The interpretation that matters

Every AMPK output moves *against* the direction the drawn activation edges predict, and every AMPK
target of repression moves *with* de-repression. A drawn activation edge whose kinase subunit rises
while its target falls is the transcriptional signature of a pathway that has been switched off
**post-translationally** — AMPK activity is set by Thr172 phosphorylation and the AMP/ATP ratio, not
by subunit transcript abundance.

This is the single most important interpretive claim in the project, and it has two consequences that
must be owned rather than discovered by an examiner:

1. **It explains why Figure 3 is negative.** If the pathway is regulated at the protein level, single-
   cell *transcript* co-expression is not expected to confirm it. The reference paper's pathway
   (actin cytoskeleton, RAC/PAK) happened to be transcriptionally visible; this one is not. Figure 3's
   negative result is therefore a prediction of Figure 2, not a contradiction of it.
2. **It dictates the validation design.** The in vivo experiment must measure p-AMPKα Thr172 / total
   AMPKα by Western blot, with metformin as a pharmacological positive control. Transcript assays
   alone cannot test this mechanism. That is recorded in `docs/invivo_validation_plan.md`.

---

## 4. Why the analysis was done this way

**Why AMPK and not one of the other 20 maps.** The selection is scripted (`scripts/16_pathway_selection.R`),
not eyeballed, and scores six criteria: ORA, phenotype-set enrichment, pathway activity (`limma::fry`),
GSEA, coherent drawn edges, and network centrality. This is a deliberate improvement on the reference,
which chose *regulation of actin cytoskeleton* — a term sitting around 17th of ~40 in its own
Figure 2B — because it matched the phenotype by eye.

**Be ready to concede one point here.** On the current 8-gene candidate set, AMPK ranks **1st on ORA
and 1st on GSEA**, but **3rd on the six-criterion composite**, behind Insulin resistance (hsa04931)
and Insulin signalling (hsa04910). `docs/pathway_selection.md` §2.1 still reports the pre-M11b
ranking in which AMPK was 1st overall; that table is stale and is flagged in
`docs/reference_paper_review.md` §6. The defensible position: AMPK wins outright on the reference's
own criterion (ORA, its Fig 2B), wins on GSEA, is the only map hitting all four annotated candidates,
and is the top-ranked *signalling* map — hsa04931 is a KEGG **disease** map, not a signalling map.
Say this before being asked; do not let the stale table be found first.

**Why the KEGG panels rest on 4 genes.** Only IRS2, PPARGC1A, CCND1 and IGF1 appear in any KEGG
pathway. The other four candidates (VSNL1, SREBF2, IGFBP1, IGFBP2) are in **no** KEGG map at all —
that is a property of KEGG's coverage, not a filtering decision. Every GeneRatio on the figure
therefore has denominator 4, and the legend states it.

---

## 5. Anticipated questions

**"Pathway enrichment on 8 genes is meaningless."** It is fragile, and the figure says so. The
defence is threefold. First, the effect size is not marginal: 4/4 of the annotated candidates on one
map at adjusted P = 5.8 × 10⁻⁶ against a 6,259-gene background. Second, the enrichment is not the
claim — the *transcript pattern on the map* (panel D) is, and that pattern involves 39 differentially
expressed map genes, not 4. Third, the redundancy is shown rather than hidden: panel C exists
precisely to make the non-independence visible.

**"You report 21 significant pathways. Isn't that inflated?"** Yes, and the legend says so in those
words. All 21 are recombinations of the same four genes; 18 of them rest on a single *pair*. The
cancer maps (melanoma, glioma, prostate, breast, proteoglycans in cancer) are all CCND1 + IGF1 — an
annotation artefact of cyclin D1's presence in every proliferation map, not evidence of cancer
biology. Read the panel as **one signal seen through 21 annotation sets**.

**"The enrichment is circular — you selected on a hyperglycaemia gene set."** Partly true, and
declared. The candidate set was defined in Fig 1G partly by curated hyperglycaemia/insulin-resistance
genes, which share annotation with the GO and KEGG terms tested here. So panels A and B are **not
independent** of the selection step. Two things limit the damage: the hyperglycaemia set was built
from CTD/GO/HPO/KEGG sources spanning many pathways, not from AMPK signalling specifically; and
panel D — the transcript pattern on the map — uses all 39 DE map genes, most of which were never in
the candidate set and never in the hyperglycaemia set. The circularity is real for A and B, and does
not reach D.

**"Why is the pathview panel showing a pathway you say is inactive?"** Because that is the finding.
The panel is not asserting activation; it shows the measured fold changes on the map, and the
discordance between the subunit and its outputs is what identifies post-translational suppression.
The figure legend explicitly says the panel "shows the *transcript* state of the map".

---

## 6. What this figure does not show

- It does not show pathway *activity*. It shows transcript abundance on a pathway map. Activity
  requires phospho-protein measurement, which is the point of the validation plan.
- It does not show that AMPK signalling is uniquely dysregulated — Insulin resistance and Insulin
  signalling score comparably on the composite, and they share genes with AMPK.
- It does not show that the 21 pathways are 21 findings.
- It does not provide evidence independent of Figure 1's gene selection for panels A and B.
