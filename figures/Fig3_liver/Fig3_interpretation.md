# Figure 3 (liver) — interpretation and defence

Companion to `Fig3_legend.md`. Pair: **Ppargc1a + Ccnd1**, both nodes of the selected map hsa04152.

---

## 1. The claim this figure supports

> In mouse liver, the AMPK-axis anchor **Ppargc1a falls in cholangiocytes under hyperglycaemic
> stress, in the same direction as human T2D liver**. The pair Ppargc1a–Ccnd1 is co-expressed in that
> cell type but is **not correlated above a detection-matched background**, so the pair-level
> criterion of the reference framework is not met.

This figure reports a **partial confirmation and an explicit negative**. Both halves are the result.
Do not present it as a positive that happens to have a caveat.

---

## 2. Panel by panel

**A — the atlas.** 16,861 cells, 13 annotated cell types. B cells dominate at 43.8 %; hepatocytes are
0.41 % (see §5). **B** — canonical markers separate the types cleanly. **C** — every cell type
contains both groups, so no cluster is a single-animal artefact.

**D — Ppargc1a is down in STZ cholangiocytes, and it is significant on the honest test.**
Per-mouse t-test P = 0.011 (n = 3 vs 4 mice). Detection falls from **31.7 % to 11.4 %** of
cholangiocytes. Human direction: PPARGC1A log2FC −0.61, P = 6 × 10⁻⁴. **Mouse and human agree.**

**E — Ccnd1 is not significant on the per-mouse test.** P = 0.18, detection 72.6 % → 66.7 %. The
cell-level Wilcoxon gives P = 3.3 × 10⁻³, and that discrepancy is exactly why the per-mouse test is
the one on the panel (§3).

**F — the pair is genuinely co-localized.** 13.8 % of cholangiocytes co-express both genes, above
the framework's 2 % gate. The joint-density peak sits on the cholangiocyte cluster.

**G — and they are not correlated.** R = 0.03, raw P = 0.41, and against 2,000 detection-matched
random gene pairs the **empirical P = 0.600**. The pair is no more co-varying than two random genes
with the same detection rates.

---

## 3. Why the analysis was done this way

**The bracket is the per-mouse t-test, not the cell-level Wilcoxon — and this is the single most
defensible choice in the figure.** Cells within an animal share genotype, cage, dissection, capture
lane and batch. They are not independent replicates of a Control-versus-STZ contrast. A cell-level
test on 636 cells answers "how many cells were sequenced", not "how different are the groups". The
reference paper reported cell-level four-star P-values from **one mouse per group** — a design in
which no per-animal test is even possible (`docs/reference_paper_review.md` §3.1). This project has
3 vs 4 mice, so the per-animal test exists and is used. Both numbers are tabulated in the legend;
`SIG_UNIT <- "cell"` in `scripts/35_fig3_paper_style.R` switches the panel back for anyone who wants
to see the anticonservative version.

Expect the follow-up: *"Isn't n = 3 vs 4 underpowered?"* Yes. That is why Ccnd1 reads `ns` here at
P = 0.18 despite a cell-level P of 3.3 × 10⁻³. Under-powering produces false negatives, not false
positives — the direction of the error is conservative, and the conservative direction is the right
one to err in when the alternative is pseudoreplication.

**The correlation is judged against a detection-matched null, not a raw P.** Cell-level correlations
in droplet data are inflated by shared detection rates: two genes that are both often-zero will
correlate simply because their zeros coincide. Comparing R against 2,000 random pairs matched on
detection rate removes that artefact. Reporting only the raw Pearson P would have made this pair look
better (P = 0.41 is still null, but a larger R with a small raw P would not have been caught).

**All of D, E and G are computed inside cholangiocytes only.** The atlas-wide number is different and
weaker (Ppargc1a per-mouse P = 0.63 over all 16,861 cells) and is not what the panels show. Cell-type
restriction is the reference's own design — its panels are neutrophil-only.

---

## 4. The negative is predicted, not anomalous

Figure 2 panel D shows that on hsa04152 the AMPK subunit transcript **rises** while every output it
activates **falls**. That is the signature of post-translational suppression: AMPK activity is set by
Thr172 phosphorylation, not by transcript abundance. **A pathway regulated at the protein level is
not expected to show transcript co-expression between its nodes.** Figure 3's negative is therefore a
prediction of Figure 2, and the two figures are consistent.

`docs/pathway_selection.md` §3.3 gives three structural reasons, all of which apply here:

1. **Post-translational control** — as above, shown directly by the discordant edges.
2. **Missing parenchyma** — this atlas contains only **69 hepatocytes** (0–5 per control mouse). The
   cells that carry hepatic gluconeogenesis and lipogenesis, where the AMPK axis would be most
   visible, are essentially absent. The compartment is missing, not negative.
3. **Model mismatch** — GSE244475 is **STZ**, i.e. β-cell ablation and insulin *deficiency*
   (C57BL/6J males, 3 Control vs 4 STZ). The human cohort is **obese, insulin-resistant** T2D. The
   lipogenic arm is expected to invert under this mismatch, and it does: Srebf1, Fasn and Fbp1 fall
   in STZ cholangiocytes although they rise in human T2D liver.

---

## 5. Anticipated questions

**"Your key cell type is cholangiocytes. Why not hepatocytes?"** Because there are 69 hepatocytes in
the atlas, 0–5 per control mouse — too few for any per-mouse test. Cholangiocytes were selected as
the cell type in which the pair could actually be evaluated, and that selection is stated. This is a
limitation of the public dataset, not a choice made to obtain a result. Whole-tissue qPCR is the
appropriate assay for the hepatocyte compartment and is in the validation plan.

**"Ccnd1 is `ns` here but significant in kidney. Isn't that inconsistent?"** No — it is
tissue-specific, which is what a cell-type-resolved analysis is for. In kidney endothelium Ccnd1 is
up at per-mouse P = 2.9 × 10⁻⁴ in the human direction. In liver cholangiocytes it does not move.
Reporting both is the point of running two tissues.

**"You changed the pair from the version that passed."** Yes, deliberately, and the reasoning is in
`figures/README.md`. The pair *Cdkn1a–Ccnd1* clears all three gates in kidney endothelium (empirical
P = 0.009), but **Cdkn1a is not a node of hsa04152** — verified against the KGML by Entrez ID. Using
it would mean the single-cell figure drew from a different pathway than the one Figure 2 selects,
which is precisely the chain break the reference avoids. The on-pathway pair is shown even though it
fails; the off-pathway pair is reported in `docs/pathway_selection.md` §3.2b as a downstream readout.
Choosing the figure that fails, over the figure that passes, because only one of them is on-pathway,
is the strongest single defence available here.

**"A null result from an underpowered experiment is uninformative."** For gate 1 (changed) that
caution applies, and Ccnd1's `ns` should not be over-read. For gate 3 (correlated) it does not: the
correlation is computed on **636 cells**, not 7 mice, and R = 0.03 against a background whose 95th
percentile is 0.153 is a well-estimated null, not an absence of power.

---

## 6. What this figure does not show

- It does not show that the AMPK axis is inactive in mouse liver — only that its two anchors do not
  co-vary transcriptionally in the one cell type where both are measurable.
- It does not show anything about hepatocytes.
- It does not test AMPK activity, which requires p-AMPKα Thr172 by Western blot.
- It does not establish that the mouse model reproduces human obese T2D; STZ does not, and the
  lipogenic arm inverts accordingly.
