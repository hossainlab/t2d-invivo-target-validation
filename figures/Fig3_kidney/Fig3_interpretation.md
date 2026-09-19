# Figure 3 (kidney) — interpretation and defence

Companion to `Fig3_legend.md`. Pair: **Ppargc1a + Ccnd1**, both nodes of the selected map hsa04152.

---

## 1. The claim this figure supports

> In mouse kidney, the AMPK-axis output node **Ccnd1 rises specifically in endothelium under
> hyperglycaemic stress, in the same direction as human T2D liver (log2FC +0.80)**. The second
> anchor, Ppargc1a, is **not expressed** in that cell type, so the pair fails the co-localization gate
> before the correlation criterion can be reached.

This figure reports **one clear confirmation and one gate failure with a mechanical cause**. The gate
failure is an expression-coverage fact, not a biological negative, and the two must not be conflated.

---

## 2. Panel by panel

**A — the atlas.** 6,224 cells, 14 annotated cell types, endothelium the most abundant at 46.2 %.
**B** — markers separate the types. **C** — both groups are present in every cluster.

**D — Ppargc1a does not move, because it is barely there.** Per-mouse P = 0.96; detection **0.5 % in
control and 0.4 % in STZ** endothelial cells. This is an absence of measurement, not an absence of
effect. Read the detection rate before reading the P-value.

**E — Ccnd1 is up in STZ endothelium, strongly and on the honest test.** Per-mouse t-test
P = 2.9 × 10⁻⁴ (n = 3 vs 4 mice); detection rises from **27.1 % to 40.4 %**. Human direction: CCND1
log2FC +0.80, P = 7.5 × 10⁻⁴, and it replicates in the external cohort GSE23343 (P = 0.016). **Mouse
and human agree, in the cell type the candidate score also localises to.**

**F — the pair is not co-localized.** Only **0.17 %** of endothelial cells co-express both genes,
far below the framework's 2 % gate. The joint-density map is dominated by the Ccnd1 signal alone.

**G — the correlation is therefore not interpretable.** R = −0.01, empirical P = 0.810. With one
partner detected in 0.4 % of cells, there is almost no cell in which both are measured. The panel is
shown for completeness and because the reference's figure has a G; the legend states that F fails
before G is reached.

---

## 3. Why the analysis was done this way

**The bracket is the per-mouse t-test.** Cells within an animal are not independent replicates of a
Control-versus-STZ contrast; a cell-level P measures sequencing depth. The reference reported
cell-level four-star P-values from **one mouse per group**, where no per-animal test is possible
(`docs/reference_paper_review.md` §3.1). Here there are 3 vs 4 mice, so the per-animal test exists
and is used. Note that in this tissue the honest test is *not* the weaker one — Ccnd1 reaches
P = 2.9 × 10⁻⁴ per mouse, so the result survives the conservative analysis rather than depending on
the anticonservative one. That is worth saying explicitly: the same rule that turned Ccnd1 `ns` in
liver leaves it highly significant here.

**The correlation is judged against a detection-matched null.** 2,000 random gene pairs matched on
detection rate. This guards against the droplet-data artefact in which two often-zero genes correlate
because their zeros coincide. It also means the null in panel G is correctly calibrated for a gene
detected in 0.4 % of cells — which is why the empirical P (0.810) is uninformative rather than
falsely reassuring.

**D, E and G are computed inside endothelium only.** The atlas-wide comparison is a different
question and gives a different answer (Ppargc1a per-mouse P = 0.048 over all 6,224 cells, driven by
other cell types). Cell-type restriction matches the reference's own neutrophil-only design.

---

## 4. Reading the two tissues together

Neither tissue satisfies the full three-gate criterion, but they fail at **different gates and for
different reasons**, and each contributes one confirmed anchor:

| | Liver (cholangiocyte) | Kidney (endothelium) |
|---|---|---|
| Ppargc1a changed, human direction | **yes**, per-mouse P = 0.011 | no — detected in 0.4 % of cells |
| Ccnd1 changed, human direction | no, per-mouse P = 0.18 | **yes**, per-mouse P = 2.9 × 10⁻⁴ |
| Co-localized (≥ 2 %) | **yes**, 13.8 % | **no**, 0.17 % |
| Correlated vs matched null | **no**, empirical P = 0.600 | not interpretable |

The defensible summary: **both anchors of the AMPK axis are confirmed in the mouse atlas, in the
human direction, each in its own tissue and cell type — but they are never confirmed together in the
same cells.** That is a weaker statement than the reference makes about RAC1/PAK1, and it is the
statement the data supports.

`docs/pathway_selection.md` §3.3 gives the structural reasons: post-translational control of AMPK
(shown directly in Fig 2D), missing parenchyma (this atlas has only 150/174 proximal-tubule cells,
the compartment where renal AMPK signalling would be most visible), and model mismatch — GSE244475 is
**STZ**, insulin-*deficient*, against a human cohort that is obese and insulin-*resistant*.

---

## 5. Anticipated questions

**"Panel G is empty of signal. Why show it?"** Because removing it would hide the gate failure. The
reference's framework has four pair-level panels (D, E, F, G) and this project reproduces all four,
including the one that fails. Showing G with its 0.17 % co-expression stated is more honest than
dropping the panel and reporting only the panels that worked.

**"Ppargc1a isn't expressed in endothelium — so why test it there?"** Because the cell type was
chosen by the pathway-score analysis (the human-up candidate score rises specifically in endothelium,
P = 0.0018, BH 0.011) and by Ccnd1's localisation, before the pair was evaluated. Testing the pair in
the cell type the score nominates, and reporting that one partner is absent there, is the correct
order of operations. Choosing the cell type *after* seeing which one gave the best correlation would
be the error.

**"So is CCND1 confirmed or not?"** Confirmed as a single gene, in this tissue: up in STZ endothelium
at per-mouse P = 2.9 × 10⁻⁴, in the human direction, replicating in an external human cohort. Not
confirmed as half of a co-expressed pair — that requires its partner, and its partner is absent here.
`docs/target_selection.md` places CCND1 in Tier 1 on exactly this basis: expression and cell-type
evidence, explicitly *not* causal evidence (Figure 5 finds CCND1 null at P = 0.55).

**"Isn't the STZ model simply wrong for T2D?"** It is a mismatch, stated as one. STZ models
hyperglycaemic stress, which is the component this project's p53/arrest and AMPK-output findings
depend on; it does not model obesity-driven insulin resistance, which is why the lipogenic arm
inverts. The appropriate claim is about hyperglycaemic stress, not about T2D as a whole, and the
in vivo validation plan uses a GV1 Ob-T2D model rather than STZ for that reason.

---

## 6. What this figure does not show

- It does not show that Ppargc1a is absent from kidney — only that it is not detected in endothelium
  at droplet depth.
- It does not show pair-level co-regulation; the co-localization gate fails.
- It does not show that Ccnd1's rise is AMPK-dependent. That requires perturbation, and the drawn
  route (AMPK ⊣ ELAVL1 → CCND1) is a two-step inference from the KEGG map, not a measurement.
- It does not test AMPK activity.
