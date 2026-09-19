# Figure 1 — interpretation and defence

Companion to `Fig1_legend.md`. The legend states what the panels show; this states what they mean,
why the analysis was done this way, and how to answer the questions this figure invites.

---

## 1. The claim this figure supports

> In human liver, type 2 diabetes is associated with a coordinated, replicated co-expression
> programme, and the intersection of differential expression, that programme, and curated
> hyperglycaemia biology yields **8 candidate genes** for downstream work.

Note what is *not* claimed: not that these 8 genes are diagnostic, not that they are causal, and not
that the differential expression survives genome-wide multiple-testing correction at a biologically
useful effect size. Those questions belong to Figures 4 and 5, and both are answered honestly there.

---

## 2. Panel by panel

**A — batch correction worked, and it is the honest batch check.**
Before correction the two cohorts separate completely on PC1 (GSE15653 at ≈ −90, GSE64998 at ≈ +95).
After ComBat they overlap. This matters because every later panel pools the two datasets; if the
dominant axis of variation were still cohort, the module–trait correlations in E would be measuring
platform, not disease. Residual structure remains — the GSE64998 ellipse is narrow on PC1 — which is
why the per-cohort replication rule in E exists rather than relying on correction alone.

**B — the differential signal is real but modest.**
96 up, 82 down at nominal P < 0.05 and |log2FC| > 0.5. Eighty-four genes reach FDR < 0.05; thirty-one
of those also clear the fold-change cut. The volcano is wide and shallow rather than tall, which is
what a 16-vs-11 human tissue contrast looks like. The most significant genes are HBB (+1.52,
P = 2.3 × 10⁻⁶), VSNL1 (−0.57, P = 4.1 × 10⁻⁶) and CDKN1A (+0.96, P = 5.3 × 10⁻⁶); the largest fold
change among nominally significant genes is IGFBP1 at −1.66. Note that CDKN1A — the p53 arrest gene
that recurs downstream — is the fourth most significant gene in the contrast and clears FDR at 0.016.

**C — power 7 is a data-driven choice, not a convention.**
Scale-free topology fit crosses 0.80 at power 7 and plateaus; mean connectivity has collapsed to near
zero by then. The reference paper also used 7. That agreement is coincidence — two different tissues
and diseases — but it does mean the network is not built on an unusual parameter.

**D — the network has real modular structure.** Fourteen modules plus the unassigned grey bin, with
two large modules (turquoise 1,796 genes, blue 1,671) and a long tail of small ones.

**E — this is the panel that does the work, and it is stricter than the reference's.**
The strongest module is **red** (r = 0.65 with T2D, P = 2.2 × 10⁻⁴, 135 genes). The Control column is
the exact arithmetic negative of the T2D column because Control is the complement indicator — that is
a property of a two-group design, not a second finding, and the reference's panel has the same
property. Four modules are carried forward (red, yellow, tan, greenyellow) under decision **M11b**.

**F — the red module is internally coherent.** Module membership and gene significance correlate at
0.64 (P = 3.8 × 10⁻¹⁷) across its 135 genes: genes central to the module are also the genes most
associated with T2D. This is the standard WGCNA sanity check, and it passes.

**G — the intersection is 8 genes.** DEGs (178) ∩ key-module genes (418) ∩ hyperglycaemia set (618
in the expression universe) = **8**: VSNL1, SREBF2, PPARGC1A, CCND1, IGFBP1, IRS2, IGF1, IGFBP2.

---

## 3. Why the analysis was done this way

**Nominal P rather than adjusted P (decision M9).** This is the choice most likely to be challenged,
so state it first and state it plainly. 84 genes do reach FDR < 0.05 — the earlier claim that none
did was wrong and has been corrected throughout. But only 31 of those also clear |log2FC| > 0.5, and
carrying that set through the module and gene-set intersections leaves **one gene** (SREBF2,
decision R17). A downstream chain resting on a single gene cannot be interrogated. The panel
therefore uses nominal P with an explicit fold-change floor, the DEG set is labelled **exploratory**
everywhere it appears, and the adjusted-P column is reported in
`results/bulk/DEG_T2D_vs_Control_all.csv` throughout. The defensible position is not "FDR was
impossible" — it is "FDR was applied, reported, and found too restrictive to support a multi-stage
design at this sample size, so the exploratory label is carried forward rather than dropped."

**Cross-cohort replication for key modules (decision M11b).** The reference took the single module
with the strongest positive correlation. This project requires a module to reach P < 0.1 pooled *and*
hold the same direction at P < 0.1 in each cohort separately. The consequence is substantial and
should be volunteered, not conceded: three modules with strong pooled correlations were **excluded**
because they do not replicate — blue (r = −0.72 in GSE15653 vs +0.11 in GSE64998; the sign flips),
green (0.69 vs 0.23) and turquoise (0.63 vs 0.16). Blue and turquoise are the two largest modules, so
this single rule cut the key-module gene pool from **4,022 to 418**. A rule that discards 90 % of the
candidate pool is not a rule chosen to produce a result.

**Why the Dataset column in E is reported but not used as a filter.** It is computed on
ComBat-corrected expression, so it is bounded near zero by construction (max |r| = 0.08). It cannot
detect residual batch effect and is not used to. The per-cohort columns are the real test.

---

## 4. Anticipated questions

**"Your T2D and control groups differ in obesity. How do you separate diabetes from adiposity?"**
You cannot, and the figure does not claim to. All 16 T2D samples are obese and all 11 controls are
lean (decision D1b) — T2D status is *perfectly collinear* with obesity in this design. No statistical
adjustment can separate two perfectly collinear variables. The honest framing is that Figure 1
identifies genes associated with **obese T2D liver versus lean non-diabetic liver**, and that
adiposity is part of the exposure, not a nuisance removed from it. This is recorded in
`docs/decisions_log.md` under D1b and should be stated in the text rather than waited for.

**"What about sex and age?"** Sex is balanced *within* condition (T2D 10 M / 6 F; control 6 M / 5 F)
but confounded with cohort (GSE15653 is 11 F / 3 M, GSE64998 is 13 M / 0 F). The limma design is
`~ condition + dataset`, so the dataset term absorbs most of the sex imbalance. Age is the weaker
point: it is recorded for only 14 of 27 samples, and among those the T2D group is about ten years
older (mean 46.1 vs 36.2). Age is therefore *not* adjusted for, and that is a stated limitation
rather than an oversight.

**"n = 27 is small for WGCNA."** Correct, and this is why the network is not the claim. WGCNA here
is a dimension-reduction and prioritisation step, and its output is subjected to a replication
requirement (M11b) that a small-n artefact would be expected to fail. Panel F's MM–GS correlation is
an internal-consistency check, not independent evidence — both quantities come from the same 27
samples, and the legend says so.

**"Isn't panel G circular? You select on a hyperglycaemia gene set and then claim hyperglycaemia
biology."** Panel G itself is not circular — it is a definition, not a test. The circularity risk
lands on Figure 2, where the same set is tested for enrichment against overlapping annotation, and it
is declared there explicitly. The correct defence is to point to that declaration rather than deny
the overlap.

**"Why does your strongest module have r = 0.65 when the reference reports 0.83?"** Because the
replication rule removed the modules that would have scored higher on the pooled correlation alone.
Blue reaches |r| = 0.72 in one cohort and +0.11 in the other. A higher pooled correlation obtained by
ignoring that sign flip would be a worse result, not a better one.

---

## 5. What this figure does not show

- It does not show that the 8 candidates are specific to T2D rather than to obesity.
- It does not show that they are causal — Figure 5 tests that directly and finds that neither Tier 1
  target is.
- It does not show that they are diagnostic — Figure 4 tests that and finds no model separates groups
  in an independent cohort at conventional confidence.
- It does not establish the direction of any individual gene beyond the two-cohort concordance
  already reported in `DEG_T2D_vs_Control_all.csv`.

The value of the figure is that it produces a **small, replicated, interrogable** candidate set. Its
smallness is the reason every later figure can be checked, and its smallness is also the reason the
later figures are honest about what they cannot support.
