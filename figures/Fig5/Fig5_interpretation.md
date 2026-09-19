# Figure 5 — interpretation and defence

Companion to `Fig5_legend.md`. Exposures: the 9 genes of hsa04152 with usable instruments.

---

## 1. The claim this figure supports

> Within the selected KEGG map, two genes show evidence of a causal effect on type 2 diabetes:
> **SREBF1** (OR 0.903 per SD of genetically predicted expression, 95 % CI 0.872–0.935, q = 8.1 × 10⁻⁸)
> and **SIRT1** (OR 1.027, 95 % CI 1.011–1.044, q = 4.3 × 10⁻³). **Neither of the two Tier 1 targets
> is causal** — PPARGC1A P = 0.16, CCND1 P = 0.55.

This is the project's clearest positive *and* a direct contradiction of its own target nomination.
Both must be presented together. The contradiction is the more important half.

---

## 2. Panel by panel

**A — the whole screen, not a selected row.**
Nine exposures × every applicable estimator, bold where P < 0.05. SREBF1 heads the table
(IVW P = 9.0 × 10⁻⁹); SIRT1 follows with 3 of its 4 estimators significant (IVW 9.5 × 10⁻⁴, weighted
median 6.9 × 10⁻³, weighted mode 2.9 × 10⁻²). PRKAA1 and CAMKK2 show method disagreement (§4).
PPARGC1A — the best-instrumented exposure at 12 SNPs, and a Tier 1 target — is null across all four
estimators.

**B — SIRT1's instruments behave.** SNP effect on outcome rises with SNP effect on exposure, and all
four method lines agree in slope and sign. That concordance is what a clean instrument set looks like.

**C — no single SNP drives SIRT1.** Four per-SNP Wald ratios, all positive, summary rows below the
rule. No instrument sits far from the others.

**D — directionality holds throughout.** Exposure R² exceeds outcome R² for **25 of 25** genes, so
the instruments act on expression first. No exposure shows evidence of reverse causality.

---

## 3. Why the analysis was done this way

**Why nine exposures rather than one.** The reference tested a single gene (PAK1) because Finan et al.
2017 had already nominated it as druggable. There is no equivalent prior nomination here, so testing
one gene would mean choosing it after seeing the data. Every instrumented gene on the selected map is
therefore tested and the whole screen shown. This is more conservative than the reference's design,
not less.

**Why FDR is recomputed over the 9.** The analysis this figure reports is *MR of hsa04152*, so that
is the family the correction belongs to. Benjamini–Hochberg over one primary estimate per gene —
the same rule `scripts/19/20b` use. Over all 24 instrumented genes the same two exposures give
q = 2.2 × 10⁻⁷ and 1.1 × 10⁻². Both sets are in the legend. **Correcting within the map changes the
q values, not which genes clear 0.05** — say this pre-emptively, because recomputing a correction on
a subset looks like q-hacking unless the invariance is stated.

**Why SIRT1 and not SREBF1 in panels B and C.** SREBF1 has the smaller P but only **2 instruments** —
too few for MR-Egger, the weighted median or the weighted mode, so its pleiotropy and heterogeneity
diagnostics would be empty. SIRT1 carries 4 instruments and 4 estimators, Cochran Q P = 0.88, Egger
intercept P = 0.62. Featuring the testable exposure rather than the smallest P is a deliberate choice
and should be presented as one.

**Why LD clumping is the primary analysis (decision R21).** The first pass used distance pruning —
one SNP per 250 kb — which for PRKAA1 and SERPINE1 left a **single lead SNP carrying the whole
estimate**. Once the locus is properly clumped and the other independent signals enter, those
estimates attenuate to null. The earlier reading "MR nominates the AMPK arm" was an artefact of
pruning and is formally withdrawn. Volunteering a withdrawn result is stronger than being asked about
it.

---

## 4. The honest weaknesses, stated first

**SREBF1 rests on 2 instruments.** With 2 SNPs, no pleiotropy test is possible. A q of 8 × 10⁻⁸
looks decisive and overstates the robustness; the correct reading is "a precise estimate from a very
small instrument set, untestable for horizontal pleiotropy". Do not lead with this number.

**PRKAA1 and CAMKK2 show method disagreement.** PRKAA1 is null by IVW (OR 0.865, P = 0.23) but
significant by the weighted median and weighted mode (OR 0.765, P = 2.2 × 10⁻⁵); CAMKK2 is
significant only by the weighted mode. On 5 instruments each, that pattern is what a single outlying
instrument produces. Reported, not interpreted.

**SIRT1's direction is counter-intuitive.** Higher genetically predicted SIRT1 expression associates
with *higher* T2D risk (OR 1.027). Most of the literature treats SIRT1 as protective. The estimate is
reported as it came out. Possible readings — whole-blood eQTL not reflecting hepatic expression,
compensatory upregulation, or a genuine tissue-specific effect — are speculation, and the figure does
not choose between them.

**Exposures are whole-blood eQTLs, not liver.** eQTLGen is whole blood. The phenotype is hepatic. A
gene whose regulation differs between blood and liver will be mis-instrumented. The GTEx liver check
exists (3 genes with usable instruments, single-SNP Wald ratios only) and cannot be tested for
pleiotropy — it is a consistency check, not a replication.

---

## 5. Anticipated questions

**"Your own Tier 1 targets are null. Doesn't that invalidate the project?"** No, and this needs a
clear answer. MR asks whether *lifelong genetically determined expression* affects disease risk.
Figures 1–3 ask whether expression *changes with disease state* and *in which cell type*. These are
different questions with different answers, and a gene can be a disease-responsive marker or a
mechanistic node without being a causal determinant of risk. `docs/target_selection.md` nominates
PPARGC1A and CCND1 on expression and cell-type evidence and states explicitly that MR does not
support them. The chain audit records that **no gene passes all five figures**. Reporting that
plainly is the integrity of the work, not a hole in it.

**"Then why does the study nominate PPARGC1A and CCND1 at all?"** Because the claim is a *mechanism
supported by expression and cell-type evidence*, not a validated causal or diagnostic target. That is
the sentence to use. A causal claim would require MR support, and it is absent.

**"Should SREBF1 or SIRT1 become the target instead?"** They are supported on this one axis only —
no cell-type confirmation, not in the Fig 1 candidate set. `docs/target_selection.md` lists both
under "not targets" for that reason. Promoting a gene on a single line of evidence is the error this
project is built to avoid, and switching targets after seeing the MR result would be exactly that
error.

**"216 cases would make this uninterpretable."** That is the reference paper's problem, not this
one's — its outcome GWAS has 216 RA cases (`docs/reference_paper_review.md` §3.2). This analysis uses
GCST006867 (Xue et al. 2018), **62,892 T2D cases and 596,424 controls**, European. Raise the
comparison unprompted; a roughly 290-fold larger case count is a genuine methodological advantage
over the framework being reproduced.

**"Why not test all 122 genes on the map?"** Only 9 have instruments meeting the criteria — cis-eQTL
P < 5 × 10⁻⁸, F > 10, surviving LD clumping. The other 113 have no usable instrument. That is an
instrument-availability constraint, not a selection.

---

## 6. What this figure does not show

- It does not show that SREBF1 or SIRT1 are druggable, or that modulating them would alter disease.
- It does not show that PPARGC1A and CCND1 are *not* mechanistically involved — only that their
  genetically predicted expression does not affect T2D risk.
- It does not establish tissue-specific causal effects; whole-blood eQTLs are the exposure.
- It does not support SREBF1's estimate against horizontal pleiotropy, which 2 instruments cannot test.
