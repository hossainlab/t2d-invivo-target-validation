# Critical review of the reference paper

**Paper.** B. Xu, Xiangqian Xue, Peijun Cheng, Xiaona Guo, Leilei Liu, Yingting Yang, Jianfei Sun,
Ning Wu. *Sidaxue Ameliorates Rheumatoid Arthritis by Inhibiting Synovial Neutrophil Migration:
Insights from an Integrated Multi-omics Analysis.* **Phytomedicine 154 (2026) 158050.**

**Why this review exists.** This project reproduces the paper's analytical framework panel for panel
(`figures/README.md`). Copying a framework means inheriting its weaknesses unless they are named.
This is an appraisal of the source, and a note on which of its problems this project has already
avoided and which it still shares.

**Scope.** The bioinformatics half (Figs 1–4) is reviewed in detail, because that is what this
project mirrors. The chemistry and wet-lab half (Figs 5–11) is assessed only at the level of design.

---

## 1. Verdict in one paragraph

The wet-lab half is competent and, on its own terms, convincing: a dose-ranged, randomised,
ethics-approved CIA rat study with a positive control, backed by in vitro work and a pharmacological
rescue. The bioinformatics half that precedes it is much weaker than its presentation suggests. Three
findings in particular — the single-cell comparison, the Mendelian randomisation, and the link
between the two — rest on data that cannot support the confidence with which they are stated. The
paper's own limitations section discusses only docking, formula complexity and the absence of genetic
knockdown; **not one of the bioinformatics weaknesses is acknowledged anywhere in the text.** The
conclusions are probably right, but they are right because of the wet-lab work, not because of the
omics that nominated the target.

---

## 2. What the paper does well

- **A genuinely complete chain.** Bulk → WGCNA → enrichment → single cell → MR → chemistry → docking
  → in vitro → in vivo. Very few papers carry a target the whole way. The framework is worth copying
  even where the execution is not.
- **The in vivo design is sound.** 48 SPF male Sprague-Dawley rats, six groups of eight, randomised,
  ethics approval (Guizhou Medical University 2200095), three SX doses plus a tripterygium glycosides
  positive control. Dose-response is what makes the in vivo result credible.
- **A rescue experiment.** Fingolimod hydrochloride is used to test PAK1 involvement rather than only
  correlating with it.
- **MR diagnostics are run, not skipped.** Five estimators, Cochran's Q, MR-Egger intercept, Steiger
  filtering, F-statistic > 10, LD clumping at r² < 0.001 / 10 Mb. The machinery is right even though
  the input is thin (§3.2).
- **Molecular dynamics in triplicate,** 100 ns each. Many docking papers report a single run.
- **Batch correction across three bulk datasets** (sva), with PCA shown before and after.

---

## 3. Serious problems

### 3.1 The single-cell analysis is n = 1 per group

The scRNA-seq data (GSE221704) is **three mice: one CIA, one vehicle control, one CIA treated with
soyasaponin Bb.** One animal per arm.

Figure 3D–E reports Rac1 at P < 0.001 and Pak1 at P < 0.0001 between groups, by Wilcoxon test on
cells with BH adjustment. With one mouse per group, cells are not independent replicates of the
biological contrast; every cell from the CIA animal shares that animal's genotype, cage, batch,
dissection and capture lane. The test measures how many cells were sequenced, not how different the
groups are. Those four-star P-values would not survive any per-animal test, and no per-animal test is
possible at n = 1.

This is the most consequential flaw in the paper, because Fig. 3 is the step that promotes RAC1/PAK1
from "enriched pathway" to "the mechanism".

### 3.2 The MR rests on 216 cases and 3 SNPs

Table 1 gives the RA outcome GWAS (`ebi-a-GCST90018873`) as **216 cases and 409,001 controls**. The
PAK1 exposure is a single gene's cis-eQTL (`eqtl-a-ENSG00000149269`), instrumented by **3 SNPs**.

The headline result is IVW **P = 0.039, OR 1.578 (95% CI 1.024–2.434)** — a confidence interval whose
lower bound is 1.024, i.e. barely clear of the null. With 216 cases the outcome effect estimates are
very imprecise, and a three-instrument IVW is one outlying SNP away from crossing 1.0. Figure 4A
shows this directly: MR-Egger gives OR 0.634 (95% CI 0.082–4.922), pointing the *opposite* way with
an interval spanning two orders of magnitude, and none of the other four estimators reaches P < 0.05.
The paper's rule — "if the P-value of the results obtained by the IVW method was less than 0.05, it
was regarded as having a causal relationship" — means the entire causal claim is one estimator at
P = 0.039 while four others disagree.

Reporting "no heterogeneity (Q, P > 0.05)" and "no pleiotropy (Egger intercept, P > 0.05)" as
reassurance is misleading at this size: with 3 instruments both tests have almost no power, so
passing them is uninformative rather than supportive.

Either 216 is the real case count, in which case the MR is not interpretable, or it is a
transcription error, in which case a core table of the paper is wrong. Neither is good.

### 3.3 Species and tissue both change under the reader's feet

- **Bulk:** human RA **synovial tissue** (GSE1919, GSE55235, GSE77298; 31 RA, 22 HC).
- **Single cell:** mouse **femoral bone marrow** (GSE221704).
- **MR:** human whole-blood cis-eQTL → human RA risk.
- **In vivo:** rat CIA synovium and paw.

The conclusion states RAC1/PAK1 activation "is involved in neutrophil migration **within RA
synovium**". The single-cell evidence for that sentence comes from bone marrow, which is where
neutrophils are *made*, not where they migrate *to*. Neutrophil Rac1/Pak1 expression in marrow is not
evidence about synovial migration, and the paper never addresses the gap — it is not raised as a
limitation, and the figure legend says "bone marrow tissue" while the conclusion says "synovium".

### 3.4 The treated arm is a different drug

The single-cell "Treat" group is a CIA mouse given **soyasaponin Bb** — not Sidaxue. The results
text reads "decreased significantly after treatment (P < 0.0001), suggesting that pharmacological
intervention targeting RAC1/PAK1 signaling may ameliorate RA". That sentence is defensible only as a
statement about soyasaponin Bb. Placed inside a Sidaxue paper, immediately before the Sidaxue
experiments, it reads as support for Sidaxue. The compound's name appears once, in Methods, and never
again.

### 3.5 Selection steps are arbitrary and unvalidated

- **GeneCards relevance score > 8** selects the 751-gene phenotype set. GeneCards relevance is a
  text-mining aggregate, not a curated ontology; the threshold is not justified and no sensitivity
  analysis is shown.
- **"The module showing the strongest positive correlation with the RA disease group was selected."**
  One module, chosen on the pooled correlation, with no requirement that it replicate in each of the
  three constituent datasets. Module–trait correlations are not corrected for the number of modules
  tested.
- **Regulation of actin cytoskeleton was not the top KEGG term.** It sits around 17th of ~40 in
  Fig. 2B. It was chosen because it matched the phenotype. That is a legitimate move, but it is a
  judgement call presented as a result — the text says KEGG "indicated strong enrichment in
  regulation of actin cytoskeleton" without noting the sixteen terms above it.
- **RAC1/PAK1 within that pathway** were picked by eye from the pathview diagram (Fig. 2D), with no
  stated rule.

Four consecutive judgement calls, each reasonable, none pre-registered, none tested for sensitivity.
Any one of them changes the target.

### 3.6 Circularity between Figs 1 and 2

The 159-gene set is defined partly by a neutrophil-migration gene list, and is then tested for
enrichment in GO terms about leukocyte migration and chemotaxis (Fig. 2A). Terms related to migration
*must* come out on top; the gene set was built from migration genes. This is reported as a finding
rather than as a property of the construction.

### 3.7 The limitations section omits the bioinformatics entirely

Three limitations are given: docking is predictive not conclusive; the formula is multi-component;
no genetic knockdown. All three concern the wet-lab half. Nothing about n = 1 single-cell, 216 GWAS
cases, the synovium/bone-marrow mismatch, the soyasaponin substitution, or the unvalidated selection
thresholds. A reader who trusts the limitations section will substantially over-read Figs 1–4.

---

## 4. What this means for this project

Where the framework was copied, these are the points of inheritance.

| Their problem | This project |
|---|---|
| scRNA n = 1 per group, cell-level P-values | **Avoided.** 3 vs 4 mice; the bracket on Fig 3D/E is the per-mouse t-test, with the cell-level Wilcoxon reported in the legend as the anticonservative number it is. |
| MR on 216 cases, 3 SNPs, one estimator | **Avoided.** Outcome is GCST006867 (~74k T2D cases); 9 on-map exposures, all estimators shown, FDR over the map. |
| Single best module, no replication | **Avoided.** Decision M11b requires the T2D association to replicate in both cohorts separately, which cut key modules from 7 to 4 and module genes from 4,022 to 418. |
| Pathway chosen by eye | **Partly avoided.** Script 16 scores six criteria. But `docs/pathway_selection.md` still reports a pre-M11b ranking in which AMPK is rank 1; on current data it is rank 3 (rank 1 on ORA and GSEA). **Still outstanding.** |
| Axis genes chosen by eye from the pathview | **Avoided.** Decision R23 restricts to KGML edges with a candidate-set node; script 25 checks reachability. |
| Circular enrichment (phenotype set → enrichment) | **Shared.** Fig 2's candidate set is defined partly by a hyperglycaemia gene set and then tested against overlapping annotation. Stated in `figures/Fig2/Fig2_legend.md`, not designed away. |
| Tissue/species mismatch | **Shared, and stated.** Human liver bulk → mouse liver/kidney scRNA, plus an STZ vs obese-T2D model mismatch. `docs/pathway_selection.md` §3.3 gives three structural reasons; the figure legends repeat them. |
| Limitations omit the bioinformatics | **Avoided.** Every figure legend carries an explicit deviations section. |

**The uncomfortable comparison.** The reference's framework produced a positive result at every
stage. This project's produces negatives at several: no pair clears the single-cell gates on the
selected map, and nothing generalises in the classifier. That difference is mostly *not* because the
biology is weaker — it is because the reference's gates were evaluated on data too small to fail.
Running the same framework on adequately powered data is what turns four-star P-values into
empirical P = 0.60.

---

## 5. Should this framework still be the template?

Yes, with the structure kept and the thresholds replaced — which is what has been done. The chain
itself (candidate set → pathway → cell type → causality → validation) is a reasonable way to organise
a target-nomination study, and matching it panel for panel makes the work legible to reviewers who
know the source.

What should **not** be carried over is the inference style: single-estimator causal claims, cell-level
statistics on one animal, thresholds chosen post hoc and reported as results, and a limitations
section that protects the computational work from scrutiny.

One practical consequence for writing: this project cannot claim a causal druggable target the way
the reference does. `figures/Fig5/Fig5_legend.md` says so directly. That is a more honest position
than the reference's, and it should be presented as a strength of the analysis rather than apologised
for.

---

## 6. Outstanding item raised by this review

`docs/pathway_selection.md` §2.1 reports the pre-M11b pathway ranking (AMPK rank 1, composite 3.33)
while `results/pathway_selection/selection_ranked.csv` now puts AMPK **3rd** (composite 5.00) behind
Insulin resistance and Insulin signalling, both of which entered after the candidate set shrank from
19 genes to 8. AMPK remains **rank 1 on ORA** (adj. P 5.8e-6 against 6.8e-4 for the next map) and
rank 1 on GSEA, and it is the only map hitting all four KEGG-mapped candidates. Note also that
hsa04931 is a KEGG *disease* map rather than a signalling map, which the doc's phrase "rank 1 among
the signalling maps" implies was meant to be excluded — but script 16 applies no such filter.

The selection is defensible; the document justifying it is out of date. It should be regenerated
against current data before submission.
