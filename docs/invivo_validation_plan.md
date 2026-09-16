# In vivo validation plan — GV1 in the Ob-T2D mouse

**Status:** pre-specified, written before any animal data are seen (2026-09-16).
**Rationale and evidence:** `docs/target_selection.md`, `docs/pathway_selection.md`, `docs/results_scRNA.md`.
**Model:** HFD/STZ obese type 2 diabetic mice, Kyung Hee University lab. **Groups:** Control, Model (Ob-T2D), GV1 (*Bacillus velezensis*), Metformin. **Tissues:** liver and kidney.

> This document fixes the endpoints, the sample size, the sample handling and the analysis **before** the experiment. Anything decided after seeing the data is exploratory and must be labelled as such in the manuscript.

---

## 1. What this experiment has to decide

The computational work produced two internally coherent readings of the same data, and it cannot separate them with 16 T2D and 11 lean human livers (`docs/target_selection.md` §2.6):

| Reading | Genes | Holds at | What it predicts in vivo |
|---|---|---|---|
| **A — mechanistic** | AMPK → PGC-1α suppressed; cyclin D1 and p21 de-repressed | nominal P to FDR 0.10 | **p-AMPKα Thr172/AMPKα falls in Model and is restored by GV1**; PGC-1α down/restored; cyclin D1 and p21 up/reversed |
| **B — statistically strict** | VSNL1, LZTFL1, GAS6, SREBF2 | FDR < 0.05 | These four change in Model and are reversed by GV1, **with no AMPK involvement** |

Both are testable here, and testing only A would make the experiment unfalsifiable. Reading B costs four extra qPCR wells.

---

## 2. Endpoint hierarchy

Tested in order. A failure at one level does not invalidate the levels below, but the levels below are then reported as exploratory.

| Level | Endpoint | Contrast(s) | Test | α |
|---|---|---|---|---|
| **Primary** | **p-AMPKα Thr172 / total AMPKα, liver (Western blot)** | Model vs Control | two-sided t-test (Welch) | 0.05 |
| **Key secondary (gatekept by the primary)** | same ratio | GV1 vs Model; Metformin vs Model | Dunnett vs Model | 0.05 family-wise |
| Secondary 1 | PGC-1α protein (liver), cyclin D1 protein (kidney), p21 protein (kidney) | Model vs Control; GV1 vs Model | Dunnett | BH across the 3 targets |
| Secondary 2 | qPCR panel, reading A: *Ppargc1a, Sirt1, Irs2, Igf1, Srebf1, Fasn, Fbp1, Pfkfb3, Cdkn1a, Ccnd1, Phlda3, Zmat3, Bax* | Model vs Control; GV1 vs Model | Dunnett | BH within panel |
| Secondary 3 | qPCR panel, **reading B**: *Vsnl1, Lztfl1, Gas6, Srebf2* | same | Dunnett | BH within panel |
| Secondary 4 | cell-type localisation: PGC-1α × CK19 (liver cholangiocytes); cyclin D1 × CD31 and p21 × CD31 (kidney endothelium) | Model vs Control; GV1 vs Model | blinded scoring, Dunnett | BH across the 3 stains |
| Exploratory | SIRT1, SREBP-1c (mature), FASN, p-p53 Ser18/p53, MDM2, histology, correlation of p-AMPK with glycemia | — | — | reported unadjusted, labelled exploratory |

**Assay-validity control, declared in advance:** metformin must raise the liver p-AMPKα/AMPKα ratio relative to Model. **If it does not, the blot has failed and the primary endpoint is uninterpretable — not negative.** Re-run before drawing any conclusion about Model vs Control.

---

## 3. Sample size

Two-sided t-test, α = 0.05, 80 % power, effect sizes in units of pooled SD (Cohen's d).

| n per group | Minimum detectable d | d with Holm across 3 contrasts | Power at d = 1.0 | d = 1.5 | d = 2.0 |
|---|---|---|---|---|---|
| 6 | 1.80 | 2.19 | 0.35 | 0.65 | 0.88 |
| 8 | 1.51 | 1.81 | 0.46 | 0.80 | 0.96 |
| **10** | **1.32** | 1.57 | 0.56 | 0.89 | 0.99 |
| 12 | 1.20 | 1.41 | 0.65 | 0.94 | 1.00 |
| 15 | 1.06 | 1.25 | 0.75 | 0.98 | 1.00 |
| 20 | 0.91 | 1.06 | 0.87 | 1.00 | 1.00 |

n required for 80 % power: d = 0.8 → 26; d = 1.0 → 17; d = 1.2 → 12; d = 1.5 → 9; d = 2.0 → 6.

**Recommendation: n = 10–12 per group for the WB cohort.** Published p-AMPK differences between diabetic and control rodent liver are typically large (d ≈ 1.5–2), so n = 10 gives ~0.9 power for the primary endpoint while still detecting a d ≈ 1.3 GV1 effect. At n = 6 the study can only detect d ≥ 1.8, and a null result would be uninformative — this must not be the design.

Allow ~15 % attrition in the STZ model; randomise **12 per group at allocation** to end with ≥ 10 analysable.

---

## 4. Randomisation, blinding, and confound control

1. **Allocation:** randomise cages, not individual mice, to treatment; record the randomisation list before dosing.
2. **Harvest order:** randomise across groups (e.g. Control, GV1, Model, Metformin, repeat), never group-by-group. Time-to-freeze is the main confounder for a phospho-endpoint and must not track with group.
3. **Record time from cervical dislocation to freeze for every animal** and include it as a covariate if it varies by more than ±30 s.
4. **Blinding:** the person running densitometry, qPCR analysis and IHC scoring is blinded to group; samples are relabelled by a third party.
5. **Confirm the model before analysing anything else:** fasting glucose, body weight, HbA1c per animal. A Model group that is not hyperglycemic invalidates the whole comparison.

---

## 5. Sample handling (critical for the primary endpoint)

p-AMPK Thr172 is lost within minutes of ischemia. The primary endpoint depends entirely on this section.

| Step | Requirement |
|---|---|
| Harvest | Excise liver/kidney and **clamp-freeze in liquid N₂ within 30 seconds**; record the time |
| Storage | −80 °C; never thaw/refreeze before lysis |
| Lysis | Ice-cold RIPA with **protease *and* phosphatase inhibitor cocktails** (Na₃VO₄, NaF, β-glycerophosphate or a commercial equivalent), added fresh on the day |
| Loading | Equal protein (BCA), same gel for all four groups of a given target; no cross-gel comparison for the primary ratio |
| Normalisation | p-AMPKα to **total AMPKα on the same membrane** (strip/re-probe or parallel gel), not to a housekeeper |
| Portioning at harvest | Three aliquots per organ: (a) snap-frozen for WB, (b) RNAlater or snap-frozen for RNA, (c) 10 % NBF → paraffin for IHC/IF |

---

## 6. Assay panel

### 6.1 Western blot

| Target | Purpose | Notes |
|---|---|---|
| **p-AMPKα Thr172; total AMPKα** | **primary endpoint (ratio)** | rabbit mAb, mouse-reactive; the pair must come from the same supplier/clone family |
| PGC-1α | Tier 1 target | liver; validate the band (~90–110 kDa, multiple isoforms reported) |
| Cyclin D1 | Tier 1 target | kidney primarily |
| p21 (CDKN1A) | Tier 2 readout | kidney and liver |
| SIRT1 | Tier 3 mechanism | exploratory |
| SREBP-1c (precursor + mature), FASN | Tier 3 / lipogenic arm | exploratory; expect the STZ-type discordance discussed in `pathway_selection.md` §3.3 |
| p-p53 Ser18 / total p53, MDM2 | arrest arm control | exploratory |
| β-actin or GAPDH | loading | **not** used to normalise the phospho/total ratio |

Confirm catalogue numbers and mouse reactivity before ordering; run a positive-control lysate for p-AMPK (e.g. AICAR- or metformin-treated hepatocytes) on the first blot.

### 6.2 qPCR

Reference genes: **two**, validated for stability in diabetic mouse liver/kidney (e.g. *Ppia*, *Rplp0*, *Tbp*); report geometric mean normalisation and the stability statistic. ΔΔCt with efficiency correction; primers must span an exon–exon junction and show 90–110 % efficiency and a single melt peak.

| Panel | Genes |
|---|---|
| Reading A — AMPK arm | *Prkaa1, Prkaa2, Ppargc1a, Sirt1, Irs2, Igf1, Fbp1, Pfkfb3, Srebf1, Fasn* |
| Reading A — arrest arm | *Cdkn1a, Ccnd1, Phlda3, Zmat3, Bax* |
| **Reading B — strict alternative** | ***Vsnl1, Lztfl1, Gas6, Srebf2*** |
| Reference | *Ppia, Rplp0* (or lab-validated equivalents) |

Primer sequences are **not** specified here: design them against the current RefSeq for each mouse gene and validate in-house rather than copying from literature.

### 6.3 IHC / IF (the panel that distinguishes this study from a tissue-lysate paper)

| Tissue | Co-stain | Why |
|---|---|---|
| Liver | **PGC-1α × CK19** | the single-cell data localise the Ppargc1a change to cholangiocytes |
| Kidney | **cyclin D1 × CD31** | Ccnd1 rises specifically in endothelium (pseudobulk P 1.2e-6) |
| Kidney | **p21 × CD31** | the arrest readout in the same compartment |

Score ≥ 5 fields per animal, blinded; report the fraction of marker-positive cells that are also target-positive, not whole-section intensity.

---

## 7. Pre-specified analysis

1. **Model check first.** Confirm Model vs Control differ in fasting glucose. If not, stop and report.
2. **Primary:** Welch t-test, Model vs Control, liver p-AMPKα/AMPKα. Report effect size with 95 % CI, not only P.
3. **Assay-validity gate:** Metformin vs Model on the same ratio must be positive.
4. **Key secondary:** Dunnett against Model for GV1 and Metformin.
5. **Secondary panels:** Dunnett, BH within each panel; panels are not pooled.
6. **Reading A vs B:** report both panels side by side. Pre-committed interpretation in §8.
7. **Correlation:** p-AMPKα/AMPKα against fasting glucose across all animals (Spearman), as a dose-response check.
8. **Transparency:** every animal, every blot, every Ct in the supplement; no exclusions without a documented, pre-specified criterion (e.g. failed model induction).

---

## 8. Interpretation, committed in advance

| Outcome | Conclusion |
|---|---|
| p-AMPK ratio down in Model, restored by GV1, metformin control positive | Reading **A** supported; GV1 acts through AMPK. Strongest result available from this design |
| p-AMPK ratio down in Model, **not** restored by GV1, but PGC-1α/cyclin D1/p21 reverse | The axis is engaged but GV1 acts downstream or in parallel; report as such, do not claim AMPK activation |
| p-AMPK ratio unchanged, metformin control positive | Reading A **not** supported in this model. The human finding remains, the mouse mechanism does not transfer — a publishable negative that the threshold analysis already anticipated |
| Metformin control negative | Assay failure. No conclusion about any group |
| Reading B genes move, reading A genes do not | The FDR-strict interpretation wins; the manuscript's mechanism section is rewritten around VSNL1/LZTFL1/GAS6/SREBF2 |
| Nothing moves | Model, handling or power problem — check the glucose data and the time-to-freeze log before interpreting |

---

## 9. Before the first animal is dosed — checklist

- [ ] n fixed at 12/group allocated (≥ 10 analysable); randomisation list generated and filed
- [ ] harvest order randomised across groups; time-to-freeze log sheet prepared
- [ ] phosphatase inhibitors in stock; p-AMPK positive-control lysate available
- [ ] antibodies confirmed mouse-reactive; p-AMPK and total AMPK from the same clone family
- [ ] qPCR primers designed and efficiency-tested, **including the four reading-B genes**
- [ ] blinding scheme agreed, relabelling person identified
- [ ] this document circulated to the KHU team and version-frozen

---

## 10. What this experiment cannot do

- It cannot establish causality in humans; the MR analysis was null for every candidate (`docs/results_MR.md`).
- It cannot validate a diagnostic signature; the 19-gene panel did not transport to an independent human cohort (`docs/results_ML.md`).
- It tests whether the mechanism proposed from human liver transcriptomes operates in one mouse model of Ob-T2D, and whether GV1 moves it. That is the claim the paper should make.
