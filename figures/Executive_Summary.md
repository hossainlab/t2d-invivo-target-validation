# Executive Summary & Central Conclusion
**Project:** Integrated Multi-Omics Target Validation for In Vivo Type 2 Diabetes (T2D) Study  
**Framework:** Adapted from *Xu et al., Phytomedicine 154 (2026) 158050*  
**Date:** September 2026  

---

## 1. Executive Summary

This study establishes a rigorous, end-to-end multi-omics target nomination and validation pipeline for Type 2 Diabetes (T2D). By integrating **human liver bulk transcriptomics (WGCNA)**, **KEGG pathway topology (KGML)**, **mouse single-cell RNA sequencing (liver & kidney atlases)**, **machine learning consensus modeling**, and **Mendelian Randomization (MR)**, we identify the primary mechanistic regulatory axis driving diabetic tissue stress and nominate prioritized targets for in vivo pharmacological validation.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                                 THE 5-STAGE ANALYTIC FUNNEL                            │
│                                                                                        │
│  [Stage 1: Human Bulk]    1,648 DEGs ∩ WGCNA Key Modules ∩ Hyperglycemia Set           │
│                                           │ (8 Core Candidates)                        │
│  [Stage 2: Pathway/KGML]  Top KEGG Map: AMPK Signaling (hsa04152; ORA q = 5.8e-6)      │
│                                           │ (Candidate-Anchored Axis)                  │
│  [Stage 3: scRNA Atlas]   Cross-Species Cell Confirmation (Liver Cholang. & Kidney EC) │
│                                           │ (Spatial Co-Localization Gate)             │
│  [Stage 4: ML Consensus]  Feature Ranking (LASSO / Random Forest / SVM Consensus)      │
│                                           │ (Multi-Variable Predictability)            │
│  [Stage 5: Human Genetics] cis-eQTL Mendelian Randomization (Outcome: 74k T2D Cases)   │
│                                           │                                            │
│  [DELIVERABLE]            Tiered Target Portfolio for In Vivo Validation (GV1 vs Met)  │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 2. Central Biological Conclusion

### Suppression of the AMPK Signaling Axis Output
The central finding across human and rodent diabetic datasets is the **functional suppression of AMPK pathway activity** under chronic hyperglycemic/insulin-resistant stress:

> **Suppressed AMPK Output:**
> * **PPARGC1A (PGC-1α) ↓** → Loss of mitochondrial biogenesis and hepatic metabolic control.
> * **CCND1 (Cyclin D1) ↑** → De-repressed cellular stress and aberrant proliferation signal.

### Key Mechanistic Insights:
1. **Post-Translational Control Signature:**  
   The human KEGG network displays discordant edge dynamics: the AMPK kinase subunit (*PRKAA2*) is transcriptionally maintained while its activated downstream targets fall. AMPK activity is governed by **Thr172 phosphorylation** and cofactor ratios, not mRNA abundance.
2. **Epithelial & Endothelial Compartmentalization:**  
   * In **liver**, metabolic suppression (*Ppargc1a* downregulation) localizes specifically to epithelial **cholangiocytes** (per-mouse *P* = 0.011, detection dropping from 31.7% to 11.4%).
   * In **kidney**, stress/cell-cycle de-repression (*Ccnd1* upregulation) localizes to **microvascular endothelial cells** (pseudobulk *P* = 1.2 × 10⁻⁶).
3. **Cross-Species Concordance:**  
   The directional shifts of both primary anchors (*PPARGC1A* down, *CCND1* up) replicate across human clinical cohorts and mouse diabetic models.

---

## 3. Core Target Tiers & Evidence Hierarchy

Targets are classified mechanically into prioritized tiers based on evidence convergence across all five analytical stages:

```
                               ┌──────────────────────────────────────────────┐
                               │       UPSTREAM MECHANISM (Tier 3)            │
                               │  p-AMPKα (PRKAA1/2) · SIRT1 · STK11 (LKB1)   │
                               │  (Phosphorylation / Post-translational node) │
                               └──────────────────────┬───────────────────────┘
                                                      │
                       ┌──────────────────────────────┴──────────────────────────────┐
                       ▼                                                             ▼
         ┌───────────────────────────┐                                 ┌───────────────────────────┐
         │     PRIMARY TARGET 1      │                                 │     PRIMARY TARGET 2      │
         │    PPARGC1A (PGC-1α)      │                                 │      CCND1 (Cyclin D1)    │
         │  Metabolic Co-activator   │                                 │    Pathway Output Node    │
         │   Human T2D: log2FC -0.61 │                                 │   Human T2D: log2FC +0.80 │
         │   Mouse Liver Cholang: ↓  │                                 │    Mouse Endothelium: ↑   │
         │     (ML Rank 3 / 19)      │                                 │    (GSE23343 Replicated)  │
         └───────────────────────────┘                                 └───────────────────────────┘
                       │                                                             │
                       ▼                                                             ▼
         ┌───────────────────────────┐                                 ┌───────────────────────────┐
         │     READOUT: GLYCOLYSIS   │                                 │     READOUT: ARREST/p53   │
         │           PFKFB3          │                                 │     CDKN1A (p21) · BAX    │
         └───────────────────────────┘                                 └───────────────────────────┘
```

### Target Classification Table

| Target | Tier | Role & Mechanism | Key Supporting Evidence | Planned In Vivo Assay |
| :--- | :--- | :--- | :--- | :--- |
| **`PPARGC1A`** *(PGC-1α)* | **Tier 1 (Primary)** | Master transcriptional coactivator of mitochondrial respiration & gluconeogenesis. | • Human liver: log₂FC = −0.61, *P* = 6.0 × 10⁻⁴<br>• Mouse cholangiocytes: *P* = 0.011, cell BH *P* = 2.0 × 10⁻⁸<br>• ML Consensus: Rank 3 / 19 | • RT-qPCR & Western blot (liver)<br>• IF co-staining: PGC-1α × CK19 |
| **`CCND1`** *(Cyclin D1)* | **Tier 1 (Primary)** | De-repressed downstream output of AMPK ⊣ HuR/ELAVL1 axis. | • Human liver: log₂FC = +0.80, *P* = 7.5 × 10⁻⁴<br>• Replicated in human cohort GSE23343 (*P* = 0.016, AUC = 0.81)<br>• Mouse kidney endothelium: *P* = 1.2 × 10⁻⁶ | • RT-qPCR & Western blot (kidney & liver)<br>• IF co-staining: Cyclin D1 × CD31 |
| **`CDKN1A`** *(p21)* | **Tier 2 (Readout)** | Cell-cycle arrest / senescence effector downstream of metabolic stress. | • Human liver DEG: log₂FC = +0.96, FDR = 0.016<br>• Elevated in mouse cholangiocytes (FDR = 0.009) & kidney endothelium (FDR = 0.008) | • RT-qPCR & Western blot<br>• IF co-staining with CD31 in kidney |
| **`PFKFB3`** | **Tier 2 (Readout)** | Inducible glycolytic regulator connected to AMPK. | • Human liver: FDR = 0.042<br>• Mouse liver LSEC (*P* = 0.049) & kidney LoH (BH *P* = 0.03) | • RT-qPCR & Western blot (liver & kidney) |
| **`PRKAA1/2`** *(AMPKα)* | **Tier 3 (Master Node)** | Central metabolic sensor and upstream serine/threonine kinase. | • Upstream driver on KEGG hsa04152; acts via Thr172 phosphorylation. | • Western blot for **p-AMPKα (Thr172) / total AMPKα**<br>• Metformin positive control |

---

## 4. Multi-Modal Evidence Matrix

| Candidate Gene | Human Bulk DE | WGCNA Module | KEGG AMPK Map | Mouse scRNA Validation | Human Replication | ML Consensus | Causal MR (GWAS) | Target Tier |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **PPARGC1A** | Down (*P* = 6.0 × 10⁻⁴) | Key Module | Direct Step (+1) | Down (Cholang. *P* = 0.011) | Concordant | Rank 3 | Null (*P* = 0.16) | **Tier 1** |
| **CCND1** | Up (*P* = 7.5 × 10⁻⁴) | Key Module | 2 Steps via HuR (−1) | Up (Kidney EC *P* = 1.2 × 10⁻⁶) | Replicated (*P* = 0.016) | Rank 14 | Null (*P* = 0.55) | **Tier 1** |
| **CDKN1A** | Up (FDR = 0.016) | Key Module | p53 / Arrest Map | Up (Cholang. & Kidney EC) | Concordant | — | Null | **Tier 2** |
| **PFKFB3** | Down (FDR = 0.042) | Turquoise | Direct Step (+1) | Down (LSEC *P* = 0.049) | Concordant | — | Null | **Tier 2** |
| **IRS2** | Down (*P* = 0.011) | Key Module | Input Node | Minor Cells only | Concordant | Rank 6 | Null | **Tier 2** |
| **IGF1** | Down (*P* = 0.021) | Key Module | Input Node | Sparse (Hepatocyte gap) | Concordant | Rank 8 | Null | **Tier 2b** |
| **PRKAA2** | Up (*P* = 0.017) | Network Node | Master Kinase | Post-translational | Concordant | — | Null | **Tier 3** |
| **SREBF1** | Down (*P* = 0.041) | Network Node | Lipogenic Output | Inverted (STZ model) | — | — | Causal (Blood, *P* = 9.0 × 10⁻⁹) | **Exploratory** |
| **VSNL1 / LZTFL1** | DE (FDR < 0.05) | Key Module | Off-Map | Not DE in mouse | Weak | Rank 1 & 2 | Null | **Non-Target** |

---

## 5. In Vivo Experimental Validation & Pharmacological Rescue Plan

The findings directly define the experimental design and endpoints for the in vivo animal validation study:

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                              IN VIVO VALIDATION DESIGN (GV1 vs MET)                     │
│                                                                                        │
│  [Groups]       1. Normal Control (Saline)                                             │
│                 2. Ob-T2D Disease Model (High-Fat Diet / STZ)                          │
│                 3. Ob-T2D + GV1 Low Dose                                               │
│                 4. Ob-T2D + GV1 Medium Dose                                            │
│                 5. Ob-T2D + GV1 High Dose                                              │
│                 6. Ob-T2D + Metformin (Positive Control, 200 mg/kg)                    │
│                                                                                        │
│  [Primary]      1. Glycemic & Metabolic Efficacy: OGTT, Fasting Glucose, HbA1c, HOMA-IR│
│  [Molecular]    2. Hepatic AMPK Restoration: Western blot p-AMPKα (Thr172) / total     │
│                 3. Mitochondrial Rescue: PGC-1α expression (RT-qPCR & Western blot)    │
│                 4. Proliferation / Stress Normalization: Cyclin D1 & p21 (qPCR & WB)   │
│  [Histology]    5. Immunofluorescence: PGC-1α × CK19 (liver); Cyclin D1 × CD31 (kidney)│
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### Hypothesized Therapeutic Mechanism of Action:
Administration of the active formulation (**GV1**) is hypothesized to mirror the therapeutic action of **Metformin** by:
1. Re-activating hepatic **AMPK phosphorylation (Thr172)**.
2. Restoring downstream **PGC-1α** expression and mitochondrial function.
3. Normalizing aberrant **Cyclin D1** and **p21** levels, resolving chronic tissue stress in liver and kidney.
