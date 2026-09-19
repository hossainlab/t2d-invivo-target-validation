#set page(
  paper: "a4",
  margin: (top: 2cm, bottom: 2cm, left: 2cm, right: 2cm),
  header: locate(loc => {
    if loc.page() > 1 {
      text(size: 8.5pt, fill: rgb("64748b"), font: "Segoe UI")[
        *Executive Summary & Central Conclusion* | T2D In Vivo Target Validation
        #h(1fr)
        September 2026
      ]
    }
  }),
  footer: locate(loc => {
    text(size: 8.5pt, fill: rgb("64748b"), font: "Segoe UI")[
      Hossain Lab · In Vivo Target Validation
      #h(1fr)
      Page #loc.page() of #counter(page).final(loc).at(0)
    ]
  })
)

#set text(
  font: ("Segoe UI", "Arial", "Roboto", "Helvetica"),
  size: 10pt,
  fill: rgb("1e293b"),
  lang: "en"
)

#set par(justify: true, leading: 0.65em)

// --- Custom Theme Colors ---
#let primary = rgb("1e40af")     // Deep Blue
#let secondary = rgb("0369a1")   // Ocean Blue
#let dark = rgb("0f172a")        // Slate Dark
#let light-bg = rgb("f8fafc")    // Off white
#let border-col = rgb("cbd5e1")  // Slate border
#let accent-bg = rgb("eff6ff")   // Soft blue tint
#let tier1-bg = rgb("dcfce7")    // Soft green
#let tier1-text = rgb("15803d")
#let tier2-bg = rgb("fef3c7")    // Soft amber
#let tier2-text = rgb("b45309")
#let tier3-bg = rgb("f3e8ff")    // Soft purple
#let tier3-text = rgb("7e22ce")

// --- Badges ---
#let badge(body, bg: accent-bg, text-color: primary) = {
  box(
    fill: bg,
    radius: 4pt,
    inset: (x: 5pt, y: 2.5pt),
    outset: 0pt,
    baseline: 0%,
    text(size: 8pt, weight: "bold", fill: text-color, body)
  )
}

// --- Header Block ---
#block(
  width: 100%,
  stroke: (bottom: 2pt + primary),
  inset: (bottom: 12pt),
  margin: (bottom: 14pt),
  [
    #grid(
      columns: (1fr, auto),
      gutter: 10pt,
      [
        #text(size: 20pt, weight: "bold", fill: dark)[Executive Summary & Central Conclusion]
        
        #v(3pt)
        #text(size: 11pt, weight: "medium", fill: secondary)[
          Integrated Multi-Omics Target Validation for In Vivo Type 2 Diabetes (T2D) Study
        ]
      ],
      [
        #align(right)[
          #badge("CONFIDENTIAL & PRE-CLINICAL", bg: rgb("fee2e2"), text-color: rgb("b91c1c"))
          #v(4pt)
          #text(size: 8.5pt, fill: rgb("64748b"))[
            *Framework:* Xu et al. (2026)\
            *Date:* September 2026
          ]
        ]
      ]
    )
  ]
)

// --- Section 1: Executive Summary ---
== 1. Strategic Overview & Analytic Funnel

This study establishes a rigorous, end-to-end multi-omics target nomination and validation pipeline for Type 2 Diabetes (T2D). By integrating *human liver bulk transcriptomics (WGCNA)*, *KEGG pathway topology (KGML)*, *mouse single-cell RNA sequencing (liver & kidney atlases)*, *machine learning consensus modeling*, and *Mendelian Randomization (MR)*, we identify the primary mechanistic regulatory axis driving diabetic tissue stress and nominate prioritized targets for in vivo pharmacological validation.

#v(4pt)

#align(center)[
  #block(
    fill: light-bg,
    stroke: 1pt + border-col,
    radius: 6pt,
    inset: 10pt,
    width: 100%,
    [
      #grid(
        columns: (1fr, auto, 1.1fr, auto, 1.1fr, auto, 1fr, auto, 1fr),
        align: center + horizon,
        gutter: 4pt,
        [
          #text(weight: "bold", size: 8.5pt, fill: primary)[Stage 1: Bulk]\
          #text(size: 7.5pt, fill: rgb("475569"))[1,648 DEGs ∩\ WGCNA ∩ Set\ *(8 Candidates)*]
        ],
        [#text(fill: rgb("94a3b8"), size: 11pt)[$arrow.r$]],
        [
          #text(weight: "bold", size: 8.5pt, fill: primary)[Stage 2: Topology]\
          #text(size: 7.5pt, fill: rgb("475569"))[AMPK Signaling\ (hsa04152)\ *(Anchored Axis)*]
        ],
        [#text(fill: rgb("94a3b8"), size: 11pt)[$arrow.r$]],
        [
          #text(weight: "bold", size: 8.5pt, fill: primary)[Stage 3: scRNA]\
          #text(size: 7.5pt, fill: rgb("475569"))[Liver Cholang. &\ Kidney Endoth.\ *(Co-Localization)*]
        ],
        [#text(fill: rgb("94a3b8"), size: 11pt)[$arrow.r$]],
        [
          #text(weight: "bold", size: 8.5pt, fill: primary)[Stage 4: ML]\
          #text(size: 7.5pt, fill: rgb("475569"))[LASSO / RF / SVM\ *(Multi-Variable\ Predictability)*]
        ],
        [#text(fill: rgb("94a3b8"), size: 11pt)[$arrow.r$]],
        [
          #text(weight: "bold", size: 8.5pt, fill: primary)[Stage 5: MR]\
          #text(size: 7.5pt, fill: rgb("475569"))[cis-eQTL MR\ 74k T2D Cases\ *(Genetics)*]
        ]
      )
      #v(6pt)
      #line(length: 100%, stroke: 0.5pt + rgb("e2e8f0"))
      #v(3pt)
      #text(weight: "bold", size: 8.5pt, fill: dark)[
        #text(fill: rgb("15803d"))[Deliverable:] Tiered Target Portfolio for In Vivo Animal Validation (GV1 Formulation vs. Metformin)
      ]
    ]
  )
]

// --- Section 2: Central Biological Conclusion ---
== 2. The Central Biological Conclusion

=== Suppression of the AMPK Signaling Axis Output
The central biological finding across human clinical cohorts and rodent diabetic models is the *functional suppression of AMPK pathway output activity* under chronic hyperglycemic/insulin-resistant stress:

#v(2pt)
#align(center)[
  #block(
    fill: accent-bg,
    stroke: 1pt + rgb("bfdbfe"),
    radius: 6pt,
    inset: (x: 14pt, y: 9pt),
    [
      #text(size: 10.5pt, weight: "bold", fill: primary)[
        $ "Suppressed AMPK Activity" ==> cases(
          bold("PPARGC1A") " (PGC-1"alpha")" arrow.b quad &"Loss of mitochondrial biogenesis & metabolic control",
          bold("CCND1") " (Cyclin D1)" arrow.t quad &"De-repressed cellular stress & aberrant proliferation signal"
        ) $
      ]
    ]
  )
]
#v(2pt)

*Key Mechanistic Insights:*
1. *Post-Translational Regulatory Control:* The human KEGG network displays discordant edge dynamics: the AMPK kinase subunit (*PRKAA2*) is transcriptionally maintained while its activated downstream effectors fall. AMPK activity is governed by *Thr172 phosphorylation* and cellular energy charge, not transcript abundance.
2. *Epithelial & Microvascular Compartmentalization:*
   - In *liver*, metabolic suppression (*Ppargc1a* downregulation) localizes specifically to epithelial *cholangiocytes* (per-mouse $P = 0.011$, detection dropping from 31.7% to 11.4%).
   - In *kidney*, stress/cell-cycle de-repression (*Ccnd1* upregulation) localizes to *microvascular endothelial cells* (pseudobulk $P = 1.2 times 10^(-6)$).
3. *Cross-Species Concordance:* The directional shifts of both primary anchors (*PPARGC1A* down, *CCND1* up) replicate across human clinical cohorts and mouse diabetic models.

#pagebreak()

// --- Section 3: Target Hierarchy & Evidence Hierarchy ---
== 3. Core Target Tiers & Evidence Hierarchy

Targets are classified mechanically into prioritized tiers based on evidence convergence across all five analytical stages:

#v(4pt)
#align(center)[
  #block(
    fill: light-bg,
    stroke: 1pt + border-col,
    radius: 6pt,
    inset: 10pt,
    width: 100%,
    [
      #grid(
        columns: (1fr),
        align: center,
        gutter: 6pt,
        [
          #block(fill: tier3-bg, stroke: 1pt + rgb("d8b4fe"), radius: 4pt, inset: 6pt, width: 70%)[
            #text(weight: "bold", size: 9pt, fill: tier3-text)[UPSTREAM MECHANISM (Tier 3)]\
            #text(size: 8pt, fill: dark)[*p-AMPKα (PRKAA1/2) · SIRT1 · STK11 (LKB1)*\ (Phosphorylation / Post-Translational Node)]
          ]
        ],
        [#text(fill: rgb("94a3b8"), size: 10pt)[$arrow.b$]],
        [
          #grid(
            columns: (1fr, 1fr),
            gutter: 14pt,
            [
              #block(fill: tier1-bg, stroke: 1pt + rgb("86efac"), radius: 4pt, inset: 6pt, width: 100%)[
                #text(weight: "bold", size: 9pt, fill: tier1-text)[PRIMARY TARGET 1 (Tier 1)]\
                #text(weight: "bold", size: 8.5pt, fill: dark)[PPARGC1A (PGC-1α)]\
                #text(size: 7.5pt, fill: rgb("475569"))[Master Metabolic Co-activator\ Human T2D: $log_2"FC" = -0.61$\ Mouse Cholangiocytes: $P = 0.011$ (ML Rank 3)]
              ]
            ],
            [
              #block(fill: tier1-bg, stroke: 1pt + rgb("86efac"), radius: 4pt, inset: 6pt, width: 100%)[
                #text(weight: "bold", size: 9pt, fill: tier1-text)[PRIMARY TARGET 2 (Tier 1)]\
                #text(weight: "bold", size: 8.5pt, fill: dark)[CCND1 (Cyclin D1)]\
                #text(size: 7.5pt, fill: rgb("475569"))[Pathway Output / Stress Readout\ Human T2D: $log_2"FC" = +0.80$\ Mouse Kidney EC: $P = 1.2 times 10^(-6)$ (Replicated)]
              ]
            ]
          )
        ],
        [#text(fill: rgb("94a3b8"), size: 10pt)[$arrow.b$]],
        [
          #grid(
            columns: (1fr, 1fr),
            gutter: 14pt,
            [
              #block(fill: tier2-bg, stroke: 1pt + rgb("fde68a"), radius: 4pt, inset: 5pt, width: 100%)[
                #text(weight: "bold", size: 8pt, fill: tier2-text)[READOUT: GLYCOLYSIS (Tier 2)]\
                #text(size: 8pt, fill: dark)[*PFKFB3* (FDR = 0.042)]
              ]
            ],
            [
              #block(fill: tier2-bg, stroke: 1pt + rgb("fde68a"), radius: 4pt, inset: 5pt, width: 100%)[
                #text(weight: "bold", size: 8pt, fill: tier2-text)[READOUT: CELL-CYCLE ARREST (Tier 2)]\
                #text(size: 8pt, fill: dark)[*CDKN1A (p21) · BAX · PHLDA3*]
              ]
            ]
          )
        ]
      )
    ]
  )
]

#v(4pt)

=== Target Classification Table
#table(
  columns: (1.2fr, 0.9fr, 1.8fr, 2.3fr, 1.8fr),
  stroke: (x, y) => if y == 0 { (bottom: 1.5pt + primary) } else { 0.5pt + border-col },
  fill: (col, row) => if row == 0 { light-bg } else if calc.even(row) { rgb("fcfdfd") } else { white },
  align: (col, row) => (left + horizon),
  
  table.header(
    [*Target*], [*Tier*], [*Role & Mechanism*], [*Key Supporting Evidence*], [*Planned In Vivo Assay*]
  ),
  
  [*PPARGC1A*\ (PGC-1α)],
  badge("Tier 1 (Primary)", bg: tier1-bg, text-color: tier1-text),
  [Master transcriptional coactivator of mitochondrial biogenesis & respiration.],
  [• Human liver: $log_2"FC" = -0.61, P = 6.0 times 10^(-4)$\ • Mouse cholangiocytes: $P = 0.011$, cell BH $P = 2.0 times 10^(-8)$\ • ML Consensus: Rank 3 / 19],
  [• RT-qPCR & Western blot (liver)\ • IF co-stain: PGC-1α × CK19],
  
  [*CCND1*\ (Cyclin D1)],
  badge("Tier 1 (Primary)", bg: tier1-bg, text-color: tier1-text),
  [De-repressed downstream output of AMPK ⊣ HuR/ELAVL1 axis.],
  [• Human liver: $log_2"FC" = +0.80, P = 7.5 times 10^(-4)$\ • Replicated in cohort GSE23343 ($P = 0.016$, AUC = 0.81)\ • Mouse kidney EC: $P = 1.2 times 10^(-6)$],
  [• RT-qPCR & Western blot (kidney & liver)\ • IF co-stain: Cyclin D1 × CD31],
  
  [*CDKN1A*\ (p21)],
  badge("Tier 2 (Readout)", bg: tier2-bg, text-color: tier2-text),
  [Cell-cycle arrest / senescence effector downstream of stress.],
  [• Human liver DEG: $log_2"FC" = +0.96, "FDR" = 0.016$\ • Elevated in mouse cholangiocytes ($"FDR" = 0.009$) & kidney EC ($"FDR" = 0.008$)],
  [• RT-qPCR & Western blot\ • IF co-staining with CD31 in kidney],
  
  [*PFKFB3*],
  badge("Tier 2 (Readout)", bg: tier2-bg, text-color: tier2-text),
  [Inducible glycolytic regulator connected to AMPK.],
  [• Human liver: $"FDR" = 0.042$\ • Mouse liver LSEC ($P = 0.049$) & kidney LoH (BH $P = 0.03$)],
  [• RT-qPCR & Western blot (liver & kidney)],
  
  [*PRKAA1/2*\ (AMPKα)],
  badge("Tier 3 (Master)", bg: tier3-bg, text-color: tier3-text),
  [Central metabolic sensor and upstream master kinase.],
  [• Upstream driver on KEGG hsa04152; acts via Thr172 phosphorylation.],
  [• Western blot for *p-AMPKα (Thr172) / total AMPKα*\ • Metformin control]
)

#v(8pt)

// --- Section 4: Multi-Modal Evidence Matrix ---
== 4. Multi-Modal Evidence Matrix

#table(
  columns: (1.3fr, 1.2fr, 1fr, 1.2fr, 1.4fr, 1.1fr, 1fr, 1.2fr, 1fr),
  stroke: (x, y) => if y == 0 { (bottom: 1.5pt + primary) } else { 0.5pt + border-col },
  fill: (col, row) => if row == 0 { light-bg } else if calc.even(row) { rgb("fcfdfd") } else { white },
  align: (col, row) => (left + horizon),
  
  table.header(
    [*Gene*], [*Human Bulk*], [*WGCNA*], [*KEGG Map*], [*Mouse scRNA*], [*Replication*], [*ML Rank*], [*cis-MR*], [*Tier*]
  ),
  
  [*PPARGC1A*], [Down ($P = 6.0 times 10^(-4)$)], [Key Module], [Direct (+1)], [Down (Cholang. $P=0.011$)], [Concordant], [Rank 3], [Null ($P = 0.16$)], badge("Tier 1", bg: tier1-bg, text-color: tier1-text),
  [*CCND1*], [Up ($P = 7.5 times 10^(-4)$)], [Key Module], [2-Step (−1)], [Up (Kidney EC $P=1.2 times 10^(-6)$)], [Replicated], [Rank 14], [Null ($P = 0.55$)], badge("Tier 1", bg: tier1-bg, text-color: tier1-text),
  [*CDKN1A*], [Up ($"FDR" = 0.016$)], [Key Module], [p53 Map], [Up (Cholang. & EC)], [Concordant], [—], [Null], badge("Tier 2", bg: tier2-bg, text-color: tier2-text),
  [*PFKFB3*], [Down ($"FDR" = 0.042$)], [Turquoise], [Direct (+1)], [Down (LSEC $P=0.049$)], [Concordant], [—], [Null], badge("Tier 2", bg: tier2-bg, text-color: tier2-text),
  [*IRS2*], [Down ($P = 0.011$)], [Key Module], [Input Node], [Minor cells only], [Concordant], [Rank 6], [Null], badge("Tier 2", bg: tier2-bg, text-color: tier2-text),
  [*IGF1*], [Down ($P = 0.021$)], [Key Module], [Input Node], [Sparse (Hepatocyte gap)], [Concordant], [Rank 8], [Null], badge("Tier 2b", bg: rgb("f1f5f9"), text-color: rgb("475569")),
  [*PRKAA2*], [Up ($P = 0.017$)], [Network Node], [Master Kinase], [Post-translational], [Concordant], [—], [Null], badge("Tier 3", bg: tier3-bg, text-color: tier3-text),
  [*SREBF1*], [Down ($P = 0.041$)], [Network Node], [Lipogenic], [Inverted (STZ model)], [—], [—], [Causal ($P=9 times 10^(-9)$)], badge("Exploratory", bg: rgb("f1f5f9"), text-color: rgb("475569")),
  [*VSNL1/LZTFL1*], [DE ($"FDR" < 0.05$)], [Key Module], [Off-Map], [Not DE in mouse], [Weak], [Rank 1 & 2], [Null], badge("Non-Target", bg: rgb("fee2e2"), text-color: rgb("991b1b"))
)

#pagebreak()

// --- Section 5: In Vivo Validation Plan ---
== 5. In Vivo Experimental Validation & Pharmacological Rescue Plan

The multi-omics analysis directly informs the experimental design, dosing groups, and molecular endpoints for the in vivo pre-clinical validation study in the Ob-T2D rodent model:

#v(4pt)

#block(
  fill: light-bg,
  stroke: 1pt + border-col,
  radius: 6pt,
  inset: 12pt,
  width: 100%,
  [
    #grid(
      columns: (1fr, 1.2fr),
      gutter: 14pt,
      [
        #text(weight: "bold", size: 9.5pt, fill: primary)[Experimental Groups (n = 8 per group)]
        #v(3pt)
        1. *Normal Control:* Standard diet + Vehicle (Saline).
        2. *Ob-T2D Disease Model:* High-Fat Diet (HFD) + STZ.
        3. *Ob-T2D + GV1 Low Dose:* Test formulation (Low).
        4. *Ob-T2D + GV1 Medium Dose:* Test formulation (Med).
        5. *Ob-T2D + GV1 High Dose:* Test formulation (High).
        6. *Ob-T2D + Metformin:* Positive control (200 mg/kg).
      ],
      [
        #text(weight: "bold", size: 9.5pt, fill: primary)[Molecular & Phenotypic Endpoints]
        #v(3pt)
        - *Primary Metabolic Efficacy:* Oral Glucose Tolerance Test (OGTT), Fasting Blood Glucose, HbA1c, HOMA-IR.
        - *Hepatic AMPK Activation:* Western blot for *p-AMPKα (Thr172) / total AMPKα*.
        - *Mitochondrial Function:* *PGC-1α* (RT-qPCR & Western blot).
        - *Cellular Stress / Proliferation:* *Cyclin D1 & p21* (qPCR & WB).
        - *Immunofluorescence:* PGC-1α × CK19 (liver cholangiocytes); Cyclin D1 × CD31 (kidney microvascular endothelium).
      ]
    )
  ]
)

#v(8pt)

=== Hypothesized Therapeutic Mechanism of Action

Administration of the active herbal/natural formulation (*GV1*) is hypothesized to mirror the therapeutic action of *Metformin* through multi-target pathway restoration:

#v(4pt)

#align(center)[
  #block(
    fill: accent-bg,
    stroke: 1pt + rgb("bfdbfe"),
    radius: 6pt,
    inset: 10pt,
    width: 100%,
    [
      #grid(
        columns: (1fr, auto, 1fr, auto, 1fr),
        align: center + horizon,
        gutter: 6pt,
        [
          #block(fill: white, stroke: 0.5pt + rgb("93c5fd"), radius: 4pt, inset: 6pt)[
            #text(weight: "bold", size: 8.5pt, fill: primary)[1. AMPK Activation]\
            #text(size: 7.5pt, fill: dark)[Restoration of Thr172 phosphorylation in liver & kidney]
          ]
        ],
        [#text(fill: primary, size: 12pt)[$arrow.r$]],
        [
          #block(fill: white, stroke: 0.5pt + rgb("93c5fd"), radius: 4pt, inset: 6pt)[
            #text(weight: "bold", size: 8.5pt, fill: primary)[2. PGC-1α Upregulation]\
            #text(size: 7.5pt, fill: dark)[Mitochondrial biogenesis & metabolic rescue]
          ]
        ],
        [#text(fill: primary, size: 12pt)[$arrow.r$]],
        [
          #block(fill: white, stroke: 0.5pt + rgb("93c5fd"), radius: 4pt, inset: 6pt)[
            #text(weight: "bold", size: 8.5pt, fill: primary)[3. Stress Normalization]\
            #text(size: 7.5pt, fill: dark)[Suppression of Cyclin D1 & p21 to baseline levels]
          ]
        ]
      )
    ]
  )
]

#v(10pt)
#line(length: 100%, stroke: 0.5pt + rgb("cbd5e1"))
#v(4pt)
#text(size: 8pt, fill: rgb("64748b"))[
  *Document Provenance:* This report synthesizes outputs from `24_target_scoring.R`, `16_pathway_selection.R`, `35_fig3_paper_style.R`, `21_ML_classifier.R`, and `19_mendelian_randomization.R`. All supporting data tables, figures, and analysis scripts are cataloged in `docs/target_selection.md` and `docs/pathway_selection.md`.
]
