# Figures: Reference-Framework Panels Only

This folder holds only the panels that mirror the analytical framework of the reference paper (Xu et al., *Phytomedicine* 2026, Figs 1–3). Everything else lives in `results/`:

| Location | Contents |
|---|---|
| `results/supplementary_figures/` | QC, diagnostic and sensitivity plots (bulk QC, scRNA demultiplexing/QC, pseudobulk DE counts, validation, prioritization) |
| `results/figure_exports/` | TIFF submission rasters, serialized panel intermediates (`fig*_panels/`), figure caches |
| `archive/superseded_figures/` | Outputs from superseded designs (3-group bulk contrast, 12-gene DEG run) |

## Publication contract

Everything in `figures/` obeys one contract, defined in `scripts/figure_style.R` and applied by the
figure scripts 27–31:

- **PDF only.** No PNG. TIFF submission rasters and any serialized intermediates go to
  `results/figure_exports/`, never here.
- **Drawn at final print size** in millimetres (single column 89 mm, 1.5 column 120 mm, double
  180 mm). Panels are never drawn large and scaled down — that is what lands text at 3 pt in a proof.
- **Arial 5–7 pt**, fonts embedded, text left live and selectable in the PDF; only dense point clouds
  are rasterised (600 dpi).
- **No titles, subtitles or explanatory prose inside the axes.** Every statistic that would sit in a
  panel title lives in that figure's `Fig<n>_legend.md`, which also carries the caveats.
- **Perceptually uniform, colour-vision-safe colour maps** (viridis for continuous, a blue–white–red
  divergent for correlations). No rainbow/jet.
- Panel tags are Nature-style bold lowercase; set `TAG_CASE <- "upper"` in `figure_style.R` for Cell.

Rebuild any figure with `Rscript scripts/<script> <stage>`, where `stage` is a single panel letter,
`panels`, `assemble`, `legend` or `all`. The scripts **re-draw only** — they read the result tables the
analysis scripts already wrote, so no number changes.

| Figure | Script | Legend |
|---|---|---|
| Fig 1 | `scripts/28_pub_fig1.R` | `Fig1_bulk_DEG_WGCNA/Fig1_legend.md` |
| Fig 2 | `scripts/29_pub_fig2.R` | `Fig2_enrichment/Fig2_legend.md` |
| Fig 3 atlas, per tissue | `scripts/27_fig3_composite.R <tissue>` | `Fig3_composite_<tissue>/Fig3_legend.md` |
| Fig 3 axis, per tissue | `scripts/32_pub_fig3_axis.R <tissue>` | `Fig3_axis_<tissue>/Fig3_axis_legend.md` |
| Fig 4 | `scripts/30_pub_fig4.R` | `Fig4_ML/Fig4_legend.md` |
| Fig 5 | `scripts/31_pub_fig5.R` | `Fig5_MR/Fig5_legend.md` |

## The reference-paper-style set (scripts 33-37)

There are **two complete figure sets in this folder**, and they are alternatives:

| | Contract set | Reference-paper-style set |
|---|---|---|
| Folders | `Fig1_bulk_DEG_WGCNA/`, `Fig2_enrichment/`, `Fig3_composite_<tissue>/`, `Fig3_axis_<tissue>/`, `Fig4_ML/`, `Fig5_MR/` | `Fig1_paper_style/`, `Fig2_paper_style/`, `Fig3_paper_style_<tissue>/`, `Fig4_paper_style/`, `Fig5_paper_style/` |
| Style | `scripts/figure_style.R`: Nature/Cell, lowercase tags, viridis, no in-panel titles | the reference paper's own look: boxed panels, bold centred panel titles, upper-case tags, its palettes |
| Scripts | 27-32 | 33-37 |

**Ship one set or the other, never both.** Both re-draw from the same result tables, so no number
differs between them; only the drawing does. Each paper-style folder carries its own
`*_paper_style_legend.md` with the caption written in the reference's wording and an explicit list of
the places this data forces a deviation from it.

| Figure | Script | Reference analogue |
|---|---|---|
| Fig 1 | `scripts/33_fig1_paper_style.R` | Xu et al. Fig. 1, panel for panel (A-G) |
| Fig 2 | `scripts/34_fig2_paper_style.R` | Xu et al. Fig. 2, panel for panel (A-D) |
| Fig 3 | `scripts/35_fig3_paper_style.R <tissue>` | Xu et al. Fig. 3, panel for panel (A-G) |
| Fig 4 | `scripts/37_fig4_paper_style.R` | **none** - the reference runs no classifier; visual language only |
| Fig 5 | `scripts/36_fig5_paper_style.R` | Xu et al. Fig. 4, panel for panel (A-D) |

The deviations that change what a reader should conclude, rather than just how it looks:

- **Fig 1B** uses nominal P. 84 genes reach FDR < 0.05, but only 31 also clear |log2FC| > 0.5 and the
  strict downstream rules then leave one intersecting gene (R17).
- **Fig 1E** ends on MEgrey, as the reference does. Grey is recomputed for display only; script 03
  keeps the unassigned bin out of the key-module rule.
- **Fig 2A** has no CC column: over-representation of 8 genes returns no cellular-component term.
- **Fig 2B, 2C** rest on the 4 candidates that map into KEGG at all, so the pathway list is one signal
  seen through 21 annotation sets.
- **Fig 3D, 3E, 3G** are computed inside the key cell type, not the whole atlas.
- **Fig 5** does **not** nominate a causal druggable target, which is the opposite of the reference's
  Fig. 4 result. Panel A therefore carries three exposures rather than one.

Analysis scripts never write into `figures/`. Scripts 01-05, 08, 15b, 18, 19, 20b, 21 and 23 send their
draft panels to `results/supplementary_figures/`, because Windows filenames are case-insensitive and a
draft `Fig1C_soft_threshold.pdf` silently overwrites the publication `Fig1c_soft_threshold.pdf`. The one
exception is the KEGG pathview diagram, which no figure script redraws.

**Fig 3 is now four figures**, two per tissue. `Fig3_composite_<tissue>/` is the atlas (UMAP, markers,
composition, the pair, co-expression density, correlation) and `Fig3_axis_<tissue>/` is the pathway-axis
evidence (AMPK-PGC-1a and p53 arrest: gene-level effects, per-cell-type scores, and each pair against a
detection-matched null). The old `Fig3_scRNA/` set is superseded and moved to
`archive/superseded_figures/Fig3_scRNA/`.

**Note on filenames:** panel files are lowercase (`Fig1a_…`). Windows is case-insensitive, so a new
lowercase name silently overwrites an old uppercase one when the rest of the name matches — check
before deleting anything that differs only in case.

## Fig 1: Bulk DEG + WGCNA (`Fig1_bulk_DEG_WGCNA/`)

Key modules follow **decision M11b**: T2D p < 0.1 pooled AND the same direction at p < 0.1 in each
cohort separately. That gives 4 modules (red, yellow, tan, greenyellow), 418 module genes and an
**8-gene** three-way candidate set. See `Fig1_legend.md` and `docs/decisions_log.md`.

| Panel | Reference paper | This project |
|---|---|---|
| a | PCA before/after batch correction | PCA before/after ComBat (GSE15653 + GSE64998), 16 T2D vs 11 lean |
| b | Volcano of DEGs | 96 up / 82 down at P < 0.05 and \|log2FC\| > 0.5 (M9; 84 genes reach FDR < 0.05, only 31 of them at \|log2FC\| > 0.5) |
| c | Scale independence / mean connectivity | Same; soft power 7, as in the paper |
| d | Module dendrogram | Same |
| e | Module–trait heatmap | T2D pooled **and per cohort**, HbA1c, Dataset†. Key modules in bold |
| f | MM vs GS in the key module | MM vs GS in all 4 key modules |
| g | Venn: DEG ∩ module ∩ phenotype set | DEGs ∩ WGCNA key-module genes ∩ hyperglycaemia set = **8** |

† The Dataset column is computed on ComBat-corrected expression, so it is bounded near zero by
construction and is reported, not used as a filter. The per-cohort columns are what M11b tests.

## Fig 2: Functional enrichment (`Fig2_enrichment/`)

| Panel | Reference paper | This project |
|---|---|---|
| a | GO enrichment, coloured by ontology | GO over-representation of the 8 candidates, faceted by ontology, colour = q |
| b | KEGG bubble plot | KEGG over-representation; top hit **AMPK signalling**, q = 5.8 × 10⁻⁶ |
| c | Pathway interaction network, top 30 KEGG | Same: top 30 KEGG terms, edges at Jaccard ≥ 0.5 |
| d | pathview diagram | GSEA Hallmark on the full ranking (the pathview raster is in `Fig2D_pathview/`) |

## Fig 3: scRNA-seq, STZ vs Control

Four figures, two per tissue. Selected map **AMPK signalling** (script 16); reported axis is
candidate-anchored (R23): **AMPKα2 (PRKAA2) → PGC-1α (PPARGC1A)** with **CCND1** as output, and the
downstream **Cdkn1a (p21) – Ccnd1** arrest readout.

### `Fig3_composite_<tissue>/` — atlas (script 27)

Liver 16,861 cells / 13 cell types; kidney 6,224 / 14. `stage` = `cache`, `a`…`g`, `panels`,
`assemble`, `legend`, `all`. Only `cache` touches the Seurat object.

| Panel | Contents |
|---|---|
| a | UMAP by cell type, legend gives per-type proportion, corner arrows in place of axes |
| b | Canonical markers (script-08 sets, 3 per present cell type) |
| c | Control/STZ fraction per cell type + cells per cell type (log10) |
| d, e | *Cdkn1a* and *Ccnd1* in the key cell type, Control vs STZ; white points are per-mouse means, the % under each violin is the detection rate |
| f | Nebulosa joint density for *Cdkn1a*+*Ccnd1* |
| g | Cell-level Pearson in the key cell type against the detection-matched null |

Brackets in **d**,**e** report the **per-mouse test** (n = 3 vs 4 mice), not the cell-level Wilcoxon;
cells within a mouse are not independent replicates. Both are in `Fig3_legend.md`; `SIG_UNIT <- "cell"`
switches the panels back. Panels d, e and g are all restricted to the key cell type, as in the
reference framework.

### `Fig3_axis_<tissue>/` — pathway-axis evidence (script 32)

Replaces the old `Fig3_scRNA/` D–J panels, redrawn from the tables scripts 18 and 23 wrote.

| Panel | Contents |
|---|---|
| a | AMPK–PGC-1α axis genes in the key cell type: pseudobulk log2FC, colour = agreement with the human direction, shape = significance, `*` = Fig 1g candidate |
| b | UCell score of the candidates on the AMPK map, per cell type |
| c | *Ppargc1a*–*Ccnd1* cell-level R against 2,000 detection-matched random pairs |
| d, e, f | As a, b, c for the p53 arrest axis and *Cdkn1a*–*Ccnd1* |

Neither pair is co-expressed above chance (liver: *Ppargc1a*–*Ccnd1* empirical P = 0.600,
*Cdkn1a*–*Ccnd1* P = 0.992). The axis claim rests on the gene-level concordance in **a** and **d**.

## Fig 4: Machine learning (`Fig4_ML/`)

Not a reference-paper panel — an addition (script 21, `docs/results_ML.md`, R22), built so the small-n
failure mode is visible.

| Panel | Contents |
|---|---|
| a | Panel-fixed CV vs nested CV vs label-permutation null, for LASSO / RF / SVM |
| b | ROC in the independent cohort GSE23343 (10 T2D / 7 NGT), three models plus a signature score |
| c | LASSO selection frequency and RF permutation importance |
| d | External predicted probability by true group |
| e | Candidate genes z-scored in both cohorts |

Read as a negative result: external AUC 0.43–0.64, every CI spanning 0.5.

## Fig 5: Mendelian randomization (`Fig5_MR/`)

Exposure: eQTLGen whole-blood cis-eQTLs. Outcome: T2D (GCST006867). **All panels are the LD-clumped
primary analysis** (R21). The instrumented set is the Fig 1g candidates plus the p53 and AMPK axis
genes, not the candidate set alone. Details in `docs/results_MR.md`.

| Panel | Reference paper | This project |
|---|---|---|
| a | Forest of the causal estimate | All 21 genes with usable instruments, OR per SD of predicted expression |
| b | Radial / SNP plot | SNP effect on expression vs on T2D for the 4 best-instrumented genes, with the IVW slope |
| c | Per-SNP forest | Leave-one-out IVW |
| d | — (added) | LD-clumped vs distance-pruned ORs; why the first pass over-called |
| e | — (added) | GTEx liver instruments (kidney cortex yields none) |
| f | Steiger directionality | Exposure vs outcome variance explained, per gene |

Of the 8 candidates only PPARGC1A, CCND1, IRS2 and IGFBP2 have instruments, and none is significant.
