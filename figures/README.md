# Figures

One figure set, drawn to match the reference paper: **Xu et al., *Phytomedicine* 154 (2026) 158050**.
Six folders, one per figure. Nothing else belongs here.

| Folder | Figure | Reference analogue | Script |
|---|---|---|---|
| `Fig1/` | Bulk DEG + WGCNA | Xu et al. Fig. 1, panel for panel (A–G) | `scripts/33_fig1_paper_style.R` |
| `Fig2/` | Functional enrichment | Xu et al. Fig. 2, panel for panel (A–D) | `scripts/34_fig2_paper_style.R` |
| `Fig3_liver/`, `Fig3_kidney/` | scRNA-seq atlas | Xu et al. Fig. 3, panel for panel (A–G) | `scripts/35_fig3_paper_style.R <tissue>` |
| `Fig4/` | Machine-learning classifier | **none** — the reference runs no classifier; visual language only | `scripts/37_fig4_paper_style.R` |
| `Fig5/` | Mendelian randomisation | Xu et al. Fig. 4, panel for panel (A–D) | `scripts/36_fig5_paper_style.R` |

Each folder holds the assembled figure (`Fig<n>.pdf`), one PDF per panel, and `Fig<n>_legend.md` —
the caption in the reference's own wording, plus an explicit list of where this data forces a
deviation from it. Read the legend before quoting any number from a panel.

Rebuild with `Rscript scripts/<script> <stage>`, where `stage` is a panel letter, `panels`,
`assemble`, `legend` or `all`. The scripts **re-draw only** — they read the result tables the analysis
scripts already wrote, so rebuilding changes no number.

```
Rscript scripts/33_fig1_paper_style.R all
Rscript scripts/34_fig2_paper_style.R all
Rscript scripts/35_fig3_paper_style.R liver all
Rscript scripts/35_fig3_paper_style.R kidney all
Rscript scripts/37_fig4_paper_style.R all
Rscript scripts/36_fig5_paper_style.R all
```

TIFF submission rasters and serialized panel intermediates go to `results/figure_exports/`, never
here. Analysis scripts never write into `figures/`.

## The chain, and where it holds

The reference's logic is: DEGs ∩ module ∩ phenotype set → enrichment → **one selected KEGG map** →
single-cell and MR on genes **taken from that map**. This set follows the same chain.

| Stage | Reference | This project |
|---|---|---|
| Candidate genes | 159 | **8** (DEG ∩ key modules ∩ hyperglycaemia set) |
| Selected map | regulation of actin cytoskeleton | **AMPK signalling (hsa04152)** |
| Genes of the set on that map | RAC1, PAK1 | **IRS2, PPARGC1A, CCND1, IGF1** |
| Single-cell pair | Rac1 + Pak1 | **Ppargc1a + Ccnd1** |
| ML panel | — | the 4 map genes above |
| MR exposures | PAK1 | the **9 map genes with usable instruments** |

Every downstream step draws from hsa04152. The four candidates that are not on it (VSNL1, SREBF2,
IGFBP1, IGFBP2) appear in **no** KEGG pathway at all and are carried in Fig. 1 only.

## What the results actually say

Read these before writing any text around the figures.

- **Fig. 1B uses nominal P.** 84 genes reach FDR < 0.05, but only 31 also clear |log2FC| > 0.5, and
  the strict downstream rules then leave a single intersecting gene (SREBF2, decision R17).
- **Fig. 2A has no CC column.** Over-representation of 8 genes returns BP and MF terms only.
- **Fig. 2B/2C rest on 4 genes.** The 21 pathways at adjusted P < 0.05 are recombinations of IRS2,
  PPARGC1A, CCND1 and IGF1 — one signal seen through 21 annotation sets, not 21 findings.
- **Fig. 3 is negative, on purpose.** The on-map pair is co-localized in liver cholangiocytes (13.8 %
  of cells) but not correlated (R 0.03, empirical P 0.60). In kidney endothelium it fails the 2 %
  co-localization gate outright at 0.17 %, because *Ppargc1a* is barely expressed there. The only
  pair in either atlas clearing all three of the reference's gates is the off-map p53 readout
  *Cdkn1a*–*Ccnd1* (kidney endothelium, empirical P 0.009); it is not part of this set because
  *Cdkn1a* is not a node of hsa04152.
- **Fig. 4 is negative, but restriction helps.** External AUC in GSE23343 is 0.71 for LASSO on the
  4-gene map panel, against 0.43 on the unrestricted 8-gene candidate panel. Every confidence
  interval still includes 0.5.
- **Fig. 5 is the one clear positive.** The selected map contains two causal genes at FDR < 0.05:
  **SREBF1** (OR 0.903, P 9.0e-9, q 8.1e-8) and **SIRT1** (OR 1.027, P 9.5e-4, q 4.3e-3). FDR is
  Benjamini–Hochberg over the 9 map genes, one primary estimate each. Neither Tier 1 target is
  causal (PPARGC1A P 0.16, CCND1 P 0.55), and SREBF1 rests on only 2 instruments, so its pleiotropy
  cannot be tested — which is why panels B and C feature SIRT1.

## Variants

The scripts can still draw the alternatives; they write to suffixed folders so the shipped set is
never overwritten.

```
# Fig 3 with the off-map p53 pair (the one that passes all three gates, in kidney)
Rscript scripts/27_fig3_composite.R kidney cache "Cdkn1a,Ccnd1"
Rscript scripts/35_fig3_paper_style.R kidney all "Cdkn1a,Ccnd1"

# Fig 4 on the unrestricted 8-gene candidate panel
Rscript scripts/21_ML_classifier.R          # writes results/ml/
Rscript scripts/37_fig4_paper_style.R all ""
```

The Nature/Cell-contract set (scripts 27–32, lowercase tags, viridis, no in-panel titles) is no
longer built. Those scripts still exist and still run; their output folders were removed in favour of
this one set. Recover them from git history if they are ever needed again.

## Elsewhere

| Location | Contents |
|---|---|
| `results/supplementary_figures/` | QC, diagnostic and sensitivity plots |
| `results/figure_exports/` | TIFF submission rasters, panel intermediates, figure caches |
| `results/pathway_selection/pathview/` | the KEGG pathview raster and KGML for hsa04152 (rendered by KEGG, not drawn by a figure script) |
| `archive/superseded_figures/` | outputs from superseded designs |

**Note on filenames:** Windows is case-insensitive, so a new lowercase name silently overwrites an old
uppercase one when the rest matches. Check before deleting anything that differs only in case.

**Do not open `results/bulk/enrich_KEGG.csv` in Excel.** It reads `GeneRatio` values such as `4/4` as
dates and writes back `4-Apr`, which makes every ratio unparseable. This has happened twice.
`scripts/34_fig2_paper_style.R` now refuses to run on a damaged file rather than silently dropping
every point; the repair is `GeneRatio = Count/4` for all 81 rows.
