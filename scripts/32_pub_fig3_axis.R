# 32_pub_fig3_axis.R - the pathway-axis evidence panels of Fig 3, redrawn to the publication contract.
#
# Usage: Rscript scripts/32_pub_fig3_axis.R <tissue> <stage>
#   tissue : liver | kidney
#   stage  : A | B | C | D | E | F | panels | assemble | legend | all
#
# Replaces the old figures/Fig3_scRNA/Fig3{D-G,H-J}_*.pdf, which were drawn on 300-430 mm canvases with
# in-panel titles. This script re-DRAWS only, from the tables scripts 18 and 23 already wrote:
#   a  AMPK-PGC-1a axis genes in the key cell type            <- 23_anchored_gene_stats_<tissue>.csv
#   b  candidate-gene score on the AMPK map, per cell type    <- 23_anchored_score_<tissue>.csv
#   c  anchored pair correlation against a matched null       <- 23_anchored_correlation_<tissue>.csv
#   d  p53 arrest axis genes in the key cell type             <- 18_axis_gene_stats_<tissue>.csv
#   e  p53 target-programme score, per cell type              <- 18_axis_score_<tissue>.csv
#   f  p53 pair correlation against a matched null            <- 18_axis_correlation_<tissue>.csv
#
# The old panels showed per-cell violins. Those need the Seurat object; the stored tables carry the
# pseudobulk effect sizes and tests, which is what the claim actually rests on, so the gene panels are
# drawn as effect-size plots with the cross-species concordance shown explicitly.
#
# The joint-density co-localization panels (old Fig3F/Fig3I) are not reproduced here: the composite
# figure from script 27 already carries the Cdkn1a+Ccnd1 density as its panel f, and the anchored pair
# adds only its co-expression percentage, which is annotated in panel c.
#
# Outputs: figures/Fig3_axis_<tissue>/Fig3ax<a-f>_*.pdf, Fig3_axis_<tissue>.pdf, Fig3_axis_legend.md
#          results/figure_exports/Fig3_axis_<tissue>.tiff, results/figure_exports/fig3ax_panels_<tissue>/

source("scripts/figure_style.R")

args   <- commandArgs(trailingOnly = TRUE)
tissue <- tolower(if (length(args) >= 1) args[1] else "liver")
stage  <- if (length(args) >= 2) args[2] else "all"
stopifnot(tissue %in% c("liver", "kidney"))
set.seed(20260918)

fig_dir   <- file.path("figures", paste0("Fig3_axis_", tissue))
panel_dir <- file.path(EXPORT_DIR, paste0("fig3ax_panels_", tissue))
psd       <- "results/pathway_selection"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
keep <- function(p, n, w, h) save_fig(p, n, w, h, dir = fig_dir, keep_dir = panel_dir)

rd <- function(f) fread(file.path(psd, sprintf(f, tissue)))

# ---------------------------------------------------------------------------------------------------
# gene effect panel: pseudobulk log2FC in the key cell type, with cross-species concordance
# ---------------------------------------------------------------------------------------------------
gene_panel <- function(d, mark_candidates = FALSE) {
  d <- copy(d)
  d[, concord := fifelse(concordant, "Same direction as human", "Discordant")]
  d[, sig := fifelse(pb_FDR < 0.05, "FDR < 0.05",
             fifelse(pb_P < 0.05, "P < 0.05", "n.s."))]
  d[, lab := gene]
  if (mark_candidates && "fig1_candidate" %in% names(d))
    d[fig1_candidate == TRUE, lab := paste0(gene, "*")]
  d[, lab := factor(lab, levels = d[order(pb_logFC), lab])]

  ggplot(d, aes(pb_logFC, lab)) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_segment(aes(x = 0, xend = pb_logFC, yend = lab), colour = "grey80", linewidth = 0.25) +
    geom_point(aes(colour = concord, shape = sig), size = 1.3) +
    scale_colour_manual(values = c(`Same direction as human` = pal_qual[1],
                                   Discordant = col_dir[["Up"]]), name = NULL) +
    scale_shape_manual(values = c(`FDR < 0.05` = 16, `P < 0.05` = 17, `n.s.` = 1), name = NULL) +
    scale_x_continuous(expand = expansion(mult = 0.12)) +
    labs(x = expression("Pseudobulk log"[2]*" fold change, STZ vs Control"), y = NULL) +
    theme_pub() +
    theme(axis.text.y = element_text(face = "italic", size = BASE - 1),
          legend.position = "top", legend.justification = "left",
          legend.box = "vertical", legend.spacing.y = unit(0.5, "mm"),
          legend.margin = margin(0, 0, 1, 0))
}

# ---------------------------------------------------------------------------------------------------
# score panel: per-cell-type difference in signature score
# ---------------------------------------------------------------------------------------------------
score_panel <- function(d, xlab) {
  d <- copy(d)
  d <- d[is.finite(delta)]
  d[, ct := factor(pretty_ct_local(cell_type), levels = pretty_ct_local(cell_type[order(delta)]))]
  d[, sig := fifelse(BH < 0.05, "FDR < 0.05", fifelse(p < 0.05, "P < 0.05", "n.s."))]
  ggplot(d, aes(delta, ct)) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_segment(aes(x = 0, xend = delta, yend = ct), colour = "grey80", linewidth = 0.25) +
    geom_point(aes(colour = sig), size = 1.3) +
    scale_colour_manual(values = c(`FDR < 0.05` = col_dir[["Up"]], `P < 0.05` = pal_qual[2],
                                   `n.s.` = "grey60"), name = NULL) +
    scale_x_continuous(expand = expansion(mult = 0.12)) +
    labs(x = xlab, y = NULL) +
    theme_pub() +
    theme(legend.position = "top", legend.justification = "left",
          legend.margin = margin(0, 0, 1, 0))
}
pretty_ct_local <- function(x) {
  x <- as.character(x)
  x <- sub("^HSC_Fibroblast$", "HSC/fibroblast", x)
  x <- sub("^Macrophage_other$", "Macrophage (other)", x)
  x <- sub("^Fibroblast_Pericyte$", "Fibroblast/pericyte", x)
  x <- sub("^Tubular_epithelium_other$", "Tubular epithelium (other)", x)
  x <- sub("^CD_IC$", "CD-IC", x)
  gsub("_", " ", x)
}

# ---------------------------------------------------------------------------------------------------
# correlation-vs-null panel
# ---------------------------------------------------------------------------------------------------
null_panel <- function(cr) {
  pair <- cr$pair[1]
  d <- data.table(
    what = factor(c("Observed", "Matched random pairs"),
                  levels = c("Matched random pairs", "Observed")),
    r    = c(cr$cell_R[1], cr$bg_mean_R[1]),
    hi   = c(NA_real_, cr$bg_q95[1]))
  lab <- sprintf('italic(R)*" = %s"', sub("-", "−", formatC(cr$cell_R[1], format = "f", digits = 3)))
  sub <- sprintf('"empirical "*italic(P)*" = %s"', formatC(cr$empirical_p[1], format = "f", digits = 3))
  ggplot(d, aes(r, what)) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_segment(aes(x = 0, xend = r, yend = what), colour = "grey80", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = r, xmax = hi), height = 0.12, linewidth = 0.3, colour = "grey45",
                   na.rm = TRUE) +
    geom_point(aes(colour = what), size = 1.6) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.7, parse = TRUE, label = lab,
             size = (BASE - 1) * ppt, family = FONT) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.07, vjust = 3.4, parse = TRUE, label = sub,
             size = (BASE - 2) * ppt, family = FONT) +
    scale_colour_manual(values = c(Observed = col_dir[["Up"]],
                                   `Matched random pairs` = "grey55"), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.25, 0.25))) +
    labs(x = sprintf("Cell-level Pearson R, %s", sub("-", "–", pair)), y = NULL) +
    theme_pub() +
    theme(axis.text.y = element_text(size = BASE - 1))
}

# ===================================================================================================
if (stage %in% c("A", "panels", "all"))
  keep(gene_panel(rd("23_anchored_gene_stats_%s.csv"), mark_candidates = TRUE),
       "Fig3axa_anchored_genes", 89, 62)

if (stage %in% c("B", "panels", "all")) {
  d <- rd("23_anchored_score_%s.csv")
  # the anchored score is computed per signature; the up-direction set is the one the axis claim uses
  if ("signature" %in% names(d) && uniqueN(d$signature) > 1) d <- d[signature == d$signature[1]]
  keep(score_panel(d, "UCell score, STZ − Control"), "Fig3axb_anchored_score", 89, 62)
}

if (stage %in% c("C", "panels", "all"))
  keep(null_panel(rd("23_anchored_correlation_%s.csv")), "Fig3axc_anchored_null", 89, 34)

if (stage %in% c("D", "panels", "all"))
  keep(gene_panel(rd("18_axis_gene_stats_%s.csv")), "Fig3axd_p53_genes", 89, 62)

if (stage %in% c("E", "panels", "all"))
  keep(score_panel(rd("18_axis_score_%s.csv"), "UCell p53-target score, STZ − Control"),
       "Fig3axe_p53_score", 89, 62)

if (stage %in% c("F", "panels", "all"))
  keep(null_panel(rd("18_axis_correlation_%s.csv")), "Fig3axf_p53_null", 89, 34)

# ===================================================================================================
if (stage %in% c("assemble", "all")) {
  need <- c("Fig3axa_anchored_genes", "Fig3axb_anchored_score", "Fig3axc_anchored_null",
            "Fig3axd_p53_genes", "Fig3axe_p53_score", "Fig3axf_p53_null")
  ps <- lapply(need, function(n) {
    f <- file.path(panel_dir, paste0(n, ".rds"))
    if (!file.exists(f)) stop("missing panel ", n, "; run stage panels first")
    readRDS(f)$plot
  })
  design <- "
AAAAAABBBBBB
AAAAAABBBBBB
AAAAAABBBBBB
CCCCCCCCCCCC
DDDDDDEEEEEE
DDDDDDEEEEEE
DDDDDDEEEEEE
FFFFFFFFFFFF
"
  comp <- ps[[1]] + ps[[2]] + ps[[3]] + ps[[4]] + ps[[5]] + ps[[6]] +
    plot_layout(design = design) +
    plot_annotation(tag_levels = list(tag_letters(6))) &
    theme(plot.tag = element_text(size = BASE + 1, face = "bold", family = FONT))
  save_fig(comp, paste0("Fig3_axis_", tissue), W_2COL, 195, dir = fig_dir, tiff = TRUE)
}

# ===================================================================================================
if (stage %in% c("legend", "panels", "all")) {
  an <- rd("23_anchored_gene_stats_%s.csv"); ac <- rd("23_anchored_correlation_%s.csv")
  pg <- rd("18_axis_gene_stats_%s.csv");     pc <- rd("18_axis_correlation_%s.csv")
  l <- c(
    sprintf("# Figure 3 (axis panels) legend - %s", tissue),
    "",
    sprintf("**Fig. 3 (continued) | Pathway-axis evidence in STZ-diabetic mouse %s.**", tissue),
    sprintf("**a**, AMPK-PGC-1a axis genes in %s, the cell type with the most human-concordant axis genes. Points are pseudobulk log2 fold changes (STZ vs Control, edgeR quasi-likelihood); colour marks agreement with the human T2D direction and shape marks significance. An asterisk marks a Fig. 1g candidate gene.",
            tolower(ac$cell_type[1])),
    "**b**, UCell score of the candidate genes lying on the AMPK map, per cell type, as the STZ minus Control difference of per-mouse means.",
    sprintf("**c**, Cell-level Pearson correlation of %s in %s against 2,000 detection-matched random gene pairs from the same cell type. Bar, the 95th percentile of the null.",
            sub("-", "–", ac$pair[1]), tolower(ac$cell_type[1])),
    sprintf("**d**,**e**,**f**, As **a**,**b**,**c** for the p53 arrest axis and the %s pair in %s.",
            sub("-", "–", pc$pair[1]), tolower(pc$cell_type[1])),
    "",
    "## Statistics to quote",
    "",
    "| Quantity | Value |",
    "|---|---|",
    sprintf("| Key cell type (anchored axis) | %s |", pretty_ct_local(ac$cell_type[1])),
    sprintf("| %s, cells | %s |", ac$pair[1], nfmt(ac$n_cells[1])),
    sprintf("| %s, co-expressing cells | %s%% |", ac$pair[1], mfmt(ac$pct_coexpr[1], 1)),
    sprintf("| %s, cell-level R | %s (%s) |", ac$pair[1], mfmt(ac$cell_R[1], 3), p_short(ac$cell_p[1])),
    sprintf("| %s, null mean / 95th pct / empirical P | %s / %s / %s |", ac$pair[1],
            mfmt(ac$bg_mean_R[1], 3), mfmt(ac$bg_q95[1], 3), mfmt(ac$empirical_p[1], 3)),
    sprintf("| %s, co-expressing cells | %s%% |", pc$pair[1], mfmt(pc$pct_coexpr[1], 1)),
    sprintf("| %s, cell-level R | %s (%s) |", pc$pair[1], mfmt(pc$cell_R[1], 3), p_short(pc$cell_p[1])),
    sprintf("| %s, null mean / 95th pct / empirical P | %s / %s / %s |", pc$pair[1],
            mfmt(pc$bg_mean_R[1], 3), mfmt(pc$bg_q95[1], 3), mfmt(pc$empirical_p[1], 3)),
    sprintf("| Axis genes concordant with human (anchored) | %d of %d |",
            an[concordant == TRUE, .N], nrow(an)),
    sprintf("| Axis genes concordant with human (p53) | %d of %d |",
            pg[concordant == TRUE, .N], nrow(pg)),
    "",
    "## Caveats that belong in the text, not the figure",
    "",
    sprintf("- Neither pair is co-expressed above chance: %s gives empirical P = %s and %s gives P = %s against detection-matched null pairs. The axis claim rests on the gene-level concordance in **a** and **d**, not on co-expression.",
            ac$pair[1], mfmt(ac$empirical_p[1], 3), pc$pair[1], mfmt(pc$empirical_p[1], 3)),
    "- Pseudobulk tests use 3 Control and 4 STZ mice; cell-type-level effects with few eligible mice are unstable, and the per-mouse cell counts are uneven (STZ_4 dominates the fibroblast compartment).",
    "- The cross-species direction check compares mouse pseudobulk against the human bulk contrast, which is an association, not a validation of the human effect."
  )
  writeLines(l, file.path(fig_dir, "Fig3_axis_legend.md"))
  message("wrote ", file.path(fig_dir, "Fig3_axis_legend.md"))
}
