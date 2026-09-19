# 35_fig3_paper_style.R - Figure 3 drawn in the visual style of the reference paper.
#
# Usage: Rscript scripts/35_fig3_paper_style.R <tissue> <stage>
#   tissue : liver | kidney
#   stage  : A | B | C | D | E | F | G | panels | assemble | legend | all   (default all)
#
# Reference: Xu et al., Phytomedicine 154 (2026) 158050, Fig. 3. This reproduces that figure's own
# visual language with this project's data, panel for panel:
#   A  UMAP of the atlas, cluster names written on the plot, "nCells: N" title, a
#      "Cell Type (Proportion)" legend
#   B  marker bubble plot, cell type on y and gene on x, size = percent expressed, colour = average
#      expression on a teal ramp, both keys laid out above the panel
#   C  left, stacked bars of the fraction of cells from each group; right, cells per cluster on a
#      log10 axis
#   D,E  violin of each pair gene by group, with a significance bracket
#   F  joint-density UMAP of the pair, titled "<gene1>+ <gene2>"
#   G  cell-level correlation of the pair in the key cell type, red fit over a grey ribbon, with the
#      R and P written in a box inside the panel
#
# This is deliberately NOT the figure_style.R contract used by scripts 27-32. It is the reference
# paper's style, kept as a separate deliverable in figures/Fig3_paper_style_<tissue>/ so the
# Nature/Cell version in figures/Fig3_composite_<tissue>/ is untouched. Pick one; do not ship both.
#
# Deviations from the reference forced by this dataset, all of them stated in the legend:
#   - Two groups, not three. The reference contrasts Ctrl / CIA / Treat; this atlas has Control and
#     STZ only, so D and E carry one bracket rather than three.
#   - The bracket reports the PER-MOUSE t-test (n = 3 vs 4 mice), not the cell-level Wilcoxon, which
#     is what script 27 does and for the same reason: cells within a mouse are not replicates. The
#     cell-level P is reported in the legend. Set SIG_UNIT <- "cell" to swap them.
#   - G is a within-cell-type correlation whose significance is judged against a background-matched
#     null (empirical P), not the raw Pearson P; both are given in the legend.
#
# Reads only results/figure_exports/fig3_composite_<tissue>_cache.rds, the cache script 27 writes, so
# it never loads the Seurat object. Run script 27 stage `cache` first if that file is missing.
#
# Outputs: figures/Fig3_paper_style_<tissue>/Fig3<A-G>_*.pdf, Fig3_paper_style_<tissue>.pdf,
#          Fig3_paper_style_legend.md, results/figure_exports/Fig3_paper_style_<tissue>.tiff

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
})
options(stringsAsFactors = FALSE)

TISSUES <- c("liver", "kidney")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || !tolower(args[1]) %in% TISSUES)
  stop("usage: Rscript scripts/35_fig3_paper_style.R <", paste(TISSUES, collapse = "|"),
       "> [stage]\n  ", if (length(args) < 1) "no tissue given" else
         paste0("unrecognised tissue: ", args[1]), call. = FALSE)
tissue <- tolower(args[1])
stage  <- if (length(args) >= 2) args[2] else "all"
set.seed(20260918)

SIG_UNIT <- "mouse"     # "mouse" (per-mouse t-test) or "cell" (cell-level Wilcoxon)

fig_dir   <- file.path("figures", paste0("Fig3_paper_style_", tissue))
panel_dir <- file.path("results/figure_exports", paste0("fig3paper_panels_", tissue))
cache_f   <- file.path("results/figure_exports", paste0("fig3_composite_", tissue, "_cache.rds"))
for (d in c(fig_dir, panel_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cache_f))
  stop("missing ", cache_f, "\n  run: Rscript scripts/27_fig3_composite.R ", tissue, " cache",
       call. = FALSE)
cc <- readRDS(cache_f)

# the paper's palette: a pastel qualitative ramp for cell types, green for the groups, a teal ramp on
# the bubble plot and a blue-yellow-red density scale
pal_ct   <- c("#F2F0A0", "#BFD7EA", "#E8C6A0", "#F2A6A2", "#A8A8A8", "#D9C6E0", "#F4C04E",
              "#BEE3B0", "#7FCDBB", "#C9B79C", "#E5A3C4", "#9ECAE1", "#CBD5A0")
cols_grp <- c(Control = "#1B7837", STZ = "#762A83")
col_bar  <- c(Control = "#C7E9C0", STZ = "#41AB5D")
col_dot_lo <- "#E4F1EE"; col_dot_hi <- "#20706B"
col_cnt  <- "#2E6E6A"

theme_baseR <- function(base = 8)
  theme_bw(base_size = base) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
        axis.ticks = element_line(colour = "black", linewidth = 0.3),
        plot.title = element_text(hjust = 0.5, face = "bold", size = base + 0.5),
        plot.margin = margin(2, 3, 2, 2))

keep <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, units = "mm",
         device = cairo_pdf)
  saveRDS(list(plot = p, width = w, height = h), file.path(panel_dir, paste0(name, ".rds")))
  message("wrote ", name, " (", w, " x ", h, " mm)")
  invisible(p)
}
nice_ct <- function(x) gsub("_", " ", x)
stars <- function(p) {
  if (is.na(p)) "ns" else if (p < 1e-4) "****" else if (p < 1e-3) "***" else
    if (p < 1e-2) "**" else if (p < 0.05) "*" else "ns"
}
# cc$gene_tests and cc$per_mouse in the cache are computed over the WHOLE atlas. D, E and G are about
# the key cell type only, so the tests are recomputed here inside that cell type, by the same rule
# script 27 uses, rather than read off the cache.
pair_stats <- function(g) {
  md <- cc$md[cell_type == cc$key]
  w  <- suppressWarnings(wilcox.test(md[[g]][md$condition == "STZ"],
                                     md[[g]][md$condition == "Control"]))
  pm <- md[, .(expr = mean(get(g)), n = .N), by = .(mouse_id, condition)][n >= 10]
  tt <- if (pm[condition == "STZ", .N] >= 2 && pm[condition == "Control", .N] >= 2)
    suppressWarnings(t.test(pm[condition == "STZ", expr], pm[condition == "Control", expr])) else NULL
  list(md = md, per_mouse = pm,
       cell_p  = w$p.value,
       mouse_p = if (is.null(tt)) NA_real_ else tt$p.value,
       det     = md[, .(pct = 100 * mean(get(g) > 0)), by = condition])
}

# the cell-type proportion table drives both the A legend and the cell-type ordering everywhere else
prop <- cc$md[, .N, by = cell_type][order(-N)]
prop[, pct := 100 * N / sum(N)]
prop[, lab := sprintf("%s (%.2f%%)", nice_ct(cell_type), pct)]
cols_ct <- setNames(rep_len(pal_ct, nrow(prop)), prop$cell_type)

# ===================================================================================================
# A - UMAP of the atlas, as in the paper
# ===================================================================================================
if (stage %in% c("A", "panels", "all")) {
  md <- copy(cc$md)
  md[, cell_type := factor(cell_type, levels = prop$cell_type)]
  cen <- md[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = cell_type]

  pA <- ggplot(md, aes(UMAP_1, UMAP_2, colour = cell_type)) +
    geom_point(size = 0.12, alpha = 0.75, stroke = 0) +
    geom_text_repel(data = cen, aes(label = nice_ct(cell_type)), colour = "black", size = 2,
                    seed = 1, min.segment.length = 0.25, segment.size = 0.15,
                    segment.colour = "grey45", box.padding = 0.2, max.overlaps = Inf) +
    scale_colour_manual(values = cols_ct, breaks = prop$cell_type, labels = prop$lab,
                        name = "Cell Type (Proportion)",
                        guide = guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
    labs(title = sprintf("nCells: %d", cc$n_cells), x = "UMAP_1", y = "UMAP_2") +
    theme_classic(base_size = 8) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 8.5),
          legend.title = element_text(size = 7), legend.text = element_text(size = 6),
          legend.key.size = unit(3, "mm"))
  keep(pA, "Fig3A_umap", 120, 88)
}

# ===================================================================================================
# B - marker bubble plot, as in the paper
# ===================================================================================================
if (stage %in% c("B", "panels", "all")) {
  d <- copy(cc$dot)
  d[, cell_type := factor(cell_type, levels = rev(prop$cell_type))]
  d[, gene := factor(gene, levels = sort(unique(gene)))]

  pB <- ggplot(d[pct_expr > 0], aes(gene, cell_type)) +
    geom_point(aes(size = pct_expr, colour = avg_expr)) +
    scale_colour_gradient(low = col_dot_lo, high = col_dot_hi, name = "Average Expression",
                          guide = guide_colourbar(barheight = unit(1.8, "mm"),
                                                  barwidth = unit(16, "mm"),
                                                  title.position = "top")) +
    scale_size_continuous(range = c(0.15, 1.9), name = "Percent Expressed",
                          breaks = c(0, 25, 50, 75, 100), limits = c(0, 100),
                          guide = guide_legend(title.position = "top", nrow = 1)) +
    scale_y_discrete(labels = nice_ct) +
    labs(x = NULL, y = NULL) +
    theme_baseR() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3.6),
          axis.text.y = element_text(size = 5),
          legend.position = "top", legend.box = "horizontal",
          legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
          legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -2, 0))
  keep(pB, "Fig3B_markers", 150, 62)
}

# ===================================================================================================
# C - group composition and cells per cluster, as in the paper
# ===================================================================================================
if (stage %in% c("C", "panels", "all")) {
  comp <- cc$md[, .N, by = .(cell_type, condition)]
  comp[, frac := N / sum(N), by = cell_type]
  comp[, cell_type := factor(cell_type, levels = rev(prop$cell_type))]
  comp[, condition := factor(condition, levels = names(col_bar))]
  cnt <- copy(prop)[, cell_type := factor(cell_type, levels = rev(prop$cell_type))]

  pC1 <- ggplot(comp, aes(frac, cell_type, fill = condition)) +
    # reverse = TRUE so the groups stack left to right in legend order, as the reference does;
    # position_stack() otherwise lays the last factor level down first
    geom_col(width = 0.72, colour = "black", linewidth = 0.12,
             position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = col_bar, name = "Group") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    scale_y_discrete(labels = nice_ct) +
    labs(x = "Fraction of cells in each group", y = NULL) +
    theme_baseR() +
    theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 5.5),
          legend.position = "top", legend.title = element_text(size = 6),
          legend.text = element_text(size = 5.5), legend.key.size = unit(2.6, "mm"),
          legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -2, 0))
  pC2 <- ggplot(cnt, aes(N, cell_type)) +
    geom_col(width = 0.72, fill = col_cnt) +
    scale_x_log10(expand = expansion(mult = c(0, 0.05)),
                  labels = scales::label_number(accuracy = 1, big.mark = "")) +
    labs(x = "Cells per cluster, log10 scale", y = NULL) +
    theme_baseR() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 5.5))
  pC <- wrap_elements(full = pC1 + pC2 + plot_layout(widths = c(1.6, 1)))
  keep(pC, "Fig3C_composition", 150, 62)
}

# ===================================================================================================
# D, E - violin of each pair gene by group, as in the paper
# ===================================================================================================
build_violin <- function(g) {
  st <- pair_stats(g)
  md <- copy(st$md); pm <- copy(st$per_mouse)
  p_show <- if (SIG_UNIT == "mouse") st$mouse_p else st$cell_p
  md[, condition := factor(condition, levels = names(cols_grp))]
  pm[, condition := factor(condition, levels = names(cols_grp))]
  ymax <- max(md[[g]])

  ggplot(md, aes(condition, .data[[g]], fill = condition)) +
    geom_violin(scale = "width", trim = TRUE, colour = "black", linewidth = 0.25, width = 0.82) +
    geom_point(data = pm, aes(condition, expr), shape = 21, size = 0.9, fill = "white",
               colour = "black", stroke = 0.3, inherit.aes = FALSE,
               position = position_jitter(width = 0.09, height = 0, seed = 1)) +
    annotate("segment", x = 1, xend = 2, y = ymax * 1.07, yend = ymax * 1.07, linewidth = 0.25) +
    annotate("segment", x = 1, xend = 1, y = ymax * 1.03, yend = ymax * 1.07, linewidth = 0.25) +
    annotate("segment", x = 2, xend = 2, y = ymax * 1.03, yend = ymax * 1.07, linewidth = 0.25) +
    annotate("text", x = 1.5, y = ymax * 1.10, label = stars(p_show), vjust = 0, size = 2.4) +
    scale_fill_manual(values = cols_grp, name = "Group") +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.17)),
                       breaks = scales::breaks_pretty(4)) +
    labs(x = NULL, y = bquote(italic(.(g)) ~ "Expression")) +
    theme_baseR() +
    theme(legend.position = "top", legend.title = element_text(size = 6),
          legend.text = element_text(size = 5.5), legend.key.size = unit(2.6, "mm"),
          legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -2, 0),
          axis.text.x = element_text(size = 6))
}
if (stage %in% c("D", "panels", "all")) keep(build_violin(cc$pair[1]), "Fig3D_violin1", 48, 62)
if (stage %in% c("E", "panels", "all")) keep(build_violin(cc$pair[2]), "Fig3E_violin2", 48, 62)

# ===================================================================================================
# F - joint-density UMAP of the pair, as in the paper
# ===================================================================================================
if (stage %in% c("F", "panels", "all")) {
  d <- cc$dens[order(density)]
  pF <- ggplot(d, aes(UMAP_1, UMAP_2, colour = density)) +
    geom_point(size = 0.12, stroke = 0) +
    scale_colour_gradientn(colours = c("#2B3A8F", "#3C8DBC", "#8CC9A8", "#F4E285", "#E8362B"),
                           name = "Joint density") +
    labs(title = sprintf("%s+ %s", cc$pair[1], cc$pair[2]), x = "umap 1", y = "umap 2") +
    theme_baseR() +
    theme(legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
          legend.key.size = unit(3, "mm"))
  keep(pF, "Fig3F_density", 72, 62)
}

# ===================================================================================================
# G - cell-level correlation of the pair in the key cell type, as in the paper
# ===================================================================================================
if (stage %in% c("G", "panels", "all")) {
  md <- cc$md[cell_type == cc$key]
  cr <- cc$corr
  lab <- sprintf("R = %.2f, p %s", cr$cell_R,
                 if (cr$cell_p < 0.001) "< 0.001" else sprintf("= %.3f", cr$cell_p))
  pG <- ggplot(md, aes(.data[[cc$pair[2]]], .data[[cc$pair[1]]])) +
    geom_point(size = 0.25, colour = "black", alpha = 0.55, stroke = 0) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "red", fill = "grey70",
                alpha = 0.4, linewidth = 0.5) +
    annotate("label", x = -Inf, y = Inf, label = lab, hjust = -0.06, vjust = 1.25,
             size = 2.1, linewidth = 0.2, label.padding = unit(0.8, "mm")) +
    labs(x = bquote(italic(.(cc$pair[2]))), y = bquote(italic(.(cc$pair[1])))) +
    theme_baseR()
  keep(pG, "Fig3G_correlation", 72, 62)
}

# ===================================================================================================
# assemble - the paper's layout: A | B over C on row 1, then D E F G across row 2
# ===================================================================================================
if (stage %in% c("assemble", "all")) {
  need <- c("Fig3A_umap", "Fig3B_markers", "Fig3C_composition", "Fig3D_violin1",
            "Fig3E_violin2", "Fig3F_density", "Fig3G_correlation")
  ps <- lapply(need, function(n) {
    f <- file.path(panel_dir, paste0(n, ".rds"))
    if (!file.exists(f)) stop("missing panel ", n, "; run stage panels first")
    readRDS(f)$plot
  })
  design <- "
AAAAAABBBBBBBB
AAAAAABBBBBBBB
AAAAAACCCCCCCC
AAAAAACCCCCCCC
DDDEEEFFFFGGGG
DDDEEEFFFFGGGG
"
  comp <- ps[[1]] + ps[[2]] + ps[[3]] + ps[[4]] + ps[[5]] + ps[[6]] + ps[[7]] +
    plot_layout(design = design) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 13, face = "bold"))
  out <- file.path(fig_dir, paste0("Fig3_paper_style_", tissue))
  ggsave(paste0(out, ".pdf"), comp, width = 260, height = 190, units = "mm", device = cairo_pdf)
  tif <- file.path("results/figure_exports", paste0("Fig3_paper_style_", tissue, ".tiff"))
  unlink(tif)
  ggsave(tif, comp, width = 260, height = 190, units = "mm", dpi = 400, bg = "white",
         compression = "lzw")
  message("wrote ", out, ".pdf and ", tif)
}

# ===================================================================================================
if (stage %in% c("legend", "panels", "all")) {
  cr <- cc$corr
  st <- lapply(setNames(cc$pair, cc$pair), pair_stats)
  nm <- st[[1]]$per_mouse[, .(n = uniqueN(mouse_id)), by = condition]
  sig_word <- if (SIG_UNIT == "mouse")
    sprintf("per-mouse t-test (n = %d vs %d mice)", nm[condition == "Control", n],
            nm[condition == "STZ", n]) else "cell-level Wilcoxon test"

  l <- c(
    sprintf("# Figure 3, reference-paper style (%s) - legend", tissue),
    "",
    "Drawn to match Xu et al., *Phytomedicine* 154 (2026) 158050, Fig. 3, with this project's data.",
    sprintf("The Nature/Cell-contract version of the same figure is in `figures/Fig3_composite_%s/`.",
            tissue),
    "**Ship one or the other, not both.**",
    "",
    sprintf("**Fig. 3. Results of scRNA-seq analysis.** (A) UMAP projection illustrating the single-cell atlas of mouse %s tissue. (B) Bubble plot showing marker genes for each cell type. (C) Bar plot displaying the proportion and number of cells from each sample source. (D-E) Analysis of %s and %s expression levels in %s across experimental groups. %s; ns: P > 0.05, *: P < 0.05, **: P < 0.01, ***: P < 0.001, ****: P < 0.0001. (F) UMAP-based heatmap showing co-localization of %s and %s across the atlas. (G) Correlation between %s and %s expression levels in %s.",
            tissue, cc$pair[1], cc$pair[2], tolower(nice_ct(cc$key)),
            paste0(toupper(substring(sig_word, 1, 1)), substring(sig_word, 2)),
            cc$pair[1], cc$pair[2], cc$pair[1], cc$pair[2], tolower(nice_ct(cc$key))),
    "",
    "## Where this data forces a deviation from the reference",
    "",
    sprintf("- **Two groups, not three.** The reference contrasts Ctrl, CIA and Treat, so its D and E carry three brackets. This atlas is Control vs STZ, so each panel carries one."),
    sprintf("- **The bracket is the %s, not the cell-level Wilcoxon.** Cells within a mouse are not independent replicates, so a cell-level P is anticonservative. Both are tabulated below; the reference reported the cell-level test.",
            sig_word),
    sprintf("- **(F) and (G) do not agree, and that is the result.** The pair is co-localized enough to draw, but in %s the cell-level correlation is %.2f and the background-matched empirical P is %.3f, so the pair-level criterion of the reference framework is **not** met in this tissue.",
            tolower(nice_ct(cc$key)), cr$cell_R, cr$empirical_p),
    sprintf("- **(A) labels every cluster on the plot** as the reference does; the reference's atlas is bone marrow with 11 types, this one is %s with %d.",
            tissue, nrow(prop)),
    sprintf("- **D, E and G are computed inside %s only,** the cell type the pair was selected in. The atlas-wide test is a different and weaker number (%s over all %d cells gives per-mouse P = %s), and is not what these panels show.",
            nice_ct(cc$key), cc$pair[1], cc$n_cells,
            format(signif(cc$gene_tests[gene == cc$pair[1], mouse_p], 2))),
    "",
    "## Statistics to quote",
    "",
    "| Quantity | Value |",
    "|---|---|",
    sprintf("| Cells in the atlas | %d |", cc$n_cells),
    sprintf("| Cell types | %d |", nrow(prop)),
    sprintf("| Most abundant type | %s (%.2f%%) |", nice_ct(prop$cell_type[1]), prop$pct[1]),
    sprintf("| Key cell type for the pair | %s (%d cells, %d mice) |", nice_ct(cc$key),
            cr$n_cells, cr$n_mice),
    sprintf("| %s in %s, Control vs STZ | cell-level P = %s; %s P = %s |", cc$pair[1],
            nice_ct(cc$key), format(signif(st[[1]]$cell_p, 2), scientific = TRUE), sig_word,
            format(signif(st[[1]]$mouse_p, 2))),
    sprintf("| %s in %s, Control vs STZ | cell-level P = %s; %s P = %s |", cc$pair[2],
            nice_ct(cc$key), format(signif(st[[2]]$cell_p, 2), scientific = TRUE), sig_word,
            format(signif(st[[2]]$mouse_p, 2))),
    sprintf("| Detection rate in %s | %s %.1f%% / %.1f%%, %s %.1f%% / %.1f%% (Control / STZ) |",
            nice_ct(cc$key), cc$pair[1],
            st[[1]]$det[condition == "Control", pct], st[[1]]$det[condition == "STZ", pct],
            cc$pair[2],
            st[[2]]$det[condition == "Control", pct], st[[2]]$det[condition == "STZ", pct]),
    sprintf("| Cells co-expressing the pair | %.1f%% |", cr$pct_coexpr),
    sprintf("| Pair correlation in %s | R = %.2f, P = %s, empirical P = %.3f |",
            nice_ct(cc$key), cr$cell_R, format(signif(cr$cell_p, 2)), cr$empirical_p)
  )
  writeLines(l, file.path(fig_dir, "Fig3_paper_style_legend.md"))
  message("wrote ", file.path(fig_dir, "Fig3_paper_style_legend.md"))
}
