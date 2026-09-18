# 33_fig1_paper_style.R - Figure 1 drawn in the visual style of the reference paper.
#
# Usage: Rscript scripts/33_fig1_paper_style.R <stage>
#   stage : A | B | C | D | E | F | G | panels | assemble | legend | all
#
# Reference: Xu et al., Phytomedicine 154 (2026) 158050, Fig. 1. This reproduces that figure's own visual
# language with this project's data, panel for panel:
#   A  PCA before / after batch correction, coloured and shaped by dataset, with confidence ellipses
#   B  volcano, Down blue / Not grey / Up red, legend titled "Sig", no gene labels
#   C  WGCNA base-R scale independence and mean connectivity, powers drawn as red numerals
#   D  gene dendrogram with a single "Merged dynamic" colour band
#   E  WGCNA labeledHeatmap, module colour squares, Control and T2D columns, r over (p)
#   F  MM vs GS for the single strongest module, open circles, threshold lines
#   G  ggvenn with counts and percentages, set names in the circle colours
#
# This is deliberately NOT the figure_style.R contract used by scripts 27-32. It is the reference
# paper's style, kept as a separate deliverable in figures/Fig1_paper_style/ so the Nature/Cell version
# in figures/Fig1_bulk_DEG_WGCNA/ is untouched. Pick one for submission; do not ship both.
#
# Deviations from the reference forced by this dataset, all of them stated in Fig1_paper_style_legend.md:
#   - B uses nominal P, not adjusted P. No gene reaches FDR < 0.05 here (decision M9), so an adjusted-P
#     volcano would show zero significant genes. The paper had 1,648 DEGs at adjusted P < 0.05.
#   - C draws its reference line at the 0.80 fit actually used here; the paper's line is at 0.90.
#   - E has 14 modules plus grey, not seven, and the strongest is red (r = 0.65) not blue (r = 0.83).
#   - F shows the red module and the M10 hub thresholds (|kME| > 0.7, |GS| > 0.2), not 0.7 / 0.7.
#
# Outputs: figures/Fig1_paper_style/Fig1<A-G>_*.pdf, Fig1_paper_style.pdf, Fig1_paper_style_legend.md
#          results/figure_exports/Fig1_paper_style.tiff

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(WGCNA)
  library(ggvenn); library(ggplotify)
})
options(stringsAsFactors = FALSE)

args  <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) >= 1) args[1] else "all"
set.seed(20260918)

fig_dir   <- "figures/Fig1_paper_style"
panel_dir <- file.path("results/figure_exports", "fig1paper_panels")
bulk      <- "results/bulk"
for (d in c(fig_dir, panel_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LFC <- 0.5; PCUT <- 0.05
FIT_LINE <- 0.80          # the paper's line is 0.90; this project selects power at 0.80
MM_CUT <- 0.7; GS_CUT <- 0.2   # M10 hub rule; the paper used 0.7 / 0.7

# the paper's palette: dataset colours in A, Down/Not/Up in B, Venn fills in G
col_ds  <- c("#F8766D", "#00BA38", "#619CFF")
col_sig <- c(Down = "#3B7DD8", Not = "grey75", Up = "#E8362B")

keep <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, units = "mm",
         device = cairo_pdf)
  saveRDS(list(plot = p, width = w, height = h), file.path(panel_dir, paste0(name, ".rds")))
  message("wrote ", name, " (", w, " x ", h, " mm)")
  invisible(p)
}
key_modules <- function()
  strsplit(fread(file.path(bulk, "intersection_summary.csv"))[item == "key_modules", value], ";")[[1]]

# ===================================================================================================
# A - PCA before / after batch correction, with confidence ellipses, as in the paper
# ===================================================================================================
if (stage %in% c("A", "panels", "all")) {
  meta <- fread(file.path(bulk, "meta_bulk.csv"))
  # drawn as one faceted object, not two plots side by side: patchwork assigns one tag per plot, and a
  # two-plot panel A would consume the tags A and B and push every later panel along by one
  pcs <- rbindlist(lapply(c(`Before batch correction` = "expr_merged_raw.rds",
                            `After batch correction`  = "expr_merged_combat.rds"), function(f) {
    e  <- readRDS(file.path(bulk, f))
    m  <- meta[match(colnames(e), GSM)]
    pv <- prcomp(t(e), scale. = TRUE)
    data.table(PC1 = pv$x[, 1], PC2 = pv$x[, 2], Type = m$dataset)
  }), idcol = "stg")
  pcs[, stg := factor(stg, levels = c("Before batch correction", "After batch correction"))]

  pA <- ggplot(pcs, aes(PC1, PC2, colour = Type, shape = Type, fill = Type)) +
    stat_ellipse(aes(group = Type), type = "norm", level = 0.95, linewidth = 0.4,
                 geom = "polygon", alpha = 0.12, show.legend = FALSE) +
    geom_point(size = 1.5) +
    facet_wrap(~ stg, scales = "free") +
    scale_colour_manual(values = col_ds[1:2]) +
    scale_fill_manual(values = col_ds[1:2]) +
    scale_shape_manual(values = c(16, 17)) +
    labs(x = "PC1", y = "PC2") +
    theme_classic(base_size = 8) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 9),
          legend.title = element_text(size = 8), legend.text = element_text(size = 7),
          legend.key.size = unit(3.5, "mm"),
          panel.spacing.x = unit(5, "mm"))
  keep(pA, "Fig1A_PCA", 170, 62)
}

# ===================================================================================================
# B - volcano, paper colour scheme and legend, no gene labels
# ===================================================================================================
if (stage %in% c("B", "panels", "all")) {
  tab <- fread(file.path(bulk, "DEG_T2D_vs_Control_all.csv"))
  tab[, Sig := fifelse(P.Value < PCUT & logFC >  LFC, "Up",
              fifelse(P.Value < PCUT & logFC < -LFC, "Down", "Not"))]
  tab[, Sig := factor(Sig, levels = c("Down", "Not", "Up"))]
  pB <- ggplot(tab[order(Sig)], aes(logFC, -log10(P.Value), colour = Sig)) +
    geom_point(size = 0.7, alpha = 0.85) +
    scale_colour_manual(values = col_sig, name = "Sig") +
    labs(x = "logFC", y = expression(-log[10]*"(P.Value)")) +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(),
          legend.title = element_text(size = 8), legend.text = element_text(size = 7),
          legend.key.size = unit(3.5, "mm"))
  keep(pB, "Fig1B_volcano", 95, 72)
}

# ===================================================================================================
# C - WGCNA soft threshold, base-R with red numerals, as in the paper
# ===================================================================================================
if (stage %in% c("C", "panels", "all")) {
  fi <- fread(file.path(bulk, "WGCNA_soft_threshold.csv"))
  pC <- as.ggplot(function() {
    par(mfrow = c(1, 2), oma = c(0, 0, 1.8, 0), mar = c(4.0, 4.8, 2.6, 0.8),
        mgp = c(2.9, 0.6, 0), cex = 0.55)
    plot(fi$Power, fi$R2, type = "n", xlab = "Soft Threshold (power)",
         ylab = "Scale Free Topology Fit, signed R^2", main = "Scale independence")
    text(fi$Power, fi$R2, labels = fi$Power, col = "red", cex = 0.85)
    abline(h = FIT_LINE, col = "red")
    plot(fi$Power, fi$mean.k., type = "n", xlab = "Soft Threshold (power)",
         ylab = "Mean Connectivity", main = "Mean connectivity")
    text(fi$Power, fi$mean.k., labels = fi$Power, col = "red", cex = 0.85)
  })
  keep(pC, "Fig1C_soft_threshold", 150, 72)
}

# ===================================================================================================
# D - gene dendrogram with a single "Merged dynamic" band, as in the paper
# ===================================================================================================
if (stage %in% c("D", "panels", "all")) {
  nt <- readRDS(file.path(bulk, "WGCNA_net.rds"))
  pD <- as.ggplot(function() {
    par(oma = c(0, 0, 1.8, 0))
    plotDendroAndColors(nt$net$dendrograms[[1]], nt$moduleColors[nt$net$blockGenes[[1]]],
                        groupLabels = "Merged dynamic", dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = "Gene dendrogram and module colors", xlab = "", sub = "",
                        cex.colorLabels = 0.55, cex.axis = 0.55, cex.lab = 0.65, cex.main = 0.75,
                        marAll = c(0.6, 5.6, 2.6, 0.4))
  })
  keep(pD, "Fig1D_dendrogram", 120, 72)
}

# ===================================================================================================
# E - labeledHeatmap, Control and T2D columns, module colour squares, as in the paper
# ===================================================================================================
if (stage %in% c("E", "panels", "all")) {
  mt <- fread(file.path(bulk, "WGCNA_module_trait.csv"))
  setorder(mt, -r_T2D)
  # Control is the complement indicator of T2D: r is the exact negative, p identical
  mat <- cbind(Control = -mt$r_T2D, T2D = mt$r_T2D)
  pmat <- cbind(mt$p_T2D, mt$p_T2D)
  rownames(mat) <- paste0("ME", mt$module)
  txt <- paste0(signif(mat, 2), "\n(", signif(pmat, 1), ")")
  dim(txt) <- dim(mat)
  pE <- as.ggplot(function() {
    par(oma = c(0, 0, 3.2, 0), mar = c(3.2, 11.5, 3.4, 3.4), cex = 0.55)
    labeledHeatmap(Matrix = mat, xLabels = c("Control", "T2D"),
                   yLabels = rownames(mat), ySymbols = rownames(mat),
                   colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = txt,
                   setStdMargins = FALSE, cex.text = 0.5, cex.lab = 0.62, zlim = c(-1, 1),
                   main = "Module-trait relationships")
  })
  keep(pE, "Fig1E_module_trait", 90, 118)
}

# ===================================================================================================
# F - MM vs GS in the single strongest module, as in the paper's blue-module panel
# ===================================================================================================
if (stage %in% c("F", "panels", "all")) {
  gt <- fread(file.path(bulk, "WGCNA_gene_module_GS_MM.csv"))
  mt <- fread(file.path(bulk, "WGCNA_module_trait.csv"))
  km <- key_modules()
  best <- mt[module %in% km][which.max(abs(r_T2D)), module]
  g <- gt[module == best]
  mm <- abs(g[[paste0("MM_", best)]]); gs <- abs(g$GS_T2D)
  ct <- suppressWarnings(cor.test(mm, gs))
  pt_col <- if (best %in% c("red", "yellow", "white", "lightyellow")) "grey25" else best
  d <- data.table(mm = mm, gs = gs)
  pF <- ggplot(d, aes(mm, gs)) +
    geom_point(shape = 1, size = 1, stroke = 0.35, colour = pt_col) +
    geom_vline(xintercept = MM_CUT, colour = "red", linewidth = 0.4) +
    geom_hline(yintercept = GS_CUT, colour = "red", linewidth = 0.4) +
    labs(title = sprintf("Module membership vs. gene significance\ncor=%.2f, p=%s",
                         ct$estimate, format(signif(ct$p.value, 2), scientific = TRUE)),
         x = sprintf("Module Membership in %s module", best),
         y = "Gene significance for T2D") +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 8.5, face = "bold"))
  keep(pF, "Fig1F_MM_GS", 90, 82)
}

# ===================================================================================================
# G - ggvenn with counts and percentages, set names in the circle colours
# ===================================================================================================
if (stage %in% c("G", "panels", "all")) {
  deg <- fread(file.path(bulk, "DEG_T2D_vs_Control_all.csv"))
  gt  <- fread(file.path(bulk, "WGCNA_gene_module_GS_MM.csv"))
  gs  <- fread(file.path(bulk, "geneset_hyperglycemia.csv"))
  sets <- list(DEGs = deg[P.Value < PCUT & abs(logFC) > LFC, gene],
               WGCNA = gt[module %in% key_modules(), gene],
               HG = intersect(gs$gene, deg$gene))
  pG <- ggvenn(sets, fill_color = c("#E8837D", "#5BB3E8", "#F4A64B"), stroke_size = 0.4,
               show_percentage = TRUE, set_name_size = 3.2, text_size = 2.6,
               fill_alpha = 0.55) +
    coord_fixed(clip = "off") +
    theme(plot.margin = margin(6, 6, 6, 6))
  keep(pG, "Fig1G_venn", 95, 95)
}

# ===================================================================================================
# assemble - the paper's layout: A|B on row 1, C|D on row 2, E|F|G on row 3
# ===================================================================================================
if (stage %in% c("assemble", "all")) {
  need <- c("Fig1A_PCA", "Fig1B_volcano", "Fig1C_soft_threshold", "Fig1D_dendrogram",
            "Fig1E_module_trait", "Fig1F_MM_GS", "Fig1G_venn")
  ps <- lapply(need, function(n) {
    f <- file.path(panel_dir, paste0(n, ".rds"))
    if (!file.exists(f)) stop("missing panel ", n, "; run stage panels first")
    readRDS(f)$plot
  })
  design <- "
AAAAAAAAABBBBB
AAAAAAAAABBBBB
CCCCCCCCDDDDDD
CCCCCCCCDDDDDD
CCCCCCCCDDDDDD
EEEEEFFFFFGGGG
EEEEEFFFFFGGGG
EEEEEFFFFFGGGG
EEEEEFFFFFGGGG
"
  comp <- ps[[1]] + ps[[2]] + ps[[3]] + ps[[4]] + ps[[5]] + ps[[6]] + ps[[7]] +
    plot_layout(design = design) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 13, face = "bold"))
  out <- file.path(fig_dir, "Fig1_paper_style")
  ggsave(paste0(out, ".pdf"), comp, width = 260, height = 250, units = "mm", device = cairo_pdf)
  tif <- file.path("results/figure_exports", "Fig1_paper_style.tiff")
  unlink(tif)
  ggsave(tif, comp, width = 260, height = 250, units = "mm", dpi = 400, bg = "white",
         compression = "lzw")
  message("wrote ", out, ".pdf and ", tif)
}

# ===================================================================================================
if (stage %in% c("legend", "panels", "all")) {
  deg  <- fread(file.path(bulk, "DEG_T2D_vs_Control_all.csv"))
  mt   <- fread(file.path(bulk, "WGCNA_module_trait.csv"))
  isum <- fread(file.path(bulk, "intersection_summary.csv"))
  km   <- key_modules()
  best <- mt[module %in% km][which.max(abs(r_T2D)), module]
  nup  <- deg[P.Value < PCUT & logFC >  LFC, .N]
  ndn  <- deg[P.Value < PCUT & logFC < -LFC, .N]
  pwr  <- readRDS(file.path(bulk, "WGCNA_net.rds"))$power

  l <- c(
    "# Figure 1, reference-paper style - legend",
    "",
    "Drawn to match Xu et al., *Phytomedicine* 154 (2026) 158050, Fig. 1, with this project's data.",
    "The Nature/Cell-contract version of the same figure is in `figures/Fig1_bulk_DEG_WGCNA/`.",
    "**Ship one or the other, not both.**",
    "",
    sprintf("**Fig. 1. Screening of DEGs related to T2D in human liver.** (A) PCA plot before and after batch effect correction. (B) Volcano plot of DEGs. Blue indicates significantly downregulated genes, red indicates significantly upregulated genes, and grey indicates non-significant genes. (C) Scale independence and mean connectivity in the weighted gene co-expression network. At a soft threshold of %d, the network exhibits high scale independence and mean connectivity close to zero. (D) Different co-expression modules within the weighted gene co-expression network. (E) Heatmap showing the correlation of %d modules with T2D case and control groups. (F) Correlation between MM and GS for all genes in the %s module. (G) Venn diagram illustrating the intersection among DEGs, key-module genes, and hyperglycaemia-related genes.",
            pwr, nrow(mt), best),
    "",
    "## Where this data forces a deviation from the reference",
    "",
    sprintf("- **(B) uses nominal P, not adjusted P.** The paper reported 1,648 DEGs at adjusted P < 0.05. Here **no gene reaches FDR < 0.05** (decision M9), so an adjusted-P volcano would be empty. The panel shows nominal P < %.2f and |log2FC| > %.1f: %d up, %d down.",
            PCUT, LFC, nup, ndn),
    sprintf("- **(C) reference line at %.2f,** not the paper's 0.90; %.2f is the fit this project selects power at. Power %d happens to match the paper.",
            FIT_LINE, FIT_LINE, pwr),
    sprintf("- **(E) has %d modules plus grey, not seven,** and the strongest is **%s** (r = %s) rather than the paper's blue (r = 0.83). Control is the complement of T2D, so its column is the exact negative with an identical P, exactly as in the reference.",
            nrow(mt), best, format(round(mt[module == best, r_T2D], 2), nsmall = 2)),
    sprintf("- **(F) shows the %s module** and this project's hub thresholds (|kME| > %.1f, |GS| > %.1f, decision M10), not the paper's 0.7 / 0.7.",
            best, MM_CUT, GS_CUT),
    sprintf("- **(G)** the third set is hyperglycaemia-related genes (%s in the expression universe), standing in for the paper's neutrophil-migration set. The three-way intersection is **%s genes**.",
            isum[item == "geneset_in_universe", value], isum[item == "triple_intersection", value]),
    "- Key modules follow decision **M11b** (T2D P < 0.1 pooled and replicated in both cohorts), which is stricter than the paper's single-best-module rule.",
    "",
    "## Statistics to quote",
    "",
    "| Quantity | Value |",
    "|---|---|",
    "| Samples | 16 T2D, 11 lean control (GSE15653 + GSE64998) |",
    sprintf("| DEGs (nominal P < %.2f, \\|log2FC\\| > %.1f) | %d up, %d down |", PCUT, LFC, nup, ndn),
    sprintf("| Genes at FDR < 0.05 | %d |", deg[adj.P.Val < 0.05, .N]),
    sprintf("| Soft-thresholding power | %d |", pwr),
    sprintf("| Modules (excluding grey) | %d |", nrow(mt)),
    sprintf("| Key modules | %s |", paste(km, collapse = ", ")),
    sprintf("| Strongest module | %s, r = %s, P = %s |", best,
            format(round(mt[module == best, r_T2D], 2), nsmall = 2),
            format(signif(mt[module == best, p_T2D], 2), scientific = TRUE)),
    sprintf("| Three-way intersection | %s genes |", isum[item == "triple_intersection", value])
  )
  writeLines(l, file.path(fig_dir, "Fig1_paper_style_legend.md"))
  message("wrote ", file.path(fig_dir, "Fig1_paper_style_legend.md"))
}
