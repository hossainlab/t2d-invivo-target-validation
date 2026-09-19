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
# This is the shipped Figure 1. It follows the reference paper's visual language, not the
# figure_style.R Nature/Cell contract that scripts 27-32 use; that set is no longer built.
#
# Deviations from the reference forced by this dataset, all of them stated in Fig1_legend.md:
#   - B uses nominal P, not adjusted P (decision M9). 84 genes do reach FDR < 0.05, but only 31 of them
#     also clear |log2FC| > 0.5, and the paper's strict downstream rules then leave one intersecting
#     gene (SREBF2, decision R17). The paper had 1,648 DEGs at adjusted P < 0.05.
#   - C draws its reference line at the 0.80 fit actually used here; the paper's line is at 0.90.
#   - E has 14 modules plus the unassigned grey bin, not seven, and the strongest is red (r = 0.65) not
#     blue (r = 0.83). Grey is recomputed here for display only; script 03 keeps it out of the rules.
#   - F shows the red module and the M10 hub thresholds (|kME| > 0.7, |GS| > 0.2), not 0.7 / 0.7. Its
#     threshold lines are black dashed because the module colour is red.
#
# Outputs: figures/Fig1/Fig1<A-G>_*.pdf, Fig1.pdf, Fig1_legend.md
#          results/figure_exports/Fig1.tiff

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(WGCNA)
  library(ggvenn); library(ggplotify)
})
options(stringsAsFactors = FALSE)

args  <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) >= 1) args[1] else "all"
set.seed(20260918)

fig_dir   <- "figures/Fig1"
panel_dir <- file.path("results/figure_exports", "fig1_panels_ref")
bulk      <- "results/bulk"
for (d in c(fig_dir, panel_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LFC <- 0.5; PCUT <- 0.05
FIT_LINE <- 0.80          # the paper's line is 0.90; this project selects power at 0.80
MM_CUT <- 0.7; GS_CUT <- 0.2   # M10 hub rule; the paper used 0.7 / 0.7

# the paper's palette: dataset colours in A, Down/Not/Up in B, Venn fills in G
col_ds  <- c("#F8766D", "#00BA38", "#619CFF")
col_sig <- c(Down = "#3B7DD8", Not = "grey75", Up = "#E8362B")

# The reference draws C, D and E with base-R WGCNA calls. C and E are redrawn here in ggplot to the same
# visual specification (black box, ticks out, no grid, bold centred title) because base-R panels captured
# by as.ggplot() do not rescale with the patchwork cell and lose their axis titles in the composite.
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
  # drawn in ggplot, not base R. The base-R version reproduced the paper's look at one device size only:
  # as.ggplot() captures at the size patchwork gives the cell, and the margins set in par() then clipped
  # both x-axis titles and the left y-axis title out of the composite. ggplot keeps the same red-numeral
  # look and sizes its own margins, so the panel survives assembly.
  pC1 <- ggplot(fi, aes(Power, R2, label = Power)) +
    geom_hline(yintercept = FIT_LINE, colour = "red", linewidth = 0.4) +
    geom_text(colour = "red", size = 2.1) +
    labs(title = "Scale independence", x = "Soft Threshold (power)",
         y = expression("Scale Free Topology Model Fit, signed R"^2)) +
    theme_baseR()
  pC2 <- ggplot(fi, aes(Power, mean.k., label = Power)) +
    geom_text(colour = "red", size = 2.1) +
    labs(title = "Mean connectivity", x = "Soft Threshold (power)", y = "Mean Connectivity") +
    theme_baseR()
  # wrap_elements so the two sub-plots take a single patchwork tag ("C") instead of one each
  pC <- wrap_elements(full = pC1 + pC2)
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
  # The reference heatmap ends on MEgrey. Script 03 drops grey from WGCNA_module_trait.csv because grey
  # is the unassigned-gene bin, not a module, so it must not enter the key-module rule. It is recomputed
  # here for display only, by exactly the rule script 03 uses for every other row: Pearson r of the
  # module eigengene against the T2D indicator, p from corPvalueStudent at n = samples.
  grey <- local({
    nt <- readRDS(file.path(bulk, "WGCNA_net.rds"))
    if (!"grey" %in% nt$moduleColors) return(NULL)
    dat <- t(readRDS(file.path(bulk, "expr_merged_combat.rds")))[, nt$genes, drop = FALSE]
    t2d <- fread(file.path(bulk, "meta_bulk.csv"))[match(rownames(dat), GSM), T2D]
    me  <- moduleEigengenes(dat, nt$moduleColors)$eigengenes[["MEgrey"]]
    r   <- cor(me, t2d, use = "p")
    data.table(module = "grey", size = sum(nt$moduleColors == "grey"),
               r_T2D = as.numeric(r), p_T2D = as.numeric(corPvalueStudent(r, nrow(dat))))
  })
  mt <- rbind(mt, grey, fill = TRUE)   # grey last, as in the reference, not in the r order

  # Control is the complement indicator of T2D: r is the exact negative, p identical
  d <- rbindlist(list(
    data.table(module = mt$module, x = 1, r = -mt$r_T2D, p = mt$p_T2D),
    data.table(module = mt$module, x = 2, r =  mt$r_T2D, p = mt$p_T2D)))
  d[, y := factor(paste0("ME", module), levels = rev(paste0("ME", mt$module)))]
  d[, lab := paste0(signif(r, 2), "\n(", signif(p, 1), ")")]
  sw <- data.table(y = factor(paste0("ME", mt$module), levels = levels(d$y)), module = mt$module)

  pE <- ggplot(d, aes(x, y)) +
    geom_tile(aes(fill = r), width = 1, height = 1) +
    # fixed fill, passed as a parameter so the module swatches bypass the r colour scale
    geom_tile(data = sw, aes(x = 0.34, y = y), fill = sw$module, width = 0.22, height = 0.9,
              colour = "black", linewidth = 0.15, inherit.aes = FALSE) +
    geom_text(aes(label = lab), size = 1.75, lineheight = 0.95) +
    scale_fill_gradientn(colours = blueWhiteRed(50), limits = c(-1, 1), name = NULL,
                         guide = guide_colourbar(barwidth = unit(2.6, "mm"),
                                                 barheight = unit(34, "mm"), ticks = FALSE)) +
    scale_x_continuous(breaks = c(1, 2), labels = c("Control", "T2D"),
                       limits = c(0.1, 2.5), expand = c(0, 0)) +
    labs(title = "Module-trait relationships", x = NULL, y = NULL) +
    theme_baseR() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 6.5),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 6),
          legend.margin = margin(0, 0, 0, 1))
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
  # the reference draws the points in the module's own colour (blue module, blue points); keep that,
  # and only fall back to grey for module names too pale to read as open circles on white
  pt_col <- if (best %in% c("white", "lightyellow", "lightcyan", "ivory")) "grey25" else best
  # the reference draws the two threshold lines in red over a blue module. The strongest module here is
  # red, so red lines would disappear into the points; they go black and dashed for a red-family module.
  reddish <- best %in% c("red", "darkred", "salmon", "orangered", "magenta", "pink")
  ln_col  <- if (reddish) "black" else "red"
  ln_type <- if (reddish) "dashed" else "solid"
  d <- data.table(mm = mm, gs = gs)
  pF <- ggplot(d, aes(mm, gs)) +
    geom_point(shape = 1, size = 1, stroke = 0.35, colour = pt_col) +
    geom_vline(xintercept = MM_CUT, colour = ln_col, linetype = ln_type, linewidth = 0.4) +
    geom_hline(yintercept = GS_CUT, colour = ln_col, linetype = ln_type, linewidth = 0.4) +
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
  # set names carry the circle colour, as in the reference. text_size is held down because the three-way
  # region here holds 8 genes: at the reference's size its count and percentage ran into the 12 and 26
  # labels in the two lower lens regions.
  pG <- ggvenn(sets, fill_color = c("#E8837D", "#5BB3E8", "#F4A64B"), stroke_size = 0.4,
               show_percentage = TRUE, digits = 1, set_name_size = 3.2, text_size = 1.85,
               set_name_color = c("#C0392B", "#2471A3", "#CA7A11"),
               fill_alpha = 0.55) +
    coord_fixed(clip = "off") +
    theme(plot.margin = margin(2, 4, 2, 4))
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
EEEEFFFFFGGGGG
EEEEFFFFFGGGGG
EEEEFFFFFGGGGG
EEEEFFFFFGGGGG
"
  comp <- ps[[1]] + ps[[2]] + ps[[3]] + ps[[4]] + ps[[5]] + ps[[6]] + ps[[7]] +
    plot_layout(design = design) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 13, face = "bold"))
  out <- file.path(fig_dir, "Fig1")
  ggsave(paste0(out, ".pdf"), comp, width = 260, height = 250, units = "mm", device = cairo_pdf)
  tif <- file.path("results/figure_exports", "Fig1.tiff")
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
  nfdr <- deg[adj.P.Val < 0.05, .N]
  nfdr_lfc <- deg[adj.P.Val < 0.05 & abs(logFC) > LFC, .N]
  pwr  <- readRDS(file.path(bulk, "WGCNA_net.rds"))$power

  l <- c(
    "# Figure 1, reference-paper style - legend",
    "",
    "Drawn to match Xu et al., *Phytomedicine* 154 (2026) 158050, Fig. 1, with this project's data.",
    "",
    sprintf("**Fig. 1. Screening of DEGs related to T2D in human liver.** (A) PCA plot before and after batch effect correction. (B) Volcano plot of DEGs. Blue indicates significantly downregulated genes, red indicates significantly upregulated genes, and grey indicates non-significant genes. (C) Scale independence and mean connectivity in the weighted gene co-expression network. At a soft threshold of %d, the network exhibits high scale independence and mean connectivity close to zero. (D) Different co-expression modules within the weighted gene co-expression network. (E) Heatmap showing the correlation of %d modules and the unassigned grey bin with T2D case and control groups. (F) Correlation between MM and GS for all genes in the %s module. (G) Venn diagram illustrating the intersection among DEGs, key-module genes, and hyperglycaemia-related genes.",
            pwr, nrow(mt), best),
    "",
    "## Where this data forces a deviation from the reference",
    "",
    sprintf("- **(B) uses nominal P, not adjusted P** (decision M9). The paper reported 1,648 DEGs at adjusted P < 0.05. Here %d genes reach FDR < 0.05 and only %d of those also clear |log2FC| > %.1f; carrying that set through the paper's strict module and gene-set rules leaves a single intersecting gene (SREBF2, decision R17), which no downstream analysis can rest on. The panel therefore shows nominal P < %.2f and |log2FC| > %.1f: %d up, %d down, labelled exploratory. The adjusted-P column is reported in `DEG_T2D_vs_Control_all.csv` throughout.",
            nfdr, nfdr_lfc, LFC, PCUT, LFC, nup, ndn),
    sprintf("- **(C) reference line at %.2f,** not the paper's 0.90; %.2f is the fit this project selects power at. Power %d happens to match the paper.",
            FIT_LINE, FIT_LINE, pwr),
    sprintf("- **(E) has %d modules plus the unassigned grey bin, not seven,** and the strongest is **%s** (r = %s) rather than the paper's blue (r = 0.83). Control is the complement of T2D, so its column is the exact negative with an identical P, exactly as in the reference. Grey is drawn last, as in the reference, and is display-only: it is the unassigned-gene bin, so script 03 excludes it from the key-module rule.",
            nrow(mt), best, format(round(mt[module == best, r_T2D], 2), nsmall = 2)),
    sprintf("- **(F) shows the %s module** and this project's hub thresholds (|kME| > %.1f, |GS| > %.1f, decision M10), not the paper's 0.7 / 0.7. Points take the module colour as in the reference; because that colour is %s here, the two threshold lines are drawn black and dashed rather than the reference's red, which would be invisible against the points.",
            best, MM_CUT, GS_CUT, best),
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
    sprintf("| Genes at FDR < 0.05 | %d |", nfdr),
    sprintf("| Genes at FDR < 0.05 and \\|log2FC\\| > %.1f | %d |", LFC, nfdr_lfc),
    sprintf("| Soft-thresholding power | %d |", pwr),
    sprintf("| Modules (excluding grey) | %d |", nrow(mt)),
    sprintf("| Key modules | %s |", paste(km, collapse = ", ")),
    sprintf("| Strongest module | %s, n = %d genes, r = %s, P = %s |", best,
            mt[module == best, size],
            format(round(mt[module == best, r_T2D], 2), nsmall = 2),
            format(signif(mt[module == best, p_T2D], 2), scientific = TRUE)),
    local({
      g  <- fread(file.path(bulk, "WGCNA_gene_module_GS_MM.csv"))[module == best]
      ct <- suppressWarnings(cor.test(abs(g[[paste0("MM_", best)]]), abs(g$GS_T2D)))
      sprintf("| MM vs GS in the %s module (panel F) | cor = %.2f, P = %s |", best, ct$estimate,
              format(signif(ct$p.value, 2), scientific = TRUE))
    }),
    sprintf("| Three-way intersection | %s genes |", isum[item == "triple_intersection", value])
  )
  writeLines(l, file.path(fig_dir, "Fig1_legend.md"))
  message("wrote ", file.path(fig_dir, "Fig1_legend.md"))
}
