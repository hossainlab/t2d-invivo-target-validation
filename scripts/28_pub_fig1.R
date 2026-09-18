# 28_pub_fig1.R - Figure 1 (bulk DEG + WGCNA) redrawn to the publication contract in figure_style.R.
#
# Usage: Rscript scripts/28_pub_fig1.R <stage>
#   stage : A | B | C | D | E | F | G | panels | assemble | legend | all
#
# This script re-DRAWS; it does not re-ANALYSE. Every panel is rebuilt from the tables scripts 01-04
# already wrote, so the numbers are identical to the originals:
#   a  PCA before/after ComBat                <- expr_merged_raw.rds, expr_merged_combat.rds, meta_bulk.csv
#   b  volcano, T2D vs Control                <- DEG_T2D_vs_Control_all.csv
#   c  soft-threshold selection               <- WGCNA_soft_threshold.csv
#   d  gene dendrogram + modules              <- WGCNA_net.rds
#   e  module-trait correlations              <- WGCNA_module_trait.csv
#   f  module membership vs gene significance <- WGCNA_gene_module_GS_MM.csv
#   g  DEG / module / gene-set overlap        <- DEG + WGCNA + geneset_hyperglycemia.csv
#
# Panel e note: the Dataset column is computed on ComBat-corrected expression, so it is bounded near
# zero by construction and is not an independent batch check. It is marked with a dagger and the
# caveat is written into Fig1_legend.md rather than left for a reader to assume.
#
# Outputs: figures/Fig1_bulk_DEG_WGCNA/Fig1<a-g>_*.pdf, Fig1_composite.pdf, Fig1_legend.md
#          results/figure_exports/Fig1_composite.tiff, results/figure_exports/fig1_panels/*.rds

source("scripts/figure_style.R")
suppressPackageStartupMessages({ library(ggrepel) })

args  <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) >= 1) args[1] else "all"
set.seed(20260918)

fig_dir   <- "figures/Fig1_bulk_DEG_WGCNA"
panel_dir <- file.path(EXPORT_DIR, "fig1_panels")
bulk      <- "results/bulk"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

LFC <- 0.5; PCUT <- 0.05
keep <- function(p, n, w, h) save_fig(p, n, w, h, dir = fig_dir, keep_dir = panel_dir)
key_modules <- function()
  strsplit(fread(file.path(bulk, "intersection_summary.csv"))[item == "key_modules", value], ";")[[1]]

# ===================================================================================================
# a - PCA before and after ComBat
# ===================================================================================================
if (stage %in% c("A", "panels", "all")) {
  meta <- fread(file.path(bulk, "meta_bulk.csv"))
  pcs <- rbindlist(lapply(c(Uncorrected = "expr_merged_raw.rds", ComBat = "expr_merged_combat.rds"),
                          function(f) {
    e  <- readRDS(file.path(bulk, f))
    m  <- meta[match(colnames(e), GSM)]
    pv <- prcomp(t(e), scale. = TRUE)
    v  <- 100 * pv$sdev^2 / sum(pv$sdev^2)
    data.table(PC1 = pv$x[, 1], PC2 = pv$x[, 2], dataset = m$dataset, condition = m$condition,
               v1 = v[1], v2 = v[2])
  }), idcol = "stg")
  pcs[, stg := factor(stg, levels = c("Uncorrected", "ComBat"))]
  # variance explained goes inside the panel; as part of the strip label it overruns a 45 mm facet
  vlab <- pcs[, .(lab = sprintf("PC1 %.0f%% · PC2 %.0f%%", v1[1], v2[1])), by = stg]

  pA <- ggplot(pcs, aes(PC1, PC2, colour = dataset, shape = condition)) +
    geom_point(size = 1.1, stroke = 0.4) +
    geom_text(data = vlab, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
              hjust = -0.08, vjust = 1.5, size = (BASE - 2) * ppt, family = FONT,
              colour = "grey30") +
    facet_wrap(~ stg, scales = "free") +
    scale_colour_manual(values = c(GSE15653 = pal_qual[1], GSE64998 = pal_qual[2]), name = NULL) +
    scale_shape_manual(values = c(Control = 1, T2D = 16), name = NULL) +
    labs(x = "PC1", y = "PC2") +
    theme_pub() +
    theme(legend.position = "bottom", legend.box = "horizontal",
          legend.spacing.x = unit(2, "mm"), panel.spacing.x = unit(3, "mm"))
  keep(pA, "Fig1a_PCA_ComBat", 89, 58)
}

# ===================================================================================================
# b - volcano
# ===================================================================================================
if (stage %in% c("B", "panels", "all")) {
  tab <- fread(file.path(bulk, "DEG_T2D_vs_Control_all.csv"))
  tab[, call := fifelse(P.Value < PCUT & logFC > LFC, "Up",
                fifelse(P.Value < PCUT & logFC < -LFC, "Down", "NS"))]
  tab[, call := factor(call, levels = c("Up", "Down", "NS"))]
  nup <- tab[call == "Up", .N]; ndn <- tab[call == "Down", .N]
  lab <- rbind(tab[call == "Up"][order(P.Value)][1:10], tab[call == "Down"][order(P.Value)][1:10])

  pB <- ggplot(tab, aes(logFC, -log10(P.Value))) +
    raster_pts(geom_point(aes(colour = call), size = 0.35, alpha = 0.7, stroke = 0)) +
    geom_vline(xintercept = c(-LFC, LFC), linetype = 2, linewidth = 0.25, colour = "grey55") +
    geom_hline(yintercept = -log10(PCUT), linetype = 2, linewidth = 0.25, colour = "grey55") +
    geom_text_repel(data = lab, aes(label = gene), size = (BASE - 2) * ppt, family = FONT,
                    fontface = "italic", max.overlaps = 30, min.segment.length = 0.2,
                    segment.size = 0.2, segment.colour = "grey50", box.padding = 0.18) +
    scale_colour_manual(values = col_dir, breaks = c("Up", "Down"),
                        labels = c(sprintf("Up (%d)", nup), sprintf("Down (%d)", ndn)), name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 1.4, alpha = 1))) +
    labs(x = expression(log[2]~"fold change"), y = expression(-log[10]~italic(P))) +
    theme_pub() +
    theme(legend.position = "top", legend.justification = "right",
          legend.margin = margin(0, 0, 1, 0))
  keep(pB, "Fig1b_volcano", 89, 66)
}

# ===================================================================================================
# c - soft threshold
# ===================================================================================================
if (stage %in% c("C", "panels", "all")) {
  fi  <- fread(file.path(bulk, "WGCNA_soft_threshold.csv"))
  pwr <- readRDS(file.path(bulk, "WGCNA_net.rds"))$power
  d   <- melt(fi[, .(Power, `Scale-free fit` = R2, `Mean connectivity` = mean.k.)],
              id.vars = "Power", variable.name = "metric", value.name = "value")
  hl  <- data.table(metric = factor("Scale-free fit", levels = levels(d$metric)), y = 0.8)
  # computed into the data, not referenced from aes(): the saved plot is re-rendered by `assemble`
  # in a session where `pwr` does not exist
  d[, selected := Power == pwr]

  pC <- ggplot(d, aes(Power, value)) +
    geom_hline(data = hl, aes(yintercept = y), linetype = 2, linewidth = 0.25, colour = "grey55") +
    geom_line(linewidth = 0.3, colour = "grey60") +
    geom_point(aes(colour = selected), size = 1) +
    geom_text_repel(data = d[selected == TRUE], aes(label = Power), size = (BASE - 2) * ppt,
                    family = FONT, colour = col_dir[["Up"]], min.segment.length = 0,
                    segment.size = 0.2, nudge_y = 0.03, nudge_x = 1.5) +
    facet_wrap(~ metric, scales = "free_y") +
    scale_colour_manual(values = c(`TRUE` = col_dir[["Up"]], `FALSE` = "grey45"), guide = "none") +
    labs(x = "Soft-thresholding power", y = NULL) +
    theme_pub() +
    theme(panel.spacing.x = unit(3, "mm"))
  keep(pC, "Fig1c_soft_threshold", 89, 50)
}

# ===================================================================================================
# d - gene dendrogram with module colours (base graphics)
# ===================================================================================================
if (stage %in% c("D", "panels", "all")) {
  suppressPackageStartupMessages(library(WGCNA))
  nt <- readRDS(file.path(bulk, "WGCNA_net.rds"))
  open_pdf("Fig1d_dendrogram", 89, 52, fig_dir)
  plotDendroAndColors(nt$net$dendrograms[[1]], nt$moduleColors[nt$net$blockGenes[[1]]],
                      groupLabels = "Module", dendroLabels = FALSE, hang = 0.03,
                      addGuide = FALSE, main = "", ylab = "Co-expression distance",
                      cex.colorLabels = 0.85, cex.axis = 0.85, cex.lab = 1,
                      marAll = c(0.4, 4.4, 0.4, 0.3))
  dev.off()
  message("wrote Fig1d_dendrogram (89 x 52 mm)")
}

# ===================================================================================================
# e - module-trait correlation heatmap
# ===================================================================================================
if (stage %in% c("E", "panels", "all")) {
  mt <- fread(file.path(bulk, "WGCNA_module_trait.csv"))
  km <- key_modules()
  ds <- "Dataset†"   # dagger: bounded near zero by ComBat, not an independent batch test
  # the per-cohort columns are the ones the key-module rule actually uses (decision M11b), so they
  # belong in the panel; without them a reader cannot see why blue is excluded and greenyellow is not
  cols <- c("T2D (all)", "T2D\nGSE15653", "T2D\nGSE64998", "HbA1c", ds)
  rr <- melt(setnames(mt[, .(module, r_T2D, r_T2D_GSE15653, r_T2D_GSE64998,
                             r_HbA1c_GSE15653, r_Dataset)], c("module", cols)),
             id.vars = "module", variable.name = "trait", value.name = "r")
  pp <- melt(setnames(mt[, .(module, p_T2D, p_T2D_GSE15653, p_T2D_GSE64998,
                             p_HbA1c_GSE15653, p_Dataset)], c("module", cols)),
             id.vars = "module", variable.name = "trait", value.name = "p")
  d <- merge(rr, pp, by = c("module", "trait"))
  d[, trait := factor(trait, levels = cols)]
  d[, module := factor(module, levels = mt[order(r_T2D), module])]
  d[, lab := sprintf("%s\n%s", mfmt(r), formatC(p, format = "g", digits = 1))]
  d[, dark := abs(r) > 0.6]
  d[, key := module %in% km]

  pE <- ggplot(d, aes(trait, module, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = lab, colour = dark), size = (BASE - 2.5) * ppt, family = FONT,
              lineheight = 0.9) +
    scale_fill_gradientn(colours = pal_div, limits = c(-1, 1), breaks = c(-1, 0, 1),
                         name = "Pearson r",
                         guide = guide_colourbar(barheight = unit(24, "mm"),
                                                 barwidth = unit(2.2, "mm"))) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none") +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_pub_grid() +
    theme(axis.ticks = element_blank(),
          axis.text.x = element_text(size = BASE - 2, lineheight = 0.9),
          axis.text.y = element_text(face = ifelse(levels(d$module) %in% km, "bold", "plain")))
  keep(pE, "Fig1e_module_trait", 92, 150)
}

# ===================================================================================================
# f - module membership vs gene significance, for the key modules
# ===================================================================================================
if (stage %in% c("F", "panels", "all")) {
  gt <- fread(file.path(bulk, "WGCNA_gene_module_GS_MM.csv"))
  km <- key_modules()
  d <- rbindlist(lapply(km, function(m) {
    g <- gt[module == m]
    data.table(module = m, MM = abs(g[[paste0("MM_", m)]]), GS = abs(g$GS_T2D))
  }))
  d[, module := factor(module, levels = km)]
  ann <- d[, {
    ct <- suppressWarnings(cor.test(MM, GS))
    .(lab = paste0("italic(r)*' = ", formatC(ct$estimate, format = "f", digits = 2), "'"))
  }, by = module]

  pF <- ggplot(d, aes(MM, GS)) +
    raster_pts(geom_point(size = 0.3, alpha = 0.35, colour = "grey35", stroke = 0)) +
    geom_smooth(method = "lm", formula = y ~ x, colour = col_dir[["Up"]], fill = col_dir[["Up"]],
                alpha = 0.15, linewidth = 0.4) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab), parse = TRUE, hjust = -0.12,
              vjust = 1.5, size = (BASE - 2) * ppt, family = FONT, inherit.aes = FALSE) +
    facet_wrap(~ module, ncol = 4) +
    scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
                       expand = expansion(mult = 0.12)) +
    labs(x = "Module membership |MM|", y = "Gene significance |GS| for T2D") +
    theme_pub() +
    theme(panel.spacing.x = unit(2.6, "mm"), panel.spacing.y = unit(1.8, "mm"))
  keep(pF, "Fig1f_MM_GS", 120, 62)
}

# ===================================================================================================
# g - DEG / key-module / hyperglycaemia gene-set overlap
# ===================================================================================================
if (stage %in% c("G", "panels", "all")) {
  deg <- fread(file.path(bulk, "DEG_T2D_vs_Control_all.csv"))
  gt  <- fread(file.path(bulk, "WGCNA_gene_module_GS_MM.csv"))
  gs  <- fread(file.path(bulk, "geneset_hyperglycemia.csv"))
  A <- deg[P.Value < PCUT & abs(logFC) > LFC, gene]   # DEGs
  B <- gt[module %in% key_modules(), gene]            # key-module genes
  C <- intersect(gs$gene, deg$gene)                   # hyperglycaemia set, restricted to the universe

  cnt <- c(A = length(setdiff(A, union(B, C))), B = length(setdiff(B, union(A, C))),
           C = length(setdiff(C, union(A, B))), AB = length(setdiff(intersect(A, B), C)),
           AC = length(setdiff(intersect(A, C), B)), BC = length(setdiff(intersect(B, C), A)),
           ABC = length(intersect(intersect(A, B), C)))

  # three circles on a fixed layout with the counts placed in the seven regions; ggvenn's defaults
  # are far too heavy for a 7 pt panel
  th <- seq(0, 2 * pi, length.out = 200)
  ctr <- data.table(set = c("A", "B", "C"), x = c(-0.55, 0.55, 0), y = c(0.32, 0.32, -0.62))
  circ <- ctr[, .(x = x + cos(th), y = y + sin(th)), by = set]
  lab <- data.table(x = c(-1.15, 1.15, 0, 0, -0.62, 0.62, 0),
                    y = c(0.62, 0.62, -1.05, 0.62, -0.32, -0.32, 0.02),
                    n = cnt[c("A", "B", "C", "AB", "AC", "BC", "ABC")],
                    bold = c(rep(FALSE, 6), TRUE))
  nm <- data.table(x = c(-1.5, 1.5, 0), y = c(1.62, 1.62, -1.95), set = c("A", "B", "C"),
                   lab = c(sprintf("DEGs\n(%s)", nfmt(length(A))),
                           sprintf("WGCNA key-module\ngenes (%s)", nfmt(length(B))),
                           sprintf("Hyperglycaemia set\n(%s)", nfmt(length(C)))))

  pG <- ggplot() +
    geom_polygon(data = circ, aes(x, y, fill = set), alpha = 0.32, colour = "grey25",
                 linewidth = 0.3) +
    geom_text(data = lab, aes(x, y, label = nfmt(n), fontface = ifelse(bold, "bold", "plain")),
              size = (BASE - 1) * ppt, family = FONT) +
    geom_text(data = nm, aes(x, y, label = lab), size = (BASE - 1.5) * ppt, family = FONT,
              lineheight = 0.95, colour = "grey15") +
    scale_fill_manual(values = c(A = pal_qual[2], B = pal_qual[1], C = pal_qual[6]), guide = "none") +
    coord_equal(xlim = c(-2.6, 2.6), ylim = c(-2.4, 2.1), clip = "off") +
    theme_void_pub()
  keep(pG, "Fig1g_venn", 89, 70)
  fwrite(data.table(region = names(cnt), genes = as.integer(cnt)),
         file.path(EXPORT_DIR, "fig1g_venn_counts.csv"))
}

# ===================================================================================================
# assemble - d is base graphics and cannot join a patchwork, so it is placed at layout time
# ===================================================================================================
if (stage %in% c("assemble", "all")) {
  need <- c("Fig1a_PCA_ComBat", "Fig1b_volcano", "Fig1c_soft_threshold", "Fig1e_module_trait",
            "Fig1f_MM_GS", "Fig1g_venn")
  ps <- lapply(need, function(n) {
    f <- file.path(panel_dir, paste0(n, ".rds"))
    if (!file.exists(f)) stop("missing panel ", n, "; run stage panels first")
    readRDS(f)$plot
  })
  names(ps) <- c("a", "b", "c", "e", "f", "g")
  L <- tag_letters(7)
  # patchwork maps design letters to the ORDER plots are added, not to panel names: A is the first
  # plot in the expression below, B the second, and so on. The visible tags are set separately.
  design <- "
AAAAAAADDDDDDD
AAAAAAADDDDDDD
BBBBBBBDDDDDDD
BBBBBBBDDDDDDD
CCCCCCCDDDDDDD
CCCCCCCDDDDDDD
EEEEEEEFFFFFFF
EEEEEEEFFFFFFF
"
  comp <- ps$a + ps$b + ps$c + ps$e + ps$f + ps$g +
    plot_layout(design = design) +
    plot_annotation(tag_levels = list(c(L[1], L[2], L[3], L[5], L[6], L[7]))) &
    theme(plot.tag = element_text(size = BASE + 1, face = "bold", family = FONT))
  save_fig(comp, "Fig1_composite", W_2COL, 205, dir = fig_dir, tiff = TRUE)
}

# ===================================================================================================
# legend
# ===================================================================================================
if (stage %in% c("legend", "panels", "all")) {
  deg <- fread(file.path(bulk, "DEG_T2D_vs_Control_all.csv"))
  mt  <- fread(file.path(bulk, "WGCNA_module_trait.csv"))
  isum <- fread(file.path(bulk, "intersection_summary.csv"))
  nup <- deg[P.Value < PCUT & logFC > LFC, .N]; ndn <- deg[P.Value < PCUT & logFC < -LFC, .N]
  km <- key_modules()
  pwr <- readRDS(file.path(bulk, "WGCNA_net.rds"))$power

  l <- c(
    "# Figure 1 legend - bulk liver transcriptome: differential expression and co-expression modules",
    "",
    "**Fig. 1 | Differentially expressed genes and T2D-associated co-expression modules in human liver.**",
    "**a**, Principal components of the merged GSE15653 + GSE64998 liver arrays before and after ComBat batch correction (27 samples: 16 T2D, 11 lean control).",
    sprintf("**b**, Volcano plot of T2D vs control (limma, ~condition + dataset). Dashed lines, |log2 fold change| > %.1f and nominal P < %.2f; %d genes up and %d down. The 10 most significant genes in each direction are labelled.",
            LFC, PCUT, nup, ndn),
    sprintf("**c**, Scale-free topology fit and mean connectivity against soft-thresholding power; power %d (red) was the lowest reaching a fit of 0.8 (dashed).", pwr),
    "**d**, Gene dendrogram from the signed-hybrid network with the assigned module colours beneath.",
    sprintf("**e**, Pearson correlation between each module eigengene and T2D status (all 27 samples and within each cohort separately), HbA1c (GSE15653 only, n = 14) and dataset of origin. Each cell gives r above P. Key modules, in bold, are those with P < 0.1 for T2D pooled AND the same direction at P < 0.1 in both cohorts (decision M11b): %s.",
            paste(km, collapse = ", ")),
    "**f**, Gene significance for T2D against module membership within each key module; line, linear fit.",
    sprintf("**g**, Overlap of the %d DEGs, the %s genes in the WGCNA key modules and the %s hyperglycaemia-associated genes present in the expression universe; the %s genes shared by all three (bold) are the candidate set carried into Fig. 2.",
            nup + ndn, nfmt(isum[item == "module_genes", as.integer(value)]),
            nfmt(isum[item == "geneset_in_universe", as.integer(value)]),
            isum[item == "triple_intersection", value]),
    "",
    "## Caveats that belong in the text, not the figure",
    "",
    "- † **e**, the Dataset column is computed on ComBat-corrected expression. ComBat removes the batch mean by construction, so this column is bounded near zero (max |r| = 0.08) and is *not* an independent test for residual batch effect. It is reported, not used as a filter; the per-cohort columns are what the key-module rule tests, and the honest batch check is **a**.",
    "- Three modules with a strong pooled T2D correlation are excluded because they do not replicate: blue (r = −0.72 in GSE15653 vs +0.11 in GSE64998, i.e. the sign flips), green (0.69 vs 0.23) and turquoise (0.63 vs 0.16). blue and turquoise are also the two largest modules, so this is what cuts the key-module gene pool from 4,022 to 418.",
    sprintf("- No gene reaches FDR < 0.05 for T2D vs control; **b** therefore uses nominal P < %.2f and the DEG set is exploratory (decision M9).", PCUT),
    "- **e**, module-trait P values are uncorrected across 14 modules x 3 traits.",
    "- **f**, MM and GS are both computed from the same 27 samples, so the correlation is not independent evidence of module relevance.",
    "",
    "## Statistics to quote",
    "",
    "| Quantity | Value |",
    "|---|---|",
    sprintf("| Samples | 16 T2D, 11 lean control (GSE15653 + GSE64998) |"),
    sprintf("| Genes tested | %s |", nfmt(nrow(deg))),
    sprintf("| DEGs (P < %.2f, |log2FC| > %.1f) | %d up, %d down |", PCUT, LFC, nup, ndn),
    sprintf("| Genes at FDR < 0.05 | %d |", deg[adj.P.Val < 0.05, .N]),
    sprintf("| Soft-thresholding power | %d |", pwr),
    sprintf("| Modules (excluding grey) | %d |", nrow(mt)),
    sprintf("| Key modules (T2D P < 0.1) | %s |", paste(km, collapse = ", ")),
    sprintf("| Strongest module-trait correlation | %s, r = %s, %s |",
            mt[which.max(abs(r_T2D)), module], mfmt(mt[which.max(abs(r_T2D)), r_T2D]),
            p_short(mt[which.max(abs(r_T2D)), p_T2D])),
    sprintf("| Three-way intersection | %s genes |", isum[item == "triple_intersection", value])
  )
  writeLines(l, file.path(fig_dir, "Fig1_legend.md"))
  message("wrote ", file.path(fig_dir, "Fig1_legend.md"))
}

if (stage %in% c("assemble", "all", "panels")) {
  writeLines(capture.output(sessionInfo()), file.path(EXPORT_DIR, "sessionInfo_28.txt"))
}
