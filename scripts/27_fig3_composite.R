# 27_fig3_composite.R - publication composite for the scRNA-seq analysis (Fig 3), built panel by panel.
#
# Usage: Rscript scripts/27_fig3_composite.R <tissue> <stage>
#   tissue : liver | kidney
#   stage  : cache | A | B | C | D | E | F | G | panels | assemble | legend | all
#
# Panels:
#   a  UMAP by cell type, legend with per-type proportion
#   b  canonical-marker dot plot (size = % expressed, colour = mean log-normalized expression)
#   c  Control/STZ composition per cell type + cells per cell type (log10)
#   d  violin of gene 1, Control vs STZ
#   e  violin of gene 2, Control vs STZ
#   f  joint co-expression density of the pair on the UMAP
#   g  cell-level correlation of the pair in the key cell type
#
# Figure conventions (Nature/Cell): drawn at final print size (double column, 180 mm), Arial 5-7 pt,
# no in-panel titles or subtitles, no explanatory prose inside the axes. Every statistic that would
# otherwise sit in a panel title is written to Fig3_legend.md by the `legend` stage, which is the text
# the figure legend should carry. Colour maps are perceptually uniform (viridis); no rainbow/jet.
#
# The pair is the p53 arrest readout CDKN1A-CCND1 (decision R23, scripts 18/23); the key cell type and the
# background-matched empirical p come from results/pathway_selection/18_axis_correlation_<tissue>.csv.
#
# Significance brackets in d/e report the PER-MOUSE test (n = 3 vs 4 mice), not the cell-level Wilcoxon:
# cells within a mouse are not independent replicates, so the cell-level p is anticonservative. It is
# still reported, in the legend. Set SIG_UNIT <- "cell" to put the cell-level stars on the panels instead.
#
# Stage `cache` is the only stage that touches the Seurat object; it writes a small cache that every panel
# stage reads, so panels can be rebuilt one at a time without reloading ~200 MB.
#
# Inputs:  results/sc/<tissue>_annotated.rds, results/pathway_selection/18_axis_correlation_<tissue>.csv
# Outputs: results/figure_exports/fig3_composite_<tissue>_cache.rds
#          figures/Fig3_composite_<tissue>/Fig3<panel>_*.pdf (one file per panel)
#          figures/Fig3_composite_<tissue>/Fig3_composite_<tissue>.pdf
#          figures/Fig3_composite_<tissue>/Fig3_legend.md
#          results/figure_exports/Fig3_composite_<tissue>.tiff (submission raster)
#          results/figure_exports/fig3_panels_<tissue>/*.rds (intermediates for `assemble`)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

args   <- commandArgs(trailingOnly = TRUE)
tissue <- tolower(if (length(args) >= 1) args[1] else "liver")
stage  <- if (length(args) >= 2) args[2] else "all"
stopifnot(tissue %in% c("liver", "kidney"))
set.seed(20260918)

pair      <- c("Cdkn1a", "Ccnd1")
fig_dir   <- file.path("figures", paste0("Fig3_composite_", tissue))
cache_dir <- "results/figure_exports"
cache_f   <- file.path(cache_dir, paste0("fig3_composite_", tissue, "_cache.rds"))
panel_dir <- file.path(cache_dir, paste0("fig3_panels_", tissue))   # intermediates, not deliverables
tiff_dir  <- cache_dir                                              # TIFF exports live with the others
for (d in c(fig_dir, cache_dir, panel_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------------------------------
# figure constants
# ---------------------------------------------------------------------------------------------------
FONT      <- if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else "sans"
BASE      <- 7          # pt; Nature body size for figure text is 5-7 pt
TAG_CASE  <- "lower"    # "lower" = Nature (bold a, b, c), "upper" = Cell (A, B, C)
SIG_UNIT  <- "mouse"    # "mouse" (per-mouse t-test) or "cell" (cell-level Wilcoxon)
W_FULL    <- 180        # mm, double-column width

# Control / STZ, and a 13-colour qualitative set that survives greyscale and CVD simulation
cols_grp <- c(Control = "#4E79A7", STZ = "#E15759")
pal_ct   <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1",
              "#FF9DA7", "#9C755F", "#BAB0AC", "#86BCB6", "#D37295", "#8CD17D")

ppt <- 1 / .pt          # ggplot geom text `size` is mm; multiply pt by this

# ---------------------------------------------------------------------------------------------------
# canonical marker sets (the sets script 08 annotated with), in the order they are drawn in panel b
# ---------------------------------------------------------------------------------------------------
shared_sets <- list(
  B_cell        = c("Cd79a", "Cd79b", "Ms4a1", "Cd19"),
  Plasma_cell   = c("Jchain", "Mzb1", "Xbp1", "Igkc"),
  T_cell        = c("Cd3e", "Cd3d", "Trac", "Cd8a", "Cd4"),
  NK_cell       = c("Nkg7", "Klrb1c", "Ncr1", "Gzma"),
  Neutrophil    = c("S100a8", "S100a9", "Retnlg", "Csf3r"),
  Monocyte      = c("Ly6c2", "Ccr2", "Plac8", "Chil3"),
  cDC           = c("Xcr1", "Clec9a", "Cd209a", "Flt3"),
  pDC           = c("Siglech", "Bst2", "Ccr9"),
  Proliferating = c("Mki67", "Top2a", "Stmn1")
)
marker_sets <- list(
  liver = c(list(
    Hepatocyte       = c("Alb", "Apoa1", "Apoc3", "Cyp2e1", "Ttr"),
    Cholangiocyte    = c("Krt19", "Epcam", "Spp1", "Sox9"),
    LSEC             = c("Stab2", "Clec4g", "Dnase1l3", "Kdr"),
    Kupffer_cell     = c("Clec4f", "Vsig4", "Cd5l", "Marco"),
    Macrophage_other = c("C1qa", "C1qb", "Adgre1", "Cx3cr1"),
    HSC_Fibroblast   = c("Dcn", "Reln", "Lrat", "Colec11")
  ), shared_sets),
  kidney = c(list(
    PT                       = c("Lrp2", "Slc34a1", "Kap", "Miox"),
    LoH                      = c("Umod", "Slc12a1"),
    Distal_tubule            = c("Atp1b1", "Atp1a1", "Wfdc2"),
    CD_IC                    = c("Atp6v1g3", "Atp6v0d2"),
    Urothelium               = c("Krt7", "Sprr1a", "Foxq1"),
    Tubular_epithelium_other = c("Cdh16", "Pkhd1", "Bicc1"),
    Endothelial              = c("Emcn", "Plvap", "Pecam1", "Ehd3"),
    Fibroblast_Pericyte      = c("Pdgfrb", "Col1a1", "Dcn", "Rgs5"),
    Macrophage               = c("C1qa", "C1qb", "Adgre1", "Cd68")
  ), shared_sets)
)[[tissue]]

# ---------------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------------
# journal-readable cell type names; the underscore labels stay in the data and the result tables
pretty_ct <- function(x) {
  x <- as.character(x)
  x <- sub("^HSC_Fibroblast$",           "HSC/fibroblast",             x)
  x <- sub("^Macrophage_other$",         "Macrophage (other)",         x)
  x <- sub("^Fibroblast_Pericyte$",      "Fibroblast/pericyte",        x)
  x <- sub("^Tubular_epithelium_other$", "Tubular epithelium (other)", x)
  x <- sub("^CD_IC$",                    "CD-IC",                      x)
  x <- sub("^CD_PC$",                    "CD-PC",                      x)
  gsub("_", " ", x)
}
# minus sign and thousands separator as a journal copy-editor would set them
nfmt  <- function(n) formatC(n, big.mark = ",", format = "d")
mfmt  <- function(x, d = 2) sub("-", "−", formatC(x, format = "f", digits = d))
sup_digits <- function(n) {
  g <- c("⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹")
  paste(g[as.integer(strsplit(as.character(abs(n)), "")[[1]]) + 1], collapse = "")
}
p_short <- function(p) {
  if (is.na(p)) return("n.d.")
  if (p < 2.2e-16) return("P < 2.2 × 10⁻¹⁶")
  if (p < 0.001) {
    e <- floor(log10(p))
    return(sprintf("P = %s × 10⁻%s", formatC(p / 10^e, format = "f", digits = 1),
                   sup_digits(-e)))
  }
  # two significant figures, without the trailing zeros formatC(flag = "#") would leave behind
  sprintf("P = %s", formatC(p, format = "f", digits = min(3, max(2, -floor(log10(p)) + 1))))
}
stars <- function(p) {
  if (is.na(p)) "n.s." else if (p < 1e-4) "****" else if (p < 1e-3) "***" else
    if (p < 1e-2) "**" else if (p < 0.05) "*" else "n.s."
}

theme_pub <- function(base = BASE) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(
      axis.line         = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks        = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(1.2, "pt"),
      axis.text         = element_text(size = base - 1, colour = "black"),
      axis.title        = element_text(size = base, colour = "black"),
      legend.title      = element_text(size = base - 1, colour = "black"),
      legend.text       = element_text(size = base - 1, colour = "black"),
      legend.key.size   = unit(2.6, "mm"),
      legend.margin     = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "pt"),
      legend.background = element_blank(),
      plot.title        = element_blank(),   # no panel titles: everything goes in the legend
      plot.subtitle     = element_blank(),
      plot.caption      = element_blank(),
      plot.margin       = margin(2, 2, 2, 2),
      plot.tag          = element_text(size = base + 1, face = "bold", family = FONT, hjust = 0, vjust = 1)
    )
}
theme_umap <- function(base = BASE) {
  theme_pub(base) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank(), axis.title = element_blank())
}
# small corner arrows in place of full UMAP axes
umap_axes <- function(d, frac = 0.16, base = BASE) {
  rx <- diff(range(d$UMAP_1)); ry <- diff(range(d$UMAP_2))
  x0 <- min(d$UMAP_1) - 0.02 * rx; y0 <- min(d$UMAP_2) - 0.06 * ry
  ar <- arrow(length = unit(1.1, "mm"), type = "closed")
  list(
    annotate("segment", x = x0, xend = x0 + frac * rx, y = y0, yend = y0,
             arrow = ar, linewidth = 0.3, colour = "black"),
    annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + frac * ry,
             arrow = ar, linewidth = 0.3, colour = "black"),
    annotate("text", x = x0 + 0.5 * frac * rx, y = y0 - 0.035 * ry, label = "UMAP 1",
             size = (base - 2) * ppt, family = FONT, vjust = 1, hjust = 0.5),
    annotate("text", x = x0 - 0.025 * rx, y = y0 + 0.5 * frac * ry, label = "UMAP 2",
             size = (base - 2) * ppt, family = FONT, angle = 90, vjust = 0, hjust = 0.5),
    scale_x_continuous(expand = expansion(mult = c(0.09, 0.02))),
    scale_y_continuous(expand = expansion(mult = c(0.09, 0.02)))
  )
}
raster_pts <- function(g) if (requireNamespace("ggrastr", quietly = TRUE)) ggrastr::rasterise(g, dpi = 600) else g

# widths/heights are millimetres at final print size. figures/ holds deliverable PDFs only; the
# serialized ggplot that `assemble` re-reads is an intermediate and goes to results/figure_exports/.
save_panel <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, units = "mm",
         device = cairo_pdf)
  saveRDS(list(plot = p, width = w, height = h), file.path(panel_dir, paste0(name, ".rds")))
  message("wrote ", name)
  invisible(p)
}
get_cache <- function() {
  if (!file.exists(cache_f)) stop("cache missing; run stage `cache` first")
  cc <- readRDS(cache_f)
  cc$md[, ct := factor(pretty_ct(cell_type), levels = pretty_ct(cc$ct_levels))]
  cc$ct_pretty <- pretty_ct(cc$ct_levels)
  cc$key_pretty <- pretty_ct(cc$key)
  cc
}

# ===================================================================================================
# stage: cache - the only stage that loads the Seurat object
# ===================================================================================================
if (stage %in% c("cache", "all")) {
  suppressPackageStartupMessages({ library(Seurat); library(Nebulosa); library(Matrix) })
  so <- readRDS(file.path("results/sc", paste0(tissue, "_annotated.rds")))
  DefaultAssay(so) <- "RNA"

  md <- data.table(cell = colnames(so),
                   cell_type = as.character(so$cell_type),
                   condition = factor(as.character(so$condition), levels = c("Control", "STZ")),
                   mouse_id  = as.character(so$mouse_id))
  um <- Embeddings(so, "umap")
  md[, `:=`(UMAP_1 = um[, 1], UMAP_2 = um[, 2])]

  ct_levels <- md[, .N, by = cell_type][order(-N), cell_type]
  md[, cell_type := factor(cell_type, levels = ct_levels)]

  dat <- LayerData(so, assay = "RNA", layer = "data")

  # ---- panel b input: mean expression and % expressed per cell type ----
  genes_b <- unique(unlist(marker_sets)); genes_b <- genes_b[genes_b %in% rownames(dat)]
  ct_f <- factor(md$cell_type, levels = ct_levels)
  mm   <- Matrix::sparse.model.matrix(~ 0 + ct_f); colnames(mm) <- levels(ct_f)
  n_ct <- colSums(mm)
  sub  <- dat[genes_b, , drop = FALSE]
  avg  <- as.matrix(sub %*% mm) / rep(n_ct, each = length(genes_b))
  pct  <- as.matrix((sub > 0) %*% mm) / rep(n_ct, each = length(genes_b)) * 100
  dot  <- melt(as.data.table(avg, keep.rownames = "gene"), id.vars = "gene",
               variable.name = "cell_type", value.name = "avg_expr")
  dot  <- merge(dot, melt(as.data.table(pct, keep.rownames = "gene"), id.vars = "gene",
                          variable.name = "cell_type", value.name = "pct_expr"),
                by = c("gene", "cell_type"))
  gene_order <- unique(unlist(lapply(marker_sets, function(g) intersect(g, genes_b))))
  dot[, gene := factor(gene, levels = gene_order)]
  dot[, cell_type := factor(as.character(cell_type), levels = ct_levels)]

  # ---- panels d/e/f/g input: pair expression per cell ----
  pr <- intersect(pair, rownames(dat))
  stopifnot(length(pr) == 2)
  ex <- as.data.table(as.matrix(t(dat[pr, , drop = FALSE])))
  md <- cbind(md, ex)

  # cell-level Wilcoxon per gene (whole tissue) + per-mouse t-test on the mouse means
  gene_tests <- rbindlist(lapply(pr, function(g) {
    w  <- suppressWarnings(wilcox.test(md[[g]][md$condition == "STZ"], md[[g]][md$condition == "Control"]))
    pm <- md[, .(m = mean(get(g))), by = .(mouse_id, condition)]
    tt <- suppressWarnings(t.test(pm[condition == "STZ", m], pm[condition == "Control", m]))
    data.table(gene = g, cell_p = w$p.value, mouse_p = tt$p.value,
               mean_ctl = mean(md[[g]][md$condition == "Control"]),
               mean_stz = mean(md[[g]][md$condition == "STZ"]),
               pct_ctl = 100 * mean(md[[g]][md$condition == "Control"] > 0),
               pct_stz = 100 * mean(md[[g]][md$condition == "STZ"] > 0))
  }))
  per_mouse <- melt(md[, c("mouse_id", "condition", pr), with = FALSE],
                    id.vars = c("mouse_id", "condition"), variable.name = "gene", value.name = "expr"
                    )[, .(expr = mean(expr)), by = .(gene, mouse_id, condition)]

  # ---- panel f input: Nebulosa joint density ----
  dens <- NULL
  pf <- try(plot_density(so, features = pr, joint = TRUE, reduction = "umap"), silent = TRUE)
  if (!inherits(pf, "try-error")) {
    dd <- as.data.table(pf[[length(pf)]]$data)
    vcol <- names(dd)[vapply(dd, is.numeric, TRUE)]
    vcol <- setdiff(vcol, grep("^(umap|UMAP)", names(dd), value = TRUE))
    dens <- data.table(UMAP_1 = md$UMAP_1, UMAP_2 = md$UMAP_2, density = dd[[vcol[1]]])
  } else {
    message("Nebulosa failed: ", as.character(pf))
  }

  # ---- panel g input: key cell type + background statistics already computed by script 18 ----
  corr_f <- file.path("results/pathway_selection", paste0("18_axis_correlation_", tissue, ".csv"))
  corr   <- if (file.exists(corr_f)) fread(corr_f) else NULL
  key    <- if (!is.null(corr)) corr$cell_type[1] else as.character(ct_levels[1])

  saveRDS(list(tissue = tissue, pair = pr, md = md, ct_levels = ct_levels, dot = dot,
               gene_tests = gene_tests, per_mouse = per_mouse, dens = dens,
               corr = corr, key = key, n_cells = nrow(md)),
          cache_f)
  message("cache written: ", cache_f, " (", nrow(md), " cells, ", length(ct_levels), " cell types)")
  rm(so, dat); gc(verbose = FALSE)
}

# ===================================================================================================
# panel a - UMAP by cell type
# ===================================================================================================
if (stage %in% c("A", "panels", "all")) {
  cc <- get_cache(); md <- cc$md
  prop <- md[, .N, by = ct][order(-N)]
  prop[, lab := sprintf("%s  %.1f%%", ct, 100 * N / sum(N))]
  cols_ct <- setNames(pal_ct[seq_len(nrow(prop))], as.character(prop$ct))
  labpos <- md[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = ct]

  pA <- ggplot(md, aes(UMAP_1, UMAP_2, colour = ct)) +
    raster_pts(geom_point(size = 0.18, alpha = 0.8, stroke = 0)) +
    ggrepel::geom_text_repel(data = labpos, aes(label = ct), colour = "black",
                             size = (BASE - 1.5) * ppt, family = FONT,
                             segment.colour = "grey45", segment.size = 0.2,
                             min.segment.length = 0.2, max.overlaps = Inf, box.padding = 0.22,
                             bg.colour = "white", bg.r = 0.14, show.legend = FALSE) +
    scale_colour_manual(values = cols_ct, breaks = as.character(prop$ct), labels = prop$lab,
                        name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 1.4, alpha = 1), ncol = 1,
                                 keyheight = unit(2.9, "mm"))) +
    umap_axes(md) + theme_umap()
  save_panel(pA, "Fig3a_umap_celltype", 88, 62)
}

# ===================================================================================================
# panel b - canonical marker dot plot
# ===================================================================================================
if (stage %in% c("B", "panels", "all")) {
  cc <- get_cache(); dot <- copy(cc$dot)
  # only the marker sets whose cell type is actually present, and at most 3 markers each: 45+ columns
  # cannot be set legibly in half a double-column width
  sets <- marker_sets[names(marker_sets) %in% as.character(cc$ct_levels)]
  keep <- unique(unlist(lapply(sets, function(g) head(intersect(g, levels(dot$gene)), 3))))
  dot <- dot[gene %in% keep]
  dot[, gene := factor(as.character(gene), levels = keep)]
  dot[, ct := factor(pretty_ct(cell_type), levels = rev(cc$ct_pretty))]
  pB <- ggplot(dot[pct_expr > 0], aes(gene, ct)) +
    geom_point(aes(size = pct_expr, colour = avg_expr)) +
    scale_size_continuous(range = c(0.1, 2.1), name = "% expressed",
                          breaks = c(25, 50, 75, 100), limits = c(0, 100)) +
    scale_colour_viridis_c(option = "viridis", name = "Mean expression",
                           guide = guide_colourbar(barheight = unit(2.2, "mm"),
                                                   barwidth = unit(12, "mm"),
                                                   title.position = "top", title.hjust = 0)) +
    guides(size = guide_legend(title.position = "top", nrow = 1, title.hjust = 0,
                               override.aes = list(colour = "grey30"))) +
    scale_x_discrete(expand = expansion(add = 0.7)) +
    labs(x = NULL, y = NULL) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = BASE - 2,
                                     face = "italic"),
          axis.text.y = element_text(size = BASE - 1),
          legend.position = "bottom", legend.box = "horizontal",
          legend.justification = "center", legend.box.just = "bottom",
          legend.spacing.x = unit(2, "mm"), legend.key.size = unit(2.1, "mm"),
          legend.text = element_text(size = BASE - 2),
          legend.title = element_text(size = BASE - 1.5),
          plot.margin = margin(2, 5, 2, 2))
  save_panel(pB, "Fig3b_marker_dotplot", 92, 62)
}

# ===================================================================================================
# panel c - composition by group + cells per cell type
# ===================================================================================================
if (stage %in% c("C", "panels", "all")) {
  cc <- get_cache(); md <- cc$md
  lev  <- rev(cc$ct_pretty)
  comp <- md[, .N, by = .(ct, condition)][, frac := N / sum(N), by = ct]
  comp[, ct := factor(as.character(ct), levels = lev)]
  tot  <- md[, .N, by = ct][, ct := factor(as.character(ct), levels = lev)]

  pC1 <- ggplot(comp, aes(frac, ct, fill = condition)) +
    geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols_grp, name = NULL) +
    scale_x_continuous(labels = function(x) paste0(x * 100), expand = c(0, 0),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = "Cells (%)", y = NULL) +
    theme_pub() +
    theme(legend.position = c(0.5, 1.07), legend.direction = "horizontal",
          legend.key.size = unit(2.4, "mm"), plot.margin = margin(9, 7, 2, 2))
  pC2 <- ggplot(tot, aes(N, ct)) +
    geom_col(fill = "grey45", width = 0.7) +
    scale_x_log10(expand = expansion(mult = c(0, 0.14)), limits = c(1, NA),
                  breaks = c(1, 100, 10000),
                  labels = c("1", expression(10^2), expression(10^4))) +
    labs(x = "Cells (n)", y = NULL) +
    theme_pub() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(), plot.margin = margin(9, 3, 2, 7))
  pC <- (pC1 | pC2) + plot_layout(widths = c(1, 0.55))
  save_panel(pC, "Fig3c_composition", 95, 62)
}

# ===================================================================================================
# panels d and e - pair expression, Control vs STZ
# ===================================================================================================
# Statistics for one gene of the pair, restricted to the key cell type. The reference framework
# (Xu et al. Fig. 3D-E) compares the pair between groups WITHIN the key cell type, not across the whole
# atlas, and scripts 18/23 do the same; an atlas-wide violin mixes cell types and is not comparable to
# panel g, which is key-cell-type only.
pair_stats <- function(cc, g) {
  md <- cc$md[cell_type == cc$key]
  w  <- suppressWarnings(wilcox.test(md[[g]][md$condition == "STZ"],
                                     md[[g]][md$condition == "Control"]))
  pm <- md[, .(expr = mean(get(g)), n = .N), by = .(mouse_id, condition)][n >= 10]
  tt <- if (pm[condition == "STZ", .N] >= 2 && pm[condition == "Control", .N] >= 2)
    suppressWarnings(t.test(pm[condition == "STZ", expr], pm[condition == "Control", expr])) else NULL
  list(md = md, per_mouse = pm,
       cell_p  = w$p.value,
       mouse_p = if (is.null(tt)) NA_real_ else tt$p.value,
       det     = md[, .(pct = 100 * mean(get(g) > 0)), by = condition],
       n_cells = nrow(md))
}
build_violin <- function(cc, g) {
  st <- pair_stats(cc, g)
  md <- st$md
  pm <- st$per_mouse
  p_show <- if (SIG_UNIT == "mouse") st$mouse_p else st$cell_p
  ymax <- max(md[[g]])
  # most cells do not detect these genes, so the detection rate carries as much of the signal as the
  # violin shape does; it is annotated rather than left as white space
  det <- st$det
  ggplot(md, aes(condition, .data[[g]], fill = condition)) +
    geom_violin(scale = "width", trim = TRUE, colour = "black", linewidth = 0.25, width = 0.82) +
    geom_point(data = pm, aes(condition, expr), shape = 21, size = 0.9, fill = "white",
               colour = "black", stroke = 0.3, inherit.aes = FALSE,
               position = position_jitter(width = 0.09, height = 0)) +
    geom_text(data = det, aes(condition, y = -ymax * 0.055,
                              label = sprintf("%.1f%%", pct)), inherit.aes = FALSE,
              size = (BASE - 2.5) * ppt, family = FONT, colour = "grey25") +
    annotate("segment", x = 1, xend = 2, y = ymax * 1.07, yend = ymax * 1.07, linewidth = 0.25) +
    annotate("segment", x = 1, xend = 1, y = ymax * 1.03, yend = ymax * 1.07, linewidth = 0.25) +
    annotate("segment", x = 2, xend = 2, y = ymax * 1.03, yend = ymax * 1.07, linewidth = 0.25) +
    annotate("text", x = 1.5, y = ymax * 1.10, label = stars(p_show), vjust = 0,
             size = (BASE - 1) * ppt, family = FONT) +
    scale_fill_manual(values = cols_grp, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.17)),
                       breaks = scales::breaks_pretty(4)) +
    labs(x = NULL, y = bquote(italic(.(g)) ~ "expression")) +
    theme_pub()
}
if (stage %in% c("D", "panels", "all")) {
  cc <- get_cache(); save_panel(build_violin(cc, cc$pair[1]), "Fig3d_violin_pair1", 42, 55)
}
if (stage %in% c("E", "panels", "all")) {
  cc <- get_cache(); save_panel(build_violin(cc, cc$pair[2]), "Fig3e_violin_pair2", 42, 55)
}

# ===================================================================================================
# panel f - joint co-expression density
# ===================================================================================================
if (stage %in% c("F", "panels", "all")) {
  cc <- get_cache()
  if (is.null(cc$dens)) stop("no density in cache")
  d <- copy(cc$dens)
  d[, rel := density / max(density)]
  thr <- 0.02
  pF <- ggplot() +
    raster_pts(geom_point(data = d[rel <= thr], aes(UMAP_1, UMAP_2), colour = "grey88",
                          size = 0.18, stroke = 0)) +
    raster_pts(geom_point(data = d[rel > thr][order(rel)], aes(UMAP_1, UMAP_2, colour = rel),
                          size = 0.22, stroke = 0)) +
    scale_colour_viridis_c(option = "viridis", limits = c(0, 1), breaks = c(0, 0.5, 1),
                           name = "Joint density",
                           guide = guide_colourbar(barheight = unit(2.2, "mm"),
                                                   barwidth = unit(14, "mm"),
                                                   title.position = "top", ticks.colour = "white")) +
    umap_axes(d) + theme_umap() +
    theme(legend.position = c(0.76, 0.10), legend.direction = "horizontal",
          legend.title = element_text(size = BASE - 1, hjust = 0.5))
  save_panel(pF, "Fig3f_coexpression_density", 85, 62)
}

# ===================================================================================================
# panel g - cell-level correlation of the pair in the key cell type
# ===================================================================================================
if (stage %in% c("G", "panels", "all")) {
  cc <- get_cache(); md <- cc$md; g1 <- cc$pair[1]; g2 <- cc$pair[2]
  ek <- md[cell_type == cc$key]
  ct <- suppressWarnings(cor.test(ek[[g1]], ek[[g2]], method = "pearson"))
  # the numbers are quoted, not parsed as plotmath arithmetic: unparsed they would lose the trailing
  # zero (0.010 -> 0.01) and gain a space after the unary minus
  lab <- sprintf('italic(R)*" = %s, "*italic(P)*" = %s"',
                 sub("-", "−", formatC(unname(ct$estimate), format = "f", digits = 2)),
                 formatC(ct$p.value, format = "f", digits = 3))
  pG <- ggplot(ek, aes(.data[[g2]], .data[[g1]])) +
    raster_pts(geom_jitter(width = 0.03, height = 0.03, alpha = 0.28, size = 0.35,
                           colour = "grey35", stroke = 0)) +
    geom_smooth(method = "lm", formula = y ~ x, colour = cols_grp[["STZ"]],
                fill = cols_grp[["STZ"]], alpha = 0.16, linewidth = 0.5) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.12, vjust = 1.8, parse = TRUE, label = lab,
             size = (BASE - 1) * ppt, family = FONT) +
    labs(x = bquote(italic(.(g2)) ~ "expression"), y = bquote(italic(.(g1)) ~ "expression")) +
    theme_pub()
  save_panel(pG, "Fig3g_pair_correlation", 90, 55)
}

# ===================================================================================================
# stage: assemble
# ===================================================================================================
if (stage %in% c("assemble", "all")) {
  need <- c("Fig3a_umap_celltype", "Fig3b_marker_dotplot", "Fig3c_composition",
            "Fig3d_violin_pair1", "Fig3e_violin_pair2", "Fig3f_coexpression_density",
            "Fig3g_pair_correlation")
  ps <- lapply(need, function(n) {
    f <- file.path(panel_dir, paste0(n, ".rds"))
    if (!file.exists(f)) stop("missing panel ", n, "; run stage panels first")
    readRDS(f)$plot
  })
  design <- "
AAAAAAAAABBBBBBBBB
AAAAAAAAABBBBBBBBB
CCCCCCCCCFFFFFFFFF
CCCCCCCCCFFFFFFFFF
DDDEEEGGGGGGGGGGGG
DDDEEEGGGGGGGGGGGG
"
  # panel c is itself a two-plot patchwork, so it occupies two slots in the tag sequence; the blank
  # keeps the visible tags reading a-g.
  lets <- if (TAG_CASE == "lower") letters[1:7] else LETTERS[1:7]
  tags <- list(c(lets[1], lets[2], lets[3], "", lets[4], lets[5], lets[6], lets[7]))
  comp <- ps[[1]] + ps[[2]] + ps[[3]] + ps[[4]] + ps[[5]] + ps[[6]] + ps[[7]] +
    plot_layout(design = design) +
    plot_annotation(tag_levels = tags,
                    theme = theme(plot.margin = margin(2, 2, 2, 2))) &
    theme(plot.tag = element_text(size = BASE + 1, face = "bold", family = FONT))
  w <- W_FULL; h <- 186
  out <- file.path(fig_dir, paste0("Fig3_composite_", tissue))
  ggsave(paste0(out, ".pdf"), comp, width = w, height = h, units = "mm", device = cairo_pdf)
  tif <- file.path(tiff_dir, paste0("Fig3_composite_", tissue, ".tiff"))
  ggsave(tif, comp, width = w, height = h, units = "mm", dpi = 600, bg = "white",
         compression = "lzw")
  message("wrote ", out, ".pdf and ", tif, " (", w, " x ", h, " mm)")
}

# ===================================================================================================
# stage: legend - every number that was removed from the panels
# ===================================================================================================
if (stage %in% c("legend", "panels", "all")) {
  cc <- get_cache(); md <- cc$md; g1 <- cc$pair[1]; g2 <- cc$pair[2]
  gt <- rbindlist(lapply(cc$pair, function(g) {
    st <- pair_stats(cc, g)
    data.table(gene = g, cell_p = st$cell_p, mouse_p = st$mouse_p,
               mean_ctl = mean(st$md[[g]][st$md$condition == "Control"]),
               mean_stz = mean(st$md[[g]][st$md$condition == "STZ"]),
               pct_ctl = st$det[condition == "Control", pct],
               pct_stz = st$det[condition == "STZ", pct],
               n_cells = st$n_cells)
  }))
  cr <- cc$corr
  ek <- md[cell_type == cc$key]
  ctst <- suppressWarnings(cor.test(ek[[g1]], ek[[g2]], method = "pearson"))
  n_ct <- md[, .N, by = ct][order(-N)]
  sig_word <- if (SIG_UNIT == "mouse") "per-mouse t-test (n = 3 vs 4 mice)" else
    "cell-level Wilcoxon rank-sum test"
  ct_plural <- paste0(tolower(cc$key_pretty), if (grepl("[^s]$", cc$key_pretty)) "s" else "")

  l <- c(
    sprintf("# Figure 3 legend - single-cell atlas of %s, STZ vs Control", cc$tissue),
    "",
    sprintf("**Fig. 3 | Single-cell transcriptomes of STZ-diabetic and control mouse %s.**", cc$tissue),
    sprintf("**a**, UMAP of %s cells from 7 mice (3 Control, 4 STZ) coloured by cell type; the legend gives the proportion of cells in each type.",
            nfmt(cc$n_cells)),
    "**b**, Canonical markers used for annotation. Dot size, percentage of cells in the type with non-zero expression; colour, mean log-normalized expression.",
    "**c**, Left, proportion of each cell type contributed by Control and STZ mice. Right, number of cells per type (log scale).",
    sprintf("**d**,**e**, *%s* (**d**) and *%s* (**e**) expression in %s (%s cells), Control vs STZ. Violins are width-scaled over all cells of that type; white points are per-mouse means; the percentage under each violin is the fraction of cells with non-zero expression. Brackets give the %s.",
            g1, g2, ct_plural, nfmt(gt[1, n_cells]), sig_word),
    sprintf("**f**, Joint kernel density of *%s* and *%s* co-expression (Nebulosa) on the UMAP of **a**; grey, cells below 2%% of the maximum density.",
            g1, g2),
    sprintf("**g**, *%s* against *%s* in %s, the cell type with the strongest concordance with the bulk axis; line, linear fit with 95%% confidence band.",
            g1, g2, ct_plural),
    "",
    "## Statistics to quote",
    "",
    "| Quantity | Value |",
    "|---|---|",
    sprintf("| Cells (total) | %s |", nfmt(cc$n_cells)),
    sprintf("| Cell types | %d |", nrow(n_ct)),
    sprintf("| Mice | 3 Control, 4 STZ |"),
    sprintf("| %s, mean expression (Control / STZ) | %s / %s |", g1,
            mfmt(gt[gene == g1, mean_ctl], 3), mfmt(gt[gene == g1, mean_stz], 3)),
    sprintf("| %s, cells detected (Control / STZ) | %s%% / %s%% |", g1,
            mfmt(gt[gene == g1, pct_ctl], 1), mfmt(gt[gene == g1, pct_stz], 1)),
    sprintf("| Panels d,e cell type | %s (n = %s cells) |", cc$key_pretty, nfmt(gt[1, n_cells])),
    sprintf("| %s, cell-level Wilcoxon | %s |", g1, p_short(gt[gene == g1, cell_p])),
    sprintf("| %s, per-mouse t-test (n = 3 vs 4) | %s |", g1, p_short(gt[gene == g1, mouse_p])),
    sprintf("| %s, mean expression (Control / STZ) | %s / %s |", g2,
            mfmt(gt[gene == g2, mean_ctl], 3), mfmt(gt[gene == g2, mean_stz], 3)),
    sprintf("| %s, cells detected (Control / STZ) | %s%% / %s%% |", g2,
            mfmt(gt[gene == g2, pct_ctl], 1), mfmt(gt[gene == g2, pct_stz], 1)),
    sprintf("| %s, cell-level Wilcoxon | %s |", g2, p_short(gt[gene == g2, cell_p])),
    sprintf("| %s, per-mouse t-test (n = 3 vs 4) | %s |", g2, p_short(gt[gene == g2, mouse_p])),
    sprintf("| Panel g cell type | %s (n = %s cells) |", cc$key_pretty, nfmt(nrow(ek))),
    sprintf("| Panel g Pearson R | %s (%s) |", mfmt(unname(ctst$estimate)), p_short(ctst$p.value)),
    if (!is.null(cr)) sprintf("| Panel g, %% cells co-expressing both | %s%% |", mfmt(cr$pct_coexpr[1], 1)),
    if (!is.null(cr)) sprintf("| Panel g background null (detection-matched random pairs) | mean R = %s, 95th percentile = %s, empirical P = %s |",
                              mfmt(cr$bg_mean_R[1]), mfmt(cr$bg_q95[1]), mfmt(cr$empirical_p[1], 3)),
    "",
    "## Caveats that belong in the text, not the figure",
    "",
    sprintf("- The cell-level Wilcoxon in **d**,**e** treats cells as independent replicates, which they are not; the brackets therefore report the per-mouse test. The two disagree in both directions here: %s gives %s across cells and %s across mice, %s gives %s across cells and %s across mice. With 3 vs 4 mice the per-mouse test is the valid one but has very little power, so neither should be read as strong evidence.",
            g1, p_short(gt[gene == g1, cell_p]), p_short(gt[gene == g1, mouse_p]),
            g2, p_short(gt[gene == g2, cell_p]), p_short(gt[gene == g2, mouse_p])),
    if (!is.null(cr)) sprintf("- The %s-%s correlation in **g** is negative (R = %s) and does not exceed a detection-matched background (empirical P = %s), i.e. the two genes are not co-expressed within %s above chance.",
                              g1, g2, mfmt(unname(ctst$estimate)), mfmt(cr$empirical_p[1], 3),
                              ct_plural),
    "- Cell numbers are uneven between mice (STZ_4 contributes most fibroblasts); proportions in **c** should not be read as differential abundance without the propeller test in `results/sc/<tissue>_propeller_DA.csv`.",
    "",
    "## Cells per type",
    "",
    "| Cell type | Cells | % |",
    "|---|---|---|",
    sprintf("| %s | %s | %s |", n_ct$ct, nfmt(n_ct$N), mfmt(100 * n_ct$N / sum(n_ct$N), 2))
  )
  writeLines(Filter(Negate(is.null), l), file.path(fig_dir, "Fig3_legend.md"))
  message("wrote ", file.path(fig_dir, "Fig3_legend.md"))
}

if (stage %in% c("all", "panels", "assemble")) {
  writeLines(capture.output(sessionInfo()), file.path(cache_dir, "sessionInfo_27.txt"))
}
