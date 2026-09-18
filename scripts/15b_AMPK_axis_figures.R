# 15b_AMPK_axis_figures.R — Fig 3D-G analogue of the reference paper for the selected pathway
# (AMPK signalling, hsa04152; decisions R17/R18). Run after 15a_AMPK_axis_screen.R.
#
# Paper (Xu et al. 2026) -> this project
#   Fig 3D-E  Rac1/Pak1 expression in neutrophils by group (cell-level Wilcoxon, BH)
#             -> 3D: AMPK axis genes (pathview nodes) in the key cell type; cell-level Wilcoxon BH (paper) and
#                    per-mouse pseudobulk edgeR P (primary, script 09) are both shown
#             -> 3E: AMPK pathway scores (UCell) per cell type: mouse orthologs of AMPK-map genes up / down in
#                    human T2D liver (bulk P < 0.05); per-mouse means, Welch t-test
#   Fig 3F    Rac1+Pak1 joint density (highest in neutrophils) -> Sirt1+Ppargc1a joint density (Nebulosa) and
#             % co-expressing cells per cell type
#   Fig 3G    Rac1-Pak1 Pearson in neutrophils -> Sirt1-Ppargc1a: cell-level Pearson (paper) and per-mouse
#             pseudobulk Pearson
#
# Key cell type (3D) = top row of AMPK_screen_celltype_summary.csv among pseudobulk-tested cell types
#   (most genes concordant in both tests, then cell-level only).
# Pair cell type (3G) = key cell type if >= 5% of its cells co-express the pair; otherwise the cell type with the
#   highest co-expression among those with >= 10 cells per group (paper: correlation where co-localization is highest).
# Pair = SIRT1 -> PGC-1a edge of the AMPK map; both suppressed in human T2D (analogue of RAC1 -> PAK1).
#
# Figures: figures/Fig3_scRNA/Fig3{D,E,F,G}_AMPK_*_<tissue>.pdf
# Tables:  results/pathway_selection/AMPK_score_<tissue>.csv, AMPK_pair_correlation_<tissue>.csv

suppressPackageStartupMessages({
  library(Seurat)
  library(Nebulosa)
  library(UCell)
  library(edgeR)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})
set.seed(20260914)
fig_main <- "results/supplementary_figures/pathway_selection"; out_dir <- "results/pathway_selection"
  # draft panels; figures/ is owned by the 27-32 publication figure scripts (Fig 3 is now
  # Fig3_composite_<tissue> from 27 and Fig3_axis_<tissue> from 32). Writing here from an analysis
  # script silently overwrites them, because Windows filenames are case-insensitive.

dir.create(fig_main, recursive = TRUE, showWarnings = FALSE)
cols <- c(Control = "#3B7DD8", STZ = "#C8322F")
axis_genes <- c("Prkaa1", "Prkaa2", "Stk11", "Sirt1", "Ppargc1a", "Irs2", "Igf1", "Srebf1", "Fasn", "Fbp1", "Ccnd1")
pair <- c("Sirt1", "Ppargc1a")
fmt_p <- function(p) ifelse(is.na(p), "NA", formatC(p, format = "g", digits = 2))

pw   <- fread(file.path(out_dir, "AMPK_pathway_genes.csv"))
summ <- fread(file.path(out_dir, "AMPK_screen_celltype_summary.csv"))
sigs_all <- list(AMPK_T2D_up   = pw[human_DE == TRUE & human_logFC > 0 & !is.na(mouse), mouse],
                 AMPK_T2D_down = pw[human_DE == TRUE & human_logFC < 0 & !is.na(mouse), mouse])

for (tt in c("liver", "kidney")) {
  so  <- readRDS(file.path("results/sc", paste0(tt, "_annotated.rds")))
  scr <- fread(file.path(out_dir, paste0("AMPK_screen_", tt, ".csv")))
  key <- summ[summ$tissue == tt & pb_tested > 0][1, cell_type]
  message(tt, ": key cell type = ", key)

  # ---- Fig 3D: axis genes in key cell type ---------------------------------------------
  st <- scr[cell_type == key & mouse %in% axis_genes & pmax(pct_STZ, pct_Control) >= 0.05]
  genes_d <- axis_genes[axis_genes %in% st$mouse]
  sub <- subset(so, subset = cell_type == key)
  ex <- as.data.table(FetchData(sub, vars = c(genes_d, "condition", "mouse_id"), layer = "data"))
  long <- melt(ex, id.vars = c("condition", "mouse_id"), variable.name = "mouse", value.name = "expr")
  mm <- long[, .(expr = mean(expr), n = .N), by = .(mouse, condition, mouse_id)][n >= 10]
  st[, lab := sprintf("%s: human %s, STZ %s (%s)\ncell BH %s | pseudobulk P %s", mouse,
                      ifelse(human_logFC > 0, "up", "down"), ifelse(cell_log2FC > 0, "up", "down"),
                      ifelse(sign(cell_log2FC) == sign(human_logFC), "same", "opposite"), fmt_p(cell_BH), fmt_p(pb_P))]
  long[, gene_lab := factor(st$lab[match(mouse, st$mouse)], levels = st$lab[match(genes_d, st$mouse)])]
  mm[, gene_lab := factor(st$lab[match(mouse, st$mouse)], levels = levels(long$gene_lab))]
  p3d <- ggplot(long, aes(condition, expr, fill = condition)) +
    geom_violin(scale = "width", alpha = 0.45, colour = NA) +
    geom_point(data = mm, shape = 21, size = 2.4, colour = "black", position = position_jitter(width = 0.08, height = 0)) +
    facet_wrap(~gene_lab, scales = "free_y", nrow = 2) + scale_fill_manual(values = cols) +
    labs(x = NULL, y = "Log-normalized expression (points = mouse means)",
         title = sprintf("%s %s: AMPK axis genes, STZ vs Control", tt, key),
         subtitle = sprintf("Cells: Control %d, STZ %d. Cell-level Wilcoxon BH (paper) and per-mouse pseudobulk edgeR P",
                            sum(sub$condition == "Control"), sum(sub$condition == "STZ"))) +
    theme_bw(base_size = 10) + theme(legend.position = "none", strip.text = element_text(size = 7.5))
  ggsave(file.path(fig_main, paste0("Fig3D_AMPK_axis_expression_", tt, ".pdf")), p3d,
         width = 2.6 * ceiling(length(genes_d) / 2) + 1, height = 6.5)

  # ---- Fig 3E: AMPK pathway scores per cell type -----------------------------------------
  sigs <- lapply(sigs_all, intersect, rownames(so)); sigs <- sigs[lengths(sigs) >= 5]
  so <- AddModuleScore_UCell(so, features = sigs, assay = "RNA", slot = "counts", ncores = 4, name = "_UCell")
  sc <- as.data.table(so@meta.data)[, c("cell_type", "condition", "mouse_id", paste0(names(sigs), "_UCell")), with = FALSE]
  sl <- melt(sc, id.vars = c("cell_type", "condition", "mouse_id"), variable.name = "score", value.name = "value")
  pm <- sl[, .(value = mean(value), n = .N), by = .(cell_type, condition, mouse_id, score)][n >= 20]
  tst <- pm[, {
    a <- value[condition == "STZ"]; b <- value[condition == "Control"]
    p <- if (length(a) >= 2 && length(b) >= 2) tryCatch(t.test(a, b)$p.value, error = function(e) NA_real_) else NA_real_
    .(n_STZ = length(a), n_Control = length(b), diff = mean(a) - mean(b), p_welch = p)
  }, by = .(cell_type, score)]
  fwrite(tst[order(score, p_welch)], file.path(out_dir, paste0("AMPK_score_", tt, ".csv")))
  tst[, lab := ifelse(is.na(p_welch), "", ifelse(p_welch < 0.05, sprintf("p=%s", fmt_p(p_welch)), "ns"))]
  ymax <- pm[, .(y = max(value) * 1.08), by = .(cell_type, score)]
  tst <- merge(tst, ymax, by = c("cell_type", "score"))
  pm[, score := sub("_UCell$", "", score)]; tst[, score := sub("_UCell$", "", score)]
  p3e <- ggplot(pm, aes(cell_type, value, colour = condition)) +
    geom_point(position = position_dodge(width = 0.6), size = 2) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, position = position_dodge(width = 0.6), linewidth = 0.3) +
    geom_text(data = tst, aes(cell_type, y, label = lab), inherit.aes = FALSE, size = 2.8) +
    facet_wrap(~score, ncol = 1, scales = "free_y") + scale_colour_manual(values = cols) +
    labs(x = NULL, y = "Mean UCell score per mouse (>= 20 cells)",
         title = paste(tt, "- AMPK pathway genes changed in human T2D liver (mouse orthologs)"),
         subtitle = sprintf("Up set %d genes, down set %d genes; Welch t-test on mouse means",
                            length(sigs$AMPK_T2D_up), length(sigs$AMPK_T2D_down))) +
    theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
  ggsave(file.path(fig_main, paste0("Fig3E_AMPK_pathway_score_", tt, ".pdf")), p3e, width = 11, height = 7)

  # ---- Fig 3F: co-localization of the pair ------------------------------------------------
  dat <- LayerData(so, assay = "RNA", layer = "data")
  both <- dat[pair[1], ] > 0 & dat[pair[2], ] > 0
  md <- as.data.table(so@meta.data)[, both := both]
  co <- md[, .(pct = 100 * mean(both), n_ctl = sum(condition == "Control"), n_stz = sum(condition == "STZ")), by = cell_type]
  eligible <- co[n_ctl >= 10 & n_stz >= 10]
  pair_ct <- if (co[cell_type == key, pct] >= 5) key else eligible[order(-pct)][1, cell_type]
  co[, highlight := cell_type == pair_ct]
  pf <- plot_density(so, features = pair, joint = TRUE, reduction = "umap")
  pf <- if (inherits(pf, "patchwork")) pf else wrap_plots(pf)
  pb <- ggplot(co, aes(reorder(cell_type, pct), pct, fill = highlight)) + geom_col() + coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#C8322F", `FALSE` = "grey70"), guide = "none") +
    labs(x = NULL, y = sprintf("%% cells co-expressing %s & %s", pair[1], pair[2]), title = paste(tt, "- AMPK axis")) +
    theme_bw(base_size = 10)
  ggsave(file.path(fig_main, paste0("Fig3F_AMPK_colocalization_", tt, ".pdf")), pf | pb, width = 17, height = 4.5)

  # ---- Fig 3G: correlation of the pair ------------------------------------------------------
  subg <- subset(so, subset = cell_type == pair_ct)
  dg <- LayerData(subg, assay = "RNA", layer = "data")
  x <- dg[pair[1], ]; y <- dg[pair[2], ]
  pc <- cor.test(x, y, method = "pearson")
  bo <- x > 0 & y > 0
  sp <- if (sum(bo) >= 10) cor.test(x[bo], y[bo], method = "spearman", exact = FALSE) else NULL
  cnt <- LayerData(subg, assay = "RNA", layer = "counts")
  grp <- factor(subg$mouse_id); mmx <- sparse.model.matrix(~0 + grp); colnames(mmx) <- levels(grp)
  keep <- colnames(mmx)[colSums(mmx) >= 20]
  pbm <- as.matrix(cnt %*% mmx)[, keep, drop = FALSE]
  lcpm <- cpm(DGEList(pbm), log = TRUE, prior.count = 1)
  pp <- if (length(keep) >= 5) cor.test(lcpm[pair[1], ], lcpm[pair[2], ]) else NULL
  cond <- as.character(subg$condition[match(keep, subg$mouse_id)])
  ctab <- data.table(tissue = tt, cell_type = pair_ct, gene1 = pair[1], gene2 = pair[2], n_cells = length(x),
                     cell_pearson_r = unname(pc$estimate), cell_pearson_p = pc$p.value,
                     cells_coexpressing = sum(bo), pct_coexpressing = round(100 * mean(bo), 2),
                     coexpr_spearman_rho = if (is.null(sp)) NA_real_ else unname(sp$estimate),
                     coexpr_spearman_p = if (is.null(sp)) NA_real_ else sp$p.value,
                     n_mice = length(keep), pseudobulk_r = if (is.null(pp)) NA_real_ else unname(pp$estimate),
                     pseudobulk_p = if (is.null(pp)) NA_real_ else pp$p.value)
  fwrite(ctab, file.path(out_dir, paste0("AMPK_pair_correlation_", tt, ".csv")))
  print(ctab)
  dc <- data.table(x = x, y = y, condition = subg$condition)
  g1 <- ggplot(dc, aes(x, y)) + geom_jitter(aes(colour = condition), width = 0.03, height = 0.03, size = 0.5, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.6) + scale_colour_manual(values = cols) +
    labs(x = pair[1], y = pair[2], title = sprintf("%s %s: all cells (paper style)", tt, pair_ct),
         subtitle = sprintf("Pearson R = %.2f, p = %s (%d cells; %.1f%% co-express)", ctab$cell_pearson_r,
                            fmt_p(ctab$cell_pearson_p), ctab$n_cells, ctab$pct_coexpressing)) +
    theme_bw(base_size = 11)
  g2 <- if (!is.null(pp)) {
    ggplot(data.table(x = lcpm[pair[1], ], y = lcpm[pair[2], ], condition = cond), aes(x, y)) +
      geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.6) + geom_point(aes(colour = condition), size = 3) +
      scale_colour_manual(values = cols) +
      labs(x = paste(pair[1], "log-CPM"), y = paste(pair[2], "log-CPM"), title = "Per-mouse pseudobulk",
           subtitle = sprintf("Pearson r = %.2f, p = %s (n = %d mice)", ctab$pseudobulk_r, fmt_p(ctab$pseudobulk_p), ctab$n_mice)) +
      theme_bw(base_size = 11)
  } else plot_spacer() + labs(title = "Pseudobulk: < 5 mice with >= 20 cells")
  ggsave(file.path(fig_main, paste0("Fig3G_AMPK_correlation_", tt, ".pdf")), g1 | g2, width = 11, height = 4.8)
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_15b.txt"))
