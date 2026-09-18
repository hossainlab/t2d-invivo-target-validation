# 18_selected_axis_figures.R — Fig 3D-G for the axis that survives the confirmation screen (script 17).
#
# Screen result (17): none of the DE-DE edges drawn on the shortlisted KEGG maps satisfies the reference paper's
# three Fig 3 criteria in this mouse atlas - the AMPK output arm is regulated post-translationally (the
# AMPK -> PPARGC1A edge is discordant in the human bulk data: PRKAA2 up, PPARGC1A down), and every edge that
# clears the correlation test does so with R ~ 0.1, i.e. on cell number rather than effect size.
# The one axis that does satisfy all three criteria is the p53 cell-cycle-arrest arm, CDKN1A (p21) - CCND1
# (cyclin D1): 3rd-ranked KEGG term on the intersect genes, the top Hallmark GSEA term on the full human ranking
# (NES 1.92, padj 1.5e-5), CDKN1A the only FDR-significant human DEG in that map, CCND1 in the 19-gene candidate
# set and replicated in GSE23343.
#
# Every correlation and co-localization statistic here is additionally referred to a background null of random
# gene pairs matched on detection rate in the same cell type, because with thousands of cells a cell-level
# Pearson R of 0.1 is "significant" without being meaningful (the reference paper reported R = 0.34).
#
# Figures: figures/Fig3_scRNA/Fig3{H,I,J}_p53_partner_*_<tissue>.pdf (the arrest readout downstream of the
#          anchored AMPK axis in script 23); the per-cell-type score plot goes to supplementary
# Tables:  results/pathway_selection/18_axis_gene_stats_<tissue>.csv
#          results/pathway_selection/18_axis_score_<tissue>.csv
#          results/pathway_selection/18_axis_correlation_<tissue>.csv
#          results/pathway_selection/18_axis_background_<tissue>.csv

suppressPackageStartupMessages({
  library(Seurat); library(Nebulosa); library(UCell); library(data.table)
  library(ggplot2); library(patchwork); library(Matrix)
})
set.seed(20260916)
fig_main <- "results/supplementary_figures/pathway_selection"; out_dir <- "results/pathway_selection"
  # draft panels; figures/ is owned by the 27-32 publication figure scripts (Fig 3 is now
  # Fig3_composite_<tissue> from 27 and Fig3_axis_<tissue> from 32). Writing here from an analysis
  # script silently overwrites them, because Windows filenames are case-insensitive.

dir.create(fig_main, recursive = TRUE, showWarnings = FALSE)
cols <- c(Control = "#3B7DD8", STZ = "#C8322F")

pair       <- c("Cdkn1a", "Ccnd1")                                     # CDKN1A - CCND1
# Curated p53 arrest / target set. Provenance matters here and is easy to lose: of these eight genes
# only CCND1 is a Fig 1g candidate. CDKN1A is a bulk DEG but sits in the turquoise module, which fails
# the M11b replication rule, and it is not in the hyperglycaemia gene set - so it never reaches the
# candidate intersection. This axis is secondary and data-driven (decisions R17, R23), NOT a Fig 1
# result, and the fig1_candidate column written below is what keeps that visible downstream.
axis_genes <- c("Cdkn1a", "Ccnd1", "Phlda3", "Zmat3", "Bax", "Serpine1", "Mdm2", "Trp53")
human_of   <- c(Cdkn1a = "CDKN1A", Ccnd1 = "CCND1", Phlda3 = "PHLDA3", Zmat3 = "ZMAT3",
                Bax = "BAX", Serpine1 = "SERPINE1", Mdm2 = "MDM2", Trp53 = "TP53")
fig1_candidates <- fread("results/bulk/intersect_genes.csv")$gene      # never typed in; see script 23
n_bg <- 2000L                                                          # background gene pairs per cell type
fmt_p <- function(p) ifelse(is.na(p), "NA", formatC(p, format = "g", digits = 2))

deg <- fread("results/bulk/DEG_T2D_vs_Control_all.csv")
hum <- deg[gene %in% human_of, .(human = gene, human_logFC = logFC, human_P = P.Value, human_FDR = adj.P.Val)]

for (tt in c("liver", "kidney")) {
  so <- readRDS(file.path("results/sc", paste0(tt, "_annotated.rds")))
  pb <- fread(file.path("results/sc", paste0(tt, "_pseudobulk_DE_all.csv")))
  md <- as.data.table(so@meta.data)
  genes <- intersect(axis_genes, rownames(so))

  eligible <- md[, .(n_ctl = sum(condition == "Control"), n_stz = sum(condition == "STZ")), by = cell_type][
    n_ctl >= 20 & n_stz >= 20]
  pb_ax <- pb[gene %in% genes & cell_type %in% eligible$cell_type]
  pb_ax[, human := human_of[gene]]
  pb_ax <- merge(pb_ax, hum, by = "human")
  pb_ax[, concordant := PValue < 0.05 & sign(logFC) == sign(human_logFC)]
  key <- pb_ax[, .(n = sum(concordant)), by = cell_type][order(-n)][1, cell_type]
  message(tt, ": key cell type = ", key, " (", pb_ax[cell_type == key & concordant == TRUE, .N], " axis genes concordant)")

  # ---- Fig 3D: axis genes in the key cell type ----------------------------------------
  sub <- subset(so, subset = cell_type == key)
  Idents(sub) <- "condition"
  fm <- FindMarkers(sub, ident.1 = "STZ", ident.2 = "Control", features = genes, test.use = "wilcox",
                    logfc.threshold = 0, min.pct = 0, verbose = FALSE)
  fm <- as.data.table(fm, keep.rownames = "gene")[, .(gene, cell_log2FC = avg_log2FC, cell_p = p_val,
                                                      pct_STZ = pct.1, pct_Control = pct.2)]
  fm[, cell_BH := p.adjust(cell_p, "BH")]
  st <- merge(fm, pb_ax[cell_type == key, .(gene, pb_logFC = logFC, pb_P = PValue, pb_FDR = FDR,
                                            human, human_logFC, human_P, human_FDR, concordant)], by = "gene")
  st[, fig1_candidate := human %in% fig1_candidates]
  fwrite(st[order(pb_P)], file.path(out_dir, paste0("18_axis_gene_stats_", tt, ".csv")))

  ex <- as.data.table(FetchData(sub, vars = c(genes, "condition", "mouse_id"), layer = "data"))
  long <- melt(ex, id.vars = c("condition", "mouse_id"), variable.name = "gene", value.name = "expr")
  pm <- long[, .(expr = mean(expr), n = .N), by = .(gene, condition, mouse_id)][n >= 10]
  st[, lab := sprintf("%s (human %s)\ncell BH %s | pseudobulk P %s", gene,
                      ifelse(human_logFC > 0, "up", "down"), fmt_p(cell_BH), fmt_p(pb_P))]
  long <- merge(long, st[, .(gene, lab)], by = "gene"); pm <- merge(pm, st[, .(gene, lab)], by = "gene")
  p3d <- ggplot(long, aes(condition, expr, fill = condition)) +
    geom_violin(scale = "width", trim = TRUE, colour = NA, alpha = 0.7) +
    geom_point(data = pm, aes(condition, expr), shape = 21, fill = "white", size = 2) +
    facet_wrap(~ lab, scales = "free_y", nrow = 2) +
    scale_fill_manual(values = cols, guide = "none") +
    labs(x = NULL, y = "normalized expression",
         title = paste0("p53 arrest axis in ", key, " (", tt, "); points = per-mouse means")) +
    theme_bw(base_size = 9)
  ggsave(file.path(fig_main, paste0("Fig3H_p53_partner_expression_", tt, ".pdf")), p3d, width = 11, height = 5.5)

  # ---- Fig 3E: p53 target score per cell type, against matched random sets ------------
  sig <- list(p53_targets = intersect(c("Cdkn1a", "Phlda3", "Zmat3", "Bax", "Ccng1", "Mdm2", "Aen", "Serpine1",
                                        "Trp53inp1", "Eda2r", "Ccnd1"), rownames(so)))
  so <- AddModuleScore_UCell(so, features = sig, name = "")
  md2 <- as.data.table(so@meta.data)[, .(cell_type, condition, mouse_id, score = p53_targets)]
  per_mouse <- md2[, .(score = mean(score), n = .N), by = .(cell_type, condition, mouse_id)][n >= 10]
  tst <- per_mouse[, {
    a <- score[condition == "STZ"]; b <- score[condition == "Control"]
    if (length(a) >= 2 && length(b) >= 2) {
      t <- t.test(a, b); .(delta = mean(a) - mean(b), p = t$p.value, n_stz = length(a), n_ctl = length(b))
    } else .(delta = NA_real_, p = NA_real_, n_stz = length(a), n_ctl = length(b))
  }, by = cell_type][order(p)]
  tst[, BH := p.adjust(p, "BH")]
  fwrite(tst, file.path(out_dir, paste0("18_axis_score_", tt, ".csv")))
  per_mouse <- merge(per_mouse, tst[, .(cell_type, p)], by = "cell_type")
  per_mouse[, lab := paste0(cell_type, "\np = ", fmt_p(p))]
  p3e <- ggplot(per_mouse, aes(condition, score, fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) + geom_point(shape = 21, size = 2, fill = "white") +
    facet_wrap(~ lab, scales = "free_y") + scale_fill_manual(values = cols, guide = "none") +
    labs(x = NULL, y = "UCell p53-target score", title = paste0("p53 target program per cell type (", tt, ")")) +
    theme_bw(base_size = 8)
  ggsave(file.path("results/supplementary_figures/pathway_selection/p53_axis", paste0("Fig3E_p53_score_", tt, ".pdf")), p3e, width = 11, height = 7)

  # ---- Fig 3F: co-localization ---------------------------------------------------------
  if (all(pair %in% rownames(so))) {
    pf <- plot_density(so, features = pair, joint = TRUE, reduction = "umap")
    co <- as.data.table(FetchData(so, vars = pair, layer = "data"))
    co[, cell_type := md$cell_type]
    cop <- co[, .(pct = 100 * mean(get(pair[1]) > 0 & get(pair[2]) > 0), n = .N), by = cell_type][order(-pct)]
    pbar <- ggplot(cop, aes(pct, factor(cell_type, levels = rev(cell_type)))) +
      geom_col(fill = "#8E5CC8") +
      labs(x = paste0("% cells co-expressing ", pair[1], " + ", pair[2]), y = NULL) + theme_bw(base_size = 9)
    ggsave(file.path(fig_main, paste0("Fig3I_p53_partner_colocalization_", tt, ".pdf")), pf | pbar, width = 17, height = 4.5)
  }

  # ---- Fig 3G: correlation, with a detection-matched background null ------------------
  e <- as.data.table(FetchData(so, vars = pair, layer = "data"))
  e[, `:=`(cell_type = md$cell_type, condition = md$condition, mouse_id = md$mouse_id)]
  ek <- e[cell_type == key]
  ct_test <- suppressWarnings(cor.test(ek[[pair[1]]], ek[[pair[2]]], method = "pearson"))
  cocells <- which(ek[[pair[1]]] > 0 & ek[[pair[2]]] > 0)
  sp <- if (length(cocells) >= 10) suppressWarnings(cor.test(ek[[pair[1]]][cocells], ek[[pair[2]]][cocells],
                                                            method = "spearman")) else NULL
  pmk <- ek[, .(a = mean(get(pair[1])), b = mean(get(pair[2])), n = .N), by = .(mouse_id, condition)][n >= 10]
  pmt <- if (nrow(pmk) >= 4) cor.test(pmk$a, pmk$b, method = "pearson") else NULL

  # background: random expressed pairs in the same cell type, matched on detection rate (+/- 5 percentage points)
  cells_k <- which(md$cell_type == key)
  cnt <- GetAssayData(so, layer = "data")[, cells_k, drop = FALSE]
  det <- Matrix::rowMeans(cnt > 0)
  det <- det[det >= 0.02 & det <= 0.9]
  d1 <- mean(ek[[pair[1]]] > 0); d2 <- mean(ek[[pair[2]]] > 0)
  pool1 <- names(det)[abs(det - d1) <= 0.05]; pool2 <- names(det)[abs(det - d2) <= 0.05]
  pool1 <- setdiff(pool1, pair); pool2 <- setdiff(pool2, pair)
  bg <- if (length(pool1) >= 10 && length(pool2) >= 10) {
    g1 <- sample(pool1, n_bg, replace = TRUE); g2 <- sample(pool2, n_bg, replace = TRUE)
    keep <- g1 != g2; g1 <- g1[keep]; g2 <- g2[keep]
    m <- as.matrix(cnt[unique(c(g1, g2)), , drop = FALSE])
    vapply(seq_along(g1), function(i) suppressWarnings(cor(m[g1[i], ], m[g2[i], ])), 0)
  } else numeric(0)
  bg <- bg[is.finite(bg)]
  emp_p <- if (length(bg)) mean(bg >= unname(ct_test$estimate)) else NA_real_
  corr <- data.table(tissue = tt, cell_type = key, pair = paste(pair, collapse = "-"),
                     n_cells = nrow(ek), pct_coexpr = 100 * length(cocells) / nrow(ek),
                     cell_R = unname(ct_test$estimate), cell_p = ct_test$p.value,
                     coexpr_rho = if (is.null(sp)) NA_real_ else unname(sp$estimate),
                     coexpr_p = if (is.null(sp)) NA_real_ else sp$p.value,
                     pb_r = if (is.null(pmt)) NA_real_ else unname(pmt$estimate),
                     pb_p = if (is.null(pmt)) NA_real_ else pmt$p.value, n_mice = nrow(pmk),
                     bg_n = length(bg), bg_mean_R = mean(bg), bg_q95 = if (length(bg)) quantile(bg, 0.95) else NA_real_,
                     empirical_p = emp_p)
  fwrite(corr, file.path(out_dir, paste0("18_axis_correlation_", tt, ".csv")))
  if (length(bg)) fwrite(data.table(R = bg), file.path(out_dir, paste0("18_axis_background_", tt, ".csv")))
  print(corr)

  g1p <- ggplot(ek, aes(get(pair[1]), get(pair[2]))) +
    geom_jitter(width = 0.02, height = 0.02, alpha = 0.25, size = 0.6, colour = "grey35") +
    geom_smooth(method = "lm", formula = y ~ x, colour = "#C8322F", se = TRUE) +
    labs(x = pair[1], y = pair[2],
         title = sprintf("%s, cell level: R = %.2f (p %s)", key, corr$cell_R, fmt_p(corr$cell_p)),
         subtitle = sprintf("background R (matched random pairs): mean %.2f, 95%% %.2f; empirical p = %s",
                            corr$bg_mean_R, corr$bg_q95, fmt_p(corr$empirical_p))) +
    theme_bw(base_size = 9)
  g2p <- ggplot(pmk, aes(a, b, fill = condition)) + geom_point(shape = 21, size = 3) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "grey30", se = FALSE) +
    scale_fill_manual(values = cols) +
    labs(x = paste0("mean ", pair[1]), y = paste0("mean ", pair[2]),
         title = sprintf("per-mouse: r = %.2f (p %s), n = %d", corr$pb_r, fmt_p(corr$pb_p), corr$n_mice)) +
    theme_bw(base_size = 9)
  g3p <- if (length(bg)) ggplot(data.table(R = bg), aes(R)) +
    geom_histogram(bins = 60, fill = "grey75", colour = NA) +
    geom_vline(xintercept = corr$cell_R, colour = "#C8322F", linewidth = 1) +
    labs(x = "cell-level Pearson R", y = "matched random pairs",
         title = "axis pair vs background") + theme_bw(base_size = 9) else plot_spacer()
  ggsave(file.path(fig_main, paste0("Fig3J_p53_partner_correlation_", tt, ".pdf")), g1p | g2p | g3p, width = 15, height = 4.6)
  rm(so); gc(verbose = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_18.txt"))
