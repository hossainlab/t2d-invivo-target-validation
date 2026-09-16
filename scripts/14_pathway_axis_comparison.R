# 14_pathway_axis_comparison.R — pathway selection following the reference paper's logic
# (Xu et al. 2026: top KEGG pathway of the intersect genes -> pathview-activated axis -> the axis genes are
#  changed, co-localized and correlated in the key cell type in scRNA -> downstream readouts from literature).
#
# Compares two candidate axes on identical criteria in mouse STZ vs Control scRNA (liver, kidney):
#   AMPK signaling (hsa04152; top KEGG of the 19 intersect genes)
#   p53 signaling  (hsa04115; CDKN1A/CCND1 axis)
# For each tissue x cell type (pseudobulk-eligible types):
#   1. per-gene STZ change (pseudobulk, script 09) and direction vs human T2D liver (script 02)
#   2. pairwise correlation within each axis: per-mouse pseudobulk log-CPM Pearson (n = 7) and
#      co-expressing-cell Spearman; co-expression fraction (co-localization metric)
#   3. axis score per cell type: concordant genes (P < 0.05), best-correlated concordant pair
# Also draws the p53 pathview map for comparison with the AMPK map (Fig 2D).
#
# Outputs: results/pathway_selection/<tissue>_axis_gene_evidence.csv, <tissue>_axis_pair_correlation.csv,
#          axis_summary.csv
# Figures: results/supplementary_figures/pathway_selection/{coloc,corr}_<axis>_<tissue>.pdf, pathview_p53/

suppressPackageStartupMessages({
  library(Seurat)
  library(Nebulosa)
  library(edgeR)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(org.Hs.eg.db)
  library(pathview)
})
set.seed(20260914)
out_dir <- normalizePath("results/pathway_selection", mustWork = FALSE)
fig_dir <- normalizePath("results/supplementary_figures/pathway_selection", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Axis genes: intersect/pathway genes plus the core signalling nodes used as protein readouts
axes <- list(
  # AMPK: intersect genes (Irs2, Igf1, Ppargc1a, Fbp1, Ccnd1) + core nodes + nodes coloured in the pathview map
  # (Fig 2D: SIRT1/PGC-1a/FOXO down; SREBP-1c/FAS up)
  AMPK = c("Prkaa1", "Prkaa2", "Stk11", "Camkk2", "Sirt1", "Ppargc1a", "Foxo1", "Irs2", "Igf1", "Fbp1",
           "Srebf1", "Fasn", "Ccnd1"),
  p53  = c("Trp53", "Mdm2", "Cdkn1a", "Ccnd1", "Bax", "Zmat3", "Phlda3", "Serpine1", "Igf1")
)
human_of <- function(m) ifelse(m == "Trp53", "TP53", toupper(m))
bulk <- fread("results/bulk/DEG_T2D_vs_Control_all.csv")
hdir <- setNames(bulk$logFC, bulk$gene); hP <- setNames(bulk$P.Value, bulk$gene)

summary_rows <- list()
for (tissue in c("liver", "kidney")) {
  so <- readRDS(file.path("results/sc", paste0(tissue, "_annotated.rds")))
  pb <- fread(file.path("results/sc", paste0(tissue, "_pseudobulk_DE_all.csv")))
  ex <- fread(file.path("results/sc", paste0(tissue, "_pseudobulk_DE_all_excl_STZ_4.csv")))
  cts <- unique(pb$cell_type)
  counts <- LayerData(so, assay = "RNA", layer = "counts"); dat <- LayerData(so, assay = "RNA", layer = "data")
  md <- as.data.table(so@meta.data, keep.rownames = "cell")

  # ---- 1. gene-level evidence --------------------------------------------------------
  all_genes <- intersect(unique(unlist(axes)), rownames(so))
  gene_ev <- rbindlist(lapply(names(axes), function(ax) {
    g <- intersect(axes[[ax]], rownames(so))
    x <- merge(data.table(axis = ax, gene = g), pb[gene %in% g, .(gene, cell_type, logFC, PValue, FDR, pct_cells_detected)],
               by = "gene", all.x = TRUE, allow.cartesian = TRUE)
    x <- merge(x, ex[gene %in% g, .(gene, cell_type, logFC_exSTZ4 = logFC, P_exSTZ4 = PValue)], by = c("gene", "cell_type"), all.x = TRUE)
    x[, human := human_of(gene)]
    x[, `:=`(human_logFC = hdir[human], human_P = hP[human])]
    x[, concordant := !is.na(logFC) & !is.na(human_logFC) & sign(logFC) == sign(human_logFC)]
    x[, concordant_P05 := concordant & PValue < 0.05]
    x[, robust_exSTZ4 := concordant_P05 & !is.na(P_exSTZ4) & P_exSTZ4 < 0.05 & sign(logFC_exSTZ4) == sign(logFC)]
    x
  }))
  fwrite(gene_ev, file.path(out_dir, paste0(tissue, "_axis_gene_evidence.csv")))

  # ---- 2. pairwise correlation within each axis per cell type ------------------------
  pair_rows <- list()
  for (ct in cts) {
    cells <- md[cell_type == ct, cell]
    grp <- factor(md[cell_type == ct, mouse_id])
    mm <- sparse.model.matrix(~0 + grp); colnames(mm) <- levels(grp)
    keep_m <- colnames(mm)[colSums(mm) >= 20]
    pbm <- as.matrix(counts[all_genes, cells, drop = FALSE] %*% mm)[, keep_m, drop = FALSE]
    libs <- colSums(counts[, cells, drop = FALSE] %*% mm)[keep_m]
    lcpm <- log2(t(t(pbm) / libs) * 1e6 + 1)
    for (ax in names(axes)) {
      g <- intersect(axes[[ax]], all_genes)
      prs <- t(combn(g, 2))
      for (i in seq_len(nrow(prs))) {
        g1 <- prs[i, 1]; g2 <- prs[i, 2]
        pr <- if (ncol(lcpm) >= 5 && sd(lcpm[g1, ]) > 0 && sd(lcpm[g2, ]) > 0) cor.test(lcpm[g1, ], lcpm[g2, ]) else NULL
        both <- dat[g1, cells] > 0 & dat[g2, cells] > 0
        sp <- if (sum(both) >= 10) suppressWarnings(cor.test(dat[g1, cells][both], dat[g2, cells][both], method = "spearman", exact = FALSE)) else NULL
        pair_rows[[length(pair_rows) + 1]] <- data.table(
          tissue, cell_type = ct, axis = ax, gene1 = g1, gene2 = g2, n_mice = ncol(lcpm),
          pseudobulk_r = if (is.null(pr)) NA_real_ else unname(pr$estimate), pseudobulk_p = if (is.null(pr)) NA_real_ else pr$p.value,
          cells_coexpressing = sum(both), pct_coexpressing = round(100 * mean(both), 2),
          cell_rho = if (is.null(sp)) NA_real_ else unname(sp$estimate), cell_p = if (is.null(sp)) NA_real_ else sp$p.value)
      }
    }
  }
  pairs <- rbindlist(pair_rows)
  conc <- gene_ev[concordant_P05 == TRUE, .(gene, cell_type, axis)]
  pairs[, both_concordant := mapply(function(a, b, c, x) nrow(conc[cell_type == c & axis == x & gene %in% c(a, b)]) == 2,
                                    gene1, gene2, cell_type, axis)]
  fwrite(pairs[order(axis, cell_type, pseudobulk_p)], file.path(out_dir, paste0(tissue, "_axis_pair_correlation.csv")))

  # ---- 3. axis summary per cell type --------------------------------------------------
  for (ax in names(axes)) for (ct in cts) {
    ge <- gene_ev[axis == ax & cell_type == ct]
    pc <- pairs[axis == ax & cell_type == ct & both_concordant == TRUE][order(pseudobulk_p)]
    summary_rows[[length(summary_rows) + 1]] <- data.table(
      tissue, axis = ax, cell_type = ct, genes_tested = ge[!is.na(PValue), .N],
      concordant_P05 = ge[concordant_P05 == TRUE, .N], concordant_FDR05 = ge[concordant_P05 == TRUE & FDR < 0.05, .N],
      robust_exSTZ4 = ge[robust_exSTZ4 == TRUE, .N],
      concordant_genes = paste(ge[concordant_P05 == TRUE][order(PValue), gene], collapse = ";"),
      best_concordant_pair = if (nrow(pc)) paste(pc$gene1[1], pc$gene2[1], sep = "-") else NA_character_,
      best_pair_pseudobulk_r = if (nrow(pc)) pc$pseudobulk_r[1] else NA_real_,
      best_pair_pseudobulk_p = if (nrow(pc)) pc$pseudobulk_p[1] else NA_real_,
      best_pair_cell_rho = if (nrow(pc)) pc$cell_rho[1] else NA_real_,
      best_pair_cell_p = if (nrow(pc)) pc$cell_p[1] else NA_real_,
      best_pair_pct_coexpr = if (nrow(pc)) pc$pct_coexpressing[1] else NA_real_)
  }

  # ---- 4. figures: best concordant pair per axis (tissue-wide best cell type) ----------
  sm <- rbindlist(summary_rows)[tissue == get("tissue", envir = environment())]
  for (ax in names(axes)) {
    best <- sm[axis == ax & !is.na(best_concordant_pair)][order(-concordant_P05, best_pair_pseudobulk_p)][1]
    if (!nrow(best) || is.na(best$best_concordant_pair)) { message(tissue, " ", ax, ": no concordant pair"); next }
    gp <- strsplit(best$best_concordant_pair, "-")[[1]]; ct <- best$cell_type
    pden <- tryCatch(plot_density(so, features = gp, joint = TRUE, reduction = "umap"),
                     error = function(e) FeaturePlot(so, features = gp, blend = TRUE, raster = TRUE))
    pct_by_ct <- md[, .(pct_coexpr = 100 * mean(dat[gp[1], cell] > 0 & dat[gp[2], cell] > 0)), by = cell_type][order(-pct_coexpr)]
    pbar <- ggplot(pct_by_ct, aes(reorder(cell_type, pct_coexpr), pct_coexpr, fill = cell_type == ct)) + geom_col() + coord_flip() +
      scale_fill_manual(values = c(`TRUE` = "#C8322F", `FALSE` = "grey70"), guide = "none") +
      labs(x = NULL, y = paste0("% cells co-expressing ", gp[1], " & ", gp[2]), title = paste(tissue, ax, "axis")) + theme_bw()
    ggsave(file.path(fig_dir, paste0("coloc_", ax, "_", tissue, ".pdf")), wrap_plots(pden) | pbar, width = 16, height = 5)

    cells <- md[cell_type == ct, cell]; grp <- factor(md[cell_type == ct, mouse_id])
    mm <- sparse.model.matrix(~0 + grp); colnames(mm) <- levels(grp)
    pbm <- as.matrix(counts[gp, cells, drop = FALSE] %*% mm); libs <- colSums(counts[, cells, drop = FALSE] %*% mm)
    lc <- log2(t(t(pbm) / libs) * 1e6 + 1)
    df <- data.table(mouse = colnames(lc), x = lc[gp[1], ], y = lc[gp[2], ])
    df[, condition := sub("_[0-9]+$", "", mouse)]
    pcor <- ggplot(df, aes(x, y)) + geom_smooth(method = "lm", colour = "black", linewidth = 0.6) +
      geom_point(aes(colour = condition), size = 3) + scale_colour_manual(values = c(Control = "#3B7DD8", STZ = "#C8322F")) +
      labs(x = paste(gp[1], "log-CPM"), y = paste(gp[2], "log-CPM"), title = paste0(tissue, " ", ct, ": ", ax, " axis"),
           subtitle = sprintf("Pearson r = %.2f, p = %.2g (n = %d mice)", best$best_pair_pseudobulk_r, best$best_pair_pseudobulk_p, ncol(lc))) +
      theme_bw()
    ggsave(file.path(fig_dir, paste0("corr_", ax, "_", tissue, ".pdf")), pcor, width = 6, height = 5)
  }
  rm(so, counts, dat); gc()
}
summ <- rbindlist(summary_rows)
setorder(summ, tissue, axis, -concordant_P05, best_pair_pseudobulk_p)
fwrite(summ, file.path(out_dir, "axis_summary.csv"))
print(summ[concordant_P05 > 0])

# ---- 5. p53 pathview map for comparison with Fig 2D (AMPK) ----------------------------
fc <- bulk[, .(gene, logFC)][, entrez := suppressMessages(mapIds(org.Hs.eg.db, gene, "ENTREZID", "SYMBOL"))][!is.na(entrez)]
pv <- file.path(fig_dir, "pathview_p53"); dir.create(pv, showWarnings = FALSE)
old <- getwd()
tryCatch({
  setwd(pv)
  pathview(gene.data = setNames(fc$logFC, fc$entrez), pathway.id = "04115", species = "hsa",
           out.suffix = "T2D_vs_Control_logFC", limit = list(gene = 1.5), kegg.native = TRUE)
}, error = function(e) message("pathview p53 failed: ", conditionMessage(e)), finally = setwd(old))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_14.txt"))
