# 13_sc_colocalization.R — Fig 3F-G analogue of the reference paper (co-localization + correlation of the
# candidate genes in the key cell type; paper: Rac1/Pak1 in neutrophils)
# Usage: Rscript scripts/13_sc_colocalization.R liver|kidney
# Run after 10_integration_bulk_sc.R.
#
# Steps:
#   1. Candidates = bulk main set (results/bulk/intersect_genes.csv) -> 1:1 mouse orthologs present in the data.
#   2. Focus cell type = the cell type with the most candidates changed in STZ in the bulk (T2D) direction
#      (pseudobulk P < 0.05); ties broken by cell number.
#   3. Fig 3F: joint kernel-density UMAP (Nebulosa) of the candidates across the atlas.
#   4. Fig 3G: correlation of candidate pairs in the focus cell type:
#        primary   = Pearson on per-mouse pseudobulk log-CPM (n = 7 mice; avoids dropout-driven cell-level r)
#        secondary = Spearman on cells expressing both genes
#
# Outputs: results/integration/<tissue>_candidate_correlation.csv, <tissue>_focus_celltype.csv
# Figures: figures/Fig3_scRNA/Fig3F_colocalization_<tissue>.pdf, Fig3G_correlation_<tissue>.pdf

suppressPackageStartupMessages({
  library(Seurat)
  library(Nebulosa)
  library(babelgene)
  library(edgeR)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})
set.seed(20260914)
tissue <- tolower(commandArgs(trailingOnly = TRUE)[1])
if (is.na(tissue) || !tissue %in% c("liver", "kidney"))
  stop("usage: Rscript scripts/13_sc_colocalization.R <liver|kidney>\n",
       if (is.na(tissue)) "  no tissue argument was given."
       else paste0("  got '", tissue, "'."), call. = FALSE)
out_dir <- "results/integration"; fig_main <- "results/supplementary_figures/candidates_19gene"  # main Fig 3F-G now from 15b (AMPK axis)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE); dir.create(fig_main, recursive = TRUE, showWarnings = FALSE)

ints <- fread("results/bulk/intersect_genes.csv")
so   <- readRDS(file.path("results/sc", paste0(tissue, "_annotated.rds")))
pb   <- fread(file.path("results/sc", paste0(tissue, "_pseudobulk_DE_all.csv")))

# ---------------- 1. Candidate orthologs ----------------------------------------------
o <- as.data.table(orthologs(genes = ints$gene, species = "mouse", human = TRUE))
o <- unique(o[, .(human = human_symbol, mouse = symbol)])[, n_h := .N, by = human][, n_m := .N, by = mouse][n_h == 1 & n_m == 1]
cand <- merge(ints[, .(human = gene, logFC_human = logFC)], o[, .(human, mouse)], by = "human")[mouse %in% rownames(so)]
message(tissue, ": candidates with 1:1 mouse ortholog in data: ", paste(cand$mouse, collapse = ", "),
        " | without: ", paste(setdiff(ints$gene, cand$human), collapse = ", "))
if (nrow(cand) < 1) stop("No candidate ortholog detected in ", tissue)

# ---------------- 2. Focus cell type ---------------------------------------------------
ev <- merge(cand, pb[, .(mouse = gene, cell_type, logFC, PValue, FDR)], by = "mouse")
ev[, concordant_nominal := PValue < 0.05 & sign(logFC) == sign(logFC_human)]
n_cells <- as.data.table(table(cell_type = so$cell_type))
focus_tab <- merge(ev[, .(n_concordant = sum(concordant_nominal), genes = paste(mouse[concordant_nominal], collapse = ";")),
                      by = cell_type], setnames(n_cells, "N", "n_cells"), by = "cell_type")
setorder(focus_tab, -n_concordant, -n_cells)
fwrite(focus_tab, file.path(out_dir, paste0(tissue, "_focus_celltype.csv")))
focus <- focus_tab$cell_type[1]
message("Focus cell type: ", focus, " (", focus_tab$n_concordant[1], " concordant candidates: ", focus_tab$genes[1], ")")

# ---------------- 3. Fig 3F: joint density ------------------------------------------
# Genes to plot: candidates changed in the focus cell type in the bulk direction (strongest first),
# then remaining candidates by bulk P value; at most 4
conc_focus <- ev[cell_type == focus & concordant_nominal == TRUE][order(PValue), mouse]
bulk_order <- merge(cand, ints[, .(human = gene, P_bulk = P.Value)], by = "human")[order(P_bulk), mouse]
genes_plot <- head(unique(c(conc_focus, bulk_order)), 4)
message("Genes plotted: ", paste(genes_plot, collapse = ", "))
pf <- tryCatch({
  p <- plot_density(so, features = genes_plot, joint = length(genes_plot) > 1, reduction = "umap")
  if (inherits(p, "patchwork")) p else wrap_plots(p)
}, error = function(e) {
  message("Nebulosa failed (", conditionMessage(e), "); using FeaturePlot")
  FeaturePlot(so, features = genes_plot, order = TRUE, raster = TRUE)
})
p_ct <- DimPlot(so, group.by = "cell_type", label = TRUE, repel = TRUE, raster = TRUE) + NoLegend() +
  ggtitle(paste(tissue, "- focus:", focus))
ggsave(file.path(fig_main, paste0("Fig3F_colocalization_", tissue, ".pdf")), p_ct | pf,
       width = 6 + 4 * min(length(genes_plot) + 1, 3), height = 5 * ceiling((length(genes_plot) + 1) / 3))

# ---------------- 4. Fig 3G: correlation in focus cell type --------------------------
sub <- subset(so, subset = cell_type == focus)
cnt <- LayerData(sub, assay = "RNA", layer = "counts")
grp <- factor(sub$mouse_id); mm <- sparse.model.matrix(~0 + grp); colnames(mm) <- levels(grp)
pbm <- as.matrix(cnt %*% mm)
lcpm <- cpm(DGEList(pbm), log = TRUE, prior.count = 1)
cond <- setNames(as.character(sub$condition[match(colnames(pbm), sub$mouse_id)]), colnames(pbm))
dat  <- LayerData(sub, assay = "RNA", layer = "data")

pairs <- if (length(genes_plot) >= 2) t(combn(genes_plot, 2)) else matrix(character(0), ncol = 2)
cor_tab <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  g1 <- pairs[i, 1]; g2 <- pairs[i, 2]
  pr <- cor.test(lcpm[g1, ], lcpm[g2, ], method = "pearson")
  both <- dat[g1, ] > 0 & dat[g2, ] > 0
  sp <- if (sum(both) >= 10) cor.test(dat[g1, both], dat[g2, both], method = "spearman", exact = FALSE) else NULL
  data.table(tissue, cell_type = focus, gene1 = g1, gene2 = g2,
             pseudobulk_pearson_r = unname(pr$estimate), pseudobulk_p = pr$p.value, n_mice = ncol(lcpm),
             cells_coexpressing = sum(both), pct_coexpressing = round(100 * mean(both), 2),
             cell_spearman_rho = if (is.null(sp)) NA_real_ else unname(sp$estimate),
             cell_spearman_p = if (is.null(sp)) NA_real_ else sp$p.value)
}))
fwrite(cor_tab, file.path(out_dir, paste0(tissue, "_candidate_correlation.csv")))
print(cor_tab)

if (nrow(cor_tab)) {
  top <- cor_tab[order(pseudobulk_p)][1]
  df <- data.table(mouse_id = colnames(lcpm), x = lcpm[top$gene1, ], y = lcpm[top$gene2, ], condition = cond[colnames(lcpm)])
  pg1 <- ggplot(df, aes(x, y)) + geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.6) +
    geom_point(aes(colour = condition), size = 3) +
    scale_colour_manual(values = c(Control = "#3B7DD8", STZ = "#C8322F")) +
    labs(x = paste(top$gene1, "log-CPM"), y = paste(top$gene2, "log-CPM"),
         title = paste0(focus, ": per-mouse pseudobulk"),
         subtitle = sprintf("Pearson r = %.2f, p = %.2g (n = %d mice)", top$pseudobulk_pearson_r, top$pseudobulk_p, top$n_mice)) +
    theme_bw(base_size = 12)
  both <- dat[top$gene1, ] > 0 & dat[top$gene2, ] > 0
  dc <- data.table(x = dat[top$gene1, both], y = dat[top$gene2, both], condition = sub$condition[both])
  pg2 <- ggplot(dc, aes(x, y, colour = condition)) + geom_point(size = 0.6, alpha = 0.5) +
    geom_smooth(aes(group = 1), method = "lm", colour = "black", linewidth = 0.6) +
    scale_colour_manual(values = c(Control = "#3B7DD8", STZ = "#C8322F")) +
    labs(x = top$gene1, y = top$gene2, title = "Co-expressing cells",
         subtitle = if (is.na(top$cell_spearman_rho)) "< 10 co-expressing cells"
                    else sprintf("Spearman rho = %.2f, p = %.2g (%d cells)", top$cell_spearman_rho, top$cell_spearman_p, top$cells_coexpressing)) +
    theme_bw(base_size = 12)
  ggsave(file.path(fig_main, paste0("Fig3G_correlation_", tissue, ".pdf")), pg1 | pg2, width = 11, height = 5)
} else message("Fewer than 2 candidate orthologs: correlation panel (Fig3G) skipped")

writeLines(capture.output(sessionInfo()), file.path(out_dir, paste0("sessionInfo_13_", tissue, ".txt")))
