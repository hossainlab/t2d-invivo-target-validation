# 10_integration_bulk_sc.R — Phases 7.4-7.5 and 8 of docs/analysis_plan.md
# Usage: Rscript scripts/10_integration_bulk_sc.R liver|kidney
# Human liver bulk intersect genes -> mouse 1:1 orthologs -> (a) UCell signature score by cell type,
# (b) per-gene concordance with STZ pseudobulk DE, (c) cell-type specificity (tau), (d) tiers.
#
# Inputs: results/bulk/intersect_genes.csv, results/bulk/DEG_T2D_vs_Control_all.csv,
#         results/sc/<tissue>_annotated.rds, results/sc/<tissue>_pseudobulk_DE_all.csv
# Outputs: results/integration/<tissue>_orthologs.csv, <tissue>_signature_score_test.csv,
#          <tissue>_candidate_concordance.csv, <tissue>_logFC_correlation.csv
# Figures: figures/integration/10_<tissue>_signature_score.pdf, 10_<tissue>_concordance_heatmap.pdf,
#          10_<tissue>_candidates_dotplot.pdf

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(babelgene)
  library(data.table)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(Matrix)
})
set.seed(20260914)
tissue <- tolower(commandArgs(trailingOnly = TRUE)[1])
bulk_dir <- "results/bulk"; sc_dir <- "results/sc"
out_dir <- "results/integration"; fig_dir <- "results/supplementary_figures/integration"
fig_main <- "results/supplementary_figures/candidates_19gene"  # 19-gene candidate panels (Fig 3D-G main panels now from 15b, AMPK axis)
dir.create(fig_main, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

ints <- fread(file.path(bulk_dir, "intersect_genes.csv"))
deg  <- fread(file.path(bulk_dir, "DEG_T2D_vs_Control_all.csv"))
so   <- readRDS(file.path(sc_dir, paste0(tissue, "_annotated.rds")))
pb   <- fread(file.path(sc_dir, paste0(tissue, "_pseudobulk_DE_all.csv")))
val_file <- "results/validation/GSE23343_limma_all.csv"  # independent human liver (script 11)
val  <- if (file.exists(val_file)) fread(val_file)[, .(gene, logFC_GSE23343 = logFC, P_GSE23343 = P.Value)] else NULL

# ---------------- 1:1 orthologs --------------------------------------------------------
map_1to1 <- function(hs) {
  o <- as.data.table(orthologs(genes = hs, species = "mouse", human = TRUE))
  o <- o[, .(human = human_symbol, mouse = symbol, support_n)]
  o <- unique(o)
  o[, n_h := .N, by = human][, n_m := .N, by = mouse]
  o[n_h == 1 & n_m == 1, .(human, mouse, support_n)]
}
sig_col <- "deg_call"
orth_all <- map_1to1(unique(c(ints$gene, deg[get(sig_col) != "NS", gene])))
orth_all <- orth_all[mouse %in% rownames(so)]
fwrite(orth_all, file.path(out_dir, paste0(tissue, "_orthologs.csv")))
message("Intersect genes: ", nrow(ints), "; with 1:1 mouse ortholog present in data: ",
        sum(ints$gene %in% orth_all$human))

# ---------------- (a) Signature scores ----------------------------------------------
to_mouse <- function(h) orth_all[human %in% h, mouse]
sigs <- list(
  intersect_up   = to_mouse(ints[logFC > 0, gene]),
  intersect_down = to_mouse(ints[logFC < 0, gene]),
  T2D_up_top100   = to_mouse(deg[get(sig_col) == "Up"][order(P.Value)][1:min(.N, 100), gene]),
  T2D_down_top100 = to_mouse(deg[get(sig_col) == "Down"][order(P.Value)][1:min(.N, 100), gene])
)
sigs <- sigs[lengths(sigs) >= 5]
message("Signature sizes (>= 5 mouse orthologs kept): ",
        if (length(sigs)) paste(names(sigs), lengths(sigs), sep = "=", collapse = ", ") else "none")
if (length(sigs)) {
  so <- AddModuleScore_UCell(so, features = sigs, assay = "RNA", slot = "counts", ncores = 4, name = "_UCell")
}
score_cols <- paste0(names(sigs), "_UCell")

md <- as.data.table(so@meta.data, keep.rownames = "cell")
if (length(score_cols)) {
per_mouse <- md[, lapply(.SD, mean), by = .(cell_type, mouse_id, condition), .SDcols = score_cols]
fwrite(per_mouse, file.path(out_dir, paste0(tissue, "_signature_score_per_mouse.csv")))
sig_test <- rbindlist(lapply(score_cols, function(sc) {
  per_mouse[, {
    if (uniqueN(condition) == 2 && min(table(condition)) >= 2) {
      tt <- t.test(get(sc) ~ condition)
      .(signature = sc, mean_Control = mean(get(sc)[condition == "Control"]),
        mean_STZ = mean(get(sc)[condition == "STZ"]), p = tt$p.value, n = .N)
    } else NULL
  }, by = cell_type]
}))
sig_test[, diff := mean_STZ - mean_Control]
sig_test[, padj := p.adjust(p, "BH"), by = signature]
fwrite(sig_test[order(signature, p)], file.path(out_dir, paste0(tissue, "_signature_score_test.csv")))
print(sig_test[order(p)][1:min(.N, 20)])

plong <- melt(per_mouse, id.vars = c("cell_type", "mouse_id", "condition"), variable.name = "signature", value.name = "score")
ps <- ggplot(plong, aes(cell_type, score, fill = condition)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), width = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 1) +
  facet_wrap(~signature, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(Control = "#9DB9E0", STZ = "#E39A98")) +
  labs(x = NULL, y = "Mean UCell score per mouse", title = paste(tissue, "- human Ob-T2D liver signatures")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_main, paste0("Fig3E_signature_score_by_celltype_", tissue, ".pdf")), ps, width = 11, height = 3 * length(sigs))
} else message("No signature with >= 5 orthologs: signature scoring (Fig3E) skipped")

# ---------------- (b) Global logFC correlation human vs mouse ------------------------
cor_tab <- rbindlist(lapply(unique(pb$cell_type), function(ct) {
  x <- merge(orth_all, deg[, .(human = gene, logFC_h = logFC, sig_h = get(sig_col))], by = "human")
  x <- merge(x, pb[cell_type == ct, .(mouse = gene, logFC_m = logFC, FDR_m = FDR)], by = "mouse")
  xd <- x[sig_h != "NS"]
  data.table(cell_type = ct, n_DEG_tested = nrow(xd),
             spearman_DEGs = if (nrow(xd) > 10) cor(xd$logFC_h, xd$logFC_m, method = "spearman") else NA,
             sign_agree_DEGs = if (nrow(xd)) mean(sign(xd$logFC_h) == sign(xd$logFC_m)) else NA,
             sign_test_p = if (nrow(xd)) binom.test(sum(sign(xd$logFC_h) == sign(xd$logFC_m)), nrow(xd))$p.value else NA)
}))
fwrite(cor_tab[order(sign_test_p)], file.path(out_dir, paste0(tissue, "_logFC_correlation.csv")))
print(cor_tab[order(sign_test_p)])

# ---------------- (c) specificity (tau) on pseudobulk mean log CPM by cell type ------
cnt <- LayerData(so, assay = "RNA", layer = "counts")
ct_f <- factor(md$cell_type)
mm <- sparse.model.matrix(~0 + ct_f); colnames(mm) <- levels(ct_f)
agg <- as.matrix(cnt %*% mm)
lcpm <- log2(t(t(agg) / colSums(agg)) * 1e6 + 1)
tau <- apply(lcpm, 1, function(x) if (max(x) == 0) NA else sum(1 - x / max(x)) / (length(x) - 1))
top_ct <- colnames(lcpm)[apply(lcpm, 1, which.max)]

# ---------------- (d) Candidate concordance & tiers --------------------------------
cand <- merge(ints, orth_all[, .(gene = human, mouse)], by = "gene", all.x = TRUE)
pbc <- pb[, .(mouse = gene, cell_type, logFC_m = logFC, FDR_m = FDR, P_m = PValue, sig_m = sig,
              deseq2_concordant, pct_cells_detected)]
long <- merge(cand[!is.na(mouse), .(gene, mouse, logFC_h = logFC)], pbc, by = "mouse", allow.cartesian = TRUE)
long[, concordant := sign(logFC_h) == sign(logFC_m)]
best <- long[order(!(sig_m != "NS" & concordant), FDR_m), .SD[1], by = gene]
best <- best[, .(gene, best_cell_type = cell_type, best_logFC_m = logFC_m, best_FDR_m = FDR_m,
                 best_sig_concordant = sig_m != "NS" & concordant, best_nominal_concordant = P_m < 0.05 & concordant,
                 best_deseq2_concordant = deseq2_concordant)]
n_ct_sig <- long[, .(n_celltypes_sig_concordant = sum(sig_m != "NS" & concordant),
                     celltypes_sig_concordant = paste(cell_type[sig_m != "NS" & concordant], collapse = ";")), by = gene]
cand <- merge(cand, best, by = "gene", all.x = TRUE)
cand <- merge(cand, n_ct_sig, by = "gene", all.x = TRUE)
cand[, tau := tau[mouse]]
cand[, top_expressing_cell_type := top_ct[match(mouse, rownames(lcpm))]]
cand[, tier := fifelse(best_sig_concordant %in% TRUE, "Tier1",
               fifelse(best_nominal_concordant %in% TRUE, "Tier1b_nominal",
               fifelse((hub %in% TRUE | in_geneset %in% TRUE) & tau > 0.5, "Tier2", "Other")))]
cand[is.na(mouse), tier := "No_1to1_ortholog"]
if (!is.null(val)) {
  cand <- merge(cand, val, by = "gene", all.x = TRUE)
  cand[, GSE23343_concordant := !is.na(logFC_GSE23343) & sign(logFC_GSE23343) == sign(logFC)]
}
setorder(cand, tier, best_FDR_m)
fwrite(cand, file.path(out_dir, paste0(tissue, "_candidate_concordance.csv")))
print(cand[, .N, by = tier])

# ---------------- Figures ---------------------------------------------------------
show <- cand[tier %in% c("Tier1", "Tier1b_nominal", "Tier2")][1:min(.N, 60)]
if (nrow(show) >= 2) {
  mat <- dcast(long[gene %in% show$gene], gene ~ cell_type, value.var = "logFC_m")
  m <- as.matrix(mat[, -1]); rownames(m) <- mat$gene
  star <- dcast(long[gene %in% show$gene], gene ~ cell_type, value.var = "FDR_m")
  s <- as.matrix(star[, -1]); rownames(s) <- star$gene; s <- s[rownames(m), colnames(m), drop = FALSE]
  ha <- rowAnnotation(human_logFC = anno_barplot(show[match(rownames(m), gene), logFC],
                                                 gp = gpar(fill = "grey40")),
                      Tier = show[match(rownames(m), gene), tier])
  pdf(file.path(fig_dir, paste0("10_", tissue, "_concordance_heatmap.pdf")), width = 10, height = max(5, 0.22 * nrow(m) + 3))
  draw(Heatmap(m, name = "mouse\nlog2FC", col = colorRamp2(c(-2, 0, 2), c("#3B7DD8", "white", "#C8322F")), cluster_rows = FALSE, cluster_columns = FALSE,
               right_annotation = ha, na_col = "grey90", row_names_gp = gpar(fontsize = 8),
               cell_fun = function(j, i, x, y, w, h, fill) if (!is.na(s[i, j]) && s[i, j] < 0.05) grid.text("*", x, y),
               column_title = paste(tissue, "- STZ vs Control pseudobulk logFC for human bulk candidates (* FDR<0.05)")))
  dev.off()

  so$ct_cond <- paste(so$cell_type, so$condition, sep = " | ")
  pdp <- DotPlot(so, features = unique(show[tier != "Other"][1:min(.N, 25), mouse]), group.by = "ct_cond") +
    RotatedAxis() + ggtitle(paste(tissue, "- top candidates by cell type and condition"))
  ggsave(file.path(fig_main, paste0("Fig3D_candidate_expression_", tissue, ".pdf")), pdp, width = 14, height = 9)
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, paste0("sessionInfo_10_", tissue, ".txt")))
