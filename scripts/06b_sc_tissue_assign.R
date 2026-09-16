# 06b_sc_tissue_assign.R — Phase 6 (step 3), stage 2 of tissue assignment (see docs/decisions_log.md S2)
# Input: stage-1 outputs of 06a_sc_demux_embed.R
#
# Rules (fixed before looking at DE results):
#   A. OVERRIDE clusters (tissue-restricted parenchyma; hashtag unreliable there):
#        Kidney if cluster mean kidney-tubular score >= 0.3, >= 10x median cluster score, and > 3x hepatocyte score
#        Liver  if cluster mean hepatocyte score     >= 0.5, >= 10x median cluster score, and > 3x kidney score
#   B. Training set = confident HTO singlets (HTODemux singlet, MULTIseq same tissue,
#      HTO margin >= lane median of singlets), outside override clusters.
#   C. kNN (k = 30, Harmony 30 dims). 5-fold CV on the training set -> per-cluster accuracy among
#      high-confidence predictions (pmax >= 0.8). A cluster is "rescue-eligible" if that accuracy
#      >= 0.90 with >= 30 CV cells.
#   D. Final tissue:
#        override cluster                                         -> RNA tissue   ("RNA_override")
#        HTO singlet and kNN prob(HTO tissue) >= 0.2              -> HTO tissue   ("HTO")
#        HTO singlet and kNN prob(HTO tissue) <  0.2              -> dropped      ("dropped_HTO_RNA_conflict")
#        HTO negative, kNN pmax >= 0.8, rescue-eligible cluster   -> kNN tissue   ("kNN_rescue")
#        HTO negative, kNN pmax >= 0.8, non-eligible cluster      -> dropped      ("dropped_rescue_unreliable")
#        HTO doublet                                              -> dropped      ("dropped_HTO_doublet")
#        otherwise                                                -> dropped      ("dropped_ambiguous")
#
# Outputs: results/sc/demux_cell_assignment.csv.gz, demux_override_clusters.csv, demux_knn_cv_cluster.csv,
#          demux_summary.csv, demux_source_by_tissue.csv, raw_<tissue>.rds, qc_metrics_<tissue>.csv
# Figures: figures/sc/06b_assignment_umap.png, 06b_composition.pdf

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
set.seed(20260914)
raw_dir <- "data/single/GSE244475_RAW"
out_dir <- "results/sc"; fig_dir <- "results/supplementary_figures/sc_qc"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
tissues <- c("Heart", "Kidney", "Liver", "Spleen")

md  <- fread(file.path(out_dir, "demux_stage1_meta.csv.gz"), colClasses = c(cluster = "character"))
emb <- readRDS(file.path(out_dir, "demux_stage1_harmony.rds"))
stopifnot(identical(rownames(emb), md$cell))
cls <- fread(file.path(out_dir, "demux_cluster_summary.csv"), colClasses = c(cluster = "character"))

# ---------------- A. Override clusters ------------------------------------------------
med_k <- median(cls$score_Kidney); med_l <- median(cls$score_Liver)
cls[, override_tissue := fifelse(score_Kidney >= 0.3 & score_Kidney >= 10 * med_k & score_Kidney > 3 * score_Liver, "Kidney",
                         fifelse(score_Liver >= 0.5 & score_Liver >= 10 * med_l & score_Liver > 3 * score_Kidney, "Liver",
                                 NA_character_))]
override <- cls[!is.na(override_tissue), .(cluster, n, override_tissue, score_Kidney, score_Liver, pct_negative,
                                           singlet_Heart, singlet_Kidney, singlet_Liver, singlet_Spleen, top_markers)]
message("Median cluster kidney score ", signif(med_k, 3), "; hepatocyte score ", signif(med_l, 3))
print(override)
fwrite(override, file.path(out_dir, "demux_override_clusters.csv"))
md[, override_tissue := cls$override_tissue[match(cluster, cls$cluster)]]

# ---------------- B. Training set -----------------------------------------------------
md[, margin_med := median(HTO_margin[hto_global == "Singlet"]), by = lane]
md[, confident := hto_global == "Singlet" & multi_id == hto_tissue & HTO_margin >= margin_med & is.na(override_tissue)]
train_idx <- which(md$confident)
message("Training cells: ", length(train_idx)); print(md[train_idx, .N, by = hto_tissue])

knn_prob <- function(query_idx, ref_idx, k = 30) {
  nn <- RANN::nn2(emb[ref_idx, , drop = FALSE], emb[query_idx, , drop = FALSE], k = k)$nn.idx
  lab <- match(md$hto_tissue[ref_idx], tissues)
  lab_mat <- matrix(lab[nn], nrow = nrow(nn))
  p <- sapply(seq_along(tissues), function(j) rowMeans(lab_mat == j)); colnames(p) <- tissues; p
}

# ---------------- C. CV and rescue-eligible clusters -----------------------------------
folds <- sample(rep(1:5, length.out = length(train_idx)))
cv <- rbindlist(lapply(1:5, function(f) {
  q <- train_idx[folds == f]; r <- train_idx[folds != f]
  p <- knn_prob(q, r)
  data.table(idx = q, true = md$hto_tissue[q], pred = tissues[max.col(p, ties.method = "first")], pmax = apply(p, 1, max))
}))
cv[, cluster := md$cluster[idx]]
cv_overall <- rbind(
  cv[, .(level = "overall", n = .N, acc = mean(true == pred), frac_hi = mean(pmax >= 0.8), acc_hi = mean((true == pred)[pmax >= 0.8]))],
  cv[, .(n = .N, acc = mean(true == pred), frac_hi = mean(pmax >= 0.8), acc_hi = mean((true == pred)[pmax >= 0.8])), by = .(level = paste0("true_", true))]
)
print(cv_overall)
cv_cl <- cv[pmax >= 0.8, .(n_cv_hi = .N, acc_hi = mean(true == pred)), by = cluster]
cv_cl <- merge(cls[, .(cluster, n, top_markers)], cv_cl, by = "cluster", all.x = TRUE)
cv_cl[, rescue_eligible := !is.na(acc_hi) & n_cv_hi >= 30 & acc_hi >= 0.90]
setorder(cv_cl, -rescue_eligible, acc_hi)
fwrite(rbind(cv_overall[, .(level, n, acc, frac_hi, acc_hi)], fill = TRUE), file.path(out_dir, "demux_knn_cv.csv"))
fwrite(cv_cl, file.path(out_dir, "demux_knn_cv_cluster.csv"))
print(cv_cl[, .(cluster, n, n_cv_hi, acc_hi = round(acc_hi, 3), rescue_eligible, top = substr(top_markers, 1, 60))])
md[, rescue_eligible := cv_cl$rescue_eligible[match(cluster, cv_cl$cluster)] %in% TRUE]

# ---------------- D. Final tissue -----------------------------------------------------
P <- knn_prob(seq_len(nrow(md)), train_idx)
md[, knn_tissue := tissues[max.col(P, ties.method = "first")]]
md[, knn_pmax := apply(P, 1, max)]
md[, knn_p_hto := P[cbind(seq_len(.N), match(hto_tissue, tissues))]]

md[, `:=`(tissue = NA_character_, tissue_source = NA_character_)]
md[!is.na(override_tissue), `:=`(tissue = override_tissue, tissue_source = "RNA_override")]
md[is.na(tissue_source) & hto_global == "Singlet" & knn_p_hto >= 0.2, `:=`(tissue = hto_tissue, tissue_source = "HTO")]
md[is.na(tissue_source) & hto_global == "Singlet", tissue_source := "dropped_HTO_RNA_conflict"]
md[is.na(tissue_source) & hto_global == "Negative" & knn_pmax >= 0.8 & rescue_eligible,
   `:=`(tissue = knn_tissue, tissue_source = "kNN_rescue")]
md[is.na(tissue_source) & hto_global == "Negative" & knn_pmax >= 0.8, tissue_source := "dropped_rescue_unreliable"]
md[is.na(tissue_source) & hto_global == "Doublet", tissue_source := "dropped_HTO_doublet"]
md[is.na(tissue_source), tissue_source := "dropped_ambiguous"]

fwrite(md[, .(cell, mouse_id, condition, cluster, hto_global, hto_tissue, multi_id, HTO_margin, confident,
              override_tissue, rescue_eligible, knn_tissue, knn_pmax, knn_p_hto, tissue, tissue_source,
              score_Kidney, score_Liver)],
       file.path(out_dir, "demux_cell_assignment.csv.gz"))

summ <- dcast(md[!is.na(tissue), .N, by = .(mouse_id, condition, tissue)], mouse_id + condition ~ tissue, value.var = "N", fill = 0)
src  <- dcast(md[, .N, by = .(mouse_id, tissue_source)], mouse_id ~ tissue_source, value.var = "N", fill = 0)
summ <- merge(summ, src, by = "mouse_id")
summ[, cells_light_qc := md[, .N, by = mouse_id][match(summ$mouse_id, mouse_id), N]]
setorder(summ, condition, mouse_id)
print(summ)
fwrite(summ, file.path(out_dir, "demux_summary.csv"))
src_t <- dcast(md[!is.na(tissue), .N, by = .(tissue, condition, tissue_source)], tissue + condition ~ tissue_source,
               value.var = "N", fill = 0)
print(src_t)
fwrite(src_t, file.path(out_dir, "demux_source_by_tissue.csv"))

# ---------------- Figures ------------------------------------------------------------
g1 <- ggplot(md, aes(UMAP_1, UMAP_2, colour = tissue)) + geom_point(size = 0.05, alpha = 0.4) +
  guides(colour = guide_legend(override.aes = list(size = 3))) + theme_void() + ggtitle("Final tissue (NA = dropped)")
g2 <- ggplot(md, aes(UMAP_1, UMAP_2, colour = tissue_source)) + geom_point(size = 0.05, alpha = 0.4) +
  guides(colour = guide_legend(override.aes = list(size = 3))) + theme_void() + ggtitle("Assignment source")
ggsave(file.path(fig_dir, "06b_assignment_umap.png"), g1 | g2, width = 16, height = 7, dpi = 120)
comp <- melt(summ[, c("mouse_id", intersect(tissues, names(summ))), with = FALSE], id.vars = "mouse_id",
             variable.name = "tissue", value.name = "cells")
pc <- ggplot(comp, aes(mouse_id, cells, fill = tissue)) + geom_col() +
  labs(x = NULL, y = "Assigned cells", title = "GSE244475 hybrid tissue assignment") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "06b_composition.pdf"), pc, width = 7, height = 4.5)

# ---------------- Per-tissue raw objects (re-read counts lane by lane) --------------------
keep <- md[!is.na(tissue)]
samples <- fread(file.path(out_dir, "sc_samples.csv"))
parts <- setNames(vector("list", length(tissues)), tissues)
for (i in seq_len(nrow(samples))) {
  s <- samples[i]
  pre <- file.path(raw_dir, paste0(s$GSM, "_", s$lane, "_"))
  mat <- ReadMtx(mtx = paste0(pre, "matrix.mtx.gz"), cells = paste0(pre, "barcodes.tsv.gz"),
                 features = paste0(pre, "features.tsv.gz"), feature.column = 2, unique.features = TRUE)
  colnames(mat) <- paste0(s$mouse_id, "_", colnames(mat))
  mat <- mat[!grepl("^HTO-", rownames(mat)), , drop = FALSE]
  for (tis in tissues) {
    cells <- keep[mouse_id == s$mouse_id & tissue == tis, cell]
    parts[[tis]][[s$mouse_id]] <- mat[, cells, drop = FALSE]
  }
  rm(mat); gc()
}
meta_cols <- c("cell", "GSM", "lane", "condition", "mouse_id", "hto_global", "hto_tissue", "HTO_margin",
               "cluster", "tissue", "tissue_source", "knn_pmax", "score_Kidney", "score_Liver")
for (tis in tissues) {
  cnt <- do.call(cbind, parts[[tis]])
  mdt <- as.data.frame(keep[tissue == tis, ..meta_cols]); rownames(mdt) <- mdt$cell; mdt$cell <- NULL
  setnames(mdt, "cluster", "demux_cluster")
  st <- CreateSeuratObject(counts = cnt, meta.data = mdt[colnames(cnt), ])
  st[["percent.mt"]]   <- PercentageFeatureSet(st, pattern = "^mt-")
  st[["percent.hb"]]   <- PercentageFeatureSet(st, features = intersect(c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt"), rownames(st)))
  st[["percent.ribo"]] <- PercentageFeatureSet(st, pattern = "^Rp[sl]")
  saveRDS(st, file.path(out_dir, paste0("raw_", tolower(tis), ".rds")))
  qc_tab <- as.data.table(st@meta.data)[, .(cells = .N,
      from_HTO = sum(tissue_source == "HTO"), from_kNN = sum(tissue_source == "kNN_rescue"),
      from_override = sum(tissue_source == "RNA_override"),
      median_nFeature = as.numeric(median(nFeature_RNA)), median_nCount = as.numeric(median(nCount_RNA)),
      median_mt = median(percent.mt), pct_mt_lt10 = round(100 * mean(percent.mt < 10), 1),
      median_hb = median(percent.hb)), by = .(mouse_id, condition)][order(condition, mouse_id)]
  fwrite(qc_tab, file.path(out_dir, paste0("qc_metrics_", tolower(tis), ".csv")))
  message(tis, ": ", ncol(st), " cells"); print(qc_tab)
  rm(st, cnt); gc()
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_06b.txt"))
