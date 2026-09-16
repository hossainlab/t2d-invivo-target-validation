# 06_sc_load_demux.R — Phase 6 (steps 1-3) of docs/analysis_plan.md
# GSE244475: 7 10x lanes (3 Control, 4 STZ), each hashtag-multiplexed (Heart/Kidney/Liver/Spleen).
# Comparison downstream: STZ vs Control.
#
# WHY HYBRID TISSUE ASSIGNMENT (see docs/decisions_log.md S1/S2):
#   HTO signal is weak and uneven. A cluster-level check (results/sc/diag_cluster_hto_purity.csv) showed
#   that liver-restricted clusters (LSEC, Kupffer, hepatocytes) and kidney endothelium carry the correct
#   hashtag (85-96%), but kidney tubular epithelium is mostly HTO-Negative and its few "singlet" calls
#   are spread over all 4 tissues; stromal cells are also largely HTO-Negative.
#
# Tissue assignment:
#   1. Light QC as in the source paper (Braithwaite 2024): 300 <= nFeature <= 5000, nCount <= 25000,
#      percent.mt <= 25.
#   2. HTODemux (+ MULTIseqDemux for concordance) per lane. HTOs selected by NAME (row order differs).
#   3. RNA embedding of all lanes (PCA on 2000 HVGs; Harmony by lane), Louvain clusters.
#   4. OVERRIDE clusters: clusters dominated by a tissue-restricted parenchymal program
#      (kidney tubular epithelium; hepatocytes) are assigned that tissue from RNA.
#   5. Training set = confident HTO singlets (HTODemux singlet, MULTIseq same tissue, HTO margin
#      >= lane median of singlets) outside override clusters.
#   6. kNN (k = 30) in Harmony space: tissue probability for every non-override cell.
#      Accuracy is estimated by 5-fold CV on the training set (reported by cluster).
#   7. Final tissue:
#        override cluster                           -> RNA tissue           (source = "RNA_override")
#        HTODemux singlet, kNN prob(HTO tissue)>=0.2 -> HTO tissue           (source = "HTO")
#        HTODemux singlet, kNN prob(HTO tissue)<0.2  -> dropped (conflict)
#        HTODemux negative / low-margin, kNN max>=0.8 -> kNN tissue           (source = "kNN_rescue")
#        HTODemux doublet or ambiguous              -> dropped
#
# Outputs: results/sc/demux_summary.csv, demux_cell_assignment.csv.gz, demux_knn_cv.csv,
#          demux_override_clusters.csv, raw_<tissue>.rds (counts + metadata), qc_metrics_<tissue>.csv
# Figures: figures/sc/06_*.pdf

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(harmony)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
set.seed(20260914)
options(future.globals.maxSize = 16 * 1024^3)
raw_dir <- "data/single/GSE244475_RAW"
out_dir <- "results/sc"; fig_dir <- "figures/sc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

samples <- data.table(
  GSM       = c("GSM7817709", "GSM7817710", "GSM7817711", "GSM7817712", "GSM7817713", "GSM7817714", "GSM7817715"),
  lane      = c("T0076_6", "T0076_9", "T0076_11", "T0072", "T0076_8", "T0076_12", "T0076_15"),
  condition = c("Control", "Control", "Control", "STZ", "STZ", "STZ", "STZ"),
  replicate = c(1, 2, 3, 1, 2, 3, 4)
)
samples[, mouse_id := paste0(condition, "_", replicate)]
hto_names <- c("HTO-Heart", "HTO-Kidney", "HTO-Liver", "HTO-Spleen")
tissues   <- sub("^HTO-", "", hto_names)

# ---------------- 1-2. Load, light QC, HTO demultiplex per lane ------------------------
load_lane <- function(i) {
  s   <- samples[i]
  pre <- file.path(raw_dir, paste0(s$GSM, "_", s$lane, "_"))
  mat <- ReadMtx(mtx = paste0(pre, "matrix.mtx.gz"), cells = paste0(pre, "barcodes.tsv.gz"),
                 features = paste0(pre, "features.tsv.gz"), feature.column = 2, unique.features = TRUE)
  colnames(mat) <- paste0(s$mouse_id, "_", colnames(mat))
  stopifnot(all(hto_names %in% rownames(mat)))
  hto <- mat[hto_names, , drop = FALSE]
  rna <- mat[!rownames(mat) %in% hto_names, , drop = FALSE]
  rm(mat)

  obj <- CreateSeuratObject(counts = rna, project = s$mouse_id)
  n_barcodes <- ncol(obj)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj <- subset(obj, subset = nFeature_RNA >= 300 & nFeature_RNA <= 5000 & nCount_RNA <= 25000 & percent.mt <= 25)
  obj[["HTO"]] <- CreateAssayObject(counts = hto[, colnames(obj)])
  obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR", margin = 2, verbose = FALSE)
  obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.99, verbose = FALSE)
  obj <- MULTIseqDemux(obj, assay = "HTO", autoThresh = TRUE, verbose = FALSE)

  Idents(obj) <- "HTO_maxID"
  pr <- RidgePlot(obj, assay = "HTO", features = rownames(obj[["HTO"]]), ncol = 2) +
    plot_annotation(title = paste(s$GSM, s$mouse_id))
  ggsave(file.path(fig_dir, paste0("06_HTO_ridge_", s$GSM, ".pdf")), pr, width = 10, height = 8)

  obj$GSM <- s$GSM; obj$lane <- s$lane; obj$condition <- s$condition; obj$mouse_id <- s$mouse_id
  obj$hto_global  <- as.character(obj$HTO_classification.global)
  obj$hto_tissue  <- sub("^HTO-", "", as.character(obj$HTO_maxID))
  obj$multi_id    <- sub("^HTO-", "", as.character(obj$MULTI_ID))
  obj$n_barcodes_lane <- n_barcodes
  obj[["HTO"]] <- NULL
  keep_cols <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "GSM", "lane", "condition",
                 "mouse_id", "hto_global", "hto_tissue", "multi_id", "HTO_margin", "n_barcodes_lane")
  obj@meta.data <- obj@meta.data[, keep_cols]
  obj
}

objs <- lapply(seq_len(nrow(samples)), function(i) { message("Lane ", samples$GSM[i]); o <- load_lane(i); gc(); o })
so <- merge(objs[[1]], objs[-1]); rm(objs); gc()
so <- JoinLayers(so)
md <- as.data.table(so@meta.data, keep.rownames = "cell")
message("Cells after light QC: ", nrow(md))

# ---------------- 3. RNA embedding ---------------------------------------------------
so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = 30, verbose = FALSE)
so <- RunHarmony(so, group.by.vars = "lane", reduction.use = "pca", dims.use = 1:30,
                 reduction.save = "harmony", verbose = FALSE)
so <- FindNeighbors(so, reduction = "harmony", dims = 1:30, verbose = FALSE)
so <- FindClusters(so, resolution = 1, verbose = FALSE)
so <- RunUMAP(so, reduction = "harmony", dims = 1:30, verbose = FALSE)
md[, cluster := as.character(so$seurat_clusters)]

# ---------------- 4. Override clusters (tissue-restricted parenchymal programs) -------
programs <- list(
  Kidney = c("Lrp2", "Slc34a1", "Kap", "Miox", "Umod", "Slc12a1", "Slc12a3", "Aqp2", "Atp6v1g3",
             "Fxyd2", "Cdh16", "Pax8"),
  Liver  = c("Alb", "Apoa1", "Apoc3", "Ttr", "Ahsg", "Serpina1c", "Fabp1", "Apoa2")
)
lognorm <- LayerData(so, assay = "RNA", layer = "data")
prog_score <- sapply(programs, function(g) Matrix::colMeans(lognorm[intersect(g, rownames(lognorm)), , drop = FALSE]))
md[, `:=`(score_Kidney = prog_score[cell, "Kidney"], score_Liver = prog_score[cell, "Liver"])]
cl_scores <- md[, .(n = .N, Kidney = mean(score_Kidney), Liver = mean(score_Liver),
                    pct_hto_negative = round(100 * mean(hto_global == "Negative"), 1)), by = cluster]
# Rule: cluster mean program score >= 1 (log-normalized) and >= 5x the median cluster score
override <- rbindlist(lapply(names(programs), function(tis) {
  x <- cl_scores[, .(cluster, n, score = get(tis), pct_hto_negative)]
  x[, `:=`(tissue = tis, median_all = median(score))]
  x[score >= 1 & score >= 5 * median_all]
}))
override <- override[order(-score)][!duplicated(cluster)][order(tissue, -score)]  # keep highest-scoring tissue per cluster
print(override)
fwrite(cl_scores[order(-Kidney)], file.path(out_dir, "demux_cluster_program_scores.csv"))
fwrite(override, file.path(out_dir, "demux_override_clusters.csv"))
md[, override_tissue := override$tissue[match(cluster, override$cluster)]]

# ---------------- 5. Training set -----------------------------------------------------
md[hto_global == "Singlet", margin_med := median(HTO_margin), by = lane]
md[, confident := hto_global == "Singlet" & multi_id == hto_tissue & HTO_margin >= margin_med & is.na(override_tissue)]
md[is.na(confident), confident := FALSE]
message("Confident HTO singlets for training: ", md[confident == TRUE, .N])
print(md[confident == TRUE, .N, by = hto_tissue])

# ---------------- 6. kNN tissue probabilities -----------------------------------------
emb <- Embeddings(so, "harmony")[md$cell, 1:30]
train_idx <- which(md$confident)
knn_prob <- function(query_idx, ref_idx, k = 30) {
  nn <- RANN::nn2(emb[ref_idx, , drop = FALSE], emb[query_idx, , drop = FALSE], k = k)$nn.idx
  lab <- factor(md$hto_tissue[ref_idx], levels = tissues)
  lab_mat <- matrix(as.integer(lab)[nn], nrow = nrow(nn))
  sapply(seq_along(tissues), function(j) rowMeans(lab_mat == j))
}

# 5-fold CV on training set
folds <- sample(rep(1:5, length.out = length(train_idx)))
cv <- rbindlist(lapply(1:5, function(f) {
  q <- train_idx[folds == f]; r <- train_idx[folds != f]
  p <- knn_prob(q, r); colnames(p) <- tissues
  data.table(cell_idx = q, true = md$hto_tissue[q], pred = tissues[max.col(p, ties.method = "first")],
             pmax = apply(p, 1, max))
}))
cv[, cluster := md$cluster[cell_idx]]
cv_summary <- rbind(
  cv[, .(level = "overall", n = .N, accuracy = mean(true == pred),
         frac_pmax_ge_0.8 = mean(pmax >= 0.8), accuracy_pmax_ge_0.8 = mean((true == pred)[pmax >= 0.8]))],
  cv[, .(level = paste0("true_", true), n = .N, accuracy = mean(true == pred),
         frac_pmax_ge_0.8 = mean(pmax >= 0.8), accuracy_pmax_ge_0.8 = mean((true == pred)[pmax >= 0.8])), by = true][, -"true"],
  cv[, .(level = paste0("cluster_", cluster), n = .N, accuracy = mean(true == pred),
         frac_pmax_ge_0.8 = mean(pmax >= 0.8), accuracy_pmax_ge_0.8 = mean((true == pred)[pmax >= 0.8])), by = cluster][, -"cluster"]
)
print(cv_summary[1:5])
fwrite(cv_summary, file.path(out_dir, "demux_knn_cv.csv"))

all_idx <- seq_len(nrow(md))
P <- knn_prob(all_idx, train_idx); colnames(P) <- tissues
md[, knn_tissue := tissues[max.col(P, ties.method = "first")]]
md[, knn_pmax := apply(P, 1, max)]
md[, knn_p_hto := P[cbind(all_idx, match(hto_tissue, tissues))]]

# ---------------- 7. Final tissue -----------------------------------------------------
md[, `:=`(tissue = NA_character_, tissue_source = NA_character_)]
md[!is.na(override_tissue), `:=`(tissue = override_tissue, tissue_source = "RNA_override")]
md[is.na(tissue) & hto_global == "Singlet" & knn_p_hto >= 0.2, `:=`(tissue = hto_tissue, tissue_source = "HTO")]
md[is.na(tissue) & hto_global == "Singlet" & knn_p_hto < 0.2, tissue_source := "dropped_HTO_RNA_conflict"]
md[is.na(tissue) & hto_global == "Negative" & knn_pmax >= 0.8, `:=`(tissue = knn_tissue, tissue_source = "kNN_rescue")]
md[is.na(tissue) & is.na(tissue_source) & hto_global == "Doublet", tissue_source := "dropped_HTO_doublet"]
md[is.na(tissue) & is.na(tissue_source), tissue_source := "dropped_ambiguous"]

fwrite(md[, .(cell, mouse_id, condition, cluster, hto_global, hto_tissue, multi_id, HTO_margin, confident,
              override_tissue, knn_tissue, knn_pmax, knn_p_hto, tissue, tissue_source, score_Kidney, score_Liver)],
       file.path(out_dir, "demux_cell_assignment.csv.gz"))

summ <- dcast(md[!is.na(tissue), .N, by = .(mouse_id, condition, tissue)], mouse_id + condition ~ tissue,
              value.var = "N", fill = 0)
src  <- dcast(md[, .N, by = .(mouse_id, tissue_source)], mouse_id ~ tissue_source, value.var = "N", fill = 0)
lane_info <- unique(md[, .(mouse_id, GSM, barcodes = n_barcodes_lane)])
summ <- Reduce(function(a, b) merge(a, b, by = "mouse_id"), list(lane_info, summ, src))
summ[, cells_light_qc := md[, .N, by = mouse_id][match(summ$mouse_id, mouse_id), N]]
print(summ)
fwrite(summ, file.path(out_dir, "demux_summary.csv"))

# Source by tissue (how much of each tissue came from rescue/override)
src_t <- dcast(md[!is.na(tissue), .N, by = .(tissue, condition, tissue_source)], tissue + condition ~ tissue_source,
               value.var = "N", fill = 0)
print(src_t)
fwrite(src_t, file.path(out_dir, "demux_source_by_tissue.csv"))

# ---------------- Figures ------------------------------------------------------------
so$hto_global <- md$hto_global; so$tissue <- md$tissue; so$tissue_source <- md$tissue_source
so$score_Kidney <- md$score_Kidney
p1 <- DimPlot(so, group.by = "hto_global", raster = TRUE) + ggtitle("HTODemux class")
p2 <- DimPlot(so, group.by = "tissue", raster = TRUE) + ggtitle("Final tissue (NA = dropped)")
p3 <- DimPlot(so, group.by = "tissue_source", raster = TRUE) + ggtitle("Assignment source")
p4 <- FeaturePlot(so, "score_Kidney", raster = TRUE) + ggtitle("Kidney tubular program")
ggsave(file.path(fig_dir, "06_demux_umap.pdf"), (p1 | p2) / (p3 | p4), width = 15, height = 12)

comp <- melt(summ[, c("mouse_id", intersect(tissues, names(summ))), with = FALSE], id.vars = "mouse_id",
             variable.name = "tissue", value.name = "cells")
pc <- ggplot(comp, aes(mouse_id, cells, fill = tissue)) + geom_col() +
  labs(x = NULL, y = "Assigned cells", title = "GSE244475 hybrid tissue assignment") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "06_demux_composition.pdf"), pc, width = 7, height = 4.5)

# ---------------- Split by tissue ---------------------------------------------------
so@meta.data <- cbind(so@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "GSM", "lane",
                                       "condition", "mouse_id")],
                      as.data.frame(md[, .(hto_global, hto_tissue, HTO_margin, tissue, tissue_source,
                                           knn_pmax, demux_cluster = cluster)]))
so[["percent.hb"]]   <- PercentageFeatureSet(so, features = intersect(c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt"), rownames(so)))
so[["percent.ribo"]] <- PercentageFeatureSet(so, pattern = "^Rp[sl]")
counts_all <- LayerData(so, assay = "RNA", layer = "counts")
meta_all <- so@meta.data
rm(so, lognorm); gc()

for (tis in tissues) {
  cells <- rownames(meta_all)[!is.na(meta_all$tissue) & meta_all$tissue == tis]
  st <- CreateSeuratObject(counts = counts_all[, cells], meta.data = meta_all[cells, ])
  saveRDS(st, file.path(out_dir, paste0("raw_", tolower(tis), ".rds")))
  qc_tab <- as.data.table(st@meta.data)[, .(cells = .N,
      from_HTO = sum(tissue_source == "HTO"), from_kNN = sum(tissue_source == "kNN_rescue"),
      from_override = sum(tissue_source == "RNA_override"),
      median_nFeature = as.numeric(median(nFeature_RNA)), median_nCount = as.numeric(median(nCount_RNA)),
      median_mt = median(percent.mt), p90_mt = as.numeric(quantile(percent.mt, 0.9)),
      pct_mt_lt10 = round(100 * mean(percent.mt < 10), 1), median_hb = median(percent.hb)),
    by = .(mouse_id, condition)][order(condition, mouse_id)]
  fwrite(qc_tab, file.path(out_dir, paste0("qc_metrics_", tolower(tis), ".csv")))
  message(tis, ": ", ncol(st), " cells"); print(qc_tab)
  rm(st); gc()
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_06.txt"))
