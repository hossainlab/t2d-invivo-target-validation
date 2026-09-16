# 06a_sc_demux_embed.R — Phase 6 (steps 1-2), stage 1 of tissue assignment
# GSE244475: 7 10x lanes (3 Control, 4 STZ), each hashtag-multiplexed (Heart/Kidney/Liver/Spleen).
#
# Stage 1 (this script, slow, run once):
#   - light QC as in the source paper (Braithwaite 2024): 300 <= nFeature <= 5000, nCount <= 25000, mt <= 25%
#   - HTODemux + MULTIseqDemux per lane (HTOs matched by NAME; row order differs between lanes)
#   - joint RNA embedding (2000 HVG, PCA 30, Harmony by lane), Louvain clusters (res 1), UMAP
#   - tissue-restricted parenchymal program scores (kidney tubular segments, hepatocyte)
#   - cluster markers and cluster-level HTO composition
# Stage 2 (06b_sc_tissue_assign.R): assignment rules + per-tissue objects.
#
# Outputs: results/sc/demux_stage1_meta.csv.gz, demux_stage1_harmony.rds, demux_cluster_markers.csv,
#          demux_cluster_summary.csv
# Figures: figures/sc/06a_HTO_ridge_<GSM>.pdf, 06a_umap_overview.png

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(harmony)
  library(presto)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
set.seed(20260914)
options(future.globals.maxSize = 16 * 1024^3)
raw_dir <- "data/single/GSE244475_RAW"
out_dir <- "results/sc"; fig_dir <- "results/supplementary_figures/sc_qc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

samples <- data.table(
  GSM       = c("GSM7817709", "GSM7817710", "GSM7817711", "GSM7817712", "GSM7817713", "GSM7817714", "GSM7817715"),
  lane      = c("T0076_6", "T0076_9", "T0076_11", "T0072", "T0076_8", "T0076_12", "T0076_15"),
  condition = c("Control", "Control", "Control", "STZ", "STZ", "STZ", "STZ"),
  replicate = c(1, 2, 3, 1, 2, 3, 4)
)
samples[, mouse_id := paste0(condition, "_", replicate)]
fwrite(samples, file.path(out_dir, "sc_samples.csv"))
hto_names <- c("HTO-Heart", "HTO-Kidney", "HTO-Liver", "HTO-Spleen")
tissues   <- sub("^HTO-", "", hto_names)

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
  ggsave(file.path(fig_dir, paste0("06a_HTO_ridge_", s$GSM, ".pdf")), pr, width = 10, height = 8)

  obj$GSM <- s$GSM; obj$lane <- s$lane; obj$condition <- s$condition; obj$mouse_id <- s$mouse_id
  obj$hto_global <- as.character(obj$HTO_classification.global)
  obj$hto_tissue <- sub("^HTO-", "", as.character(obj$HTO_maxID))
  obj$multi_id   <- sub("^HTO-", "", as.character(obj$MULTI_ID))
  obj$n_barcodes_lane <- n_barcodes
  obj[["HTO"]] <- NULL
  obj@meta.data <- obj@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt", "GSM", "lane", "condition",
                                     "mouse_id", "hto_global", "hto_tissue", "multi_id", "HTO_margin",
                                     "n_barcodes_lane")]
  obj
}

objs <- lapply(seq_len(nrow(samples)), function(i) { message("Lane ", samples$GSM[i]); o <- load_lane(i); gc(); o })
so <- merge(objs[[1]], objs[-1]); rm(objs); gc()
so <- JoinLayers(so)
message("Cells after light QC: ", ncol(so))

so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = 30, verbose = FALSE)
so <- RunHarmony(so, group.by.vars = "lane", reduction.use = "pca", dims.use = 1:30,
                 reduction.save = "harmony", verbose = FALSE)
so <- FindNeighbors(so, reduction = "harmony", dims = 1:30, verbose = FALSE)
so <- FindClusters(so, resolution = 1, verbose = FALSE)
so <- RunUMAP(so, reduction = "harmony", dims = 1:30, verbose = FALSE)

# ---------------- Parenchymal program scores (per kidney segment; max = kidney tubular score) ----
programs <- list(
  Kidney_PT      = c("Lrp2", "Slc34a1", "Kap", "Miox"),
  Kidney_LoH_DCT = c("Umod", "Slc12a1", "Slc12a3"),
  Kidney_CD      = c("Aqp2", "Atp6v1g3", "Atp6v0d2"),
  Kidney_pan     = c("Fxyd2", "Cdh16", "Pax8"),
  Hepatocyte     = c("Alb", "Apoa1", "Apoc3", "Ttr", "Ahsg", "Serpina1c", "Apoa2")
)
lognorm <- LayerData(so, assay = "RNA", layer = "data")
ps <- sapply(programs, function(g) Matrix::colMeans(lognorm[intersect(g, rownames(lognorm)), , drop = FALSE]))

md <- as.data.table(so@meta.data, keep.rownames = "cell")
md[, cluster := as.character(so$seurat_clusters)]
for (nm in colnames(ps)) md[, paste0("prog_", nm) := ps[, nm]]
md[, score_Kidney := apply(ps[, grep("^Kidney", colnames(ps)), drop = FALSE], 1, max)]
md[, score_Liver := ps[, "Hepatocyte"]]
um <- Embeddings(so, "umap"); md[, `:=`(UMAP_1 = um[, 1], UMAP_2 = um[, 2])]

# ---------------- Cluster markers and HTO composition ------------------------------------
mk <- as.data.table(wilcoxauc(so, group_by = "seurat_clusters", seurat_assay = "RNA"))
mk <- mk[padj < 0.05 & logFC > 0.5][order(group, -auc)][, head(.SD, 20), by = group]
fwrite(mk, file.path(out_dir, "demux_cluster_markers.csv"))

cl <- md[, c(list(n = .N,
                  pct_singlet = round(100 * mean(hto_global == "Singlet"), 1),
                  pct_negative = round(100 * mean(hto_global == "Negative"), 1),
                  pct_doublet = round(100 * mean(hto_global == "Doublet"), 1),
                  score_Kidney = mean(score_Kidney), score_Liver = mean(score_Liver)),
             setNames(lapply(tissues, function(t) sum(hto_global == "Singlet" & hto_tissue == t)), paste0("singlet_", tissues))),
         by = cluster]
top <- mk[, .(top_markers = paste(head(feature, 12), collapse = ",")), by = .(cluster = as.character(group))]
cl <- merge(cl, top, by = "cluster", all.x = TRUE)
setorder(cl, -score_Kidney)
fwrite(cl, file.path(out_dir, "demux_cluster_summary.csv"))
print(cl[, .(cluster, n, pct_negative, score_Kidney = round(score_Kidney, 2), score_Liver = round(score_Liver, 2),
             singlet_Heart, singlet_Kidney, singlet_Liver, singlet_Spleen, top = substr(top_markers, 1, 70))])

# ---------------- Save stage-1 outputs --------------------------------------------------
fwrite(md, file.path(out_dir, "demux_stage1_meta.csv.gz"))
saveRDS(Embeddings(so, "harmony")[md$cell, 1:30], file.path(out_dir, "demux_stage1_harmony.rds"))

p1 <- ggplot(md, aes(UMAP_1, UMAP_2, colour = hto_global)) + geom_point(size = 0.05, alpha = 0.4) +
  guides(colour = guide_legend(override.aes = list(size = 3))) + theme_void() + ggtitle("HTODemux class")
p2 <- ggplot(md, aes(UMAP_1, UMAP_2, colour = score_Kidney)) + geom_point(size = 0.05, alpha = 0.4) +
  scale_colour_viridis_c() + theme_void() + ggtitle("Kidney tubular program (max segment)")
cent <- md[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = cluster]
p3 <- ggplot(md, aes(UMAP_1, UMAP_2, colour = cluster)) + geom_point(size = 0.05, alpha = 0.4, show.legend = FALSE) +
  geom_text(data = cent, aes(label = cluster), colour = "black", size = 3) + theme_void() + ggtitle("Clusters")
p4 <- ggplot(md[hto_global == "Singlet"], aes(UMAP_1, UMAP_2, colour = hto_tissue)) + geom_point(size = 0.05, alpha = 0.4) +
  guides(colour = guide_legend(override.aes = list(size = 3))) + theme_void() + ggtitle("HTO singlet tissue")
ggsave(file.path(fig_dir, "06a_umap_overview.png"), (p1 | p2) / (p3 | p4), width = 14, height = 12, dpi = 120)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_06a.txt"))
