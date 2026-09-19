# 07_sc_qc_cluster.R — Phase 6 (steps 4-8) of docs/analysis_plan.md
# Usage: Rscript scripts/07_sc_qc_cluster.R liver|kidney
#
# QC (reference-paper rules kept where valid): genes in >= 3 cells, nFeature >= 200, percent.hb < 3;
#   percent.mt: tissue MAD rule (median + 3 MAD), never stricter than the paper's 10%, capped at liver 20 /
#   kidney 30 (mito-rich hepatocytes and tubular cells); cells passing the paper's <10% rule reported as sensitivity.
#   Upper nCount/nFeature: median + 5 MAD (log10 scale).
# Ambient RNA: DecontX per mouse (filtered matrices only). Doublets: scDblFinder per mouse.
# LogNormalize, 2000 HVG, PCA; Harmony by mouse if it improves mixing; UMAP on 20 PCs (as paper);
# Louvain over resolutions 0.2-1.2; presto markers; SingleR (ImmGen + MouseRNAseq).
#
# Outputs: results/sc/<tissue>_processed.rds, <tissue>_qc_summary.csv, <tissue>_markers_res*.csv,
#          <tissue>_singler_cluster.csv, <tissue>_integration_check.csv
# Figures: figures/sc/07_<tissue>_*.pdf

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(celda)
  library(scDblFinder)
  library(harmony)
  library(presto)
  library(SingleR)
  library(celldex)
  library(clustree)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})
set.seed(20260914)
options(future.globals.maxSize = 8 * 1024^3)

TISSUES <- c("liver", "kidney", "heart", "spleen")
tissue  <- tolower(commandArgs(trailingOnly = TRUE)[1])
# a bare stopifnot() here reports "tissue %in% c(...) is not TRUE", which reads as though the tissue
# names were wrong when the real cause is almost always a missing argument
if (is.na(tissue) || !tissue %in% TISSUES)
  stop("usage: Rscript scripts/07_sc_qc_cluster.R <", paste(TISSUES, collapse = "|"), ">\n",
       if (is.na(tissue)) "  no tissue argument was given."
       else paste0("  got '", tissue, "', which is not one of the four tissues."),
       call. = FALSE)
out_dir <- "results/sc"; fig_dir <- "results/supplementary_figures/sc_qc"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
mt_cap <- c(liver = 20, kidney = 30, heart = 30, spleen = 10)[[tissue]]

so <- readRDS(file.path(out_dir, paste0("raw_", tissue, ".rds")))
n0 <- ncol(so)
md <- as.data.table(so@meta.data, keep.rownames = "cell")

# ---------------- QC thresholds ------------------------------------------------------
mad_hi <- function(x, k) median(x) + k * mad(x)
thr <- list(
  nFeature_min = 200,
  hb_max       = 3,
  mt_max       = min(max(mad_hi(md$percent.mt, 3), 10), mt_cap),  # floor = reference paper's 10%
  nCount_max   = 10^mad_hi(log10(md$nCount_RNA), 5),
  nFeature_max = 10^mad_hi(log10(md$nFeature_RNA), 5)
)
md[, pass_qc := nFeature_RNA >= thr$nFeature_min & percent.hb < thr$hb_max & percent.mt <= thr$mt_max &
                nCount_RNA <= thr$nCount_max & nFeature_RNA <= thr$nFeature_max]
md[, pass_paper_mt10 := nFeature_RNA >= 200 & percent.hb < 3 & percent.mt < 10]
print(unlist(thr))

p_qc <- ggplot(md, aes(log10(nCount_RNA), percent.mt, colour = pass_qc)) +
  geom_point(size = 0.2, alpha = 0.3) +
  geom_hline(yintercept = 10, linetype = 3) + geom_hline(yintercept = thr$mt_max, linetype = 2) +
  facet_wrap(~mouse_id, nrow = 2) + scale_colour_manual(values = c(`TRUE` = "grey30", `FALSE` = "red")) +
  labs(title = paste(tissue, "QC: dashed = MAD mito cutoff, dotted = paper 10%")) + theme_bw()
ggsave(file.path(fig_dir, paste0("07_", tissue, "_qc_scatter.pdf")), p_qc, width = 12, height = 6)

so <- subset(so, cells = md[pass_qc == TRUE, cell])
counts <- LayerData(so, assay = "RNA", layer = "counts")
counts <- counts[rowSums(counts > 0) >= 3, ]
n_qc <- ncol(counts)

# ---------------- DecontX (per mouse) ------------------------------------------------
sce <- SingleCellExperiment(assays = list(counts = counts), colData = so@meta.data[colnames(counts), ])
sce <- decontX(sce, batch = sce$mouse_id, seed = 20260914, verbose = FALSE)
contam <- sce$decontX_contamination
dec <- round(decontXcounts(sce))
dec <- as(dec, "dgCMatrix")

# ---------------- scDblFinder (per mouse, raw counts) --------------------------------
sce_d <- SingleCellExperiment(assays = list(counts = counts), colData = so@meta.data[colnames(counts), ])
sce_d <- scDblFinder(sce_d, samples = "mouse_id", BPPARAM = BiocParallel::SerialParam())
dbl <- sce_d$scDblFinder.class
rm(sce, sce_d); gc()

keep <- dbl == "singlet" & contam < 0.5
so <- CreateSeuratObject(counts = dec[, keep], meta.data = so@meta.data[colnames(dec)[keep], ])
so[["RAW"]] <- CreateAssay5Object(counts = counts[, keep])
so$decontX_contamination <- contam[keep]
DefaultAssay(so) <- "RNA"

qc_sum <- md[, .(barcodes = .N, pass_qc = sum(pass_qc), pass_paper_mt10 = sum(pass_paper_mt10)), by = .(mouse_id, condition)]
fin <- as.data.table(so@meta.data)[, .(final = .N, median_contam = median(decontX_contamination)), by = mouse_id]
dbl_tab <- data.table(mouse_id = sub("_[^_]+$", "", colnames(counts)),  # cell names are <mouse_id>_<barcode>
                      doublet = dbl == "doublet", high_contam = contam >= 0.5)[
  , .(doublets = sum(doublet), high_contam = sum(high_contam & !doublet)), by = mouse_id]
qc_sum <- merge(merge(qc_sum, dbl_tab, by = "mouse_id"), fin, by = "mouse_id")
qc_sum[, mt_threshold := round(thr$mt_max, 2)]
print(qc_sum)
fwrite(qc_sum, file.path(out_dir, paste0(tissue, "_qc_summary.csv")))
message(tissue, ": ", n0, " demuxed singlets -> ", n_qc, " pass QC -> ", ncol(so), " after doublet/ambient filtering")

# ---------------- Normalize, PCA, integration check ----------------------------------
so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = 50, verbose = FALSE)
so <- RunHarmony(so, group.by.vars = "mouse_id", reduction.use = "pca", dims.use = 1:30,
                 reduction.save = "harmony", verbose = FALSE)

mixing_entropy <- function(emb, labels, k = 30) {
  nn <- RANN::nn2(emb, k = k + 1)$nn.idx[, -1]
  lab <- as.integer(factor(labels)); L <- length(unique(lab))
  ent <- apply(nn, 1, function(i) { p <- tabulate(lab[i], L) / k; p <- p[p > 0]; -sum(p * log(p)) / log(L) })
  mean(ent)
}
ent_pca <- mixing_entropy(Embeddings(so, "pca")[, 1:20], so$mouse_id)
ent_har <- mixing_entropy(Embeddings(so, "harmony")[, 1:20], so$mouse_id)
cond_pca <- mixing_entropy(Embeddings(so, "pca")[, 1:20], so$condition)
use_red <- if (ent_har - ent_pca > 0.05) "harmony" else "pca"
fwrite(data.table(tissue, mouse_entropy_pca = ent_pca, mouse_entropy_harmony = ent_har,
                  condition_entropy_pca = cond_pca, reduction_used = use_red),
       file.path(out_dir, paste0(tissue, "_integration_check.csv")))
message("Mouse mixing entropy PCA=", round(ent_pca, 3), " Harmony=", round(ent_har, 3), " -> using ", use_red)

so <- RunUMAP(so, reduction = use_red, dims = 1:20, verbose = FALSE)
so <- FindNeighbors(so, reduction = use_red, dims = 1:20, verbose = FALSE)
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
so <- FindClusters(so, resolution = resolutions, verbose = FALSE)

ggsave(file.path(fig_dir, paste0("07_", tissue, "_clustree.pdf")),
       clustree(so, prefix = "RNA_snn_res."), width = 10, height = 10)
pu <- (DimPlot(so, group.by = "RNA_snn_res.0.6", label = TRUE) + NoLegend()) |
      DimPlot(so, group.by = "mouse_id", shuffle = TRUE) | DimPlot(so, group.by = "condition", shuffle = TRUE)
ggsave(file.path(fig_dir, paste0("07_", tissue, "_umap_overview.pdf")), pu, width = 18, height = 5.5)

# ---------------- Markers -----------------------------------------------------------
for (r in c(0.4, 0.6, 1.0)) {
  col <- paste0("RNA_snn_res.", r)
  mk <- as.data.table(wilcoxauc(so, group_by = col, seurat_assay = "RNA"))
  mk <- mk[padj < 0.05 & logFC > 0.25][order(group, -auc)][, head(.SD, 30), by = group]
  fwrite(mk, file.path(out_dir, paste0(tissue, "_markers_res", r, ".csv")))
}

canon <- list(
  liver = c("Alb", "Apoa1", "Cyp2e1", "Cyp2f2", "Kdr", "Stab2", "Clec4g", "Clec4f", "Vsig4", "Adgre1",
            "Ly6c2", "Ccr2", "Dcn", "Reln", "Lrat", "Pdgfrb", "Krt19", "Epcam", "Cd79a", "Cd3e",
            "Nkg7", "S100a8", "Siglech", "Xcr1", "Mki67"),
  kidney = c("Lrp2", "Slc34a1", "Slc5a2", "Slc22a6", "Umod", "Slc12a1", "Slc12a3", "Aqp2", "Atp6v1g3",
             "Nphs1", "Nphs2", "Emcn", "Plvap", "Pecam1", "Pdgfrb", "Col1a1", "Acta2", "C1qa",
             "Adgre1", "Ly6c2", "Cd79a", "Cd3e", "Nkg7", "S100a8", "Mki67"),
  heart = c("Ttn", "Myh6", "Pecam1", "Fabp4", "Col1a1", "Pdgfra", "Rgs5", "C1qa", "Cd3e", "Cd79a"),
  spleen = c("Cd79a", "Cd3e", "Nkg7", "C1qa", "Hba-a1", "S100a8", "Siglech", "Xcr1")
)[[tissue]]
pd <- DotPlot(so, features = intersect(canon, rownames(so)), group.by = "RNA_snn_res.0.6") +
  RotatedAxis() + ggtitle(paste(tissue, "canonical markers (res 0.6)"))
ggsave(file.path(fig_dir, paste0("07_", tissue, "_canonical_dotplot.pdf")), pd, width = 13, height = 7)

# ---------------- SingleR ------------------------------------------------------------
refs <- list(ImmGen = celldex::ImmGenData(), MouseRNAseq = celldex::MouseRNAseqData())
avg <- AggregateExpression(so, group.by = "RNA_snn_res.0.6", assays = "RNA", return.seurat = FALSE)$RNA
colnames(avg) <- sub("^g", "", colnames(avg))  # Seurat prefixes numeric group names with "g"
avg <- log1p(t(t(avg) / colSums(avg)) * 1e4)
sr_tab <- data.table(cluster = colnames(avg))
for (nm in names(refs)) {
  pc <- SingleR(test = avg, ref = refs[[nm]], labels = refs[[nm]]$label.main)
  sr_tab[, paste0(nm, "_label") := pc$labels]
  sc <- SingleR(test = LayerData(so, "RNA", layer = "data"), ref = refs[[nm]],
                labels = refs[[nm]]$label.main, BPPARAM = BiocParallel::SerialParam())
  so[[paste0("SingleR_", nm)]] <- sc$pruned.labels
}
sizes <- as.data.table(table(cluster = so$RNA_snn_res.0.6))
sr_tab <- merge(sr_tab, setnames(sizes, "N", "n_cells"), by = "cluster")
comp <- dcast(as.data.table(so@meta.data)[, .N, by = .(cluster = RNA_snn_res.0.6, condition)],
              cluster ~ condition, value.var = "N", fill = 0)
sr_tab <- merge(sr_tab, comp, by = "cluster")
fwrite(sr_tab[order(as.integer(cluster))], file.path(out_dir, paste0(tissue, "_singler_cluster.csv")))
print(sr_tab[order(as.integer(cluster))])

saveRDS(so, file.path(out_dir, paste0(tissue, "_processed.rds")))
writeLines(capture.output(sessionInfo()), file.path(out_dir, paste0("sessionInfo_07_", tissue, ".txt")))
