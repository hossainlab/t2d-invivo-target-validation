# 08_sc_annotate.R — Phase 6 (step 8) of docs/analysis_plan.md
# Usage: Rscript scripts/08_sc_annotate.R liver|kidney [resolution, default 0.6]
#
# Cluster-level annotation, as in the reference paper (canonical markers; CellMarker 2.0 / PanglaoDB /
# Braithwaite 2024), with an automated first pass and a manual override file:
#   1. For each marker set: per-cluster average expression of each marker (log-normalized), z-scored
#      across clusters, averaged over the set's genes -> set score per cluster.
#   2. Label = best-scoring set. Confidence gap = best - second best.
#      Flags: "low_confidence" (gap < 0.3), "mixed_lineage" (>= 2 sets from different lineages score > 1).
#   3. SingleR (from 07) majority label per cluster reported for cross-checking.
#   4. Manual overrides: docs/annotation/<tissue>_overrides.csv with columns cluster,cell_type
#      (use cell_type = "Exclude" to drop a cluster, e.g. doublets or low quality).
#   5. Broad classes for the intake's priority populations: Epithelial, Endothelial, Stromal,
#      Myeloid, Lymphoid, Other.
#
# Inputs:  results/sc/<tissue>_processed.rds (from 07)
# Outputs: results/sc/<tissue>_cluster_annotation_auto.csv, <tissue>_annotated.rds,
#          <tissue>_celltype_counts.csv
# Figures: figures/sc/08_<tissue>_umap_celltype.pdf, 08_<tissue>_marker_dotplot.pdf (Fig 3B analogue),
#          08_<tissue>_composition.pdf (Fig 3C analogue)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
set.seed(20260914)
args <- commandArgs(trailingOnly = TRUE)
tissue <- tolower(args[1]); stopifnot(tissue %in% c("liver", "kidney"))
res <- if (length(args) >= 2) args[2] else "0.6"
out_dir <- "results/sc"; fig_dir <- "figures/Fig3_scRNA"; ann_dir <- "docs/annotation"  # Fig 3A-C analogue
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ann_dir, recursive = TRUE, showWarnings = FALSE)

shared <- list(
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
sets <- list(
  liver = c(list(
    Hepatocyte        = c("Alb", "Apoa1", "Apoc3", "Cyp2e1", "Ttr", "Serpina1c"),
    LSEC              = c("Stab2", "Clec4g", "Dnase1l3", "Fcgr2b", "Kdr"),
    Endothelial_vasc  = c("Vwf", "Efnb2", "Ly6a", "Rspo3"),
    Kupffer_cell      = c("Clec4f", "Vsig4", "Cd5l", "Marco"),
    Macrophage_other  = c("C1qa", "C1qb", "Adgre1", "Cx3cr1"),
    HSC_Fibroblast    = c("Dcn", "Reln", "Lrat", "Pdgfrb", "Colec11", "Ecm1"),
    Cholangiocyte     = c("Krt19", "Epcam", "Spp1", "Sox9")
  ), shared),
  kidney = c(list(
    PT                = c("Lrp2", "Slc34a1", "Kap", "Miox", "Slc22a6"),
    LoH               = c("Umod", "Slc12a1"),
    DCT               = c("Slc12a3", "Pvalb"),
    CD_PC             = c("Aqp2", "Hsd11b2"),
    CD_IC             = c("Atp6v1g3", "Atp6v0d2"),
    Podocyte          = c("Nphs1", "Nphs2", "Podxl"),
    Endothelial       = c("Emcn", "Plvap", "Pecam1", "Kdr", "Ehd3"),
    Fibroblast_Pericyte = c("Pdgfrb", "Col1a1", "Dcn", "Rgs5"),
    Macrophage        = c("C1qa", "C1qb", "Adgre1", "Cd68")
  ), shared)
)[[tissue]]
lineage <- c(Hepatocyte = "Epithelial", Cholangiocyte = "Epithelial", PT = "Epithelial", LoH = "Epithelial",
             DCT = "Epithelial", CD_PC = "Epithelial", CD_IC = "Epithelial", Podocyte = "Epithelial",
             Distal_tubule = "Epithelial", Urothelium = "Epithelial", Tubular_epithelium_other = "Epithelial",  # manual override labels
             LSEC = "Endothelial", Endothelial_vasc = "Endothelial", Endothelial = "Endothelial",
             HSC_Fibroblast = "Stromal", Fibroblast_Pericyte = "Stromal",
             Kupffer_cell = "Myeloid", Macrophage_other = "Myeloid", Macrophage = "Myeloid", Monocyte = "Myeloid",
             Neutrophil = "Myeloid", cDC = "Myeloid", pDC = "Myeloid",
             B_cell = "Lymphoid", Plasma_cell = "Lymphoid", T_cell = "Lymphoid", NK_cell = "Lymphoid",
             Proliferating = "Other")

so <- readRDS(file.path(out_dir, paste0(tissue, "_processed.rds")))
cl_col <- paste0("RNA_snn_res.", res); stopifnot(cl_col %in% colnames(so@meta.data))
so$cluster <- as.character(so[[cl_col, drop = TRUE]])

genes <- intersect(unique(unlist(sets)), rownames(so))
# per-cluster mean of log-normalized expression (data layer)
dat <- LayerData(so, "RNA", layer = "data")[genes, , drop = FALSE]
cl_f <- factor(so$cluster)
mm <- Matrix::sparse.model.matrix(~0 + cl_f); colnames(mm) <- levels(cl_f)
avg <- as.matrix(dat %*% mm) / rep(colSums(mm), each = nrow(dat))
z <- t(scale(t(avg))); z[is.na(z)] <- 0

set_scores <- sapply(sets, function(g) { g <- intersect(g, rownames(z)); if (!length(g)) rep(NA, ncol(z)) else colMeans(z[g, , drop = FALSE]) })
rownames(set_scores) <- colnames(z)

ann <- rbindlist(lapply(rownames(set_scores), function(k) {
  s <- sort(set_scores[k, ], decreasing = TRUE)
  hi <- names(s)[s > 1]
  data.table(cluster = k, n_cells = sum(so$cluster == k),
             auto_label = names(s)[1], score1 = round(s[1], 2),
             second = names(s)[2], score2 = round(s[2], 2), gap = round(s[1] - s[2], 2),
             mixed_lineage = length(unique(lineage[hi])) >= 2,
             high_sets = paste(hi, collapse = ";"))
}))
ann[, low_confidence := gap < 0.3]

for (ref in c("SingleR_ImmGen", "SingleR_MouseRNAseq")) {
  if (ref %in% colnames(so@meta.data)) {
    maj <- as.data.table(so@meta.data)[, .(lab = names(which.max(table(get(ref)))),
                                           frac = round(max(table(get(ref))) / .N, 2)), by = cluster]
    ann <- merge(ann, maj[, setNames(.SD, c("cluster", paste0(ref, "_major"), paste0(ref, "_frac")))], by = "cluster", all.x = TRUE)
  }
}
mk_file <- file.path(out_dir, paste0(tissue, "_markers_res", res, ".csv"))
if (file.exists(mk_file)) {
  mk <- fread(mk_file)
  ann <- merge(ann, mk[, .(top_markers = paste(head(feature, 12), collapse = ",")), by = .(cluster = as.character(group))],
               by = "cluster", all.x = TRUE)
}

# ---------------- Manual overrides -------------------------------------------------
ann[, cell_type := auto_label]
ov_file <- file.path(ann_dir, paste0(tissue, "_overrides.csv"))
if (file.exists(ov_file)) {
  ov <- fread(ov_file, colClasses = c(cluster = "character"))
  ann[ov, on = "cluster", cell_type := i.cell_type]
  message("Applied ", nrow(ov), " manual overrides from ", ov_file)
} else {
  message("No override file (", ov_file, "); automated labels used. Review ", tissue, "_cluster_annotation_auto.csv")
}
ann[, cell_class := fifelse(cell_type == "Exclude", "Exclude", unname(lineage[cell_type]))]
ann[is.na(cell_class), cell_class := "Other"]
setorder(ann, cell_type, cluster)
fwrite(ann, file.path(out_dir, paste0(tissue, "_cluster_annotation_auto.csv")))
print(ann[, .(cluster, n_cells, auto_label, score1, second, gap, low_confidence, mixed_lineage, cell_type,
              top = substr(top_markers, 1, 50))])

so$cell_type  <- ann$cell_type[match(so$cluster, ann$cluster)]
so$cell_class <- ann$cell_class[match(so$cluster, ann$cluster)]
so$annotation_flag <- fifelse(ann$mixed_lineage[match(so$cluster, ann$cluster)], "mixed_lineage",
                      fifelse(ann$low_confidence[match(so$cluster, ann$cluster)], "low_confidence", "ok"))
n_excl <- sum(so$cell_type == "Exclude")
if (n_excl) so <- subset(so, subset = cell_type != "Exclude")
message(tissue, ": ", ncol(so), " cells annotated (", n_excl, " excluded)")

counts <- dcast(as.data.table(so@meta.data)[, .N, by = .(cell_type, mouse_id)], cell_type ~ mouse_id, value.var = "N", fill = 0)
fwrite(counts, file.path(out_dir, paste0(tissue, "_celltype_counts.csv")))
print(counts)

# ---------------- Figures -----------------------------------------------------------
p_umap <- DimPlot(so, group.by = "cell_type", label = TRUE, repel = TRUE, raster = TRUE) + NoLegend() +
  ggtitle(paste(tissue, "cell types")) |
  DimPlot(so, group.by = "condition", shuffle = TRUE, raster = TRUE)
ggsave(file.path(fig_dir, paste0("Fig3A_umap_celltype_", tissue, ".pdf")), p_umap, width = 14, height = 6)

top_genes <- unique(unlist(lapply(sets[unique(so$cell_type)[unique(so$cell_type) %in% names(sets)]], head, 3)))
p_dot <- DotPlot(so, features = intersect(top_genes, rownames(so)), group.by = "cell_type") + RotatedAxis() +
  ggtitle(paste(tissue, "marker genes by cell type"))
ggsave(file.path(fig_dir, paste0("Fig3B_marker_bubble_", tissue, ".pdf")), p_dot, width = 14, height = 6)

prop <- as.data.table(so@meta.data)[, .N, by = .(mouse_id, condition, cell_type)][, frac := N / sum(N), by = mouse_id]
p_comp <- ggplot(prop, aes(mouse_id, frac, fill = cell_type)) + geom_col() +
  facet_grid(~condition, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Fraction of cells", title = paste(tissue, "composition by mouse")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, paste0("Fig3C_cell_proportion_", tissue, ".pdf")), p_comp, width = 8, height = 5)

saveRDS(so, file.path(out_dir, paste0(tissue, "_annotated.rds")))
writeLines(capture.output(sessionInfo()), file.path(out_dir, paste0("sessionInfo_08_", tissue, ".txt")))
