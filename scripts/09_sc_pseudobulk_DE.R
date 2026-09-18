# 09_sc_pseudobulk_DE.R — Phase 7 (steps 1-2) of docs/analysis_plan.md
# Usage: Rscript scripts/09_sc_pseudobulk_DE.R liver|kidney [exclude_mouse_ids, comma-separated]
#   e.g. Rscript scripts/09_sc_pseudobulk_DE.R kidney STZ_4  -> outputs get suffix "_excl_STZ_4"
# Input: results/sc/<tissue>_annotated.rds with meta columns cell_type, mouse_id, condition
#
# Pseudobulk DE (primary): sum DecontX counts per mouse x cell type; cell types with >= 20 cells
#   in >= 3 mice per condition; edgeR quasi-likelihood ~condition (STZ vs Control).
#   Check: DESeq2 Wald on same pseudobulk. Thresholds: FDR < 0.05, |log2FC| > 0.5.
# Differential abundance: speckle::propeller (logit transform).
#
# Outputs: results/sc/<tissue>_pseudobulk_DE_all.csv, <tissue>_pseudobulk_DE_summary.csv,
#          <tissue>_propeller_DA.csv
# Figures: figures/sc/09_<tissue>_DE_counts.pdf, 09_<tissue>_composition.pdf

suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(DESeq2)
  library(speckle)
  library(data.table)
  library(ggplot2)
  library(Matrix)
})
set.seed(20260914)
args <- commandArgs(trailingOnly = TRUE)
tissue <- tolower(args[1])
if (is.na(tissue) || !tissue %in% c("liver", "kidney"))
  stop("usage: Rscript scripts/09_sc_pseudobulk_DE.R <liver|kidney> [mice_to_exclude,comma,separated]\n",
       if (is.na(tissue)) "  no tissue argument was given."
       else paste0("  got '", tissue, "'."), call. = FALSE)
exclude <- if (length(args) >= 2 && nzchar(args[2])) strsplit(args[2], ",")[[1]] else character(0)
sfx <- if (length(exclude)) paste0("_excl_", paste(exclude, collapse = "-")) else ""
out_dir <- "results/sc"; fig_dir <- "results/supplementary_figures/sc_de"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
LFC <- 0.5; FDR <- 0.05; MIN_CELLS <- 20; MIN_MICE <- 3

so <- readRDS(file.path(out_dir, paste0(tissue, "_annotated.rds")))
if (length(exclude)) {
  so <- subset(so, cells = colnames(so)[!so$mouse_id %in% exclude])
  message("Sensitivity analysis: excluded ", paste(exclude, collapse = ", "), "; ", ncol(so), " cells remain")
}
stopifnot(all(c("cell_type", "mouse_id", "condition") %in% colnames(so@meta.data)))
md <- as.data.table(so@meta.data, keep.rownames = "cell")
md[, condition := factor(condition, levels = c("Control", "STZ"))]

# ---------------- Eligible cell types ---------------------------------------------
ncell <- md[, .N, by = .(cell_type, mouse_id, condition)]
elig <- ncell[N >= MIN_CELLS, .(mice = uniqueN(mouse_id)), by = .(cell_type, condition)]
elig <- dcast(elig, cell_type ~ condition, value.var = "mice", fill = 0)
for (cc in c("Control", "STZ")) if (!cc %in% names(elig)) elig[, (cc) := 0]
eligible <- elig[Control >= MIN_MICE & STZ >= MIN_MICE, cell_type]
message("Eligible cell types: ", paste(eligible, collapse = ", "))
fwrite(merge(dcast(ncell, cell_type ~ mouse_id, value.var = "N", fill = 0), elig, by = "cell_type"),
       file.path(out_dir, paste0(tissue, "_celltype_counts_per_mouse", sfx, ".csv")))

counts <- LayerData(so, assay = "RNA", layer = "counts")

run_ct <- function(ct) {
  cells_ct <- md[cell_type == ct]
  keep_m   <- cells_ct[, .N, by = mouse_id][N >= MIN_CELLS, mouse_id]
  cells_ct <- cells_ct[mouse_id %in% keep_m]
  grp <- factor(cells_ct$mouse_id)
  mm  <- sparse.model.matrix(~0 + grp); colnames(mm) <- levels(grp)
  pb  <- as.matrix(counts[, cells_ct$cell] %*% mm)
  cd  <- unique(cells_ct[, .(mouse_id, condition)])[match(colnames(pb), mouse_id)]
  cd[, n_cells := as.integer(table(grp)[mouse_id])]

  # edgeR QL
  y <- DGEList(pb, group = cd$condition)
  keep <- filterByExpr(y, group = cd$condition, min.count = 10)
  y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  design <- model.matrix(~condition, data = cd)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "conditionSTZ")
  tt <- as.data.table(topTags(qlf, n = Inf)$table, keep.rownames = "gene")

  # DESeq2 check
  dds <- DESeqDataSetFromMatrix(pb[keep, ], colData = as.data.frame(cd), design = ~condition)
  dds <- DESeq(dds, quiet = TRUE, fitType = "local")
  dr  <- as.data.table(as.data.frame(results(dds, name = "condition_STZ_vs_Control")), keep.rownames = "gene")

  tt <- merge(tt, dr[, .(gene, log2FC_DESeq2 = log2FoldChange, padj_DESeq2 = padj)], by = "gene", all.x = TRUE)
  # mean CPM per condition & detection rate
  cpm_m <- cpm(y, log = FALSE)
  tt[, cpm_Control := rowMeans(cpm_m[gene, cd$condition == "Control", drop = FALSE])]
  tt[, cpm_STZ     := rowMeans(cpm_m[gene, cd$condition == "STZ", drop = FALSE])]
  det <- rowMeans(counts[tt$gene, cells_ct$cell] > 0)
  tt[, pct_cells_detected := round(100 * det[gene], 2)]
  tt[, `:=`(cell_type = ct, n_mice_Control = sum(cd$condition == "Control"), n_mice_STZ = sum(cd$condition == "STZ"),
            n_cells = nrow(cells_ct))]
  tt[, sig := fifelse(FDR < 0.05 & logFC > LFC, "Up", fifelse(FDR < 0.05 & logFC < -LFC, "Down", "NS"))]
  tt[, deseq2_concordant := !is.na(padj_DESeq2) & padj_DESeq2 < 0.05 & sign(log2FC_DESeq2) == sign(logFC)]
  tt
}

de_all <- rbindlist(lapply(eligible, function(ct) {
  message("  DE: ", ct); tryCatch(run_ct(ct), error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
}), fill = TRUE)
setorder(de_all, cell_type, PValue)
fwrite(de_all, file.path(out_dir, paste0(tissue, "_pseudobulk_DE_all", sfx, ".csv")))

de_sum <- de_all[, .(genes_tested = .N, up = sum(sig == "Up"), down = sum(sig == "Down"),
                     sig_deseq2_concordant = sum(sig != "NS" & deseq2_concordant),
                     n_cells = n_cells[1], mice_ctrl = n_mice_Control[1], mice_stz = n_mice_STZ[1]),
                 by = cell_type][order(-(up + down))]
print(de_sum)
fwrite(de_sum, file.path(out_dir, paste0(tissue, "_pseudobulk_DE_summary", sfx, ".csv")))

p <- ggplot(melt(de_sum, id.vars = "cell_type", measure.vars = c("up", "down")),
            aes(reorder(cell_type, value), ifelse(variable == "up", value, -value), fill = variable)) +
  geom_col() + coord_flip() + scale_fill_manual(values = c(up = "#C8322F", down = "#3B7DD8")) +
  labs(x = NULL, y = "Pseudobulk DEGs (STZ vs Control)", title = paste(tissue, "- FDR<0.05, |log2FC|>0.5")) +
  theme_bw()
ggsave(file.path(fig_dir, paste0("09_", tissue, "_DE_counts", sfx, ".pdf")), p, width = 7, height = 5)

# ---------------- Differential abundance (propeller) ---------------------------------
da <- propeller(clusters = md$cell_type, sample = md$mouse_id, group = md$condition, transform = "logit")
da <- as.data.table(da, keep.rownames = "cell_type")
fwrite(da, file.path(out_dir, paste0(tissue, "_propeller_DA", sfx, ".csv")))
print(da)

prop <- md[, .N, by = .(mouse_id, condition, cell_type)][, frac := N / sum(N), by = mouse_id]
pc <- ggplot(prop, aes(mouse_id, frac, fill = cell_type)) + geom_col() +
  facet_grid(~condition, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Fraction of cells", title = paste(tissue, "cell-type composition")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, paste0("09_", tissue, "_composition", sfx, ".pdf")), pc, width = 8, height = 5)

writeLines(capture.output(sessionInfo()), file.path(out_dir, paste0("sessionInfo_09_", tissue, sfx, ".txt")))
