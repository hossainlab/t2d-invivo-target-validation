# 01_bulk_preprocess.R — Phase 1 of docs/analysis_plan.md
# Human liver microarrays GSE15653 (GPL96, HG-U133A) + GSE64998 (GPL11532, HuGene 1.1 ST)
# Sample selection (D1b: 27 samples, T2D vs lean controls) -> raw CEL of selected samples -> RMA per platform
# -> gene-level -> merge on common genes (no expression filter, M8) -> ComBat -> QC/outliers
#
# Outputs (results/bulk/):
#   meta_bulk.csv, expr_<GSE>_gene.rds, expr_merged_raw.rds, expr_merged_combat.rds,
#   qc_outlier_flags.csv
# Figures: results/supplementary_figures/bulk/ 01_boxplots_rma.pdf, 01_rle_nuse.pdf, 01_sample_tree.pdf;
#          figures/Fig1_bulk_DEG_WGCNA/Fig1A_PCA_before_after_ComBat.pdf; TIFF in results/figure_exports/

suppressPackageStartupMessages({
  library(oligo)
  library(AnnotationDbi)
  library(hgu133a.db)
  library(hugene11sttranscriptcluster.db)
  library(WGCNA)
  library(sva)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260914)
proj    <- normalizePath(".")
out_dir <- file.path(proj, "results", "bulk")
fig_dir  <- file.path(proj, "results", "supplementary_figures", "bulk")  # QC plots
fig_main <- file.path(proj, "results", "supplementary_figures", "bulk")  # draft panels; figures/ is
# owned by scripts/28_pub_fig1.R, which redraws Fig 1 to the publication contract from the tables below.
# Writing here from an analysis script would clobber those panels (Windows filenames are case-insensitive).
exp_dir  <- file.path(proj, "results", "figure_exports")                  # TIFF exports
for (d in c(fig_main, exp_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- c("GSE15653", "GSE64998")

# ------------------------------------------------------------------
# 1. Sample metadata from GEO SOFT (data/meta/<GSE>_gsm_soft.txt)
# ------------------------------------------------------------------
parse_soft <- function(path) {
  lines  <- readLines(path, warn = FALSE)
  starts <- grep("^\\^SAMPLE = ", lines)
  ends   <- c(starts[-1] - 1, length(lines))
  rows <- lapply(seq_along(starts), function(i) {
    blk   <- lines[starts[i]:ends[i]]
    gsm   <- sub("^\\^SAMPLE = ", "", blk[1])
    title <- sub("^!Sample_title = ", "", grep("^!Sample_title = ", blk, value = TRUE))
    ch    <- sub("^!Sample_characteristics_ch1 = ", "",
                 grep("^!Sample_characteristics_ch1 = ", blk, value = TRUE))
    kv <- setNames(trimws(sub("^[^:]*:", "", ch)), trimws(sub(":.*$", "", ch)))
    data.frame(GSM = gsm, title = title, as.list(kv), check.names = FALSE)
  })
  rbindlist(rows, fill = TRUE)
}

assign_group <- function(title, gse) {
  t <- tolower(title)
  if (gse == "GSE15653") {
    g <- ifelse(grepl("^liver_lean", t), "Lean",
         ifelse(grepl("^liver_obese_nodm", t), "ObeseND",
         ifelse(grepl("^liver_obese_dm", t), "ObeseT2D", NA)))
  } else {
    g <- ifelse(grepl("^non-obese", t), "Lean",
         ifelse(grepl("^nd-obese", t), "ObeseND",
         ifelse(grepl("^t2d-obese", t), "ObeseT2D", NA)))
  }
  g
}

meta_list <- lapply(datasets, function(gse) {
  m <- parse_soft(file.path(proj, "data", "meta", paste0(gse, "_gsm_soft.txt")))
  num <- function(x) suppressWarnings(as.numeric(x))
  data.table(
    GSM       = m$GSM,
    title     = m$title,
    dataset   = gse,
    group3    = assign_group(m$title, gse),
    dm_control = ifelse(grepl("well-controlled", m$title), "well",
                 ifelse(grepl("poorly-controlled", m$title), "poor", NA)),
    sex       = tolower(m$gender),
    age       = if ("age" %in% names(m)) num(m$age) else NA_real_,
    BMI       = if ("body mass index (kg/m2)" %in% names(m)) num(m[["body mass index (kg/m2)"]]) else NA_real_,
    HbA1c     = if ("hemoglobin a1c" %in% names(m)) num(m[["hemoglobin a1c"]]) else NA_real_,
    fasting_glucose = if ("fasting glucose" %in% names(m)) num(m[["fasting glucose"]]) else NA_real_,
    fasting_insulin = if ("fasting insulin" %in% names(m)) num(m[["fasting insulin"]]) else NA_real_
  )
})
meta <- rbindlist(meta_list, fill = TRUE)
stopifnot(!anyNA(meta$group3))
meta[, group3 := factor(group3, levels = c("Lean", "ObeseND", "ObeseT2D"))]
meta[, T2D := as.integer(group3 == "ObeseT2D")]
meta[, obese := as.integer(group3 != "Lean")]
# Condition: T2D (ObeseT2D) vs Control (Lean). Obese non-diabetic samples are removed by the D1b selection below.
meta[, condition := factor(ifelse(group3 == "ObeseT2D", "T2D", "Control"), levels = c("Control", "T2D"))]
# Sample selection (user decision D1b, 2026-09-14): T2D vs lean/non-obese controls only; obese non-diabetic excluded
keep_gsm <- c(
  "GSM391693", "GSM391694", "GSM391695", "GSM391696", "GSM391697",                       # GSE15653 Control (lean)
  "GSM391702", "GSM391703", "GSM391704", "GSM391705", "GSM391706",                       # GSE15653 T2D (well-controlled)
  "GSM391707", "GSM391708", "GSM391709", "GSM391710",                                    # GSE15653 T2D (poorly controlled)
  "GSM1585585", "GSM1585586", "GSM1585587", "GSM1585588", "GSM1585589", "GSM1585590", "GSM1585591",  # GSE64998 T2D
  "GSM1585592", "GSM1585593", "GSM1585594", "GSM1585595", "GSM1585596", "GSM1585597")   # GSE64998 Control (non-obese)
stopifnot(all(keep_gsm %in% meta$GSM))
meta <- meta[GSM %in% keep_gsm]
meta[, group3 := droplevels(group3)]
stopifnot(nrow(meta) == 27, all(meta$group3 %in% c("Lean", "ObeseT2D")))
print(meta[, .N, by = .(dataset, condition, group3)])

# ------------------------------------------------------------------
# 2. Read CEL files (decompress .CEL.gz to a working folder once)
# ------------------------------------------------------------------
gunzip_to <- function(src, dest_dir) {
  dest <- file.path(dest_dir, sub("\\.gz$", "", basename(src)))
  if (!file.exists(dest)) {
    zin <- gzfile(src, "rb"); zout <- file(dest, "wb")
    repeat {
      buf <- readBin(zin, raw(), 1e7)
      if (!length(buf)) break
      writeBin(buf, zout)
    }
    close(zin); close(zout)
  }
  dest
}

read_cels <- function(gse) {
  gz  <- list.files(file.path(proj, "data", "bulk", paste0(gse, "_RAW")),
                    pattern = "\\.CEL\\.gz$", full.names = TRUE, ignore.case = TRUE)
  tmp <- file.path(proj, "data", "bulk", paste0(gse, "_CEL"))
  dir.create(tmp, showWarnings = FALSE)
  gz  <- gz[sub("^(GSM[0-9]+).*$", "\\1", basename(gz)) %in% meta$GSM]  # selected samples only (D1b)
  stopifnot(length(gz) == sum(meta$dataset == gse))
  cel <- vapply(gz, gunzip_to, character(1), dest_dir = tmp)
  raw <- read.celfiles(cel)
  sampleNames(raw) <- sub("^(GSM[0-9]+).*$", "\\1", basename(cel))
  raw
}

raw15653 <- read_cels("GSE15653")
raw64998 <- read_cels("GSE64998")

# ------------------------------------------------------------------
# 3. Array-level QC: RLE / NUSE from probe-level model
# ------------------------------------------------------------------
plm15653 <- fitProbeLevelModel(raw15653)
plm64998 <- fitProbeLevelModel(raw64998)

rle_nuse_stats <- function(plm, gse) {
  r <- RLE(plm, type = "values"); n <- NUSE(plm, type = "values")
  data.table(GSM = sub("^(GSM[0-9]+).*$", "\\1", colnames(r)), dataset = gse,
             RLE_median  = apply(r, 2, median, na.rm = TRUE),
             RLE_IQR     = apply(r, 2, IQR, na.rm = TRUE),
             NUSE_median = apply(n, 2, median, na.rm = TRUE))
}
arr_qc <- rbind(rle_nuse_stats(plm15653, "GSE15653"), rle_nuse_stats(plm64998, "GSE64998"))

pdf(file.path(fig_dir, "01_rle_nuse.pdf"), width = 12, height = 8)
par(mfrow = c(2, 2), mar = c(7, 4, 3, 1))
RLE(plm15653, las = 2, main = "GSE15653 RLE");  NUSE(plm15653, las = 2, main = "GSE15653 NUSE")
RLE(plm64998, las = 2, main = "GSE64998 RLE");  NUSE(plm64998, las = 2, main = "GSE64998 NUSE")
dev.off()

# ------------------------------------------------------------------
# 4. RMA + annotation + gene-level collapse (MaxMean)
# ------------------------------------------------------------------
eset15653 <- rma(raw15653)
eset64998 <- rma(raw64998, target = "core")

to_gene <- function(eset, annot_db, gse) {
  ex  <- exprs(eset)
  ids <- rownames(ex)
  sym <- suppressMessages(AnnotationDbi::mapIds(annot_db, keys = ids, column = "SYMBOL",
                                               keytype = "PROBEID", multiVals = "list"))
  n_sym <- lengths(sym)
  keep  <- n_sym == 1 & !grepl("^AFFX", ids)
  sym1  <- unlist(sym[keep])
  keep_ids <- names(sym1)
  ex <- ex[keep_ids, , drop = FALSE]
  sym1 <- sym1[!is.na(sym1)]
  ex <- ex[names(sym1), , drop = FALSE]
  cr <- collapseRows(ex, rowGroup = unname(sym1), rowID = names(sym1), method = "MaxMean")
  g  <- cr$datETcollapsed
  message(gse, ": ", length(ids), " probes -> ", length(sym1), " uniquely annotated -> ",
          nrow(g), " genes")
  g
}

g15653 <- to_gene(eset15653, hgu133a.db, "GSE15653")
g64998 <- to_gene(eset64998, hugene11sttranscriptcluster.db, "GSE64998")

pdf(file.path(fig_dir, "01_boxplots_rma.pdf"), width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(7, 4, 3, 1))
boxplot(g15653, las = 2, outline = FALSE, main = "GSE15653 RMA (gene level)")
boxplot(g64998, las = 2, outline = FALSE, main = "GSE64998 RMA (gene level)")
dev.off()

saveRDS(g15653, file.path(out_dir, "expr_GSE15653_gene.rds"))
saveRDS(g64998, file.path(out_dir, "expr_GSE64998_gene.rds"))

# ------------------------------------------------------------------
# 5. Expression filter: gene above the dataset's 25th percentile intensity in
#    >= smallest-group-size samples, required in BOTH datasets; then merge
# ------------------------------------------------------------------
expr_filter <- function(g, gse) {
  thr  <- quantile(g, 0.25)
  kmin <- min(table(meta[dataset == gse, condition]))
  keep <- rowSums(g > thr) >= kmin
  message(gse, ": threshold ", round(thr, 2), ", min samples ", kmin, ", kept ", sum(keep), "/", nrow(g))
  rownames(g)[keep]
}
common <- intersect(rownames(g15653), rownames(g64998))
# No expression filter (user decision M8; the reference paper and the earlier project do not filter).
# The filter used before removed 25 nominal DEGs, e.g. SLC16A4. expr_filter() is kept for reference only.
keep_genes <- common
message("Common genes kept (no expression filter): ", length(keep_genes))

expr_raw <- cbind(g15653[keep_genes, ], g64998[keep_genes, ])
meta <- meta[match(colnames(expr_raw), GSM)]
stopifnot(identical(meta$GSM, colnames(expr_raw)))

# ------------------------------------------------------------------
# 6. ComBat (batch = dataset, protect group3) — for PCA/WGCNA/heatmaps/ML only
# ------------------------------------------------------------------
combat_fun <- function(ex, md) {
  ComBat(dat = ex, batch = md$dataset, mod = model.matrix(~condition, data = md), par.prior = TRUE)
}
expr_combat <- combat_fun(expr_raw, meta)

# ------------------------------------------------------------------
# 7. Outlier detection (remove only if flagged by >= 2 of 3 criteria)
#    (a) PCA distance: within-dataset PC1-PC2 Euclidean distance > median + 3*MAD
#    (b) WGCNA standardized connectivity Z.k < -2.5 (on ComBat data)
#    (c) array QC: NUSE median > 1.05 or RLE IQR > dataset median + 3*MAD
# ------------------------------------------------------------------
pca_flag <- rbindlist(lapply(datasets, function(gse) {
  s  <- meta[dataset == gse, GSM]
  pc <- prcomp(t(expr_raw[, s]), scale. = FALSE)$x[, 1:2]
  d  <- sqrt(rowSums(sweep(pc, 2, colMeans(pc))^2))
  data.table(GSM = s, pca_dist = d, flag_pca = d > median(d) + 3 * mad(d))
}))

A    <- adjacency(expr_combat, type = "distance")
k    <- colSums(A) - 1
Zk   <- (k - mean(k)) / sd(k)
zk_flag <- data.table(GSM = names(Zk), Zk = Zk, flag_zk = Zk < -2.5)

arr_qc[, flag_array := NUSE_median > 1.05 |
         RLE_IQR > median(RLE_IQR) + 3 * mad(RLE_IQR), by = dataset]

flags <- Reduce(function(a, b) merge(a, b, by = "GSM"), list(pca_flag, zk_flag, arr_qc))
flags[, n_flags := flag_pca + flag_zk + flag_array]
flags <- merge(meta[, .(GSM, title, dataset, condition, group3)], flags, by = c("GSM", "dataset"))
fwrite(flags, file.path(out_dir, "qc_outlier_flags.csv"))
print(flags[n_flags > 0])

outliers <- flags[n_flags >= 2, GSM]
message("Outliers removed (>=2 flags): ", if (length(outliers)) paste(outliers, collapse = ", ") else "none")

hc <- hclust(dist(t(expr_combat)), method = "average")
pdf(file.path(fig_dir, "01_sample_tree.pdf"), width = 12, height = 6)
plotDendroAndColors(hc,
  colors = cbind(condition = labels2colors(as.integer(meta$condition)),
                 dataset = labels2colors(as.integer(factor(meta$dataset)) + 3),
                 Zk = ifelse(Zk[meta$GSM] < -2.5, "red", "white")),
  groupLabels = c("Condition", "Dataset", "Zk < -2.5"), main = "Sample clustering (ComBat)")
dev.off()

if (length(outliers)) {
  keep_s   <- setdiff(colnames(expr_raw), outliers)
  expr_raw <- expr_raw[, keep_s]
  meta     <- meta[match(keep_s, GSM)]
  expr_combat <- combat_fun(expr_raw, meta)
}
meta[, outlier_removed := FALSE]

# ------------------------------------------------------------------
# 8. PCA before / after ComBat (Fig 1A analogue)
# ------------------------------------------------------------------
pca_df <- function(ex, label) {
  p  <- prcomp(t(ex))
  ve <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
  d  <- data.table(GSM = rownames(p$x), PC1 = p$x[, 1], PC2 = p$x[, 2])
  d  <- merge(d, meta[, .(GSM, dataset, condition, group3)], by = "GSM")
  list(d = d, ve = ve, label = label)
}
plot_pca <- function(o) {
  ggplot(o$d, aes(PC1, PC2, colour = condition, shape = dataset)) +
    geom_point(size = 3, alpha = 0.9) +
    stat_ellipse(aes(group = dataset, linetype = dataset), colour = "grey40", level = 0.9) +
    scale_colour_manual(values = c(Control = "#3B7DD8", T2D = "#C8322F")) +
    labs(title = o$label, x = paste0("PC1 (", o$ve[1], "%)"), y = paste0("PC2 (", o$ve[2], "%)"),
         colour = "Condition", shape = "Dataset", linetype = "Dataset") +
    theme_bw(base_size = 12)
}
p_before <- plot_pca(pca_df(expr_raw, "Before batch correction"))
p_after  <- plot_pca(pca_df(expr_combat, "After ComBat"))
fig1a <- p_before + p_after + plot_layout(guides = "collect")
ggsave(file.path(fig_main, "Fig1A_PCA_before_after_ComBat.pdf"), fig1a, width = 11, height = 4.8)
unlink(file.path(exp_dir, "Fig1A_PCA_before_after_ComBat.tiff")); ggsave(file.path(exp_dir, "Fig1A_PCA_before_after_ComBat.tiff"), fig1a, width = 11, height = 4.8,
       dpi = 300, compression = "lzw")

# ------------------------------------------------------------------
# 9. Save
# ------------------------------------------------------------------
fwrite(meta, file.path(out_dir, "meta_bulk.csv"))
saveRDS(expr_raw,    file.path(out_dir, "expr_merged_raw.rds"))
saveRDS(expr_combat, file.path(out_dir, "expr_merged_combat.rds"))
fwrite(data.table(gene = rownames(expr_combat), expr_combat), file.path(out_dir, "expr_merged_combat.csv"))

message("Final: ", nrow(expr_raw), " genes x ", ncol(expr_raw), " samples")
print(meta[, .N, by = .(dataset, condition, group3)])
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_01.txt"))
