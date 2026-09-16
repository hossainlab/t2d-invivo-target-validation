# 11_external_validation_GSE23343.R — Phase 10 of docs/analysis_plan.md (human liver)
# Independent validation of bulk candidates in GSE23343 (human liver biopsies, T2D vs normal glucose
# tolerance; Affymetrix HG-U133 Plus 2.0, GPL570).
#
# Tests (discovery = T2D vs Control, GSE15653 + GSE64998):
#   1. limma T2D vs NGT for all genes; direction concordance of the 38 intersect genes and of all discovery DEGs
#      (binomial sign test vs 50%).
#   2. Per-gene AUC (pROC) for the 38 intersect genes, oriented by the discovery direction.
#   3. Signature score = mean z(up genes) - mean z(down genes); AUC and Wilcoxon test.
#
# Inputs:  results/bulk/intersect_genes.csv, results/bulk/DEG_T2D_vs_Control_all.csv
# Outputs: results/validation/GSE23343_meta.csv, GSE23343_limma_all.csv, GSE23343_candidates.csv,
#          GSE23343_summary.csv
# Figures: figures/validation/11_GSE23343_candidates_boxplot.pdf, 11_GSE23343_signature.pdf

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(limma)
  library(WGCNA)
  library(pROC)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
set.seed(20260914)
bulk_dir <- "results/bulk"; out_dir <- "results/validation"; fig_dir <- "results/supplementary_figures/validation"
geo_dir <- "data/validation"
for (d in c(out_dir, fig_dir, geo_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---------------- Load GSE23343 --------------------------------------------------------
gse <- getGEO("GSE23343", GSEMatrix = TRUE, destdir = geo_dir, AnnotGPL = FALSE)
eset <- gse[[1]]
pd <- as.data.table(pData(eset), keep.rownames = "GSM")
ch_cols <- grep("^characteristics|^title|^source_name|description", names(pd), value = TRUE)
txt <- tolower(apply(pd[, ..ch_cols], 1, paste, collapse = " | "))
pd[, group := fifelse(grepl("type 2 diabet|t2d|diabetic", txt) & !grepl("non-?diabet|normal glucose", txt), "T2D",
              fifelse(grepl("normal glucose|ngt|normoglyc|non-?diabet|control", txt), "NGT", NA_character_))]
print(pd[, .(GSM, title, group)])
stopifnot(!anyNA(pd$group), all(c("T2D", "NGT") %in% pd$group))
fwrite(pd[, c("GSM", "group", ch_cols), with = FALSE], file.path(out_dir, "GSE23343_meta.csv"))

# Expression: raw CEL + RMA (same pipeline as discovery) when data/validation/GSE23343_RAW.tar exists.
# The GEO series-matrix values for this study are pre-processed and centred near 0
# (range about -4 to 2), so they are used only as a fallback.
cel_tar <- file.path(geo_dir, "GSE23343_RAW.tar")
if (file.exists(cel_tar) && file.size(cel_tar) > 1e6) {
  suppressPackageStartupMessages({ library(oligo); library(hgu133plus2.db) })
  cel_dir <- file.path(geo_dir, "GSE23343_CEL"); dir.create(cel_dir, showWarnings = FALSE)
  untar(cel_tar, exdir = cel_dir)
  gz <- list.files(cel_dir, pattern = "\\.CEL\\.gz$", full.names = TRUE, ignore.case = TRUE)
  for (f in gz) {
    dest <- sub("\\.gz$", "", f, ignore.case = TRUE)
    if (!file.exists(dest)) {
      zin <- gzfile(f, "rb"); zout <- file(dest, "wb")
      repeat { buf <- readBin(zin, raw(), 1e7); if (!length(buf)) break; writeBin(buf, zout) }
      close(zin); close(zout)
    }
  }
  cels <- list.files(cel_dir, pattern = "\\.CEL$", full.names = TRUE, ignore.case = TRUE)
  raw <- read.celfiles(cels)
  sampleNames(raw) <- sub("^(GSM[0-9]+).*$", "\\1", basename(cels))
  ex <- exprs(rma(raw))
  sym <- suppressMessages(AnnotationDbi::mapIds(hgu133plus2.db, rownames(ex), "SYMBOL", "PROBEID", multiVals = "list"))
  keep <- lengths(sym) == 1 & !grepl("^AFFX", rownames(ex))
  sym1 <- unlist(sym[keep]); sym1 <- sym1[!is.na(sym1)]
  cr <- collapseRows(ex[names(sym1), ], rowGroup = unname(sym1), rowID = names(sym1), method = "MaxMean")
  expr_source <- "raw CEL, oligo RMA, hgu133plus2.db, MaxMean"
} else {
  ex <- exprs(eset)
  qx <- quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
  if (qx[5] > 100 || (qx[6] - qx[1] > 50 && qx[2] > 0)) {  # GEO2R log2 heuristic
    ex[ex <= 0] <- NA; ex <- log2(ex); message("Applied log2 transform")
  }
  ex <- normalizeBetweenArrays(ex, method = "quantile")
  fd <- as.data.table(fData(eset), keep.rownames = "probe")
  sym_col <- grep("^gene.?symbol$|^symbol$", names(fd), ignore.case = TRUE, value = TRUE)[1]
  stopifnot(!is.na(sym_col))
  fd[, symbol := trimws(sub(" ///.*$", "", get(sym_col)))]
  fd <- fd[symbol != "" & !is.na(symbol) & !grepl("^AFFX", probe)]
  ex <- ex[fd$probe, , drop = FALSE]
  ok <- rowSums(is.na(ex)) == 0
  cr <- collapseRows(ex[ok, ], rowGroup = fd$symbol[ok], rowID = fd$probe[ok], method = "MaxMean")
  expr_source <- "GEO series matrix (pre-processed values)"
}
gx <- cr$datETcollapsed
message("Expression source: ", expr_source)
pd <- pd[match(colnames(gx), GSM)]
message("GSE23343: ", nrow(gx), " genes x ", ncol(gx), " samples; ", paste(names(table(pd$group)), table(pd$group), collapse = ", "))

# ---------------- limma T2D vs NGT ----------------------------------------------------
pd[, group := factor(group, levels = c("NGT", "T2D"))]
fit <- eBayes(lmFit(gx, model.matrix(~group, data = pd)), trend = TRUE)
val <- data.table(gene = rownames(gx), topTable(fit, coef = "groupT2D", number = Inf, sort.by = "none"))
fwrite(val[order(P.Value)], file.path(out_dir, "GSE23343_limma_all.csv"))

# ---------------- Candidates ----------------------------------------------------------
disc <- fread(file.path(bulk_dir, "DEG_T2D_vs_Control_all.csv"))
ints <- fread(file.path(bulk_dir, "intersect_genes.csv"))
degs <- disc[deg_call != "NS"]

auc_gene <- function(g, dir) {
  if (!g %in% rownames(gx)) return(NA_real_)
  x <- gx[g, ]; if (dir < 0) x <- -x
  as.numeric(auc(roc(pd$group, x, levels = c("NGT", "T2D"), direction = "<", quiet = TRUE)))
}
cand <- merge(ints[, .(gene, logFC_discovery = logFC, P_discovery = P.Value, module, hub)],
              val[, .(gene, logFC_val = logFC, P_val = P.Value)], by = "gene", all.x = TRUE)
cand[, present := !is.na(logFC_val)]
cand[, concordant := present & sign(logFC_val) == sign(logFC_discovery)]
cand[, replicated := concordant & P_val < 0.05]
cand[, AUC_val := mapply(auc_gene, gene, sign(logFC_discovery))]
setorder(cand, -replicated, P_val)
fwrite(cand, file.path(out_dir, "GSE23343_candidates.csv"))
print(cand)

sign_test <- function(tab) {
  x <- merge(tab[, .(gene, d = logFC)], val[, .(gene, v = logFC)], by = "gene")
  k <- sum(sign(x$d) == sign(x$v)); n <- nrow(x)
  list(n = n, concordant = k, frac = k / n, p = if (n) binom.test(k, n, 0.5, alternative = "greater")$p.value else NA)
}
st_int <- sign_test(ints); st_deg <- sign_test(degs)

# ---------------- Signature score ----------------------------------------------------
zs <- t(scale(t(gx)))
sig_score <- function(up, down) {
  up <- intersect(up, rownames(zs)); down <- intersect(down, rownames(zs))
  (if (length(up)) colMeans(zs[up, , drop = FALSE]) else 0) - (if (length(down)) colMeans(zs[down, , drop = FALSE]) else 0)
}
pd[, score_intersect := sig_score(ints[logFC > 0, gene], ints[logFC < 0, gene])]
pd[, score_DEGs    := sig_score(degs[logFC > 0, gene], degs[logFC < 0, gene])]
sig_stats <- rbindlist(lapply(c("score_intersect", "score_DEGs"), function(s) {
  r <- roc(pd$group, pd[[s]], levels = c("NGT", "T2D"), direction = "<", quiet = TRUE)
  ci <- as.numeric(ci.auc(r))
  data.table(signature = s, AUC = as.numeric(auc(r)), AUC_CI_low = ci[1], AUC_CI_high = ci[3],
             wilcox_p = wilcox.test(pd[[s]] ~ pd$group)$p.value)
}))

summ <- data.table(
  item = c("samples_NGT", "samples_T2D", "genes",
           "intersect_present", "intersect_concordant", "intersect_sign_test_p", "intersect_replicated_P05",
           "DEGs_present", "DEGs_concordant", "DEGs_sign_test_p",
           "signature_intersect_AUC", "signature_intersect_AUC_CI", "signature_intersect_wilcox_p",
           "signature_DEGs_AUC", "signature_DEGs_AUC_CI", "signature_DEGs_wilcox_p"),
  value = c(sum(pd$group == "NGT"), sum(pd$group == "T2D"), nrow(gx),
            st_int$n, st_int$concordant, signif(st_int$p, 3), cand[replicated == TRUE, .N],
            st_deg$n, st_deg$concordant, signif(st_deg$p, 3),
            round(sig_stats$AUC[1], 3), sprintf("%.2f-%.2f", sig_stats$AUC_CI_low[1], sig_stats$AUC_CI_high[1]), signif(sig_stats$wilcox_p[1], 3),
            round(sig_stats$AUC[2], 3), sprintf("%.2f-%.2f", sig_stats$AUC_CI_low[2], sig_stats$AUC_CI_high[2]), signif(sig_stats$wilcox_p[2], 3))
)
summ <- rbind(data.table(item = "expression_source", value = expr_source), summ)
summ <- rbind(summ, data.table(item = "genomewide_logFC_spearman_vs_discovery",
  value = round(cor(merge(disc[, .(gene, d = logFC)], val[, .(gene, v = logFC)], by = "gene")[, .(d, v)], method = "spearman")[1, 2], 3)))
print(summ)
fwrite(summ, file.path(out_dir, "GSE23343_summary.csv"))

# ---------------- Figures ------------------------------------------------------------
show <- cand[present == TRUE][order(P_val)][1:min(.N, 20), gene]
long <- melt(data.table(GSM = colnames(gx), t(gx[show, , drop = FALSE])), id.vars = "GSM", variable.name = "gene", value.name = "expr")
long <- merge(long, pd[, .(GSM, group)], by = "GSM")
long[, gene := factor(gene, levels = show)]
pb <- ggplot(long, aes(group, expr, fill = group)) + geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.12, size = 1) + facet_wrap(~gene, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c(NGT = "#3B7DD8", T2D = "#C8322F"), guide = "none") +
  labs(x = NULL, y = "log2 expression", title = "GSE23343 (independent liver): candidate genes, T2D vs NGT") + theme_bw(base_size = 10)
ggsave(file.path(fig_dir, "11_GSE23343_candidates_boxplot.pdf"), pb, width = 11, height = 2.3 * ceiling(length(show) / 5) + 1)

ps <- lapply(c("score_intersect", "score_DEGs"), function(s) {
  a <- sig_stats[signature == s]
  ggplot(pd, aes(group, .data[[s]], fill = group)) + geom_boxplot(outlier.shape = NA, width = 0.6) + geom_jitter(width = 0.1) +
    scale_fill_manual(values = c(NGT = "#3B7DD8", T2D = "#C8322F"), guide = "none") +
    labs(x = NULL, y = "Signature score", title = s, subtitle = sprintf("AUC %.2f (%.2f-%.2f), Wilcoxon p = %.2g", a$AUC, a$AUC_CI_low, a$AUC_CI_high, a$wilcox_p)) +
    theme_bw()
})
ggsave(file.path(fig_dir, "11_GSE23343_signature.pdf"), wrap_plots(ps, nrow = 1), width = 9, height = 4.5)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_11.txt"))
