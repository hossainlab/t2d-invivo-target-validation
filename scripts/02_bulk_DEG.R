# 02_bulk_DEG.R — Phase 2 of docs/analysis_plan.md
# limma DEG on merged human liver arrays: T2D vs Control (user decision 2026-09-14).
#   T2D     = obese type 2 diabetic (GSE15653 Obese_DM well + poorly controlled; GSE64998 T2D-obese), n = 16
#   Control = lean (GSE15653) / non-obese metabolically healthy (GSE64998),                         n = 11
#   Obese non-diabetic samples are excluded in script 01 (user decision D1b); 27 samples.
# Primary model: ~condition + dataset on the uncorrected merged matrix.
# Checks: paper-style ComBat -> limma; + sex covariate; per-dataset.
# Thresholds as in the reference paper: |log2FC| > 0.5 and BH adj.P < 0.05.
# DEG definition for downstream steps (user decision M9): nominal P < 0.05 and |log2FC| > 0.5,
# because no gene reaches FDR < 0.05 for T2D vs Control. Reported as exploratory.
#
# Outputs: results/bulk/DEG_T2D_vs_Control_all.csv, DEG_summary.csv, DEG_concordance.csv,
#          HbA1c_assoc_GSE15653.csv
# Figures: figures/bulk/Fig1B_volcano_T2D_vs_Control.pdf/.tiff, 02_heatmap_top50_T2D_vs_Control.pdf

suppressPackageStartupMessages({
  library(limma)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
})

set.seed(20260914)
out_dir <- "results/bulk"; fig_dir <- "results/supplementary_figures/bulk"
fig_main <- "results/supplementary_figures/bulk"; exp_dir <- "results/figure_exports"  # draft panels;
# figures/Fig1 is owned by scripts/28_pub_fig1.R
for (d in c(fig_dir, fig_main, exp_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
LFC <- 0.5; FDR <- 0.05

expr_raw    <- readRDS(file.path(out_dir, "expr_merged_raw.rds"))
expr_combat <- readRDS(file.path(out_dir, "expr_merged_combat.rds"))
meta <- fread(file.path(out_dir, "meta_bulk.csv"))
meta <- meta[match(colnames(expr_raw), GSM)]
stopifnot(identical(meta$GSM, colnames(expr_raw)), identical(colnames(expr_combat), meta$GSM))
meta[, condition := factor(condition, levels = c("Control", "T2D"))]
meta[, group3 := factor(group3, levels = c("Lean", "ObeseT2D"))]  # obese non-diabetic excluded in 01 (D1b)
meta[, dataset := factor(dataset)]
meta[, sex := factor(sex)]
print(meta[, .N, by = .(dataset, condition, group3)][order(dataset, condition)])

fit_coef <- function(ex, design, coef) {
  fit <- eBayes(lmFit(ex, design), trend = TRUE)
  tt  <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  tt  <- data.table(gene = rownames(tt), tt)
  tt[, sig := fifelse(adj.P.Val < FDR & logFC > LFC, "Up",
              fifelse(adj.P.Val < FDR & logFC < -LFC, "Down", "NS"))]
  tt
}

# ---------------- Primary: ~condition + dataset -------------------------------------
res_primary <- fit_coef(expr_raw, model.matrix(~condition + dataset, data = meta), "conditionT2D")

# ---------------- Paper-style: ComBat (protected condition) -> ~condition -----------
res_paper <- fit_coef(expr_combat, model.matrix(~condition, data = meta), "conditionT2D")

# ---------------- Sensitivity: + sex ------------------------------------------------
res_sex <- fit_coef(expr_raw, model.matrix(~condition + dataset + sex, data = meta), "conditionT2D")

# ---------------- Sensitivity: per dataset ------------------------------------------
res_ds <- lapply(split(meta$GSM, meta$dataset), function(s) {
  md <- droplevels(meta[match(s, GSM)])
  fit_coef(expr_raw[, s], model.matrix(~condition, data = md), "conditionT2D")
})


# ---------------- DEG definition used downstream (user decision M6) ---------------------
P_DEG <- 0.05; LFC_DEG <- 0.5
res_primary[, deg_call := fifelse(P.Value < P_DEG & logFC > LFC_DEG, "Up",
                           fifelse(P.Value < P_DEG & logFC < -LFC_DEG, "Down", "NS"))]
message("T2D vs Control: ", res_primary[sig != "NS", .N], " genes at FDR<0.05 & |log2FC|>0.5; ",
        res_primary[deg_call != "NS", .N], " DEGs at nominal P<0.05 & |log2FC|>0.5 (used downstream)")

# ---------------- Merge result table --------------------------------------------------
tab <- copy(res_primary)
tab <- merge(tab, res_paper[, .(gene, logFC_paper = logFC, P_paper = P.Value, adjP_paper = adj.P.Val, sig_paper = sig)], by = "gene")
tab <- merge(tab, res_sex[, .(gene, logFC_sex = logFC, P_sex = P.Value, adjP_sex = adj.P.Val)], by = "gene")
tab <- merge(tab, res_ds$GSE15653[, .(gene, logFC_GSE15653 = logFC, P_GSE15653 = P.Value)], by = "gene")
tab <- merge(tab, res_ds$GSE64998[, .(gene, logFC_GSE64998 = logFC, P_GSE64998 = P.Value)], by = "gene")
tab[, consistent_both_datasets := sign(logFC_GSE15653) == sign(logFC) & sign(logFC_GSE64998) == sign(logFC)]
setorder(tab, P.Value)
fwrite(tab, file.path(out_dir, "DEG_T2D_vs_Control_all.csv"))

# ---------------- Summary & concordance ---------------------------------------------
count_sig <- function(x) c(up = x[sig == "Up", .N], down = x[sig == "Down", .N])
summ <- rbindlist(list(
  data.table(analysis = "primary ~condition+dataset", t(count_sig(res_primary))),
  data.table(analysis = "paper-style ComBat->limma",  t(count_sig(res_paper))),
  data.table(analysis = "sensitivity +sex",           t(count_sig(res_sex))),
  data.table(analysis = "GSE15653 only",              t(count_sig(res_ds$GSE15653))),
  data.table(analysis = "GSE64998 only",              t(count_sig(res_ds$GSE64998))),
  data.table(analysis = "DEG definition used downstream: P<0.05 & |log2FC|>0.5",
             up = res_primary[deg_call == "Up", .N], down = res_primary[deg_call == "Down", .N])
))
print(summ)
fwrite(summ, file.path(out_dir, "DEG_summary.csv"))

deg <- tab[deg_call != "NS"]
jac <- function(a, b) if (length(union(a, b))) length(intersect(a, b)) / length(union(a, b)) else NA_real_
conc <- data.table(
  n_primary_DEG = nrow(deg),
  jaccard_primary_vs_paper = jac(deg$gene, tab[P_paper < P_DEG & abs(logFC_paper) > LFC_DEG, gene]),
  frac_DEG_nominal_in_sex_model = if (nrow(deg)) mean(deg$P_sex < P_DEG & sign(deg$logFC_sex) == sign(deg$logFC)) else NA,
  sign_agree_GSE15653 = if (nrow(deg)) mean(sign(deg$logFC_GSE15653) == sign(deg$logFC)) else NA,
  sign_agree_GSE64998 = if (nrow(deg)) mean(sign(deg$logFC_GSE64998) == sign(deg$logFC)) else NA,
  cor_logFC_between_datasets_allgenes = cor(tab$logFC_GSE15653, tab$logFC_GSE64998)
)
print(conc)
fwrite(conc, file.path(out_dir, "DEG_concordance.csv"))

# ---------------- HbA1c association (GSE15653 only) ----------------------------------
m15 <- meta[dataset == "GSE15653" & !is.na(HbA1c)]
fit_h <- eBayes(lmFit(expr_raw[, m15$GSM], model.matrix(~HbA1c, data = m15)), trend = TRUE)
hb <- data.table(gene = rownames(expr_raw), topTable(fit_h, coef = "HbA1c", number = Inf, sort.by = "none"))
fwrite(hb[order(P.Value)], file.path(out_dir, "HbA1c_assoc_GSE15653.csv"))

# ---------------- Volcano (Fig 1B) ---------------------------------------------------
lab <- rbind(tab[deg_call == "Up"][1:min(.N, 12)], tab[deg_call == "Down"][1:min(.N, 12)])
g1 <- ggplot(tab, aes(logFC, -log10(P.Value), colour = deg_call)) +
  geom_point(size = 1.1, alpha = 0.7) +
  geom_vline(xintercept = c(-LFC_DEG, LFC_DEG), linetype = 2, colour = "grey50") +
  geom_hline(yintercept = -log10(P_DEG), linetype = 2, colour = "grey50") +
  geom_text_repel(data = lab, aes(label = gene), size = 3, max.overlaps = 40, show.legend = FALSE) +
  scale_colour_manual(values = c(Up = "#C8322F", Down = "#3B7DD8", NS = "grey75"),
                      labels = c(Up = paste0("Up (", tab[deg_call == "Up", .N], ")"),
                                 Down = paste0("Down (", tab[deg_call == "Down", .N], ")"), NS = "NS")) +
  labs(title = "T2D vs Control (human liver, GSE15653 + GSE64998)",
       subtitle = "limma ~condition + dataset; nominal P < 0.05 and |log2FC| > 0.5 (no gene at FDR < 0.05)",
       x = expression(log[2]~fold~change), y = expression(-log[10]~P), colour = NULL) +
  theme_bw(base_size = 12)
ggsave(file.path(fig_main, "Fig1B_volcano_T2D_vs_Control.pdf"), g1, width = 6.5, height = 5.5)
unlink(file.path(exp_dir, "Fig1B_volcano_T2D_vs_Control.tiff")); ggsave(file.path(exp_dir, "Fig1B_volcano_T2D_vs_Control.tiff"), g1, width = 6.5, height = 5.5, dpi = 300, compression = "lzw")

# ---------------- Heatmap top 25 up + 25 down (ComBat values) -----------------------
sig_col <- "deg_call"
top <- rbind(tab[get(sig_col) == "Up"][1:min(.N, 25)], tab[get(sig_col) == "Down"][1:min(.N, 25)])
if (nrow(top) >= 2) {
  s  <- meta[order(condition, dataset)]
  hm <- t(scale(t(expr_combat[top$gene, s$GSM])))
  ha <- HeatmapAnnotation(Condition = s$condition, Dataset = s$dataset,
          col = list(Condition = c(Control = "#3B7DD8", T2D = "#C8322F"),
                     Dataset = c(GSE15653 = "grey30", GSE64998 = "grey75")))
  pdf(file.path(fig_dir, "02_heatmap_top50_T2D_vs_Control.pdf"), width = 8, height = 9)
  draw(Heatmap(hm, name = "z", top_annotation = ha, cluster_columns = FALSE,
               column_split = s$condition, show_column_names = FALSE, row_names_gp = gpar(fontsize = 7),
               col = colorRamp2(c(-2, 0, 2), c("#3B7DD8", "white", "#C8322F"))))
  dev.off()
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_02.txt"))
