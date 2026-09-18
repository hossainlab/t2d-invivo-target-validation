# 30_pub_fig4.R - Figure 4 (machine-learning classifier) redrawn to the publication contract.
#
# Usage: Rscript scripts/30_pub_fig4.R <stage>
#   stage : A | B | C | D | E | panels | assemble | legend | all
#
# Re-draws only; every number comes from what script 21 wrote:
#   a  CV AUC: panel-fixed vs nested vs permutation null  <- cv_auc_panel_fixed/nested, permutation_null
#   b  ROC in the independent cohort GSE23343             <- external_predictions.csv (+ external_metrics)
#   c  LASSO selection stability and RF importance        <- feature_ranking.csv
#   d  external predicted probability by true group       <- external_predictions.csv
#   e  candidate-gene expression in both cohorts          <- expr_merged_combat.rds, GSE23343_expr_gene.rds
#
# The point of this figure is that the optimistic estimate and the honest one disagree, so panel a is
# drawn so the three distributions can be compared directly rather than in three separate panels.
#
# Outputs: figures/Fig4_ML/Fig4<a-e>_*.pdf, Fig4_composite.pdf, Fig4_legend.md
#          results/figure_exports/Fig4_composite.tiff, results/figure_exports/fig4_panels/*.rds

source("scripts/figure_style.R")
suppressPackageStartupMessages({ library(pROC) })

args  <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) >= 1) args[1] else "all"
set.seed(20260918)

fig_dir   <- "figures/Fig4_ML"
panel_dir <- file.path(EXPORT_DIR, "fig4_panels")
ml        <- "results/ml"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
keep <- function(p, n, w, h) save_fig(p, n, w, h, dir = fig_dir, keep_dir = panel_dir)

MODELS <- c("LASSO", "RF", "SVM")
col_mod <- setNames(pal_qual[c(1, 3, 5)], MODELS)

# ===================================================================================================
# a - cross-validated AUC: optimistic, honest, and the null
# ===================================================================================================
if (stage %in% c("A", "panels", "all")) {
  pf <- fread(file.path(ml, "cv_auc_panel_fixed.csv"))[, .(model, auc, scheme = "Panel-fixed CV")]
  nx <- fread(file.path(ml, "cv_auc_nested.csv"))[, .(model, auc, scheme = "Nested CV")]
  nl <- fread(file.path(ml, "permutation_null.csv"))[, .(model = "Null", auc, scheme = "Permuted labels")]
  d  <- rbindlist(list(pf, nx, nl))
  d[, scheme := factor(scheme, levels = c("Panel-fixed CV", "Nested CV", "Permuted labels"))]
  d[, model := factor(model, levels = c(MODELS, "Null"))]

  pA <- ggplot(d, aes(model, auc, fill = scheme)) +
    geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.25, colour = "grey55") +
    geom_boxplot(outlier.size = 0.25, outlier.stroke = 0, linewidth = 0.25, width = 0.68,
                 position = position_dodge(preserve = "single")) +
    scale_fill_manual(values = c(`Panel-fixed CV` = "#E8B4A0", `Nested CV` = pal_qual[1],
                                 `Permuted labels` = "grey78"), name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = NULL, y = "Cross-validated AUC") +
    theme_pub() +
    theme(legend.position = "top", legend.justification = "left")
  keep(pA, "Fig4a_cv_auc", 89, 60)
}

# ===================================================================================================
# b - ROC in the independent cohort
# ===================================================================================================
if (stage %in% c("B", "panels", "all")) {
  pr <- fread(file.path(ml, "external_predictions.csv"))
  mx <- fread(file.path(ml, "external_metrics.csv"))
  truth <- factor(pr$truth)
  roc_d <- rbindlist(lapply(c(MODELS, "Score"), function(m) {
    r <- suppressMessages(pROC::roc(truth, pr[[m]], quiet = TRUE))
    data.table(model = m, fpr = 1 - r$specificities, tpr = r$sensitivities)
  }))
  auctab <- mx[order(-AUC), .(model, auclab = sprintf("%s  %.2f", model, AUC))]
  roc_d <- merge(roc_d, auctab, by = "model")
  roc_d[, auclab := factor(auclab, levels = auctab$auclab)]

  pB <- ggplot(roc_d[order(fpr)], aes(fpr, tpr, colour = auclab)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.25, colour = "grey55") +
    geom_step(linewidth = 0.45) +
    scale_colour_manual(values = setNames(pal_qual[c(1, 3, 5, 9)], auctab$auclab), name = "AUC") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
    coord_equal() +
    labs(x = "1 − specificity", y = "Sensitivity") +
    theme_pub() +
    theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
          legend.key.size = unit(2.2, "mm"))
  keep(pB, "Fig4b_external_roc", 70, 68)
}

# ===================================================================================================
# c - feature ranking: LASSO stability and RF importance
# ===================================================================================================
if (stage %in% c("C", "panels", "all")) {
  fr <- fread(file.path(ml, "feature_ranking.csv"))
  d <- melt(fr[, .(gene, `LASSO frequency` = lasso_freq, `RF importance` = rf_importance)],
            id.vars = "gene", variable.name = "metric", value.name = "value")
  ord <- fr[order(lasso_freq + rank(rf_importance) / 100), gene]
  d[, gene := factor(gene, levels = ord)]

  pC <- ggplot(d, aes(value, gene)) +
    geom_segment(aes(x = 0, xend = value, yend = gene), colour = "grey80", linewidth = 0.25) +
    geom_point(aes(colour = metric), size = 1.1) +
    facet_wrap(~ metric, scales = "free_x") +
    scale_colour_manual(values = setNames(pal_qual[c(1, 3)], levels(d$metric)), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.10))) +
    labs(x = NULL, y = NULL) +
    theme_pub() +
    theme(axis.text.y = element_text(face = "italic", size = BASE - 2),
          panel.spacing.x = unit(3, "mm"))
  keep(pC, "Fig4c_feature_ranking", 89, 78)
}

# ===================================================================================================
# d - predicted probability in the independent cohort, by true group
# ===================================================================================================
if (stage %in% c("D", "panels", "all")) {
  pr <- fread(file.path(ml, "external_predictions.csv"))
  d  <- melt(pr, id.vars = c("GSM", "truth"), measure.vars = c(MODELS, "Score"),
             variable.name = "model", value.name = "p")
  d[, truth := factor(truth)]

  pD <- ggplot(d, aes(truth, p, fill = truth)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.25, width = 0.62, alpha = 0.55) +
    geom_point(shape = 21, size = 0.9, stroke = 0.25, colour = "grey20",
               position = position_jitter(width = 0.14, height = 0)) +
    facet_wrap(~ model, nrow = 1) +
    scale_fill_manual(values = unname(col_grp[c("Control", "T2D")]), guide = "none") +
    labs(x = NULL, y = "Predicted probability / score") +
    theme_pub() +
    theme(panel.spacing.x = unit(2, "mm"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  keep(pD, "Fig4d_external_predictions", 89, 62)
}

# ===================================================================================================
# e - candidate-gene expression in both cohorts
# ===================================================================================================
if (stage %in% c("E", "panels", "all")) {
  fr <- fread(file.path(ml, "feature_ranking.csv"))
  genes <- fr[order(-consensus)][1:min(.N, 19), gene]

  zmat <- function(e, keep_genes) {
    g <- intersect(keep_genes, rownames(e))
    m <- t(scale(t(e[g, , drop = FALSE])))
    m[is.na(m)] <- 0
    m
  }
  e1 <- readRDS("results/bulk/expr_merged_combat.rds")
  m1 <- fread("results/bulk/meta_bulk.csv")[match(colnames(e1), GSM)]
  z1 <- zmat(e1, genes)
  d1 <- melt(as.data.table(z1, keep.rownames = "gene"), id.vars = "gene",
             variable.name = "GSM", value.name = "z")
  d1[, `:=`(group = m1$condition[match(GSM, m1$GSM)], cohort = "Discovery")]

  e2 <- readRDS("results/validation/GSE23343_expr_gene.rds")
  m2 <- fread("results/validation/GSE23343_meta.csv")
  gcol <- intersect(c("condition", "group", "T2D"), names(m2))[1]
  z2 <- zmat(e2, genes)
  d2 <- melt(as.data.table(z2, keep.rownames = "gene"), id.vars = "gene",
             variable.name = "GSM", value.name = "z")
  d2[, `:=`(group = as.character(m2[[gcol]][match(GSM, m2$GSM)]), cohort = "GSE23343")]

  d <- rbind(d1, d2)
  d[, group := fifelse(grepl("T2D|diab|1", group, ignore.case = TRUE), "T2D", "Control")]
  d[, gene := factor(gene, levels = rev(genes))]
  d[, z := pmax(pmin(z, 2), -2)]
  setorder(d, cohort, group, GSM)
  d[, GSM := factor(GSM, levels = unique(GSM))]

  pE <- ggplot(d, aes(GSM, gene, fill = z)) +
    geom_tile() +
    facet_grid(~ cohort + group, scales = "free_x", space = "free_x") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_gradientn(colours = pal_div, limits = c(-2, 2), breaks = c(-2, 0, 2),
                         name = "z",
                         guide = guide_colourbar(barheight = unit(16, "mm"),
                                                 barwidth = unit(2, "mm"))) +
    labs(x = NULL, y = NULL) +
    theme_pub_grid() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(face = "italic", size = BASE - 2),
          panel.spacing.x = unit(0.8, "mm"),
          strip.text = element_text(size = BASE - 2.5, margin = margin(0.5, 0.5, 1, 0.5)))
  keep(pE, "Fig4e_gene_heatmap", 120, 72)
}

# ===================================================================================================
# assemble
# ===================================================================================================
if (stage %in% c("assemble", "all")) {
  need <- c("Fig4a_cv_auc", "Fig4b_external_roc", "Fig4c_feature_ranking",
            "Fig4d_external_predictions", "Fig4e_gene_heatmap")
  ps <- lapply(need, function(n) {
    f <- file.path(panel_dir, paste0(n, ".rds"))
    if (!file.exists(f)) stop("missing panel ", n, "; run stage panels first")
    readRDS(f)$plot
  })
  design <- "
AAAAAAABBBBB
AAAAAAABBBBB
CCCCCDDDDDDD
CCCCCDDDDDDD
CCCCCEEEEEEE
CCCCCEEEEEEE
"
  comp <- ps[[1]] + ps[[2]] + ps[[3]] + ps[[4]] + ps[[5]] +
    plot_layout(design = design) +
    plot_annotation(tag_levels = list(tag_letters(5))) &
    theme(plot.tag = element_text(size = BASE + 1, face = "bold", family = FONT))
  save_fig(comp, "Fig4_composite", W_2COL, 180, dir = fig_dir, tiff = TRUE)
}

# ===================================================================================================
# legend
# ===================================================================================================
if (stage %in% c("legend", "panels", "all")) {
  pf <- fread(file.path(ml, "cv_auc_panel_fixed.csv"))
  nx <- fread(file.path(ml, "cv_auc_nested.csv"))
  nl <- fread(file.path(ml, "permutation_null.csv"))
  mx <- fread(file.path(ml, "external_metrics.csv"))
  fr <- fread(file.path(ml, "feature_ranking.csv"))
  med <- function(d) d[, .(m = median(auc)), by = model]

  l <- c(
    "# Figure 4 legend - classification of T2D from the candidate panel",
    "",
    "**Fig. 4 | A 19-gene candidate panel does not generalise to an independent cohort.**",
    "**a**, Cross-validated AUC for LASSO, random forest and SVM under two schemes: panel-fixed CV, in which the gene panel selected on the full data is held constant across folds, and nested CV, in which selection is repeated inside each fold. Grey, AUC from 200 label permutations. Boxes, median and interquartile range across CV repeats.",
    sprintf("**b**, ROC curves in the independent cohort GSE23343 (%d samples); the legend gives AUC. Score is a parameter-free signature score.",
            nrow(fread(file.path(ml, "external_predictions.csv")))),
    "**c**, LASSO selection frequency across folds and random-forest permutation importance for each candidate gene.",
    "**d**, Predicted probability (or signature score) in GSE23343 by true group; points are samples.",
    "**e**, Candidate-gene expression, z-scored within cohort, in the discovery data and in GSE23343.",
    "",
    "## Caveats that belong in the text, not the figure",
    "",
    sprintf("- The gap between panel-fixed and nested CV in **a** is the selection bias: panel-fixed median AUC %s vs nested %s (LASSO). Only the nested estimate is honest, and it is close to the permutation null.",
            mfmt(med(pf)[model == "LASSO", m]), mfmt(med(nx)[model == "LASSO", m])),
    sprintf("- No model separates groups in the independent cohort at conventional confidence: best AUC %s (%s), 95%% CI %s-%s, and every CI includes 0.5.",
            mfmt(mx[which.max(AUC), AUC]), mx[which.max(AUC), model],
            mfmt(mx[which.max(AUC), AUC_lo]), mfmt(mx[which.max(AUC), AUC_hi])),
    "- With 16 vs 11 samples in discovery and 10 vs 7 in validation, all of these intervals are wide; the figure should be read as a negative result, not a ranking of models.",
    "- **e** is z-scored within cohort, so it shows relative pattern, not comparable absolute expression.",
    "",
    "## Statistics to quote",
    "",
    "| Quantity | Value |",
    "|---|---|",
    sprintf("| Permutation null AUC (median) | %s |", mfmt(median(nl$auc))),
    sprintf("| Panel-fixed CV AUC, median | %s |",
            paste(sprintf("%s %s", med(pf)$model, mfmt(med(pf)$m)), collapse = "; ")),
    sprintf("| Nested CV AUC, median | %s |",
            paste(sprintf("%s %s", med(nx)$model, mfmt(med(nx)$m)), collapse = "; ")),
    sprintf("| External AUC (95%% CI) | %s |",
            paste(sprintf("%s %s (%s-%s)", mx$model, mfmt(mx$AUC), mfmt(mx$AUC_lo), mfmt(mx$AUC_hi)),
                  collapse = "; ")),
    sprintf("| Most stable LASSO feature | %s (selected in %s%% of folds) |",
            fr[which.max(lasso_freq), gene], mfmt(100 * fr[which.max(lasso_freq), lasso_freq], 0)),
    sprintf("| Top RF feature | %s |", fr[which.max(rf_importance), gene])
  )
  writeLines(l, file.path(fig_dir, "Fig4_legend.md"))
  message("wrote ", file.path(fig_dir, "Fig4_legend.md"))
}
