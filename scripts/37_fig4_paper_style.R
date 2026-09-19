# 37_fig4_paper_style.R - Figure 4 (machine-learning classifier) drawn in the visual style of the
# reference paper.
#
# Usage: Rscript scripts/37_fig4_paper_style.R <stage>
#   stage : A | B | C | D | E | panels | assemble | legend | all   (default all)
#
# Reference: Xu et al., Phytomedicine 154 (2026) 158050. Unlike Figures 1, 2, 3 and 5 of this project,
# this figure has NO panel-for-panel source: the reference paper runs no classifier, and its Figures 5
# to 11 are chemistry and wet-lab validation that this study does not have. What is matched here is the
# reference's visual language, so the shipped set reads as one figure set:
#   - boxed panels, ticks out, no grid, bold centred panel titles
#   - the reference's two-group colours, the same green and purple used in Figure 3
#   - a red-white-blue expression heatmap, as in the Figure 1 module-trait panel
#   - upper-case bold panel tags
#
# Panels, unchanged in content from script 30 - every number comes from what script 21 wrote:
#   A  CV AUC: panel-fixed vs nested vs the permutation null
#   B  ROC in the independent cohort GSE23343, AUC in a box inside the panel
#   C  LASSO selection frequency and random-forest importance per candidate gene
#   D  external predicted probability / score by true group
#   E  candidate-gene expression, z-scored within cohort, in both cohorts
#
# The point of the figure is that the optimistic estimate and the honest one disagree and that nothing
# generalises to the independent cohort. That is a negative result, and the legend says so plainly.
#
# This is the shipped Figure 4. It follows the reference paper's visual language, not the
# figure_style.R Nature/Cell contract that scripts 27-32 use; that set is no longer built.
#
# Outputs: figures/Fig4/Fig4<A-E>_*.pdf, Fig4.pdf, Fig4_legend.md
#          results/figure_exports/Fig4.tiff

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(pROC)
})
options(stringsAsFactors = FALSE)

args  <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) >= 1) args[1] else "all"
# Optional second argument selects a pathway-restricted ML run, e.g. "hsa04152" for the panel
# restricted to the selected KEGG map. It reads results/ml_<map>/ and writes its own figure folder.
# The shipped panel is the Fig 1 candidates that are nodes of the selected KEGG map. Pass "" as the
# second argument to draw the unrestricted 8-gene panel instead; it then gets its own folder.
MAP_DEFAULT <- "hsa04152"
PANEL_MAP <- if (length(args) >= 2) args[2] else MAP_DEFAULT
tag       <- if (nzchar(PANEL_MAP)) paste0("_", PANEL_MAP) else ""
dir_tag   <- if (identical(PANEL_MAP, MAP_DEFAULT)) "" else if (nzchar(PANEL_MAP)) tag else "_allcandidates"
set.seed(20260918)

fig_dir   <- paste0("figures/Fig4", dir_tag)
panel_dir <- file.path("results/figure_exports", paste0("fig4_panels_ref", dir_tag))
mld       <- paste0("results/ml", tag)
if (!dir.exists(mld))
  stop("missing ", mld, "
  run: Rscript scripts/21_ML_classifier.R ", PANEL_MAP, call. = FALSE)
bulk      <- "results/bulk"
for (dd in c(fig_dir, panel_dir)) dir.create(dd, recursive = TRUE, showWarnings = FALSE)

MODELS   <- c("LASSO", "RF", "SVM", "Score")
cols_grp <- c(Control = "#1B7837", T2D = "#762A83")          # as in Figure 3, paper style
cols_sch <- c(`Panel-fixed CV` = "#E8362B", `Nested CV` = "#3B7DD8", `Permuted labels` = "grey70")
cols_mod <- c(LASSO = "#3B7DD8", RF = "#E8362B", SVM = "#2E9E5B", Score = "#8E5CC8")
col_lo   <- "#3B7DD8"; col_mid <- "#FFFFFF"; col_hi <- "#E8362B"   # the Fig 1E ramp

theme_baseR <- function(base = 8)
  theme_bw(base_size = base) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
        axis.ticks = element_line(colour = "black", linewidth = 0.3),
        plot.title = element_text(hjust = 0.5, face = "bold", size = base + 0.5),
        strip.background = element_rect(fill = "grey92", colour = "black", linewidth = 0.4),
        strip.text = element_text(size = 6.5),
        plot.margin = margin(2, 3, 2, 2))

keep <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, units = "mm",
         device = cairo_pdf)
  saveRDS(list(plot = p, width = w, height = h), file.path(panel_dir, paste0(name, ".rds")))
  message("wrote ", name, " (", w, " x ", h, " mm)")
  invisible(p)
}

# ===================================================================================================
# A - cross-validated AUC under the two schemes, against the permutation null
# ===================================================================================================
if (stage %in% c("A", "panels", "all")) {
  pf <- fread(file.path(mld, "cv_auc_panel_fixed.csv"))[, .(model, auc, scheme = "Panel-fixed CV")]
  ns <- fread(file.path(mld, "cv_auc_nested.csv"))[, .(model, auc, scheme = "Nested CV")]
  nl <- fread(file.path(mld, "permutation_null.csv"))[, .(model = "Null", auc,
                                                          scheme = "Permuted labels")]
  d <- rbindlist(list(pf, ns, nl))
  d[, model := factor(model, levels = c("LASSO", "RF", "SVM", "Null"))]
  d[, scheme := factor(scheme, levels = names(cols_sch))]

  pA <- ggplot(d, aes(model, auc, fill = scheme)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3, colour = "grey30") +
    geom_boxplot(outlier.size = 0.25, linewidth = 0.25, width = 0.68,
                 position = position_dodge2(preserve = "single")) +
    scale_fill_manual(values = cols_sch, name = "Scheme") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(title = "Cross-validated AUC", x = NULL, y = "AUC") +
    theme_baseR() +
    theme(legend.position = "top", legend.title = element_text(size = 6),
          legend.text = element_text(size = 5.5), legend.key.size = unit(2.6, "mm"),
          legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -2, 0))
  keep(pA, "Fig4A_cv_auc", 95, 72)
}

# ===================================================================================================
# B - ROC in the independent cohort, with the AUC key boxed inside the panel
# ===================================================================================================
if (stage %in% c("B", "panels", "all")) {
  ep <- fread(file.path(mld, "external_predictions.csv"))
  em <- fread(file.path(mld, "external_metrics.csv"))
  roc_dt <- rbindlist(lapply(MODELS, function(m) {
    r <- suppressMessages(roc(ep$truth, ep[[m]], levels = c("Control", "T2D"), direction = "<"))
    data.table(model = m, fpr = 1 - r$specificities, tpr = r$sensitivities)
  }))
  roc_dt[, model := factor(model, levels = MODELS)]
  lab <- em[match(MODELS, model), sprintf("%s  %.2f", model, AUC)]

  pB <- ggplot(roc_dt, aes(fpr, tpr, colour = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey55") +
    geom_step(linewidth = 0.45) +
    scale_colour_manual(values = cols_mod, name = "AUC", labels = lab) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(title = "GSE23343 (independent)", x = "1 - specificity", y = "Sensitivity") +
    theme_baseR() +
    theme(legend.position = "inside", legend.position.inside = c(0.98, 0.02),
          legend.justification = c(1, 0),
          legend.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.2),
          legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
          legend.key.size = unit(2.6, "mm"), legend.margin = margin(1, 2, 1, 2))
  keep(pB, "Fig4B_external_roc", 78, 72)
}

# ===================================================================================================
# C - LASSO selection frequency and random-forest importance per gene
# ===================================================================================================
if (stage %in% c("C", "panels", "all")) {
  fr <- fread(file.path(mld, "feature_ranking.csv"))[order(consensus)]
  fr[, gene := factor(gene, levels = rev(gene))]
  d <- melt(fr, id.vars = "gene", measure.vars = c("lasso_freq", "rf_importance"),
            variable.name = "metric", value.name = "value")
  d[, metric := factor(metric, levels = c("lasso_freq", "rf_importance"),
                       labels = c("LASSO selection frequency", "RF permutation importance"))]

  pC <- ggplot(d, aes(value, gene, fill = metric)) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_col(width = 0.7, colour = "black", linewidth = 0.12) +
    facet_wrap(~ metric, scales = "free_x") +
    scale_fill_manual(values = setNames(c(cols_mod[["LASSO"]], cols_mod[["RF"]]), levels(d$metric)),
                      guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.08)), n.breaks = 4) +
    labs(x = NULL, y = NULL) +
    theme_baseR() +
    theme(axis.text.y = element_text(size = 6, face = "italic"),
          axis.text.x = element_text(size = 5.5), panel.spacing.x = unit(3, "mm"))
  keep(pC, "Fig4C_feature_ranking", 95, 72)
}

# ===================================================================================================
# D - external predicted probability / score by true group
# ===================================================================================================
if (stage %in% c("D", "panels", "all")) {
  ep <- fread(file.path(mld, "external_predictions.csv"))
  d  <- melt(ep, id.vars = c("GSM", "truth"), measure.vars = MODELS,
             variable.name = "model", value.name = "value")
  d[, truth := factor(truth, levels = names(cols_grp))]
  d[, model := factor(model, levels = MODELS)]

  pD <- ggplot(d, aes(truth, value, fill = truth)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.25, width = 0.62, alpha = 0.55) +
    geom_point(shape = 21, size = 0.8, stroke = 0.25, colour = "black",
               position = position_jitter(width = 0.12, height = 0, seed = 1)) +
    facet_wrap(~ model, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = cols_grp, name = "Group") +
    labs(title = "GSE23343 predictions", x = NULL, y = "Predicted probability / score") +
    theme_baseR() +
    theme(axis.text.x = element_text(size = 5.5, angle = 30, hjust = 1),
          legend.position = "top", legend.title = element_text(size = 6),
          legend.text = element_text(size = 5.5), legend.key.size = unit(2.6, "mm"),
          legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -2, 0))
  keep(pD, "Fig4D_external_predictions", 95, 72)
}

# ===================================================================================================
# E - candidate-gene expression in both cohorts, z-scored within cohort
# ===================================================================================================
if (stage %in% c("E", "panels", "all")) {
  fr   <- fread(file.path(mld, "feature_ranking.csv"))[order(consensus)]
  meta <- fread(file.path(bulk, "meta_bulk.csv"))
  zc <- function(e, genes) {
    e <- e[rownames(e) %in% genes, , drop = FALSE]
    t(scale(t(e)))
  }
  disc <- zc(readRDS(file.path(bulk, "expr_merged_combat.rds")), fr$gene)
  extf <- file.path(bulk, "..", "validation", "GSE23343_expr_gene.rds")
  extf <- if (file.exists(extf)) extf else file.path(bulk, "GSE23343_expr_gene.rds")
  ext  <- if (file.exists(extf)) zc(readRDS(extf), fr$gene) else NULL

  long <- function(m, cohort, grp) {
    d <- as.data.table(as.table(m))
    setnames(d, c("gene", "sample", "z"))
    d[, `:=`(cohort = cohort, group = grp[match(sample, names(grp))])]
    d
  }
  g1 <- setNames(meta$condition, meta$GSM)
  d  <- long(disc, "Discovery", g1)
  if (!is.null(ext)) {
    ev <- fread(file.path(mld, "external_predictions.csv"))
    g2 <- setNames(ev$truth, ev$GSM)
    d  <- rbind(d, long(ext, "GSE23343", g2))
  }
  d <- d[!is.na(group)]
  d[, gene := factor(gene, levels = rev(fr$gene))]
  d[, facet := factor(paste(cohort, group), levels = c("Discovery Control", "Discovery T2D",
                                                       "GSE23343 Control", "GSE23343 T2D"))]
  d <- d[!is.na(facet)]
  d[, z := pmax(pmin(z, 2), -2)]                 # clip, so a single outlier does not eat the ramp

  pE <- ggplot(d, aes(sample, gene, fill = z)) +
    geom_tile() +
    facet_wrap(~ facet, nrow = 1, scales = "free_x") +
    scale_fill_gradient2(low = col_lo, mid = col_mid, high = col_hi, midpoint = 0,
                         limits = c(-2, 2), name = "z",
                         guide = guide_colourbar(barwidth = unit(2.4, "mm"),
                                                 barheight = unit(16, "mm"))) +
    labs(x = NULL, y = NULL) +
    theme_baseR() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 6, face = "italic"),
          legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
          panel.spacing.x = unit(1.2, "mm"))
  keep(pE, "Fig4E_gene_heatmap", 190, 62)
}

# ===================================================================================================
# assemble
# ===================================================================================================
if (stage %in% c("assemble", "all")) {
  need <- c("Fig4A_cv_auc", "Fig4B_external_roc", "Fig4C_feature_ranking",
            "Fig4D_external_predictions", "Fig4E_gene_heatmap")
  ps <- lapply(need, function(n) {
    f <- file.path(panel_dir, paste0(n, ".rds"))
    if (!file.exists(f)) stop("missing panel ", n, "; run stage panels first")
    readRDS(f)$plot
  })
  design <- "
AAAAABBBB
AAAAABBBB
CCCCDDDDD
CCCCDDDDD
EEEEEEEEE
EEEEEEEEE
"
  comp <- ps[[1]] + ps[[2]] + ps[[3]] + ps[[4]] + ps[[5]] +
    plot_layout(design = design) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 13, face = "bold"))
  out <- file.path(fig_dir, paste0("Fig4", dir_tag))
  ggsave(paste0(out, ".pdf"), comp, width = 230, height = 205, units = "mm", device = cairo_pdf)
  tif <- file.path("results/figure_exports", paste0("Fig4", dir_tag, ".tiff"))
  unlink(tif)
  ggsave(tif, comp, width = 230, height = 205, units = "mm", dpi = 400, bg = "white",
         compression = "lzw")
  message("wrote ", out, ".pdf and ", tif)
}

# ===================================================================================================
if (stage %in% c("legend", "panels", "all")) {
  sm <- fread(file.path(mld, "ml_summary.csv"))
  em <- fread(file.path(mld, "external_metrics.csv"))
  fr <- fread(file.path(mld, "feature_ranking.csv"))[order(consensus)]
  nl <- fread(file.path(mld, "permutation_null.csv"))
  best <- em[which.max(AUC)]
  pf <- sm[analysis %like% "Panel-fixed" & model == "LASSO"]
  ns <- sm[analysis %like% "Nested" & model == "LASSO"]

  l <- c(
    "# Figure 4, reference-paper style - legend",
    "",
    "Drawn in the visual language of Xu et al., *Phytomedicine* 154 (2026) 158050.",
    "",
    "## This figure has no panel-for-panel source in the reference",
    "",
    "The reference paper runs no classifier, and its Figures 5 to 11 are chemistry and wet-lab validation that this study does not have. Figures 1, 2, 3 and 5 of this project each reproduce a reference figure panel for panel; this one reproduces only the reference's **visual language**, so that the shipped set reads as one figure set.",
    "",
    if (nzchar(PANEL_MAP))
      sprintf("**The panel is restricted to the selected KEGG map (%s).** Of the 8 Fig 1 candidate genes, %d are nodes of that map: %s. The other four (VSNL1, SREBF2, IGFBP1, IGFBP2) appear in no KEGG pathway at all and are excluded, so the classifier is tested on the same genes the pathway, single-cell and MR figures are about. The restriction is applied inside the nested-CV selection rule as well, not only to the fixed panel.",
              PANEL_MAP, nrow(fr), paste(sort(fr$gene), collapse = ", "))
    else
      "The panel is the 8-gene Fig 1 candidate set. Nothing about the analysis changes: every number is what script 21 wrote.",
    "",
    sprintf("**Fig. 4. Classification of T2D from the %s.** (A) Cross-validated AUC for LASSO, random forest and SVM under two schemes: panel-fixed CV, in which the gene panel selected on the full data is held constant across folds, and nested CV, in which selection is repeated inside each fold. Grey, AUC from %d label permutations; dashed line, chance. (B) ROC curves in the independent cohort GSE23343 (%d samples); the key gives AUC. Score is a parameter-free signature score. (C) LASSO selection frequency across folds and random-forest permutation importance for each candidate gene. (D) Predicted probability, or signature score, in GSE23343 by true group; points are samples. (E) Candidate-gene expression, z-scored within cohort and clipped at +/-2, in the discovery data and in GSE23343.",
            if (nzchar(PANEL_MAP)) sprintf("%d-gene %s panel", nrow(fr), PANEL_MAP)
            else "candidate panel",
            nrow(nl), nrow(fread(file.path(mld, "external_predictions.csv")))),
    "",
    "## Read this with the figure",
    "",
    sprintf("- **The two CV estimates in (A) disagree, and only the lower one is honest.** Panel-fixed CV gives LASSO a mean AUC of %.2f; nested CV, which repeats gene selection inside every fold, gives %.2f. The gap is selection bias, not model quality.",
            pf$mean_AUC, ns$mean_AUC),
    local({
      any_sep <- em[AUC_lo > 0.5, .N] > 0
      below   <- em[AUC < 0.5]
      sprintf("- **%s** The best external AUC is %.2f (%s), 95%% CI %.2f-%.2f%s.%s This is %s and should be reported as one.",
              if (any_sep) "One model separates the groups in the independent cohort."
              else "No model separates the groups in the independent cohort at conventional confidence.",
              best$AUC, best$model, best$AUC_lo, best$AUC_hi,
              if (any_sep) "" else ", and every confidence interval includes 0.5",
              if (nrow(below))
                sprintf(" %s is below chance at %s.",
                        paste(below$model, collapse = " and "),
                        paste(sprintf("%.2f", below$AUC), collapse = " and "))
              else "",
              if (any_sep) "a weak positive result" else "a negative result")
    }),
    sprintf("- **(C) ranks the genes the models lean on.** %s and %s are selected in %.0f%% and %.0f%% of LASSO folds. With %d features and %d discovery samples, that stability says the panel is small and internally consistent, not that these are validated markers.",
            fr$gene[1], fr$gene[2], 100 * fr$lasso_freq[1], 100 * fr$lasso_freq[2],
            nrow(fr), 27L),
    "- With 16 vs 11 samples in discovery and 10 vs 7 in validation, every interval on this figure is wide. Read it as a negative result, not as a ranking of models.",
    "- **(E) is z-scored within cohort**, so it shows relative pattern; absolute expression is not comparable across the two cohorts.",
    "",
    "## Statistics to quote",
    "",
    "| Quantity | Value |",
    "|---|---|",
    sprintf("| Panel-fixed CV AUC (optimistic) | %s |",
            paste(sm[analysis %like% "Panel-fixed", sprintf("%s %.2f", model, mean_AUC)],
                  collapse = ", ")),
    sprintf("| Nested CV AUC (honest) | %s |",
            paste(sm[analysis %like% "Nested", sprintf("%s %.2f", model, mean_AUC)],
                  collapse = ", ")),
    sprintf("| Permutation null AUC | %.2f (sd %.2f), %d permutations |", mean(nl$auc), sd(nl$auc),
            nrow(nl)),
    sprintf("| External AUC, GSE23343 | %s |",
            paste(em[, sprintf("%s %.2f", model, AUC)], collapse = ", ")),
    sprintf("| Best external model | %s, AUC %.2f (95%% CI %.2f-%.2f) |", best$model, best$AUC,
            best$AUC_lo, best$AUC_hi),
    sprintf("| Most stable features | %s |",
            paste(sprintf("%s (LASSO %.0f%%)", fr$gene[1:min(3, nrow(fr))],
                          100 * fr$lasso_freq[1:min(3, nrow(fr))]), collapse = ", ")),
    if (nzchar(PANEL_MAP) && file.exists("results/ml/external_metrics.csv")) {
      base <- fread("results/ml/external_metrics.csv")
      sprintf("| Same models on the unrestricted 8-gene panel | %s |",
              paste(base[, sprintf("%s %.2f", model, AUC)], collapse = ", "))
    } else NULL
  )
  l <- l[!vapply(l, is.null, TRUE)]
  writeLines(l, file.path(fig_dir, "Fig4_legend.md"))
  message("wrote ", file.path(fig_dir, "Fig4_legend.md"))
}
