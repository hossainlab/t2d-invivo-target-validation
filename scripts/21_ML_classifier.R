# 21_ML_classifier.R — machine-learning evaluation of the candidate panel as a classifier of type 2 diabetes,
# and a consensus ranking of the genes (Fig 4).
#
# This step is NOT in the reference paper; it is an addition, and it is built so that the usual failure mode of
# small-n signature papers is visible rather than hidden:
#   1. Panel-fixed CV        — repeated stratified 5-fold CV on the 19-gene panel. This is what most papers
#                              report, and it is optimistic: the panel was selected on all 27 samples.
#   2. Honest nested CV      — the selection rule (limma DEG P < 0.05 and |logFC| > 0.5, intersected with the
#                              hyperglycemia gene set) is re-run inside every training fold, so no test sample
#                              ever influences feature choice. This is the defensible estimate.
#   3. Permutation null      — label-shuffled CV, to show what AUC this design produces by chance at n = 27.
#   4. External validation   — trained on all 27 discovery samples, tested on GSE23343 (10 T2D vs 7 NGT),
#                              an independent human liver cohort on a different platform.
# Cross-dataset prediction uses per-dataset gene-wise z-scoring, because discovery is ComBat-corrected
# GPL96 + HuGene and the external set is hgu133plus2 RMA.
#
# Models: LASSO logistic regression (glmnet), random forest (ranger), SVM with RBF kernel (e1071).
# Ranking: LASSO selection frequency, random-forest permutation importance, SVM-RFE rank.
#
# Outputs: results/ml/cv_auc_panel_fixed.csv, cv_auc_nested.csv, permutation_null.csv,
#          external_predictions.csv, external_metrics.csv, feature_ranking.csv, ml_summary.csv
# Figures: figures/Fig4_ML/Fig4A..F_*.pdf and the assembled figures/Fig4_ML/Fig4_ML.pdf

suppressPackageStartupMessages({
  library(data.table); library(glmnet); library(ranger); library(e1071); library(pROC)
  library(limma); library(ggplot2); library(patchwork); library(WGCNA)
})
set.seed(20260916)
out_dir <- "results/ml"; fig_dir <- "figures/Fig4_ML"
for (d in c(out_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
N_REPEATS <- 20L; K <- 5L; N_PERM <- 200L
cols_model <- c(LASSO = "#C8322F", RF = "#3B7DD8", SVM = "#E0A100")
cols_plot  <- c(cols_model, Score = "#2E7D5B")

# ---- 1. Discovery data -------------------------------------------------------------------
expr <- readRDS("results/bulk/expr_merged_combat.rds")
meta <- fread("results/bulk/meta_bulk.csv")[match(colnames(expr), GSM)]
stopifnot(identical(meta$GSM, colnames(expr)))
y <- factor(meta$condition, levels = c("Control", "T2D"))
panel <- fread("results/bulk/intersect_genes.csv")$gene
hyper <- fread("results/bulk/geneset_hyperglycemia.csv")$gene
message("Discovery: ", ncol(expr), " samples (", sum(y == "T2D"), " T2D / ", sum(y == "Control"),
        " lean); panel genes measured: ", length(intersect(panel, rownames(expr))))

# ---- 2. External cohort (GSE23343), cached ------------------------------------------------
ext_file <- "results/validation/GSE23343_expr_gene.rds"
if (file.exists(ext_file)) {
  gx <- readRDS(ext_file)
} else {
  suppressPackageStartupMessages({ library(oligo); library(hgu133plus2.db) })
  cel_dir <- "data/validation/GSE23343_CEL"
  cels <- list.files(cel_dir, pattern = "\\.CEL$", full.names = TRUE, ignore.case = TRUE)
  raw <- read.celfiles(cels)
  sampleNames(raw) <- sub("^(GSM[0-9]+).*$", "\\1", basename(cels))
  ex <- exprs(rma(raw))
  sym <- suppressMessages(AnnotationDbi::mapIds(hgu133plus2.db, rownames(ex), "SYMBOL", "PROBEID", multiVals = "list"))
  keep <- lengths(sym) == 1 & !grepl("^AFFX", rownames(ex))
  sym1 <- unlist(sym[keep]); sym1 <- sym1[!is.na(sym1)]
  cr <- collapseRows(ex[names(sym1), ], rowGroup = unname(sym1), rowID = names(sym1), method = "MaxMean")
  gx <- cr$datETcollapsed
  saveRDS(gx, ext_file)
}
ext_meta <- fread("results/validation/GSE23343_meta.csv")[match(colnames(gx), GSM)]
y_ext <- factor(ifelse(ext_meta$group == "T2D", "T2D", "Control"), levels = c("Control", "T2D"))
message("External GSE23343: ", ncol(gx), " samples (", sum(y_ext == "T2D"), " T2D / ", sum(y_ext == "Control"), " NGT)")

panel_use <- intersect(intersect(panel, rownames(expr)), rownames(gx))
message("Panel genes present in both cohorts: ", length(panel_use), " (", paste(panel_use, collapse = ", "), ")")

zscore <- function(m) t(scale(t(m)))          # gene-wise, within a dataset
Xd <- t(zscore(expr[panel_use, , drop = FALSE]))
Xe <- t(zscore(gx[panel_use, , drop = FALSE]))

# Direction-based signature score: mean z of the human-up genes minus mean z of the human-down genes.
# It has no fitted parameters, so it transports across platforms better than a trained classifier; the
# directions come from the discovery contrast, so its discovery AUC is optimistic and only the external
# value is a clean test.
deg_dir <- fread("results/bulk/DEG_T2D_vs_Control_all.csv")[gene %in% panel_use, .(gene, logFC)]
up_g <- deg_dir[logFC > 0, gene]; dn_g <- deg_dir[logFC < 0, gene]
sig_score <- function(X) rowMeans(X[, intersect(up_g, colnames(X)), drop = FALSE]) -
                         rowMeans(X[, intersect(dn_g, colnames(X)), drop = FALSE])
message("Signature score: ", length(up_g), " up genes, ", length(dn_g), " down genes")

# ---- 3. Model wrappers --------------------------------------------------------------------
fit_predict <- function(Xtr, ytr, Xte, model) {
  if (model == "LASSO") {
    nf <- min(5, min(table(ytr)))
    cv <- cv.glmnet(as.matrix(Xtr), ytr, family = "binomial", alpha = 1, nfolds = nf, type.measure = "deviance")
    as.numeric(predict(cv, newx = as.matrix(Xte), s = "lambda.min", type = "response"))
  } else if (model == "RF") {
    rf <- ranger(x = as.data.frame(Xtr), y = ytr, probability = TRUE, num.trees = 1000,
                 mtry = max(1, floor(sqrt(ncol(Xtr)))), min.node.size = 2)
    predict(rf, data = as.data.frame(Xte))$predictions[, "T2D"]
  } else {
    sv <- svm(x = as.matrix(Xtr), y = ytr, kernel = "radial", probability = TRUE, scale = FALSE)
    attr(predict(sv, as.matrix(Xte), probability = TRUE), "probabilities")[, "T2D"]
  }
}
folds_of <- function(y, k, seed) {
  set.seed(seed)
  idx <- integer(length(y))
  for (lvl in levels(y)) {
    w <- which(y == lvl); idx[w] <- sample(rep_len(seq_len(k), length(w)))
  }
  idx
}
auc_of <- function(obs, prob) as.numeric(pROC::auc(pROC::roc(obs, prob, levels = c("Control", "T2D"),
                                                             direction = "<", quiet = TRUE)))

# ---- 4. Panel-fixed repeated CV (the optimistic estimate) ---------------------------------
res_fixed <- list(); oof <- list()
for (r in seq_len(N_REPEATS)) {
  fo <- folds_of(y, K, 1000 + r)
  for (m in names(cols_model)) {
    p <- numeric(length(y))
    for (k in seq_len(K)) {
      tr <- fo != k; te <- fo == k
      p[te] <- fit_predict(Xd[tr, , drop = FALSE], y[tr], Xd[te, , drop = FALSE], m)
    }
    res_fixed[[length(res_fixed) + 1L]] <- data.table(repeat_id = r, model = m, auc = auc_of(y, p))
    if (r == 1L) oof[[m]] <- p
  }
}
cv_fixed <- rbindlist(res_fixed)
fwrite(cv_fixed, file.path(out_dir, "cv_auc_panel_fixed.csv"))
print(cv_fixed[, .(mean_AUC = round(mean(auc), 3), sd = round(sd(auc), 3),
                   lo = round(quantile(auc, 0.025), 3), hi = round(quantile(auc, 0.975), 3)), by = model])

# ---- 5. Honest nested CV: selection re-run inside each training fold ----------------------
select_in_fold <- function(idx_tr) {
  ex <- expr[, idx_tr, drop = FALSE]; md <- meta[idx_tr]
  des <- model.matrix(~ condition + dataset, data = md)
  fit <- eBayes(lmFit(ex, des), trend = TRUE)
  tt <- as.data.table(topTable(fit, coef = "conditionT2D", number = Inf), keep.rownames = "gene")
  g <- tt[P.Value < 0.05 & abs(logFC) > 0.5, gene]
  g <- intersect(g, hyper)
  if (length(g) < 3) g <- head(intersect(tt[order(P.Value), gene], hyper), 10)
  intersect(g, rownames(gx))                 # keep only genes that would also be portable
}
res_nested <- list()
for (r in seq_len(10L)) {
  fo <- folds_of(y, K, 5000 + r)
  for (m in names(cols_model)) {
    p <- numeric(length(y)); nsel <- integer(K)
    for (k in seq_len(K)) {
      tr <- which(fo != k); te <- which(fo == k)
      g <- select_in_fold(tr); nsel[k] <- length(g)
      if (length(g) < 2) { p[te] <- 0.5; next }
      Xtr <- t(zscore(expr[g, tr, drop = FALSE])); Xte <- t(zscore(expr[g, te, drop = FALSE]))
      p[te] <- fit_predict(Xtr, y[tr], Xte, m)
    }
    res_nested[[length(res_nested) + 1L]] <-
      data.table(repeat_id = r, model = m, auc = auc_of(y, p), mean_n_features = mean(nsel))
  }
}
cv_nested <- rbindlist(res_nested)
fwrite(cv_nested, file.path(out_dir, "cv_auc_nested.csv"))
print(cv_nested[, .(mean_AUC = round(mean(auc), 3), sd = round(sd(auc), 3),
                    mean_features = round(mean(mean_n_features), 1)), by = model])

# ---- 6. Permutation null (LASSO, panel fixed) ----------------------------------------------
perm <- vapply(seq_len(N_PERM), function(i) {
  set.seed(20000 + i); yp <- factor(sample(as.character(y)), levels = levels(y))
  fo <- folds_of(yp, K, 30000 + i); p <- numeric(length(yp))
  for (k in seq_len(K)) {
    tr <- fo != k; te <- fo == k
    p[te] <- tryCatch(fit_predict(Xd[tr, , drop = FALSE], yp[tr], Xd[te, , drop = FALSE], "LASSO"),
                      error = function(e) rep(0.5, sum(te)))
  }
  auc_of(yp, p)
}, 0)
fwrite(data.table(auc = perm), file.path(out_dir, "permutation_null.csv"))
obs_lasso <- cv_fixed[model == "LASSO", mean(auc)]
perm_p <- (sum(perm >= obs_lasso) + 1) / (N_PERM + 1)
message("Permutation null: mean AUC ", round(mean(perm), 3), ", 95th pct ", round(quantile(perm, 0.95), 3),
        "; observed ", round(obs_lasso, 3), ", empirical P ", signif(perm_p, 3))

# ---- 7. External validation ------------------------------------------------------------------
ext_pred <- data.table(GSM = rownames(Xe), truth = y_ext)
ext_metrics <- list()
for (m in names(cols_model)) {
  set.seed(77)
  pr <- fit_predict(Xd, y, Xe, m)
  ext_pred[[m]] <- pr
  ro <- pROC::roc(y_ext, pr, levels = c("Control", "T2D"), direction = "<", quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(ro))
  cut <- 0.5
  tp <- sum(pr >= cut & y_ext == "T2D"); tn <- sum(pr < cut & y_ext == "Control")
  ext_metrics[[m]] <- data.table(model = m, AUC = as.numeric(pROC::auc(ro)), AUC_lo = ci[1], AUC_hi = ci[3],
                                 accuracy = (tp + tn) / length(y_ext),
                                 sensitivity = tp / sum(y_ext == "T2D"),
                                 specificity = tn / sum(y_ext == "Control"),
                                 p_vs_chance = pROC::roc.test(ro, pROC::roc(y_ext, rep(0.5, length(y_ext)),
                                   levels = c("Control", "T2D"), direction = "<", quiet = TRUE))$p.value)
}
pr_s <- sig_score(Xe)
ext_pred[["Score"]] <- pr_s
ro_s <- pROC::roc(y_ext, pr_s, levels = c("Control", "T2D"), direction = "<", quiet = TRUE)
ci_s <- as.numeric(pROC::ci.auc(ro_s))
ext_metrics[["Score"]] <- data.table(model = "Score", AUC = as.numeric(pROC::auc(ro_s)),
  AUC_lo = ci_s[1], AUC_hi = ci_s[3],
  accuracy = mean((pr_s > median(sig_score(Xd))) == (y_ext == "T2D")),
  sensitivity = mean(pr_s[y_ext == "T2D"] > median(sig_score(Xd))),
  specificity = mean(pr_s[y_ext == "Control"] <= median(sig_score(Xd))),
  p_vs_chance = NA_real_)
ext_metrics <- rbindlist(ext_metrics)
message("Discovery AUC of the signature score (optimistic, directions from these samples): ",
        round(auc_of(y, sig_score(Xd)), 3))
fwrite(ext_pred, file.path(out_dir, "external_predictions.csv"))
fwrite(ext_metrics, file.path(out_dir, "external_metrics.csv"))
print(ext_metrics[, .(model, AUC = round(AUC, 3), lo = round(AUC_lo, 3), hi = round(AUC_hi, 3),
                      acc = round(accuracy, 2), sens = round(sensitivity, 2), spec = round(specificity, 2))])

# ---- 8. Feature ranking ------------------------------------------------------------------------
sel_freq <- setNames(numeric(length(panel_use)), panel_use)
for (r in seq_len(N_REPEATS)) {
  fo <- folds_of(y, K, 1000 + r)
  for (k in seq_len(K)) {
    tr <- fo != k
    cv <- cv.glmnet(as.matrix(Xd[tr, , drop = FALSE]), y[tr], family = "binomial", alpha = 1,
                    nfolds = min(5, min(table(y[tr]))))
    co <- as.matrix(coef(cv, s = "lambda.min"))[-1, 1]
    sel_freq[names(co)[co != 0]] <- sel_freq[names(co)[co != 0]] + 1
  }
}
sel_freq <- sel_freq / (N_REPEATS * K)
rf_full <- ranger(x = as.data.frame(Xd), y = y, probability = TRUE, num.trees = 5000,
                  importance = "permutation", mtry = max(1, floor(sqrt(ncol(Xd)))), min.node.size = 2)
rf_imp <- ranger::importance(rf_full)
rfe_rank <- local({                                  # SVM-RFE: drop the smallest-weight gene each round
  keep <- panel_use; rank <- integer(0)
  while (length(keep) > 1) {
    sv <- svm(x = as.matrix(Xd[, keep, drop = FALSE]), y = y, kernel = "linear", scale = FALSE)
    w <- abs(drop(t(sv$coefs) %*% sv$SV))
    worst <- keep[which.min(w)]
    rank <- c(worst, rank); keep <- setdiff(keep, worst)
  }
  r <- c(keep, rank)                                  # best first
  setNames(seq_along(r), r)
})
rankdt <- data.table(gene = panel_use, lasso_freq = sel_freq[panel_use],
                     rf_importance = rf_imp[panel_use], svmrfe_rank = rfe_rank[panel_use])
rankdt[, consensus := (frank(-lasso_freq) + frank(-rf_importance) + frank(svmrfe_rank)) / 3]
setorder(rankdt, consensus)
fwrite(rankdt, file.path(out_dir, "feature_ranking.csv"))
print(rankdt)

# ---- 9. Figure 4 ---------------------------------------------------------------------------------
cvA <- melt(rbind(cv_fixed[, .(model, auc, design = "Panel-fixed CV")],
                  cv_nested[, .(model, auc, design = "Nested CV (selection inside folds)")],
                  data.table(model = "LASSO", auc = perm, design = "Permutation null")),
            id.vars = c("model", "design"), measure.vars = "auc")
cvA[, design := factor(design, levels = c("Panel-fixed CV", "Nested CV (selection inside folds)", "Permutation null"))]
pA <- ggplot(cvA, aes(design, value, fill = model)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
  geom_boxplot(outlier.size = 0.4, position = position_dodge(preserve = "single")) +
  scale_fill_manual(values = cols_model, name = NULL) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "AUC", title = "A  Cross-validated AUC: optimistic vs honest vs chance") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 12, hjust = 1))

rocs <- lapply(names(cols_plot), function(m) {
  ro <- pROC::roc(y_ext, ext_pred[[m]], levels = c("Control", "T2D"), direction = "<", quiet = TRUE)
  data.table(model = m, fpr = 1 - ro$specificities, tpr = ro$sensitivities)
})
pB <- ggplot(rbindlist(rocs), aes(fpr, tpr, colour = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_step(linewidth = 0.8) + scale_colour_manual(values = cols_plot, name = NULL) +
  labs(x = "1 - specificity", y = "sensitivity",
       title = "B  External cohort GSE23343 (10 T2D / 7 NGT)",
       subtitle = paste(ext_metrics[, sprintf("%s AUC %.2f", model, AUC)], collapse = "   ")) +
  theme_bw(base_size = 9)

rk <- copy(rankdt)[, gene := factor(gene, levels = rev(gene))]
pC <- ggplot(rk, aes(lasso_freq, gene)) + geom_col(fill = "#C8322F") +
  labs(x = "LASSO selection frequency (100 folds)", y = NULL, title = "C  Stability of gene selection") +
  theme_bw(base_size = 9)
pD <- ggplot(rk, aes(rf_importance, gene)) + geom_col(fill = "#3B7DD8") +
  labs(x = "random-forest permutation importance", y = NULL, title = "D  Random-forest importance") +
  theme_bw(base_size = 9)

epred <- melt(ext_pred, id.vars = c("GSM", "truth"), variable.name = "model", value.name = "prob")
pE <- ggplot(epred, aes(truth, prob, fill = truth)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_point(shape = 21, size = 1.6, position = position_jitter(width = 0.12)) +
  facet_wrap(~ model) + scale_fill_manual(values = c(Control = "#3B7DD8", T2D = "#C8322F"), guide = "none") +
  labs(x = NULL, y = "predicted P(T2D)  /  signature score", title = "E  External predictions by true group") +
  theme_bw(base_size = 9)

top_g <- head(rankdt$gene, 8)
hm <- rbind(
  data.table(cohort = "Discovery", melt(data.table(sample = rownames(Xd), group = as.character(y),
                                                   Xd[, top_g, drop = FALSE]),
                                        id.vars = c("sample", "group"), variable.name = "gene", value.name = "z")),
  data.table(cohort = "GSE23343", melt(data.table(sample = rownames(Xe), group = as.character(y_ext),
                                                  Xe[, top_g, drop = FALSE]),
                                       id.vars = c("sample", "group"), variable.name = "gene", value.name = "z")))
hm[, sample := factor(sample, levels = unique(sample[order(cohort, group)]))]
pF <- ggplot(hm, aes(sample, gene, fill = pmin(pmax(z, -2), 2))) + geom_tile() +
  facet_grid(. ~ cohort + group, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#3B7DD8", mid = "white", high = "#C8322F", midpoint = 0, name = "z") +
  labs(x = NULL, y = NULL, title = "F  Top consensus genes in both cohorts") +
  theme_bw(base_size = 8) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

for (nm in c("A", "B", "C", "D", "E", "F")) {
  p <- get(paste0("p", nm))
  ggsave(file.path(fig_dir, paste0("Fig4", nm, "_", c(A = "cv_auc", B = "external_roc", C = "lasso_stability",
         D = "rf_importance", E = "external_predictions", F = "gene_heatmap")[[nm]], ".pdf")),
         p, width = if (nm %in% c("A", "F")) 9 else 6.5, height = if (nm == "F") 4 else 4.5)
}
fig4 <- (pA | pB) / (pC | pD) / (pE | pF) + plot_annotation(
  title = "Fig 4  Machine-learning evaluation of the candidate panel",
  subtitle = sprintf("Discovery: %d samples (%d T2D / %d lean), %d-gene panel. External: GSE23343 (%d T2D / %d NGT). Permutation P = %s",
                     ncol(expr), sum(y == "T2D"), sum(y == "Control"), length(panel_use),
                     sum(y_ext == "T2D"), sum(y_ext == "Control"), signif(perm_p, 3)))
ggsave(file.path(fig_dir, "Fig4_ML.pdf"), fig4, width = 15, height = 15)

summ <- rbind(
  cv_fixed[, .(analysis = "Panel-fixed CV (optimistic)", mean_AUC = mean(auc), sd = sd(auc)), by = model],
  cv_nested[, .(analysis = "Nested CV (honest)", mean_AUC = mean(auc), sd = sd(auc)), by = model],
  ext_metrics[, .(model, analysis = "External GSE23343", mean_AUC = AUC, sd = NA_real_)],
  data.table(model = "LASSO", analysis = "Permutation null", mean_AUC = mean(perm), sd = sd(perm)))
summ[, permutation_p := ifelse(analysis == "Panel-fixed CV (optimistic)" & model == "LASSO", perm_p, NA_real_)]
fwrite(summ, file.path(out_dir, "ml_summary.csv"))
print(summ)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_21.txt"))
