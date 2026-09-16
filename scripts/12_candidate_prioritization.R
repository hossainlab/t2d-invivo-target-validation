# 12_candidate_prioritization.R — Phase 12 of docs/analysis_plan.md
# Weighted evidence score for candidates (weights pre-specified in the plan, before scRNA results).
#
# Candidate pool:
#   A. Bulk main set (DEG ∩ WGCNA key modules; results/bulk/intersect_genes.csv) with mouse evidence
#      from liver and kidney (results/integration/<tissue>_candidate_concordance.csv, script 10).
#   B. scRNA-first genes: STZ vs Control pseudobulk DEGs (FDR < 0.05, |log2FC| > 0.5, DESeq2-concordant)
#      with a 1:1 human ortholog present in the bulk data, even if not in the bulk main set
#      (plan Tier 3; kidney candidates depend on this route).
#
# Criteria (score 0-3; weight):
#   bulk_DE (2)          P < 0.001 -> 3, < 0.01 -> 2, < 0.05 -> 1; 0 if direction differs between discovery datasets
#   sc_STZ (3)           concordant with bulk direction (pool A) or bulk logFC sign (pool B):
#                        FDR < 0.05 & DESeq2 concordant -> 3; FDR < 0.05 -> 2; P < 0.05 -> 1
#   specificity (1)      tau > 0.8 -> 3, > 0.6 -> 2, > 0.4 -> 1
#   network (1)          WGCNA hub -> 3; |kME| > 0.7 -> 2; key-module member -> 1
#   external_val (2)     GSE23343 concordant & P < 0.05 -> 3; concordant & P < 0.2 -> 2; concordant -> 1
#   clinical (1)         HbA1c association in GSE15653 in the T2D direction: P < 0.05 -> 3; P < 0.2 -> 1
#   protein (2), literature (1), druggability (1): manual scores from docs/annotation/candidate_manual_scores.csv
#                        (columns gene, protein_detectability, literature, druggability; 0-3); 0 + flag if absent
#   feasibility (1)      best-cell-type detection: > 30% cells -> 3, > 10% -> 2, > 3% -> 1
#
# Outputs: results/integration/candidate_scores_<tissue>.csv, candidate_scores_top.csv
# Figure:  figures/integration/12_candidate_score_heatmap_<tissue>.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
})
bulk_dir <- "results/bulk"; sc_dir <- "results/sc"; int_dir <- "results/integration"; val_dir <- "results/validation"
fig_dir <- "results/supplementary_figures/prioritization"; dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

weights <- c(bulk_DE = 2, sc_STZ = 3, specificity = 1, network = 1, external_val = 2, clinical = 1,
             protein = 2, literature = 1, druggability = 1, feasibility = 1)

deg  <- fread(file.path(bulk_dir, "DEG_T2D_vs_Control_all.csv"))
wg   <- fread(file.path(bulk_dir, "WGCNA_gene_module_GS_MM.csv"))
ints <- fread(file.path(bulk_dir, "intersect_genes.csv"))
hb   <- fread(file.path(bulk_dir, "HbA1c_assoc_GSE15653.csv"))[, .(gene, hb_logFC = logFC, hb_P = P.Value)]
val  <- if (file.exists(file.path(val_dir, "GSE23343_limma_all.csv")))
  fread(file.path(val_dir, "GSE23343_limma_all.csv"))[, .(gene, val_logFC = logFC, val_P = P.Value)] else data.table(gene = character())
man_file <- "docs/annotation/candidate_manual_scores.csv"
man <- if (file.exists(man_file)) fread(man_file) else NULL
# Key modules: same rule as 03/04 (user decision M11), not only modules of the few intersect genes
mt_key <- fread(file.path(bulk_dir, "WGCNA_module_trait.csv"))
key_modules <- mt_key[module != "grey" & p_T2D < 0.1 & abs(r_T2D) > abs(r_Dataset), module]

score_tissue <- function(tissue) {
  conc_file <- file.path(int_dir, paste0(tissue, "_candidate_concordance.csv"))
  pb_file   <- file.path(sc_dir, paste0(tissue, "_pseudobulk_DE_all.csv"))
  orth_file <- file.path(int_dir, paste0(tissue, "_orthologs.csv"))
  if (!all(file.exists(c(conc_file, pb_file)))) { message(tissue, ": integration outputs missing, skipped"); return(NULL) }
  conc <- fread(conc_file); pb <- fread(pb_file)

  # Pool B: scRNA-first genes with human ortholog in bulk universe
  o <- as.data.table(babelgene::orthologs(genes = unique(pb[sig != "NS", gene]), species = "mouse", human = FALSE))
  o <- unique(o[, .(human = human_symbol, mouse = symbol)])[, n_h := .N, by = human][, n_m := .N, by = mouse][n_h == 1 & n_m == 1]
  sc_first <- pb[sig != "NS" & deseq2_concordant == TRUE, .SD[order(FDR)][1], by = gene]
  sc_first <- merge(sc_first, o[, .(gene = mouse, human)], by = "gene")[human %in% deg$gene & !human %in% conc$gene]

  # Unified long table of mouse evidence (all cell types) for pools A and B
  pool <- unique(rbind(conc[!is.na(mouse), .(gene, mouse, pool = "A_bulk_main")],
                       sc_first[, .(gene = human, mouse = gene, pool = "B_scRNA_first")]))
  pool <- merge(pool, deg[, .(gene, logFC, P.Value, logFC_GSE15653, logFC_GSE64998)], by = "gene")
  ev <- merge(pool, pb[, .(mouse = gene, cell_type, logFC_m = logFC, FDR_m = FDR, P_m = PValue, sig_m = sig,
                           deseq2_concordant, pct_cells_detected)], by = "mouse", allow.cartesian = TRUE)
  ev[, concordant := sign(logFC_m) == sign(logFC)]
  ev[, sc_score := fifelse(concordant & sig_m != "NS" & deseq2_concordant, 3,
                   fifelse(concordant & sig_m != "NS", 2, fifelse(concordant & P_m < 0.05, 1, 0)))]
  best <- ev[order(-sc_score, FDR_m), .SD[1], by = gene][, .(gene, sc_STZ = sc_score, best_cell_type = cell_type,
                                                              best_logFC_mouse = logFC_m, best_FDR_mouse = FDR_m,
                                                              best_pct_detected = pct_cells_detected)]
  s <- merge(unique(pool[, .(gene, mouse, pool, logFC, P.Value, logFC_GSE15653, logFC_GSE64998)]), best, by = "gene", all.x = TRUE)
  s <- merge(s, conc[, .(gene, tau)], by = "gene", all.x = TRUE)
  s <- merge(s, wg[, .(gene, module, kME_own, hub)], by = "gene", all.x = TRUE)
  s <- merge(s, val, by = "gene", all.x = TRUE)
  s <- merge(s, hb, by = "gene", all.x = TRUE)

  s[, bulk_DE := fifelse(P.Value < 0.001, 3, fifelse(P.Value < 0.01, 2, fifelse(P.Value < 0.05, 1, 0)))]
  s[sign(logFC_GSE15653) != sign(logFC) | sign(logFC_GSE64998) != sign(logFC), bulk_DE := 0]
  s[is.na(sc_STZ), sc_STZ := 0]
  s[, specificity := fifelse(is.na(tau), 0, fifelse(tau > 0.8, 3, fifelse(tau > 0.6, 2, fifelse(tau > 0.4, 1, 0))))]
  s[, network := fifelse(hub %in% TRUE, 3, fifelse(module %in% key_modules & abs(kME_own) > 0.7, 2,
                 fifelse(module %in% key_modules, 1, 0)))]
  s[, external_val := fifelse(is.na(val_logFC) | sign(val_logFC) != sign(logFC), 0,
                      fifelse(val_P < 0.05, 3, fifelse(val_P < 0.2, 2, 1)))]
  s[, clinical := fifelse(is.na(hb_logFC) | sign(hb_logFC) != sign(logFC), 0, fifelse(hb_P < 0.05, 3, fifelse(hb_P < 0.2, 1, 0)))]
  s[, feasibility := fifelse(is.na(best_pct_detected), 0, fifelse(best_pct_detected > 30, 3,
                     fifelse(best_pct_detected > 10, 2, fifelse(best_pct_detected > 3, 1, 0))))]
  if (!is.null(man)) {
    s <- merge(s, man[, .(gene, protein = protein_detectability, literature, druggability)], by = "gene", all.x = TRUE)
  } else s[, `:=`(protein = NA_real_, literature = NA_real_, druggability = NA_real_)]
  s[, manual_scores_missing := is.na(protein) | is.na(literature) | is.na(druggability)]
  for (k in c("protein", "literature", "druggability")) s[is.na(get(k)), (k) := 0]

  crit <- names(weights)
  s[, score := as.numeric(as.matrix(.SD) %*% weights[crit]), .SDcols = crit]
  s[, max_score := sum(3 * weights)]
  s[, score_pct := round(100 * score / max_score, 1)]
  s[, tissue := tissue]
  setorder(s, -score)
  fwrite(s, file.path(int_dir, paste0("candidate_scores_", tissue, ".csv")))

  top <- s[1:min(.N, 25)]
  if (nrow(top) >= 2) {
    m <- as.matrix(top[, ..crit]); rownames(m) <- paste0(top$gene, " (", top$pool, ")")
    pdf(file.path(fig_dir, paste0("12_candidate_score_heatmap_", tissue, ".pdf")), width = 9, height = max(4, 0.28 * nrow(m) + 2))
    draw(Heatmap(m, name = "score", col = colorRamp2(c(0, 3), c("white", "#C8322F")), cluster_rows = FALSE, cluster_columns = FALSE,
                 column_labels = paste0(crit, " (w", weights[crit], ")"), row_names_gp = gpar(fontsize = 8),
                 right_annotation = rowAnnotation(total = anno_barplot(top$score)),
                 cell_fun = function(j, i, x, y, w, h, f) grid.text(m[i, j], x, y, gp = gpar(fontsize = 7)),
                 column_title = paste(tissue, "- candidate evidence scores")))
    dev.off()
  }
  message(tissue, ": ", nrow(s), " candidates scored (pool A ", s[pool == "A_bulk_main", .N], ", pool B ", s[pool == "B_scRNA_first", .N], ")")
  s
}

res <- rbindlist(lapply(c("liver", "kidney"), score_tissue), fill = TRUE)
if (nrow(res)) {
  top <- res[order(-score)][, head(.SD, 15), by = tissue]
  fwrite(top, file.path(int_dir, "candidate_scores_top.csv"))
  print(top[, .(tissue, gene, pool, score, score_pct, bulk_DE, sc_STZ, best_cell_type, specificity, network,
                external_val, clinical, feasibility, manual_scores_missing)])
  if (any(res$manual_scores_missing)) message("Manual scores (protein detectability, literature, druggability) missing for ",
                                              sum(res$manual_scores_missing), " rows; fill ", man_file, " and rerun.")
}
