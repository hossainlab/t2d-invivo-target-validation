# 24_target_scoring.R — assemble the evidence per gene and assign target tiers (docs/target_selection.md).
#
# Nothing new is computed here: every column comes from a stage that was run and reported earlier. The point
# is that the tier rule is applied mechanically to a table, so the target list is reproducible and can be
# re-derived if any upstream stage changes.
#
# Evidence axes (columns):
#   human_logFC/P/FDR   discovery contrast, obese T2D vs lean liver (script 02)
#   module              member of a WGCNA module associated with T2D (p < 0.1; script 03)
#   phenotype_set       member of the curated hyperglycemia / insulin-resistance set (script 04)
#   candidate19         the Fig 1 three-way intersection (script 04)
#   on_map / axis_anchor  on the selected KEGG map hsa04152 / node of the candidate-anchored axis (script 16)
#   sc_cell_type, sc_P  cell type where the mouse ortholog moves in the human direction (scripts 17, 18, 23)
#   ml_rank             consensus machine-learning rank over the 19-gene panel (script 21)
#   mr_OR, mr_p, mr_FDR LD-clumped cis-MR on type 2 diabetes (script 20b)
#   replicated, val_P   independent human liver cohort GSE23343 (script 11)
#   protein_readout     whether the mechanism is visible at transcript level or needs protein/phospho assay
#
# Tier rule, fixed before reading the table:
#   Tier 1  in the Fig 1 candidate set AND on the selected map AND confirmed in the mouse atlas
#           AND (independently replicated in human OR top-5 of the machine-learning consensus)
#   Tier 2  confirmed in the mouse atlas and human-DE, but missing one of the Tier 1 anchors
#           (used as readout / biomarker rather than as a target)
#   Tier 2b candidate on the selected map whose mouse ortholog cannot be tested here, because the cell type
#           that expresses it (hepatocyte, proximal tubule) is absent from the atlas - measure in whole tissue
#   Tier 3  mechanistically required node that transcript data cannot test (phospho-regulated), carried on
#           pathway logic plus literature, to be measured at protein level
#   Not a target  evidence from a single axis only (classifier-only, or MR-only with a direction that
#           contradicts the tissue biology)
#
# Outputs: results/targets/target_tiers.csv, results/targets/target_evidence_long.csv

suppressPackageStartupMessages({ library(data.table) })
out_dir <- "results/targets"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

deg  <- fread("results/bulk/DEG_T2D_vs_Control_all.csv")
cons <- fread("results/consistency/gene_stage_matrix.csv")
ml   <- fread("results/ml/feature_ranking.csv")
mr   <- fread("results/mr/mr_results_blood_ldclumped.csv")[method %in% c("IVW", "Wald ratio")]
val  <- fread("results/validation/GSE23343_candidates.csv")
anch <- fread("results/pathway_selection/selection_axis_anchored.csv")[pathway == "hsa04152"]
sc23 <- rbindlist(lapply(c("liver", "kidney"), function(t)
  fread(paste0("results/pathway_selection/23_anchored_gene_stats_", t, ".csv"))[, tissue := t]), fill = TRUE)
sc18 <- rbindlist(lapply(c("liver", "kidney"), function(t)
  fread(paste0("results/pathway_selection/18_axis_gene_stats_", t, ".csv"))[, tissue := t]), fill = TRUE)
key_ct <- c(liver = "Cholangiocyte", kidney = "Endothelial")   # key cell type per tissue; the screen table carries its own

# The axis-gene tables (18, 23) only cover the axis itself; the map-wide screen (15a) covers every gene of
# hsa04152 in every eligible cell type, and must be included or map genes look untested when they are not.
scr <- rbindlist(lapply(c("liver", "kidney"), function(t)
  fread(paste0("results/pathway_selection/AMPK_screen_", t, ".csv"))), fill = TRUE)
sc <- rbind(sc23[, .(human, tissue, cell_type = key_ct[tissue], pb_logFC, pb_P, cell_BH, concordant)],
            sc18[, .(human, tissue, cell_type = key_ct[tissue], pb_logFC, pb_P, cell_BH, concordant)],
            scr[!is.na(pb_P), .(human, tissue, cell_type, pb_logFC, pb_P, cell_BH,
                                concordant = pb_concordant %in% TRUE)])
sc <- sc[concordant == TRUE][order(pb_P)]
sc_best <- sc[, .SD[1], by = human][, .(gene = human, sc_tissue = tissue, sc_cell_type = cell_type,
                                        sc_logFC = pb_logFC, sc_P = pb_P, sc_cell_BH = cell_BH)]

tab <- cons[, .(gene, candidate19, module, phenotype_set, on_map = AMPK_map, axis_anchor = AMPK_axis,
                DEG_FDR, human_logFC = logFC, human_P = P, human_FDR = FDR)]
tab <- merge(tab, sc_best, by = "gene", all.x = TRUE)
tab <- merge(tab, ml[, .(gene, ml_rank = frank(consensus), lasso_freq)], by = "gene", all.x = TRUE)
tab <- merge(tab, mr[, .(gene, mr_OR = OR, mr_p = p, mr_FDR = FDR)], by = "gene", all.x = TRUE)
tab <- merge(tab, val[, .(gene, replicated, val_logFC = logFC_val, val_P = P_val, val_AUC = AUC_val)],
             by = "gene", all.x = TRUE)
tab[is.na(replicated), replicated := FALSE]
tab[, sc_confirmed := !is.na(sc_P)]

# genes whose mechanism is post-translational: nodes of the anchored axis that are kinases/deacetylases, plus
# the kinase subunits of the selected map
phospho <- c("PRKAA1", "PRKAA2", "STK11", "CAMKK2", "SIRT1", "TP53")
tab[, protein_readout := gene %in% phospho]

tab[, tier := fifelse(
  candidate19 & on_map & sc_confirmed & (replicated | (!is.na(ml_rank) & ml_rank <= 5)), "Tier 1",
  fifelse(sc_confirmed & human_P < 0.05, "Tier 2",
  fifelse(candidate19 & on_map & !sc_confirmed, "Tier 2b (not testable in this atlas)",
  fifelse(protein_readout & on_map, "Tier 3 (protein-level)", "Not a target"))))]

setorder(tab, tier, -candidate19, human_P)
fwrite(tab, file.path(out_dir, "target_tiers.csv"))

cat("\n=== Tier 1 ===\n")
print(tab[tier == "Tier 1", .(gene, human_logFC = round(human_logFC, 2), human_P = signif(human_P, 2),
                              axis_anchor, sc_cell_type, sc_P = signif(sc_P, 2), ml_rank,
                              replicated, val_P = signif(val_P, 2), mr_p = signif(mr_p, 2))])
cat("\n=== Tier 2 ===\n")
print(tab[tier == "Tier 2", .(gene, candidate19, human_logFC = round(human_logFC, 2),
                              human_FDR = signif(human_FDR, 2), sc_cell_type, sc_P = signif(sc_P, 2),
                              ml_rank, mr_p = signif(mr_p, 2))])
cat("\n=== Tier 2b: anchored candidates the mouse atlas cannot test ===\n")
print(tab[tier == "Tier 2b (not testable in this atlas)",
          .(gene, human_logFC = round(human_logFC, 2), human_P = signif(human_P, 2), axis_anchor,
            ml_rank, mr_p = signif(mr_p, 2), replicated)])

cat("\n=== Tier 3 (protein-level) ===\n")
print(tab[tier == "Tier 3 (protein-level)", .(gene, human_logFC = round(human_logFC, 2),
                                              human_P = signif(human_P, 2), on_map, axis_anchor,
                                              mr_p = signif(mr_p, 2))])
cat("\n=== Not a target (single-axis evidence) ===\n")
print(tab[tier == "Not a target" & (!is.na(ml_rank) | !is.na(mr_FDR)),
          .(gene, candidate19, ml_rank, lasso_freq, mr_OR = round(mr_OR, 3),
            mr_FDR = signif(mr_FDR, 2), sc_confirmed, replicated)][order(ml_rank)][seq_len(min(12, .N))])

long <- melt(tab[, .(gene, tier, candidate19, module, phenotype_set, on_map, axis_anchor, sc_confirmed,
                     replicated, ml_top5 = !is.na(ml_rank) & ml_rank <= 5, mr_sig = !is.na(mr_FDR) & mr_FDR < 0.05)],
             id.vars = c("gene", "tier"), variable.name = "evidence", value.name = "supported")
fwrite(long, file.path(out_dir, "target_evidence_long.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_24.txt"))
