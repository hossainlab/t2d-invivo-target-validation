# 22_consistency_audit.R — does one gene set actually flow through the whole study?
#
# The manuscript chain is: Fig 1 (bulk DEG + WGCNA + phenotype set) -> Fig 2 (pathway selection) ->
# pathway axis -> Fig 3 (single-cell confirmation) -> Fig 4 (machine learning) -> Fig 5 (Mendelian
# randomization). This script builds one gene x stage matrix so that every break in that chain is explicit
# and countable rather than argued in prose.
#
# Stages (columns):
#   S1  DEG            nominal P < 0.05 and |logFC| > 0.5 in the human bulk contrast
#   S1b DEG_FDR        adj.P < 0.05
#   S2  module         member of a WGCNA module with T2D p < 0.1
#   S3  phenotype_set  member of the curated hyperglycemia / insulin-resistance set
#   S4  candidate19    the Fig 1 three-way intersection actually carried forward
#   S5  AMPK_map       gene on KEGG hsa04152, the selected map
#   S6  AMPK_axis      node of the primary coherent DE-DE component of that map
#   S7  p53_axis       the axis that survived single-cell confirmation (CDKN1A, CCND1 + PHLDA3, ZMAT3, BAX)
#   S8  sc_confirmed   pseudobulk-concordant with the human direction in kidney endothelium or cholangiocytes
#   S9  ML_top7        top 7 of the consensus machine-learning ranking
#   S10 MR_FDR05       FDR < 0.05 in the LD-clumped cis-MR
#   S11 replicated     replicates in the independent bulk cohort GSE23343
#
# Outputs: results/consistency/gene_stage_matrix.csv, stage_overlap_summary.csv
# Figure:  results/supplementary_figures/consistency/22_chain_consistency.pdf

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
out_dir <- "results/consistency"
fig_dir <- "results/supplementary_figures/consistency"
for (d in c(out_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

deg   <- fread("results/bulk/DEG_T2D_vs_Control_all.csv")
mm    <- fread("results/bulk/WGCNA_gene_module_GS_MM.csv")
mt    <- fread("results/bulk/WGCNA_module_trait.csv")
hyper <- fread("results/bulk/geneset_hyperglycemia.csv")$gene
ints  <- fread("results/bulk/intersect_genes.csv")$gene
pwg   <- fread("results/pathway_selection/AMPK_pathway_genes.csv")      # hsa04152 membership
comps <- fread("results/pathway_selection/17_axis_components.csv")
ml    <- fread("results/ml/feature_ranking.csv")
mr    <- fread("results/mr/mr_results_blood_ldclumped.csv")[method %in% c("IVW", "Wald ratio")]
val   <- fread("results/validation/GSE23343_candidates.csv")
scg   <- fread("results/pathway_selection/17_gene_confirmation.csv")
sc18  <- rbindlist(lapply(c("kidney", "liver"), function(t)
  fread(paste0("results/pathway_selection/18_axis_gene_stats_", t, ".csv"))[, tissue := t]), fill = TRUE)

key_modules <- mt[p_T2D < 0.1, module]
ampk_axis_unanchored <- strsplit(comps[pathway == "hsa04152"][order(-n_nodes, median_worst_p)][1, nodes], ";")[[1]]
anch <- fread("results/pathway_selection/selection_axis_anchored.csv")[pathway == "hsa04152"]
ampk_axis <- unique(c(anch$gene1, anch$gene2, "CCND1"))   # reported, candidate-anchored axis (R23)
p53_axis    <- c("CDKN1A", "CCND1", "PHLDA3", "ZMAT3", "BAX")
ml_top      <- head(ml[order(consensus), gene], 7)

# genes confirmed in the mouse single-cell data (human symbols)
sc_ok <- unique(c(toupper(sc18[concordant == TRUE, gene]),
                  unique(scg[concordant == TRUE, human])))

universe <- sort(unique(c(ints, ampk_axis, ampk_axis_unanchored, p53_axis, ml_top,
                          mr[FDR < 0.05, gene], pwg[human_DE == TRUE, human])))

tab <- data.table(gene = universe)
tab[, DEG           := gene %in% deg[P.Value < 0.05 & abs(logFC) > 0.5, gene]]
tab[, DEG_FDR       := gene %in% deg[adj.P.Val < 0.05, gene]]
tab[, module        := gene %in% mm[module %in% key_modules, gene]]
tab[, phenotype_set := gene %in% hyper]
tab[, candidate19   := gene %in% ints]
tab[, AMPK_map      := gene %in% pwg$human]
tab[, AMPK_axis     := gene %in% ampk_axis]
tab[, axis_unanchored := gene %in% ampk_axis_unanchored]
tab[, p53_axis      := gene %in% p53_axis]
tab[, sc_confirmed  := gene %in% sc_ok]
tab[, ML_top7       := gene %in% ml_top]
tab[, MR_FDR05      := gene %in% mr[FDR < 0.05, gene]]
tab[, replicated    := gene %in% val[replicated == TRUE, gene]]
tab <- merge(tab, deg[, .(gene, logFC, P = P.Value, FDR = adj.P.Val)], by = "gene", all.x = TRUE)
tab <- merge(tab, ml[, .(gene, ml_rank = frank(consensus))], by = "gene", all.x = TRUE)
tab <- merge(tab, mr[, .(gene, mr_OR = OR, mr_p = p, mr_FDR = FDR)], by = "gene", all.x = TRUE)
stage_cols <- c("DEG", "DEG_FDR", "module", "phenotype_set", "candidate19", "AMPK_map", "AMPK_axis",
                "p53_axis", "sc_confirmed", "ML_top7", "MR_FDR05", "replicated")
tab[, n_stages := rowSums(as.matrix(.SD)), .SDcols = stage_cols]
setorder(tab, -n_stages, gene)
fwrite(tab, file.path(out_dir, "gene_stage_matrix.csv"))

cat("\n=== Genes by number of stages passed ===\n")
print(tab[, .(gene, n_stages, candidate19, AMPK_map, AMPK_axis, p53_axis, sc_confirmed, ML_top7, MR_FDR05,
              replicated)][seq_len(min(25, .N))])

# ---- pairwise overlap between the sets that define each figure ---------------------------
sets <- list(`Fig1 candidate 19` = ints,
             `Fig2 axis (unanchored, superseded)` = ampk_axis_unanchored,
             `Fig2 selected map (hsa04152)` = pwg[human_DE == TRUE, human],
             `Fig2 axis (anchored, reported)` = ampk_axis,
             `Fig3 confirmed axis` = p53_axis,
             `Fig4 ML top 7` = ml_top,
             `Fig5 MR FDR<0.05` = mr[FDR < 0.05, gene])
ov <- rbindlist(lapply(names(sets), function(a) rbindlist(lapply(names(sets), function(b)
  data.table(set_a = a, set_b = b, n_a = length(sets[[a]]), n_b = length(sets[[b]]),
             overlap = length(intersect(sets[[a]], sets[[b]])),
             genes = paste(intersect(sets[[a]], sets[[b]]), collapse = ";"))))))
fwrite(ov, file.path(out_dir, "stage_overlap_summary.csv"))
cat("\n=== Overlap between the gene sets that define each figure ===\n")
print(dcast(ov, set_a ~ set_b, value.var = "overlap"))

cat("\n=== Genes carried by every stage of the chain ===\n")
print(tab[candidate19 == TRUE & AMPK_map == TRUE & sc_confirmed == TRUE & ML_top7 == TRUE & MR_FDR05 == TRUE, gene])
cat("(empty means no single gene survives all five figures)\n")
cat("\nCandidate19 AND on the selected map: ", paste(intersect(ints, pwg$human), collapse = ", "), "\n")
cat("Candidate19 AND single-cell confirmed: ", paste(intersect(ints, sc_ok), collapse = ", "), "\n")
cat("Candidate19 AND ML top 7: ", paste(intersect(ints, ml_top), collapse = ", "), "\n")
cat("Candidate19 AND MR FDR<0.05: ", paste(intersect(ints, mr[FDR < 0.05, gene]), collapse = ", "), "\n")

# ---- figure -------------------------------------------------------------------------------
long <- melt(tab[n_stages >= 2, c("gene", stage_cols, "n_stages"), with = FALSE],
             id.vars = c("gene", "n_stages"), variable.name = "stage", value.name = "member")
long[, gene := factor(gene, levels = tab[n_stages >= 2][order(n_stages), gene])]
long[, stage := factor(stage, levels = stage_cols)]
p <- ggplot(long, aes(stage, gene, fill = member)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(`TRUE` = "#C8322F", `FALSE` = "grey92"), guide = "none") +
  labs(x = NULL, y = NULL,
       title = "Consistency of the analysis chain",
       subtitle = "Red = the gene is carried by that stage. A consistent chain would show continuous rows.") +
  theme_bw(base_size = 8) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "22_chain_consistency.pdf"), p, width = 7.5,
       height = 0.16 * long[, uniqueN(gene)] + 2, limitsize = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_22.txt"))
