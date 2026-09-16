# 17_axis_confirmation_scRNA.R — systematic scRNA confirmation of the bulk-shortlisted pathways and their axes
# (paper Fig 3D-G analogue, made exhaustive and multiplicity-controlled).
#
# The reference paper tested ONE pair (RAC1 -> PAK1) in ONE cell type and reported that it changed, co-localized
# and correlated. Its bulk and single-cell data came from the same disease and tissue, so one pathway could be
# carried straight through. Here the human bulk (obese T2D liver) and the mouse atlas (STZ liver + kidney) differ
# in species, model and tissue, so a two-stage design is used and declared:
#   stage 1 (script 16, human bulk only): shortlist = KEGG maps enriched in the intersect genes, scored on six
#           paper-derived criteria; nothing from the mouse data enters this stage;
#   stage 2 (this script, mouse only):    every axis candidate of every shortlisted map is put through the
#           paper's three Fig 3 criteria - changed, co-localized, correlated - in every eligible cell type,
#           with BH correction across all tested pairs within a tissue.
# The full matrix is reported, so the final choice is transparent rather than post hoc.
#
# Fixed definitions (set before looking at the mouse data):
#   axis component   = connected component of the coherent DE-DE edge graph of a map (script 16);
#                      a map's primary axis is the component with the most nodes (ties: lowest median worst_p).
#   discordant edge  = drawn edge whose two DE nodes move against the edge sign; marks nodes whose regulation is
#                      post-translational and therefore invisible to transcript data.
#   testable pair    = both 1:1 mouse orthologs detected in >= 5% / >= 2% of that cell type's cells.
#   changed          = per-mouse pseudobulk (script 09) P < 0.05 with the same sign as the human logFC.
#   co-localized     = >= 2% of the cell type's cells co-express the pair.
#   correlated       = cell-level Pearson BH < 0.05 with R > 0 (the paper's Fig 3G test).
#
# Outputs: results/pathway_selection/17_axis_components.csv
#          results/pathway_selection/17_gene_confirmation.csv
#          results/pathway_selection/17_edge_confirmation_<tissue>.csv
#          results/pathway_selection/17_pathway_confirmation.csv
#          results/supplementary_figures/pathway_selection/17_edge_confirmation.pdf

suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(babelgene); library(igraph)
  library(ggplot2); library(Matrix)
})
set.seed(20260916)
out_dir <- "results/pathway_selection"
fig_dir <- "results/supplementary_figures/pathway_selection"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

short <- fread(file.path(out_dir, "selection_ranked.csv"))
deg   <- fread("results/bulk/DEG_T2D_vs_Control_all.csv")
edges <- fread(file.path(out_dir, "selection_axis_edges.csv"))[pathway %in% short$pathway]
edges <- merge(edges, short[, .(pathway, pw_name = name, composite)], by = "pathway")
message("Shortlisted maps: ", nrow(short), "; DE-DE edges on them: ", nrow(edges))

# ---- 1. Axis components per map (human bulk only) -------------------------------------
comp_list <- lapply(split(edges[coherent == TRUE], by = "pathway"), function(e) {
  if (!nrow(e)) return(NULL)
  g <- graph_from_data_frame(e[, .(gene1, gene2)], directed = FALSE)
  memb <- components(g)$membership
  e[, component := memb[gene1]]
  e[, .(n_nodes = uniqueN(c(gene1, gene2)), n_edges = .N, median_worst_p = median(worst_p),
        n_pheno_nodes = uniqueN(c(gene1[in_pheno > 0], gene2[in_pheno > 0])),
        nodes = paste(sort(unique(c(gene1, gene2))), collapse = ";")),
    by = .(pathway, pw_name, component)]
})
comps <- rbindlist(comp_list)[order(pathway, -n_nodes, median_worst_p)]
comps[, rank_in_map := seq_len(.N), by = pathway]
fwrite(comps, file.path(out_dir, "17_axis_components.csv"))
print(comps[rank_in_map == 1, .(pathway, pw_name, n_nodes, n_edges, n_pheno_nodes, nodes)])

# ---- 2. Human -> mouse orthologs for every gene on the shortlisted maps ---------------
kegg_genes <- unique(c(edges$gene1, edges$gene2))
pw_all <- fread(file.path(out_dir, "selection_candidates.csv"))
o <- as.data.table(orthologs(genes = kegg_genes, species = "mouse", human = TRUE))
o <- unique(o[, .(human = human_symbol, mouse = symbol)])[, n_h := .N, by = human][, n_m := .N, by = mouse][n_h == 1 & n_m == 1]
m_of <- setNames(o$mouse, o$human)
edges[, `:=`(m1 = m_of[gene1], m2 = m_of[gene2])]
edges_t <- edges[!is.na(m1) & !is.na(m2)]
message("Edges with 1:1 orthologs for both nodes: ", nrow(edges_t), " of ", nrow(edges))

hum_dir <- setNames(deg$logFC, deg$gene)

# ---- 3. Confirmation per tissue --------------------------------------------------------
gene_all <- list(); edge_all <- list()
for (tt in c("liver", "kidney")) {
  so <- readRDS(file.path("results/sc", paste0(tt, "_annotated.rds")))
  pb <- fread(file.path("results/sc", paste0(tt, "_pseudobulk_DE_all.csv")))
  md <- as.data.table(so@meta.data)
  eligible <- md[, .(n_ctl = sum(condition == "Control"), n_stz = sum(condition == "STZ")), by = cell_type][
    n_ctl >= 20 & n_stz >= 20][order(-(n_ctl + n_stz))]

  # ---- 3a. Gene level: pseudobulk concordance with the human direction ----------------
  gmap <- data.table(human = names(m_of), mouse = unname(m_of))
  gmap[, human_logFC := hum_dir[human]]
  gmap <- gmap[!is.na(human_logFC) & human %in% deg[P.Value < 0.05, gene]]
  gd <- merge(pb[, .(mouse = gene, cell_type, pb_logFC = logFC, pb_P = PValue, pb_FDR = FDR)], gmap, by = "mouse")
  gd[, concordant := pb_P < 0.05 & sign(pb_logFC) == sign(human_logFC)]
  gd[, tissue := tt]
  gene_all[[tt]] <- gd

  # ---- 3b. Edge level: co-localization and correlation --------------------------------
  genes_need <- intersect(unique(c(edges_t$m1, edges_t$m2)), rownames(so))
  ex <- as.data.table(FetchData(so, vars = genes_need, layer = "data"))
  ex[, `:=`(cell_type = md$cell_type, condition = md$condition, mouse_id = md$mouse_id)]
  ep <- unique(edges_t[m1 %in% genes_need & m2 %in% genes_need, .(m1, m2)])
  rows <- list()
  for (ct in eligible$cell_type) {
    e <- ex[cell_type == ct]
    det <- vapply(genes_need, function(gn) mean(e[[gn]] > 0), 0)
    for (k in seq_len(nrow(ep))) {
      a <- ep$m1[k]; b <- ep$m2[k]
      if (max(det[a], det[b]) < 0.05 || min(det[a], det[b]) < 0.02) next
      x <- e[[a]]; y <- e[[b]]
      if (sd(x) == 0 || sd(y) == 0) next
      ctt <- suppressWarnings(cor.test(x, y, method = "pearson"))
      co  <- which(x > 0 & y > 0)
      sp  <- if (length(co) >= 10) suppressWarnings(cor.test(x[co], y[co], method = "spearman")) else NULL
      pm  <- e[, .(a = mean(get(a)), b = mean(get(b)), n = .N), by = .(mouse_id, condition)][n >= 10]
      pmt <- if (nrow(pm) >= 4) suppressWarnings(cor.test(pm$a, pm$b, method = "pearson")) else NULL
      rows[[length(rows) + 1L]] <- data.table(
        tissue = tt, cell_type = ct, mouse1 = a, mouse2 = b, n_cells = nrow(e),
        det1 = det[a], det2 = det[b], pct_coexpr = length(co) / nrow(e),
        cell_R = unname(ctt$estimate), cell_p = ctt$p.value,
        coexpr_rho = if (is.null(sp)) NA_real_ else unname(sp$estimate),
        coexpr_p   = if (is.null(sp)) NA_real_ else sp$p.value,
        pb_r = if (is.null(pmt)) NA_real_ else unname(pmt$estimate),
        pb_p = if (is.null(pmt)) NA_real_ else pmt$p.value, n_mice = nrow(pm))
    }
  }
  r <- rbindlist(rows)
  if (nrow(r)) {
    r[, cell_BH := p.adjust(cell_p, "BH")]
    r[, pb_BH   := p.adjust(pb_p, "BH")]
    r <- merge(r, unique(edges_t[, .(mouse1 = m1, mouse2 = m2, pathway, pw_name,
                                     gene1, gene2, kind = ifelse(coherent, "coherent", "discordant"))]),
               by = c("mouse1", "mouse2"), allow.cartesian = TRUE)
    fwrite(r[order(cell_p)], file.path(out_dir, paste0("17_edge_confirmation_", tt, ".csv")))
    edge_all[[tt]] <- r
  }
  rm(so, ex); gc(verbose = FALSE)
}

gene_dt <- rbindlist(gene_all)
fwrite(gene_dt, file.path(out_dir, "17_gene_confirmation.csv"))
edge_dt <- rbindlist(edge_all)

# ---- 4. Per-map confirmation summary ---------------------------------------------------
# genes of each map that are human-DE, with the best cell type per map/tissue
map_genes <- unique(rbind(edges_t[, .(pathway, pw_name, human = gene1)], edges_t[, .(pathway, pw_name, human = gene2)]))
gsumm <- merge(gene_dt, map_genes, by = "human", allow.cartesian = TRUE)[
  , .(genes_tested = uniqueN(mouse), genes_concordant = sum(concordant, na.rm = TRUE)),
  by = .(pathway, pw_name, tissue, cell_type)][order(pathway, -genes_concordant)]
gbest <- gsumm[, .SD[1], by = .(pathway, pw_name)]

esumm <- if (nrow(edge_dt)) edge_dt[kind == "coherent", .(
  edges_tested   = uniqueN(paste(mouse1, mouse2)),
  edges_coloc    = uniqueN(paste(mouse1, mouse2)[pct_coexpr >= 0.02]),
  edges_corr     = uniqueN(paste(mouse1, mouse2)[cell_BH < 0.05 & cell_R > 0 & pct_coexpr >= 0.02]),
  best_edge      = paste(mouse1[which.min(cell_BH)], mouse2[which.min(cell_BH)], sep = "-"),
  best_edge_type = cell_type[which.min(cell_BH)],
  best_edge_R    = round(cell_R[which.min(cell_BH)], 2),
  best_edge_BH   = signif(cell_BH[which.min(cell_BH)], 2)), by = .(pathway, pw_name)] else data.table()

final <- Reduce(function(a, b) merge(a, b, by = c("pathway", "pw_name"), all = TRUE),
                list(short[, .(pathway, pw_name = name, composite, ora_padj, pheno_p, fry_p, gsea_p, n_coherent)],
                     gbest[, .(pathway, pw_name, best_tissue = tissue, best_cell_type = cell_type,
                               genes_tested, genes_concordant)],
                     esumm))
setorder(final, composite)
fwrite(final, file.path(out_dir, "17_pathway_confirmation.csv"))
print(final)

if (nrow(edge_dt)) {
  cat("\n--- Pairs meeting all three paper criteria in a cell type (changed nodes + coloc + correlated) ---\n")
  chg <- gene_dt[concordant == TRUE, paste(tissue, cell_type, mouse)]
  ok <- edge_dt[kind == "coherent" & pct_coexpr >= 0.02 & cell_BH < 0.05 & cell_R > 0]
  ok[, both_changed := paste(tissue, cell_type, mouse1) %in% chg & paste(tissue, cell_type, mouse2) %in% chg]
  ok[, one_changed  := paste(tissue, cell_type, mouse1) %in% chg | paste(tissue, cell_type, mouse2) %in% chg]
  print(ok[order(-both_changed, -one_changed, cell_BH),
           .(pw_name, tissue, cell_type, mouse1, mouse2, pct_coexpr = round(pct_coexpr, 3),
             cell_R = round(cell_R, 2), cell_BH = signif(cell_BH, 2), pb_r = round(pb_r, 2),
             both_changed, one_changed)][seq_len(min(30, .N))])

  p <- ggplot(edge_dt[!is.na(cell_R) & kind == "coherent"],
              aes(paste(mouse1, mouse2, sep = "-"), cell_type, fill = cell_R, size = pct_coexpr)) +
    geom_point(shape = 21, colour = "grey30") +
    facet_grid(tissue ~ pw_name, scales = "free", space = "free") +
    scale_fill_gradient2(low = "#3B7DD8", mid = "white", high = "#C8322F", midpoint = 0) +
    labs(x = NULL, y = NULL, fill = "cell-level R", size = "% co-expressing",
         title = "scRNA confirmation of the shortlisted KEGG map axes") +
    theme_bw(base_size = 7) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
  ggsave(file.path(fig_dir, "17_edge_confirmation.pdf"), p, width = 18, height = 9, limitsize = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_17.txt"))
