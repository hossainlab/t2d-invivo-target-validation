# 16_pathway_selection.R — reproducible KEGG pathway + signalling-axis selection (paper step 3-5 analogue)
#
# Reference paper (Xu et al., Phytomedicine 2026) selected its pathway by:
#   (i)   GO enrichment of the intersect genes pointing at the phenotype (Fig 2A),
#   (ii)  KEGG enrichment of the same genes; the chosen term matched the phenotype, not the top rank (Fig 2B),
#   (iii) a pathway interaction network of the top-30 KEGG terms (Fig 2C),
#   (iv)  pathview showing that one axis inside the map was activated (Fig 2D; RAC1 -> PAK1, an edge on the map).
# Steps (i)-(iii) were narrative in the paper. This script turns each of them into an explicit, scored criterion,
# and replaces the visual reading of the pathview map (iv) by parsing the KEGG KGML relations: an "axis" is a
# drawn edge between two differentially expressed nodes whose logFC signs agree with the sign of the edge.
#
# All criteria here use the human bulk data ONLY (no mouse scRNA), so the scRNA analysis in script 17 stays an
# independent confirmation step, exactly as in the paper.
#
# Outputs: results/pathway_selection/selection_candidates.csv   (all criteria per pathway)
#          results/pathway_selection/selection_axis_edges.csv   (all DE-DE KGML edges of the screened maps)
#          results/pathway_selection/selection_ranked.csv       (composite ranking of the ORA-significant maps)
#          results/pathway_selection/kgml/                      (cached KGML)
#          results/supplementary_figures/pathway_selection/16_selection_criteria.pdf

suppressPackageStartupMessages({
  library(data.table); library(limma); library(fgsea); library(clusterProfiler)
  library(org.Hs.eg.db); library(xml2); library(ggplot2); library(pathview)
})
# clusterProfiler reaches KEGG through yulab.utils::yread, which calls readLines() on the URL as a
# bare string. That goes through file() -> url() with R's default URL method, which on this
# Windows box hangs on https and surfaces as "cannot read from connection" after the timeout.
# curl and readLines(url(...)) both fetch the same table in ~3 s, so this is the method, not the
# network. Forcing libcurl makes the bare-string form work.
options(url.method = "libcurl", timeout = max(300, getOption("timeout")))

set.seed(20260916)
out_dir  <- "results/pathway_selection"
fig_dir  <- "results/supplementary_figures/pathway_selection"
kgml_dir <- file.path(out_dir, "kgml")
for (d in c(out_dir, fig_dir, kgml_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

deg   <- fread("results/bulk/DEG_T2D_vs_Control_all.csv")
ints  <- fread("results/bulk/intersect_genes.csv")$gene
hyper <- fread("results/bulk/geneset_hyperglycemia.csv")$gene
ora   <- fread("results/bulk/enrich_KEGG.csv")
expr  <- readRDS("results/bulk/expr_merged_raw.rds")
meta  <- fread("results/bulk/meta_bulk.csv")
meta  <- meta[match(colnames(expr), GSM)]
stopifnot(identical(meta$GSM, colnames(expr)))
meta[, condition := factor(condition, levels = c("Control", "T2D"))]

measured <- deg$gene
setkey(deg, gene)

# ---- 1. KEGG pathway gene sets --------------------------------------------------------
kg  <- download_KEGG("hsa")
p2g <- as.data.table(kg$KEGGPATHID2EXTID); setnames(p2g, c("pathway", "entrez"))
p2n <- as.data.table(kg$KEGGPATHID2NAME);  setnames(p2n, c("pathway", "name"))
p2g[, symbol := suppressMessages(mapIds(org.Hs.eg.db, entrez, "SYMBOL", "ENTREZID"))]
p2g <- p2g[!is.na(symbol)]
if (!grepl("^hsa", p2g$pathway[1])) { p2g[, pathway := paste0("hsa", pathway)]; p2n[, pathway := paste0("hsa", pathway)] }

# Keep mechanistic maps: drop global/overview maps (hsa01xxx) and disease maps (hsa05xxx); the paper's own
# choice (regulation of actin cytoskeleton, hsa04810) is a mechanism map.
cand  <- p2g[!grepl("^hsa(01|05)", pathway) & symbol %in% measured]
sizes <- cand[, .(n_measured = uniqueN(symbol)), by = pathway][n_measured >= 15 & n_measured <= 500]
cand  <- cand[pathway %in% sizes$pathway]
sets  <- lapply(split(cand$symbol, cand$pathway), unique)
message("Candidate mechanistic KEGG maps with 15-500 measured genes: ", length(sets))

# ---- 2. Criterion A: ORA on the intersect genes (paper Fig 2B) ------------------------
ora_dt <- ora[, .(pathway = ID, ora_p = pvalue, ora_padj = p.adjust, ora_count = Count)]

# ---- 3. Criterion B: phenotype match (paper judged this by eye) -----------------------
# Objective stand-in: enrichment of the curated hyperglycemia / insulin-resistance gene set (script 04),
# which is this study's phenotype definition, inside the pathway.
n_hyper_meas <- length(intersect(hyper, measured))
pheno <- rbindlist(lapply(names(sets), function(p) data.table(
  pathway     = p,
  n_measured  = length(sets[[p]]),
  n_pheno     = length(intersect(sets[[p]], hyper)),
  n_intersect = length(intersect(sets[[p]], ints)),
  pheno_frac  = length(intersect(sets[[p]], hyper)) / length(sets[[p]]))))
pheno[, pheno_p := phyper(n_pheno - 1L, n_hyper_meas, length(measured) - n_hyper_meas, n_measured, lower.tail = FALSE)]

# ---- 4. Criterion C: pathway-level differential activity in the full bulk data --------
ranks <- sort(setNames(deg$t, deg$gene), decreasing = TRUE)
gs <- as.data.table(fgsea(sets, ranks, minSize = 15, maxSize = 500, eps = 0))[
  , .(pathway, gsea_NES = NES, gsea_p = pval, gsea_padj = padj)]

design <- model.matrix(~ condition + dataset, data = meta)   # same design as the primary DEG fit (script 02)
idx <- lapply(sets, function(g) which(rownames(expr) %in% g))
idx <- idx[vapply(idx, length, 1L) >= 15]
fr  <- as.data.table(fry(expr, index = idx, design = design, contrast = "conditionT2D"), keep.rownames = "pathway")[
  , .(pathway, fry_dir = Direction, fry_p = PValue, fry_padj = FDR, fry_mixed_p = PValue.Mixed)]

# ---- 5. Criterion D: DE content and directional coherence ----------------------------
de_stats <- rbindlist(lapply(names(sets), function(p) {
  d <- deg[J(sets[[p]]), nomatch = 0L]
  s <- d[P.Value < 0.05]
  data.table(pathway = p, n_DE = nrow(s), n_DE_FDR = d[adj.P.Val < 0.05, .N],
             de_frac = nrow(s) / nrow(d),
             coherence = if (nrow(s)) max(mean(s$logFC > 0), mean(s$logFC < 0)) else NA_real_,
             mean_absLFC = mean(abs(d$logFC)))
}))

# ---- 6. Criterion E: centrality in the KEGG term network (paper Fig 2C) --------------
# Jaccard similarity between the measured gene sets of the top-30 ORA maps; degree counted at J >= 0.1.
top_terms <- ora_dt[order(ora_p)][pathway %in% names(sets)][seq_len(min(30, .N)), pathway]
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
deg_cent <- rbindlist(lapply(top_terms, function(p) data.table(
  pathway = p,
  degree  = sum(vapply(setdiff(top_terms, p), function(q) jac(sets[[p]], sets[[q]]) >= 0.1, TRUE)))))

# ---- 7. Criterion F: coherent DE-DE edges drawn on the KEGG map (paper Fig 2D) -------
act_sub <- c("activation", "expression", "phosphorylation", "indirect effect", "binding/association")
inh_sub <- c("inhibition", "repression", "dephosphorylation")

edges_for <- function(pid) {
  f <- file.path(kgml_dir, paste0(pid, ".xml"))
  if (!file.exists(f)) {
    ok <- tryCatch({ download.kegg(pathway.id = sub("^hsa", "", pid), species = "hsa",
                                   kegg.dir = kgml_dir, file.type = "xml"); TRUE },
                   error = function(e) FALSE)
    if (!ok || !file.exists(f)) return(NULL)
  }
  x <- tryCatch(read_xml(f), error = function(e) NULL); if (is.null(x)) return(NULL)
  en <- xml_find_all(x, "//entry[@type='gene']")
  rel <- xml_find_all(x, "//relation")
  if (!length(en) || !length(rel)) return(NULL)
  emap <- setNames(lapply(xml_attr(en, "name"), function(s) sub("^hsa:", "", strsplit(s, " ")[[1]])), xml_attr(en, "id"))
  rbindlist(lapply(rel, function(r) {
    e1 <- xml_attr(r, "entry1"); e2 <- xml_attr(r, "entry2")
    if (is.null(emap[[e1]]) || is.null(emap[[e2]])) return(NULL)
    sub <- xml_attr(xml_find_all(r, "./subtype"), "name")
    sgn <- if (any(sub %in% inh_sub)) -1L else if (any(sub %in% act_sub)) 1L else NA_integer_
    if (is.na(sgn)) return(NULL)
    g1 <- unique(na.omit(suppressMessages(mapIds(org.Hs.eg.db, emap[[e1]], "SYMBOL", "ENTREZID"))))
    g2 <- unique(na.omit(suppressMessages(mapIds(org.Hs.eg.db, emap[[e2]], "SYMBOL", "ENTREZID"))))
    d1 <- deg[J(g1), nomatch = 0L][P.Value < 0.05]
    d2 <- deg[J(g2), nomatch = 0L][P.Value < 0.05]
    if (!nrow(d1) || !nrow(d2)) return(NULL)
    CJ(i = seq_len(nrow(d1)), j = seq_len(nrow(d2)))[, .(
      pathway = pid, rel_type = xml_attr(r, "type"), subtype = paste(sub, collapse = ";"), edge_sign = sgn,
      gene1 = d1$gene[i], lfc1 = d1$logFC[i], p1 = d1$P.Value[i], fdr1 = d1$adj.P.Val[i],
      gene2 = d2$gene[j], lfc2 = d2$logFC[j], p2 = d2$P.Value[j], fdr2 = d2$adj.P.Val[j])]
  }))
}

screen_ids <- unique(c(top_terms, ora_dt[ora_padj < 0.05, pathway],
                       gs[gsea_padj < 0.05, pathway], fr[fry_padj < 0.05, pathway]))
screen_ids <- intersect(screen_ids, names(sets))
message("Parsing KGML for ", length(screen_ids), " maps")
edges <- rbindlist(lapply(screen_ids, function(p) { e <- edges_for(p); if (!is.null(e) && nrow(e)) e else NULL }))
if (nrow(edges)) {
  edges <- edges[gene1 != gene2]
  edges[, coherent     := sign(lfc1) * sign(lfc2) == edge_sign]
  edges[, both_FDR     := fdr1 < 0.05 & fdr2 < 0.05]
  edges[, in_intersect := (gene1 %in% ints) + (gene2 %in% ints)]
  edges[, in_pheno     := (gene1 %in% hyper) + (gene2 %in% hyper)]
  edges[, worst_p      := pmax(p1, p2)]
  edges <- unique(edges, by = c("pathway", "gene1", "gene2", "edge_sign"))
  fwrite(edges[order(pathway, worst_p)], file.path(out_dir, "selection_axis_edges.csv"))
}
edge_stats <- if (nrow(edges)) edges[, .(n_edges = .N, n_coherent = sum(coherent),
                                         n_coherent_FDR = sum(coherent & both_FDR),
                                         n_coherent_pheno = sum(coherent & in_pheno > 0)), by = pathway] else
  data.table(pathway = character(), n_edges = integer(), n_coherent = integer(),
             n_coherent_FDR = integer(), n_coherent_pheno = integer())

# ---- 7b. Candidate-anchored axis (decision R23) ----------------------------------------
# The unrestricted edge rule can select an axis that shares no gene with the Fig 1 candidate set, which
# breaks the chain (script 22). The reported axis is therefore restricted to edges with at least one node
# in the three-way candidate set, so the axis stays inside the funnel that selected the pathway.
# Discordant edges are kept here: a drawn activation edge whose two DE nodes move apart is the
# transcript-level signature of post-translational control, which is the study's mechanism.
if (nrow(edges)) {
  anch <- edges[in_intersect > 0]
  anch[, anchor_gene := ifelse(gene1 %in% ints, gene1, gene2)]
  fwrite(anch[order(pathway, worst_p)], file.path(out_dir, "selection_axis_anchored.csv"))
  message("Candidate-anchored DE-DE edges: ", nrow(anch), " on ", uniqueN(anch$pathway), " maps")
  print(anch[pathway %in% ora_dt[ora_padj < 0.05, pathway],
             .(pathway, gene1, gene2, subtype, edge_sign, coherent, anchor_gene,
               worst_p = signif(worst_p, 2))][order(pathway, worst_p)])
}

# ---- 8. Assemble and rank -------------------------------------------------------------
tab <- Reduce(function(a, b) merge(a, b, by = "pathway", all.x = TRUE),
              list(pheno, p2n, ora_dt, gs, fr, de_stats, deg_cent, edge_stats))
for (cl in c("n_edges", "n_coherent", "n_coherent_FDR", "n_coherent_pheno", "degree"))
  tab[is.na(get(cl)), (cl) := 0L]
fwrite(tab[order(ora_p)], file.path(out_dir, "selection_candidates.csv"))

short <- tab[!is.na(ora_padj) & ora_padj < 0.05]   # the paper's step-3 pool: maps enriched in the intersect genes
rk <- function(x, decreasing = FALSE) frank(if (decreasing) -x else x, ties.method = "average", na.last = "keep")
short[, score_ora   := rk(ora_p)]
short[, score_pheno := rk(pheno_p)]
short[, score_act   := rk(fry_p)]
short[, score_gsea  := rk(gsea_p)]
short[, score_axis  := rk(n_coherent, decreasing = TRUE)]
short[, score_net   := rk(degree, decreasing = TRUE)]
sc <- c("score_ora", "score_pheno", "score_act", "score_gsea", "score_axis", "score_net")
short[, composite := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = sc]
setorder(short, composite)
fwrite(short, file.path(out_dir, "selection_ranked.csv"))
print(short[, .(pathway, name, n_measured, n_DE, n_DE_FDR, ora_padj, n_pheno, pheno_p, fry_dir, fry_p,
                gsea_NES, gsea_p, degree, n_coherent, n_coherent_FDR, composite)])

# ---- 9. Criteria figure ---------------------------------------------------------------
pl <- melt(short[seq_len(min(10, .N)), .(name = factor(name, levels = rev(name)),
                                         `ORA (intersect genes)`     = -log10(ora_p),
                                         `Phenotype-set enrichment`  = -log10(pheno_p),
                                         `Pathway activity (fry)`    = -log10(fry_p),
                                         `GSEA (full ranking)`       = -log10(gsea_p),
                                         `Coherent DE-DE map edges`  = as.numeric(n_coherent),
                                         `Network degree`            = as.numeric(degree))],
           id.vars = "name", variable.name = "criterion", value.name = "value")
p <- ggplot(pl, aes(value, name)) + geom_col(fill = "#3B7DD8") +
  facet_wrap(~ criterion, scales = "free_x", nrow = 2) +
  labs(x = NULL, y = NULL, title = "Pathway selection criteria (human bulk data only)") + theme_bw(base_size = 9)
ggsave(file.path(fig_dir, "16_selection_criteria.pdf"), p, width = 12, height = 6)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_16.txt"))
