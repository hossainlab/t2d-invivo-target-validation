# 25_map_path_provenance.R — where exactly does each candidate sit on the selected KEGG map, and does the
# drawn route predict the direction we observe?
#
# "PPARGC1A and CCND1 are AMPK genes" is only meaningful if the map actually connects them to AMPK and if the
# sign of that connection matches the human data. This script parses the hsa04152 KGML into a signed directed
# graph, finds the shortest signed path from the AMPK heterotrimer to every Fig 1 candidate on the map, and
# checks the composite sign against the observed logFC under one hypothesis: AMPK output is suppressed in
# obese T2D liver.
#
#   net sign of a path = product of the edge signs (activation/expression = +1, inhibition = -1)
#   prediction under AMPK suppression = -1 * net sign   (less AMPK output => target moves opposite to the path)
#
# Outputs: results/pathway_selection/25_map_paths.csv

suppressPackageStartupMessages({ library(data.table); library(xml2); library(org.Hs.eg.db); library(igraph) })
out_dir <- "results/pathway_selection"
kgml <- file.path(out_dir, "kgml", "hsa04152.xml")
deg  <- fread("results/bulk/DEG_T2D_vs_Control_all.csv"); setkey(deg, gene)
ints <- fread("results/bulk/intersect_genes.csv")$gene

act_sub <- c("activation", "expression", "phosphorylation", "indirect effect", "binding/association")
inh_sub <- c("inhibition", "repression", "dephosphorylation")

x  <- read_xml(kgml)
en <- xml_find_all(x, "//entry[@type='gene']")
sym_of <- setNames(lapply(xml_attr(en, "name"), function(nm)
  unique(na.omit(suppressMessages(mapIds(org.Hs.eg.db, sub("^hsa:", "", strsplit(nm, " ")[[1]]),
                                         "SYMBOL", "ENTREZID"))))), xml_attr(en, "id"))
node_label <- vapply(sym_of, function(s) if (length(s) > 3) paste0(s[1], "+", length(s) - 1L) else paste(s, collapse = "/"), "")

rel <- xml_find_all(x, "//relation")
ed <- rbindlist(lapply(rel, function(r) {
  e1 <- xml_attr(r, "entry1"); e2 <- xml_attr(r, "entry2")
  if (is.null(sym_of[[e1]]) || is.null(sym_of[[e2]])) return(NULL)
  sub <- xml_attr(xml_find_all(r, "./subtype"), "name")
  sgn <- if (any(sub %in% inh_sub)) -1L else if (any(sub %in% act_sub)) 1L else NA_integer_
  if (is.na(sgn)) return(NULL)
  data.table(from = e1, to = e2, sign = sgn, subtype = paste(sub, collapse = ";"))
}))
ed <- unique(ed, by = c("from", "to", "sign"))
g <- graph_from_data_frame(ed[, .(from, to)], directed = TRUE)
E(g)$sign <- ed$sign; E(g)$subtype <- ed$subtype

# the AMPK heterotrimer entry: the node whose symbol set contains the catalytic subunits
ampk_id <- names(sym_of)[vapply(sym_of, function(s) all(c("PRKAA1", "PRKAA2") %in% s), TRUE)][1]
stopifnot(!is.na(ampk_id))
message("AMPK complex entry id ", ampk_id, ": ", node_label[[ampk_id]])

targets <- unique(unlist(sym_of))          # every gene on the map, not only the candidates
res <- rbindlist(lapply(targets, function(tg) {
  ids <- names(sym_of)[vapply(sym_of, function(s) tg %in% s, TRUE)]
  ids <- intersect(ids, V(g)$name)
  if (!length(ids)) return(data.table(gene = tg, on_map = tg %in% unlist(sym_of), path = NA_character_))
  best <- NULL
  for (id in ids) {
    sp <- suppressWarnings(shortest_paths(g, from = ampk_id, to = id, output = "both"))
    if (!length(sp$vpath[[1]])) next
    v <- sp$vpath[[1]]; e <- sp$epath[[1]]
    if (is.null(best) || length(v) < length(best$v)) best <- list(v = v, e = e)
  }
  if (is.null(best)) return(data.table(gene = tg, on_map = TRUE, path = "no directed path from AMPK"))
  net <- prod(E(g)$sign[best$e])
  lab <- paste(vapply(names(best$v), function(n) node_label[[n]], ""),
               collapse = " -> ")
  steps <- paste(ifelse(E(g)$sign[best$e] > 0, "+", "-"), collapse = "")
  d <- deg[J(tg), nomatch = 0L]
  data.table(gene = tg, on_map = TRUE, n_steps = length(best$e), path = lab, edge_signs = steps,
             subtypes = paste(E(g)$subtype[best$e], collapse = " | "), net_sign = net,
             predicted_if_AMPK_suppressed = ifelse(-net > 0, "up", "down"),
             observed_logFC = if (nrow(d)) round(d$logFC, 2) else NA_real_,
             observed_P = if (nrow(d)) signif(d$P.Value, 2) else NA_real_)
}), fill = TRUE)
res[, direction_matches := !is.na(observed_logFC) &
      ((predicted_if_AMPK_suppressed == "up" & observed_logFC > 0) |
       (predicted_if_AMPK_suppressed == "down" & observed_logFC < 0))]
res[, fig1_candidate := gene %in% ints]

# ---- attach the downstream evidence, so the attrition can be counted ---------------------
pw   <- fread(file.path(out_dir, "AMPK_pathway_genes.csv"))
res  <- merge(res, pw[, .(gene = human, measured = !is.na(human_P), human_DE, human_FDR, mouse)],
              by = "gene", all.x = TRUE)
sc23 <- rbindlist(lapply(c("liver", "kidney"), function(t)
  fread(paste0(out_dir, "/23_anchored_gene_stats_", t, ".csv"))[, tissue := t]), fill = TRUE)
scr  <- rbindlist(lapply(c("liver", "kidney"), function(t)
  fread(paste0(out_dir, "/AMPK_screen_", t, ".csv"))), fill = TRUE)
sc_ok <- unique(c(sc23[concordant == TRUE, human], scr[pb_concordant == TRUE, human]))
res[, sc_confirmed := gene %in% sc_ok]
ml  <- fread("results/ml/feature_ranking.csv")
mr  <- fread("results/mr/mr_results_blood_ldclumped.csv")[method %in% c("IVW", "Wald ratio")]
val <- fread("results/validation/GSE23343_candidates.csv")
res <- merge(res, ml[, .(gene, ml_rank = frank(consensus))], by = "gene", all.x = TRUE)
res <- merge(res, mr[, .(gene, mr_FDR = FDR)], by = "gene", all.x = TRUE)
res <- merge(res, val[, .(gene, replicated)], by = "gene", all.x = TRUE)
res[is.na(replicated), replicated := FALSE]

cat("
=== Attrition over the AMPK map ===
")
att <- data.table(
  step = c("1. genes on KEGG hsa04152",
           "2. measured in the human bulk data",
           "3. differentially expressed (P < 0.05, |logFC| > 0.5)",
           "4. + reachable from the AMPK node by a drawn directed path",
           "5. + drawn sign predicts the observed direction",
           "6. + confirmed in the mouse atlas (pseudobulk, human direction)",
           "7. + inside the Fig 1 candidate set (Tier 1 rule)"),
  n = c(nrow(res),
        res[measured == TRUE, .N],
        res[human_DE == TRUE, .N],
        res[human_DE == TRUE & !is.na(n_steps), .N],
        res[human_DE == TRUE & !is.na(n_steps) & direction_matches == TRUE, .N],
        res[human_DE == TRUE & !is.na(n_steps) & direction_matches == TRUE & sc_confirmed == TRUE, .N],
        res[human_DE == TRUE & !is.na(n_steps) & direction_matches == TRUE & sc_confirmed == TRUE &
            fig1_candidate == TRUE, .N]))
print(att)
cat("
Genes surviving step 5 (DE + reachable + sign matches):
")
print(res[human_DE == TRUE & !is.na(n_steps) & direction_matches == TRUE,
          .(gene, n_steps, path, observed_logFC, human_FDR = signif(human_FDR, 2), sc_confirmed,
            fig1_candidate, ml_rank, replicated)][order(n_steps, -sc_confirmed)])
cat("
FDR-significant map genes that are NOT Fig 1 candidates (would the gate exclude a good target?):
")
print(res[human_FDR < 0.05 & fig1_candidate == FALSE,
          .(gene, observed_logFC, human_FDR = signif(human_FDR, 2), n_steps, direction_matches,
            sc_confirmed, ml_rank, replicated)])
setorder(res, -fig1_candidate, n_steps, gene)
fwrite(res, file.path(out_dir, "25_map_paths.csv"))
print(res[!is.na(n_steps), .(gene, fig1_candidate, n_steps, path, edge_signs, net_sign,
                             predicted_if_AMPK_suppressed, observed_logFC, observed_P, direction_matches)])
cat("\nCandidates on the map with a directed path from AMPK: ", res[!is.na(n_steps) & fig1_candidate, .N], "\n")
cat("of which the drawn sign predicts the observed direction: ",
    res[!is.na(n_steps) & fig1_candidate & direction_matches, .N], "\n")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_25.txt"))
