# 23_anchored_axis_figures.R — Fig 3D-G for the candidate-anchored axis (decision R23).
#
# Script 22 showed the chain broke at the axis step: the unrestricted KGML rule picked an axis
# (CAMKK2/PRKAA2/GYS2/PFKFB4/MLYCD) that shares no gene with the Fig 1 candidate set. Script 16 now also
# reports candidate-anchored edges, and on the selected map (hsa04152) those are:
#     PRKAA2 -> PPARGC1A   activation;phosphorylation, DISCORDANT (subunit up, target down)
#     IRS2   -> PIK3CA     activation, coherent (both down)
#     IRS2   -> PIK3R2 / PIK3CD, discordant
# with CCND1 as the de-repressed output node of the same map, and SIRT1 -> PPARGC1A coherent on hsa04211.
# The reported axis is therefore AMPK(alpha2) -> PGC-1alpha with CCND1 as output; both anchors (PPARGC1A,
# CCND1) are Fig 1 candidates, both are on the selected map, and both are confirmed in the mouse atlas.
#
# Same three paper criteria as script 18, same background null. The p53 arrest pair (Cdkn1a-Ccnd1) remains
# the downstream readout and keeps its own panels (Fig 3H-I) from script 18.
#
# Outputs: results/pathway_selection/23_anchored_gene_stats_<tissue>.csv
#          results/pathway_selection/23_anchored_score_<tissue>.csv
#          results/pathway_selection/23_anchored_correlation_<tissue>.csv
# Figures: figures/Fig3_scRNA/Fig3{D,E,F,G}_anchored_*_<tissue>.pdf

suppressPackageStartupMessages({
  library(Seurat); library(Nebulosa); library(UCell); library(data.table)
  library(ggplot2); library(patchwork); library(Matrix)
  library(xml2); library(babelgene); library(org.Hs.eg.db); library(AnnotationDbi)
})
set.seed(20260916)
fig_main <- "results/supplementary_figures/pathway_selection"; out_dir <- "results/pathway_selection"
  # draft panels; figures/ is owned by the 27-32 publication figure scripts (Fig 3 is now
  # Fig3_composite_<tissue> from 27 and Fig3_axis_<tissue> from 32). Writing here from an analysis
  # script silently overwrites them, because Windows filenames are case-insensitive.

dir.create(fig_main, recursive = TRUE, showWarnings = FALSE)
cols <- c(Control = "#3B7DD8", STZ = "#C8322F")

SELECTED_MAP <- "hsa04152"     # AMPK signalling; the map script 16 selects (decision R23)
pair         <- c("Ppargc1a", "Ccnd1")   # the two Fig 1 candidates that anchor the axis
n_bg         <- 2000L
fmt_p <- function(p) ifelse(is.na(p), "NA", formatC(p, format = "g", digits = 2))

# ---- gene sets, derived rather than typed in ------------------------------------------------
# The Fig 1 candidate set is read from the Fig 1g intersection every run. It used to be a literal
# vector here, which went stale the moment decision M11b cut the set from 19 genes to 8: the list
# still named FBP1, a turquoise-module gene that is no longer a candidate. Anything that marks
# "this is a Fig 1 candidate" has to come from the file that defines them.
fig1_candidates <- fread("results/bulk/intersect_genes.csv")$gene

# Map membership comes from the KGML of the selected map, so "lies on the selected map" is checked
# against KEGG rather than asserted.
map_xml <- file.path("results/pathway_selection/kgml", paste0(SELECTED_MAP, ".xml"))
if (!file.exists(map_xml)) map_xml <- file.path("results/pathway_selection/pathview",
                                                paste0(SELECTED_MAP, ".xml"))
stopifnot(file.exists(map_xml))
map_ids <- unique(unlist(strsplit(
  xml2::xml_attr(xml2::xml_find_all(xml2::read_xml(map_xml), "//entry[@type='gene']"), "name"), " ")))
map_human <- unique(na.omit(AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db, sub("^hsa:", "", map_ids[grepl("^hsa:", map_ids)]),
  "SYMBOL", "ENTREZID", multiVals = "first")))

candidates <- intersect(fig1_candidates, map_human)   # Fig 1 candidates that lie on the selected map
stopifnot(length(candidates) > 0)

# The axis panel also needs the upstream regulators and the output node. These are NOT candidates and
# cannot be: PRKAA2, STK11 and SIRT1 are not differentially expressed in bulk, so they can never reach
# the Fig 1g intersection. Including them is the point - testing the axis only with the genes that
# selected the map would be circular. The set is curated, but every member is asserted to be on the
# map, and any candidate on the map is added automatically so this cannot go stale again.
axis_human <- union(c("PRKAA1", "PRKAA2", "STK11", "CAMKK2", "SIRT1"), candidates)
missing_on_map <- setdiff(axis_human, map_human)
if (length(missing_on_map))
  stop("not on ", SELECTED_MAP, ": ", paste(missing_on_map, collapse = ", "))

orth <- as.data.table(babelgene::orthologs(genes = axis_human, species = "mouse", human = TRUE))
orth <- unique(orth[, .(human = human_symbol, mouse = symbol)])[
  , n_h := .N, by = human][, n_m := .N, by = mouse][n_h == 1 & n_m == 1]
human_of   <- setNames(orth$human, orth$mouse)
axis_genes <- orth$mouse
message(SELECTED_MAP, ": ", length(candidates), " Fig 1 candidates on the map (",
        paste(candidates, collapse = ", "), "); axis genes tested: ",
        paste(axis_genes, collapse = ", "))

deg <- fread("results/bulk/DEG_T2D_vs_Control_all.csv")
hum <- deg[gene %in% human_of, .(human = gene, human_logFC = logFC, human_P = P.Value, human_FDR = adj.P.Val)]

for (tt in c("liver", "kidney")) {
  so <- readRDS(file.path("results/sc", paste0(tt, "_annotated.rds")))
  pb <- fread(file.path("results/sc", paste0(tt, "_pseudobulk_DE_all.csv")))
  md <- as.data.table(so@meta.data)
  genes <- intersect(axis_genes, rownames(so))

  eligible <- md[, .(n_ctl = sum(condition == "Control"), n_stz = sum(condition == "STZ")), by = cell_type][
    n_ctl >= 20 & n_stz >= 20]
  pb_ax <- pb[gene %in% genes & cell_type %in% eligible$cell_type]
  pb_ax[, human := human_of[gene]]
  pb_ax <- merge(pb_ax, hum, by = "human")
  pb_ax[, concordant := PValue < 0.05 & sign(logFC) == sign(human_logFC)]
  # key cell type: most anchored (Fig 1 candidate) genes concordant, then most axis genes concordant
  rankct <- pb_ax[, .(n_anchor = sum(concordant & human %in% candidates), n_all = sum(concordant)),
                  by = cell_type][order(-n_anchor, -n_all)]
  key <- rankct[1, cell_type]
  message(tt, ": key cell type = ", key, " (", rankct[1, n_anchor], " candidate genes, ",
          rankct[1, n_all], " axis genes concordant)")

  # ---- Fig 3D ---------------------------------------------------------------------------
  sub <- subset(so, subset = cell_type == key)
  Idents(sub) <- "condition"
  fm <- FindMarkers(sub, ident.1 = "STZ", ident.2 = "Control", features = genes, test.use = "wilcox",
                    logfc.threshold = 0, min.pct = 0, verbose = FALSE)
  fm <- as.data.table(fm, keep.rownames = "gene")[, .(gene, cell_log2FC = avg_log2FC, cell_p = p_val,
                                                      pct_STZ = pct.1, pct_Control = pct.2)]
  fm[, cell_BH := p.adjust(cell_p, "BH")]
  st <- merge(fm, pb_ax[cell_type == key, .(gene, pb_logFC = logFC, pb_P = PValue, pb_FDR = FDR,
                                            human, human_logFC, human_P, human_FDR, concordant)], by = "gene")
  st[, fig1_candidate := human %in% candidates]
  fwrite(st[order(pb_P)], file.path(out_dir, paste0("23_anchored_gene_stats_", tt, ".csv")))

  ex <- as.data.table(FetchData(sub, vars = c(genes, "condition", "mouse_id"), layer = "data"))
  long <- melt(ex, id.vars = c("condition", "mouse_id"), variable.name = "gene", value.name = "expr")
  pm <- long[, .(expr = mean(expr), n = .N), by = .(gene, condition, mouse_id)][n >= 10]
  st[, lab := sprintf("%s%s (human %s)\ncell BH %s | pseudobulk P %s", gene,
                      ifelse(fig1_candidate, "*", ""), ifelse(human_logFC > 0, "up", "down"),
                      fmt_p(cell_BH), fmt_p(pb_P))]
  long <- merge(long, st[, .(gene, lab)], by = "gene"); pm <- merge(pm, st[, .(gene, lab)], by = "gene")
  p3d <- ggplot(long, aes(condition, expr, fill = condition)) +
    geom_violin(scale = "width", trim = TRUE, colour = NA, alpha = 0.7) +
    geom_point(data = pm, aes(condition, expr), shape = 21, fill = "white", size = 2) +
    facet_wrap(~ lab, scales = "free_y", nrow = 2) +
    scale_fill_manual(values = cols, guide = "none") +
    labs(x = NULL, y = "normalized expression",
         title = paste0("AMPK-PGC-1a axis in ", key, " (", tt, "); * = Fig 1 candidate gene")) +
    theme_bw(base_size = 9)
  ggsave(file.path(fig_main, paste0("Fig3D_anchored_axis_expression_", tt, ".pdf")), p3d, width = 12, height = 5.5)

  # ---- Fig 3E: score of the Fig 1 candidates on the map, by human direction --------------
  up_h <- deg[gene %in% candidates & logFC > 0, gene]; dn_h <- deg[gene %in% candidates & logFC < 0, gene]
  m_up <- intersect(names(human_of)[human_of %in% up_h], rownames(so))
  m_dn <- intersect(names(human_of)[human_of %in% dn_h], rownames(so))
  sig <- list(candidates_T2D_up = m_up, candidates_T2D_down = m_dn)
  sig <- sig[lengths(sig) > 0]
  so <- AddModuleScore_UCell(so, features = sig, name = "")
  md2 <- as.data.table(so@meta.data)[, c("cell_type", "condition", "mouse_id", names(sig)), with = FALSE]
  pmelt <- melt(md2, id.vars = c("cell_type", "condition", "mouse_id"), variable.name = "signature",
                value.name = "score")
  per_mouse <- pmelt[, .(score = mean(score), n = .N), by = .(signature, cell_type, condition, mouse_id)][n >= 10]
  tst <- per_mouse[, {
    a <- score[condition == "STZ"]; b <- score[condition == "Control"]
    if (length(a) >= 2 && length(b) >= 2) { t <- t.test(a, b); .(delta = mean(a) - mean(b), p = t$p.value) }
    else .(delta = NA_real_, p = NA_real_)
  }, by = .(signature, cell_type)][order(p)]
  tst[, BH := p.adjust(p, "BH"), by = signature]
  fwrite(tst, file.path(out_dir, paste0("23_anchored_score_", tt, ".csv")))
  pmm <- merge(per_mouse, tst[, .(signature, cell_type, p)], by = c("signature", "cell_type"))
  pmm[, lab := paste0(cell_type, "\np = ", fmt_p(p))]
  p3e <- ggplot(pmm, aes(condition, score, fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) + geom_point(shape = 21, size = 1.8, fill = "white") +
    facet_grid(signature ~ lab, scales = "free_y") + scale_fill_manual(values = cols, guide = "none") +
    labs(x = NULL, y = "UCell score", title = paste0("Fig 1 candidate genes on the AMPK map, by cell type (", tt, ")")) +
    theme_bw(base_size = 7)
  ggsave(file.path(fig_main, paste0("Fig3E_anchored_score_", tt, ".pdf")), p3e, width = 12, height = 6)

  # ---- Fig 3F: co-localization of the anchored pair ---------------------------------------
  if (all(pair %in% rownames(so))) {
    pf <- plot_density(so, features = pair, joint = TRUE, reduction = "umap")
    co <- as.data.table(FetchData(so, vars = pair, layer = "data")); co[, cell_type := md$cell_type]
    cop <- co[, .(pct = 100 * mean(get(pair[1]) > 0 & get(pair[2]) > 0), n = .N), by = cell_type][order(-pct)]
    pbar <- ggplot(cop, aes(pct, factor(cell_type, levels = rev(cell_type)))) + geom_col(fill = "#8E5CC8") +
      labs(x = paste0("% cells co-expressing ", pair[1], " + ", pair[2]), y = NULL) + theme_bw(base_size = 9)
    ggsave(file.path(fig_main, paste0("Fig3F_anchored_colocalization_", tt, ".pdf")), pf | pbar,
           width = 17, height = 4.5)
  }

  # ---- Fig 3G: correlation of the anchored pair against a matched background --------------
  e <- as.data.table(FetchData(so, vars = pair, layer = "data"))
  e[, `:=`(cell_type = md$cell_type, condition = md$condition, mouse_id = md$mouse_id)]
  ek <- e[cell_type == key]
  ct_test <- suppressWarnings(cor.test(ek[[pair[1]]], ek[[pair[2]]], method = "pearson"))
  cocells <- which(ek[[pair[1]]] > 0 & ek[[pair[2]]] > 0)
  sp <- if (length(cocells) >= 10) suppressWarnings(cor.test(ek[[pair[1]]][cocells], ek[[pair[2]]][cocells],
                                                            method = "spearman")) else NULL
  pmk <- ek[, .(a = mean(get(pair[1])), b = mean(get(pair[2])), n = .N), by = .(mouse_id, condition)][n >= 10]
  pmt <- if (nrow(pmk) >= 4) cor.test(pmk$a, pmk$b, method = "pearson") else NULL

  cells_k <- which(md$cell_type == key)
  cnt <- GetAssayData(so, layer = "data")[, cells_k, drop = FALSE]
  det <- Matrix::rowMeans(cnt > 0); det <- det[det >= 0.02 & det <= 0.9]
  d1 <- mean(ek[[pair[1]]] > 0); d2 <- mean(ek[[pair[2]]] > 0)
  pool1 <- setdiff(names(det)[abs(det - d1) <= 0.05], pair)
  pool2 <- setdiff(names(det)[abs(det - d2) <= 0.05], pair)
  bg <- if (length(pool1) >= 10 && length(pool2) >= 10) {
    g1 <- sample(pool1, n_bg, TRUE); g2 <- sample(pool2, n_bg, TRUE)
    keep <- g1 != g2; g1 <- g1[keep]; g2 <- g2[keep]
    m <- as.matrix(cnt[unique(c(g1, g2)), , drop = FALSE])
    vapply(seq_along(g1), function(i) suppressWarnings(cor(m[g1[i], ], m[g2[i], ])), 0)
  } else numeric(0)
  bg <- bg[is.finite(bg)]
  emp_p <- if (length(bg)) mean(bg >= unname(ct_test$estimate)) else NA_real_
  corr <- data.table(tissue = tt, cell_type = key, pair = paste(pair, collapse = "-"), n_cells = nrow(ek),
                     pct_coexpr = 100 * length(cocells) / nrow(ek),
                     cell_R = unname(ct_test$estimate), cell_p = ct_test$p.value,
                     coexpr_rho = if (is.null(sp)) NA_real_ else unname(sp$estimate),
                     coexpr_p = if (is.null(sp)) NA_real_ else sp$p.value,
                     pb_r = if (is.null(pmt)) NA_real_ else unname(pmt$estimate),
                     pb_p = if (is.null(pmt)) NA_real_ else pmt$p.value, n_mice = nrow(pmk),
                     bg_mean_R = mean(bg), bg_q95 = if (length(bg)) quantile(bg, 0.95) else NA_real_,
                     empirical_p = emp_p)
  fwrite(corr, file.path(out_dir, paste0("23_anchored_correlation_", tt, ".csv")))
  print(corr)

  g1p <- ggplot(ek, aes(get(pair[1]), get(pair[2]))) +
    geom_jitter(width = 0.02, height = 0.02, alpha = 0.25, size = 0.6, colour = "grey35") +
    geom_smooth(method = "lm", formula = y ~ x, colour = "#C8322F", se = TRUE) +
    labs(x = pair[1], y = pair[2],
         title = sprintf("%s, cell level: R = %.2f (p %s)", key, corr$cell_R, fmt_p(corr$cell_p)),
         subtitle = sprintf("background R: mean %.2f, 95%% %.2f; empirical p = %s",
                            corr$bg_mean_R, corr$bg_q95, fmt_p(corr$empirical_p))) +
    theme_bw(base_size = 9)
  g2p <- ggplot(pmk, aes(a, b, fill = condition)) + geom_point(shape = 21, size = 3) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "grey30", se = FALSE) +
    scale_fill_manual(values = cols) +
    labs(x = paste0("mean ", pair[1]), y = paste0("mean ", pair[2]),
         title = sprintf("per-mouse: r = %.2f (p %s), n = %d", corr$pb_r, fmt_p(corr$pb_p), corr$n_mice)) +
    theme_bw(base_size = 9)
  g3p <- if (length(bg)) ggplot(data.table(R = bg), aes(R)) +
    geom_histogram(bins = 60, fill = "grey75", colour = NA) +
    geom_vline(xintercept = corr$cell_R, colour = "#C8322F", linewidth = 1) +
    labs(x = "cell-level Pearson R", y = "matched random pairs", title = "axis pair vs background") +
    theme_bw(base_size = 9) else plot_spacer()
  ggsave(file.path(fig_main, paste0("Fig3G_anchored_correlation_", tt, ".pdf")), g1p | g2p | g3p,
         width = 15, height = 4.6)
  rm(so); gc(verbose = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_23.txt"))
