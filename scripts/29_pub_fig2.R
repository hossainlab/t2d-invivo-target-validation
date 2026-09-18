# 29_pub_fig2.R - Figure 2 (functional enrichment) redrawn to the publication contract.
#
# Usage: Rscript scripts/29_pub_fig2.R <stage>
#   stage : A | B | C | D | panels | assemble | legend | all
#
# Re-draws only; the enrichment itself stays as script 05 computed it:
#   a  GO over-representation of the candidate set   <- enrich_GO.csv
#   b  KEGG over-representation                      <- enrich_KEGG.csv
#   c  term-similarity network of the enriched sets  <- enrich_GO.csv / enrich_KEGG.csv geneID overlap
#   d  GSEA on the full ranking, Hallmark            <- GSEA_T2DvsControl_Hallmark.csv
#
# The KEGG pathway diagram (Fig2D_pathview/) is rendered by KEGG itself and is left untouched; it is a
# supplied raster, not something this script can restyle.
#
# Outputs: figures/Fig2_enrichment/Fig2<a-d>_*.pdf, Fig2_composite.pdf, Fig2_legend.md
#          results/figure_exports/Fig2_composite.tiff, results/figure_exports/fig2_panels/*.rds

source("scripts/figure_style.R")
suppressPackageStartupMessages({ library(igraph); library(ggraph) })

args  <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) >= 1) args[1] else "all"
set.seed(20260918)

fig_dir   <- "figures/Fig2_enrichment"
panel_dir <- file.path(EXPORT_DIR, "fig2_panels")
bulk      <- "results/bulk"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
keep <- function(p, n, w, h) save_fig(p, n, w, h, dir = fig_dir, keep_dir = panel_dir)

TOP <- 12   # terms shown per ontology / database

ratio_num <- function(x) vapply(strsplit(x, "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]), 0)

# ===================================================================================================
# a - GO over-representation
# ===================================================================================================
if (stage %in% c("A", "panels", "all")) {
  go <- fread(file.path(bulk, "enrich_GO.csv"))
  go[, gene_ratio := ratio_num(GeneRatio)]
  d <- go[order(pvalue)][, head(.SD, TOP), by = ONTOLOGY]
  d[, term := factor(wrap_term(tidy_term(Description), 38),
                     levels = rev(wrap_term(tidy_term(Description), 38)))]
  d[, ONTOLOGY := factor(ONTOLOGY, levels = c("BP", "CC", "MF"),
                         labels = c("Biological process", "Cellular component", "Molecular function"))]

  pA <- ggplot(d, aes(gene_ratio, term)) +
    geom_segment(aes(x = 0, xend = gene_ratio, yend = term), colour = "grey80", linewidth = 0.25) +
    geom_point(aes(size = Count, colour = -log10(p.adjust))) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_colour_viridis_c(option = "viridis", direction = -1,
                           name = expression(-log[10]~italic(q)),
                           guide = guide_colourbar(barheight = unit(14, "mm"),
                                                   barwidth = unit(2, "mm"))) +
    scale_size_continuous(range = c(0.6, 2.4), name = "Genes", breaks = scales::breaks_pretty(3)) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(x = "Gene ratio", y = NULL) +
    theme_pub() +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90),
          axis.text.y = element_text(size = BASE - 2, lineheight = 0.9),
          panel.spacing.y = unit(1.5, "mm"))
  keep(pA, "Fig2a_GO", 100, 118)
}

# ===================================================================================================
# b - KEGG over-representation
# ===================================================================================================
if (stage %in% c("B", "panels", "all")) {
  kg <- fread(file.path(bulk, "enrich_KEGG.csv"))
  kg[, gene_ratio := ratio_num(GeneRatio)]
  d <- kg[p.adjust < 0.05][order(pvalue)][1:min(.N, TOP)]
  d[, term := factor(wrap_term(tidy_term(Description), 34),
                     levels = rev(wrap_term(tidy_term(Description), 34)))]

  pB <- ggplot(d, aes(gene_ratio, term)) +
    geom_segment(aes(x = 0, xend = gene_ratio, yend = term), colour = "grey80", linewidth = 0.25) +
    geom_point(aes(size = Count, colour = -log10(p.adjust))) +
    scale_colour_viridis_c(option = "viridis", direction = -1,
                           name = expression(-log[10]~italic(q)),
                           guide = guide_colourbar(barheight = unit(14, "mm"),
                                                   barwidth = unit(2, "mm"))) +
    scale_size_continuous(range = c(0.8, 2.6), name = "Genes", breaks = scales::breaks_pretty(3)) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(x = "Gene ratio", y = NULL) +
    theme_pub() +
    theme(axis.text.y = element_text(size = BASE - 2, lineheight = 0.9))
  keep(pB, "Fig2b_KEGG", 89, 62)
}

# ===================================================================================================
# c - term-similarity network (the enrichment map, rebuilt from gene overlap)
# ===================================================================================================
if (stage %in% c("C", "panels", "all")) {
  # the reference framework (Xu et al. Fig. 2C) builds this network from the top 30 enriched KEGG
  # terms, not from GO; only 9 KEGG pathways reach q < 0.05 here, so the top 30 by nominal P are used
  # and significance is carried by node colour
  NTERM <- 30
  go <- fread(file.path(bulk, "enrich_KEGG.csv"))[order(pvalue)][1:NTERM]
  sets <- strsplit(go$geneID, "/")
  names(sets) <- tidy_term(go$Description)
  n <- length(sets)
  jac <- matrix(0, n, n, dimnames = list(names(sets), names(sets)))
  for (i in seq_len(n)) for (j in seq_len(n)) if (i < j) {
    jac[i, j] <- length(intersect(sets[[i]], sets[[j]])) / length(union(sets[[i]], sets[[j]]))
  }
  ed <- as.data.table(which(jac >= 0.5, arr.ind = TRUE))
  ed[, `:=`(from = rownames(jac)[row], to = colnames(jac)[col], w = jac[cbind(row, col)])]
  vt <- data.frame(name = names(sets), q = -log10(go$p.adjust), size = go$Count)
  # terms that share nothing with any other term carry no information in a network and only push the
  # connected component into a corner; they are listed in panel a instead
  connected <- unique(c(ed$from, ed$to))
  vt <- vt[vt$name %in% connected, ]
  g <- graph_from_data_frame(ed[, .(from, to, w)], directed = FALSE, vertices = vt)
  lay <- create_layout(g, layout = "fr")

  pC <- ggraph(lay) +
    geom_edge_link(aes(width = w), colour = "grey78", alpha = 0.8, show.legend = FALSE) +
    scale_edge_width(range = c(0.15, 0.8)) +
    geom_node_point(aes(size = size, colour = q)) +
    geom_node_text(aes(label = wrap_term(name, 18)), size = (BASE - 3) * ppt, family = FONT,
                   repel = TRUE, lineheight = 0.88, max.overlaps = Inf, segment.size = 0.15,
                   segment.colour = "grey60", box.padding = 0.45, force = 8, max.iter = 20000) +
    scale_colour_viridis_c(option = "viridis", direction = -1, breaks = scales::breaks_pretty(3),
                           name = expression(-log[10]~italic(q)),
                           guide = guide_colourbar(barheight = unit(2, "mm"),
                                                   barwidth = unit(14, "mm"),
                                                   title.position = "top")) +
    scale_size_continuous(range = c(1, 3.4), name = "Genes",
                          guide = guide_legend(title.position = "top", nrow = 1)) +
    coord_equal(clip = "off") +
    theme_void_pub() +
    theme(legend.position = "bottom", legend.box = "horizontal",
          legend.spacing.x = unit(3, "mm"))
  keep(pC, "Fig2c_kegg_network", 110, 95)
}

# ===================================================================================================
# d - GSEA, Hallmark
# ===================================================================================================
if (stage %in% c("D", "panels", "all")) {
  gs <- fread(file.path(bulk, "GSEA_T2DvsControl_Hallmark.csv"))
  d <- gs[padj < 0.25][order(NES)]
  if (nrow(d) > 20) d <- rbind(head(d, 10), tail(d, 10))
  d[, term := factor(wrap_term(tidy_term(pathway), 34), levels = wrap_term(tidy_term(pathway), 34))]
  d[, dir := fifelse(NES > 0, "Up in T2D", "Down in T2D")]

  pD <- ggplot(d, aes(NES, term, fill = dir)) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_col(width = 0.68) +
    geom_point(aes(size = size), shape = 21, colour = "grey20", stroke = 0.2, show.legend = TRUE) +
    scale_fill_manual(values = c(`Up in T2D` = col_dir[["Up"]], `Down in T2D` = col_dir[["Down"]]),
                      name = NULL) +
    scale_size_continuous(range = c(0.6, 2.2), name = "Set size", breaks = scales::breaks_pretty(3),
                          guide = guide_legend(nrow = 1)) +
    labs(x = "Normalised enrichment score", y = NULL) +
    theme_pub() +
    theme(axis.text.y = element_text(size = BASE - 2, lineheight = 0.9),
          legend.position = "bottom", legend.box = "vertical",
          legend.spacing.y = unit(0.5, "mm"), legend.justification = "left")
  keep(pD, "Fig2d_GSEA_hallmark", 89, 78)
}

# ===================================================================================================
# assemble
# ===================================================================================================
if (stage %in% c("assemble", "all")) {
  need <- c("Fig2a_GO", "Fig2b_KEGG", "Fig2c_kegg_network", "Fig2d_GSEA_hallmark")
  ps <- lapply(need, function(n) {
    f <- file.path(panel_dir, paste0(n, ".rds"))
    if (!file.exists(f)) stop("missing panel ", n, "; run stage panels first")
    readRDS(f)$plot
  })
  design <- "
AAAAAACCCCCC
AAAAAACCCCCC
AAAAAACCCCCC
AAAAAABBBBBB
AAAAAABBBBBB
DDDDDDBBBBBB
DDDDDDBBBBBB
"
  comp <- ps[[1]] + ps[[2]] + ps[[3]] + ps[[4]] +
    plot_layout(design = design) +
    plot_annotation(tag_levels = list(tag_letters(4)[c(1, 2, 3, 4)])) &
    theme(plot.tag = element_text(size = BASE + 1, face = "bold", family = FONT))
  save_fig(comp, "Fig2_composite", W_2COL, 200, dir = fig_dir, tiff = TRUE)
}

# ===================================================================================================
# legend
# ===================================================================================================
if (stage %in% c("legend", "panels", "all")) {
  go <- fread(file.path(bulk, "enrich_GO.csv"))
  kg <- fread(file.path(bulk, "enrich_KEGG.csv"))
  gs <- fread(file.path(bulk, "GSEA_T2DvsControl_Hallmark.csv"))
  isum <- fread(file.path(bulk, "intersection_summary.csv"))
  ncand <- isum[item == "triple_intersection", value]

  l <- c(
    "# Figure 2 legend - functional enrichment of the candidate set",
    "",
    "**Fig. 2 | Pathways enriched in the T2D candidate genes.**",
    sprintf("**a**, GO over-representation of the %s candidate genes from Fig. 1g, top %d terms per ontology by P. Point size, genes from the candidate set annotated to the term; colour, BH-adjusted q.",
            ncand, TOP),
    sprintf("**b**, KEGG over-representation, all %d pathways at q < 0.05.", kg[p.adjust < 0.05, .N]),
    "**c**, Pathway interaction network of the 30 most enriched KEGG terms; nodes are pathways, edges join pathways sharing at least half of their candidate genes (Jaccard >= 0.5). Node colour, BH-adjusted q; size, candidate genes in the pathway.",
    sprintf("**d**, GSEA of the full T2D vs control ranking against the Hallmark collection, sets at q < 0.25 (%d of 50). Bars, normalised enrichment score; points, set size.",
            gs[padj < 0.25, .N]),
    "",
    "## Caveats that belong in the text, not the figure",
    "",
    sprintf("- **a**,**b** test %s genes against a %s-gene background; with a set this small, an over-representation q value is fragile and a single gene changes several terms.",
            ncand, sub(".*/", "", kg$BgRatio[1])),
    sprintf("- **The enriched KEGG pathways are not independent findings.** All %d pathways at q < 0.05 are driven by the same four genes, IRS2, PPARGC1A, CCND1 and IGF1, recombined: AMPK signalling has all four, and %d of the %d pathways rest on a single pair. Apelin signalling and alcoholic liver disease are both PPARGC1A/CCND1; melanoma, glioma, prostate cancer, breast cancer, endocrine resistance, focal adhesion, integrin signalling and proteoglycans in cancer are all CCND1/IGF1, identically. The list should be read as one signal seen through %d annotation sets, not as %d separate pathway results, and the cancer maps in particular are annotation artefacts of CCND1 rather than evidence of cancer biology.",
            kg[p.adjust < 0.05, .N],
            kg[p.adjust < 0.05][vapply(strsplit(geneID, "/"), length, 1L) == 2, .N],
            kg[p.adjust < 0.05, .N], kg[p.adjust < 0.05, .N], kg[p.adjust < 0.05, .N]),
    "- The candidate set was defined partly by a hyperglycaemia gene set (Fig. 1g), which shares annotation with the GO and KEGG terms tested here, so **a** and **b** are not independent of that selection step.",
    "- **c** edges are gene-overlap similarity between pathways, not curated pathway-pathway relationships; with only 19 candidate genes most KEGG pathways share the same few genes, which is why the overlap threshold is high.",
    "",
    "## Statistics to quote",
    "",
    "| Quantity | Value |",
    "|---|---|",
    sprintf("| Candidate genes tested | %s |", ncand),
    sprintf("| GO terms at q < 0.05 | %d |", go[p.adjust < 0.05, .N]),
    sprintf("| Top GO term | %s (%s) |", tidy_term(go[order(pvalue)][1, Description]),
            p_short(go[order(pvalue)][1, pvalue])),
    sprintf("| KEGG pathways at q < 0.05 | %d |", kg[p.adjust < 0.05, .N]),
    sprintf("| Top KEGG pathway | %s (%s) |", kg[order(pvalue)][1, Description],
            p_short(kg[order(pvalue)][1, pvalue])),
    sprintf("| Hallmark sets at q < 0.25 | %d of 50 |", gs[padj < 0.25, .N]),
    sprintf("| Strongest Hallmark NES | %s, NES = %s (q = %s) |", tidy_term(gs[which.max(abs(NES)), pathway]),
            mfmt(gs[which.max(abs(NES)), NES]), mfmt(gs[which.max(abs(NES)), padj], 3))
  )
  writeLines(l, file.path(fig_dir, "Fig2_legend.md"))
  message("wrote ", file.path(fig_dir, "Fig2_legend.md"))
}
