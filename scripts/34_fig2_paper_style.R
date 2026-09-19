# 34_fig2_paper_style.R - Figure 2 drawn in the visual style of the reference paper.
#
# Usage: Rscript scripts/34_fig2_paper_style.R <stage>
#   stage : A | B | C | D | panels | assemble | legend | all
#
# Reference: Xu et al., Phytomedicine 154 (2026) 158050, Fig. 2. This reproduces that figure's own
# visual language with this project's data, panel for panel:
#   A  GO enrichment bars, one column per ontology, BP blue / CC yellow / MF purple. Top row, the
#      -log10 adjusted P bar with the term written inside it; bottom row, the same bar with a Rich
#      Factor line and points on a secondary top axis
#   B  KEGG bubble chart, x = GeneRatio, size = Count, colour = p.adjust on a red-to-blue ramp, with
#      the selected pathway boxed in red
#   C  pathway interaction network of the top 30 KEGG terms, nodes coloured by p.adjust and sized by
#      gene count, grey edges, selected pathway boxed in red
#   D  the pathview rendering of the selected KEGG map, as supplied by script 05
#
# This is deliberately NOT the figure_style.R contract used by scripts 27-32. It is the reference
# paper's style, kept as a separate deliverable in figures/Fig2_paper_style/ so the Nature/Cell version
# in figures/Fig2_enrichment/ is untouched. Pick one for submission; do not ship both.
#
# Deviations from the reference forced by this dataset, all of them stated in Fig2_paper_style_legend.md:
#   - A has no CC column. Over-representation of the 8 candidate genes returns BP and MF terms only;
#     no cellular-component term is returned at all, so the paper's yellow column has nothing to draw.
#   - B and C rest on 4 genes. Only 4 of the 8 candidates (IRS2, PPARGC1A, CCND1, IGF1) map into any
#     KEGG pathway, so every GeneRatio here has denominator 4 and the pathways are recombinations of
#     the same four genes rather than independent results.
#   - D shows AMPK signalling (hsa04152), the map selected in docs/pathway_selection.md, standing in
#     for the paper's regulation of actin cytoskeleton.
#
# Outputs: figures/Fig2_paper_style/Fig2<A-D>_*.pdf, Fig2_paper_style.pdf, Fig2_paper_style_legend.md
#          results/figure_exports/Fig2_paper_style.tiff

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(igraph); library(ggraph); library(ggrepel); library(png); library(grid)
})
options(stringsAsFactors = FALSE)

args  <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) >= 1) args[1] else "all"
set.seed(20260918)

fig_dir   <- "figures/Fig2_paper_style"
panel_dir <- file.path("results/figure_exports", "fig2paper_panels")
bulk      <- "results/bulk"
for (d in c(fig_dir, panel_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

TOP_GO   <- 10          # terms per ontology in panel A, as in the reference
TOP_KEGG <- 30          # bubbles in B and nodes in C, as in the reference
SEL_ID   <- "hsa04152"  # the selected map, boxed in B and C and drawn in D
PATHVIEW <- "figures/Fig2_enrichment/Fig2D_pathview/hsa04152.T2D_vs_Control_logFC.png"

# the paper's palette: one hue per ontology in A, a red-to-blue p.adjust ramp in B and C
col_onto <- c(BP = "#9ECAE1", CC = "#FDBE85", MF = "#BCBDDC")
col_line <- c(BP = "#08306B", CC = "#E6550D", MF = "#54278F")
col_q_lo <- "#E8362B"   # most significant
col_q_hi <- "#3B7DD8"   # least significant

theme_baseR <- function(base = 8)
  theme_bw(base_size = base) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
        axis.ticks = element_line(colour = "black", linewidth = 0.3),
        plot.title = element_text(hjust = 0.5, face = "bold", size = base + 0.5),
        plot.margin = margin(2, 3, 2, 2))

keep <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, units = "mm",
         device = cairo_pdf)
  saveRDS(list(plot = p, width = w, height = h), file.path(panel_dir, paste0(name, ".rds")))
  message("wrote ", name, " (", w, " x ", h, " mm)")
  invisible(p)
}

wrap_term <- function(x, width = 42)
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), "", USE.NAMES = FALSE)

ratio_num <- function(x) {
  z <- strsplit(as.character(x), "/", fixed = TRUE)
  out <- vapply(z, function(p) if (length(p) == 2L) as.numeric(p[1]) / as.numeric(p[2]) else NA_real_, 0)
  # enrich_KEGG.csv has been damaged once by a round trip through Excel, which turned "4/4" into
  # "4-Apr" and left every GeneRatio unparseable and every bubble silently dropped. Fail loudly.
  if (anyNA(out)) stop("GeneRatio will not parse (", paste(utils::head(unique(x[is.na(out)]), 3),
                       collapse = ", "), "). Re-run scripts/05_bulk_enrichment.R; do not open the ",
                       "enrichment CSVs in Excel.")
  out
}

# ===================================================================================================
# A - GO bars, one column per ontology, as in the paper
# ===================================================================================================
if (stage %in% c("A", "panels", "all")) {
  go <- fread(file.path(bulk, "enrich_GO.csv"))[p.adjust < 0.05]
  d  <- go[order(p.adjust), head(.SD, TOP_GO), by = ONTOLOGY]
  d[, `:=`(negl = -log10(p.adjust), term = wrap_term(Description, 46))]
  onts <- intersect(c("BP", "CC", "MF"), unique(d$ONTOLOGY))

  # one column per ontology; two rows, as the reference draws it
  mk_top <- function(o) {
    x <- d[ONTOLOGY == o][order(negl)]
    x[, term := factor(term, levels = term)]
    ggplot(x, aes(negl, term)) +
      geom_col(fill = col_onto[[o]], width = 0.78) +
      geom_text(aes(x = max(negl) * 0.012, label = term), hjust = 0, size = 1.55,
                lineheight = 0.85) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
      labs(title = o, x = expression(-log[10]*"(adjust.p)"), y = NULL) +
      theme_baseR() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  # bottom row: the same bars with the Rich Factor line on a secondary top axis
  mk_bot <- function(o) {
    x <- d[ONTOLOGY == o][order(negl)]
    x[, term := factor(term, levels = term)]
    sc <- max(x$negl) / max(x$RichFactor)   # map Rich Factor onto the bar axis
    ggplot(x, aes(y = term)) +
      geom_col(aes(x = negl), fill = col_onto[[o]], width = 0.78) +
      geom_line(aes(x = RichFactor * sc, group = 1), colour = col_line[[o]], linewidth = 0.3) +
      geom_point(aes(x = RichFactor * sc), colour = col_line[[o]], size = 0.7) +
      scale_x_continuous(
        expand = expansion(mult = c(0, 0.04)),
        sec.axis = sec_axis(~ . / sc, name = "Rich Factor")) +
      labs(x = expression(-log[10]*"(adjust.p)"), y = NULL) +
      theme_baseR() +
      theme(axis.text.y = element_text(size = 4.6, lineheight = 0.85))
  }

  pA <- wrap_elements(full =
    wrap_plots(lapply(onts, mk_top), nrow = 1) /
    wrap_plots(lapply(onts, mk_bot), nrow = 1))
  keep(pA, "Fig2A_GO", 150, 108)
}

# ===================================================================================================
# B - KEGG bubble chart, as in the paper
# ===================================================================================================
if (stage %in% c("B", "panels", "all")) {
  kg <- fread(file.path(bulk, "enrich_KEGG.csv"))
  kg[, gene_ratio := ratio_num(GeneRatio)]
  d <- kg[order(pvalue)][1:min(.N, TOP_KEGG)]
  d[, term := factor(Description, levels = rev(Description))]
  # the reference boxes its chosen pathway in red; ours is the map selected in pathway_selection.md.
  # The box is drawn after the points, on the discrete axis by row index, so it does not train the
  # y scale from a one-row subset and silently reorder the panel.
  sel_row <- which(levels(d$term) == d[ID == SEL_ID, Description])
  sel_x   <- d[ID == SEL_ID, gene_ratio]
  ylab_col  <- ifelse(seq_along(levels(d$term)) == sel_row, "red", "black")
  ylab_face <- ifelse(seq_along(levels(d$term)) == sel_row, "bold", "plain")

  pB <- ggplot(d, aes(gene_ratio, term)) +
    geom_point(aes(size = Count, colour = p.adjust)) +
    annotate("rect", xmin = sel_x - 0.075, xmax = sel_x + 0.075,
             ymin = sel_row - 0.48, ymax = sel_row + 0.48,
             fill = NA, colour = "red", linewidth = 0.35) +
    scale_colour_gradient(low = col_q_lo, high = col_q_hi, name = "p.adjust",
                          guide = guide_colourbar(barwidth = unit(2.4, "mm"),
                                                  barheight = unit(18, "mm"))) +
    scale_size_continuous(range = c(0.9, 3.4), name = "Count", breaks = scales::breaks_pretty(4)) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.10))) +
    labs(x = "GeneRatio", y = NULL) +
    theme_baseR() +
    theme(axis.text.y = element_text(size = 5.2, colour = ylab_col, face = ylab_face),
          legend.title = element_text(size = 6.5), legend.text = element_text(size = 5.8),
          legend.key.size = unit(3, "mm"))
  keep(pB, "Fig2B_KEGG", 110, 108)
}

# ===================================================================================================
# C - pathway interaction network of the top 30 KEGG terms, as in the paper
# ===================================================================================================
if (stage %in% c("C", "panels", "all")) {
  kg <- fread(file.path(bulk, "enrich_KEGG.csv"))[order(pvalue)][1:TOP_KEGG]
  sets <- strsplit(kg$geneID, "/"); names(sets) <- kg$Description
  n <- length(sets)
  jac <- matrix(0, n, n, dimnames = list(names(sets), names(sets)))
  for (i in seq_len(n)) for (j in seq_len(n)) if (i < j)
    jac[i, j] <- length(intersect(sets[[i]], sets[[j]])) / length(union(sets[[i]], sets[[j]]))
  ed <- as.data.table(which(jac >= 0.5, arr.ind = TRUE))
  ed[, `:=`(from = rownames(jac)[row], to = colnames(jac)[col], w = jac[cbind(row, col)])]

  vt <- data.frame(name = names(sets), q = kg$p.adjust, size = kg$Count,
                   sel = kg$ID == SEL_ID)
  g   <- graph_from_data_frame(ed[, .(from, to, w)], directed = FALSE, vertices = vt)
  lay <- create_layout(g, layout = "fr")

  pC <- ggraph(lay) +
    geom_edge_link(aes(width = w), colour = "grey78", alpha = 0.85, show.legend = FALSE) +
    scale_edge_width(range = c(0.15, 0.7)) +
    geom_node_point(aes(size = size, colour = q)) +
    # label colour and weight are passed as plain vectors in node order, not mapped: the node points
    # already own the colour scale, and a second colour mapping would collide with it
    # plain repelled text with a white halo, as the reference draws it. geom_node_label() puts a
    # bordered box round every term, which the reference has only on its selected pathway; here the
    # selection is carried by red bold text instead.
    geom_node_text(aes(label = wrap_term(name, 20)),
                   colour = ifelse(lay$sel, "red", "black"),
                   fontface = ifelse(lay$sel, "bold", "plain"),
                   bg.colour = "white", bg.r = 0.12,
                   size = 1.4, lineheight = 0.85, repel = TRUE, max.overlaps = Inf,
                   segment.size = 0.12, segment.colour = "grey55",
                   box.padding = 0.34, point.padding = 0.16, force = 14, force_pull = 0.4,
                   max.iter = 60000, show.legend = FALSE) +
    scale_colour_gradient(low = col_q_lo, high = col_q_hi, name = "p.adjust",
                          guide = guide_colourbar(barwidth = unit(2.4, "mm"),
                                                  barheight = unit(16, "mm"))) +
    scale_size_continuous(range = c(1.2, 4), name = "number of genes",
                          breaks = scales::breaks_pretty(4)) +
    # no coord_equal: the force-directed layout is wider than it is tall, and locking the aspect
    # left the node cloud as a thin band across the middle of the composite cell
    coord_cartesian(clip = "off") +
    theme_void(base_size = 8) +
    theme(legend.title = element_text(size = 6.5), legend.text = element_text(size = 5.8),
          legend.key.size = unit(3, "mm"), plot.margin = margin(4, 4, 4, 4))
  keep(pC, "Fig2C_kegg_network", 140, 125)
}

# ===================================================================================================
# D - the pathview map of the selected pathway, as in the paper
# ===================================================================================================
if (stage %in% c("D", "panels", "all")) {
  if (!file.exists(PATHVIEW))
    stop("missing ", PATHVIEW, "; run scripts/05_bulk_enrichment.R to render the pathview map")
  img <- readPNG(PATHVIEW)
  pD  <- wrap_elements(full = rasterGrob(img, interpolate = TRUE))
  keep(pD, "Fig2D_pathview", 130, 118)
}

# ===================================================================================================
# assemble - the paper's layout: A|B on row 1, C|D on row 2
# ===================================================================================================
if (stage %in% c("assemble", "all")) {
  need <- c("Fig2A_GO", "Fig2B_KEGG", "Fig2C_kegg_network", "Fig2D_pathview")
  ps <- lapply(need, function(n) {
    f <- file.path(panel_dir, paste0(n, ".rds"))
    if (!file.exists(f)) stop("missing panel ", n, "; run stage panels first")
    readRDS(f)$plot
  })
  design <- "
AAAAAAAABBBBBB
AAAAAAAABBBBBB
AAAAAAAABBBBBB
CCCCCCCDDDDDDD
CCCCCCCDDDDDDD
CCCCCCCDDDDDDD
"
  comp <- ps[[1]] + ps[[2]] + ps[[3]] + ps[[4]] +
    plot_layout(design = design) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 13, face = "bold"))
  out <- file.path(fig_dir, "Fig2_paper_style")
  ggsave(paste0(out, ".pdf"), comp, width = 260, height = 230, units = "mm", device = cairo_pdf)
  tif <- file.path("results/figure_exports", "Fig2_paper_style.tiff")
  unlink(tif)
  ggsave(tif, comp, width = 260, height = 230, units = "mm", dpi = 400, bg = "white",
         compression = "lzw")
  message("wrote ", out, ".pdf and ", tif)
}

# ===================================================================================================
if (stage %in% c("legend", "panels", "all")) {
  go   <- fread(file.path(bulk, "enrich_GO.csv"))
  kg   <- fread(file.path(bulk, "enrich_KEGG.csv"))
  isum <- fread(file.path(bulk, "intersection_summary.csv"))
  ncand <- isum[item == "triple_intersection", value]
  ot    <- sort(unique(go[p.adjust < 0.05, ONTOLOGY]))
  onts  <- if (length(ot) > 1) paste(paste(utils::head(ot, -1), collapse = ", "), "and",
                                     utils::tail(ot, 1)) else ot
  sel   <- kg[ID == SEL_ID]
  kegg_n <- kg[p.adjust < 0.05, .N]
  npair  <- kg[p.adjust < 0.05][vapply(strsplit(geneID, "/"), length, 1L) == 2, .N]
  kgenes <- sort(unique(unlist(strsplit(kg$geneID, "/"))))

  l <- c(
    "# Figure 2, reference-paper style - legend",
    "",
    "Drawn to match Xu et al., *Phytomedicine* 154 (2026) 158050, Fig. 2, with this project's data.",
    "The Nature/Cell-contract version of the same figure is in `figures/Fig2_enrichment/`.",
    "**Ship one or the other, not both.**",
    "",
    sprintf("**Fig. 2. Functional enrichment analysis of intersecting genes.** (A) Results of GO enrichment analysis. Blue represents biological process (BP) terms, yellow represents cellular component (CC) terms, and purple represents molecular function (MF) terms. (B) Bubble chart displaying the results of KEGG analysis. (C) Pathway interaction network of the top %d enriched KEGG terms. (D) Pathway diagram generated using the pathview package, illustrating signal activation within %s.",
            TOP_KEGG, sel$Description),
    "",
    "## Where this data forces a deviation from the reference",
    "",
    sprintf("- **(A) has no CC column.** Over-representation of the %s candidate genes returns %s terms only; not one cellular-component term is returned, so the paper's yellow column has nothing to draw. The colour key is kept as the reference states it.",
            ncand, onts),
    sprintf("- **(B) and (C) rest on %d genes, not 159.** Only %d of the %s candidates (%s) map into any KEGG pathway, so every GeneRatio on this figure has denominator %d.",
            length(kgenes), length(kgenes), ncand, paste(kgenes, collapse = ", "), length(kgenes)),
    sprintf("- **The %d pathways at adjusted P < 0.05 are not %d independent findings.** They are recombinations of those same %d genes: %s has all %d, and %d of the %d rest on a single pair. The cancer maps in the list are annotation artefacts of CCND1, not evidence of cancer biology. Read the panel as one signal seen through %d annotation sets.",
            kegg_n, kegg_n, length(kgenes), sel$Description, length(kgenes), npair, kegg_n, kegg_n),
    sprintf("- **(C) shows the top %d KEGG terms by nominal P,** as the reference does, because only %d reach adjusted P < 0.05; significance is carried by node colour. Edges join pathways sharing at least half of their candidate genes (Jaccard >= 0.5), which is gene-overlap similarity, not a curated pathway-pathway relation.",
            TOP_KEGG, kegg_n),
    sprintf("- **(D) shows %s (%s)**, the map selected in `docs/pathway_selection.md`, standing in for the paper's regulation of actin cytoskeleton. Node colour is the T2D vs control log2 fold change, so the panel shows the *transcript* state of the map, which for this pathway is the finding: the AMPK subunit transcript rises while its outputs fall.",
            sel$Description, SEL_ID),
    "- The candidate set was defined partly by a hyperglycaemia gene set (Fig. 1G), which shares annotation with the GO and KEGG terms tested here, so (A) and (B) are not independent of that selection step.",
    "",
    "## Statistics to quote",
    "",
    "| Quantity | Value |",
    "|---|---|",
    sprintf("| Candidate genes tested | %s |", ncand),
    sprintf("| Background genes | %s (GO), %s (KEGG) |", sub(".*/", "", go$BgRatio[1]),
            sub(".*/", "", kg$BgRatio[1])),
    sprintf("| GO terms at adjusted P < 0.05 | %d (%s) |", go[p.adjust < 0.05, .N], onts),
    sprintf("| Top GO term | %s, adjusted P = %s |", go[order(pvalue)][1, Description],
            format(signif(go[order(pvalue)][1, p.adjust], 2), scientific = TRUE)),
    sprintf("| KEGG pathways at adjusted P < 0.05 | %d |", kegg_n),
    sprintf("| Selected pathway | %s (%s), GeneRatio %s, adjusted P = %s |", sel$Description, SEL_ID,
            sel$GeneRatio, format(signif(sel$p.adjust, 2), scientific = TRUE)),
    sprintf("| Genes driving the KEGG result | %s |", paste(kgenes, collapse = ", "))
  )
  writeLines(l, file.path(fig_dir, "Fig2_paper_style_legend.md"))
  message("wrote ", file.path(fig_dir, "Fig2_paper_style_legend.md"))
}
