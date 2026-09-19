# 05_bulk_enrichment.R — Phase 5 of docs/analysis_plan.md
# ORA (GO BP/CC/MF, KEGG) on intersect genes with filtered-gene universe; pathway network;
# pathview of lead KEGG pathway; GSEA (fgsea) on the T2D vs Control ranked list.
#
# Outputs: results/bulk/enrich_GO.csv, enrich_KEGG.csv, enrich_DEG_up_down_GO.csv,
#          GSEA_<contrast>_<collection>.csv
# Figures: figures/bulk/Fig2A_GO.pdf, Fig2B_KEGG.pdf, Fig2C_emap.pdf; the pathview raster and its
#          KGML go to results/pathway_selection/pathview/,
#          05_GSEA_hallmark.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(fgsea)
  library(ggplot2)
  library(patchwork)
  library(pathview)  # must be attached: pathview() calls data(bods) internally
})
set.seed(20260914)
out_dir <- "results/bulk"; fig_dir <- "results/supplementary_figures/bulk"
fig_main <- "results/supplementary_figures/enrichment"  # draft panels; figures/Fig2 is owned by
# scripts/29_pub_fig2.R. The pathview diagram is the exception and still goes to figures/ (see below).
for (d in c(fig_dir, fig_main)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Remove outputs of previous runs first: files are only written when enrichment is found,
# so a null result must not leave stale files behind.
unlink(c(file.path(out_dir, c("enrich_GO.csv", "enrich_KEGG.csv", "enrich_DEG_up_down_GO.csv")),
         Sys.glob(file.path(out_dir, "GSEA_*.csv")),
         file.path(fig_main, c("Fig2A_GO.pdf", "Fig2B_KEGG.pdf", "Fig2C_emap_KEGG.pdf", "Fig2C_emap_GOBP.pdf",
                               "Fig2_GSEA_hallmark.pdf")),
         file.path(fig_dir, "05_GO_DEG_up_down.pdf")))
# the pathview diagram is a rendered asset, not a drawn panel: KEGG produces the raster and no figure
# script redraws it. It therefore lives with the other analysis outputs, and the figure script that
# needs it reads it from there. figures/ holds only what a figure script writes.
pv_main <- "results/pathway_selection"
unlink(file.path(pv_main, "pathview"), recursive = TRUE)

deg <- fread(file.path(out_dir, "DEG_T2D_vs_Control_all.csv"))
ints <- fread(file.path(out_dir, "intersect_genes.csv"))

to_entrez <- function(sym) {
  sym <- unique(sym[!is.na(sym) & nzchar(sym)])
  if (!length(sym)) return(character(0))  # empty gene list (e.g. no DEG ∩ hub genes)
  e <- suppressMessages(mapIds(org.Hs.eg.db, sym, "ENTREZID", "SYMBOL"))
  unique(na.omit(e))
}
universe_entrez <- to_entrez(deg$gene)
genes_entrez    <- to_entrez(ints$gene)
message("Intersect genes: ", nrow(ints), " (", length(genes_entrez), " Entrez); universe: ", length(universe_entrez))

# ---------------- GO ORA ---------------------------------------------------------------
run_ora <- length(genes_entrez) >= 5
if (!run_ora) message("Only ", length(genes_entrez), " intersect genes: over-representation analysis skipped; GSEA below is used")
ego <- if (run_ora) enrichGO(genes_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL",
                universe = universe_entrez, pAdjustMethod = "BH", pvalueCutoff = 0.05,
                qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = TRUE) else NULL
ego_dt <- if (!is.null(ego)) as.data.table(ego@result) else data.table(ONTOLOGY = character(), Description = character(), p.adjust = numeric())
if (nrow(ego_dt)) fwrite(ego_dt, file.path(out_dir, "enrich_GO.csv"))
ego_sig <- ego_dt[p.adjust < 0.05]
message("GO terms adj.P<0.05: ", nrow(ego_sig))

if (nrow(ego_sig)) {
  top_go <- ego_sig[order(p.adjust), head(.SD, 10), by = ONTOLOGY]
  top_go[, Description := factor(Description, levels = rev(unique(Description)))]
  pA <- ggplot(top_go, aes(-log10(p.adjust), Description, fill = ONTOLOGY)) +
    geom_col(width = 0.75) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c(BP = "#3B7DD8", CC = "#E0A100", MF = "#8E5CC8")) +
    labs(x = expression(-log[10]~adjusted~P), y = NULL, title = "GO enrichment") +
    theme_bw(base_size = 11) + theme(legend.position = "none")
  ggsave(file.path(fig_main, "Fig2A_GO.pdf"), pA, width = 8.5, height = 8)
}

# ---------------- KEGG ORA -------------------------------------------------------------
ekg <- if (!run_ora) NULL else tryCatch(enrichKEGG(genes_entrez, organism = "hsa", universe = universe_entrez,
                           pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10),
                error = function(e) { message("KEGG failed: ", conditionMessage(e)); NULL })
if (!is.null(ekg)) {
  ekg <- setReadable(ekg, org.Hs.eg.db, keyType = "ENTREZID")
  ekg_dt <- as.data.table(ekg@result)
  fwrite(ekg_dt, file.path(out_dir, "enrich_KEGG.csv"))
  message("KEGG pathways adj.P<0.05: ", ekg_dt[p.adjust < 0.05, .N])
  if (ekg_dt[p.adjust < 0.05, .N]) {
    pB <- dotplot(ekg, showCategory = 20) + ggtitle("KEGG enrichment")
    ggsave(file.path(fig_main, "Fig2B_KEGG.pdf"), pB, width = 7.5, height = 7)

    # Pathway interaction network (top 30 GO+KEGG terms by adj.P)
    ekg_sim <- pairwise_termsim(ekg)
    pC <- emapplot(ekg_sim, showCategory = min(30, ekg_dt[p.adjust < 0.05, .N])) +
      ggtitle("KEGG pathway interaction network")
    ggsave(file.path(fig_main, "Fig2C_emap_KEGG.pdf"), pC, width = 9, height = 8)

    # pathview for lead non-disease pathway, coloured by T2D vs Control logFC
    lead <- ekg_dt[p.adjust < 0.05][order(p.adjust)]
    lead <- lead[!grepl("disease|infection|cancer|carcinoma|virus", Description, ignore.case = TRUE)][1]
    if (nrow(lead) && requireNamespace("pathview", quietly = TRUE)) {
      fc <- deg[, .(gene, logFC)]
      fc[, entrez := suppressMessages(mapIds(org.Hs.eg.db, gene, "ENTREZID", "SYMBOL"))]
      fc <- fc[!is.na(entrez)]
      fc_vec <- setNames(fc$logFC, fc$entrez)
      pv_dir <- normalizePath(file.path(pv_main, "pathview"), mustWork = FALSE)
      dir.create(pv_dir, recursive = TRUE, showWarnings = FALSE)
      old <- getwd()
      # pathview writes into the working directory; always restore it, even if pathview fails
      tryCatch({
        setwd(pv_dir)
        pathview(gene.data = fc_vec, pathway.id = sub("^hsa", "", lead$ID), species = "hsa",
                 out.suffix = "T2D_vs_Control_logFC", limit = list(gene = 1.5), kegg.native = TRUE)
      }, error = function(e) message("pathview failed: ", conditionMessage(e)),
         finally = setwd(old))
      stopifnot(identical(normalizePath(getwd()), normalizePath(old)))
      message("pathview drawn for ", lead$ID, " ", lead$Description)
    }
  }
}

if (nrow(ego_sig)) {
  ego_bp <- enrichGO(genes_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                     universe = universe_entrez, pvalueCutoff = 0.05, readable = TRUE)
  if (nrow(as.data.frame(ego_bp))) {
    pC2 <- emapplot(pairwise_termsim(ego_bp), showCategory = 30) + ggtitle("GO BP term network")
    ggsave(file.path(fig_main, "Fig2C_emap_GOBP.pdf"), pC2, width = 10, height = 9)
  }
}

# ---------------- Context: all T2D vs Control DEGs up / down -------------------------------------
sig_col <- "deg_call"
deg_lists <- list(Up = to_entrez(deg[get(sig_col) == "Up", gene]), Down = to_entrez(deg[get(sig_col) == "Down", gene]))
deg_lists <- deg_lists[lengths(deg_lists) >= 5]
cc <- if (!length(deg_lists)) NULL else compareCluster(deg_lists,
                     fun = "enrichGO", OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                     universe = universe_entrez, pvalueCutoff = 0.05, readable = TRUE)
cc_dt <- if (is.null(cc)) data.table() else as.data.table(as.data.frame(cc))
if (nrow(cc_dt)) fwrite(cc_dt, file.path(out_dir, "enrich_DEG_up_down_GO.csv")) else
  message("GO BP on up/down DEGs: no enrichment at adj.P < 0.05")
if (nrow(cc_dt)) {
  ggsave(file.path(fig_dir, "05_GO_DEG_up_down.pdf"), dotplot(cc, showCategory = 12), width = 8, height = 9)
}

# ---------------- GSEA (fgsea) on moderated t ----------------------------------------
get_msig <- function(collection, subcollection = NULL) {
  f <- names(formals(msigdbr::msigdbr))
  args <- list(species = "Homo sapiens")
  if ("collection" %in% f) { args$collection <- collection; if (!is.null(subcollection)) args$subcollection <- subcollection }
  else { args$category <- collection; if (!is.null(subcollection)) args$subcategory <- subcollection }
  as.data.table(do.call(msigdbr, args))
}
collections <- list(
  Hallmark = get_msig("H"),
  Reactome = get_msig("C2", "CP:REACTOME"),
  KEGG     = tryCatch(get_msig("C2", "CP:KEGG_MEDICUS"), error = function(e) get_msig("C2", "CP:KEGG_LEGACY"))
)
gsea_plots <- list()
for (cn in c("T2DvsControl")) {
  tab   <- deg
  ranks <- sort(setNames(tab$t, tab$gene), decreasing = TRUE)
  for (nm in names(collections)) {
    pw  <- split(collections[[nm]]$gene_symbol, collections[[nm]]$gs_name)
    res <- fgsea(pw, ranks, minSize = 15, maxSize = 500, eps = 0)
    res[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
    setorder(res, padj)
    fwrite(res, file.path(out_dir, paste0("GSEA_", cn, "_", nm, ".csv")))
    message(cn, " ", nm, ": ", res[padj < 0.05, .N], " sets padj<0.05")
    if (nm == "Hallmark") {
      top <- res[padj < 0.05][order(NES)]
      if (nrow(top)) {
        top[, pathway := factor(sub("^HALLMARK_", "", pathway), levels = sub("^HALLMARK_", "", pathway))]
        gsea_plots[[cn]] <- ggplot(top, aes(NES, pathway, fill = NES > 0)) + geom_col() +
          scale_fill_manual(values = c(`TRUE` = "#C8322F", `FALSE` = "#3B7DD8"), guide = "none") +
          labs(title = paste0(cn, " Hallmark GSEA (padj < 0.05)"), y = NULL) + theme_bw(base_size = 10)
      }
    }
  }
}
if (length(gsea_plots)) {
  ggsave(file.path(fig_main, "Fig2_GSEA_hallmark.pdf"), wrap_plots(gsea_plots, nrow = 1), width = 13, height = 7)
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05.txt"))