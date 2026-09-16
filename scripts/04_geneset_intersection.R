# 04_geneset_intersection.R — Phase 4 of docs/analysis_plan.md
# Build hyperglycemia / insulin-resistance phenotype gene set (analogue of the reference's
# GeneCards neutrophil-migration set) and intersect: DEG(T2D vs Control) ∩ WGCNA key modules ∩ gene set.
#
# Gene-set sources (each recorded per gene):
#   - MSigDB GO:BP: response to glucose, cellular response to glucose stimulus, response to insulin
#   - KEGG (REST): hsa04933 AGE-RAGE in diabetic complications, hsa04931 insulin resistance,
#                  hsa04930 type II diabetes mellitus
#   - HPO: HP:0003074 Hyperglycemia, HP:0000855 Insulin resistance
#   - CTD curated (DirectEvidence): Hyperglycemia D006943, Insulin Resistance D007333
#   - GeneCards (optional manual export at data/genesets/GeneCards_hyperglycemia.csv, relevance > 8)
#
# Outputs: results/bulk/geneset_hyperglycemia.csv, intersect_genes.csv, intersection_summary.csv
# Figure:  figures/Fig1_bulk_DEG_WGCNA/Fig1G_venn.pdf (TIFF in results/figure_exports)

suppressPackageStartupMessages({
  library(data.table)
  library(msigdbr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggvenn)
  library(ggplot2)
})
out_dir <- "results/bulk"; fig_dir <- "results/supplementary_figures/bulk"; gs_dir <- "data/genesets"
fig_main <- "figures/Fig1_bulk_DEG_WGCNA"; exp_dir <- "results/figure_exports"  # framework panels / TIFF exports
for (d in c(fig_dir, fig_main, exp_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
dir.create(gs_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------- MSigDB GO:BP ----------------------------------------------------------
get_msig <- function(collection, subcollection) {
  f <- names(formals(msigdbr::msigdbr))
  if ("collection" %in% f) msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection)
  else msigdbr(species = "Homo sapiens", category = collection, subcategory = subcollection)
}

gobp <- as.data.table(get_msig("C5", "GO:BP"))
# Select by GO ID. In the installed MSigDB release, "response to glucose" (GO:0009749) and
# "cellular response to glucose stimulus" (GO:0071333) are not separate sets; their parent terms
# are used: response to monosaccharide (GO:0034284) and cellular response to carbohydrate
# stimulus (GO:0071322).
go_ids <- c("GO:0034284", "GO:0071322", "GO:0032868")
go_sets <- gobp[gs_exact_source %in% go_ids, .(gene = unique(gene_symbol)), by = .(source = gs_name)]
missing_go <- setdiff(go_ids, gobp[gs_exact_source %in% go_ids, unique(gs_exact_source)])
if (length(missing_go)) warning("GO IDs not found in this msigdbr release: ", paste(missing_go, collapse = ", "))
message("GO sets used: ", paste(unique(go_sets$source), collapse = ", "))
stopifnot(nrow(go_sets) > 0)

# ---------------- KEGG via REST ---------------------------------------------------------
kegg_ids <- c(hsa04933 = "KEGG_AGE_RAGE_diabetic_complications",
              hsa04931 = "KEGG_insulin_resistance",
              hsa04930 = "KEGG_type_II_diabetes")
kegg_cache <- file.path(gs_dir, "kegg_diabetes_pathways.csv")
if (!file.exists(kegg_cache)) {
  kegg <- rbindlist(lapply(names(kegg_ids), function(pid) {
    x <- readLines(paste0("https://rest.kegg.jp/link/hsa/", pid), warn = FALSE)
    entrez <- sub("^hsa:", "", sub("^.*\t", "", x))
    data.table(source = kegg_ids[[pid]], entrez = entrez)
  }))
  kegg[, gene := mapIds(org.Hs.eg.db, entrez, "SYMBOL", "ENTREZID")]
  fwrite(kegg, kegg_cache)
}
kegg_sets <- fread(kegg_cache, colClasses = c(entrez = "character"))[!is.na(gene), .(source, gene)]

# ---------------- HPO -------------------------------------------------------------------
hpo <- fread(file.path(gs_dir, "hpo_genes_to_phenotype.txt"))
hpo_sets <- unique(hpo[hpo_id %in% c("HP:0003074", "HP:0000855"),
                       .(source = paste0("HPO_", gsub(" ", "_", hpo_name)), gene = gene_symbol)])

# ---------------- CTD curated ----------------------------------------------------------
ctd_file <- file.path(gs_dir, "CTD_genes_diseases.tsv.gz")
ctd_sets <- data.table(source = character(), gene = character())
if (file.exists(ctd_file) && file.size(ctd_file) > 1e5) {
  # CTD_curated_genes_diseases: '#'-prefixed comment lines, 7 tab-separated columns, no header row
  con <- gzfile(ctd_file); ln <- readLines(con, warn = FALSE); close(con)
  ln  <- ln[!grepl("^#", ln) & nzchar(ln)]
  ctd <- fread(text = ln, sep = "\t", quote = "", header = FALSE, fill = TRUE,
               col.names = c("GeneSymbol", "GeneID", "DiseaseName", "DiseaseID",
                             "DirectEvidence", "OmimIDs", "PubMedIDs"))
  ctd_sets <- unique(ctd[DiseaseID %in% c("MESH:D006943", "MESH:D007333") & DirectEvidence != "",
                         .(source = paste0("CTD_", gsub(" ", "_", DiseaseName)), gene = GeneSymbol)])
} else {
  message("CTD file not available; skipping CTD source")
}

# ---------------- GeneCards (optional) --------------------------------------------------
gc_file <- file.path(gs_dir, "GeneCards_hyperglycemia.csv")
gc_sets <- data.table(source = character(), gene = character())
if (file.exists(gc_file)) {
  gc <- fread(gc_file)
  sym_col <- grep("symbol", names(gc), ignore.case = TRUE, value = TRUE)[1]
  rel_col <- grep("relevance", names(gc), ignore.case = TRUE, value = TRUE)[1]
  gc_sets <- gc[get(rel_col) > 8, .(source = "GeneCards_relevance_gt8", gene = get(sym_col))]
} else {
  message("GeneCards export not found (", gc_file, "); skipping. Export manually to add it.")
}

gs_long <- unique(rbind(go_sets, kegg_sets, hpo_sets, ctd_sets, gc_sets))
gs_long <- gs_long[!is.na(gene) & gene != ""]
gs <- gs_long[, .(n_sources = uniqueN(source), sources = paste(sort(unique(source)), collapse = ";")), by = gene]
fwrite(gs[order(-n_sources, gene)], file.path(out_dir, "geneset_hyperglycemia.csv"))
src_summary <- gs_long[, .N, by = source]
print(src_summary)
message("Hyperglycemia gene set: ", nrow(gs), " unique genes")

# ---------------- Inputs from Phases 2-3 ------------------------------------------------
deg <- fread(file.path(out_dir, "DEG_T2D_vs_Control_all.csv"))
use_fallback <- FALSE  # DEG definition fixed in 02_bulk_DEG.R (deg_call; user decision M6)
deg[, sig_use := deg_call]
deg_genes <- deg[sig_use != "NS", gene]

mt <- fread(file.path(out_dir, "WGCNA_module_trait.csv"))
cand <- mt[module != "grey" & p_T2D < 0.1 & abs(r_T2D) > abs(r_Dataset)]  # user decision M11 (same rule as 03)
key_modules <- cand[order(p_T2D), module]
wg <- fread(file.path(out_dir, "WGCNA_gene_module_GS_MM.csv"))
module_genes <- wg[module %in% key_modules, gene]

universe <- deg$gene
gs_in <- intersect(gs$gene, universe)
message("Key modules: ", paste(key_modules, collapse = ", "),
        " | DEGs: ", length(deg_genes), if (use_fallback) " (fallback)" else "",
        " | module genes: ", length(module_genes), " | gene set in universe: ", length(gs_in))

# Main candidate set (user decision M13, reference-paper rule): DEGs ∩ all genes of the T2D-associated key
# modules ∩ hyperglycemia gene set (3-way). WGCNA hub status (|kME| > 0.7 & |GS| > 0.2) is kept as annotation.
hub_genes  <- wg[hub == TRUE, gene]
triple <- Reduce(intersect, list(deg_genes, module_genes, gs_in))
double <- intersect(deg_genes, module_genes)
main_genes <- triple
main_set_type <- "DEG ∩ key-module genes ∩ hyperglycemia gene set (paper 3-way); hub status as annotation"
message("Hub genes: ", length(hub_genes), " | DEG∩hub: ", length(intersect(deg_genes, hub_genes)),
        " | DEG∩module: ", length(double), " | 3-way DEG∩module∩geneset (main set): ", length(triple))

res <- deg[gene %in% main_genes,
           .(gene, logFC, AveExpr, P.Value, adj.P.Val, sig_use,
             logFC_GSE15653, logFC_GSE64998)]
res <- merge(res, wg[, .(gene, module, GS_T2D, kME_own, hub)], by = "gene", all.x = TRUE)
res <- merge(res, gs[, .(gene, in_geneset = TRUE, geneset_sources = sources)], by = "gene", all.x = TRUE)
res[is.na(in_geneset), in_geneset := FALSE]
setorder(res, P.Value)
fwrite(res, file.path(out_dir, "intersect_genes.csv"))

# Fig 1F in figures/: MM-GS panels (written by 03 to results/supplementary_figures/bulk) for key modules
# that contain candidate genes; panels of other modules stay supplementary
unlink(Sys.glob(file.path(fig_main, "Fig1F_MM_GS_*.pdf")))
for (m in unique(na.omit(res$module))) {
  src <- file.path(fig_dir, paste0("Fig1F_MM_GS_", m, ".pdf"))
  if (file.exists(src)) file.copy(src, file.path(fig_main, basename(src)), overwrite = TRUE)
}

fwrite(data.table(
  item = c("DEG_T2D_vs_Control", "hub_genes", "key_modules", "module_genes", "geneset_total", "geneset_in_universe",
           "triple_intersection", "DEG_module_intersection", "main_set_type", "main_set_n",
           "main_set_in_geneset"),
  value = c(length(deg_genes), length(hub_genes), paste(key_modules, collapse = ";"), length(module_genes),
            nrow(gs), length(gs_in), length(triple), length(double), main_set_type, length(main_genes),
            res[in_geneset == TRUE, .N])
), file.path(out_dir, "intersection_summary.csv"))

venn <- ggvenn(list(`DEGs` = deg_genes, `WGCNA key-module genes` = module_genes, `Hyperglycemia genes` = gs_in),
               fill_color = c("#E07B39", "#3B8BC2", "#E0A100"), stroke_size = 0.5,
               show_percentage = FALSE,  # counts only
               set_name_size = 4.5, text_size = 4.5)
wrap <- function(x) paste(strwrap(x, 80), collapse = "
")
up_g <- res[logFC > 0, gene]; dn_g <- res[logFC < 0, gene]
cap <- paste(wrap(paste0("Up in T2D: ", if (length(up_g)) paste(up_g, collapse = ", ") else "none")),
             wrap(paste0("Down in T2D: ", if (length(dn_g)) paste(dn_g, collapse = ", ") else "none")),
             paste0("DEG ∩ WGCNA hub genes (annotation): ", length(intersect(deg_genes, hub_genes)), " gene(s)"), sep = "
")
venn <- venn + labs(caption = cap) + theme(plot.caption = element_text(hjust = 0.5, size = 10))
ggsave(file.path(fig_main, "Fig1G_venn.pdf"), venn, width = 6.5, height = 6, bg = "white")
unlink(file.path(exp_dir, "Fig1G_venn.tiff")); ggsave(file.path(exp_dir, "Fig1G_venn.tiff"), venn, width = 6.5, height = 6, bg = "white", dpi = 300, compression = "lzw")
