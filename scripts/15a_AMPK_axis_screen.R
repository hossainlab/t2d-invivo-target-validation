# 15a_AMPK_axis_screen.R — screen the selected pathway (AMPK signalling, hsa04152; decision R17) in mouse
# STZ vs Control scRNA, to choose the key cell type and gene pair for the Fig 3D-G analogue.
#
# Reference paper (Xu et al. 2026) Fig 3D-E: cell-level Wilcoxon (BH) of the axis genes in the key cell type.
# Here both are reported: cell-level Wilcoxon (paper style) and per-mouse pseudobulk edgeR (script 09; primary,
# because the paper had n = 1 mouse per group and we have 3 vs 4).
#
# Gene universe: all human genes on the KEGG hsa04152 map (pathview KGML from script 05, in
#                results/pathway_selection/pathview/) -> 1:1 mouse orthologs.
# "Human-DE" genes = bulk T2D vs lean P < 0.05 (script 02); direction taken from the human logFC.
#
# Outputs: results/pathway_selection/AMPK_screen_<tissue>.csv (gene x cell type)
#          results/pathway_selection/AMPK_screen_celltype_summary.csv
#          results/pathway_selection/AMPK_pathway_genes.csv (human/mouse map + human bulk stats)

suppressPackageStartupMessages({
  library(Seurat)
  library(babelgene)
  library(org.Hs.eg.db)
  library(data.table)
  library(xml2)
})
set.seed(20260914)
out_dir <- "results/pathway_selection"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. AMPK pathway genes (KEGG map) -------------------------------------------------
xml <- read_xml("results/pathway_selection/pathview/hsa04152.xml")
ids <- unique(unlist(strsplit(xml_attr(xml_find_all(xml, "//entry[@type='gene']"), "name"), " ")))
entrez <- sub("^hsa:", "", ids[grepl("^hsa:", ids)])
hum <- unique(na.omit(AnnotationDbi::mapIds(org.Hs.eg.db, entrez, "SYMBOL", "ENTREZID")))
o <- as.data.table(orthologs(genes = hum, species = "mouse", human = TRUE))
o <- unique(o[, .(human = human_symbol, mouse = symbol)])[, n_h := .N, by = human][, n_m := .N, by = mouse][n_h == 1 & n_m == 1]
bulk <- fread("results/bulk/DEG_T2D_vs_Control_all.csv")
pw <- merge(data.table(human = hum), bulk[, .(human = gene, human_logFC = logFC, human_P = P.Value, human_FDR = adj.P.Val)],
            by = "human", all.x = TRUE)
pw <- merge(pw, o[, .(human, mouse)], by = "human", all.x = TRUE)
pw[, human_DE := !is.na(human_P) & human_P < 0.05]
fwrite(pw[order(human_P)], file.path(out_dir, "AMPK_pathway_genes.csv"))
message("KEGG hsa04152: ", length(hum), " human genes; measured in bulk: ", pw[!is.na(human_P), .N],
        "; human-DE (P<0.05): ", pw[human_DE == TRUE, .N], "; with 1:1 mouse ortholog: ", pw[!is.na(mouse), .N])

# ---- 2. Screen per tissue x cell type ---------------------------------------------------
summ <- list()
for (tissue in c("liver", "kidney")) {
  so <- readRDS(file.path("results/sc", paste0(tissue, "_annotated.rds")))
  pb <- fread(file.path("results/sc", paste0(tissue, "_pseudobulk_DE_all.csv")))
  genes <- intersect(pw[!is.na(mouse), mouse], rownames(so))
  md <- so@meta.data
  res <- list()
  for (ct in sort(unique(md$cell_type))) {
    n_ctl <- sum(md$cell_type == ct & md$condition == "Control"); n_stz <- sum(md$cell_type == ct & md$condition == "STZ")
    if (n_ctl < 10 || n_stz < 10) next
    sub <- subset(so, subset = cell_type == ct)
    Idents(sub) <- "condition"
    fm <- tryCatch(FindMarkers(sub, ident.1 = "STZ", ident.2 = "Control", features = genes, test.use = "wilcox",
                               logfc.threshold = 0, min.pct = 0.05, verbose = FALSE),
                   error = function(e) NULL)
    if (is.null(fm) || !nrow(fm)) next
    fm <- as.data.table(fm, keep.rownames = "mouse")
    fm[, cell_BH := p.adjust(p_val, "BH")]
    res[[ct]] <- fm[, .(tissue, cell_type = ct, mouse, n_ctl, n_stz, cell_log2FC = avg_log2FC, pct_STZ = pct.1,
                        pct_Control = pct.2, cell_p = p_val, cell_BH)]
  }
  r <- rbindlist(res)
  r <- merge(r, pb[, .(mouse = gene, cell_type, pb_logFC = logFC, pb_P = PValue, pb_FDR = FDR)], by = c("mouse", "cell_type"), all.x = TRUE)
  r <- merge(r, pw[!is.na(mouse), .(mouse, human, human_logFC, human_P, human_DE)], by = "mouse", all.x = TRUE)
  r[, cell_concordant := human_DE & cell_BH < 0.05 & sign(cell_log2FC) == sign(human_logFC)]
  r[, pb_concordant := human_DE & !is.na(pb_P) & pb_P < 0.05 & sign(pb_logFC) == sign(human_logFC)]
  fwrite(r[order(cell_type, cell_p)], file.path(out_dir, paste0("AMPK_screen_", tissue, ".csv")))
  summ[[tissue]] <- r[, .(genes_tested = .N, human_DE_tested = sum(human_DE, na.rm = TRUE),
                          cell_concordant = sum(cell_concordant, na.rm = TRUE),
                          pb_tested = sum(!is.na(pb_P)), pb_concordant = sum(pb_concordant, na.rm = TRUE),
                          both_concordant = sum(cell_concordant & pb_concordant, na.rm = TRUE),
                          genes_both = paste(mouse[cell_concordant & pb_concordant %in% TRUE][order(cell_p[cell_concordant & pb_concordant %in% TRUE])], collapse = ";"),
                          genes_cell_only = paste(mouse[cell_concordant %in% TRUE & !(pb_concordant %in% TRUE)], collapse = ";"),
                          n_ctl = n_ctl[1], n_stz = n_stz[1]),
                      by = .(tissue, cell_type)]
}
s <- rbindlist(summ)[order(-both_concordant, -cell_concordant)]
fwrite(s, file.path(out_dir, "AMPK_screen_celltype_summary.csv"))
print(s)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_15a.txt"))
