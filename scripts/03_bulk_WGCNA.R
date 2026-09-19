# 03_bulk_WGCNA.R — Phase 3 of docs/analysis_plan.md
# Signed-hybrid WGCNA on ComBat-corrected merged liver arrays.
# Traits: T2D (T2D vs lean Control, 27 samples, D1b), dataset (confound check), HbA1c (GSE15653 only).
# Obesity/group traits are not used: with lean controls only they are identical to T2D.
# With < 30 samples the top 5000 MAD genes are used (plan Phase 3).
# Key modules (decision M11b): modules correlated with T2D at p < 0.1 that also hold in BOTH cohorts
# separately (same sign, p < 0.1 each). The superseded M11 rule also required abs(r_T2D) > abs(r_Dataset),
# which is vacuous because datExpr is ComBat-corrected; see the module-trait section.
# Hub genes: key-module genes with |kME| > 0.7 and |GS| > 0.2.
# Robustness: modulePreservation GSE15653 (ref) -> GSE64998 (test).
#
# Outputs: results/bulk/WGCNA_*.csv, WGCNA_net.rds
# Figures: figures/bulk/Fig1C_soft_threshold.pdf, Fig1D_dendrogram.pdf,
#          Fig1E_module_trait.pdf, 03_module_preservation.pdf; Fig1F_MM_GS_<module>.pdf for all key modules go to
#          results/supplementary_figures/bulk (04 copies panels of candidate modules into figures/)

suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
  library(ggplot2)
})
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 8)
set.seed(20260914)

REQUIRE_REPLICATION <- TRUE   # decision M11b; FALSE reproduces the superseded M11 module set
P_COHORT            <- 0.1    # per-cohort significance a module must reach in BOTH cohorts
out_dir <- "results/bulk"; fig_dir <- "results/supplementary_figures/bulk"
fig_main <- "results/supplementary_figures/bulk"; exp_dir <- "results/figure_exports"  # draft panels;
# figures/Fig1 is owned by scripts/28_pub_fig1.R
for (d in c(fig_dir, fig_main, exp_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

expr <- readRDS(file.path(out_dir, "expr_merged_combat.rds"))
meta <- fread(file.path(out_dir, "meta_bulk.csv"))
meta <- meta[match(colnames(expr), GSM)]
stopifnot(identical(meta$GSM, colnames(expr)))

# ---------------- Input genes ---------------------------------------------------------
datExpr <- t(expr)  # samples x genes
if (nrow(datExpr) < 30) {
  mad_g   <- apply(expr, 1, mad)
  datExpr <- datExpr[, names(sort(mad_g, decreasing = TRUE))[1:5000]]
  message("n < 30: using top 5000 MAD genes")
}
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
message("WGCNA input: ", nrow(datExpr), " samples x ", ncol(datExpr), " genes")

# ---------------- Soft threshold --------------------------------------------------------
powers <- c(1:10, seq(12, 30, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed hybrid",
                         corFnc = "bicor", corOptions = list(maxPOutliers = 0.1), verbose = 0)
fi <- data.table(sft$fitIndices)
fi[, R2 := -sign(slope) * SFT.R.sq]
fwrite(fi, file.path(out_dir, "WGCNA_soft_threshold.csv"))
power <- fi[R2 >= 0.80, min(Power)]
if (!is.finite(power)) {
  power <- if (nrow(datExpr) < 20) 9 else if (nrow(datExpr) < 30) 8 else if (nrow(datExpr) < 40) 7 else 6  # Langfelder FAQ (unsigned/signed-hybrid)
  message("Scale-free R2 >= 0.80 not reached; using guideline power = ", power)
}
message("Soft power: ", power)

pdf(file.path(fig_main, "Fig1C_soft_threshold.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plot(fi$Power, fi$R2, type = "n", xlab = "Soft threshold (power)",
     ylab = "Scale-free topology fit (signed R²)", main = "Scale independence")
text(fi$Power, fi$R2, labels = fi$Power, col = "red"); abline(h = 0.80, col = "red", lty = 2)
plot(fi$Power, fi$mean.k., type = "n", xlab = "Soft threshold (power)", ylab = "Mean connectivity",
     main = "Mean connectivity")
text(fi$Power, fi$mean.k., labels = fi$Power, col = "red")
dev.off()

# ---------------- Network -------------------------------------------------------------
net <- blockwiseModules(datExpr, power = power, networkType = "signed hybrid", TOMType = "signed",
                        corType = "bicor", maxPOutliers = 0.1,
                        minModuleSize = 30, mergeCutHeight = 0.25, deepSplit = 2,
                        pamRespectsDendro = FALSE, maxBlockSize = 20000,
                        numericLabels = TRUE, saveTOMs = FALSE, verbose = 2)

moduleColors <- labels2colors(net$colors)
saveRDS(list(net = net, moduleColors = moduleColors, power = power, genes = colnames(datExpr)),
        file.path(out_dir, "WGCNA_net.rds"))
print(table(moduleColors))

pdf(file.path(fig_main, "Fig1D_dendrogram.pdf"), width = 11, height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], "Module",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# ---------------- Module-trait ---------------------------------------------------------
MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
MEs <- MEs[, colnames(MEs) != "MEgrey", drop = FALSE]
traits <- data.frame(
  T2D = meta$T2D,
  Dataset = as.integer(meta$dataset == "GSE64998"),
  row.names = meta$GSM
)
mt_cor <- cor(MEs, traits, use = "p")
mt_p   <- corPvalueStudent(mt_cor, nrow(datExpr))

# HbA1c within GSE15653
s15 <- meta[dataset == "GSE15653" & !is.na(HbA1c), GSM]
hb_cor <- cor(MEs[s15, ], meta[match(s15, GSM), HbA1c], use = "p")
hb_p   <- corPvalueStudent(hb_cor, length(s15))
mt_cor <- cbind(mt_cor, HbA1c_GSE15653 = hb_cor[, 1])
mt_p   <- cbind(mt_p,   HbA1c_GSE15653 = hb_p[, 1])

# ---- per-cohort T2D correlation (decision M11b) ---------------------------------------
# datExpr is ComBat-corrected, so the Dataset column above is bounded near zero by construction
# (max |r| = 0.08 here, against |r| up to 1.00 on the uncorrected matrix). It therefore cannot act as
# a confound filter: the old rule abs(r_T2D) > abs(r_Dataset) excluded 0 of 14 modules. It is kept as
# a reported column, but the filter is now whether a module's T2D association holds in BOTH cohorts
# separately, which is what a residual batch or cohort effect would break.
i1 <- meta$dataset == "GSE15653"; i2 <- meta$dataset == "GSE64998"
r1 <- cor(MEs[i1, ], meta$T2D[i1], use = "p"); p1 <- corPvalueStudent(r1, sum(i1))
r2 <- cor(MEs[i2, ], meta$T2D[i2], use = "p"); p2 <- corPvalueStudent(r2, sum(i2))
mt_cor <- cbind(mt_cor, T2D_GSE15653 = r1[, 1], T2D_GSE64998 = r2[, 1])
mt_p   <- cbind(mt_p,   T2D_GSE15653 = p1[, 1], T2D_GSE64998 = p2[, 1])
message("per-cohort n: GSE15653 = ", sum(i1), ", GSE64998 = ", sum(i2))

mt_tab <- data.table(module = sub("^ME", "", rownames(mt_cor)),
                     size = as.integer(table(moduleColors)[sub("^ME", "", rownames(mt_cor))]),
                     setnames(data.table(mt_cor), paste0("r_", colnames(mt_cor))),
                     setnames(data.table(mt_p), paste0("p_", colnames(mt_p))))
mt_tab[, concordant := sign(r_T2D_GSE15653) == sign(r_T2D_GSE64998)]
mt_tab[, replicates := concordant & p_T2D_GSE15653 < P_COHORT & p_T2D_GSE64998 < P_COHORT]
fwrite(mt_tab, file.path(out_dir, "WGCNA_module_trait.csv"))

pdf(file.path(fig_main, "Fig1E_module_trait.pdf"), width = 9, height = max(5, 0.35 * nrow(mt_cor) + 2))
par(mar = c(7, 9, 3, 2))
labeledHeatmap(Matrix = mt_cor, xLabels = colnames(mt_cor), yLabels = rownames(mt_cor),
               ySymbols = rownames(mt_cor), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = paste0(signif(mt_cor, 2), "\n(", signif(mt_p, 1), ")"),
               setStdMargins = FALSE, cex.text = 0.6, zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()

# ---------------- Key modules ----------------------------------------------------------
# M11  (superseded): p_T2D < 0.1 and abs(r_T2D) > abs(r_Dataset). The second clause is vacuous on
#                    ComBat-corrected input, so this reduced to p_T2D < 0.1 alone.
# M11b (current):    p_T2D < 0.1 AND the same direction in both cohorts AND p < P_COHORT in each.
#                    Set REQUIRE_REPLICATION <- FALSE to reproduce the M11 module set.
cand_all <- mt_tab[module != "grey" & p_T2D < 0.1]
cand <- if (REQUIRE_REPLICATION) cand_all[replicates == TRUE] else cand_all
key_modules <- cand[order(p_T2D), module]
dropped <- setdiff(cand_all[order(p_T2D), module], key_modules)
message("Key modules (M11b, T2D p<0.1 + replicated): ", paste(key_modules, collapse = ", "))
if (length(dropped))
  message("Dropped for failing replication: ",
          paste(sprintf("%s (r %.2f / %.2f, p %.3g / %.3g)", dropped,
                        mt_tab[match(dropped, module), r_T2D_GSE15653],
                        mt_tab[match(dropped, module), r_T2D_GSE64998],
                        mt_tab[match(dropped, module), p_T2D_GSE15653],
                        mt_tab[match(dropped, module), p_T2D_GSE64998]), collapse = "; "))
fwrite(mt_tab[, .(module, size, r_T2D, p_T2D, r_T2D_GSE15653, p_T2D_GSE15653,
                  r_T2D_GSE64998, p_T2D_GSE64998, concordant, replicates,
                  key = module %in% key_modules)][order(p_T2D)],
       file.path(out_dir, "WGCNA_key_module_selection.csv"))

# ---------------- GS / MM -------------------------------------------------------------
GS   <- as.numeric(cor(datExpr, traits$T2D, use = "p"))
GS_p <- corPvalueStudent(GS, nrow(datExpr))
MM   <- cor(datExpr, MEs, use = "p")
gene_tab <- data.table(gene = colnames(datExpr), module = moduleColors, GS_T2D = GS, GS_p = GS_p)
for (m in key_modules) {
  gene_tab[, paste0("MM_", m) := MM[, paste0("ME", m)]]
}
gene_tab[, kME_own := MM[cbind(seq_len(.N), match(paste0("ME", module), colnames(MM)))]]
gene_tab[, hub := module %in% key_modules & abs(kME_own) > 0.7 & abs(GS_T2D) > 0.2]  # hub rule of the earlier canonical analysis (M10)
fwrite(gene_tab, file.path(out_dir, "WGCNA_gene_module_GS_MM.csv"))

for (m in key_modules) {
  g  <- gene_tab[module == m]
  ct <- cor.test(abs(g$kME_own), abs(g$GS_T2D))
  p  <- ggplot(g, aes(abs(kME_own), abs(GS_T2D))) +
    geom_point(colour = ifelse(m %in% c("white", "lightyellow", "ivory"), "grey40", m), alpha = 0.7) +
    geom_smooth(method = "lm", colour = "black", se = FALSE, linewidth = 0.6) +
    geom_vline(xintercept = 0.7, linetype = 2) + geom_hline(yintercept = 0.2, linetype = 2) +
    labs(title = paste0(m, " module (n = ", nrow(g), ")"),
         subtitle = sprintf("cor = %.2f, p = %.1e", ct$estimate, ct$p.value),
         x = "Module membership |MM|", y = "Gene significance |GS| for T2D") +
    theme_bw(base_size = 12)
  ggsave(file.path(fig_dir, paste0("Fig1F_MM_GS_", m, ".pdf")), p, width = 5.5, height = 5)
}

# ---------------- Module preservation GSE15653 -> GSE64998 ---------------------------
raw <- readRDS(file.path(out_dir, "expr_merged_raw.rds"))[colnames(datExpr), ]
s1 <- meta[dataset == "GSE15653", GSM]; s2 <- meta[dataset == "GSE64998", GSM]
multiExpr  <- list(GSE15653 = list(data = t(raw[, s1])), GSE64998 = list(data = t(raw[, s2])))
multiColor <- list(GSE15653 = moduleColors)
mp <- modulePreservation(multiExpr, multiColor, referenceNetworks = 1, nPermutations = 100,
                         networkType = "signed hybrid", corFnc = "bicor",
                         randomSeed = 20260914, quickCor = 0, verbose = 1)
ps <- mp$preservation$Z$ref.GSE15653$inColumnsAlsoPresentIn.GSE64998
pres <- data.table(module = rownames(ps), size = ps$moduleSize, Zsummary = ps$Zsummary.pres,
                   medianRank = mp$preservation$observed$ref.GSE15653$inColumnsAlsoPresentIn.GSE64998$medianRank.pres)
pres <- pres[!module %in% c("gold", "grey")]
fwrite(pres, file.path(out_dir, "WGCNA_module_preservation.csv"))
p_pres <- ggplot(pres, aes(size, Zsummary, label = module)) +
  geom_point(aes(colour = module), size = 3, show.legend = FALSE) +
  geom_text(vjust = -0.8, size = 3) + scale_colour_identity() +
  geom_hline(yintercept = c(2, 10), linetype = 2, colour = c("blue", "darkgreen")) +
  scale_x_log10() + labs(x = "Module size", y = "Zsummary (GSE64998 vs GSE15653)") + theme_bw()
ggsave(file.path(fig_dir, "03_module_preservation.pdf"), p_pres, width = 6, height = 5)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_03.txt"))
