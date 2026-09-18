# 31_pub_fig5.R - Figure 5 (Mendelian randomisation) redrawn to the publication contract.
#
# Usage: Rscript scripts/31_pub_fig5.R <stage>
#   stage : A | B | C | D | E | panels | assemble | legend | all
#
# Re-draws only; every estimate comes from what scripts 19/20 wrote. The primary analysis is the
# LD-clumped one (decision R21), so that is what panel a shows - the superseded distance-pruned pass is
# shown only in panel d, as the comparison that justifies the switch.
#   a  forest of IVW causal estimates, LD-clumped   <- mr_results_blood_ldclumped.csv
#   b  SNP-level exposure vs outcome effects        <- instruments_blood_ldclumped.csv
#   c  leave-one-out IVW                            <- mr_leaveoneout.csv
#   d  LD-clumped vs distance-pruned estimates      <- mr_blood_pruning_comparison.csv
#   e  GTEx liver instruments                       <- mr_results_tissue.csv
#
# Outputs: figures/Fig5_MR/Fig5<a-e>_*.pdf, Fig5_composite.pdf, Fig5_legend.md
#          results/figure_exports/Fig5_composite.tiff, results/figure_exports/fig5_panels/*.rds

source("scripts/figure_style.R")
suppressPackageStartupMessages({ library(ggrepel) })

args  <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) >= 1) args[1] else "all"
set.seed(20260918)

fig_dir   <- "figures/Fig5_MR"
panel_dir <- file.path(EXPORT_DIR, "fig5_panels")
mrd       <- "results/mr"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
keep <- function(p, n, w, h) save_fig(p, n, w, h, dir = fig_dir, keep_dir = panel_dir)

PRIMARY <- "IVW"
# IMPORTANT: mr_instruments.csv is the only table carrying the outcome side, and it belongs to the
# SUPERSEDED distance-pruned pass. instruments_blood_ldclumped.csv lists the SNPs that survived
# clumping but holds exposure columns only, and the two barely overlap (SIRT1: 10 vs 5 SNPs, 0 shared).
# The harmonised per-SNP data for the primary LD-clumped analysis was never written to disk, so the
# per-SNP diagnostics below CANNOT be drawn for the primary pass without re-running script 20b. They
# are therefore written to results/supplementary_figures/mr/ with distance-pruned in the file name,
# and are deliberately kept out of Fig. 5, whose panels are all the LD-clumped analysis.
supp_dir <- "results/supplementary_figures/mr"
read_instruments <- function() {
  d <- fread(file.path(mrd, "mr_instruments.csv"))
  need <- c("gene", "rsid", "exp_beta", "exp_se", "out_beta", "out_se")
  if (!all(need %in% names(d))) stop("harmonised instrument table lacks: ",
                                     paste(setdiff(need, names(d)), collapse = ", "))
  d
}

# ===================================================================================================
# a - forest of the primary IVW estimates
# ===================================================================================================
if (stage %in% c("A", "panels", "all")) {
  m <- fread(file.path(mrd, "mr_results_blood_ldclumped.csv"))[method == PRIMARY]
  m <- m[is.finite(OR) & is.finite(OR_lo) & is.finite(OR_hi)]
  m[, sig := fifelse(FDR < 0.05, "FDR < 0.05", "Not significant")]
  m[, gene := factor(gene, levels = m[order(OR), gene])]

  pA <- ggplot(m, aes(OR, gene, colour = sig)) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.25, colour = "grey55") +
    geom_errorbarh(aes(xmin = OR_lo, xmax = OR_hi), height = 0, linewidth = 0.35) +
    geom_point(size = 1.1) +
    geom_text(aes(x = Inf, label = sprintf("%d", n_snp)), hjust = 1.4, size = (BASE - 2.5) * ppt,
              family = FONT, colour = "grey35", show.legend = FALSE) +
    scale_colour_manual(values = c(`FDR < 0.05` = col_dir[["Up"]],
                                   `Not significant` = "grey55"), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.16))) +
    labs(x = "OR for T2D per s.d. of genetically predicted expression", y = NULL) +
    theme_pub() +
    theme(axis.text.y = element_text(face = "italic", size = BASE - 2),
          legend.position = "top", legend.justification = "left")
  keep(pA, "Fig5a_forest", 89, 96)
}

# ===================================================================================================
# b - SNP-level exposure vs outcome effect, for the best-instrumented genes
# ===================================================================================================
if (stage %in% c("B", "panels", "all")) {
  ins <- read_instruments()
  m   <- fread(file.path(mrd, "mr_results.csv"))[method == PRIMARY]   # distance-pruned, matches ins
  sel <- intersect(m[order(p)][n_snp >= 3, gene], ins[, .N, by = gene][N >= 3, gene])[1:4]
  sel <- sel[!is.na(sel)]
  d   <- ins[gene %in% sel]
  d[, gene := factor(gene, levels = sel)]
  sl  <- m[gene %in% sel, .(gene = factor(gene, levels = sel), b)]

  pB <- ggplot(d, aes(exp_beta, out_beta)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey75") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey75") +
    geom_errorbar(aes(ymin = out_beta - out_se, ymax = out_beta + out_se), width = 0,
                  linewidth = 0.25, colour = "grey60") +
    geom_errorbarh(aes(xmin = exp_beta - exp_se, xmax = exp_beta + exp_se), height = 0,
                   linewidth = 0.25, colour = "grey60") +
    geom_point(size = 0.7, colour = "grey20") +
    geom_abline(data = sl, aes(slope = b, intercept = 0), colour = col_dir[["Up"]],
                linewidth = 0.45) +
    facet_wrap(~ gene, scales = "free", nrow = 2) +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    scale_y_continuous(breaks = scales::breaks_pretty(3)) +
    labs(x = "SNP effect on expression", y = "SNP effect on T2D") +
    theme_pub() +
    theme(strip.text = element_text(face = "italic", size = BASE - 1),
          panel.spacing.x = unit(3.2, "mm"), panel.spacing.y = unit(2.2, "mm"))
  save_fig(pB, "MR_snp_effects_distance_pruned_SUPERSEDED", 120, 60, dir = supp_dir)
}

# ===================================================================================================
# c - leave-one-out
# ===================================================================================================
if (stage %in% c("C", "panels", "all")) {
  lo <- fread(file.path(mrd, "mr_leaveoneout.csv"))
  m  <- fread(file.path(mrd, "mr_results.csv"))[method == PRIMARY]   # distance-pruned, matches lo
  sel <- lo[, .N, by = gene][order(-N)][1:3, gene]
  d  <- lo[gene %in% sel]
  d[, gene := factor(gene, levels = sel)]
  d[, lo := b - 1.96 * se][, hi := b + 1.96 * se]
  d[, dropped := factor(dropped, levels = unique(dropped[order(b)]))]
  full <- m[gene %in% d$gene, .(gene = factor(gene, levels = levels(d$gene)), b)]

  pC <- ggplot(d, aes(b, dropped)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.25, colour = "grey55") +
    geom_vline(data = full, aes(xintercept = b), colour = col_dir[["Up"]], linewidth = 0.35) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.3, colour = "grey45") +
    geom_point(size = 0.7) +
    facet_wrap(~ gene, scales = "free", nrow = 1) +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    labs(x = "IVW estimate with the SNP removed", y = NULL) +
    theme_pub() +
    theme(strip.text = element_text(face = "italic", size = BASE - 1),
          axis.text.y = element_text(size = BASE - 3),
          panel.spacing.x = unit(3.2, "mm"))
  save_fig(pC, "MR_leave_one_out_distance_pruned_SUPERSEDED", 120, 54, dir = supp_dir)
}

# ===================================================================================================
# d - LD clumping vs distance pruning: why the first pass over-called
# ===================================================================================================
if (stage %in% c("D", "panels", "all")) {
  cp <- fread(file.path(mrd, "mr_blood_pruning_comparison.csv"))
  cp[, `:=`(or_ld = exp(b_ld), or_dist = exp(b_dist))]
  cp[, status := fifelse(FDR_dist < 0.05 & FDR_ld >= 0.05, "Lost after LD clumping",
                 fifelse(FDR_ld < 0.05, "Significant in both", "Not significant"))]
  lbl <- cp[status != "Not significant"]

  pD <- ggplot(cp, aes(or_dist, or_ld, colour = status)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.25, colour = "grey55") +
    geom_hline(yintercept = 1, linewidth = 0.2, colour = "grey80") +
    geom_vline(xintercept = 1, linewidth = 0.2, colour = "grey80") +
    geom_point(aes(size = n_snp_ld)) +
    geom_text_repel(data = lbl, aes(label = gene), size = (BASE - 2.5) * ppt, family = FONT,
                    fontface = "italic", show.legend = FALSE, min.segment.length = 0.2,
                    segment.size = 0.2, box.padding = 0.25, max.overlaps = Inf) +
    scale_colour_manual(values = c(`Lost after LD clumping` = col_dir[["Up"]],
                                   `Significant in both` = pal_qual[1],
                                   `Not significant` = "grey65"), name = NULL) +
    scale_size_continuous(range = c(0.5, 2.2), name = "SNPs", breaks = scales::breaks_pretty(3)) +
    labs(x = "OR, distance-pruned instruments", y = "OR, LD-clumped instruments") +
    theme_pub() +
    theme(legend.position = "right", legend.box = "vertical",
          legend.key.size = unit(2.2, "mm"))
  keep(pD, "Fig5b_clumping_comparison", 90, 62)
}

# ===================================================================================================
# e - tissue-specific instruments (GTEx liver)
# ===================================================================================================
if (stage %in% c("E", "panels", "all")) {
  ti <- fread(file.path(mrd, "mr_results_tissue.csv"))
  ti[, gene := factor(gene, levels = ti[order(OR), gene])]

  pE <- ggplot(ti, aes(OR, gene)) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.25, colour = "grey55") +
    geom_errorbarh(aes(xmin = OR_lo, xmax = OR_hi), height = 0, linewidth = 0.35,
                   colour = "grey45") +
    geom_point(size = 1.1, colour = pal_qual[1]) +
    geom_text(aes(x = Inf, label = sprintf("%d", n_snp)), hjust = 1.6, size = (BASE - 2.5) * ppt,
              family = FONT, colour = "grey35") +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.18))) +
    labs(x = "OR for T2D (GTEx liver instruments)", y = NULL) +
    theme_pub() +
    theme(axis.text.y = element_text(face = "italic"))
  keep(pE, "Fig5c_tissue_forest", 90, 40)
}


# ===================================================================================================
# f - Steiger directionality (reference framework Fig. 4D; previously tabular only)
# ===================================================================================================
if (stage %in% c("F", "panels", "all")) {
  st <- fread(file.path(mrd, "mr_steiger.csv"))
  d <- melt(st[, .(gene, Exposure = r2_exposure, Outcome = r2_outcome)],
            id.vars = "gene", variable.name = "trait", value.name = "r2")
  d[, gene := factor(gene, levels = st[order(r2_exposure), gene])]
  wrong <- st[correct_direction == FALSE, gene]

  pF <- ggplot(d, aes(r2, gene)) +
    geom_line(aes(group = gene), colour = "grey80", linewidth = 0.3) +
    geom_point(aes(colour = trait), size = 1) +
    scale_x_log10(labels = scales::label_log()) +
    scale_colour_manual(values = c(Exposure = pal_qual[1], Outcome = "grey55"), name = NULL) +
    labs(x = expression("Variance explained ("*italic(r)^2*", log scale)"), y = NULL) +
    theme_pub() +
    theme(axis.text.y = element_text(face = "italic", size = BASE - 2.5),
          legend.position = "top", legend.justification = "left")
  keep(pF, "Fig5d_steiger", 89, 96)
  if (length(wrong)) message("Steiger flags wrong direction for: ", paste(wrong, collapse = ", "))
}

# ===================================================================================================
# assemble
# ===================================================================================================
if (stage %in% c("assemble", "all")) {
  need <- c("Fig5a_forest", "Fig5b_clumping_comparison", "Fig5c_tissue_forest", "Fig5d_steiger")
  ps <- lapply(need, function(n) {
    f <- file.path(panel_dir, paste0(n, ".rds"))
    if (!file.exists(f)) stop("missing panel ", n, "; run stage panels first")
    readRDS(f)$plot
  })
  design <- "
AAAAAADDDDDD
AAAAAADDDDDD
AAAAAADDDDDD
BBBBBBDDDDDD
BBBBBBCCCCCC
"
  comp <- ps[[1]] + ps[[2]] + ps[[3]] + ps[[4]] +
    plot_layout(design = design) +
    plot_annotation(tag_levels = list(tag_letters(4)[c(1, 2, 3, 4)])) &
    theme(plot.tag = element_text(size = BASE + 1, face = "bold", family = FONT))
  save_fig(comp, "Fig5_composite", W_2COL, 165, dir = fig_dir, tiff = TRUE)
}

# ===================================================================================================
# legend
# ===================================================================================================
if (stage %in% c("legend", "panels", "all")) {
  m  <- fread(file.path(mrd, "mr_results_blood_ldclumped.csv"))[method == PRIMARY]
  cp <- fread(file.path(mrd, "mr_blood_pruning_comparison.csv"))
  ti <- fread(file.path(mrd, "mr_results_tissue.csv"))
  hit <- m[FDR < 0.05][order(FDR)]
  lost <- cp[FDR_dist < 0.05 & FDR_ld >= 0.05, gene]

  l <- c(
    "# Figure 5 legend - Mendelian randomisation of candidate gene expression on T2D",
    "",
    "**Fig. 5 | Cis-eQTL Mendelian randomisation does not support a causal role for most candidates.**",
    sprintf("**a**, Inverse-variance-weighted causal estimates for the %d instrumentable candidate genes, as odds ratio for T2D per standard deviation of genetically predicted whole-blood expression (eQTLGen exposures, GCST006867 outcome). Bars, 95%% confidence interval; the number at the right of each row is the instrument count. Coloured, FDR < 0.05.",
            nrow(m)),
    "**b**, SNP-level effect on expression against effect on T2D for the four best-instrumented genes; line, the IVW slope. Error bars, standard errors.",
    "**c**, Leave-one-out IVW estimates; the vertical line is the estimate using all SNPs.",
    "**d**, Estimates from LD-clumped instruments against the superseded distance-pruned pass, per gene; point size, instrument count after clumping.",
    sprintf("**c**, The %d genes with usable GTEx liver instruments (Wald ratio, single SNP each); kidney cortex yielded none.",
            nrow(ti)),
    sprintf("**d**, Steiger directionality test: variance in expression explained by the instruments (exposure) against variance explained in T2D (outcome), per gene. Exposure exceeds outcome for %d of %d genes, i.e. the instruments act on expression first.",
            fread(file.path(mrd, "mr_steiger.csv"))[correct_direction == TRUE, .N],
            nrow(fread(file.path(mrd, "mr_steiger.csv")))),
    "",
    "## Caveats that belong in the text, not the figure",
    "",
    sprintf("- Panel a is the LD-clumped analysis (decision R21). The first pass used distance pruning and over-called: %s %s significant there and not after clumping, which is what panel d shows.",
            paste(lost, collapse = ", "), if (length(lost) == 1) "was" else "were"),
    sprintf("- %d of %d genes reach FDR < 0.05, and %s.",
            nrow(hit), nrow(m),
            if (nrow(hit)) sprintf("the strongest is %s (OR %s, 95%% CI %s-%s, %s)", hit[1, gene],
                                   mfmt(hit[1, OR]), mfmt(hit[1, OR_lo]), mfmt(hit[1, OR_hi]),
                                   p_short(hit[1, p])) else "none survives"),
    "- Exposures are whole-blood eQTLs, not liver; **e** is the liver-specific check and is limited to single-SNP Wald ratios, which cannot be tested for pleiotropy.",
    "- Instrument counts are small for most genes, so horizontal pleiotropy tests (Q, Egger intercept) have little power; they are tabulated in `results/mr/mr_results_blood_ldclumped.csv` rather than plotted.",
    "",
    "## Statistics to quote",
    "",
    "| Quantity | Value |",
    "|---|---|",
    sprintf("| Genes with instruments | %d |", nrow(m)),
    sprintf("| Genes at FDR < 0.05 (LD-clumped) | %d |", nrow(hit)),
    if (nrow(hit)) sprintf("| Strongest estimate | %s, OR %s (95%% CI %s-%s), %s, FDR %s |",
                           hit[1, gene], mfmt(hit[1, OR]), mfmt(hit[1, OR_lo]), mfmt(hit[1, OR_hi]),
                           p_short(hit[1, p]), mfmt(hit[1, FDR], 3)),
    sprintf("| Lost after LD clumping | %s |",
            if (length(lost)) paste(lost, collapse = ", ") else "none"),
    sprintf("| Median instrument count | %s |", mfmt(median(m$n_snp), 0)),
    sprintf("| GTEx liver genes tested | %s |", paste(ti$gene, collapse = ", "))
  )
  writeLines(Filter(Negate(is.null), l), file.path(fig_dir, "Fig5_legend.md"))
  message("wrote ", file.path(fig_dir, "Fig5_legend.md"))
}
