# 36_fig5_paper_style.R - Figure 5 (Mendelian randomisation) drawn in the visual style of the
# reference paper's Fig. 4.
#
# Usage: Rscript scripts/36_fig5_paper_style.R <stage>
#   stage : A | B | C | D | panels | assemble | legend | all   (default all)
#
# Reference: Xu et al., Phytomedicine 154 (2026) 158050, Fig. 4. This reproduces that figure's own
# visual language with this project's data, panel for panel:
#   A  the table-forest: Exposure | Method | nSNP | forest | OR (95% CI) | P value, with the
#      significant row in bold
#   B  SNP effect on the outcome against SNP effect on the exposure, with one fitted line per MR
#      method ("MR Test" key)
#   C  forest of the per-SNP Wald ratios, with the All - MR Egger and All - IVW summary rows below a
#      rule, drawn in red as the reference does
#   D  Steiger directionality: variance explained in the outcome against variance explained in the
#      exposure, with the y = x line dashed in red
#
# This is deliberately NOT the figure_style.R contract used by scripts 27-32. It is the reference
# paper's style, kept as a separate deliverable in figures/Fig5_paper_style/ so the Nature/Cell
# version in figures/Fig5_MR/ is untouched. Pick one for submission; do not ship both.
#
# The single most important deviation, stated at the top of the legend:
#   The reference's Fig. 4 nominates a causal druggable target (PAK1, IVW P = 0.039). This analysis
#   does NOT. Neither Tier 1 target is causal for T2D (decision R21, docs/results_MR.md), and the two
#   genes that do reach FDR < 0.05, SREBF1 and SIRT1, are listed in docs/target_selection.md as
#   "not targets" because they are supported on this one axis alone. Panel A therefore carries all
#   three exposures - the two Tier 1 targets and the strongest MR hit - rather than one, so the null
#   result is not hidden behind the one gene that happens to be significant.
#
# Every estimate comes from what scripts 19/20 wrote; nothing is recomputed. The primary analysis is
# the LD-clumped one (decision R21).
#
# Inputs:  results/mr/mr_results_blood_ldclumped.csv, instruments_blood_ldclumped_harmonised.csv,
#          mr_steiger.csv
# Outputs: figures/Fig5_paper_style/Fig5<A-D>_*.pdf, Fig5_paper_style.pdf,
#          Fig5_paper_style_legend.md, results/figure_exports/Fig5_paper_style.tiff

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
})
options(stringsAsFactors = FALSE)

args  <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) >= 1) args[1] else "all"
set.seed(20260918)

fig_dir   <- "figures/Fig5_paper_style"
panel_dir <- file.path("results/figure_exports", "fig5paper_panels")
mrd       <- "results/mr"
for (dd in c(fig_dir, panel_dir)) dir.create(dd, recursive = TRUE, showWarnings = FALSE)

TIER1   <- c("PPARGC1A", "CCND1")   # docs/target_selection.md
FEATURE <- "PPARGC1A"               # the exposure drawn per-SNP in B and C: best-instrumented Tier 1
FDR_CUT <- 0.05
METHODS <- c("MR-Egger", "Weighted median", "IVW", "Weighted mode", "Wald ratio")
col_meth <- c(`MR-Egger` = "#3B7DD8", `Weighted median` = "#2E9E5B", IVW = "#E8362B",
              `Weighted mode` = "#8E5CC8", `Wald ratio` = "#E6A700")

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
p_sci <- function(p) formatC(p, format = "E", digits = 2)

res <- fread(file.path(mrd, "mr_results_blood_ldclumped.csv"))
top <- res[method == "IVW"][order(FDR)][1, gene]          # strongest MR hit, whatever it is
shown <- unique(c(TIER1, top))

# ===================================================================================================
# A - the table-forest, as in the paper
# ===================================================================================================
if (stage %in% c("A", "panels", "all")) {
  d <- res[gene %in% shown]
  d[, gene := factor(gene, levels = shown)]
  d[, method := factor(method, levels = METHODS)]
  setorder(d, gene, method)
  d[, row := .N:1]                                        # first row at the top
  d[, bold := ifelse(p < 0.05, "bold", "plain")]

  # the forest half: OR on a log axis, as the reference draws it
  xr <- range(c(d$OR_lo, d$OR_hi, 1))
  pad <- diff(log(xr)) * 0.08
  pForest <- ggplot(d, aes(OR, row, colour = method)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey30") +
    geom_errorbar(aes(xmin = OR_lo, xmax = OR_hi), orientation = "y", width = 0,
                  linewidth = 0.5) +
    geom_point(shape = 15, size = 1.5) +
    scale_colour_manual(values = col_meth, guide = "none") +
    scale_x_continuous(trans = "log", breaks = c(0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2),
                       limits = exp(log(xr) + c(-pad, pad))) +
    # the same y range as the table half, so the two align row for row under patchwork
    scale_y_continuous(limits = c(0.4, nrow(d) + 1.9), expand = c(0, 0)) +
    labs(x = "OR (95% CI) for T2D", y = NULL) +
    theme_baseR() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.border = element_blank(), axis.line.x = element_line(linewidth = 0.4))

  # the table half: four text columns laid out on a blank canvas, as the reference does
  tab <- rbindlist(list(
    d[, .(row, col = 1, lab = as.character(gene), face = "italic", hj = 0)],
    d[, .(row, col = 2, lab = as.character(method), face = bold, hj = 0)],
    d[, .(row, col = 3, lab = as.character(n_snp), face = bold, hj = 0.5)],
    d[, .(row, col = 4, lab = sprintf("%.3f (%.3f to %.3f)", OR, OR_lo, OR_hi), face = bold, hj = 0)],
    d[, .(row, col = 5, lab = p_sci(p), face = bold, hj = 0)]))
  # x and y are resolved into the data rather than referenced from aes(). aes() is evaluated when the
  # plot is DRAWN, not when it is built, and `assemble` draws a panel saved by an earlier stage: any
  # name the mapping reaches for is looked up in whatever the session holds by then. `nrow(d) + 1`
  # inside aes() silently became stage C's `d` and threw the header row off the panel.
  xpos <- c(0.00, 0.30, 0.58, 0.70, 1.06)
  nrow_d <- nrow(d)
  tab[, x := xpos[col]]
  hdr  <- data.table(x = xpos, y = nrow_d + 1,
                     lab = c("Exposure", "Method", "nSNP", "OR (95% CI)", "P value"),
                     hj = c(0, 0, 0.5, 0, 0))
  pTab <- ggplot() +
    geom_text(data = tab, aes(x = x, y = row, label = lab, fontface = face, hjust = hj),
              size = 1.9) +
    geom_text(data = hdr, aes(x = x, y = y, label = lab, hjust = hj),
              size = 2.1, fontface = "bold") +
    geom_hline(yintercept = nrow_d + 0.55, linewidth = 0.4) +
    scale_x_continuous(limits = c(-0.02, 1.30), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.4, nrow_d + 1.9), expand = c(0, 0)) +
    theme_void()

  pA <- wrap_elements(full = pTab + pForest + plot_layout(widths = c(1.85, 1)))
  keep(pA, "Fig5A_table_forest", 190, 78)
}

# ===================================================================================================
# B - SNP effect on outcome against SNP effect on exposure, one line per method
# ===================================================================================================
if (stage %in% c("B", "panels", "all")) {
  h <- fread(file.path(mrd, "instruments_blood_ldclumped_harmonised.csv"))[gene == FEATURE]
  # orient every instrument to the expression-increasing allele, the convention these slopes assume
  h[exp_beta < 0, `:=`(exp_beta = -exp_beta, out_beta = -out_beta)]
  sl <- res[gene == FEATURE & method %in% METHODS]
  sl[, method := factor(method, levels = METHODS)]
  # MR-Egger is the only method with a free intercept; the rest pass through the origin
  eg <- res[gene == FEATURE & method == "MR-Egger"]
  ints <- sl[, .(method, slope = b,
                 intercept = ifelse(method == "MR-Egger" & !is.na(eg$egger_int), eg$egger_int, 0))]

  pB <- ggplot(h, aes(exp_beta, out_beta)) +
    geom_errorbar(aes(ymin = out_beta - out_se, ymax = out_beta + out_se),
                  width = 0, linewidth = 0.25, colour = "grey55") +
    geom_errorbar(aes(xmin = exp_beta - exp_se, xmax = exp_beta + exp_se),
                  orientation = "y", width = 0, linewidth = 0.25, colour = "grey55") +
    geom_point(size = 0.8, colour = "black") +
    geom_abline(data = ints, aes(slope = slope, intercept = intercept, colour = method),
                linewidth = 0.4) +
    scale_colour_manual(values = col_meth, name = "MR Test", drop = TRUE) +
    labs(x = sprintf("SNP effect on %s", FEATURE), y = "SNP effect on type 2 diabetes") +
    theme_baseR() +
    theme(legend.position = "inside", legend.position.inside = c(0.99, 0.99),
          legend.justification = c(1, 1),
          legend.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.2),
          legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
          legend.key.size = unit(2.6, "mm"), legend.margin = margin(1, 2, 1, 2))
  keep(pB, "Fig5B_snp_effects", 85, 78)
}

# ===================================================================================================
# C - per-SNP forest with the summary rows below a rule, as in the paper
# ===================================================================================================
if (stage %in% c("C", "panels", "all")) {
  h <- fread(file.path(mrd, "instruments_blood_ldclumped_harmonised.csv"))[gene == FEATURE]
  # a single-SNP Wald ratio and its delta-method standard error, which is what the reference plots
  h[, `:=`(b = out_beta / exp_beta, se = abs(out_se / exp_beta))]
  snp <- h[, .(lab = rsid, b, lo = b - 1.96 * se, hi = b + 1.96 * se, grp = "snp")]
  setorder(snp, b)
  sm <- res[gene == FEATURE & method %in% c("MR-Egger", "IVW")]
  sm <- sm[, .(lab = paste("All -", ifelse(method == "IVW", "Inverse variance weighted", method)),
               b, lo = b - 1.96 * se, hi = b + 1.96 * se, grp = "summary")]
  fc <- rbind(snp, sm)                    # SNPs first, so the summary rows land at the bottom
  fc[, row := .N:1]

  pC <- ggplot(fc, aes(b, row, colour = grp)) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, colour = "grey30") +
    geom_hline(yintercept = nrow(sm) + 0.5, linewidth = 0.3, colour = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                  linewidth = 0.45) +
    geom_point(size = 1.1) +
    scale_colour_manual(values = c(snp = "black", summary = "#E8362B"), guide = "none") +
    scale_y_continuous(breaks = fc$row, labels = fc$lab, expand = expansion(add = 0.7)) +
    labs(x = sprintf("MR effect size for\n'%s' on 'type 2 diabetes'", FEATURE), y = NULL) +
    theme_baseR() +
    theme(axis.text.y = element_text(size = 5))
  keep(pC, "Fig5C_snp_forest", 85, 78)
}

# ===================================================================================================
# D - Steiger directionality, as in the paper
# ===================================================================================================
if (stage %in% c("D", "panels", "all")) {
  st <- fread(file.path(mrd, "mr_steiger.csv"))
  st[, flag := ifelse(gene %in% shown, gene, NA_character_)]

  pD <- ggplot(st, aes(r2_exposure, r2_outcome)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red", linewidth = 0.4) +
    geom_point(aes(colour = correct_direction), size = 1.1) +
    geom_text_repel(aes(label = flag), size = 1.9, seed = 1, na.rm = TRUE,
                    min.segment.length = 0.2, segment.size = 0.15, segment.colour = "grey55",
                    box.padding = 0.35, max.overlaps = Inf) +
    scale_colour_manual(values = c(`TRUE` = "#3B7DD8", `FALSE` = "#E8362B"),
                        name = "Correct direction") +
    # explicit decade breaks: label_log() on the default breaks prints things like 10^-2.523
    scale_x_continuous(trans = "log10", breaks = 10^(-4:0), labels = scales::label_log()) +
    scale_y_continuous(trans = "log10", breaks = 10^(-7:-3), labels = scales::label_log()) +
    labs(x = expression("Variance explained in exposure, "*R^2),
         y = expression("Variance explained in outcome, "*R^2)) +
    theme_baseR() +
    theme(legend.position = "inside", legend.position.inside = c(0.99, 0.02),
          legend.justification = c(1, 0),
          legend.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.2),
          legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
          legend.key.size = unit(2.6, "mm"), legend.margin = margin(1, 2, 1, 2))
  keep(pD, "Fig5D_steiger", 85, 78)
}

# ===================================================================================================
# assemble - the paper's layout: A across the top, B | C | D beneath
# ===================================================================================================
if (stage %in% c("assemble", "all")) {
  need <- c("Fig5A_table_forest", "Fig5B_snp_effects", "Fig5C_snp_forest", "Fig5D_steiger")
  ps <- lapply(need, function(n) {
    f <- file.path(panel_dir, paste0(n, ".rds"))
    if (!file.exists(f)) stop("missing panel ", n, "; run stage panels first")
    readRDS(f)$plot
  })
  design <- "
AAAAAAAAA
AAAAAAAAA
BBBCCCDDD
BBBCCCDDD
"
  comp <- ps[[1]] + ps[[2]] + ps[[3]] + ps[[4]] +
    plot_layout(design = design) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 13, face = "bold"))
  out <- file.path(fig_dir, "Fig5_paper_style")
  ggsave(paste0(out, ".pdf"), comp, width = 260, height = 165, units = "mm", device = cairo_pdf)
  tif <- file.path("results/figure_exports", "Fig5_paper_style.tiff")
  unlink(tif)
  ggsave(tif, comp, width = 260, height = 165, units = "mm", dpi = 400, bg = "white",
         compression = "lzw")
  message("wrote ", out, ".pdf and ", tif)
}

# ===================================================================================================
if (stage %in% c("legend", "panels", "all")) {
  st  <- fread(file.path(mrd, "mr_steiger.csv"))
  ivw <- res[method == "IVW"]
  sig <- ivw[FDR < FDR_CUT][order(FDR)]
  ft  <- res[gene == FEATURE]
  h   <- fread(file.path(mrd, "instruments_blood_ldclumped_harmonised.csv"))[gene == FEATURE]
  ccnd <- res[gene == "CCND1"]

  l <- c(
    "# Figure 5, reference-paper style - legend",
    "",
    "Drawn to match Xu et al., *Phytomedicine* 154 (2026) 158050, Fig. 4, with this project's data.",
    "The Nature/Cell-contract version of the same figure is in `figures/Fig5_MR/`.",
    "**Ship one or the other, not both.**",
    "",
    "## Read this before the caption",
    "",
    sprintf("**The reference's Fig. 4 nominates a causal druggable target. This figure does not.** PAK1 reached IVW P = 0.039 there and was carried forward. Here neither Tier 1 target is causal for T2D: %s. The %d gene(s) that do reach FDR < %.2f (%s) are listed in `docs/target_selection.md` as *not* targets, because each is supported on this axis alone. Panel A therefore carries three exposures rather than one, so the null is on the figure rather than behind it.",
            paste(sprintf("%s OR %.3f, P = %.2f", ivw[gene %in% TIER1, gene],
                          ivw[gene %in% TIER1, OR], ivw[gene %in% TIER1, p]), collapse = "; "),
            nrow(sig), FDR_CUT, paste(sig$gene, collapse = ", ")),
    "",
    sprintf("**Fig. 5. Results of MR analysis.** (A) Forest plot illustrating the causal inference between the two Tier 1 targets (%s), the strongest MR hit (%s) and type 2 diabetes. (B) SNP-level effect on the outcome against effect on the exposure for %s, with one fitted line per MR method. (C) Forest plot for each SNP of the %s analysis. (D) Results of Steiger filtering for directionality testing.",
            paste(TIER1, collapse = ", "), top, FEATURE, FEATURE),
    "",
    "## Where this data forces a deviation from the reference",
    "",
    sprintf("- **(A) has three exposures, not one,** for the reason given above. Rows are bold where P < 0.05."),
    sprintf("- **(B) and (C) feature %s,** the better-instrumented of the two Tier 1 targets (%d SNPs against %d for CCND1). Both panels are diagnostics of pleiotropy and heterogeneity, and they are worth showing on a null estimate: the reference's versions carry the same information about a positive one.",
            FEATURE, nrow(h), ccnd[method == "IVW", n_snp]),
    sprintf("- **CCND1 is the one exposure whose methods disagree.** IVW gives OR %.3f (P = %.2f) while MR-Egger gives OR %.3f (P = %.3f) with an intercept at P = %.3f, i.e. directional pleiotropy on 5 instruments. It is reported, not interpreted.",
            ccnd[method == "IVW", OR], ccnd[method == "IVW", p],
            ccnd[method == "MR-Egger", OR], ccnd[method == "MR-Egger", p],
            ccnd[method == "MR-Egger", egger_int_p]),
    sprintf("- **(D) labels the three exposures of panel A only;** all %d instrumented genes are plotted. Exposure R2 exceeds outcome R2 for %d of %d, so the instruments act on expression first.",
            nrow(st), st[correct_direction == TRUE, .N], nrow(st)),
    "- Exposures are whole-blood eQTLs (eQTLGen), the outcome is GCST006867. The primary analysis is the LD-clumped one; the earlier distance-pruned pass over-called SERPINE1, PPARGC1A and PRKAA1 and is withdrawn (decision R21).",
    "",
    "## Statistics to quote",
    "",
    "| Quantity | Value |",
    "|---|---|",
    sprintf("| Genes with usable instruments | %d |", nrow(ivw)),
    sprintf("| Genes at FDR < %.2f | %d (%s) |", FDR_CUT, nrow(sig), paste(sig$gene, collapse = ", ")),
    sprintf("| Strongest estimate | %s, OR %.3f (95%% CI %.3f-%.3f), P = %s |", sig$gene[1],
            sig$OR[1], sig$OR_lo[1], sig$OR_hi[1], p_sci(sig$p[1])),
    sprintf("| %s (Tier 1) | %d SNPs, IVW OR %.3f (95%% CI %.3f-%.3f), P = %.2f |", TIER1[1],
            ivw[gene == TIER1[1], n_snp], ivw[gene == TIER1[1], OR], ivw[gene == TIER1[1], OR_lo],
            ivw[gene == TIER1[1], OR_hi], ivw[gene == TIER1[1], p]),
    sprintf("| %s (Tier 1) | %d SNPs, IVW OR %.3f (95%% CI %.3f-%.3f), P = %.2f |", TIER1[2],
            ivw[gene == TIER1[2], n_snp], ivw[gene == TIER1[2], OR], ivw[gene == TIER1[2], OR_lo],
            ivw[gene == TIER1[2], OR_hi], ivw[gene == TIER1[2], p]),
    sprintf("| %s heterogeneity | Cochran Q = %.2f, P = %.2f |", FEATURE,
            ft[method == "IVW", Q], ft[method == "IVW", Q_p]),
    sprintf("| %s pleiotropy | MR-Egger intercept P = %.2f |", FEATURE,
            ft[method == "MR-Egger", egger_int_p]),
    sprintf("| Steiger, correct direction | %d of %d genes |",
            st[correct_direction == TRUE, .N], nrow(st))
  )
  writeLines(l, file.path(fig_dir, "Fig5_paper_style_legend.md"))
  message("wrote ", file.path(fig_dir, "Fig5_paper_style_legend.md"))
}
