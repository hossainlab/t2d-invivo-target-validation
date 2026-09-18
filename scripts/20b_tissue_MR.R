# 20b_tissue_MR.R — MR with tissue-specific instruments and with LD-clumped whole-blood instruments.
# Run after 20a_tissue_instruments.R and the outcome-lookup shell step documented in its header.
#
# Exposures: GTEx liver (n = 208) and GTEx kidney cortex (n = 73) cis-eQTLs from the eQTL Catalogue, and
#            eQTLGen whole blood (n = 31,684) re-selected with LD clumping instead of distance pruning.
# Outcome:   type 2 diabetes, GCST006867 (Xue et al. 2018), harmonised summary statistics.
#
# Outputs: results/mr/mr_results_tissue.csv, mr_results_blood_ldclumped.csv, mr_blood_pruning_comparison.csv,
#          results/mr/mr_steiger_tissue.csv
#          figures/Fig5_MR/Fig5E_tissue_forest.pdf, Fig5F_blood_clumping_comparison.pdf

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source("scripts/mr_functions.R")
set.seed(20260916)
scratch <- Sys.getenv("MR_SCRATCH", unset = file.path(tempdir(), "mr"))
out_dir <- "results/mr"; fig_dir <- "results/supplementary_figures/mr"
  # draft panels; figures/ is owned by the 27-31 publication figure scripts. Writing here from an
  # analysis script silently overwrites them, because Windows filenames are case-insensitive.

for (d in c(out_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
N_CASE <- 62892; N_CTRL <- 596424; N_BLOOD <- 31684

gwas <- fread(file.path(scratch, "t2d_instruments_stage2.tsv"))
g <- gwas[, .(rsid = hm_rsid, out_ea = hm_effect_allele, out_oa = hm_other_allele, out_beta = hm_beta,
              out_se = standard_error, out_p = p_value, eaf = hm_effect_allele_frequency)]
g <- unique(g[!is.na(out_beta) & !is.na(out_se) & out_se > 0], by = "rsid")
flip <- function(a) c(A = "T", T = "A", C = "G", G = "C")[a]

harmonise <- function(x) {
  d <- merge(x, g, by = "rsid")
  d[, same   := exp_ea == out_ea & exp_oa == out_oa]
  d[, swap   := exp_ea == out_oa & exp_oa == out_ea]
  d[, same_f := exp_ea == flip(out_ea) & exp_oa == flip(out_oa)]
  d[, swap_f := exp_ea == flip(out_oa) & exp_oa == flip(out_ea)]
  d <- d[same | swap | same_f | swap_f]
  d[swap | swap_f, `:=`(out_beta = -out_beta, eaf = 1 - eaf)]
  d[, palindromic := (exp_ea == "A" & exp_oa == "T") | (exp_ea == "T" & exp_oa == "A") |
                     (exp_ea == "C" & exp_oa == "G") | (exp_ea == "G" & exp_oa == "C")]
  d[!(palindromic & !is.na(eaf) & eaf > 0.42 & eaf < 0.58)]
}

# ---- 1. Tissue instruments ---------------------------------------------------------------
tis <- fread(file.path(out_dir, "instruments_tissue.csv"))
th <- harmonise(tis)
message("Tissue instruments harmonised: ", nrow(th), " of ", nrow(tis))
th[, F_stat := (exp_beta / exp_se)^2]

mr_tis <- rbindlist(lapply(split(th, by = c("dataset", "gene")), function(x) {
  r <- run_mr(x); r[, dataset := x$dataset[1]][, p_threshold := x$p_threshold[1]][]
}))
prim_t <- mr_tis[method %in% c("IVW", "Wald ratio")]
prim_t[, FDR := p.adjust(p, "BH"), by = dataset]
mr_tis <- merge(mr_tis, prim_t[, .(dataset, gene, FDR)], by = c("dataset", "gene"), all.x = TRUE)
setorder(mr_tis, dataset, FDR, gene, method)
fwrite(mr_tis, file.path(out_dir, "mr_results_tissue.csv"))
print(prim_t[order(p), .(dataset, gene, method, n_snp, OR = round(OR, 3), lo = round(OR_lo, 3),
                         hi = round(OR_hi, 3), p = signif(p, 3), FDR = signif(FDR, 3), thr = p_threshold)])

st_t <- th[, {
  r2_exp <- sum(2 * pmin(pmax(maf, 0.001), 0.999) * (1 - pmin(pmax(maf, 0.001), 0.999)) * exp_beta^2)
  r2_out <- sum(get_r2_bin(out_beta, out_se, N_CASE, N_CTRL))
  .(n_snp = .N, r2_exposure = r2_exp, r2_outcome = r2_out, correct_direction = r2_exp > r2_out)
}, by = .(dataset, gene)]
fwrite(st_t, file.path(out_dir, "mr_steiger_tissue.csv"))

# ---- 2. Whole blood, LD-clumped ------------------------------------------------------------
bl <- fread(file.path(out_dir, "instruments_blood_ldclumped.csv"))
bh <- harmonise(bl)
bh <- bh[!is.na(eaf)]
bh[is.na(exp_n), exp_n := N_BLOOD]
bh[, denom := sqrt(2 * eaf * (1 - eaf) * (exp_n + exp_z^2))]
bh[, `:=`(exp_beta = exp_z / denom, exp_se = 1 / denom, F_stat = exp_z^2)]
bh <- bh[F_stat >= 10]
mr_bl <- rbindlist(lapply(split(bh, by = "gene"), run_mr))
prim_b <- mr_bl[method %in% c("IVW", "Wald ratio")]
prim_b[, FDR := p.adjust(p, "BH")]
mr_bl <- merge(mr_bl, prim_b[, .(gene, FDR)], by = "gene", all.x = TRUE)
setorder(mr_bl, FDR, gene, method)
fwrite(mr_bl, file.path(out_dir, "mr_results_blood_ldclumped.csv"))
print(prim_b[order(p), .(gene, method, n_snp, OR = round(OR, 3), lo = round(OR_lo, 3), hi = round(OR_hi, 3),
                         p = signif(p, 3), FDR = signif(FDR, 3), Q_p = signif(Q_p, 2))][seq_len(min(12, .N))])

# ---- 3. Clumped vs distance-pruned (script 19) ---------------------------------------------
old <- fread(file.path(out_dir, "mr_results.csv"))[method %in% c("IVW", "Wald ratio"),
      .(gene, n_snp_dist = n_snp, b_dist = b, se_dist = se, p_dist = p, FDR_dist = FDR)]
cmp <- merge(prim_b[, .(gene, n_snp_ld = n_snp, b_ld = b, se_ld = se, p_ld = p, FDR_ld = FDR)], old, by = "gene")
fwrite(cmp[order(p_ld)], file.path(out_dir, "mr_blood_pruning_comparison.csv"))
print(cmp[order(p_ld), .(gene, n_snp_ld, n_snp_dist, OR_ld = round(exp(b_ld), 3), OR_dist = round(exp(b_dist), 3),
                         p_ld = signif(p_ld, 3), p_dist = signif(p_dist, 3), FDR_ld = signif(FDR_ld, 3))][seq_len(min(12, .N))])

# ---- 4. Figures ------------------------------------------------------------------------------
if (nrow(prim_t)) {
  pt <- prim_t[order(dataset, b)]
  pt[, lab := factor(paste0(gene, " (", n_snp, ")"), levels = rev(unique(paste0(gene, " (", n_snp, ")"))))]
  p4e <- ggplot(pt, aes(OR, lab)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = OR_lo, xmax = OR_hi), height = 0.25) +
    geom_point(aes(colour = FDR < 0.05), size = 2.2) +
    facet_wrap(~ dataset, scales = "free_y") +
    scale_colour_manual(values = c(`TRUE` = "#C8322F", `FALSE` = "grey35"), name = "FDR < 0.05") +
    scale_x_continuous(trans = "log10") +
    labs(x = "OR for type 2 diabetes per SD of genetically predicted expression", y = NULL,
         title = "cis-MR with tissue-specific instruments (eQTL Catalogue, GTEx)",
         subtitle = "LD-clumped against 1000G EUR (r2 < 0.1, 500 kb); instrument count in brackets") +
    theme_bw(base_size = 9)
  ggsave(file.path(fig_dir, "Fig5E_tissue_forest.pdf"), p4e,
         width = 11, height = 0.28 * nrow(pt) + 2.5, limitsize = FALSE)
}
if (nrow(cmp)) {
  p4f <- ggplot(cmp, aes(exp(b_dist), exp(b_ld))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbar(aes(ymin = exp(b_ld - 1.96 * se_ld), ymax = exp(b_ld + 1.96 * se_ld)), colour = "grey75") +
    geom_errorbarh(aes(xmin = exp(b_dist - 1.96 * se_dist), xmax = exp(b_dist + 1.96 * se_dist)), colour = "grey75") +
    geom_point(aes(colour = FDR_ld < 0.05), size = 2) +
    ggrepel::geom_text_repel(aes(label = gene), size = 2.6, max.overlaps = 20) +
    scale_colour_manual(values = c(`TRUE` = "#C8322F", `FALSE` = "grey35"), name = "FDR < 0.05 (LD-clumped)") +
    labs(x = "OR, distance-pruned (script 19)", y = "OR, LD-clumped (1000G EUR)",
         title = "Whole-blood MR: LD clumping vs distance pruning") + theme_bw(base_size = 9)
  ggsave(file.path(fig_dir, "Fig5F_blood_clumping_comparison.pdf"), p4f, width = 7, height = 6)
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_20b.txt"))
