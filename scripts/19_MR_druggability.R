# 19_MR_druggability.R — two-sample cis-Mendelian randomization of the candidate genes on type 2 diabetes
# (reference paper Fig 4 analogue: "is the target causal for the disease, and is it druggable?").
#
# Paper: exposure = PAK1 expression, outcome = RA (IVW, simple/weighted median, simple/weighted mode,
#        Cochran's Q, MR-Egger intercept, Steiger filtering), plus the Finan et al. 2017 druggable-genome list.
# Here:  exposure = cis-eQTL of each candidate gene in eQTLGen whole blood (n = 31,684; Vosa et al. 2021),
#        outcome  = type 2 diabetes GWAS, Xue et al. 2018 (GCST006867, 62,892 cases / 596,424 controls,
#        European), GWAS Catalog harmonised summary statistics (GRCh38, so positions match eQTLGen rsIDs).
#
# Declared deviations (Methods):
#   - Whole-blood eQTLs are used because the GTEx tissue-level eQTL API returns no data at present; the
#     instrument is therefore not liver- or kidney-specific.
#   - LD clumping uses distance pruning (keep the lowest-P SNP per 250 kb window) because no LD reference
#     panel or LD server token is available offline. Cis instruments in one locus are correlated, so
#     single-instrument genes are reported as Wald ratios and multi-instrument results are read with that
#     caveat; results are also reported for a stricter 1 Mb pruning window as a sensitivity analysis.
#   - eQTLGen publishes Z scores, so beta and se are reconstructed as beta = z / sqrt(2f(1-f)(n + z^2)) and
#     se = 1 / sqrt(2f(1-f)(n + z^2)) (Zhu et al. 2016), with f the effect-allele frequency taken from the
#     harmonised outcome file.
#
# Inputs (prepared by the shell step, see docs/pathway_selection.md):
#   <scratch>/mr/eqtlgen_cis.tsv        cis-eQTLs of the candidate genes
#   <scratch>/mr/t2d_GCST006867.h.tsv.gz  harmonised T2D summary statistics
# Outputs: results/mr/mr_instruments.csv, mr_results.csv, mr_sensitivity.csv, mr_steiger.csv
#          figures/Fig5_MR/Fig5A_forest.pdf, Fig5B_scatter_<gene>.pdf, Fig5C_leaveoneout_<gene>.pdf,
#          figures/Fig5_MR/Fig5D_funnel_<gene>.pdf

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(20260916)
scratch <- Sys.getenv("MR_SCRATCH", unset = file.path(tempdir(), "mr"))
out_dir <- "results/mr"; fig_dir <- "results/supplementary_figures/mr"
  # draft panels; figures/ is owned by the 27-31 publication figure scripts. Writing here from an
  # analysis script silently overwrites them, because Windows filenames are case-insensitive.

for (d in c(out_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

N_EXP  <- 31684        # eQTLGen cis discovery sample size (per-SNP NrSamples is used when present)
N_CASE <- 62892; N_CTRL <- 596424
PRUNE_KB <- 250e3
F_MIN <- 10

eqtl <- fread(file.path(scratch, "eqtlgen_cis.tsv"))
gwas <- fread(file.path(scratch, "t2d_instruments.tsv"))   # pre-filtered to the eQTL rsIDs
message("cis-eQTL rows: ", nrow(eqtl), " for ", uniqueN(eqtl$GeneSymbol), " genes; outcome rows: ", nrow(gwas))

# ---- 1. Harmonise ----------------------------------------------------------------------
setnames(eqtl, c("SNP", "AssessedAllele", "OtherAllele", "Zscore", "GeneSymbol", "NrSamples", "Pvalue",
                 "SNPChr", "SNPPos"),
         c("rsid", "exp_ea", "exp_oa", "exp_z", "gene", "exp_n", "exp_p", "chr", "pos"), skip_absent = TRUE)
g <- gwas[, .(rsid = hm_rsid, out_ea = hm_effect_allele, out_oa = hm_other_allele,
              out_beta = hm_beta, out_se = standard_error, out_p = p_value,
              eaf = hm_effect_allele_frequency)]
g <- g[!is.na(out_beta) & !is.na(out_se) & out_se > 0]
d <- merge(eqtl[, .(rsid, chr, pos, gene, exp_ea, exp_oa, exp_z, exp_p, exp_n)], g, by = "rsid")
message("SNPs present in both: ", nrow(d))

flip <- function(a) c(A = "T", T = "A", C = "G", G = "C")[a]
d[, same := exp_ea == out_ea & exp_oa == out_oa]
d[, swap := exp_ea == out_oa & exp_oa == out_ea]
d[, same_f := exp_ea == flip(out_ea) & exp_oa == flip(out_oa)]
d[, swap_f := exp_ea == flip(out_oa) & exp_oa == flip(out_ea)]
d <- d[same | swap | same_f | swap_f]
d[swap | swap_f, `:=`(out_beta = -out_beta, eaf = 1 - eaf)]        # align to the eQTL effect allele
d[, eaf := pmin(pmax(eaf, 0.001), 0.999)]
d <- d[!is.na(eaf)]
# palindromic SNPs with ambiguous frequency cannot be aligned safely
d[, palindromic := (exp_ea == "A" & exp_oa == "T") | (exp_ea == "T" & exp_oa == "A") |
                   (exp_ea == "C" & exp_oa == "G") | (exp_ea == "G" & exp_oa == "C")]
d <- d[!(palindromic & eaf > 0.42 & eaf < 0.58)]

# ---- 2. beta / se of the exposure from Z, n and allele frequency -----------------------
d[is.na(exp_n), exp_n := N_EXP]
d[, denom := sqrt(2 * eaf * (1 - eaf) * (exp_n + exp_z^2))]
d[, `:=`(exp_beta = exp_z / denom, exp_se = 1 / denom)]
d[, F_stat := (exp_beta / exp_se)^2]
message("Harmonised SNPs: ", nrow(d), " (", uniqueN(d$gene), " genes)")

# ---- 3. Instrument selection: genome-wide-significant cis-eQTLs, distance-pruned -------
prune <- function(x, kb) {
  x <- x[order(exp_p)]
  keep <- logical(nrow(x))
  for (i in seq_len(nrow(x))) {
    if (!any(keep & abs(x$pos[keep] - x$pos[i]) < kb & x$chr[keep] == x$chr[i])) keep[i] <- TRUE
  }
  x[keep]
}
sel <- d[exp_p < 5e-8 & F_stat >= F_MIN]
inst <- rbindlist(lapply(split(sel, by = "gene"), prune, kb = PRUNE_KB))
inst_1mb <- rbindlist(lapply(split(sel, by = "gene"), prune, kb = 1e6))
fwrite(inst[order(gene, exp_p)], file.path(out_dir, "mr_instruments.csv"))
print(inst[, .(n_snp = .N, min_p = min(exp_p), mean_F = round(mean(F_stat))), by = gene][order(-n_snp)])

# ---- 4. MR estimators -------------------------------------------------------------------
ivw <- function(b_exp, se_exp, b_out, se_out) {
  w <- 1 / se_out^2
  fit <- lm(b_out ~ 0 + b_exp, weights = w)
  b <- unname(coef(fit)[1]); se_fe <- summary(fit)$coef[1, 2] / summary(fit)$sigma
  q <- sum(w * (b_out - b * b_exp)^2); df <- length(b_exp) - 1
  se_re <- if (df > 0) se_fe * max(1, sqrt(q / df)) else se_fe
  list(b = b, se = se_re, p = 2 * pnorm(-abs(b / se_re)), Q = q, Q_df = df,
       Q_p = if (df > 0) pchisq(q, df, lower.tail = FALSE) else NA_real_)
}
egger <- function(b_exp, se_exp, b_out, se_out) {
  if (length(b_exp) < 3) return(NULL)
  s <- sign(b_exp); be <- abs(b_exp); bo <- b_out * s
  fit <- summary(lm(bo ~ be, weights = 1 / se_out^2))
  list(b = fit$coef[2, 1], se = fit$coef[2, 2] / min(1, fit$sigma), p = fit$coef[2, 4],
       int = fit$coef[1, 1], int_se = fit$coef[1, 2], int_p = fit$coef[1, 4])
}
w_median <- function(b_exp, se_exp, b_out, se_out, nboot = 1000) {
  if (length(b_exp) < 3) return(NULL)
  ratio <- b_out / b_exp; w <- (se_out / b_exp)^-2
  wmed <- function(r, w) {
    o <- order(r); r <- r[o]; w <- w[o]; cw <- cumsum(w) - 0.5 * w; cw <- cw / sum(w)
    below <- max(which(cw < 0.5))
    r[below] + (r[below + 1] - r[below]) * (0.5 - cw[below]) / (cw[below + 1] - cw[below])
  }
  b <- wmed(ratio, w)
  bs <- replicate(nboot, wmed(rnorm(length(b_out), b_out, se_out) / b_exp, w))
  list(b = b, se = sd(bs), p = 2 * pnorm(-abs(b / sd(bs))))
}
w_mode <- function(b_exp, se_exp, b_out, se_out, phi = 1, nboot = 1000) {
  if (length(b_exp) < 3) return(NULL)
  mode_of <- function(r, s) {
    w <- 1 / s^2; bw <- phi * sd(r) * length(r)^(-1/5)
    dd <- density(r, weights = w / sum(w), bw = bw); dd$x[which.max(dd$y)]
  }
  ratio <- b_out / b_exp; se_r <- abs(se_out / b_exp)
  b <- mode_of(ratio, se_r)
  bs <- replicate(nboot, mode_of(rnorm(length(b_out), b_out, se_out) / b_exp, se_r))
  list(b = b, se = sd(bs), p = 2 * pnorm(-abs(b / sd(bs))))
}
run_mr <- function(x) {
  res <- list()
  add <- function(m, r) if (!is.null(r)) res[[length(res) + 1L]] <<-
    data.table(gene = x$gene[1], method = m, n_snp = nrow(x), b = r$b, se = r$se, p = r$p,
               OR = exp(r$b), OR_lo = exp(r$b - 1.96 * r$se), OR_hi = exp(r$b + 1.96 * r$se),
               Q = if (!is.null(r$Q)) r$Q else NA_real_, Q_p = if (!is.null(r$Q_p)) r$Q_p else NA_real_,
               egger_int = if (!is.null(r$int)) r$int else NA_real_,
               egger_int_p = if (!is.null(r$int_p)) r$int_p else NA_real_)
  if (nrow(x) == 1) {
    b <- x$out_beta / x$exp_beta; se <- abs(x$out_se / x$exp_beta)
    add("Wald ratio", list(b = b, se = se, p = 2 * pnorm(-abs(b / se))))
  } else {
    add("IVW", ivw(x$exp_beta, x$exp_se, x$out_beta, x$out_se))
    add("MR-Egger", egger(x$exp_beta, x$exp_se, x$out_beta, x$out_se))
    add("Weighted median", w_median(x$exp_beta, x$exp_se, x$out_beta, x$out_se))
    add("Weighted mode", w_mode(x$exp_beta, x$exp_se, x$out_beta, x$out_se))
  }
  rbindlist(res)
}
mr <- rbindlist(lapply(split(inst, by = "gene"), run_mr))
prim <- mr[method %in% c("IVW", "Wald ratio")]
prim[, FDR := p.adjust(p, "BH")]
mr <- merge(mr, prim[, .(gene, FDR)], by = "gene", all.x = TRUE)
setorder(mr, FDR, gene, method)
fwrite(mr, file.path(out_dir, "mr_results.csv"))
print(mr[method %in% c("IVW", "Wald ratio")][order(p), .(gene, method, n_snp, OR = round(OR, 3),
      OR_lo = round(OR_lo, 3), OR_hi = round(OR_hi, 3), p = signif(p, 3), FDR = signif(FDR, 3),
      Q_p = signif(Q_p, 2))])

# ---- 4b. Cross-check against the MendelianRandomization package -------------------------
if (requireNamespace("MendelianRandomization", quietly = TRUE)) {
  MR <- asNamespace("MendelianRandomization")
  chk <- rbindlist(lapply(split(inst, by = "gene"), function(x) {
    if (nrow(x) < 3) return(NULL)
    obj <- MR$mr_input(bx = x$exp_beta, bxse = x$exp_se, by = x$out_beta, byse = x$out_se)
    a <- MR$mr_allmethods(obj, method = "main")@Values
    data.table(gene = x$gene[1], method = a[["Method"]], b = a[["Estimate"]], se = a[["Std Error"]],
               p = a[["P-value"]])
  }))
  if (nrow(chk)) fwrite(chk, file.path(out_dir, "mr_results_MRpackage.csv"))
}

# ---- 5. Sensitivity: 1 Mb pruning, and leave-one-out --------------------------------------
mr_1mb <- rbindlist(lapply(split(inst_1mb, by = "gene"), run_mr))[method %in% c("IVW", "Wald ratio")]
setnames(mr_1mb, c("b", "se", "p", "n_snp"), c("b_1mb", "se_1mb", "p_1mb", "n_snp_1mb"))
sens <- merge(prim[, .(gene, n_snp, b, se, p)], mr_1mb[, .(gene, n_snp_1mb, b_1mb, se_1mb, p_1mb)],
              by = "gene", all.x = TRUE)
loo <- rbindlist(lapply(split(inst, by = "gene"), function(x) {
  if (nrow(x) < 3) return(NULL)
  rbindlist(lapply(seq_len(nrow(x)), function(i) {
    r <- ivw(x$exp_beta[-i], x$exp_se[-i], x$out_beta[-i], x$out_se[-i])
    data.table(gene = x$gene[1], dropped = x$rsid[i], b = r$b, se = r$se, p = r$p)
  }))
}))
fwrite(sens, file.path(out_dir, "mr_sensitivity.csv"))
if (nrow(loo)) fwrite(loo, file.path(out_dir, "mr_leaveoneout.csv"))

# ---- 6. Steiger directionality ----------------------------------------------------------
get_r2_bin <- function(b, se, ncase, nctrl) {           # Lee et al. liability-scale approximation
  k <- ncase / (ncase + nctrl); n <- ncase + nctrl
  z <- b / se; r2 <- z^2 / (z^2 + n)
  r2 * (k * (1 - k)) / (dnorm(qnorm(k))^2) * 1          # scale to liability
}
st <- inst[, {
  r2_exp <- sum(2 * eaf * (1 - eaf) * exp_beta^2)
  r2_out <- sum(get_r2_bin(out_beta, out_se, N_CASE, N_CTRL))
  n1 <- exp_n[1]; n2 <- N_CASE + N_CTRL
  z1 <- atanh(sqrt(min(r2_exp, 0.999))); z2 <- atanh(sqrt(min(r2_out, 0.999)))
  se_d <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  .(n_snp = .N, r2_exposure = r2_exp, r2_outcome = r2_out,
    correct_direction = r2_exp > r2_out, steiger_p = 2 * pnorm(-abs((z1 - z2) / se_d)))
}, by = gene]
fwrite(st, file.path(out_dir, "mr_steiger.csv"))
print(st[order(-r2_exposure)])

# ---- 7. Figures (paper Fig 4 analogue) --------------------------------------------------
pf <- prim[order(b)]
pf[, gene := factor(gene, levels = gene)]
p4a <- ggplot(pf, aes(OR, gene)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = OR_lo, xmax = OR_hi), height = 0.25) +
  geom_point(aes(colour = FDR < 0.05), size = 2.4) +
  scale_colour_manual(values = c(`TRUE` = "#C8322F", `FALSE` = "grey35"), name = "FDR < 0.05") +
  scale_x_continuous(trans = "log10") +
  labs(x = "OR for type 2 diabetes per SD of genetically predicted expression", y = NULL,
       title = "cis-MR of the candidate genes on type 2 diabetes",
       subtitle = "eQTLGen whole blood -> Xue et al. 2018 (GCST006867); IVW, or Wald ratio for single-instrument genes") +
  theme_bw(base_size = 10)
ggsave(file.path(fig_dir, "Fig5A_forest.pdf"), p4a, width = 8.5, height = 0.3 * nrow(pf) + 2.5, limitsize = FALSE)

for (gn in prim[FDR < 0.05 | gene %in% c("CDKN1A", "CCND1"), gene]) {
  x <- inst[gene == gn]
  if (nrow(x) >= 2) {
    b <- prim[gene == gn, b]
    ggsave(file.path(fig_dir, paste0("Fig5B_scatter_", gn, ".pdf")),
           ggplot(x, aes(exp_beta, out_beta)) +
             geom_errorbar(aes(ymin = out_beta - out_se, ymax = out_beta + out_se), width = 0, colour = "grey70") +
             geom_errorbarh(aes(xmin = exp_beta - exp_se, xmax = exp_beta + exp_se), height = 0, colour = "grey70") +
             geom_point(size = 2) + geom_abline(slope = b, intercept = 0, colour = "#C8322F") +
             labs(x = paste0("SNP effect on ", gn, " expression"), y = "SNP effect on T2D (log OR)",
                  title = paste0(gn, ": ", nrow(x), " instruments")) + theme_bw(base_size = 10),
           width = 5, height = 4.5)
    ggsave(file.path(fig_dir, paste0("Fig5D_funnel_", gn, ".pdf")),
           ggplot(x, aes(out_beta / exp_beta, abs(exp_beta / exp_se))) + geom_point(size = 2) +
             geom_vline(xintercept = prim[gene == gn, b], colour = "#C8322F") +
             labs(x = "Wald ratio", y = "instrument strength |beta/se|", title = paste0(gn, " funnel")) +
             theme_bw(base_size = 10), width = 5, height = 4.5)
  }
  if (nrow(loo) && gn %in% loo$gene) {
    l <- loo[gene == gn][order(b)]
    l[, dropped := factor(dropped, levels = dropped)]
    ggsave(file.path(fig_dir, paste0("Fig5C_leaveoneout_", gn, ".pdf")),
           ggplot(l, aes(b, dropped)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
             geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0.2) +
             geom_point() + labs(x = "IVW log OR, SNP removed", y = NULL,
                                 title = paste0(gn, ": leave-one-out")) + theme_bw(base_size = 9),
           width = 6, height = 0.25 * nrow(l) + 2)
  }
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_19.txt"))
