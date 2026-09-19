# 26_funnel_permutation.R — is the funnel's output better than chance, and how much does it depend on the
# DEG threshold?
#
# Two validity checks that the study did not have:
#
# A. Threshold sensitivity. The Fig 1 candidate set uses nominal P < 0.05 because the paper-strict rule
#    (FDR < 0.05, single strongest module) leaves one gene. Both Tier 1 genes sit just outside FDR
#    (PPARGC1A 0.054, CCND1 0.058), so the candidate set is recomputed under a grid of DEG and module rules
#    and the fate of every Tier 1/2 gene is reported.
#
# B. Empirical FDR of the whole funnel. The disease labels are permuted within dataset (so the batch
#    structure is preserved), the entire funnel is re-run on each permutation, and the observed survivor
#    counts are compared with the null distribution:
#       step 3  DE genes on KEGG hsa04152
#       step 4  + reachable from the AMPK node by a drawn directed path   (graph fixed, from script 25)
#       step 5  + the drawn sign predicts the observed direction
#       step 6  + in the fixed set of genes confirmed in the mouse atlas  (mouse data are not permuted)
#       step 7  + inside the permutation's own Fig 1 candidate set
#    Step 5 is the informative one: under the null roughly half of the reachable DE genes should agree with
#    the map by chance, so the observed fraction is the test of whether the map explains the data.
#
# Outputs: results/consistency/26_threshold_sensitivity.csv
#          results/consistency/26_funnel_permutation.csv   (per-permutation counts)
#          results/consistency/26_funnel_summary.csv       (observed vs null, empirical P)
# Figure:  results/supplementary_figures/consistency/26_funnel_null.pdf

suppressPackageStartupMessages({ library(data.table); library(limma); library(ggplot2); library(patchwork) })
set.seed(20260916)
out_dir <- "results/consistency"; fig_dir <- "results/supplementary_figures/consistency"
for (d in c(out_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
N_PERM <- 1000L

expr <- readRDS("results/bulk/expr_merged_raw.rds")
meta <- fread("results/bulk/meta_bulk.csv")[match(colnames(expr), GSM)]
meta[, condition := factor(condition, levels = c("Control", "T2D"))]
stopifnot(identical(meta$GSM, colnames(expr)))
mm    <- fread("results/bulk/WGCNA_gene_module_GS_MM.csv")
mt    <- fread("results/bulk/WGCNA_module_trait.csv")
hyper <- fread("results/bulk/geneset_hyperglycemia.csv")$gene
paths <- fread("results/pathway_selection/25_map_paths.csv")
tiers <- fread("results/targets/target_tiers.csv")

# Key modules are decided once, by script 03 (decision M11b: T2D p < 0.1 AND replicated in both
# cohorts). Reading that decision here instead of re-deriving it keeps every downstream step on the
# same module set; the rule used to be copied into four scripts, which is how the vacuous
# abs(r_T2D) > abs(r_Dataset) clause survived unnoticed.
read_key_modules <- function(dir = "results/bulk") {
  f <- file.path(dir, "WGCNA_key_module_selection.csv")
  if (!file.exists(f)) stop("missing ", f, "; run scripts/03_bulk_WGCNA.R first")
  fread(f)[key == TRUE][order(p_T2D), module]
}
key_modules <- read_key_modules()
strongest   <- mt[order(-r_T2D)][1, module]
module_genes <- list(`7 key modules` = mm[module %in% key_modules, gene],
                     `strongest module only` = mm[module == strongest, gene])
map_genes   <- paths[on_map == TRUE, gene]
reachable   <- paths[!is.na(n_steps), .(gene, predicted = predicted_if_AMPK_suppressed)]
sc_confirmed <- tiers[!is.na(sc_P), gene]

fit_deg <- function(cond) {
  des <- model.matrix(~ cond + meta$dataset)
  colnames(des) <- make.names(colnames(des))
  fit <- eBayes(lmFit(expr, des), trend = TRUE)
  as.data.table(topTable(fit, coef = 2, number = Inf), keep.rownames = "gene")
}
obs <- fit_deg(meta$condition)

# ---- A. Threshold sensitivity -------------------------------------------------------------
deg_rules <- list(
  `P < 0.05 & |logFC| > 0.5 (used)` = function(d) d[P.Value < 0.05 & abs(logFC) > 0.5, gene],
  `P < 0.01 & |logFC| > 0.5`        = function(d) d[P.Value < 0.01 & abs(logFC) > 0.5, gene],
  `FDR < 0.10 & |logFC| > 0.5`      = function(d) d[adj.P.Val < 0.10 & abs(logFC) > 0.5, gene],
  `FDR < 0.05 & |logFC| > 0.5`      = function(d) d[adj.P.Val < 0.05 & abs(logFC) > 0.5, gene],
  `FDR < 0.05 (no logFC cutoff)`    = function(d) d[adj.P.Val < 0.05, gene])
watch <- c("PPARGC1A", "CCND1", "PFKFB3", "CDKN1A", "IRS2", "VSNL1", "LZTFL1", "GAS6", "SREBF2")
sens <- rbindlist(lapply(names(deg_rules), function(dn) rbindlist(lapply(names(module_genes), function(mn) {
  cand <- Reduce(intersect, list(deg_rules[[dn]](obs), module_genes[[mn]], hyper))
  on_m <- intersect(cand, map_genes)
  data.table(deg_rule = dn, module_rule = mn, n_candidates = length(cand),
             n_on_AMPK_map = length(on_m),
             n_reachable = length(intersect(cand, reachable$gene)),
             kept = paste(intersect(watch, cand), collapse = ";"),
             candidates = paste(cand, collapse = ";"))
}))))
fwrite(sens, file.path(out_dir, "26_threshold_sensitivity.csv"))
print(sens[, .(deg_rule, module_rule, n_candidates, n_on_AMPK_map, n_reachable, kept)])

# ---- B. Funnel permutation ------------------------------------------------------------------
funnel_counts <- function(d, cand_genes) {
  de <- d[P.Value < 0.05 & abs(logFC) > 0.5]
  s3 <- intersect(de$gene, map_genes)
  s4 <- intersect(s3, reachable$gene)
  if (!length(s4)) return(c(s3 = length(s3), s4 = 0, s5 = 0, s6 = 0, s7 = 0))
  r <- merge(de[gene %in% s4, .(gene, logFC)], reachable, by = "gene")
  r[, match := (predicted == "up" & logFC > 0) | (predicted == "down" & logFC < 0)]
  s5 <- r[match == TRUE, gene]
  s6 <- intersect(s5, sc_confirmed)
  s7 <- intersect(s6, cand_genes)
  c(s3 = length(s3), s4 = length(s4), s5 = length(s5), s6 = length(s6), s7 = length(s7))
}
obs_cand <- Reduce(intersect, list(deg_rules[[1]](obs), module_genes[[1]], hyper))
obs_counts <- funnel_counts(obs, obs_cand)
obs_frac5 <- if (obs_counts["s4"] > 0) obs_counts["s5"] / obs_counts["s4"] else NA_real_
message("Observed: ", paste(names(obs_counts), obs_counts, sep = "=", collapse = ", "),
        "; sign-match fraction ", round(obs_frac5, 3))

perm <- rbindlist(lapply(seq_len(N_PERM), function(i) {
  set.seed(40000 + i)
  cond <- meta$condition
  for (ds in unique(meta$dataset)) {                 # permute within dataset: batch structure preserved
    idx <- which(meta$dataset == ds)
    cond[idx] <- sample(cond[idx])
  }
  d <- fit_deg(cond)
  cand <- Reduce(intersect, list(d[P.Value < 0.05 & abs(logFC) > 0.5, gene], module_genes[[1]], hyper))
  as.data.table(as.list(c(perm = i, funnel_counts(d, cand), n_cand = length(cand))))
}))
fwrite(perm, file.path(out_dir, "26_funnel_permutation.csv"))

summ <- rbindlist(lapply(c("s3", "s4", "s5", "s6", "s7"), function(s) data.table(
  step = s, observed = obs_counts[[s]], null_mean = mean(perm[[s]]), null_q95 = quantile(perm[[s]], 0.95),
  empirical_p = (sum(perm[[s]] >= obs_counts[[s]]) + 1) / (N_PERM + 1))))
perm[, frac5 := fifelse(s4 > 0, s5 / s4, NA_real_)]
summ <- rbind(summ, data.table(step = "sign-match fraction (s5/s4)", observed = obs_frac5,
                               null_mean = mean(perm$frac5, na.rm = TRUE),
                               null_q95 = quantile(perm$frac5, 0.95, na.rm = TRUE),
                               empirical_p = (sum(perm$frac5 >= obs_frac5, na.rm = TRUE) + 1) / (N_PERM + 1)))
summ <- rbind(summ, data.table(step = "Fig 1 candidate-set size", observed = length(obs_cand),
                               null_mean = mean(perm$n_cand), null_q95 = quantile(perm$n_cand, 0.95),
                               empirical_p = (sum(perm$n_cand >= length(obs_cand)) + 1) / (N_PERM + 1)))
fwrite(summ, file.path(out_dir, "26_funnel_summary.csv"))
print(summ)

labs <- c(s3 = "step 3: DE on the AMPK map", s4 = "step 4: + reachable from AMPK",
          s5 = "step 5: + sign matches", s6 = "step 6: + confirmed in mouse", s7 = "step 7: + Fig 1 candidate")
pl <- melt(perm[, c("s3", "s4", "s5", "s6", "s7")], measure.vars = names(labs),
           variable.name = "step", value.name = "n")
pl[, step := factor(labs[as.character(step)], levels = labs)]
obs_dt <- data.table(step = factor(labs, levels = labs), n = as.numeric(obs_counts[names(labs)]))
p1 <- ggplot(pl, aes(n)) + geom_histogram(binwidth = 1, fill = "grey75", colour = NA) +
  geom_vline(data = obs_dt, aes(xintercept = n), colour = "#C8322F", linewidth = 0.9) +
  facet_wrap(~ step, scales = "free", nrow = 1) +
  labs(x = "genes surviving the step", y = paste0("permutations (n = ", N_PERM, ")"),
       title = "Funnel under permuted disease labels (red = observed)") + theme_bw(base_size = 8)
p2 <- ggplot(perm[!is.na(frac5)], aes(frac5)) + geom_histogram(bins = 30, fill = "grey75", colour = NA) +
  geom_vline(xintercept = obs_frac5, colour = "#C8322F", linewidth = 0.9) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey40") +
  labs(x = "fraction of reachable DE genes whose direction the map predicts", y = "permutations",
       title = "Does the map explain the directions? (dashed = coin flip)") + theme_bw(base_size = 8)
ggsave(file.path(fig_dir, "26_funnel_null.pdf"), p1 / p2, width = 13, height = 7)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_26.txt"))
