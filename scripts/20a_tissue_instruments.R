# 20a_tissue_instruments.R — tissue-specific cis-eQTL instruments + LD clumping against 1000G EUR.
#
# Answers the two limitations of script 19: whole-blood instruments, and distance pruning in place of LD
# clumping. Both are now done without any access token:
#   - tissue eQTLs: eQTL Catalogue tabix-indexed summary statistics, queried remotely by region
#       GTEx liver        QTD000266 (n = 208)
#       GTEx kidney cortex QTD000261 (n = 73)
#     (the GTEx portal's own eQTL API returns empty results, and the eQTL Catalogue REST API was retired
#      with HTTP 410, so the FTP tabix files are the remaining route)
#   - LD clumping: Ensembl REST /ld/human/<rsid>/1000GENOMES:phase_3:EUR, r2 < 0.1 within 500 kb.
# The same clumping is also applied to the eQTLGen whole-blood instruments of script 19, so the blood result
# can be reported LD-clumped rather than distance-pruned.
#
# After this script, fetch the outcome rows for the selected instruments:
#   gzip -dc <scratch>/mr/t2d_GCST006867.h.tsv.gz \
#     | awk -F'\t' 'NR==FNR{r[$1];next} FNR==1{print;next} ($2 in r)' <scratch>/mr/rsids_stage2.txt - \
#     > <scratch>/mr/t2d_instruments_stage2.tsv
# then run 20b_tissue_MR.R.
#
# Outputs: results/mr/gene_coords.csv
#          results/mr/tissue_eqtl_<dataset>.csv         all cis-eQTLs of the candidate genes in that tissue
#          results/mr/instruments_tissue.csv            clumped tissue instruments
#          results/mr/instruments_blood_ldclumped.csv   clumped eQTLGen instruments
#          <scratch>/mr/rsids_stage2.txt

suppressPackageStartupMessages({
  library(data.table); library(Rsamtools); library(GenomicRanges); library(IRanges)
  library(httr); library(jsonlite)
})
source("scripts/mr_functions.R")
set.seed(20260916)
scratch <- Sys.getenv("MR_SCRATCH", unset = file.path(tempdir(), "mr"))
out_dir <- "results/mr"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

P_STRICT <- 5e-8; P_RELAXED <- 1e-5; F_MIN <- 10; CIS_KB <- 250e3; LD_R2 <- 0.1
datasets <- data.table(
  dataset = c("GTEx_liver", "GTEx_kidney_cortex"),
  n       = c(208L, 73L),
  url     = c("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000266/QTD000266.all.tsv.gz",
              "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000261/QTD000261.all.tsv.gz"))
eqtl_cols <- c("molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples", "maf",
               "pvalue", "beta", "se", "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id",
               "median_tpm", "rsid")

genes <- readLines(file.path(scratch, "genes.txt"))
genes <- genes[nzchar(genes)]

# ---- 1. GRCh38 coordinates -------------------------------------------------------------
coord_file <- file.path(out_dir, "gene_coords.csv")
if (file.exists(coord_file)) {
  coords <- fread(coord_file)
} else {
  coords <- rbindlist(lapply(genes, function(g) {
    r <- tryCatch(fromJSON(content(GET(
      sprintf("https://rest.ensembl.org/lookup/symbol/homo_sapiens/%s?expand=0", g),
      add_headers(`Content-Type` = "application/json"), timeout(60)), "text", encoding = "UTF-8")),
      error = function(e) NULL)
    Sys.sleep(0.12)
    if (is.null(r) || is.null(r$id)) return(NULL)
    data.table(gene = g, ensg = r$id, chr = as.character(r$seq_region_name),
               start = as.integer(r$start), end = as.integer(r$end))
  }))
  fwrite(coords, coord_file)
}
message("Gene coordinates resolved: ", nrow(coords), " / ", length(genes))

# ---- 2. Pull cis-eQTLs from each tissue dataset -----------------------------------------
pull_tissue <- function(url, cd) {
  tf <- TabixFile(url)
  rbindlist(lapply(seq_len(nrow(cd)), function(i) {
    gr <- GRanges(cd$chr[i], IRanges(max(1, cd$start[i] - CIS_KB), cd$end[i] + CIS_KB))
    raw <- tryCatch(scanTabix(tf, param = gr)[[1]], error = function(e) character(0))
    if (!length(raw)) return(NULL)
    d <- fread(text = paste(raw, collapse = "\n"), header = FALSE, col.names = eqtl_cols,
               colClasses = list(character = c(1, 2, 4, 5, 6, 12, 16, 17, 19)))
    d <- d[gene_id == cd$ensg[i] & type == "SNP" & !is.na(rsid) & rsid != "NA" & !is.na(se) & se > 0]
    if (!nrow(d)) return(NULL)
    d[, gene := cd$gene[i]][]
  }))
}

tissue_all <- list()
for (k in seq_len(nrow(datasets))) {
  message("Querying ", datasets$dataset[k], " ...")
  d <- pull_tissue(datasets$url[k], coords)
  if (is.null(d) || !nrow(d)) { message("  no rows"); next }
  d[, dataset := datasets$dataset[k]][, n_exp := datasets$n[k]]
  fwrite(d, file.path(out_dir, paste0("tissue_eqtl_", datasets$dataset[k], ".csv")))
  message("  cis rows: ", nrow(d), " genes: ", uniqueN(d$gene),
          "; P<5e-8: ", d[pvalue < P_STRICT, .N], " in ", d[pvalue < P_STRICT, uniqueN(gene)], " genes")
  tissue_all[[datasets$dataset[k]]] <- d
}
tis <- rbindlist(tissue_all)

# ---- 3. Instrument selection per gene x dataset ------------------------------------------
setnames(tis, c("pvalue", "beta", "se", "chromosome", "position"),
         c("exp_p", "exp_beta", "exp_se", "chr", "pos"), skip_absent = TRUE)
tis[, F_stat := (exp_beta / exp_se)^2]
tis[, exp_ea := alt][, exp_oa := ref]

pick <- function(x) {
  s <- x[exp_p < P_STRICT & F_stat >= F_MIN]
  thr <- "5e-8"
  if (!nrow(s)) { s <- x[exp_p < P_RELAXED & F_stat >= F_MIN]; thr <- "1e-5 (relaxed)" }
  if (!nrow(s)) return(NULL)
  s[, p_threshold := thr][]
}
sel_tis <- rbindlist(lapply(split(tis, by = c("dataset", "gene")), pick))
message("Tissue instruments before clumping: ", nrow(sel_tis), " (",
        sel_tis[, uniqueN(paste(dataset, gene))], " gene x tissue pairs)")

# ---- 4. LD clumping (Ensembl, 1000G EUR) -------------------------------------------------
clump_group <- function(x) {
  x <- x[order(exp_p)]
  keep <- character(0); unresolved <- character(0); pool <- unique(x$rsid)
  while (length(pool)) {
    lead <- pool[1]; keep <- c(keep, lead); pool <- pool[-1]
    if (!length(pool)) break
    pr <- ld_proxies(lead, r2 = LD_R2, window_kb = 500)
    Sys.sleep(0.12)
    if (is.null(pr)) { unresolved <- c(unresolved, lead); next }
    pool <- setdiff(pool, pr[r2 >= LD_R2, rsid])
  }
  y <- x[rsid %in% keep]
  y[, ld_unresolved := rsid %in% unresolved][]
}
inst_tis <- rbindlist(lapply(split(sel_tis, by = c("dataset", "gene")), clump_group))
fwrite(inst_tis, file.path(out_dir, "instruments_tissue.csv"))
print(inst_tis[, .(n_snp = .N, min_p = signif(min(exp_p), 2), thr = p_threshold[1],
                   unresolved = sum(ld_unresolved)), by = .(dataset, gene)][order(dataset, -n_snp)])

# ---- 5. Re-clump the eQTLGen whole-blood instruments of script 19 -------------------------
blood <- fread(file.path(scratch, "eqtlgen_cis.tsv"))
setnames(blood, c("SNP", "AssessedAllele", "OtherAllele", "Zscore", "GeneSymbol", "NrSamples", "Pvalue",
                  "SNPChr", "SNPPos"),
         c("rsid", "exp_ea", "exp_oa", "exp_z", "gene", "exp_n", "exp_p", "chr", "pos"), skip_absent = TRUE)
blood_sel <- blood[exp_p < P_STRICT]
inst_blood <- rbindlist(lapply(split(blood_sel, by = "gene"), function(x) {
  x <- x[order(exp_p)][seq_len(min(.N, 400))]      # cap the clumping queue per gene
  clump_group(x)
}), fill = TRUE)
fwrite(inst_blood, file.path(out_dir, "instruments_blood_ldclumped.csv"))
print(inst_blood[, .(n_snp = .N), by = gene][order(-n_snp)][1:10])

# ---- 6. rsIDs needed from the outcome GWAS ------------------------------------------------
rs <- unique(c(inst_tis$rsid, inst_blood$rsid))
writeLines(rs, file.path(scratch, "rsids_stage2.txt"))
message("rsIDs written for the outcome lookup: ", length(rs))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_20a.txt"))
