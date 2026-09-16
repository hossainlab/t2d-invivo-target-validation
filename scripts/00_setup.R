# 00_setup.R — install packages required by docs/analysis_plan.md
# Run once: Rscript scripts/00_setup.R

options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 1200)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c(
  # Phase 1: microarray preprocessing
  "oligo", "pd.hg.u133a", "pd.hugene.1.1.st.v1",
  "hgu133a.db", "hugene11sttranscriptcluster.db",
  # Phase 5/9: enrichment, networks
  "STRINGdb",
  # Phase 6/7: single-cell
  "scDblFinder", "celda", "SingleR", "celldex", "glmGamPoi", "scuttle",
  "speckle", "UCell", "decoupleR", "miloR"
)
cran_pkgs <- c("RobustRankAggreg", "clustree", "renv")
gh_pkgs   <- c("immunogenomics/presto", "jinworks/CellChat", "powellgenomicslab/Nebulosa")

missing_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
missing_cran <- setdiff(cran_pkgs, rownames(installed.packages()))

if (length(missing_cran)) install.packages(missing_cran)
if (length(missing_bioc)) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

for (repo in gh_pkgs) {
  pkg <- basename(repo)
  if (!requireNamespace(pkg, quietly = TRUE)) {
    tryCatch(remotes::install_github(repo, upgrade = "never"),
             error = function(e) message("GitHub install failed: ", repo, " — ", conditionMessage(e)))
  }
}

all_pkgs <- c(bioc_pkgs, cran_pkgs, basename(gh_pkgs))
status <- setNames(all_pkgs %in% rownames(installed.packages()), all_pkgs)
print(status)
if (!all(status)) message("Missing: ", paste(names(status)[!status], collapse = ", "))
