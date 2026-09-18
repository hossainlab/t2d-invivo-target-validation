# figure_style.R - shared publication style for everything in figures/.
#
# Sourced by the 27-31 figure scripts. Nothing here touches data; it only fixes how figures are drawn
# and written, so that every panel in figures/ obeys the same contract:
#
#   * drawn at FINAL PRINT SIZE in millimetres (single column 89 mm, 1.5 column 120 mm, double 180 mm),
#     never drawn large and scaled down - that is what makes text land at 3 pt in a proof
#   * Arial 5-7 pt, fonts embedded, text left as live text in the PDF
#   * no titles, subtitles or explanatory prose inside the axes; that belongs in the figure legend
#   * perceptually uniform, colour-vision-safe colour maps (viridis); no rainbow/jet, no pure red/green
#   * figures/ receives PDF only. TIFF submission rasters and any serialized intermediates go to
#     results/figure_exports/
#
# Journal switches: TAG_CASE "lower" = Nature (bold a, b, c), "upper" = Cell (A, B, C).

suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(patchwork)
})

FONT     <- if (requireNamespace("systemfonts", quietly = TRUE) &&
                "Arial" %in% systemfonts::system_fonts()$family) "Arial" else "sans"
BASE     <- 7
TAG_CASE <- "lower"
W_1COL   <- 89     # mm
W_15COL  <- 120    # mm
W_2COL   <- 180    # mm
EXPORT_DIR <- "results/figure_exports"

ppt <- 1 / .pt     # geom_text `size` is mm; multiply a pt value by this

# ---- palettes -------------------------------------------------------------------------------------
# two-group contrast used across the paper (Control/lean vs T2D/STZ)
col_grp  <- c(Control = "#4E79A7", T2D = "#E15759", STZ = "#E15759", Lean = "#4E79A7")
# direction of change
col_dir  <- c(Up = "#C0392B", Down = "#2C6FAD", NS = "grey80",
              Significant = "#C0392B", "Not significant" = "grey80")
# qualitative set: muted, distinguishable in greyscale and under deuteranopia
pal_qual <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1",
              "#FF9DA7", "#9C755F", "#BAB0AC", "#86BCB6", "#D37295", "#8CD17D", "#499894",
              "#F1CE63", "#A0CBE8")
# diverging map for correlation/z heatmaps: blue-white-red, equal luminance ramp either side
pal_div  <- c("#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#EF8A62", "#B2182B")

# ---- formatters -----------------------------------------------------------------------------------
nfmt <- function(n) formatC(n, big.mark = ",", format = "d")
mfmt <- function(x, d = 2) sub("-", "−", formatC(x, format = "f", digits = d))
sup_digits <- function(n) {
  g <- c("⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹")
  paste(g[as.integer(strsplit(as.character(abs(n)), "")[[1]]) + 1], collapse = "")
}
p_short <- function(p) {
  if (is.na(p)) return("n.d.")
  if (p < 2.2e-16) return("P < 2.2 × 10⁻¹⁶")
  if (p < 0.001) {
    e <- floor(log10(p))
    return(sprintf("P = %s × 10⁻%s", formatC(p / 10^e, format = "f", digits = 1), sup_digits(-e)))
  }
  sprintf("P = %s", formatC(p, format = "f", digits = min(3, max(2, -floor(log10(p)) + 1))))
}
stars <- function(p) {
  if (is.na(p)) "n.s." else if (p < 1e-4) "****" else if (p < 1e-3) "***" else
    if (p < 1e-2) "**" else if (p < 0.05) "*" else "n.s."
}
# wrap long ontology/pathway terms so they do not run off a 89 mm panel
wrap_term <- function(x, width = 42) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), "", USE.NAMES = FALSE)
}
# sentence case for GO/KEGG/Hallmark terms, which come back SHOUTING or lower_underscored.
# Acronyms that must stay upper case are restored afterwards.
tidy_term <- function(x) {
  x <- gsub("_", " ", x)
  x <- sub("^(GOBP|GOCC|GOMF|GO|HALLMARK|KEGG|REACTOME)\\s+", "", x, ignore.case = TRUE)
  allcaps <- !grepl("[a-z]", x)               # untouched lower-case terms keep their own casing
  x[allcaps] <- tolower(x[allcaps])
  x <- paste0(toupper(substring(x, 1, 1)), substring(x, 2))
  acr <- c("dna", "rna", "mrna", "atp", "amp", "ampk", "tnfa", "tgf beta", "il2", "il6", "stat3",
           "stat5", "jak", "nfkb", "mtorc1", "myc", "e2f", "g2m", "uv", "ros", "hif 1", "foxo",
           "p53", "pi3k", "akt", "kras", "wnt", "notch", "hedgehog", "epithelial mesenchymal")
  for (a in acr) {
    x <- gsub(paste0("\\b", a, "\\b"), toupper(a), x, ignore.case = TRUE)
  }
  x <- gsub("\\bIl([26])\\b", "IL\\1", x)
  x <- gsub("\\bP53\\b", "p53", x)          # gene-product convention, not an acronym
  x <- gsub("\\bresponse dn\\b", "response (down)", x, ignore.case = TRUE)
  x
}

# ---- themes ---------------------------------------------------------------------------------------
theme_pub <- function(base = BASE) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(
      axis.line          = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks         = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length  = unit(1.2, "pt"),
      axis.text          = element_text(size = base - 1, colour = "black"),
      axis.title         = element_text(size = base, colour = "black"),
      legend.title       = element_text(size = base - 1, colour = "black"),
      legend.text        = element_text(size = base - 1, colour = "black"),
      legend.key.size    = unit(2.6, "mm"),
      legend.margin      = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "pt"),
      legend.background  = element_blank(),
      strip.background   = element_blank(),
      strip.text         = element_text(size = base, colour = "black", margin = margin(1, 1, 2, 1)),
      plot.title         = element_blank(),
      plot.subtitle      = element_blank(),
      plot.caption       = element_blank(),
      plot.margin        = margin(2, 2, 2, 2),
      plot.tag           = element_text(size = base + 1, face = "bold", family = FONT,
                                        hjust = 0, vjust = 1)
    )
}
theme_pub_grid <- function(base = BASE) {   # for heatmap/tile panels, which want a frame not axes
  theme_pub(base) +
    theme(axis.line = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3))
}
theme_void_pub <- function(base = BASE) {
  theme_pub(base) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank(), axis.title = element_blank())
}

raster_pts <- function(g, dpi = 600) {
  if (requireNamespace("ggrastr", quietly = TRUE)) ggrastr::rasterise(g, dpi = dpi) else g
}

# ---- output ---------------------------------------------------------------------------------------
# PDF into figures/, nothing else. `tiff = TRUE` additionally writes the submission raster into
# results/figure_exports/. `keep` serializes the ggplot for a later assemble step, also outside figures/.
save_fig <- function(p, name, w, h, dir, tiff = FALSE, keep_dir = NULL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, units = "mm",
         device = cairo_pdf)
  if (isTRUE(tiff)) {
    dir.create(EXPORT_DIR, recursive = TRUE, showWarnings = FALSE)
    f <- file.path(EXPORT_DIR, paste0(name, ".tiff"))
    unlink(f)
    ggsave(f, p, width = w, height = h, units = "mm", dpi = 600, bg = "white", compression = "lzw")
  }
  if (!is.null(keep_dir)) {
    dir.create(keep_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(plot = p, width = w, height = h), file.path(keep_dir, paste0(name, ".rds")))
  }
  message("wrote ", name, " (", w, " x ", h, " mm)")
  invisible(p)
}

# base-graphics panels (WGCNA dendrogram, soft threshold) need the same size and font contract
open_pdf <- function(name, w, h, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  cairo_pdf(file.path(dir, paste0(name, ".pdf")), width = w / 25.4, height = h / 25.4,
            family = FONT, pointsize = BASE)
}

tag_letters <- function(n) if (TAG_CASE == "lower") letters[1:n] else LETTERS[1:n]
