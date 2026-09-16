# mr_functions.R — shared two-sample MR estimators and helpers (scripts 19, 20).
# Estimators are implemented here rather than taken from a package so that every step is visible; they are
# cross-checked against the MendelianRandomization package in script 19 (identical point estimates).

suppressPackageStartupMessages({ library(data.table) })

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
    if (!is.finite(bw) || bw <= 0) return(median(r))
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

# distance pruning: keep the lowest-P SNP per window (fallback when no LD reference is available)
prune_distance <- function(x, kb) {
  x <- x[order(exp_p)]
  keep <- logical(nrow(x))
  for (i in seq_len(nrow(x))) {
    if (!any(keep & abs(x$pos[keep] - x$pos[i]) < kb & x$chr[keep] == x$chr[i])) keep[i] <- TRUE
  }
  x[keep]
}

# LD clumping against 1000 Genomes EUR through the Ensembl REST API (no token required).
# Returns the clumped table plus an attribute recording SNPs Ensembl could not resolve.
ld_proxies <- function(rsid, r2 = 0.1, window_kb = 500, retries = 3) {
  url <- sprintf("https://rest.ensembl.org/ld/human/%s/1000GENOMES:phase_3:EUR?window_size=%d;r2=%s",
                 rsid, window_kb, format(r2, scientific = FALSE))
  for (k in seq_len(retries)) {
    res <- tryCatch(jsonlite::fromJSON(httr::content(
      httr::GET(url, httr::add_headers(`Content-Type` = "application/json"), httr::timeout(60)), "text",
      encoding = "UTF-8")), error = function(e) NULL)
    if (!is.null(res)) {
      if (is.data.frame(res) && nrow(res)) return(data.table(rsid = res$variation2, r2 = as.numeric(res$r2)))
      return(data.table(rsid = character(), r2 = numeric()))
    }
    Sys.sleep(1 + k)
  }
  NULL
}

ld_clump <- function(x, r2 = 0.1, window_kb = 500, sleep = 0.12) {
  x <- x[order(exp_p)]
  keep <- character(0); unresolved <- character(0); pool <- x$rsid
  while (length(pool)) {
    lead <- pool[1]; keep <- c(keep, lead); pool <- pool[-1]
    pr <- ld_proxies(lead, r2 = r2, window_kb = window_kb)
    Sys.sleep(sleep)
    if (is.null(pr)) { unresolved <- c(unresolved, lead); next }
    drop <- pr[r2 >= r2_threshold_env(r2), rsid]
    pool <- setdiff(pool, drop)
  }
  out <- x[rsid %in% keep]
  attr(out, "unresolved") <- unresolved
  out
}
r2_threshold_env <- function(r2) r2

# variance explained by a binary-outcome association (used for Steiger directionality)
get_r2_bin <- function(b, se, ncase, nctrl) {
  k <- ncase / (ncase + nctrl); n <- ncase + nctrl
  z <- b / se; r2 <- z^2 / (z^2 + n)
  r2 * (k * (1 - k)) / (dnorm(qnorm(k))^2)
}
