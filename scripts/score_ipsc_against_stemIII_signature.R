#!/usr/bin/env Rscript
# score_ipsc_against_stemIII_signature.R
# Apply the stemIII-only LASSO signature (K=6 A5SS events) to the iPSC rMATS
# output. Parallel to score_ipsc_against_insT_signature.R.
#
# mu_ctrl[j] recomputed on iPSC KOLF controls per DIV.
# direction[j], beta[j], intercept come from the stemIII export unchanged.
#
# Inputs:
#   data/splicing_signature/insT_stemIII_export/stemIII_signature_export.csv
#   data/splicing_signature/insT_stemIII_export/stemIII_signature_export_meta.txt
#   data/splicing/div{1,5,9,14}/A5SS.MATS.JC.txt
# Outputs:
#   outputs/ipsc_stemIII_signature_scores.rds
#   stdout: headline metrics

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(purrr)
})

SIG_CSV <- "data/splicing_signature/insT_stemIII_export/stemIII_signature_export.csv"
META    <- "data/splicing_signature/insT_stemIII_export/stemIII_signature_export_meta.txt"
OUT_RDS <- "outputs/ipsc_stemIII_signature_scores.rds"
DIVS    <- c("div1", "div5", "div9", "div14")

RENU_LABELS <- c("ReNU8", "ReNU11", "ReNU12")
KOLF_LABELS <- c("KOLF1", "KOLF2",  "KOLF3")

sigmoid  <- function(z) 1 / (1 + exp(-z))
bb_logit <- function(ijc, sjc) {
  p <- (ijc + 0.5) / (ijc + sjc + 1)
  log(p / (1 - p))
}
mw_auc   <- function(case, ctrl) {
  case <- case[!is.na(case)]; ctrl <- ctrl[!is.na(ctrl)]
  r    <- rank(c(case, ctrl))
  (sum(r[seq_along(case)]) - length(case) * (length(case) + 1L) / 2L) /
    (length(case) * length(ctrl))
}

stopifnot("stemIII signature CSV not found"  = file.exists(SIG_CSV),
          "stemIII signature meta not found" = file.exists(META))
sig    <- read_csv(SIG_CSV, show_col_types = FALSE)
K_FULL <- nrow(sig)

# Intercept read directly from meta.txt -- never hardcode LOO output values.
INTERCEPT_1SE <- as.numeric(sub(".*: *", "",
                           grep("^intercept_at_1se", readLines(META), value = TRUE)))
stopifnot("Could not parse intercept_at_1se from meta.txt" =
            length(INTERCEPT_1SE) == 1L && is.finite(INTERCEPT_1SE))

cat(sprintf("stemIII signature loaded: %d events; sum(beta) = %.3f\n",
            nrow(sig), sum(sig$beta_at_lambda1se)))

sig <- sig %>%
  mutate(coord_key = paste(chr,
                           longExonStart_0base, longExonEnd_1based,
                           shortES_0base,       shortEE_1based,
                           flankingES_0base,    flankingEE_1based,
                           sep = ":"))
stopifnot(!anyDuplicated(sig$coord_key))

parse_vec <- function(s) do.call(rbind, lapply(strsplit(s, ","), as.integer))

load_div <- function(div) {
  path <- file.path("data/splicing", div, "A5SS.MATS.JC.txt")
  stopifnot(file.exists(path))
  d <- read_tsv(path, show_col_types = FALSE,
                col_types = cols(.default = col_character()))
  d <- d %>%
    mutate(coord_key = paste(chr, longExonStart_0base, longExonEnd,
                             shortES, shortEE, flankingES, flankingEE, sep = ":"))
  d    <- d[d$coord_key %in% sig$coord_key, , drop = FALSE]
  d    <- d[match(sig$coord_key, d$coord_key), , drop = FALSE]
  miss <- which(is.na(d$coord_key))
  if (length(miss) > 0L)
    warning(sprintf("[%s] missing %d stemIII events: %s",
                    div, length(miss),
                    paste(sig$gene_symbol[miss], collapse = ", ")))
  ijc1 <- parse_vec(d$IJC_SAMPLE_1)
  sjc1 <- parse_vec(d$SJC_SAMPLE_1)
  ijc2 <- parse_vec(d$IJC_SAMPLE_2)
  sjc2 <- parse_vec(d$SJC_SAMPLE_2)
  stopifnot(ncol(ijc1) == 3L, ncol(ijc2) == 3L, nrow(ijc1) == nrow(sig))
  list(ijc1 = ijc1, sjc1 = sjc1, ijc2 = ijc2, sjc2 = sjc2,
       hits_ok = !is.na(d$coord_key))
}

per_div <- setNames(lapply(DIVS, load_div), DIVS)

cov_ok <- sapply(per_div, function(x) sum(x$hits_ok))
cat(sprintf("Signature events recovered per DIV: %s\n",
            paste(sprintf("%s=%d/%d", names(cov_ok), cov_ok, nrow(sig)),
                  collapse = ", ")))

keep_global <- Reduce(`&`, lapply(per_div, `[[`, "hits_ok"))
n_keep <- sum(keep_global)
if (n_keep < nrow(sig)) {
  cat(sprintf("Dropping %d events not quantified in iPSC: %s\n",
              nrow(sig) - n_keep,
              paste(sig$gene_symbol[!keep_global], collapse = ", ")))
  sig     <- sig[keep_global, , drop = FALSE]
  per_div <- lapply(per_div, function(pd) {
    list(ijc1    = pd$ijc1[keep_global, , drop = FALSE],
         sjc1    = pd$sjc1[keep_global, , drop = FALSE],
         ijc2    = pd$ijc2[keep_global, , drop = FALSE],
         sjc2    = pd$sjc2[keep_global, , drop = FALSE],
         hits_ok = pd$hits_ok[keep_global])
  })
}
N_EVT_USED <- nrow(sig)

make_logit_block <- function(div, ijc, sjc, labels, group) {
  m <- t(bb_logit(ijc, sjc))
  rownames(m) <- paste0(labels, "_", div)
  list(mat  = m,
       meta = data.frame(sample = rownames(m),
                         line   = labels,
                         DIV    = div,
                         group  = group,
                         stringsAsFactors = FALSE))
}

blocks <- list()
for (div in DIVS) {
  pd <- per_div[[div]]
  blocks[[length(blocks) + 1L]] <- make_logit_block(div, pd$ijc1, pd$sjc1,
                                                     RENU_LABELS, "ReNU")
  blocks[[length(blocks) + 1L]] <- make_logit_block(div, pd$ijc2, pd$sjc2,
                                                     KOLF_LABELS, "KOLF")
}
logit_psi   <- do.call(rbind, lapply(blocks, `[[`, "mat"))
colnames(logit_psi) <- sig$gene_symbol
sample_meta <- do.call(rbind, lapply(blocks, `[[`, "meta"))
rownames(sample_meta) <- sample_meta$sample
stopifnot(nrow(logit_psi) == 24L,
          ncol(logit_psi) == N_EVT_USED,
          identical(rownames(logit_psi), sample_meta$sample))

is_kolf <- sample_meta$group == "KOLF"
mu_ctrl <- matrix(NA_real_, nrow = length(DIVS), ncol = ncol(logit_psi),
                  dimnames = list(DIVS, colnames(logit_psi)))
for (d in DIVS) {
  rows_d       <- is_kolf & sample_meta$DIV == d
  mu_ctrl[d, ] <- colMeans(logit_psi[rows_d, , drop = FALSE], na.rm = TRUE)
}
stopifnot(!anyNA(mu_ctrl))

mu_per_sample <- mu_ctrl[sample_meta$DIV, , drop = FALSE]
oriented      <- (logit_psi - mu_per_sample)
oriented      <- sweep(oriented, 2, sig$direction, FUN = "*")
uniform_sum   <- rowSums(oriented)
lasso_linpred <- INTERCEPT_1SE + as.numeric(oriented %*% sig$beta_at_lambda1se)
lasso_prob    <- sigmoid(lasso_linpred)

sample_meta$uniform_sum   <- uniform_sum
sample_meta$lasso_linpred <- lasso_linpred
sample_meta$lasso_prob    <- lasso_prob

case_us <- sample_meta$uniform_sum[sample_meta$group == "ReNU"]
ctrl_us <- sample_meta$uniform_sum[sample_meta$group == "KOLF"]
auc_us  <- mw_auc(case_us, ctrl_us)
mw_us   <- wilcox.test(case_us, ctrl_us, alternative = "greater", exact = FALSE)

case_lp <- sample_meta$lasso_prob[sample_meta$group == "ReNU"]
ctrl_lp <- sample_meta$lasso_prob[sample_meta$group == "KOLF"]
auc_lp  <- mw_auc(case_lp, ctrl_lp)
mw_lp   <- wilcox.test(case_lp, ctrl_lp, alternative = "greater", exact = FALSE)

case_mean <- colMeans(logit_psi[!is_kolf, , drop = FALSE], na.rm = TRUE)
ctrl_mean <- colMeans(logit_psi[ is_kolf, , drop = FALSE], na.rm = TRUE)
delta     <- case_mean - ctrl_mean
concord   <- sign(delta) == sign(sig$direction)
n_concord <- sum(concord)
binom     <- binom.test(n_concord, length(concord), p = 0.5,
                        alternative = "greater")

res <- list(
  signature            = sig,
  mu_ctrl_ipsc_per_div = mu_ctrl,
  logit_psi            = logit_psi,
  oriented             = oriented,
  per_sample           = sample_meta,
  per_event            = data.frame(
    rank            = sig$rank,
    gene            = sig$gene_symbol,
    direction       = sig$direction,
    beta_1se        = sig$beta_at_lambda1se,
    ctrl_mean_logit = ctrl_mean,
    case_mean_logit = case_mean,
    delta_logit     = delta,
    concordant      = concord,
    stringsAsFactors = FALSE),
  metrics = list(
    auc_uniform  = auc_us,  mw_uniform_p = mw_us$p.value,
    auc_lasso    = auc_lp,  mw_lasso_p   = mw_lp$p.value,
    n_concord    = n_concord, n_events   = length(concord),
    binom_p      = binom$p.value, binom_ci = binom$conf.int
  )
)
dir.create("outputs", showWarnings = FALSE)
saveRDS(res, OUT_RDS)

cat("\n========================================================\n")
cat("iPSC stemIII-specific LASSO signature application\n")
cat("========================================================\n")
cat(sprintf("n samples            : %d (12 ReNU vs 12 KOLF; 3 lines x 4 DIV)\n",
            nrow(sample_meta)))
cat(sprintf("K events used        : %d / %d\n", N_EVT_USED, K_FULL))
cat(sprintf("Uniform-sum score    : AUC = %.3f   MWU one-sided p = %.3g\n",
            auc_us, mw_us$p.value))
cat(sprintf("LASSO P(case)        : AUC = %.3f   MWU one-sided p = %.3g\n",
            auc_lp, mw_lp$p.value))
cat(sprintf("Direction concordance: %d/%d events  binomial one-sided p = %.3g\n",
            n_concord, length(concord), binom$p.value))
cat("\nPer-event delta-logit-PSI (ReNU mean - KOLF mean):\n")
print(res$per_event[, c("rank", "gene", "direction", "delta_logit", "concordant")],
      row.names = FALSE)
cat(sprintf("\nWritten: %s\n", OUT_RDS))
