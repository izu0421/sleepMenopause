# ===============================================================
# Mendelian Randomization: Testosterone -> Sleep latency / insomnia
#
# Exposures (Ruth et al. 2020, UKB):
#   - Total testosterone, females  (ebi-a-GCST90012112, N=230,454)
#   - Total testosterone, males    (ebi-a-GCST90012113, N=194,453)
#   - Bioavailable T, females      (ebi-a-GCST90012102, N=188,507)
#   - Bioavailable T, males        (ebi-a-GCST90012103, N=178,782)
#
# Primary outcome:
#   - Amin et al. 2016 sleep latency GWAS (non-UKB, no overlap)
#     GWAS Central HGVST1836 — download separately and set path below.
#
# Secondary outcome:
#   - ukb-b-3957 (Elsworth, sleeplessness/insomnia, combined-sex)
#     Sample overlap with exposures -> use MRlap correction.
# ===============================================================

# ---- Setup ---------------------------------------------------------------

# install.packages(c("remotes", "dplyr", "readr", "tidyr"))
# remotes::install_github("MRCIEU/TwoSampleMR")
# remotes::install_github("MRCIEU/ieugwasr")
# remotes::install_github("rondolab/MR-PRESSO")
# remotes::install_github("n-mounier/MRlap")
# install.packages("LDlinkR")

library(TwoSampleMR)
library(ieugwasr)
library(MRPRESSO)
library(MRlap)
library(LDlinkR)
library(dplyr)
library(readr)
library(tidyr)


# ===============================================================
# MR: Testosterone -> Insomnia symptom (ukb-b-3957)
#
# Exposures (Ruth et al. 2020, UKB, sex-stratified):
#   - Total testosterone, females  (ebi-a-GCST90012112, N=230,454)
#   - Total testosterone, males    (ebi-a-GCST90012113, N=194,453)
#   - Bioavailable T, females      (ebi-a-GCST90012102, N=188,507)
#   - Bioavailable T, males        (ebi-a-GCST90012103, N=178,782)
#
# Outcome:
#   - Sleeplessness / insomnia (Elsworth ukb-b-3957, combined-sex)
#
# NOTE: All analyses have sample overlap (UKB on both sides).
#       Bias direction with strong instruments is toward observational.
#       Report this as a limitation.
# ===============================================================

library(TwoSampleMR)
library(dplyr)
library(readr)

# Optional: MR-PRESSO if installed
have_presso <- requireNamespace("MRPRESSO", quietly = TRUE)
if (have_presso) library(MRPRESSO)

# OpenGWAS token: https://api.opengwas.io/profile/

# ===============================================================
# MR: Testosterone -> Insomnia (ukb-b-3957)
# Resilient version: per-exposure save, API retries, no PRESSO.
# ===============================================================

library(TwoSampleMR)
library(dplyr)
library(readr)

# Paste a FRESH token (https://api.opengwas.io/profile/)
Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJzaGluMkBjYW0uYWMudWsiLCJpYXQiOjE3NzkzOTk2MTYsImV4cCI6MTc4MDYwOTIxNn0.KfKlde4LmJdYc-JkisiFEZn4deJZknK5vIu5ik2gz3ujDYEYUWgGeOYQNoNpcxJwrGr7TGdPV41KXzTdF1oYi4y-trqKwvrSNRbLNFBMgiPYG2WdcWzqnvo30AOwbIURIUnmWRjs0hxRqL3JpKcdZNDRPKU6uFB4pHjN2KpEZyzrbgt0vI2-BExRifdQd1sTNw16tQTZEE2wfkbBsns2H0p2O7ECjqB51mPMwnBKQ3hO_qKzBmMUwd3HpymITPmmKntKVkK5DcaTZicar57qTLCF1qgY-9iFCMuImNgRkC6cF0cX3-rHaxs1AvZpJ7RCizUiyCZYZBLKpuOpvrjoPQ")

out_dir <- "results/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---- Exposures and outcome -----------------------------------------------

exposures <- tibble::tribble(
  ~id,                       ~label,
  "ebi-a-GCST90012112",      "Total T (F)",
  "ebi-a-GCST90012113",      "Total T (M)",
  "ebi-a-GCST90012102",      "Bioavailable T (F)",
  "ebi-a-GCST90012103",      "Bioavailable T (M)"
)

outcome_id    <- "ukb-b-3957"
outcome_label <- "Insomnia (ukb-b-3957)"


# ---- 1. Load or extract instruments (cached) -----------------------------

inst_cache <- file.path(out_dir, "instruments.rds")

if (file.exists(inst_cache)) {
  instruments <- readRDS(inst_cache)
  message("Loaded cached instruments from ", inst_cache)
} else {
  extract_qc <- function(id, label) {
    message("--- extracting ", label, " ---")
    inst <- extract_instruments(
      outcomes = id, p1 = 5e-8,
      clump = TRUE, r2 = 0.001, kb = 10000
    )
    inst$F_stat   <- (inst$beta.exposure / inst$se.exposure)^2
    inst$exposure <- label
    message("  ", nrow(inst), " SNPs | min F = ", round(min(inst$F_stat), 1))
    inst
  }
  instruments <- Map(extract_qc, exposures$id, exposures$label)
  names(instruments) <- exposures$label
  saveRDS(instruments, inst_cache)
}


# ---- 2. Resilient per-exposure runner ------------------------------------

run_one_robust <- function(exp_dat, exp_label, outcome_id, outcome_label,
                           max_retries = 3) {
  
  safe_lbl    <- gsub("[^A-Za-z0-9]", "_", exp_label)
  result_file <- file.path(out_dir, paste0("result_", safe_lbl, ".rds"))
  
  if (file.exists(result_file)) {
    message("  Already done: ", result_file)
    return(readRDS(result_file))
  }
  
  # Outcome extraction with retry/backoff
  out_dat <- NULL
  for (i in seq_len(max_retries)) {
    out_dat <- tryCatch(
      extract_outcome_data(snps = exp_dat$SNP, outcomes = outcome_id),
      error = function(e) {
        message("  Attempt ", i, " failed: ", e$message)
        Sys.sleep(30 * i)
        NULL
      }
    )
    if (!is.null(out_dat)) break
  }
  if (is.null(out_dat) || nrow(out_dat) == 0) {
    warning("Outcome extraction failed for ", exp_label)
    return(NULL)
  }
  out_dat$outcome <- outcome_label
  
  harm <- harmonise_data(exp_dat, out_dat, action = 2)
  message("  Harmonised SNPs kept: ", sum(harm$mr_keep), "/", nrow(harm))
  if (sum(harm$mr_keep) < 3) {
    warning("Too few SNPs for MR (<3) in ", exp_label)
    return(list(harm = harm))
  }
  
  mr_res <- mr(harm, method_list = c(
    "mr_ivw", "mr_egger_regression",
    "mr_weighted_median", "mr_weighted_mode"
  ))
  
  result <- list(
    exposure_label = exp_label,
    harm    = harm,
    mr      = mr_res,
    egger   = mr_pleiotropy_test(harm),
    het     = mr_heterogeneity(harm),
    steiger = tryCatch(directionality_test(harm), error = function(e) NULL),
    loo     = mr_leaveoneout(harm),
    sing    = mr_singlesnp(harm)
  )
  
  saveRDS(result, result_file)
  message("  Saved: ", result_file)
  result
}


# ---- 3. Run each, catching per-exposure failures -------------------------

results <- list()
for (lbl in names(instruments)) {
  message("\n=== ", lbl, " ===")
  results[[lbl]] <- tryCatch(
    run_one_robust(instruments[[lbl]], lbl, outcome_id, outcome_label),
    error = function(e) {
      message("  FAILED: ", e$message)
      NULL
    }
  )
}


# ---- 4. Compile results --------------------------------------------------

flatten <- function(field) {
  bind_rows(lapply(names(results), function(lbl) {
    r <- results[[lbl]]
    if (is.null(r) || is.null(r[[field]])) return(NULL)
    cbind(exposure_label = lbl, r[[field]])
  }))
}

all_mr    <- flatten("mr")
all_egger <- flatten("egger")
all_het   <- flatten("het")

all_mr$bonf_sig <- all_mr$method == "Inverse variance weighted" &
  all_mr$pval < (0.05 / 4)

print(all_mr %>%
        filter(method == "Inverse variance weighted") %>%
        select(exposure_label, nsnp, b, se, pval, bonf_sig))

write_tsv(all_mr,    file.path(out_dir, "mr_results.tsv"))
write_tsv(all_egger, file.path(out_dir, "egger_intercepts.tsv"))
write_tsv(all_het,   file.path(out_dir, "heterogeneity.tsv"))


# ---- 5. Plots ------------------------------------------------------------

for (lbl in names(results)) {
  r <- results[[lbl]]
  if (is.null(r) || is.null(r$mr)) next
  safe <- gsub("[^A-Za-z0-9]", "_", lbl)
  
  pdf(file.path(out_dir, paste0("scatter_", safe, ".pdf")), 6, 6)
  print(mr_scatter_plot(r$mr, r$harm)); dev.off()
  
  pdf(file.path(out_dir, paste0("funnel_",  safe, ".pdf")), 6, 6)
  print(mr_funnel_plot(r$sing)); dev.off()
}

save(results, instruments, file = file.path(out_dir, "workspace.RData"))
message("\nDone. Results in: ", normalizePath(out_dir))


### PRESSO #####


library(TwoSampleMR)
library(MRPRESSO)
library(dplyr)
library(readr)

out_dir <- "results/"
result_files <- list.files(out_dir, pattern = "^result_.*\\.rds$", full.names = TRUE)

run_presso_one <- function(result_rds, NbDistribution = 10000) {
  r        <- readRDS(result_rds)
  safe     <- gsub("[^A-Za-z0-9]", "_", r$exposure_label)
  out_file <- file.path(out_dir, paste0("presso_", safe, ".rds"))
  
  # Skip if already done
  if (file.exists(out_file)) {
    message("Already done: ", out_file)
    return(readRDS(out_file))
  }
  
  harm_filt <- r$harm[r$harm$mr_keep, ]
  message("\n=== ", r$exposure_label, " | ", nrow(harm_filt), " SNPs ===")
  message("Start: ", Sys.time())
  t0 <- Sys.time()
  
  presso <- tryCatch(
    MRPRESSO::mr_presso(
      BetaOutcome     = "beta.outcome",
      BetaExposure    = "beta.exposure",
      SdOutcome       = "se.outcome",
      SdExposure      = "se.exposure",
      OUTLIERtest     = TRUE,
      DISTORTIONtest  = TRUE,
      data            = harm_filt,
      NbDistribution  = NbDistribution,
      SignifThreshold = 0.05
    ),
    error = function(e) { message("FAILED: ", e$message); NULL }
  )
  
  message("End:   ", Sys.time(), " (elapsed: ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min)")
  
  saved <- list(
    exposure       = r$exposure_label,
    nsnps          = nrow(harm_filt),
    NbDistribution = NbDistribution,
    presso         = presso,
    harm_filt      = harm_filt,
    snp_order      = harm_filt$SNP  # so you can map outlier indices back to rsids
  )
  saveRDS(saved, out_file)
  message("Saved: ", out_file)
  saved
}

# Run sequentially
for (f in result_files) {
  run_presso_one(f, NbDistribution = 10000)
}

# Compile summary
presso_files <- list.files(out_dir, pattern = "^presso_.*\\.rds$", full.names = TRUE)

summary_table <- bind_rows(lapply(presso_files, function(f) {
  s <- readRDS(f)
  if (is.null(s$presso)) return(NULL)
  
  main      <- s$presso$`Main MR results`
  global_p  <- s$presso$`MR-PRESSO results`$`Global Test`$Pvalue
  dist_p    <- s$presso$`MR-PRESSO results`$`Distortion Test`$Pvalue
  out_idx   <- s$presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
  outliers  <- if (length(out_idx) > 0 && !identical(out_idx, "No significant outliers"))
    s$snp_order[out_idx] else character(0)
  
  tibble(
    exposure            = s$exposure,
    nsnps               = s$nsnps,
    raw_beta            = main$`Causal Estimate`[main$`MR Analysis` == "Raw"],
    raw_se              = main$Sd[main$`MR Analysis` == "Raw"],
    raw_p               = main$`P-value`[main$`MR Analysis` == "Raw"],
    corrected_beta      = main$`Causal Estimate`[main$`MR Analysis` == "Outlier-corrected"],
    corrected_se        = main$Sd[main$`MR Analysis` == "Outlier-corrected"],
    corrected_p         = main$`P-value`[main$`MR Analysis` == "Outlier-corrected"],
    global_p_pleiotropy = global_p,
    distortion_p        = dist_p,
    n_outliers          = length(outliers),
    outlier_snps        = paste(outliers, collapse = ";")
  )
}))

print(summary_table)
write_tsv(summary_table, file.path(out_dir, "mr_presso_summary.tsv"))


