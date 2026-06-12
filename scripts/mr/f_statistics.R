# ===============================================================
# F-statistic summary table for testosterone instruments
# Uses cached instruments from results/instruments.rds
# ===============================================================

library(dplyr)
library(readr)

out_dir <- "results/"

instruments <- readRDS(file.path(out_dir, "instruments.rds"))

# ---- Per-SNP F-stats and R² ---------------------------------------------

per_snp <- bind_rows(lapply(names(instruments), function(lbl) {
  inst <- instruments[[lbl]]
  inst$F_stat <- (inst$beta.exposure / inst$se.exposure)^2
  inst$maf    <- pmin(inst$eaf.exposure, 1 - inst$eaf.exposure)
  # R² per SNP assumes standardised exposure (Ruth used inverse-normalised T)
  inst$R2     <- 2 * inst$maf * (1 - inst$maf) * inst$beta.exposure^2
  tibble(
    exposure  = lbl,
    SNP       = inst$SNP,
    chr       = inst$chr.exposure,
    pos       = inst$pos.exposure,
    EA        = inst$effect_allele.exposure,
    OA        = inst$other_allele.exposure,
    EAF       = inst$eaf.exposure,
    beta      = inst$beta.exposure,
    se        = inst$se.exposure,
    p         = inst$pval.exposure,
    F_stat    = inst$F_stat,
    R2        = inst$R2
  )
}))

write_tsv(per_snp, file.path(out_dir, "per_snp_f_stats.tsv"))


# ---- Per-exposure summary -----------------------------------------------

# Sample sizes for the overall F-statistic (Ruth 2020)
exposure_N <- c(
  "Total T (F)"          = 230454,
  "Total T (M)"          = 194453,
  "Bioavailable T (F)"   = 188507,
  "Bioavailable T (M)"   = 178782
)

f_summary <- per_snp %>%
  group_by(exposure) %>%
  summarise(
    n_snps          = n(),
    mean_F          = mean(F_stat),
    median_F        = median(F_stat),
    min_F           = min(F_stat),
    max_F           = max(F_stat),
    n_weak_F_lt_10  = sum(F_stat < 10),
    total_R2_pct    = sum(R2) * 100,
    mean_per_snp_R2_pct = mean(R2) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    N = exposure_N[exposure],
    # Overall F = R²(N-k-1) / (k(1-R²))
    overall_F = (total_R2_pct/100) * (N - n_snps - 1) /
                (n_snps * (1 - total_R2_pct/100))
  ) %>%
  select(exposure, N, n_snps, total_R2_pct, overall_F,
         mean_F, median_F, min_F, max_F, n_weak_F_lt_10)

print(f_summary, width = Inf)
write_tsv(f_summary, file.path(out_dir, "f_stat_summary.tsv"))
