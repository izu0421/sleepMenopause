# ===============================================================
# Final female-only figures
#
#   - FDR computed across all 10 primary tests (2 exposures × 5 methods)
#   - Weighted mode dropped throughout
#   - Methods retained: IVW, MR-Egger, Weighted median, PRESSO (raw),
#                       PRESSO (corrected)
#   - Outputs:
#       Individual plots in plots_female_final/
#       4-panel composite (forests top, scatters bottom)
#       F-statistic table
# ===============================================================
setwd("C:/Users/suley/OneDrive/Documents/")
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

out_dir  <- "results/"
plot_dir <- "results/plots_female_final/"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

female_labels <- c("Total T (F)", "Bioavailable T (F)")
exposure_N    <- c("Total T (F)" = 230454, "Bioavailable T (F)" = 188507)

method_levels <- c("IVW", "MR-Egger", "Weighted median",
                   "PRESSO (raw)", "PRESSO (corrected)")
method_colors <- c(
  "IVW"                 = "#1f77b4",
  "MR-Egger"            = "#ff7f0e",
  "Weighted median"     = "#2ca02c",
  "PRESSO (raw)"        = "#7f7f7f",
  "PRESSO (corrected)"  = "#d62728"
)


# ---- 1. Load saved objects ----------------------------------------------

result_files <- list.files(out_dir, pattern = "^result_.*\\.rds$", full.names = TRUE)
results <- lapply(result_files, readRDS)
names(results) <- sapply(results, function(r) r$exposure_label)
results <- results[female_labels]

presso_files <- list.files(out_dir, pattern = "^presso_.*\\.rds$", full.names = TRUE)
presso <- lapply(presso_files, readRDS)
names(presso) <- sapply(presso, function(r) r$exposure)
presso <- presso[female_labels]


# ---- 2. Build MR results table (5 methods × 2 exposures) ----------------

std <- bind_rows(lapply(results, function(r) {
  r$mr %>%
    filter(method %in% c("Inverse variance weighted",
                         "MR Egger",
                         "Weighted median")) %>%
    mutate(exposure_label = r$exposure_label,
           method_clean = case_when(
             method == "Inverse variance weighted" ~ "IVW",
             method == "MR Egger"                  ~ "MR-Egger",
             method == "Weighted median"           ~ "Weighted median"
           ))
})) %>%
  select(exposure_label, method_clean, b, se, pval, nsnp)

pr <- bind_rows(lapply(presso, function(s) {
  m     <- s$presso$`Main MR results`
  n_out <- length(s$presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
  bind_rows(
    tibble(exposure_label = s$exposure,
           method_clean   = "PRESSO (raw)",
           b              = m$`Causal Estimate`[m$`MR Analysis` == "Raw"],
           se             = m$Sd[m$`MR Analysis` == "Raw"],
           pval           = m$`P-value`[m$`MR Analysis` == "Raw"],
           nsnp           = s$nsnps),
    tibble(exposure_label = s$exposure,
           method_clean   = "PRESSO (corrected)",
           b              = m$`Causal Estimate`[m$`MR Analysis` == "Outlier-corrected"],
           se             = m$Sd[m$`MR Analysis` == "Outlier-corrected"],
           pval           = m$`P-value`[m$`MR Analysis` == "Outlier-corrected"],
           nsnp           = s$nsnps - n_out)
  )
}))

mr_data <- bind_rows(std, pr) %>%
  mutate(ci_lower = b - 1.96 * se,
         ci_upper = b + 1.96 * se,
         # FDR across all 10 tests (both exposures × all 5 methods)
         fdr      = p.adjust(pval, method = "BH"))

mr_data$method_clean   <- factor(mr_data$method_clean, levels = method_levels)
mr_data$exposure_label <- factor(mr_data$exposure_label, levels = female_labels)

write_tsv(mr_data, file.path(plot_dir, "mr_results_female_with_fdr.tsv"))


# ---- 3. F-statistic table -----------------------------------------------

instruments <- readRDS(file.path(out_dir, "instruments.rds"))[female_labels]

fstat_table <- bind_rows(lapply(names(instruments), function(lbl) {
  inst    <- instruments[[lbl]]
  F_stats <- (inst$beta.exposure / inst$se.exposure)^2
  maf     <- pmin(inst$eaf.exposure, 1 - inst$eaf.exposure)
  R2      <- 2 * maf * (1 - maf) * inst$beta.exposure^2
  tot_R2  <- sum(R2)
  k       <- nrow(inst)
  N       <- exposure_N[lbl]
  tibble(
    Exposure     = lbl,
    N            = N,
    n_SNPs       = k,
    Total_R2_pct = tot_R2 * 100,
    Overall_F    = tot_R2 * (N - k - 1) / (k * (1 - tot_R2)),
    Mean_F       = mean(F_stats),
    Median_F     = median(F_stats),
    Min_F        = min(F_stats),
    Max_F        = max(F_stats),
    SNPs_F_lt_10 = sum(F_stats < 10)
  )
}))

print(fstat_table, width = Inf)
write_tsv(fstat_table, file.path(plot_dir, "fstat_table.tsv"))


# ---- 4. Forest plot per exposure ----------------------------------------

make_forest <- function(exposure_name) {
  d <- mr_data %>% filter(exposure_label == exposure_name)

  ggplot(d, aes(x = b, y = method_clean)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                   height = 0.2, color = "grey30") +
    geom_point(aes(fill = method_clean),
               shape = 22, size = 3.5, color = "black") +
    geom_text(aes(label = sprintf("p=%.2g, q=%.2g", pval, fdr),
                  x = ci_upper),
              hjust = -0.12, size = 3) +
    scale_fill_manual(values = method_colors) +
    scale_y_discrete(limits = rev(method_levels)) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.65))) +
    labs(x = expression(beta*" (95% CI)"),
         y = NULL,
         title = paste0(exposure_name, " → Insomnia")) +
    theme_bw(base_size = 11) +
    theme(legend.position    = "none",
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          plot.title         = element_text(face = "bold", size = 11))
}

forest_total <- make_forest("Total T (F)")
forest_bioav <- make_forest("Bioavailable T (F)")


# ---- 5. Scatter plot per exposure ---------------------------------------

make_scatter <- function(exposure_name) {
  r <- results[[exposure_name]]
  harm <- r$harm[r$harm$mr_keep, ]

  # Standard MR convention: align exposure effects to be positive
  flip <- harm$beta.exposure < 0
  harm$beta.outcome[flip]  <- -harm$beta.outcome[flip]
  harm$beta.exposure[flip] <- -harm$beta.exposure[flip]

  # Overlay IVW, MR-Egger, weighted median (drop weighted mode)
  mr_res <- r$mr %>%
    filter(method %in% c("Inverse variance weighted",
                         "MR Egger",
                         "Weighted median"))

  scatter_colors <- c(
    "Inverse variance weighted" = "#1f77b4",
    "MR Egger"                  = "#ff7f0e",
    "Weighted median"           = "#2ca02c"
  )
  scatter_labels <- c(
    "Inverse variance weighted" = "IVW",
    "MR Egger"                  = "MR-Egger",
    "Weighted median"           = "Weighted median"
  )

  ggplot(harm, aes(x = beta.exposure, y = beta.outcome)) +
    geom_hline(yintercept = 0, color = "grey80") +
    geom_vline(xintercept = 0, color = "grey80") +
    geom_errorbar(aes(ymin = beta.outcome - 1.96 * se.outcome,
                      ymax = beta.outcome + 1.96 * se.outcome),
                  width = 0, color = "grey75", alpha = 0.5) +
    geom_errorbarh(aes(xmin = beta.exposure - 1.96 * se.exposure,
                       xmax = beta.exposure + 1.96 * se.exposure),
                   height = 0, color = "grey75", alpha = 0.5) +
    geom_point(size = 1.2, alpha = 0.6, color = "grey30") +
    geom_abline(data = mr_res,
                aes(intercept = 0, slope = b, color = method),
                linewidth = 0.8) +
    scale_color_manual(values = scatter_colors, labels = scatter_labels) +
    labs(title = exposure_name,
         x = "SNP effect on testosterone",
         y = "SNP effect on insomnia",
         color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position      = c(0.02, 0.98),
          legend.justification  = c(0, 1),
          legend.background     = element_rect(fill = alpha("white", 0.7), color = NA),
          legend.key.size       = unit(0.4, "cm"),
          legend.text           = element_text(size = 8),
          plot.title            = element_text(face = "bold", size = 11),
          panel.grid.minor      = element_blank())
}

scatter_total <- make_scatter("Total T (F)")
scatter_bioav <- make_scatter("Bioavailable T (F)")


# ---- 6. Save individual figures -----------------------------------------

save_pair <- function(plot, name, w, h) {
  ggsave(file.path(plot_dir, paste0(name, ".pdf")), plot, width = w, height = h)
  ggsave(file.path(plot_dir, paste0(name, ".png")), plot, width = w, height = h,
         dpi = 300, bg = "white")
}

save_pair(forest_total,  "forest_TotalTF",   8, 4.5)
save_pair(forest_bioav,  "forest_BioavailTF", 8, 4.5)
save_pair(scatter_total, "scatter_TotalTF",   6, 5)
save_pair(scatter_bioav, "scatter_BioavailTF", 6, 5)


# ---- 7. 4-panel composite (forests top, scatters bottom) ----------------

composite <- (forest_total | forest_bioav) /
             (scatter_total | scatter_bioav) +
  plot_annotation(
    tag_levels = 'A',
    title      = "Testosterone → Insomnia (female only)",
    subtitle   = paste0("UK Biobank: testosterone (Ruth 2020) and insomnia symptom ",
                        "(ukb-b-3957). Sample overlap acknowledged."),
    theme      = theme(plot.title    = element_text(face = "bold", size = 14),
                       plot.subtitle = element_text(color = "grey40", size = 10))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 14))

save_pair(composite, "panel_4_female_composite", 14, 10)

message("\nDone. Outputs in: ", normalizePath(plot_dir))
