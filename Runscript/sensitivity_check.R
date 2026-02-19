## =========================================================
## run_sensitivity_analysis_diff.R
##
## Goal: 
## 1. Run Senescence vs No-senescence models under different
##    maturity thresholds (10%, 25%, ..., 90%).
## 2. Calculate the "Difference" (Relative % Change) between
##    No-senescence and the baseline Senescence model.
##    Formula: (Sen - NoSen) / Sen
## 3. Visualize how this difference changes with the threshold.
## =========================================================

rm(list = ls())

library(readxl)
library(ggplot2)
library(minpack.lm)
library(Matrix)
library(exactLTRE)
library(dplyr)
library(tidyr)
library(scales) # Needed for pseudo_log_trans

## -------------------------------
## 1. Setup
## -------------------------------
excel_file <- "Jones2014.xls"  # <-- Confirm filename
output_dir <- "output_sensitivity_diff"
thresholds_to_test <- c(0.10, 0.25, 0.50, 0.75, 0.90)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load your functions
# (Ensure these file names match your local files)
source("LuckFunctions.R")
source("senescence_functions.R") 

## -------------------------------
## 2. Main Loop
## -------------------------------
sheets <- excel_sheets(excel_file)
sensitivity_results <- list()

message("Starting sensitivity analysis on ", length(sheets), " species...")

for (sh in sheets) {
  # Optional: Print progress less frequently if many species
  message("  Processing: ", sh)
  
  df_raw <- try(read_excel(excel_file, sheet = sh, skip = 1), silent = TRUE)
  if (inherits(df_raw, "try-error") || !("x" %in% names(df_raw))) next
  df <- as.data.frame(df_raw)
  
  # Loop through thresholds
  for (thresh in thresholds_to_test) {
    
    tryCatch({
      # 1. Prep data with specific maturity threshold
      demog <- prepare_demography_data_from_df(
        dat = df, input_type = "auto", maturity_age = "auto",
        estimate_tail = TRUE, maturity_method = "logistic",
        maturity_prob = thresh, # <--- Changing Threshold
        auto_min_fx = 1e-8
      )
      
      # 2. Build Models (Only Senescence & No-senescence needed for this check)
      mpm_sen <- build_MPM_senescence(demog$ages, demog$sx, demog$fx)
      mpm_no  <- build_MPM_no_senescence(demog$ages, mpm_sen$sx, mpm_sen$fx, demog$maturity_age)
      
      # 3. Mixing & Summary
      mix_sen <- mixing_distro(mpm_sen$A, mpm_sen$F)
      mix_no  <- mixing_distro(mpm_no$A, mpm_no$F)
      
      # We create dummy matrices for the other 2 models just to satisfy the function input,
      # or we calculate manually. Calling compute_summary_table is easiest.
      summ <- compute_summary_table(
        mpm_sen$U, mpm_no$U, mpm_sen$U, mpm_sen$U, 
        mpm_sen$F, mpm_no$F, mpm_sen$F, mpm_sen$F,
        mix_sen, mix_no, mix_sen, mix_sen, repro_var = "Poisson"
      )
      
      # 4. Store Result
      summ$species <- sh
      summ$threshold <- thresh
      summ$actual_maturity_age <- demog$maturity_age 
      
      # Filter only the 2 models we care about
      summ <- summ %>% filter(model %in% c("Senescence", "No-senescence"))
      
      sensitivity_results[[length(sensitivity_results) + 1]] <- summ
      
    }, error = function(e) {
      # Skip error
    })
  }
}

## -------------------------------
## 3. Calculate Differences
## -------------------------------
all_res <- do.call(rbind, sensitivity_results)

# Pivot wide to calculate difference
# UPDATED FORMULA: (Sen - NoSen) / Sen
df_diff <- all_res %>%
  select(species, threshold, model, 
         mean_lifespan, var_lifespan, skew_lifespan,
         mean_LRO, var_LRO, skew_LRO) %>%
  pivot_wider(
    names_from = model,
    values_from = c(mean_lifespan, var_lifespan, skew_lifespan, mean_LRO, var_LRO, skew_LRO)
  ) %>%
  mutate(
    # NOTE: Modified to (Sen - NoSen) / Sen
    Diff_Mean_Lifespan = (mean_lifespan_Senescence - `mean_lifespan_No-senescence`) / mean_lifespan_Senescence,
    Diff_Var_Lifespan  = (var_lifespan_Senescence  - `var_lifespan_No-senescence`)  / var_lifespan_Senescence,
    Diff_Skew_Lifespan = (skew_lifespan_Senescence - `skew_lifespan_No-senescence`) / skew_lifespan_Senescence,
    
    Diff_Mean_LRO      = (mean_LRO_Senescence      - `mean_LRO_No-senescence`)      / mean_LRO_Senescence,
    Diff_Var_LRO       = (var_LRO_Senescence       - `var_LRO_No-senescence`)       / var_LRO_Senescence,
    Diff_Skew_LRO      = (skew_LRO_Senescence      - `skew_LRO_No-senescence`)      / skew_LRO_Senescence
  )

# Save raw calculations
write.csv(df_diff, file.path(output_dir, "sensitivity_differences_by_species.csv"), row.names = FALSE)

## -------------------------------
## 4. Visualize (Boxplots across Thresholds)
## -------------------------------

# Reshape to long format for Faceted Plotting
df_plot <- df_diff %>%
  select(species, threshold, starts_with("Diff_")) %>%
  pivot_longer(
    cols = starts_with("Diff_"),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric_Label = case_when(
      Metric == "Diff_Mean_Lifespan" ~ "Mean Lifespan",
      Metric == "Diff_Var_Lifespan"  ~ "Variance Lifespan",
      Metric == "Diff_Skew_Lifespan" ~ "Skewness Lifespan",
      Metric == "Diff_Mean_LRO"      ~ "Mean LRO",
      Metric == "Diff_Var_LRO"       ~ "Variance LRO",
      Metric == "Diff_Skew_LRO"      ~ "Skewness LRO"
    ),
    Metric_Type = ifelse(grepl("Lifespan", Metric), "Lifespan Metrics", "LRO Metrics")
  )

# Plotting Function
plot_sensitivity <- function(data, title_suffix) {
  ggplot(data, aes(x = factor(threshold), y = Value, fill = factor(threshold))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_boxplot(outlier.shape = 1, outlier.alpha = 0.3) +
    facet_wrap(~Metric_Label, scales = "free_y", ncol = 3) +
    # Use pseudo_log because Variance differences can be huge
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma = 0.1),
      breaks = c(-10000, -100, -10, -1, 0, 1),
      labels = function(x) paste0(x * 100, "%")
    ) +
    scale_fill_brewer(palette = "Blues") +
    theme_bw(base_size = 14) +
    labs(
      title = paste("Sensitivity to Maturity Threshold:", title_suffix),
      subtitle = "Y-axis shows Relative Difference: (Sen - NoSen)/Sen", # Updated Label
      x = "Maturity Probability Threshold (Logistic)",
      y = "Relative Difference (Pseudo-Log Scale)",
      fill = "Threshold"
    ) +
    theme(legend.position = "none")
}

# 1. Plot Lifespan Metrics
p_life <- plot_sensitivity(
  df_plot %>% filter(Metric_Type == "Lifespan Metrics"), 
  "Lifespan"
)
ggsave(file.path(output_dir, "Sensitivity_Lifespan_Diff.png"), p_life, width = 12, height = 6)

# 2. Plot LRO Metrics
p_lro <- plot_sensitivity(
  df_plot %>% filter(Metric_Type == "LRO Metrics"), 
  "LRO"
)
ggsave(file.path(output_dir, "Sensitivity_LRO_Diff.png"), p_lro, width = 12, height = 6)

message("Sensitivity plots saved to: ", output_dir)