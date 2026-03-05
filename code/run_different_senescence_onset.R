## =========================================================
## run_peak_vs_maturity_comparison.Rb(To be DONE)
##
## Purpose: 
## Final refined visualization for methodology comparison.
## Fixes: 
## 1. Prevents label clipping (species names now visible outside plot).
## 2. Cleans legend (Class legend shows lines only, no triangles).
## 3. Uses triangles only for plot highlighting of sensitive species.
## =========================================================

rm(list = ls())

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales) 

## -------------------------------
## 1. Setup Environment
## -------------------------------
excel_file <- "data/Jones2014.xls" 
output_dir <- "Results/sensitivity analysis"
species_plots_dir <- file.path("Results/sensitivity analysis/species_diagnostic_plots")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(species_plots_dir)) dir.create(species_plots_dir, recursive = TRUE)

# Ensure helper scripts are present
source("code/LuckFunctions.R")
source("code/senescence_functions.R") 


## -------------------------------
## 2. Main Data Processing
## -------------------------------
sheets <- excel_sheets(excel_file)
comparison_results <- list()

message("Starting analysis loop...")

for (sh in sheets) {
  
  message("\n==============================")
  message("Processing sheet: ", sh)
  message("==============================")
  
  species_name <- sh
  
  meta_raw <- try(
    read_excel(excel_file, sheet = sh, range = "E1", col_names = FALSE),
    silent = TRUE
  )
  study_type <- if (!inherits(meta_raw, "try-error") && nrow(meta_raw) > 0) as.character(meta_raw[[1,1]]) else "Unknown"
  message("  -> Data type detected: ", study_type)
  
  ## ---- Read data ----
  df_raw <- try(
    read_excel(excel_file, sheet = sh, skip = 1, col_names = TRUE),
    silent = TRUE
  )
  
  if (inherits(df_raw, "try-error")) {
    message("  -> ERROR reading sheet: ", sh)
    next
  }
  
  df <- as.data.frame(df_raw)
  
  if (!("x" %in% names(df))) {
    message("  -> SKIP: sheet ", sh, " has no 'x' column.")
    next
  }
  
  ## -------------------------------
  ## Attempt full demographic analysis
  ## -------------------------------
  res <- try({
    ## 1. Prepare demographic inputs
    demog <- prepare_demography_data_from_df(dat = df, input_type = "auto")
    
    # Cap the maximum survival rate to avoid exact 1.0
    demog$sx[demog$sx >= 1] <- 0.9999 
    
    ## Basic validity checks
    if (!any(is.finite(demog$sx))) stop("All sx values are non-finite")
    if (!any(is.finite(demog$fx))) stop("All fx values are non-finite")
    
    ## 2. Build MPMs for four senescence scenarios
    mpm_sen <- build_MPM_senescence(ages = demog$ages, sx = demog$sx, fx = demog$fx)
    mpm_no_peak <- build_MPM_no_senescence(ages = demog$ages, sx_senescence = mpm_sen$sx, fx_senescence = mpm_sen$fx, senescence_onset_age = demog$senescence_onset_age)
    mpm_no_early <- build_MPM_no_senescence(ages = demog$ages, sx_senescence = mpm_sen$sx, fx_senescence = mpm_sen$fx, senescence_onset_age = demog$early_onset)
    mpm_no_late <- build_MPM_no_senescence(ages = demog$ages, sx_senescence = mpm_sen$sx, fx_senescence = mpm_sen$fx, senescence_onset_age = demog$late_onset)
    
    ## 3. Create species vector for plotting
    species_vectors <- list(
      sen      = list(ages = demog$ages, sx = mpm_sen$sx,      fx = mpm_sen$fx),
      no_peak  = list(ages = demog$ages, sx = mpm_no_peak$sx,  fx = mpm_no_peak$fx),
      no_early = list(ages = demog$ages, sx = mpm_no_early$sx, fx = mpm_no_early$fx),
      no_late  = list(ages = demog$ages, sx = mpm_no_late$sx,  fx = mpm_no_late$fx)
    )
    
    ## 4. Summary statistics
    summary_df <- compute_sensitivity_table(
      U_sen        = mpm_sen$U,
      U_no_peak    = mpm_no_peak$U,
      U_no_early   = mpm_no_early$U,
      U_no_late    = mpm_no_late$U,
      F_sen        = mpm_sen$F,
      F_no_peak    = mpm_no_peak$F,
      F_no_early   = mpm_no_early$F,
      F_no_late    = mpm_no_late$F,
      repro_var    = "Poisson"  
    )
    
    summary_df$species <- species_name
    summary_df$sheet   <- sh
    
    summary_df 
  }) #<- end of try()
  
  
  ## 4. Save section (outside try() but inside for() )
  
  if (!inherits(res, "try-error")) {
    message("  -> OK: added species ", species_name)
    
    # 1. Save 'res' into the global list using species_name as the key
    comparison_results[[species_name]] <- res
    
    # 2. Save individual CSV for this species
    write.csv(
      res,
      file = file.path(output_dir, paste0(species_name, "_summary_stats.csv")),
      row.names = FALSE
    )
  } else {
    message("  -> ERROR on species ", species_name)
  }
  
  # Diagnostic plots for each species
  if (!is.null(species_vectors$sen) && length(species_vectors) >= 3) {
    spec_dir <- file.path(species_plots_dir, gsub(" ", "_", sh))
    if (!dir.exists(spec_dir)) dir.create(spec_dir)
    
    plot_colors <- c("Senescence" = "#BEBEBE", "No Peak" = "#56B4E9", "No Early" = "#E69F00", "No Late" = "#009E73")
    
    diag_df <- data.frame(Age = species_vectors$sen$ages, sx_sen = species_vectors$sen$sx, fx_sen = species_vectors$sen$fx,
                          sx_no_peak = species_vectors$no_peak$sx, fx_no_peak = species_vectors$no_peak$fx,
                          sx_early = species_vectors$no_early$sx, fx_early = species_vectors$no_early$fx,
                          sx_late = species_vectors$no_late$sx, fx_late = species_vectors$no_late$fx)
    
    p_sx <- ggplot(diag_df, aes(x = Age)) +
      geom_line(aes(y = sx_sen, color = "Senescence"), size = 1.2) +
      geom_line(aes(y = sx_no_peak, color = "No Peak"),size=1.2) +
      geom_line(aes(y = sx_early, color = "No Early"), size = 1.2) +
      geom_line(aes(y = sx_late, color = "No Late"), size = 1.2) + theme_bw() + labs(title = paste(sh, "sx"))
    
    p_fx <- ggplot(diag_df, aes(x = Age)) +
      geom_line(aes(y = fx_sen, color = "Senescence"), size = 1.2) +
      geom_line(aes(y = fx_no_peak, color = "No Peak"),size=1.2) +
      geom_line(aes(y = fx_early, color = "No Early"), size = 1.2) +
      geom_line(aes(y = fx_late, color = "No Late"), size = 1.2) +  theme_bw() + labs(title = paste(sh, "fx"))
    
    ggsave(file.path(spec_dir, "sx_compare.png"), p_sx, width = 6, height = 4)
    ggsave(file.path(spec_dir, "fx_compare.png"), p_fx, width = 6, height = 4)
  }
} ##<--- end of the FOR loop

## -------------------------------
## 4. Compile Results & Export
## -------------------------------
all_res <- do.call(rbind, comparison_results)
# Calculate percentage change relative to Senescence
df_relative <- all_res %>%
  select(species, model, mean_lifespan, var_lifespan, skew_lifespan, mean_LRO, var_LRO, skew_LRO) %>%
  pivot_wider(names_from = model, values_from = c(mean_lifespan, var_lifespan, skew_lifespan, mean_LRO, var_LRO, skew_LRO)) %>%
  mutate(
    ##------------------
    ##  A. Lifespan 
    ##------------------
    # 1. No-senescence-peak vs Senescence
    pct_diff_NoPeak_Mean_Lifespan = (mean_lifespan_Senescence - `mean_lifespan_No-senescence-peak`) / mean_lifespan_Senescence,
    pct_diff_NoPeak_Var_Lifespan  = (var_lifespan_Senescence - `var_lifespan_No-senescence-peak`)  / var_lifespan_Senescence,
    pct_diff_NoPeak_Skew_Lifespan = (skew_lifespan_Senescence - `skew_lifespan_No-senescence-peak`) / skew_lifespan_Senescence,
    
    # 2. Early onset vs Senescence
    pct_diff_NoEarly_Mean_Lifespan = (mean_lifespan_Senescence - `mean_lifespan_No-senescence-early`) / mean_lifespan_Senescence,
    pct_diff_NoEarly_Var_Lifespan  = (var_lifespan_Senescence - `var_lifespan_No-senescence-early`)  / var_lifespan_Senescence,
    pct_diff_NoEarly_Skew_Lifespan = (skew_lifespan_Senescence - `skew_lifespan_No-senescence-early`) / skew_lifespan_Senescence,
    
    # 3. Late onset vs Senescence
    pct_diff_NoLate_Mean_Lifespan = (mean_lifespan_Senescence - `mean_lifespan_No-senescence-late`) / mean_lifespan_Senescence,
    pct_diff_NoLate_Var_Lifespan  = (var_lifespan_Senescence - `var_lifespan_No-senescence-late`)  / var_lifespan_Senescence,
    pct_diff_NoLate_Skew_Lifespan = (skew_lifespan_Senescence - `skew_lifespan_No-senescence-late`) / skew_lifespan_Senescence,
    
    ##------------------
    ##  B. LRO
    ##------------------
    # 1. No-senescence-peak vs Senescence
    pct_diff_NoPeak_Mean_LRO = (mean_LRO_Senescence - `mean_LRO_No-senescence-peak`) / mean_LRO_Senescence,
    pct_diff_NoPeak_Var_LRO  = (var_LRO_Senescence - `var_LRO_No-senescence-peak`)  / var_LRO_Senescence,
    pct_diff_NoPeak_Skew_LRO = (skew_LRO_Senescence - `skew_LRO_No-senescence-peak`) / skew_LRO_Senescence,
    
    # 2. Early onset vs Senescence
    pct_diff_NoEarly_Mean_LRO = (mean_LRO_Senescence - `mean_LRO_No-senescence-early`) / mean_LRO_Senescence,
    pct_diff_NoEarly_Var_LRO  = (var_LRO_Senescence - `var_LRO_No-senescence-early`)  / var_LRO_Senescence,
    pct_diff_NoEarly_Skew_LRO = (skew_LRO_Senescence - `skew_LRO_No-senescence-early`) / skew_LRO_Senescence,
    
    # 3. Late onset vs Senescence
    pct_diff_NoLate_Mean_LRO = (mean_LRO_Senescence - `mean_LRO_No-senescence-late`) / mean_LRO_Senescence,
    pct_diff_NoLate_Var_LRO  = (var_LRO_Senescence - `var_LRO_No-senescence-late`)  / var_LRO_Senescence,
    pct_diff_NoLate_Skew_LRO = (skew_LRO_Senescence - `skew_LRO_No-senescence-late`) / skew_LRO_Senescence
  ) %>%
  # Extract the calculated difference columns for reshaping
  select(species, starts_with("pct_diff")) %>%
  pivot_longer(cols = -species, names_to = "comparison", values_to = "pct_change") %>%
  # Split column names exactly into 5 parts: pct, diff, Model_Code, Metric_Code, Trait_Code
  separate(comparison, into = c("dummy", "dummy2", "Model_Code", "Metric_Code", "Trait_Code"), sep = "_") %>%
  select(-dummy, -dummy2) %>%
  mutate(
    # Convert to percentage values
    pct_change = pct_change * 100,
    
    # Map back to full model names
    Model = case_when(
      Model_Code == "NoPeak"  ~ "No-senescence-peak",
      Model_Code == "NoEarly" ~ "No-senescence-early",
      Model_Code == "NoLate"  ~ "No-senescence-late"
    ),
    
    # Reconstruct full metric names (e.g., combine "Mean" and "Lifespan")
    Metric = case_when(
      Metric_Code == "Mean" ~ paste("Mean", Trait_Code),
      Metric_Code == "Var"  ~ paste("Variance", Trait_Code),
      Metric_Code == "Skew" ~ paste("Skewness", Trait_Code)
    ),
    
    # Set factor levels for correct plotting order
    Model = factor(Model, levels = c("No-senescence-peak", "No-senescence-early", "No-senescence-late")),
    Metric = factor(Metric, levels = c(
      "Mean Lifespan", "Variance Lifespan", "Skewness Lifespan", 
      "Mean LRO", "Variance LRO", "Skewness LRO"
    ))
  )

write.csv(df_relative, file.path(output_dir, "methodology_full_results.csv"), row.names = FALSE)

## -------------------------------
## 5. Summary Plots with Improved Visuals
## -------------------------------
df_plot <- df_relative %>%
  select(species, Model, Metric, Trait_Code, pct_change)

create_comparison_plot <- function(data_subset, title_suffix) {
  
  # x-axis is now mapped to 'Model' and y-axis to 'pct_change'
  ggplot(data_subset, aes(x = Model, y = pct_change)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    
    # Violin and boxplots
    geom_violin(aes(fill = Model), alpha = 0.2, color = NA) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
    
    # Facet by specific metrics (Mean, Variance, Skewness)
    facet_wrap(~Metric, scales = "free_y", ncol = 3) +
    
    # The data is already in percentage format (multiplied by 100), 
    # so no need to multiply by 100 here. Just round it and add '%'. 
    # Sigma is set to 10 for better visualization on the percentage scale.
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 10),
                       labels = function(x) paste0(round(x), "%")) +
    
    # Updated to match the current three model names
    scale_fill_manual(values = c("No-senescence-peak"  = "#56B4E9", 
                                 "No-senescence-early" = "#E69F00",
                                 "No-senescence-late"  = "#009E73")) +
    theme_bw(base_size = 12) +
    
    coord_cartesian(clip = "off") +
    
    labs(title = paste("Methodological Sensitivity:", title_suffix),
         subtitle = "Distributions of relative differences",
         y = "Relative Difference (%)", 
         x = "", 
         fill = "Model") +
    
    theme(legend.position = "bottom", 
          plot.title = element_text(face="bold"),
          # Left margin reduced to 10 since there are no left-aligned labels anymore
          plot.margin = margin(10, 10, 10, 10), 
          # Slightly tilt x-axis text to prevent overlapping
          axis.text.x = element_text(angle = 15, hjust = 1)) 
}

# Use Trait_Code instead of Metric_Type for filtering
p_life <- create_comparison_plot(df_plot %>% filter(Trait_Code == "Lifespan"), "Lifespan")
p_lro  <- create_comparison_plot(df_plot %>% filter(Trait_Code == "LRO"), "LRO")

ggsave(file.path("Results/summary figures", "Summary_Lifespan_onset_of_senescence.png"), p_life, width = 14, height = 7)
ggsave(file.path("Results/summary figures","Summary_LRO_Final_onset_of_senescence.png"), p_lro, width = 14, height = 7)

message("Done! Files saved to: ", "Results/summary figures")