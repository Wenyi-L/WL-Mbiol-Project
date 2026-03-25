## =========================================================
## run_different_senescence_onset.R
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

# CH addressed: We now process data ONCE in Loop 1, save the objects, 
# and use the pre-calculated data for plotting in Loop 2!

## -------------------------------
## 2. Main Data Processing (Loop 1)
## -------------------------------
sheets <- excel_sheets(excel_file)

# Lists to store pre-calculated data
all_onset_results <- list()
all_demog_tables <- list()

message("Starting analysis loop (Calculation only)...")

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
    ## 1. Prepare demographic inputs (Includes s0 adjustment & Age 0 removal)
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
    
    ## 3. Summary statistics
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
    
    # Bundle everything to save for plotting later
    list(
      summary = summary_df,
      demog = demog,
      mpm_sen = mpm_sen,
      mpm_no_peak = mpm_no_peak,
      mpm_no_early = mpm_no_early,
      mpm_no_late = mpm_no_late
    )
  }) #<- end of try()
  
  ## 4. Save section
  if (!inherits(res, "try-error")) {
    message("  -> OK: added species ", species_name)
    
    # Save into the global list
    all_onset_results[[species_name]] <- res
    
    # Save sx/fx to master table
    sp_demog <- data.frame(
      Species = species_name,
      Age = res$demog$ages,
      sx = res$demog$sx,
      fx = res$demog$fx
    )
    all_demog_tables[[species_name]] <- sp_demog
    
    # Save individual CSV for this species
    write.csv(
      res$summary,
      file = file.path(output_dir, paste0(species_name, "_summary_stats.csv")),
      row.names = FALSE
    )
  } else {
    message("  -> ERROR on species ", species_name)
  }
} ##<--- end of LOOP 1

## -------------------------------
## 3. Export Master Demography CSV
## -------------------------------
if(length(all_demog_tables) > 0) {
  master_demog <- do.call(rbind, all_demog_tables)
  write.csv(master_demog, file.path(output_dir, "All_Species_Onset_Demography.csv"), row.names = FALSE)
  message("\nMaster Demography CSV exported successfully!")
}

## -------------------------------
## 4. Plotting Diagnostic Plots (Loop 2)
## -------------------------------
message("\nStarting diagnostic plotting loop...")

for (sh in names(all_onset_results)) {
  
  res <- all_onset_results[[sh]]
  demog <- res$demog
  mpm_sen <- res$mpm_sen
  mpm_no_peak <- res$mpm_no_peak
  mpm_no_early <- res$mpm_no_early
  mpm_no_late <- res$mpm_no_late
  
  spec_dir <- file.path(species_plots_dir, gsub(" ", "_", sh))
  if (!dir.exists(spec_dir)) dir.create(spec_dir)
  
  plot_colors <- c("Senescence" = "#BEBEBE", "No Peak" = "#56B4E9", "No Early" = "#E69F00", "No Late" = "#009E73")
  
  diag_df <- data.frame(
    Age = demog$ages, 
    sx_sen = mpm_sen$sx, fx_sen = mpm_sen$fx,
    sx_no_peak = mpm_no_peak$sx, fx_no_peak = mpm_no_peak$fx,
    sx_early = mpm_no_early$sx, fx_early = mpm_no_early$fx,
    sx_late = mpm_no_late$sx, fx_late = mpm_no_late$fx
  )
  
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
} ##<--- end of LOOP 2

## -------------------------------
## 5. Compile Results & Export
## -------------------------------
# Extract just the summary_df from our pre-calculated list
summary_list <- lapply(all_onset_results, function(x) x$summary)
all_res <- do.call(rbind, summary_list)

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
## 6. Summary Plots with Improved Visuals
## -------------------------------
df_plot <- df_relative %>%
  select(species, Model, Metric, Trait_Code, pct_change)

create_comparison_plot <- function(data_subset, title_suffix) {
  
  ggplot(data_subset, aes(x = Model, y = pct_change)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    
    # Violin and boxplots
    geom_violin(aes(fill = Model), alpha = 0.2, color = NA) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
    
    # Facet by specific metrics (Mean, Variance, Skewness)
    facet_wrap(~Metric, scales = "free_y", ncol = 3) +
    
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
          plot.margin = margin(10, 10, 10, 10), 
          axis.text.x = element_text(angle = 15, hjust = 1)) 
}

# Use Trait_Code instead of Metric_Type for filtering
p_life <- create_comparison_plot(df_plot %>% filter(Trait_Code == "Lifespan"), "Lifespan")
p_lro  <- create_comparison_plot(df_plot %>% filter(Trait_Code == "LRO"), "LRO")

if(!dir.exists("Results/summary figures")) dir.create("Results/summary figures", recursive = TRUE)

ggsave(file.path("Results/summary figures", "Summary_Lifespan_onset_of_senescence.png"), p_life, width = 14, height = 7)
ggsave(file.path("Results/summary figures","Summary_LRO_Final_onset_of_senescence.png"), p_lro, width = 14, height = 7)

message("Done! Files saved to: Results/summary figures")