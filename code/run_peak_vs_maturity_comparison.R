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
output_dir <- "Results/summary figures"
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
  if (!is.null(species_vectors$Sen) && length(species_vectors) >= 3) {
    spec_dir <- file.path(species_plots_dir, gsub(" ", "_", sh))
    if (!dir.exists(spec_dir)) dir.create(spec_dir)
    
    plot_colors <- c("Senescence" = "#BEBEBE", "No Peak" = "#56B4E9", "No Early" = "#E69F00", "No Late" = "#009E73")
    
    diag_df <- data.frame(Age = species_vectors$sen$ages, sx_sen = species_vectors$sen$sx, fx_sen = species_vectors$sen$fx,
                          sx_no_peak = species_vectors$no_peak$sx, fx_no_peak = species_vectors$no_peak$fx,
                          sx_early = species_vectors$no_early$sx, fx_early = species_vectors$no_early$fx,
                          sx_late = species_vectors$no_late$sx, fx_no_peak = species_vectors$no_late$fx,)
    
    p_sx <- ggplot(diag_df, aes(x = Age)) +
      geom_line(aes(y = sx_sen, color = "Senescence"), size = 1.2) +
      geom_line(aes(y = sx_no_peak, color = "No Peak"),size=1.2) +
      geom_line(aes(y = sx_no_early, color = "No Early"), size = 1.2) +
      geom_line(aes(y = sx_no_late, color = "No Late"), size = 1.2) + theme_bw() + labs(title = paste(sh, "sx"))
    
    p_fx <- ggplot(diag_df, aes(x = Age)) +
      geom_line(aes(y = fx_sen, color = "Senescence"), size = 1.2) +
      geom_line(aes(y = fx_no_peak, color = "No Peak"),size=1.2) +
      geom_line(aes(y = fx_no_early, color = "No Early"), size = 1.2) +
      geom_line(aes(y = fx_no_late, color = "No Late"), size = 1.2) +  theme_bw() + labs(title = paste(sh, "fx"))
    
    ggsave(file.path(spec_dir, "sx_compare.png"), p_sx, width = 6, height = 4)
    ggsave(file.path(spec_dir, "fx_compare.png"), p_fx, width = 6, height = 4)
  }
}

## -------------------------------
## 4. Compile Results & Export
## -------------------------------
all_res <- do.call(rbind, comparison_results)
# Calculate percentage change relative to Senescence
df_relative <- all_res %>%
  select(species, model, mean_lifespan, var_lifespan,skew_lifespan,mean_LRO, var_LRO, skew_LRO) %>%
  pivot_wider(names_from = model, values_from = c(mean_lifespan, var_lifespan,skew_lifespan,mean_LRO, var_LRO, skew_LRO)) %>%
  mutate(
    # NOTE: Using backticks (`) to handle hyphens and slashes in column names
    
    # 1. No-senescence vs Senescence
    pct_diff_NoSen_Mean = (`mean_LRO_No-senescence` - mean_LRO_Senescence) / mean_LRO_Senescence,
    pct_diff_NoSen_Var  = (`var_LRO_No-senescence` - var_LRO_Senescence)  / var_LRO_Senescence,
    pct_diff_NoSen_Skew = (`skew_LRO_No-senescence` - skew_LRO_Senescence) / skew_LRO_Senescence,
    
    # 2. No-Actuarial vs Senescence
    pct_diff_NoAct_Mean = (`mean_LRO_No-actuarial/Yes-reproductive` - mean_LRO_Senescence) / mean_LRO_Senescence,
    pct_diff_NoAct_Var  = (`var_LRO_No-actuarial/Yes-reproductive`  - var_LRO_Senescence)  / var_LRO_Senescence,
    pct_diff_NoAct_Skew = (`skew_LRO_No-actuarial/Yes-reproductive` - skew_LRO_Senescence) / skew_LRO_Senescence,
    
    # 3. Yes-Actuarial (No-Repro) vs Senescence
    pct_diff_NoRep_Mean = (`mean_LRO_Yes-actuarial/No-reproductive` - mean_LRO_Senescence) / mean_LRO_Senescence,
    pct_diff_NoRep_Var  = (`var_LRO_Yes-actuarial/No-reproductive`  - var_LRO_Senescence)  / var_LRO_Senescence,
    pct_diff_NoRep_Skew = (`skew_LRO_Yes-actuarial/No-reproductive` - skew_LRO_Senescence) / skew_LRO_Senescence
  ) %>%
  select(species, starts_with("pct_diff")) %>%
  pivot_longer(cols = -species, names_to = "comparison", values_to = "pct_change") %>%
  separate(comparison, into = c("dummy", "dummy2", "Model_Code", "Metric_Code"), sep = "_") %>%
  select(-dummy, -dummy2) %>%
  mutate(
    # Convert to percentage values
    pct_change = pct_change * 100,
    # Map back to full model names
    Model = case_when(
      Model_Code == "NoSen" ~ "No-senescence",
      Model_Code == "NoAct" ~ "No-actuarial/Yes-reproductive",
      Model_Code == "NoRep" ~ "Yes-actuarial/No-reproductive"
    ),
    Metric = case_when(
      Metric_Code == "Mean" ~ "Mean LRO",
      Metric_Code == "Var"  ~ "Variance LRO",
      Metric_Code == "Skew" ~ "Skewness LRO"
    ),
    # Set factor levels for plotting order
    Model = factor(Model, levels = c("No-senescence", "No-actuarial/Yes-reproductive", "Yes-actuarial/No-reproductive")),
    Metric = factor(Metric, levels = c("Mean LRO", "Variance LRO", "Skewness LRO"))
  )

write.csv(df_diff, file.path(output_dir, "methodology_full_results.csv"), row.names = FALSE)

## -------------------------------
## 5. Summary Plots with Improved Visuals
## -------------------------------
df_plot <- df_diff %>%
  select(species, Class, Method, starts_with("Diff_")) %>%
  pivot_longer(cols = starts_with("Diff_"), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric_Type = ifelse(grepl("Lifespan", Metric), "Lifespan", "LRO"),
         Metric_Label = case_when(Metric == "Diff_Mean_Lifespan" ~ "Mean Lifespan",
                                  Metric == "Diff_Var_Lifespan"  ~ "Variance Lifespan",
                                  Metric == "Diff_Skew_Lifespan" ~ "Skewness Lifespan",
                                  Metric == "Diff_Mean_LRO"      ~ "Mean LRO",
                                  Metric == "Diff_Var_LRO"       ~ "Variance LRO",
                                  Metric == "Diff_Skew_LRO"      ~ "Skewness LRO"))

create_comparison_plot <- function(data_subset, title_suffix) {
  # Select sensitive species for labels
  df_sensitive <- data_subset %>% filter(abs(Value) > 0.5)
  df_labels_left <- df_sensitive %>% filter(Method == "Logistic_50")
  
  ggplot(data_subset, aes(x = Method, y = Value)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_violin(aes(fill = Method), alpha = 0.2, color = NA) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
    
    # Lines for Taxonomic Class (This is the source for the color legend)
    geom_line(aes(group = species, color = Class), alpha = 0.4) + 
    
    # Points - show.legend = FALSE to prevent shapes from appearing in the Class legend
    geom_point(alpha = 0.2, shape = 16, size = 1, show.legend = FALSE) +
    
    # Sensitive Points as Triangles - show.legend = FALSE to keep legend as lines
    geom_point(data = df_sensitive, aes(color = Class), shape = 17, size = 3, show.legend = FALSE) +
    
    # Species labels (Left side only)
    geom_text(data = df_labels_left, aes(label = species), 
              size = 2.5, hjust = 1.1, nudge_x = -0.05) +
    
    facet_wrap(~Metric_Label, scales = "free_y", ncol = 3) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 0.1),
                       labels = function(x) paste0(round(x * 100), "%")) +
    
    # scale_x_discrete expand: adds space on the left side for labels
    scale_x_discrete(expand = expansion(mult = c(0.4, 0.1))) +
    
    scale_fill_manual(values = c("Logistic_50" = "#56B4E9", "Peak_Fertility" = "#E69F00")) +
    theme_bw(base_size = 12) +
    
    # clip = "off": ensures text outside the panel remains visible
    coord_cartesian(clip = "off") +
    
    labs(title = paste("Methodological Sensitivity:", title_suffix),
         subtitle = "Triangles: >50% Diff | Labels: Left-side only | Legend: Lines for Taxon",
         y = "Relative Difference (Sen - NoSen)/Sen", x = "", color = "Taxonomic Class") +
    theme(legend.position = "bottom", 
          plot.title = element_text(face="bold"),
          plot.margin = margin(10, 10, 10, 60)) # Extra left margin for long names
}

p_life <- create_comparison_plot(df_plot %>% filter(Metric_Type == "Lifespan"), "Lifespan")
p_lro  <- create_comparison_plot(df_plot %>% filter(Metric_Type == "LRO"), "LRO")

ggsave(file.path(output_dir, "Summary_Lifespan_Final_Fixed.png"), p_life, width = 14, height = 7)
ggsave(file.path(output_dir, "Summary_LRO_Final_Fixed.png"), p_lro, width = 14, height = 7)

message("Done! Files saved to: ", output_dir)