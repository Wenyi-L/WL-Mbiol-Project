## =========================================================
## run_peak_vs_maturity_comparison.R
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
excel_file <- "Jones2014.xls" 
output_dir <- "output_peak_weighted" # Reverted to your preferred directory
species_plots_dir <- file.path(output_dir, "species_diagnostic_plots")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(species_plots_dir)) dir.create(species_plots_dir, recursive = TRUE)

# Ensure helper scripts are present
source("LuckFunctions.R")
source("senescence_functions.R") 

## -------------------------------
## 2. Helper: No-Senescence Weighted Mean
## -------------------------------
create_nosen_vectors <- function(ages, sx, fx, maturity_age) {
  k <- length(ages)
  lx <- numeric(k)
  lx[1] <- 1
  if (k > 1) {
    for (i in 1:(k-1)) lx[i+1] <- lx[i] * sx[i]
  }
  
  idx_adult <- which(ages >= maturity_age)
  if (length(idx_adult) == 0) idx_adult <- k
  
  w <- lx[idx_adult]
  if (sum(w) == 0) w <- rep(1, length(w))
  
  sx_adult_vals <- sx[idx_adult]
  fx_adult_vals <- fx[idx_adult]
  
  sx_mean <- weighted.mean(sx_adult_vals, w, na.rm = TRUE)
  fx_mean <- weighted.mean(fx_adult_vals, w, na.rm = TRUE)
  
  sx_nosen <- sx
  fx_nosen <- fx
  sx_nosen[idx_adult] <- sx_mean
  fx_nosen[idx_adult] <- fx_mean
  
  sx_nosen[!is.finite(sx_nosen)] <- 0
  fx_nosen[!is.finite(fx_nosen)] <- 0
  
  return(list(sx = sx_nosen, fx = fx_nosen))
}

## -------------------------------
## 3. Main Data Processing
## -------------------------------
sheets <- excel_sheets(excel_file)
comparison_results <- list()

message("Starting analysis loop...")

for (sh in sheets) {
  # Extract Taxon from D1
  class_val <- as.character(read_excel(excel_file, sheet = sh, range = "D1:D1", col_names = FALSE)[1,1])
  
  df_raw <- try(read_excel(excel_file, sheet = sh, skip = 1), silent = TRUE)
  if (inherits(df_raw, "try-error") || !("x" %in% names(df_raw))) next
  df <- as.data.frame(df_raw)
  
  species_vectors <- list()
  methods_list <- c("Logistic_50", "Peak_Fertility")
  
  for (method in methods_list) {
    tryCatch({
      if (method == "Logistic_50") {
        demog <- prepare_demography_data_from_df(df, maturity_method = "logistic", maturity_prob = 0.50)
      } else {
        temp <- prepare_demography_data_from_df(df, maturity_age = 1)
        peak_idx <- which.max(temp$fx)
        demog <- temp
        demog$maturity_age <- temp$ages[peak_idx]
      }
      
      mpm_sen <- build_MPM_senescence(demog$ages, demog$sx, demog$fx)
      nosen_vecs <- create_nosen_vectors(demog$ages, demog$sx, demog$fx, demog$maturity_age)
      mpm_no <- build_MPM_senescence(demog$ages, nosen_vecs$sx, nosen_vecs$fx)
      
      species_vectors[[method]] <- list(sx = nosen_vecs$sx, fx = nosen_vecs$fx)
      if(method == "Logistic_50") species_vectors[["Sen"]] <- list(ages = demog$ages, sx = demog$sx, fx = demog$fx)
      
      mix_sen <- mixing_distro(mpm_sen$A, mpm_sen$F)
      mix_no  <- mixing_distro(mpm_no$A, mpm_no$F)
      summ <- compute_summary_table(mpm_sen$U, mpm_no$U, mpm_sen$U, mpm_sen$U, 
                                    mpm_sen$F, mpm_no$F, mpm_sen$F, mpm_sen$F,
                                    mix_sen, mix_no, mix_sen, mix_sen, repro_var = "Poisson")
      
      summ$species <- sh
      summ$Class   <- class_val
      summ$Method  <- method
      comparison_results[[length(comparison_results) + 1]] <- summ %>% filter(model %in% c("Senescence", "No-senescence"))
    }, error = function(e) { })
  }
  
  # Diagnostic plots for each species
  if (!is.null(species_vectors$Sen) && length(species_vectors) >= 3) {
    spec_dir <- file.path(species_plots_dir, gsub(" ", "_", sh))
    if (!dir.exists(spec_dir)) dir.create(spec_dir)
    diag_df <- data.frame(Age = species_vectors$Sen$ages, sx_Raw = species_vectors$Sen$sx, fx_Raw = species_vectors$Sen$fx,
                          sx_Logi = species_vectors$Logistic_50$sx, fx_Logi = species_vectors$Logistic_50$fx,
                          sx_Peak = species_vectors$Peak_Fertility$sx, fx_Peak = species_vectors$Peak_Fertility$fx)
    
    p_sx <- ggplot(diag_df, aes(x = Age)) +
      geom_line(aes(y = sx_Raw, color = "Raw"), size = 1.2) +
      geom_line(aes(y = sx_Logi, color = "Logi_50"), linetype = "dashed") +
      geom_line(aes(y = sx_Peak, color = "Peak_Fert"), linetype = "dotted", size = 1.1) +
      scale_color_manual(values = c("black", "#56B4E9", "#E69F00")) + theme_bw() + labs(title = paste(sh, "sx"))
    
    p_fx <- ggplot(diag_df, aes(x = Age)) +
      geom_line(aes(y = fx_Raw, color = "Raw"), size = 1.2) +
      geom_line(aes(y = fx_Logi, color = "Logi_50"), linetype = "dashed") +
      geom_line(aes(y = fx_Peak, color = "Peak_Fert"), linetype = "dotted", size = 1.1) +
      scale_color_manual(values = c("black", "#56B4E9", "#E69F00")) + theme_bw() + labs(title = paste(sh, "fx"))
    
    ggsave(file.path(spec_dir, "sx_compare.png"), p_sx, width = 6, height = 4)
    ggsave(file.path(spec_dir, "fx_compare.png"), p_fx, width = 6, height = 4)
  }
}

## -------------------------------
## 4. Compile Results & Export
## -------------------------------
all_res <- do.call(rbind, comparison_results)
df_diff <- all_res %>%
  select(species, Class, Method, model, mean_lifespan, var_lifespan, skew_lifespan, mean_LRO, var_LRO, skew_LRO) %>%
  pivot_wider(names_from = model, values_from = c(mean_lifespan, var_lifespan, skew_lifespan, mean_LRO, var_LRO, skew_LRO)) %>%
  mutate(Diff_Mean_Lifespan = (mean_lifespan_Senescence - `mean_lifespan_No-senescence`) / mean_lifespan_Senescence,
         Diff_Var_Lifespan  = (var_lifespan_Senescence  - `var_lifespan_No-senescence`)  / var_lifespan_Senescence,
         Diff_Skew_Lifespan = (skew_lifespan_Senescence - `skew_lifespan_No-senescence`) / skew_lifespan_Senescence,
         Diff_Mean_LRO      = (mean_LRO_Senescence      - `mean_LRO_No-senescence`)      / mean_LRO_Senescence,
         Diff_Var_LRO       = (var_LRO_Senescence       - `var_LRO_No-senescence`)       / var_LRO_Senescence,
         Diff_Skew_LRO      = (skew_LRO_Senescence      - `skew_LRO_No-senescence`)      / skew_LRO_Senescence)

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