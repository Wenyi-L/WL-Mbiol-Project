## =========================================================
## plot_summary_grouped_comparison.R (Updated 02/03/2026)
## =========================================================

rm(list = ls())

## -------------------------------
## 0. Configuration
## -------------------------------
excel_file   <- "Jones2014.xls"
summary_file <- "output_batch/all_species_summary_stats.csv"
fig_dir      <- "summary_figures_comparison_v8"

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

## -------------------------------
## 1. Load packages
## -------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(readxl)

## -------------------------------
## 2. Metadata Extraction
## -------------------------------
if (!file.exists(excel_file)) {
  stop(paste("Error: File not found:", excel_file))
}

message("Fetching metadata from ", excel_file, "...")
sheets <- excel_sheets(excel_file)
metadata_list <- list()

for (sh in sheets) {
  class_val <- suppressMessages(
    as.character(read_excel(excel_file, sheet = sh, range = "D1:D1", col_names = FALSE)[1,1])
  )
  metadata_list[[length(metadata_list) + 1]] <- data.frame(species = sh, Class = class_val)
}
df_metadata <- do.call(rbind, metadata_list)

## -------------------------------
## 3. Data Processing & Grouping
## -------------------------------
if (!file.exists(summary_file)) {
  stop(paste("Error: File not found:", summary_file))
}

df <- read.csv(summary_file, stringsAsFactors = FALSE)

# Calculate Percentage Difference #No longer calculate absolute percentage difference
df_diff <- df_merged %>%
  filter(model %in% c("Senescence", "No-senescence")) %>%
  select(species, Group, model, mean_lifespan, var_lifespan, skew_lifespan, mean_LRO, var_LRO, skew_LRO) %>%
  pivot_longer(cols = -c(species, Group, model), names_to = "metric", values_to = "value") %>%
  pivot_wider(names_from = model, values_from = value) %>%
  mutate(
    perc_diff = ifelse(
      is.na(Senescence) | is.na(`No-senescence`) | Senescence == 0,
      NA_real_,
      100 * (Senescence - `No-senescence`) / Senescence)
  )


## -------------------------------
## 4. Plotting Function
## -------------------------------
create_comparison_plot <- function(data, title, filename) {
  
  dodge_w <- 0.8
  
  p <- ggplot(data, aes(x = metric_label, y = perc_diff, fill = Group)) +
    
    # --- Layer 1: Violin (Background Color) ---
    geom_violin(
      position = position_dodge(width = dodge_w),
      scale = "width",
      trim = FALSE,
      alpha = 0.5,       # Slightly more opaque to show the yellow/green better
      color = NA,
      aes(fill = model)
    ) +
    
    # --- Layer 2: Boxplot (Transparent Fill) ---
    geom_boxplot(
      aes(group = interaction(metric_label, Group)), 
      width = 0.2,
      position = position_dodge(width = dodge_w),
      outlier.shape = NA,
      fill = "transparent",
      color = "grey30",
      size = 0.4
    ) +
    
    # --- Layer 3: Jitter Points (Darker Colors) ---
    geom_point(
      aes(color = Group),
      position = position_jitterdodge(dodge.width = dodge_w, jitter.width = 0.2),
      size = 1.2, 
      alpha = 0.8
    ) +
    
    # Scales
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma = 1),
      breaks = c(0, 1, 10, 100, 1000, 10000),
      labels = scales::comma_format()
    ) +
    
    
    # Theme
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 14),
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
    ) +
    labs(
      title = title,
      x     = NULL, 
      y     = " % Change (Log Scale)"
    )
  
  ggsave(file.path(fig_dir, filename), p, width = 11, height = 7, dpi = 300)
}

## -------------------------------
## 6. Generate Plots
## -------------------------------

message("Generating plots with Green(Plant) and Yellow(Animal)...")

# A. Lifespan Comparison
df_LS <- df_diff %>%
  filter(metric %in% c("mean_lifespan", "var_lifespan", "skew_lifespan")) %>%
  mutate(metric_label = factor(metric, 
                               levels = c("mean_lifespan", "var_lifespan", "skew_lifespan"),
                               labels = c("Mean", "Variance", "Skewness"))) %>%
  filter(!is.na(perc_diff))

create_comparison_plot(df_LS, "Lifespan Sensitivity: Senescence VS No-senescence", "compare_lifespan_sen_nosen.png")

# B. LRO Comparison
df_LRO <- df_diff %>%
  filter(metric %in% c("mean_LRO", "var_LRO", "skew_LRO")) %>%
  mutate(metric_label = factor(metric, 
                               levels = c("mean_LRO", "var_LRO", "skew_LRO"),
                               labels = c("Mean", "Variance", "Skewness"))) %>%
  filter(!is.na(perc_diff))

create_comparison_plot(df_LRO, "LRO Sensitivity: Senescence VS No-senescence", "compare_LRO_sen_nosen.png")

message("Done! Plots saved to: ", fig_dir)