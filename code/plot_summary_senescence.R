## =========================================================
## plot_summary_grouped_comparison.R (Updated 02/03/2026)
## =========================================================

rm(list = ls())

## -------------------------------
## 0. Configuration
## -------------------------------
summary_file <- "Results/senescence analysis/all_species_summary_stats.csv"
fig_dir      <- "Results/summary figures"

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
## 2. Data Processing & Grouping
## -------------------------------
if (!file.exists(summary_file)) {
  stop(paste("Error: File not found:", summary_file))
}

df <- read.csv(summary_file, stringsAsFactors = FALSE)

# Calculate Percentage Difference #No longer calculate absolute percentage difference
df_diff <- df %>%
  filter(model %in% c("Senescence", "No-senescence")) %>%
  select(species,model, mean_lifespan, var_lifespan, skew_lifespan, mean_LRO, var_LRO, skew_LRO) %>%
  pivot_longer(cols = -c(species, model), names_to = "metric", values_to = "value") %>%
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

my_greens <- c(
  "#B4EEB4", #DarkSeaGreen2
  "#9ACD32",  # YellowGreen 
  "#3CB371" # MediumSeaGreen 
)
create_comparison_plot <- function(data, title, filename) {
  
  dodge_w <- 0.8
  
  p <- ggplot(data, aes(x = metric, y = perc_diff)) +
    
    # --- Layer 1: Boxplot (Green Fill) ---
    geom_boxplot(
      aes(fill = metric),
      width = 0.5,
      position = position_dodge(width = dodge_w),
      outlier.shape = NA,
      color = "grey30",
      alpha = 0.5,
      size = 0.4
    ) +
    
    # --- Layer 2: Jitter Points (Darker Colors) ---
    geom_point(
      aes(color = metric),
      position = position_jitterdodge(dodge.width = dodge_w, jitter.width = 0.2),
      size = 1.2, 
      alpha = 0.8
    ) +
    
    #Colour
    scale_fill_manual(values = my_greens) +
    scale_color_manual(values = my_greens) +
    
    # Scales
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma = 1),
      breaks = c(-10000, -1000, -100, -10, -1, 0, 1, 10, 100, 1000, 10000),
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

# A. Lifespan Comparison
df_LS <- df_diff %>%
  filter(metric %in% c("mean_lifespan", "var_lifespan", "skew_lifespan")) %>%
  mutate(metric = factor(metric, 
                         levels = c("mean_lifespan", "var_lifespan", "skew_lifespan"),
                         labels = c("Mean", "Variance", "Skewness"))) %>%
  filter(!is.na(perc_diff))

create_comparison_plot(df_LS, "Lifespan Sensitivity: Senescence VS No-senescence", "compare_lifespan_sen_nosen.png")

# B. LRO Comparison
df_LRO <- df_diff %>%
  filter(metric %in% c("mean_LRO", "var_LRO", "skew_LRO")) %>%
  mutate(metric = factor(metric, 
                         levels = c("mean_LRO", "var_LRO", "skew_LRO"),
                         labels = c("Mean", "Variance", "Skewness"))) %>%
  filter(!is.na(perc_diff))

create_comparison_plot(df_LRO, "LRO Sensitivity: Senescence VS No-senescence", "compare_LRO_sen_nosen.png")

message("Done! Plots saved to: ", fig_dir)