## =========================================================
## plot_methodology_comparison_v5_repel.R
##
## Purpose: 
## 1. Load CSV from D:/Nian/Mbiol/...
## 2. Fetch Taxonomic Class from Jones2014.xls (Cell D1).
## 3. Visualize with ggrepel to prevent label overlap.
## 4. Paths are fixed as requested.
## =========================================================

# Ensure ggrepel is installed
if (!require("ggrepel")) install.packages("ggrepel")

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggrepel) # Added for smart labeling

# 1. Paths (Fixed as requested)
excel_file <- "D:/Nian/Mbiol/Jones2014.xls"
input_csv  <- "D:/Nian/Mbiol/output_peak_weighted/peak_vs_50_diff.csv"
output_dir <- "D:/Nian/Mbiol/output_peak_weighted" 

if(!dir.exists(output_dir)) dir.create(output_dir)

# 2. Fetch Metadata (Class) from Excel
message("Fetching taxonomic classes from Jones2014.xls...")
sheets <- excel_sheets(excel_file)
metadata_list <- list()

for (sh in sheets) {
  # Read only the single cell D1 for Taxonomic Class
  class_val <- as.character(read_excel(excel_file, sheet = sh, range = "D1:D1", col_names = FALSE)[1,1])
  metadata_list[[length(metadata_list) + 1]] <- data.frame(species = sh, Class = class_val)
}
df_metadata <- do.call(rbind, metadata_list)

# 3. Load CSV and Merge with Metadata
df_results <- read.csv(input_csv)

# Merge: remove existing Class column in CSV (if any) to avoid duplicates like Class.x/Class.y
df_final <- df_results %>%
  select(-any_of("Class")) %>% 
  left_join(df_metadata, by = "species") %>%
  drop_na(starts_with("Diff_")) # Clean missing values

# 4. Reshape data for Faceting
df_plot <- df_final %>%
  select(species, Class, Method, starts_with("Diff_")) %>%
  pivot_longer(cols = starts_with("Diff_"), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric_Type = ifelse(grepl("Lifespan", Metric), "Lifespan", "LRO"),
         Metric_Label = case_when(Metric == "Diff_Mean_Lifespan" ~ "Mean Lifespan",
                                  Metric == "Diff_Var_Lifespan"  ~ "Variance Lifespan",
                                  Metric == "Diff_Skew_Lifespan" ~ "Skewness Lifespan",
                                  Metric == "Diff_Mean_LRO"      ~ "Mean LRO",
                                  Metric == "Diff_Var_LRO"       ~ "Variance LRO",
                                  Metric == "Diff_Skew_LRO"      ~ "Skewness LRO")) %>%
  # Set order for metrics
  mutate(Metric_Label = factor(Metric_Label, 
                               levels = c("Mean Lifespan", "Skewness Lifespan", "Variance Lifespan",
                                          "Mean LRO", "Skewness LRO", "Variance LRO")))

# 5. Professional Plotting Function with ggrepel
create_final_plot <- function(data_subset, title_suffix) {
  
  # Filter for species with > 50% relative difference
  df_sensitive <- data_subset %>% filter(abs(Value) > 0.5)
  # Only show labels for the Left side (Logistic_50)
  df_labels_left <- df_sensitive %>% filter(Method == "Logistic_50")
  
  ggplot(data_subset, aes(x = Method, y = Value)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_violin(aes(fill = Method), alpha = 0.2, color = NA) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
    
    # Class lines (Legend source)
    geom_line(aes(group = species, color = Class), alpha = 0.4) + 
    
    # Points (Hidden from legend)
    geom_point(alpha = 0.2, shape = 16, size = 1, show.legend = FALSE) +
    
    # Highlighted points (Triangles, hidden from legend)
    geom_point(data = df_sensitive, aes(color = Class), 
               shape = 17, size = 3, show.legend = FALSE) +
    
    # SMART LABELS (Updated to use ggrepel)
    geom_text_repel(data = df_labels_left, aes(label = species), 
                    size = 2.8, 
                    direction = "y",      # Stack vertically
                    nudge_x = -0.6,       # Push to the left significantly
                    hjust = 1,            # Right align text
                    segment.size = 0.2,   # Thin connecting lines
                    segment.color = "grey50",
                    min.segment.length = 0, # Always draw lines
                    max.overlaps = Inf,   # FORCE all labels to show
                    force = 2) +          # Strength of separation
    
    facet_wrap(~Metric_Label, scales = "free_y", ncol = 3) +
    
    # Formatting
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 0.1),
                       labels = function(x) paste0(round(x * 100), "%")) +
    scale_x_discrete(expand = expansion(mult = c(0.8, 0.1))) + # Large left margin for labels
    scale_fill_manual(values = c("Logistic_50" = "#56B4E9", "Peak_Fertility" = "#E69F00")) +
    
    theme_bw(base_size = 12) +
    coord_cartesian(clip = "off") + 
    
    labs(title = paste("Methodological Sensitivity:", title_suffix),
         subtitle = "Triangles: >50% Diff | Labels: Left-side only (via ggrepel) | Class: Lines only",
         y = "Relative Difference: (Sen - NoSen)/Sen", x = "", color = "Taxonomic Class") +
    
    theme(legend.position = "bottom",
          plot.title = element_text(face="bold"),
          plot.margin = margin(10, 10, 10, 100), 
          legend.box = "horizontal")
}

# 6. Generate and Save
message("Generating plots...")
p_life <- create_final_plot(df_plot %>% filter(Metric_Type == "Lifespan"), "Lifespan Metrics")
p_lro  <- create_final_plot(df_plot %>% filter(Metric_Type == "LRO"), "LRO Metrics")

ggsave(file.path(output_dir, "Summary_Lifespan_Professional_v5_Repel.png"), p_life, width = 15, height = 7)
ggsave(file.path(output_dir, "Summary_LRO_Professional_v5_Repel.png"), p_lro, width = 15, height = 7)

message("Success! The updated plots (with ggrepel) have been saved to: ", output_dir)