## =========================================================
## plot_LRO_relative_impact.R (Checked 02/03/2026)
## (Fixed Variable Names matching CSV)
## =========================================================

rm(list = ls())

## -------------------------------
## 1. Configuration and Loading
## -------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales) 

summary_file <- "Results/senescence analysis/all_species_summary_stats.csv"
output_dir   <- "Results/summary figures"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Define consistent color palette
MODEL_COLORS <- c(
  "Senescence"                      = "#0072B2", # Deep Blue (Baseline)
  "No-senescence"                   = "#D55E00", # Vermilion
  "No-actuarial/Yes-reproductive"   = "#009E73", # Bluish Green
  "Yes-actuarial/No-reproductive"   = "#CC79A7"  # Reddish Purple
)

## -------------------------------
## 2. Data Reading and Cleaning
## -------------------------------
df <- read.csv(summary_file, stringsAsFactors = FALSE)

# Ensure model factor order (Senescence first as baseline)
df$model <- factor(df$model, levels = c(
  "Senescence", 
  "No-senescence", 
  "No-actuarial/Yes-reproductive", 
  "Yes-actuarial/No-reproductive"
))

## =========================================================
## Data Preparation: Relative % Change vs Senescence
## =========================================================

# Calculate percentage change relative to Senescence
df_relative <- df %>%
  select(species, model, mean_LRO, var_LRO, skew_LRO) %>%
  pivot_wider(names_from = model, values_from = c(mean_LRO, var_LRO, skew_LRO)) %>%
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

## =========================================================
## Plotting
## =========================================================

p2 <- ggplot(df_relative, aes(x = Model, y = pct_change, fill = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") + 
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.1, fill = "transparent", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, size = 1.2, alpha = 0.4) +
  facet_wrap(~Metric, scales = "free_y") +
  # Pseudo-log scale for handling large variations
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 1),
    breaks = c(-50, 0, 10, 100, 1000, 10000, 1e5, 1e6),
    labels = scales::comma_format()
  ) +
  scale_fill_manual(values = MODEL_COLORS) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.direction = "vertical",
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(
    title = "Relative Impact of Models vs. Senescence Baseline",
    subtitle = "Values show % change relative to the Senescence model (0% = No change)",
    x = NULL,
    y = "% Change (Pseudo-Log Scale)",
    fill = "Comparison Model"
  )

# Save the plot
ggsave(file.path(output_dir, "4models_Relative_Impact_LogScale.png"), p2, width = 10, height = 7, dpi = 300)

message("Relative impact plot saved to: ", output_dir)