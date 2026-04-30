## =========================================================
## run_phyloANOVA.R 
##
## Purpose: 
## 1. Build and Plot Phylogeny (Original Tree, NO modifications!)
## 2. Clean data (remove brackets, underscores, duplicates)
## 3. Run PGLS (corPagel) - lambda estimated via ML, following 
##    Harmon (2015) PGLS tutorial: 
##    https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/
## =========================================================

# Load required libraries
library(readxl)
library(rotl)
library(ape)
library(nlme)

##-----------------------------------------
## 1. Build and Plot Phylogenetic Tree
##-----------------------------------------
sheet_names <- excel_sheets("data/Jones2014.xls")
my_taxa <- unique(gsub(" \\(.*\\)", "", sheet_names))
resolved_names <- tnrs_match_names(names = my_taxa)
my_tree <- tol_induced_subtree(resolved_names$ott_id)

my_tree$tip.label <- gsub('_ott.*', '', my_tree$tip.label) 
my_tree$tip.label <- gsub('_', ' ', my_tree$tip.label)

# Assigns relative branch lengths
my_tree_brlen <- compute.brlen(my_tree, method = "Grafen")

# Plot the tree
plot(my_tree_brlen, cex = 0.8, main = "Phylogenetic Tree of the Species")

##-----------------------------------------
## 2. Load and Clean Dataset
##-----------------------------------------
data_all <- read.csv("Results/senescence analysis/all_species_summary_stats.csv", stringsAsFactors = FALSE)

# Clean species names (Your original logic)
data_all$species <- gsub("\\s*\\(.*?\\)", "", data_all$species)
data_all$species <- gsub("_", " ", data_all$species)

##-----------------------------------------
## 3. Define PGLS Function (corPagel - estimates lambda via ML)
##-----------------------------------------
get_pgls_row <- function(metric_name, comp_model_name) {
  
  # Subset and deduplicate data
  data_sen <- data_all[data_all$model == "Senescence", ]
  data_comp <- data_all[data_all$model == comp_model_name, ]
  
  data_sen <- data_sen[!duplicated(data_sen$species), ]
  data_comp <- data_comp[!duplicated(data_comp$species), ]
  
  # Align species
  common_spp <- intersect(data_sen$species, data_comp$species)
  data_sen <- data_sen[match(common_spp, data_sen$species), ]
  data_comp <- data_comp[match(common_spp, data_comp$species), ]
  
  # Calculate differences
  diff_values <- data_sen[[metric_name]] - data_comp[[metric_name]]
  
  df_diff <- data.frame(
    species = as.character(data_sen$species),
    diff_val = as.numeric(diff_values),
    stringsAsFactors = FALSE
  )
  
  # Match with tree tips directly
  df_diff <- df_diff[df_diff$species %in% my_tree_brlen$tip.label, ]
  
  # Prune tree to match data safely
  pruned_tree <- keep.tip(my_tree_brlen, df_diff$species)
  
  # Rescale branch lengths (Harmon tutorial fix)
  tempTree <- pruned_tree
  tempTree$edge.length <- tempTree$edge.length * 1000
  
  # gls control: use optim (more robust than default nlminb for 
  # intercept-only models), allow more iterations, looser tolerance.
  gls_ctrl <- list(opt = "optim", maxIter = 200, msMaxIter = 200,
                   tolerance = 1e-6, msTol = 1e-6)
  
  # Try several starting values for lambda. tryCatch handles BOTH 
  # errors (NA/NaN/Inf) and convergence failures (false convergence).
  pgls_model <- NULL
  for (lam_init in c(0.5, 0.3, 0.1, 0.7, 0)) {
    pgls_model <- tryCatch(
      gls(diff_val ~ 1, 
          data = df_diff, 
          correlation = corPagel(value = lam_init, phy = tempTree, 
                                 form = ~species, fixed = FALSE), 
          method = "ML",
          control = gls_ctrl),
      error = function(e) NULL,
      warning = function(w) NULL
    )
    if (!is.null(pgls_model)) break
  }
  
  if (is.null(pgls_model)) {
    return(data.frame(
      Metric     = metric_name,
      Comparison = paste("Senescence vs", comp_model_name),
      Lambda     = NA,
      P_Value    = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  # Extract intercept p-value
  p_val <- summary(pgls_model)$tTable[1, 4]
  
  # Extract estimated lambda from the fitted correlation structure
  lambda_val <- coef(pgls_model$modelStruct$corStruct, unconstrained = FALSE)[[1]]
  
  return(data.frame(
    Metric     = metric_name,
    Comparison = paste("Senescence vs", comp_model_name),
    Lambda     = round(lambda_val, 4),
    P_Value    = round(p_val, 4),
    stringsAsFactors = FALSE
  ))
}

##-----------------------------------------
## 4. Run Analysis Loop & Collect Results
##-----------------------------------------
results_list <- list()

# --- Part A: Lifespan Metrics ---
lifespan_metrics <- c("mean_lifespan", "var_lifespan", "skew_lifespan")
for (m in lifespan_metrics) {
  results_list[[length(results_list) + 1]] <- get_pgls_row(m, "No-senescence")
}

# --- Part B: LRO Metrics ---
lro_metrics <- c("mean_LRO", "var_LRO", "skew_LRO")
other_models <- c("No-senescence", 
                  "No-actuarial/Yes-reproductive", 
                  "Yes-actuarial/No-reproductive")

for (m in lro_metrics) {
  for (target in other_models) {
    results_list[[length(results_list) + 1]] <- get_pgls_row(m, target)
  }
}

final_results_table <- do.call(rbind, results_list)

##-----------------------------------------
## 5. View and Export Results
##-----------------------------------------
print(final_results_table)
write.csv(final_results_table, "Results/PGLS_Analysis_Summary_Pagel.csv", row.names = FALSE)