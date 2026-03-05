## =========================================================
## run_senescence_analysis_batch_excel.R (Updated 03/03/2026)
## Batch analysis for multiple species stored in different
## sheets of a single Excel file.
##
## Updated to use:
## 1. senescence_compute.R (Data prep, MPM build, lifespandist)
## 2. senescence_plots.R   (ggplot2 visualization)
## 3. LuckFunctions.R      (Moments calculation)
## 4. distTraitCondR.R
## =========================================================

rm(list = ls())

# Load necessary libraries
library(readxl)
library(ggplot2)
library(minpack.lm)
library(Matrix)

# exactLTRE is required by LuckFunctions.R
if (!requireNamespace("exactLTRE", quietly = TRUE)) {
  stop("Package 'exactLTRE' is required. Please install or ensure it is available.")
}
library(exactLTRE)

## -------------------------------
## Configuration
## -------------------------------
excel_file <- "data/Jones2014.xls" 
output_dir <- "Results/senescence analysis"



## -------------------------------
## Load required functions (UPDATED)
## -------------------------------
# 1. Load Math/Moments functions first
source("code/LuckFunctions.R")

# 2. Load Compute functions 
source("code/senescence_functions.R")
source("code/distTraitCondR.R")
source("code/megamatrixFunctionsLRO.R")

# 3. Load Plotting functions (Survival, Fecundity, LS, LRO plots)
source("code/senescence_plot.R")

## -------------------------------
## Read all sheet names
## -------------------------------
sheets <- excel_sheets(excel_file)
message("Found sheets: ", paste(sheets, collapse = ", "))

all_results <- list() #store results


## -------------------------------
## Loop over sheets (Calculations)
## -------------------------------


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
    mpm_no <- build_MPM_no_senescence(ages = demog$ages, sx_senescence = mpm_sen$sx, fx_senescence = mpm_sen$fx, senescence_onset_age = demog$senescence_onset_age)
    mpm_noA_yesR <- build_MPM_no_actuarial_yes_reproductive(ages = demog$ages, sx_senescence = mpm_sen$sx, fx_senescence = mpm_sen$fx, senescence_onset_age = demog$senescence_onset_age)
    mpm_yesA_noR <- build_MPM_yes_actuarial_no_reproductive(ages = demog$ages, sx_senescence = mpm_sen$sx, fx_senescence = mpm_sen$fx, senescence_onset_age = demog$senescence_onset_age)
    
    ## 3. Summary statistics
    summary_df <- compute_summary_table(
      U_sen        = mpm_sen$U,
      U_no         = mpm_no$U,
      U_noA_yesR   = mpm_noA_yesR$U,
      U_yesA_noR   = mpm_yesA_noR$U,
      F_sen        = mpm_sen$F,
      F_no         = mpm_no$F,
      F_noA_yesR   = mpm_noA_yesR$F,
      F_yesA_noR   = mpm_yesA_noR$F,
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
    all_results[[species_name]] <- res
    
    # 2. Save individual CSV for this species
    write.csv(
      res,
      file = file.path(output_dir, paste0(species_name, "_summary_stats.csv")),
      row.names = FALSE
    )
  } else {
    message("  -> ERROR on species ", species_name)
  }
  
} # <--- END of the for loop! 


## -------------------------------
## Combine global summary table
## -------------------------------
if (length(all_results) == 0) {
  warning("No species produced valid summary results. No global summary created.")
} else {
  global_df <- do.call(rbind, all_results)
  write.csv(
    global_df,
    file = file.path(output_dir, "all_species_summary_stats.csv"),
    row.names = FALSE
  )
  message("\nGlobal summary saved successfully.")
}


    
    ## -------------------------------
    ## ---- Plotting ----
    ## -------------------------------
  

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
    
    # Source shield: Cap the maximum survival rate to avoid exact 1.0
    demog$sx[demog$sx >= 1] <- 0.9999 
    
    ## Basic validity checks
    if (!any(is.finite(demog$sx))) stop("All sx values are non-finite")
    if (!any(is.finite(demog$fx))) stop("All fx values are non-finite")
    
    ## 2. Build MPMs for four senescence scenarios
    mpm_sen <- build_MPM_senescence(ages = demog$ages, sx = demog$sx, fx = demog$fx)
    mpm_no <- build_MPM_no_senescence(ages = demog$ages, sx_senescence = mpm_sen$sx, fx_senescence = mpm_sen$fx, senescence_onset_age = demog$senescence_onset_age)
    mpm_noA_yesR <- build_MPM_no_actuarial_yes_reproductive(ages = demog$ages, sx_senescence = mpm_sen$sx, fx_senescence = mpm_sen$fx, senescence_onset_age = demog$senescence_onset_age)
    mpm_yesA_noR <- build_MPM_yes_actuarial_no_reproductive(ages = demog$ages, sx_senescence = mpm_sen$sx, fx_senescence = mpm_sen$fx, senescence_onset_age = demog$senescence_onset_age)
    
    if (!dir.exists(species_dir)) dir.create(species_dir, recursive = TRUE)
    
    ## 4a. Survival (senescence vs no-senescence)
    p_surv2 <- plot_survival_models(
      ages        = demog$ages,
      sx_sen      = mpm_sen$sx,
      sx_no       = mpm_no$sx,
      species_name = species_name
    )
    ggsave(
      file.path(species_dir, paste0(species_name, "_survival_sx_compare.png")),
      p_surv2, width = 7, height = 5, dpi = 300, bg = "white"
    )
    
    ## 4a-2. Fecundity (senescence vs no-senescence)
    p_fec2 <- plot_fecundity_models(
      ages         = demog$ages,
      fx_sen       = mpm_sen$fx,
      fx_no        = mpm_no$fx,
      species_name = species_name
    )
    ggsave(
      file.path(species_dir, paste0(species_name, "_fecundity_fx_compare.png")),
      p_fec2, width = 7, height = 5, dpi = 300, bg = "white"
    )
    
    ## 4b. Lifespan distributions (two-model)
    # Using exact matrix method (from senescence_compute)
    p_LS <- plot_lifespan_distributions(
      U_sen   = mpm_sen$U,
      U_no    = mpm_no$U
    )
    ggsave(
      file.path(species_dir, paste0(species_name, "_lifespan_pmf.png")),
      p_LS, width = 7, height = 5, dpi = 300, bg = "white"
    )
    
    ## 4c-0
    ## ---- Automatic maxClutchSize / maxLRO (for Plotting Grid) ----
    fx <- demog$fx
    sx <- demog$sx
    ages <- demog$ages
    
    # max clutch size
    max_fx = max(as.numeric(fx), na.rm=TRUE)
    maxClutchSize <- qpois(0.999999, max_fx)
    if (maxClutchSize < 5) maxClutchSize <- 5 # Minimum buffer
    
    # survivorship lx
    n_age <- length(sx)
    lx <- numeric(n_age)
    lx[1] <- 1
    if (n_age > 1) {
      for (i in 2:n_age) lx[i] <- lx[i-1] * sx[i-1]
    }
    
    idx_m <- min(which(fx > 0))
    if (is.infinite(idx_m)) idx_m <- 1
    
    # expected post-breeding LRO (unconditional estimate for grid sizing)
    expectLRO_post_uncond <- sum(lx[idx_m:n_age] * fx[idx_m:n_age], na.rm = TRUE)
    prob_reach_m <- lx[idx_m]
    
    # conditional expectation given survival to maturity
    expectLRO_post_cond <- if (prob_reach_m > 0) {
      expectLRO_post_uncond / prob_reach_m
    } else {
      NA_real_
    }
    
    if (is.na(expectLRO_post_cond) || expectLRO_post_cond == 0) expectLRO_post_cond <- 1
    
    # max LRO support (heuristic for plotting limit)
    maxLRO <- ceiling(3 * expectLRO_post_cond)
    if (maxLRO < 20) maxLRO <- 20
    if (maxLRO > 100) maxLRO <- 100
    
    ## 4c-1. LRO: two-model comparison (analytical iterative, post-breeding)
    # This now uses the efficient calcDistLRO_iterative inside plotting function
    p_LRO2 <- plot_LRO_distributions(
      U_sen   = mpm_sen$U,
      U_no    = mpm_no$U,
      F_sen   = mpm_sen$F,
      F_no    = mpm_no$F,
      maxClutchSize = maxClutchSize,    
      maxLRO        = maxLRO,
      include_zero  = TRUE
    )
    
    ggsave(
      file.path(species_dir, paste0(species_name, "_LRO_pmf_2models.png")),
      p_LRO2, width = 7, height = 5, dpi = 300, bg = "white"
    )
    
    ## 4c-2. LRO: four-model comparison (analytical iterative, post-breeding)
    p_LRO4 <- plot_LRO_distributions_4(
      U_sen        = mpm_sen$U,
      U_no         = mpm_no$U,
      U_noA_yesR   = mpm_noA_yesR$U,
      U_yesA_noR   = mpm_yesA_noR$U,
      F_sen        = mpm_sen$F,
      F_no         = mpm_no$F,
      F_noA_yesR   = mpm_noA_yesR$F,
      F_yesA_noR   = mpm_yesA_noR$F,
      maxClutchSize = maxClutchSize,
      maxLRO        = maxLRO,
      include_zero  = TRUE
    )
    
    ## Save four-model figure for this species
    ggsave(
      file.path(species_dir, paste0(species_name, "_LRO_pmf_4models.png")),
      p_LRO4, width = 7, height = 5, dpi = 300, bg = "white"
    )
  }) # <--- end of try()
    
  
  ## -----------------------------------
  ## Error handling
  ## -----------------------------------
  if (inherits(res, "try-error")) {
    message("  -> ERROR on species ", species_name)
    message("     ", conditionMessage(attr(res, "condition")))
    next
  }
  } # <--- end of loop
  
  
  
  
  
  
  
  
  
  
  
