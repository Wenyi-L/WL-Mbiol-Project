<<<<<<< HEAD
## =========================================================
## run_senescence_analysis_batch_excel.R (Updated 03/30/2026)
## Batch analysis for multiple species stored in different
## sheets of a single Excel file.
##
## Updated to use:
## 1. senescence_compute.R (Data prep, MPM build, lifespandist) 
## 2. senescence_plot.R    (ggplot2 visualization)
## 3. LuckFunctions.R      (Moments calculation)
## 4. distTraitCondR.R
## =========================================================

### CH: Comments from Chrissy will be tagged with a triple comment symbol and "CH"
### so that you can ctrl+F to find them (or use Edit->Find in Files to search all
### files in this directory)

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
excel_file  <- "data/Jones2014.xls" 
output_dir  <- "Results/senescence analysis"
species_dir <- file.path(output_dir, "species_plots") # Ensure directory path is defined

# Create directories if they do not exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(species_dir)) dir.create(species_dir, recursive = TRUE)

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

all_results    <- list() # Store global summary results
processed_data <- list() # Store pre-processed data and MPMs


## =========================================================
## PRE-PROCESSING LOOP
## Read data and prepare demography once (KEEPING age=0)
## =========================================================
message("\n==============================")
message("Starting Pre-processing...")
message("==============================")

for (sh in sheets) {
  
  species_name <- sh
  
  meta_raw <- try(
    read_excel(excel_file, sheet = sh, range = "E1", col_names = FALSE),
    silent = TRUE
  )
  study_type <- if (!inherits(meta_raw, "try-error") && nrow(meta_raw) > 0) as.character(meta_raw[[1,1]]) else "Unknown"
  
  ## Read data
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
  
  ## Prepare demographic inputs (Includes age=0 if present in data)
  demog <- try(prepare_demography_data_from_df(dat = df, input_type = "auto"), silent = TRUE)
  
  if (inherits(demog, "try-error")) {
    message("  -> ERROR preparing demography for: ", sh)
    next
  }
  
  # Cap the maximum survival rate to avoid exact 1.0
  demog$sx[demog$sx >= 1] <- 0.9999 
  
  ## Basic validity checks
  if (!any(is.finite(demog$sx))) {
    message("  -> SKIP: All sx values are non-finite for ", sh)
    next
  }
  if (!any(is.finite(demog$fx))) {
    message("  -> SKIP: All fx values are non-finite for ", sh)
    next
  }
  
  # Save to the processed list
  processed_data[[species_name]] <- list(
    study_type = study_type,
    demog      = demog
  )
  message("  -> Pre-processed OK: ", species_name, " (Type: ", study_type, ")")
}


## =========================================================
## LOOP 1: CALCULATIONS
## Calculate matrices, extract sx/fx, and save summaries
## =========================================================
message("\n==============================")
message("Starting Calculation Loop...")
message("==============================")

for (species_name in names(processed_data)) {
  
  message("Calculating models for: ", species_name)
  demog <- processed_data[[species_name]]$demog
  
  res <- try({
    
    ## Export individual sx and fx csv for each age (INCLUDING age=0)
    sx_fx_df <- data.frame(
      age = demog$ages,
      sx  = demog$sx,
      fx  = demog$fx
    )
    write.csv(
      sx_fx_df,
      file = file.path(output_dir, paste0(species_name, "_sx_fx.csv")),
      row.names = FALSE
    )
    
    ## ---------------------------------------------------------
    ## Filter out age=0 EXCLUSIVELY for building MPMs
    ## ---------------------------------------------------------
    build_ages <- demog$ages
    build_sx   <- demog$sx
    build_fx   <- demog$fx
    
    if (length(build_ages) > 0 && build_ages[1] == 0) {
      build_ages <- build_ages[-1]
      build_sx   <- build_sx[-1]
      build_fx   <- build_fx[-1]
    }
    
    ## Build MPMs for four senescence scenarios using FILTERED data
    mpm_sen      <- build_MPM_senescence(ages = build_ages, sx = build_sx, fx = build_fx)
    mpm_no       <- build_MPM_no_senescence(ages = build_ages, sx_senescence = mpm_sen$sx, fx_senescence = mpm_sen$fx, senescence_onset_age = demog$senescence_onset_age)
    mpm_noA_yesR <- build_MPM_no_actuarial_yes_reproductive(ages = build_ages, sx_senescence = mpm_sen$sx, fx_senescence = mpm_sen$fx, senescence_onset_age = demog$senescence_onset_age)
    mpm_yesA_noR <- build_MPM_yes_actuarial_no_reproductive(ages = build_ages, sx_senescence = mpm_sen$sx, fx_senescence = mpm_sen$fx, senescence_onset_age = demog$senescence_onset_age)
    
    ## Save built matrices back into our list so we can plot them directly in Loop 2!
    processed_data[[species_name]]$mpm_sen      <- mpm_sen
    processed_data[[species_name]]$mpm_no       <- mpm_no
    processed_data[[species_name]]$mpm_noA_yesR <- mpm_noA_yesR
    processed_data[[species_name]]$mpm_yesA_noR <- mpm_yesA_noR
    
    ## Summary statistics
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
    summary_df$sheet   <- species_name
    
    summary_df 
  }) # <- end of try()
  
  
  if (!inherits(res, "try-error")) {
    message("  -> OK: completed summary stats for ", species_name)
    
    all_results[[species_name]] <- res
    
    # Save individual CSV for this species
    write.csv(
      res,
      file = file.path(output_dir, paste0(species_name, "_summary_stats.csv")),
      row.names = FALSE
    )
  } else {
    message("  -> ERROR on calculations for species ", species_name)
    message("     ", conditionMessage(attr(res, "condition")))
  }
}

## Combine global summary table
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


## =========================================================
## LOOP 2: PLOTTING
## Use the already processed data and matrices to generate plots
## =========================================================
message("\n==============================")
message("Starting Plotting Loop...")
message("==============================")

for (species_name in names(processed_data)) {
  
  message("Plotting results for: ", species_name)
  
  # Load saved variables from the list
  demog        <- processed_data[[species_name]]$demog
  mpm_sen      <- processed_data[[species_name]]$mpm_sen
  mpm_no       <- processed_data[[species_name]]$mpm_no
  mpm_noA_yesR <- processed_data[[species_name]]$mpm_noA_yesR
  mpm_yesA_noR <- processed_data[[species_name]]$mpm_yesA_noR
  
  # Skip if matrices were not built successfully
  if (is.null(mpm_sen)) {
    message("  -> SKIP Plotting: Missing MPM data for ", species_name)
    next
  }
  
  res <- try({
    
    ## ---------------------------------------------------------
    ## Filter out age=0 for plotting to match MPM dimensions
    ## ---------------------------------------------------------
    plot_ages <- demog$ages
    plot_sx   <- demog$sx
    plot_fx   <- demog$fx
    
    if (length(plot_ages) > 0 && plot_ages[1] == 0) {
      plot_ages <- plot_ages[-1]
      plot_sx   <- plot_sx[-1]
      plot_fx   <- plot_fx[-1]
    }
    
    ## 4a. Survival (senescence vs no-senescence)
    p_surv2 <- plot_survival_models(
      ages         = plot_ages,
      sx_sen       = mpm_sen$sx,
      sx_no        = mpm_no$sx,
      species_name = species_name
    )
    ggsave(
      file.path(species_dir, paste0(species_name, "_survival_sx_compare.png")),
      p_surv2, width = 7, height = 5, dpi = 300, bg = "white"
    )
    
    ## 4a-2. Fecundity (senescence vs no-senescence)
    p_fec2 <- plot_fecundity_models(
      ages         = plot_ages,
      fx_sen       = mpm_sen$fx,
      fx_no        = mpm_no$fx,
      species_name = species_name
    )
    ggsave(
      file.path(species_dir, paste0(species_name, "_fecundity_fx_compare.png")),
      p_fec2, width = 7, height = 5, dpi = 300, bg = "white"
    )
    
    ## 4b. Lifespan distributions (two-model)
    # Using exact matrix method 
    p_LS <- plot_lifespan_distributions(
      U_sen   = mpm_sen$U,
      U_no    = mpm_no$U
    )
    ggsave(
      file.path(species_dir, paste0(species_name, "_lifespan_pmf.png")),
      p_LS, width = 7, height = 5, dpi = 300, bg = "white"
    )
    
    ## 4c-0. Automatic maxClutchSize / maxLRO (for Plotting Grid) 
    ## Note: using filtered plot_fx and plot_sx to match MPM space
    fx   <- plot_fx
    sx   <- plot_sx
    ages <- plot_ages
    
    # Max clutch size
    max_fx = max(as.numeric(fx), na.rm=TRUE)
    maxClutchSize <- qpois(0.999999, max_fx)
    if (maxClutchSize < 5) maxClutchSize <- 5 # Minimum buffer
    
    # Survivorship lx
    n_age <- length(sx)
    lx <- numeric(n_age)
    lx[1] <- 1
    if (n_age > 1) {
      for (i in 2:n_age) lx[i] <- lx[i-1] * sx[i-1]
    }
    
    idx_m <- min(which(fx > 0))
    if (is.infinite(idx_m)) idx_m <- 1
    
    # Expected post-breeding LRO (unconditional estimate for grid sizing)
    expectLRO_post_uncond <- sum(lx[idx_m:n_age] * fx[idx_m:n_age], na.rm = TRUE)
    prob_reach_m <- lx[idx_m]
    
    # Conditional expectation given survival to maturity
    expectLRO_post_cond <- if (prob_reach_m > 0) {
      expectLRO_post_uncond / prob_reach_m
    } else {
      NA_real_
    }
    
    if (is.na(expectLRO_post_cond) || expectLRO_post_cond == 0) expectLRO_post_cond <- 1
    
    # Max LRO support (heuristic for plotting limit)
    maxLRO <- ceiling(3 * expectLRO_post_cond)
    if (maxLRO < 20) maxLRO <- 20
    if (maxLRO > 100) maxLRO <- 100
    
    ## 4c-1. LRO: two-model comparison (analytical iterative, post-breeding)
    p_LRO2 <- plot_LRO_distributions(
      U_sen         = mpm_sen$U,
      U_no          = mpm_no$U,
      F_sen         = mpm_sen$F,
      F_no          = mpm_no$F,
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
      U_sen         = mpm_sen$U,
      U_no          = mpm_no$U,
      U_noA_yesR    = mpm_noA_yesR$U,
      U_yesA_noR    = mpm_yesA_noR$U,
      F_sen         = mpm_sen$F,
      F_no          = mpm_no$F,
      F_noA_yesR    = mpm_noA_yesR$F,
      F_yesA_noR    = mpm_yesA_noR$F,
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
  
  if (inherits(res, "try-error")) {
    message("  -> ERROR on plotting species ", species_name)
    message("     ", conditionMessage(attr(res, "condition")))
  } else {
    message("  -> OK: Plots generated for ", species_name)
  }
=======
## =========================================================
## run_senescence_analysis_batch_excel.R (Updated 03/30/2026)
## Batch analysis for multiple species stored in different
## sheets of a single Excel file.
##
## Updated to use:
## 1. senescence_compute.R (Data prep, MPM build, lifespandist) 
## 2. senescence_plot.R    (ggplot2 visualization)
## 3. LuckFunctions.R      (Moments calculation)
## 4. distTraitCondR.R
## =========================================================

### CH: Comments from Chrissy will be tagged with a triple comment symbol and "CH"
### so that you can ctrl+F to find them (or use Edit->Find in Files to search all
### files in this directory)

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
excel_file  <- "data/Jones2014.xls" 
output_dir  <- "Results/senescence analysis"
species_dir <- file.path(output_dir, "species_plots") # Ensure directory path is defined

# Create directories if they do not exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(species_dir)) dir.create(species_dir, recursive = TRUE)

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

all_results    <- list() # Store global summary results
processed_data <- list() # Store pre-processed data and MPMs


## =========================================================
## PRE-PROCESSING LOOP
## Read data and prepare demography once (KEEPING age=0)
## =========================================================
message("\n==============================")
message("Starting Pre-processing...")
message("==============================")

for (sh in sheets) {
  
  species_name <- sh
  
  meta_raw <- try(
    read_excel(excel_file, sheet = sh, range = "E1", col_names = FALSE),
    silent = TRUE
  )
  study_type <- if (!inherits(meta_raw, "try-error") && nrow(meta_raw) > 0) as.character(meta_raw[[1,1]]) else "Unknown"
  
  ## Read data
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
  
  ## Prepare demographic inputs (Includes age=0 if present in data)
  demog <- try(prepare_demography_data_from_df(dat = df, input_type = "auto"), silent = TRUE)
  
  if (inherits(demog, "try-error")) {
    message("  -> ERROR preparing demography for: ", sh)
    next
  }
  
  # Cap the maximum survival rate to avoid exact 1.0
  demog$sx[demog$sx >= 1] <- 0.9999 
  
  ## Basic validity checks
  if (!any(is.finite(demog$sx))) {
    message("  -> SKIP: All sx values are non-finite for ", sh)
    next
  }
  if (!any(is.finite(demog$fx))) {
    message("  -> SKIP: All fx values are non-finite for ", sh)
    next
  }
  
  # Save to the processed list
  processed_data[[species_name]] <- list(
    study_type = study_type,
    demog      = demog
  )
  message("  -> Pre-processed OK: ", species_name, " (Type: ", study_type, ")")
}


## =========================================================
## LOOP 1: CALCULATIONS
## Calculate matrices, extract sx/fx, and save summaries
## =========================================================
message("\n==============================")
message("Starting Calculation Loop...")
message("==============================")

for (species_name in names(processed_data)) {
  
  message("Calculating models for: ", species_name)
  demog <- processed_data[[species_name]]$demog
  
  res <- try({
    
    ## Export individual sx and fx csv for each age (INCLUDING age=0)
    sx_fx_df <- data.frame(
      age = demog$ages,
      sx  = demog$sx,
      fx  = demog$fx
    )
    write.csv(
      sx_fx_df,
      file = file.path(output_dir, paste0(species_name, "_sx_fx.csv")),
      row.names = FALSE
    )
    
    ## ---------------------------------------------------------
    ## Filter out age=0 EXCLUSIVELY for building MPMs
    ## ---------------------------------------------------------
    build_ages <- demog$ages
    build_sx   <- demog$sx
    build_fx   <- demog$fx
    
    if (length(build_ages) > 0 && build_ages[1] == 0) {
      build_ages <- build_ages[-1]
      build_sx   <- build_sx[-1]
      build_fx   <- build_fx[-1]
    }
    
    ## Build MPMs for four senescence scenarios using FILTERED data
    mpm_sen      <- build_MPM_senescence(ages = build_ages, sx = build_sx, fx = build_fx)
    mpm_no       <- build_MPM_no_senescence(ages = build_ages, sx_senescence = mpm_sen$sx, fx_senescence = mpm_sen$fx, senescence_onset_age = demog$senescence_onset_age)
    mpm_noA_yesR <- build_MPM_no_actuarial_yes_reproductive(ages = build_ages, sx_senescence = mpm_sen$sx, fx_senescence = mpm_sen$fx, senescence_onset_age = demog$senescence_onset_age)
    mpm_yesA_noR <- build_MPM_yes_actuarial_no_reproductive(ages = build_ages, sx_senescence = mpm_sen$sx, fx_senescence = mpm_sen$fx, senescence_onset_age = demog$senescence_onset_age)
    
    ## Save built matrices back into our list so we can plot them directly in Loop 2!
    processed_data[[species_name]]$mpm_sen      <- mpm_sen
    processed_data[[species_name]]$mpm_no       <- mpm_no
    processed_data[[species_name]]$mpm_noA_yesR <- mpm_noA_yesR
    processed_data[[species_name]]$mpm_yesA_noR <- mpm_yesA_noR
    
    ## Summary statistics
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
    summary_df$sheet   <- species_name
    
    summary_df 
  }) # <- end of try()
  
  
  if (!inherits(res, "try-error")) {
    message("  -> OK: completed summary stats for ", species_name)
    
    all_results[[species_name]] <- res
    
    # Save individual CSV for this species
    write.csv(
      res,
      file = file.path(output_dir, paste0(species_name, "_summary_stats.csv")),
      row.names = FALSE
    )
  } else {
    message("  -> ERROR on calculations for species ", species_name)
    message("     ", conditionMessage(attr(res, "condition")))
  }
}

## Combine global summary table
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


## =========================================================
## LOOP 2: PLOTTING
## Use the already processed data and matrices to generate plots
## =========================================================
message("\n==============================")
message("Starting Plotting Loop...")
message("==============================")

for (species_name in names(processed_data)) {
  
  message("Plotting results for: ", species_name)
  
  # Load saved variables from the list
  demog        <- processed_data[[species_name]]$demog
  mpm_sen      <- processed_data[[species_name]]$mpm_sen
  mpm_no       <- processed_data[[species_name]]$mpm_no
  mpm_noA_yesR <- processed_data[[species_name]]$mpm_noA_yesR
  mpm_yesA_noR <- processed_data[[species_name]]$mpm_yesA_noR
  
  # Skip if matrices were not built successfully
  if (is.null(mpm_sen)) {
    message("  -> SKIP Plotting: Missing MPM data for ", species_name)
    next
  }
  
  res <- try({
    
    ## ---------------------------------------------------------
    ## Filter out age=0 for plotting to match MPM dimensions
    ## ---------------------------------------------------------
    plot_ages <- demog$ages
    plot_sx   <- demog$sx
    plot_fx   <- demog$fx
    
    if (length(plot_ages) > 0 && plot_ages[1] == 0) {
      plot_ages <- plot_ages[-1]
      plot_sx   <- plot_sx[-1]
      plot_fx   <- plot_fx[-1]
    }
    
    ## 4a. Survival (senescence vs no-senescence)
    p_surv2 <- plot_survival_models(
      ages         = plot_ages,
      sx_sen       = mpm_sen$sx,
      sx_no        = mpm_no$sx,
      species_name = species_name
    )
    ggsave(
      file.path(species_dir, paste0(species_name, "_survival_sx_compare.png")),
      p_surv2, width = 7, height = 5, dpi = 300, bg = "white"
    )
    
    ## 4a-2. Fecundity (senescence vs no-senescence)
    p_fec2 <- plot_fecundity_models(
      ages         = plot_ages,
      fx_sen       = mpm_sen$fx,
      fx_no        = mpm_no$fx,
      species_name = species_name
    )
    ggsave(
      file.path(species_dir, paste0(species_name, "_fecundity_fx_compare.png")),
      p_fec2, width = 7, height = 5, dpi = 300, bg = "white"
    )
    
    ## 4b. Lifespan distributions (two-model)
    # Using exact matrix method 
    p_LS <- plot_lifespan_distributions(
      U_sen   = mpm_sen$U,
      U_no    = mpm_no$U
    )
    ggsave(
      file.path(species_dir, paste0(species_name, "_lifespan_pmf.png")),
      p_LS, width = 7, height = 5, dpi = 300, bg = "white"
    )
    
    ## 4c-0. Automatic maxClutchSize / maxLRO (for Plotting Grid) 
    ## Note: using filtered plot_fx and plot_sx to match MPM space
    fx   <- plot_fx
    sx   <- plot_sx
    ages <- plot_ages
    
    # Max clutch size
    max_fx = max(as.numeric(fx), na.rm=TRUE)
    maxClutchSize <- qpois(0.999999, max_fx)
    if (maxClutchSize < 5) maxClutchSize <- 5 # Minimum buffer
    
    # Survivorship lx
    n_age <- length(sx)
    lx <- numeric(n_age)
    lx[1] <- 1
    if (n_age > 1) {
      for (i in 2:n_age) lx[i] <- lx[i-1] * sx[i-1]
    }
    
    idx_m <- min(which(fx > 0))
    if (is.infinite(idx_m)) idx_m <- 1
    
    # Expected post-breeding LRO (unconditional estimate for grid sizing)
    expectLRO_post_uncond <- sum(lx[idx_m:n_age] * fx[idx_m:n_age], na.rm = TRUE)
    prob_reach_m <- lx[idx_m]
    
    # Conditional expectation given survival to maturity
    expectLRO_post_cond <- if (prob_reach_m > 0) {
      expectLRO_post_uncond / prob_reach_m
    } else {
      NA_real_
    }
    
    if (is.na(expectLRO_post_cond) || expectLRO_post_cond == 0) expectLRO_post_cond <- 1
    
    # Max LRO support (heuristic for plotting limit)
    maxLRO <- ceiling(3 * expectLRO_post_cond)
    if (maxLRO < 20) maxLRO <- 20
    if (maxLRO > 100) maxLRO <- 100
    
    ## 4c-1. LRO: two-model comparison (analytical iterative, post-breeding)
    p_LRO2 <- plot_LRO_distributions(
      U_sen         = mpm_sen$U,
      U_no          = mpm_no$U,
      F_sen         = mpm_sen$F,
      F_no          = mpm_no$F,
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
      U_sen         = mpm_sen$U,
      U_no          = mpm_no$U,
      U_noA_yesR    = mpm_noA_yesR$U,
      U_yesA_noR    = mpm_yesA_noR$U,
      F_sen         = mpm_sen$F,
      F_no          = mpm_no$F,
      F_noA_yesR    = mpm_noA_yesR$F,
      F_yesA_noR    = mpm_yesA_noR$F,
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
  
  if (inherits(res, "try-error")) {
    message("  -> ERROR on plotting species ", species_name)
    message("     ", conditionMessage(attr(res, "condition")))
  } else {
    message("  -> OK: Plots generated for ", species_name)
  }
>>>>>>> 047bb209d73b65a39fe99f339d3891d2a92d655e
}