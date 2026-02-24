## =========================================================
## run_senescence_analysis_batch_excel.R
## Batch analysis for multiple species stored in different
## sheets of a single Excel file.
##
## Updated to use:
## 1. senescence_compute.R (Data prep, MPM build, Iterative LRO)
## 2. senescence_plots.R   (ggplot2 visualization)
## 3. LuckFunctions.R      (Moments calculation)
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
excel_file <- "Jones2014.xls"  # <-- update this if needed
output_dir <- "output_batch"
estimate_tail <- TRUE
tail_k        <- 2
qcut          <- 0.999

## Maturity detection options 
maturity_method <- "logistic50"   # "logistic50" or "absolute"
auto_min_fx     <- 1e-8
auto_span       <- 0.5            # loess span for smoothing before logistic fit
plot_maturity   <- TRUE           # whether to draw/save maturity detection plot

# Function for naming maturity plots
maturity_plot_filename <- function(species_name) {
  file.path(output_dir, species_name, paste0(species_name, "_maturity_logistic50.png"))
}

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## -------------------------------
## Load required functions (UPDATED)
## -------------------------------
# 1. Load Math/Moments functions first
source("LuckFunctions.R")

# 2. Load Compute functions (Data prep, MPM, Iterative LRO, Summary Table)
#    Note: This file contains the efficient 'calcDistLRO_iterative' 
#    and the wrapper 'calcDistLROPostBreedingNoEnv'
source("senescence_functions.R")

# 3. Load Plotting functions (Survival, Fecundity, LS, LRO plots)
source("senescence_plot.R")

## -------------------------------
## Read all sheet names
## -------------------------------
sheets <- excel_sheets(excel_file)
message("Found sheets: ", paste(sheets, collapse = ", "))

all_results <- list()

## -------------------------------
## Prepare a global PDF for 4-model LRO plots
## -------------------------------
global_LRO_pdf <- file.path(output_dir, "all_species_LRO_4models.pdf")
pdf(global_LRO_pdf, width = 7, height = 5)

## -------------------------------
## Loop over sheets
## -------------------------------
for (sh in sheets) {
  
  message("\n==============================")
  message("Processing sheet: ", sh)
  message("==============================")
  
  species_name <- sh
  
  ## ---- Read data (first row contains column names) ----
  df_raw <- try(
    read_excel(
      excel_file,
      sheet     = sh,
      skip      = 1,
      col_names = TRUE
    ),
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
    ## ---- 1. Prepare demographic inputs ----
    # prepare plot path for maturity (per-species)
    species_dir <- file.path(output_dir, species_name)
    if (!dir.exists(species_dir)) dir.create(species_dir, recursive = TRUE)
    plot_path_for_species <- maturity_plot_filename(species_name)
    
    demog <- prepare_demography_data_from_df(
      dat              = df,
      input_type       = "auto",
      maturity_age     = "auto",
      estimate_tail    = estimate_tail,
      tail_k           = tail_k,
      maturity_method  = maturity_method,
      auto_min_fx      = auto_min_fx,
      auto_consecutive = 1,
      auto_prop_of_max = 0.05,
      auto_span        = auto_span,
      plot_maturity    = plot_maturity,
      plot_path        = if (isTRUE(plot_maturity)) plot_path_for_species else NULL
    )
    
    ## Basic validity checks
    if (!any(is.finite(demog$sx))) stop("All sx values are non-finite")
    if (!any(is.finite(demog$fx))) stop("All fx values are non-finite")
    
    ## ---- 2. Build MPMs for four senescence scenarios ----
    mpm_sen <- build_MPM_senescence(
      ages = demog$ages,
      sx   = demog$sx,
      fx   = demog$fx
    )
    
    mpm_no <- build_MPM_no_senescence(
      ages          = demog$ages,
      sx_senescence = mpm_sen$sx,
      fx_senescence = mpm_sen$fx,
      maturity_age  = demog$maturity_age
    )
    
    mpm_noA_yesR <- build_MPM_no_actuarial_yes_reproductive(
      ages          = demog$ages,
      sx_senescence = mpm_sen$sx,
      fx_senescence = mpm_sen$fx,
      maturity_age  = demog$maturity_age
    )
    
    mpm_yesA_noR <- build_MPM_yes_actuarial_no_reproductive(
      ages          = demog$ages,
      sx_senescence = mpm_sen$sx,
      fx_senescence = mpm_sen$fx,
      maturity_age  = demog$maturity_age
    )
    
    ## Mixing distributions (needed for LuckFunctions)
    # Note: mixing_distro is in LuckFunctions.R
    mix_sen      <- mixing_distro(mpm_sen$A,      mpm_sen$F)
    mix_no       <- mixing_distro(mpm_no$A,       mpm_no$F)
    mix_noA_yesR <- mixing_distro(mpm_noA_yesR$A, mpm_noA_yesR$F)
    mix_yesA_noR <- mixing_distro(mpm_yesA_noR$A, mpm_yesA_noR$F)
    
    ## ---- 3. Summary statistics (Using LuckFunctions via senescence_compute) ----
    # This computes moments (Mean, Var, Skew) without simulating distributions
    summary_df <- compute_summary_table(
      U_sen        = mpm_sen$U,
      U_no         = mpm_no$U,
      U_noA_yesR   = mpm_noA_yesR$U,
      U_yesA_noR   = mpm_yesA_noR$U,
      F_sen        = mpm_sen$F,
      F_no         = mpm_no$F,
      F_noA_yesR   = mpm_noA_yesR$F,
      F_yesA_noR   = mpm_yesA_noR$F,
      mix_sen      = mix_sen,
      mix_no       = mix_no,
      mix_noA_yesR = mix_noA_yesR,
      mix_yesA_noR = mix_yesA_noR,
      repro_var    = "Poisson"  # LuckFunctions option
    )
    
    summary_df$species <- species_name
    summary_df$sheet   <- sh
    
    ## -------------------------------
    ## ---- 4. Plotting routines ----
    ## -------------------------------
    species_dir <- file.path(output_dir, species_name)
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
      U_no    = mpm_no$U,
      mix_sen = mix_sen,
      mix_no  = mix_no,
      qcut    = qcut
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
    k_clutch <- 2
    maxClutchSize <- ceiling(
      max(fx + k_clutch * sqrt(fx), na.rm = TRUE)
    )
    if (maxClutchSize < 5) maxClutchSize <- 5 # Minimum buffer
    
    # survivorship lx
    n_age <- length(sx)
    lx <- numeric(n_age)
    lx[1] <- 1
    if (n_age > 1) {
      for (i in 2:n_age) lx[i] <- lx[i-1] * sx[i-1]
    }
    
    # maturity index
    mat_age <- demog$maturity_age
    idx_m <- which.min(abs(ages - mat_age))
    
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
    maxLRO <- ceiling(5 * expectLRO_post_cond)
    if (maxLRO < 20) maxLRO <- 20
    
    ## 4c-1. LRO: two-model comparison (analytical iterative, post-breeding)
    # This now uses the efficient calcDistLRO_iterative inside plotting function
    p_LRO2 <- plot_LRO_distributions(
      U_sen   = mpm_sen$U,
      U_no    = mpm_no$U,
      F_sen   = mpm_sen$F,
      F_no    = mpm_no$F,
      mix_sen = mix_sen,
      mix_no  = mix_no,
      maxClutchSize = maxClutchSize,    
      maxLRO        = maxLRO,
      include_zero  = FALSE,
      qcut          = qcut
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
      mix_sen      = mix_sen,
      mix_no       = mix_no,
      mix_noA_yesR = mix_noA_yesR,
      mix_yesA_noR = mix_yesA_noR,
      maxClutchSize = maxClutchSize,
      maxLRO        = maxLRO,
      include_zero  = FALSE,
      qcut          = qcut
    )
    
    ## Save four-model figure for this species
    ggsave(
      file.path(species_dir, paste0(species_name, "_LRO_pmf_4models.png")),
      p_LRO4, width = 7, height = 5, dpi = 300, bg = "white"
    )
    
    ## Append to global multi-page PDF
    print(p_LRO4)
    
    summary_df
  }, silent = TRUE)
  
  ## -----------------------------------
  ## Error handling
  ## -----------------------------------
  if (inherits(res, "try-error")) {
    message("  -> ERROR on species ", species_name)
    message("     ", conditionMessage(attr(res, "condition")))
    next
  }
  
  ## -----------------------------------
  ## Successful species
  ## -----------------------------------
  message("  -> OK: added species ", species_name)
  all_results[[length(all_results) + 1]] <- res
  
  ## Save per-species summary table
  write.csv(
    res,
    file = file.path(output_dir, paste0(species_name, "_summary_stats.csv")),
    row.names = FALSE
  )
}

## Close the global 4-model LRO PDF
dev.off()

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
  message("\nGlobal summary saved.")
}