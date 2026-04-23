## =========================================================
## run_senescence_analysis_single_species.R
## Single-species run:
## CSV file -> MPM -> summary & plots
## Updated to use new modular architecture (Iterative LRO)
## =========================================================

rm(list = ls())

## -------------------------
## 0. Configuration
## -------------------------
species_name  <- "Viburnum furcatum"
input_csv     <- "Viburnum furcatum.csv"   # change to your CSV
output_dir    <- "output"

## input_type:
##   "auto"   : detect counts vs rates
##   "counts" : x, Nx, noffspring
##   "rates"  : x, lx, fert.mx
input_type    <- "auto"

## maturity_age:
##   numeric : fixed
##   "auto"  : detect from fx (uses logistic50 by default here)
maturity_age  <- "auto"

estimate_tail <- TRUE
tail_k        <- 2
qcut          <- 0.999  # used only for plotting/truncation heuristics

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## -------------------------
## 1. Load packages & function files (UPDATED)
## -------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Please install ggplot2: install.packages('ggplot2')")
}
if (!requireNamespace("minpack.lm", quietly = TRUE)) {
  message("Package 'minpack.lm' not found. If logistic fitting fails, please install minpack.lm.")
}
# exactLTRE is required by LuckFunctions.R
if (!requireNamespace("exactLTRE", quietly = TRUE)) {
  stop("Package 'exactLTRE' is required. Please install or ensure it is available.")
}

library(ggplot2)
library(minpack.lm)
library(Matrix)
library(exactLTRE)

# Source updated helper function files
# 1. Load Math/Moments functions first
source("LuckFunctions.R")

# 2. Load Compute functions (Data prep, MPM build, Iterative LRO, Summary Table)
source("senescence_functions.R")

# 3. Load Plotting functions (Survival, Fecundity, LS, LRO plots)
source("senescence_plot.R")

## -------------------------
## 2. Data preparation + MPM construction
## -------------------------

## 2.1 Demography data
maturity_plot_path <- file.path(output_dir,
                                paste0(gsub("[[:punct:]\\s]+", "_", species_name),
                                       "_maturity_detection_logistic50.png"))

demog <- prepare_demography_data(
  file             = input_csv,
  input_type       = input_type,
  maturity_age     = maturity_age,
  estimate_tail    = estimate_tail,
  tail_k           = tail_k,
  maturity_method  = "logistic50",
  auto_min_fx      = 1e-8,
  auto_span        = 0.5,
  plot_maturity    = TRUE,
  plot_path        = maturity_plot_path
)

if (!any(is.finite(demog$sx))) stop("All sx values are non-finite")
if (!any(is.finite(demog$fx))) stop("All fx values are non-finite")


## 2.2 Senescence MPM
mpm_sen <- build_MPM_senescence(
  ages = demog$ages,
  sx   = demog$sx,
  fx   = demog$fx
)

## 2.3 No-senescence MPM
mpm_no <- build_MPM_no_senescence(
  ages          = demog$ages,
  sx_senescence = mpm_sen$sx,
  fx_senescence = mpm_sen$fx,
  maturity_age  = demog$maturity_age
)

## 2.3b & 2.3c
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

## 2.5 Mixing distributions (needed for LuckFunctions)
mixdist_sen      <- mixing_distro(mpm_sen$A,      mpm_sen$F)
mixdist_no       <- mixing_distro(mpm_no$A,       mpm_no$F)
mixdist_noA_yesR <- mixing_distro(mpm_noA_yesR$A, mpm_noA_yesR$F)
mixdist_yesA_noR <- mixing_distro(mpm_yesA_noR$A, mpm_yesA_noR$F)

cat("Sum(mixdist_senescence) =", sum(mixdist_sen),
    " | Sum(mixdist_no_senescence) =", sum(mixdist_no), "\n")


## -------------------------
## 3. Plotting routines
## -------------------------
# Create species output directory if it doesn't match root output
species_file_prefix <- gsub("[[:punct:]\\s]+", "_", species_name)

## 3.1 Survival curves: senescence vs no-senescence
p_surv <- plot_survival_models(
  ages        = demog$ages,
  sx_sen      = mpm_sen$sx,
  sx_no       = mpm_no$sx,
  species_name = species_name
)
print(p_surv)
ggsave(
  filename = file.path(output_dir, paste0(species_file_prefix, "_survival_sen_vs_no.png")),
  plot = p_surv, width = 7, height = 5, dpi = 300, bg = "white"
)

## 3.2 fx (fertility) plot：senescence vs no-senescence
p_fecu <- plot_fecundity_models(
  ages        = demog$ages,
  fx_sen      = mpm_sen$fx,
  fx_no       = mpm_no$fx,
  species_name = species_name
)
print(p_fecu)
ggsave(
  filename = file.path(output_dir, paste0(species_file_prefix, "_fecundity_sen_vs_no.png")),
  plot = p_fecu, width = 7, height = 5, dpi = 300, bg = "white"
)

## 3.3 Lifespan distributions (exact)
# Note: plot_lifespan_distributions uses exact_lifespan_pmf from senescence_compute.R
p_LS <- plot_lifespan_distributions(
  U_sen   = mpm_sen$U,
  U_no    = mpm_no$U,
  mix_sen = mixdist_sen,
  mix_no  = mixdist_no,
  qcut    = qcut
)
print(p_LS)
ggsave(
  filename = file.path(output_dir, paste0(species_file_prefix, "_lifespan_pmf.png")),
  plot = p_LS, width = 7, height = 5, dpi = 300, bg = "white"
)


## -------------------------
## 4. Deterministic LRO distribution (post-breeding, no env)
##    Using new Iterative Convolution method
## -------------------------

# ---- Automatic maxClutchSize / maxLRO (Configuration) ----
fx <- demog$fx
sx <- demog$sx
ages <- demog$ages

# choose k for maxClutchSize calculation (used to estimate an upper bound)
k_clutch <- 2
maxClutchSize <- ceiling(
  max(fx + k_clutch * sqrt(fx), na.rm = TRUE)
)
if (maxClutchSize < 5) maxClutchSize <- 5 # Minimum buffer

# find the index of maturity
mat_age <- demog$maturity_age
idx_m <- which.min(abs(ages - mat_age))
cat("maturity age ~", ages[idx_m], " (index =", idx_m, ")\n")

# survivorship lx for estimation
n <- length(sx)
lx <- numeric(n); lx[1] <- 1
if (n > 1) for (i in 2:n) lx[i] <- lx[i-1] * sx[i-1]

# unconditional expected LRO estimate
expectLRO_post_uncond <- sum(lx[idx_m:n] * fx[idx_m:n], na.rm = TRUE)

# conditional expectation given surviving to maturity
prob_reach_m <- lx[idx_m]
expectLRO_post_cond <- if (prob_reach_m > 0) expectLRO_post_uncond / prob_reach_m else NA_real_

if (is.na(expectLRO_post_cond) || expectLRO_post_cond == 0) expectLRO_post_cond <- 1
maxLRO <- ceiling(5 * expectLRO_post_cond)
if (maxLRO < 20) maxLRO <- 20 # Minimum buffer

cat("Estimated maxLRO for plotting:", maxLRO, "\n")

# ---- Plotting using senescence_plots.R functions ----
# These functions now use calcDistLRO_iterative internally via calcDistLROPostBreedingNoEnv

# 4.1 LRO: 2-model comparison
p_LRO_2 <- plot_LRO_distributions(
  U_sen   = mpm_sen$U,
  U_no    = mpm_no$U,
  F_sen   = mpm_sen$F,
  F_no    = mpm_no$F,
  mix_sen = mixdist_sen,
  mix_no  = mixdist_no,
  maxClutchSize = maxClutchSize,
  maxLRO        = maxLRO,
  include_zero  = FALSE,
  qcut          = qcut
)
print(p_LRO_2)
ggsave(
  filename = file.path(output_dir, paste0(species_file_prefix, "_LRO_pmf_postbreeding.png")),
  plot = p_LRO_2, width = 7, height = 5, dpi = 300, bg = "white"
)

# 4.2 LRO: 4-model comparison
p_LRO_4 <- plot_LRO_distributions_4(
  U_sen        = mpm_sen$U,
  U_no         = mpm_no$U,
  U_noA_yesR   = mpm_noA_yesR$U,
  U_yesA_noR   = mpm_yesA_noR$U,
  F_sen        = mpm_sen$F,
  F_no         = mpm_no$F,
  F_noA_yesR   = mpm_noA_yesR$F,
  F_yesA_noR   = mpm_yesA_noR$F,
  mix_sen      = mixdist_sen,
  mix_no       = mixdist_no,
  mix_noA_yesR = mixdist_noA_yesR,
  mix_yesA_noR = mixdist_yesA_noR,
  maxClutchSize = maxClutchSize,
  maxLRO        = maxLRO,
  include_zero  = FALSE,
  qcut          = qcut
)
print(p_LRO_4)
ggsave(
  filename = file.path(output_dir, paste0(species_file_prefix, "_LRO_pmf_4models_postbreeding.png")),
  plot = p_LRO_4, width = 8, height = 5, dpi = 300, bg = "white"
)


## -------------------------
## 5. Compute summary table
##    Using LuckFunctions via compute_summary_table
## -------------------------
summary_df <- compute_summary_table(
  U_sen        = mpm_sen$U,
  U_no         = mpm_no$U,
  U_noA_yesR   = mpm_noA_yesR$U,
  U_yesA_noR   = mpm_yesA_noR$U,
  F_sen        = mpm_sen$F,
  F_no         = mpm_no$F,
  F_noA_yesR   = mpm_noA_yesR$F,
  F_yesA_noR   = mpm_yesA_noR$F,
  mix_sen      = mixdist_sen,
  mix_no       = mixdist_no,
  mix_noA_yesR = mixdist_noA_yesR,
  mix_yesA_noR = mixdist_yesA_noR,
  repro_var    = "Poisson"
)

summary_df$species <- species_name
summary_df <- summary_df[, c("species", "model", "mean_lifespan", "var_lifespan", "skew_lifespan", "mean_LRO", "var_LRO", "skew_LRO")]

print(summary_df)

write.csv(
  summary_df,
  file = file.path(output_dir, paste0(species_file_prefix, "_summary_stats.csv")),
  row.names = FALSE
)

# Optional: Append to a global file if you are running multiple single runs
global_file <- file.path(output_dir, "all_species_summary_stats.csv")
if (file.exists(global_file)) {
  old_df <- read.csv(global_file, stringsAsFactors = FALSE)
  # Basic check to avoid duplicates if re-running same species
  new_df <- rbind(old_df[old_df$species != species_name, ], summary_df)
} else {
  new_df <- summary_df
}
write.csv(new_df, file = global_file, row.names = FALSE)

message("All done. Outputs saved in: ", normalizePath(output_dir))