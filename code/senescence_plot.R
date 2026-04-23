## =========================================================
## senescence_plots.R (Updated 03/03/2026)
## Content: Pure plotting functions for survival, fecundity, LS, LRO
## Note: Requires library(ggplot2) to be loaded globally.
## =========================================================


## ---------------------------------------------------------
## Survival & fecundity plotting
## ---------------------------------------------------------
plot_survivorship_curve <- function(dat, species_name) {
  if (!("sx" %in% names(dat))) stop("Column 'sx' not found in dat. Make sure prepare_demography_data() was run.")
  if (!("x" %in% names(dat))) stop("Column 'x' (age) not found in dat.")
  df_plot <- data.frame(x = as.numeric(dat$x), sx = as.numeric(dat$sx))
  df_plot <- df_plot[is.finite(df_plot$sx), , drop = FALSE]
  
  p <- ggplot(df_plot, aes(x = x, y = sx)) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 2) +
    theme_classic(base_size = 14) +
    labs(title = paste0(species_name, " Survival Curve (senescence)"), x = "Age (x)", y = "Age-specific survival probability (sx)") +
    coord_cartesian(ylim = c(0, 1))
  p
}

plot_survival_models <- function(ages, sx_sen, sx_no, species_name) {
  df_plot <- data.frame(x = c(ages, ages), sx = c(as.numeric(sx_sen), as.numeric(sx_no)), model = rep(c("Senescence", "No-senescence"), each = length(ages)))
  df_plot <- df_plot[is.finite(df_plot$sx) & is.finite(df_plot$x), , drop = FALSE]
  df_plot$model <- factor(df_plot$model, levels = c("Senescence", "No-senescence"))
  
  p <- ggplot(df_plot, aes(x = x, y = sx, color = model, linetype = model)) +
    geom_line(linewidth = 1.1, alpha = 0.9) +
    theme_classic(base_size = 14) +
    labs(title = paste0(species_name, " Survival (senescence vs no-senescence)"), x = "Age (x)", y = "Age-specific survival probability (sx)", color = NULL, linetype = NULL) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_color_manual(values = c("Senescence" = "#1F78B4", "No-senescence" = "#E31A1C")) +
    scale_linetype_manual(values = c("Senescence" = "solid", "No-senescence" = "solid")) +
    theme(legend.position = "top")
  p
}

plot_fecundity_models <- function(ages, fx_sen, fx_no, species_name) {
  df_plot <- data.frame(x = c(ages, ages), fx = c(as.numeric(fx_sen), as.numeric(fx_no)), model = rep(c("Senescence", "No-senescence"), each = length(ages)))
  df_plot <- df_plot[is.finite(df_plot$fx) & is.finite(df_plot$x), , drop = FALSE]
  df_plot$model <- factor(df_plot$model, levels = c("Senescence", "No-senescence"))
  
  p <- ggplot(df_plot, aes(x = x, y = fx, color = model, linetype = model)) +
    geom_line(linewidth = 1.1, alpha = 0.9) +
    theme_classic(base_size = 14) +
    labs(title = paste0(species_name, " Fecundity (senescence vs no-senescence)"), x = "Age (x)", y = "Age-specific fecundity (fx)", color = NULL, linetype = NULL) +
    theme(legend.position = "top") +
    scale_color_manual(values = c("Senescence" = "#1F78B4", "No-senescence" = "#E31A1C")) +
    scale_linetype_manual(values = c("Senescence" = "solid", "No-senescence" = "solid"))
  p
}



## ---------------------------------------------------------
## Lifespan distributions (No LONGER NEED mix_sen,mix_no)
## ---------------------------------------------------------
plot_lifespan_distributions <- function(U_sen, U_no) {
  LS_exact_sen <- calcDistLifespan(U_sen)
  LS_exact_no  <- calcDistLifespan(U_no)
  LS_exact_sen <- data.frame(lifespan = 1:length(LS_exact_sen), pmf = as.numeric(LS_exact_sen))
  LS_exact_no  <- data.frame(lifespan = 1:length(LS_exact_no),  pmf = as.numeric(LS_exact_no)) #convert numeric values into list
  ages_all <- sort(unique(c(LS_exact_sen$lifespan, LS_exact_no$lifespan)))
  pmf_sen_all <- rep(0, length(ages_all))
  pmf_no_all  <- rep(0, length(ages_all))
  idx_sen <- match(LS_exact_sen$lifespan, ages_all)
  idx_no  <- match(LS_exact_no$lifespan,  ages_all)
  pmf_sen_all[idx_sen] <- LS_exact_sen$pmf
  pmf_no_all[idx_no]   <- LS_exact_no$pmf
  pmf_global <- pmf_sen_all + pmf_no_all
  pmf_global <- pmf_global / sum(pmf_global)
  LS_exact_sen$model <- "Senescence"
  LS_exact_no$model  <- "No-senescence"
  LS_all <- rbind(LS_exact_sen, LS_exact_no)
  xmax_exact <- max(LS_all$lifespan, na.rm = TRUE) #limit the range of x-axis
  LS_all_vis <- subset(LS_all, lifespan <= xmax_exact & is.finite(pmf))
  
  # Ensure factor levels for consistency
  LS_all_vis$model <- factor(LS_all_vis$model, levels = c("Senescence", "No-senescence"))
  
  
  pLS <- ggplot(LS_all_vis, aes(x = lifespan, y = pmf, color = model, fill = model)) +
    geom_line(linewidth = 1.2, alpha = 0.6) +
    geom_area(alpha = 0.10, position = "identity") +
    scale_color_manual(values = c("Senescence" = "#1F78B4", "No-senescence" = "#E31A1C")) +
    scale_fill_manual(values = c("Senescence" = "#1F78B4", "No-senescence" = "#E31A1C")) +
    labs(x = "Lifespan (steps, min = 1)", y = "Probability mass (PMF)", title = "Lifespan distribution (Analytical)") +
    coord_cartesian(xlim = c(1, xmax_exact)) +
    theme_classic(base_size = 14) +
    theme(legend.title = element_blank(), legend.position = "top")
  pLS
}

## ---------------------------------------------------------
## LRO distributions (ANALYTICAL, pre-breeding) NO longer need mix_sen, mix_no
## ---------------------------------------------------------
# Single-model (2-model wrapper kept for backward compat)
plot_LRO_distributions <- function(U_sen, U_no,
                                   F_sen, F_no,
                                   maxClutchSize = 20,
                                   maxLRO = 50,
                                   include_zero = TRUE
                                   ) {
  # U_* and F_* are age-structured pre-breeding U and F (as in build_MPM_*).
  # We'll use calcDistLRONoEnv since models are no-env here.
  
  k_sen <- nrow(U_sen)
  k_no  <- nrow(U_no)
  
  # c0: newborn state distribution -- assume newborns start in first age/class
  c0_sen <- rep(0, k_sen); c0_sen[1] <- 1
  c0_no  <- rep(0, k_no);  c0_no[1] <- 1
  
  dist_sen <- tryCatch({
    calcDistLRONoEnv(U_sen, F_sen, c0_sen, maxClutchSize, maxLRO)
  }, error = function(e) {
    stop("calcDistLRONoEnv failed for Senescence model: ", conditionMessage(e))
  })
  dist_no <- tryCatch({
    calcDistLRONoEnv(U_no, F_no, c0_no, maxClutchSize, maxLRO)
  }, error = function(e) {
    stop("calcDistLRONoEnv failed for No-senescence model: ", conditionMessage(e))
  })
  
  # dist_* are numeric vectors indexed from LRO = 0..maxLRO
  # Ensure they are normalized
  dist_sen <- as.numeric(dist_sen); dist_no <- as.numeric(dist_no)
  if (sum(dist_sen) > 0) dist_sen <- dist_sen / sum(dist_sen)
  if (sum(dist_no) > 0)  dist_no  <- dist_no  / sum(dist_no)
  
  df_sen <- data.frame(LRO = 0:(length(dist_sen)-1), pmf = dist_sen, model = "Senescence")
  df_no  <- data.frame(LRO = 0:(length(dist_no)-1),  pmf = dist_no,  model = "No-senescence")
  df_all <- rbind(df_sen, df_no)
  
  xmax_all <- max(df_all$LRO, na.rm = TRUE) #limit the range of x-axis
  
  # optionally drop LRO == 0 from plot if requested
  if (!include_zero) df_all <- subset(df_all, LRO > 0)
  
  
  pLRO <- ggplot(df_all, aes(x = LRO, y = pmf, color = model, fill = model)) +
    geom_line(linewidth = 1.2, alpha = 0.7) +
    geom_area(alpha = 0.10, position = "identity") +
    scale_color_manual(values = c("Senescence" = "#1F78B4", "No-senescence" = "#E31A1C")) +
    scale_fill_manual(values = c("Senescence" = "#1F78B4", "No-senescence" = "#E31A1C")) +
    labs(x = "LRO (T)", y = "Probability mass (PMF)", title = "LRO distribution (Analytical, pre-breeding)") +
    coord_cartesian(xlim = c(0, xmax_all)) +
    theme_classic(base_size = 14) +
    theme(legend.title = element_blank(), legend.position = "top")
  
  pLRO
}


# Four-model version: uses analytical distributions
plot_LRO_distributions_4 <- function(
    U_sen, U_no, U_noA_yesR, U_yesA_noR,
    F_sen, F_no, F_noA_yesR, F_yesA_noR,
    maxClutchSize = 20,
    maxLRO = 50,
    include_zero = TRUE
) {
  
  # c0 for each model (assume newborns start in class 1)
  c0_sen <- rep(0, nrow(U_sen)); c0_sen[1] <- 1
  c0_no  <- rep(0, nrow(U_no));  c0_no[1]  <- 1
  c0_noA <- rep(0, nrow(U_noA_yesR)); c0_noA[1] <- 1
  c0_yesA <- rep(0, nrow(U_yesA_noR)); c0_yesA[1] <- 1
  
  dist_sen      <- calcDistLRONoEnv(U_sen,      F_sen,      c0_sen,    maxClutchSize, maxLRO)
  dist_no       <- calcDistLRONoEnv(U_no,       F_no,       c0_no,     maxClutchSize, maxLRO)
  dist_noA_yesR <- calcDistLRONoEnv(U_noA_yesR, F_noA_yesR, c0_noA,    maxClutchSize, maxLRO)
  dist_yesA_noR <- calcDistLRONoEnv(U_yesA_noR, F_yesA_noR, c0_yesA,   maxClutchSize, maxLRO)
  
  dist_sen      <- as.numeric(dist_sen); dist_no       <- as.numeric(dist_no)
  dist_noA_yesR <- as.numeric(dist_noA_yesR); dist_yesA_noR <- as.numeric(dist_yesA_noR)
  if (sum(dist_sen) > 0) dist_sen <- dist_sen / sum(dist_sen)
  if (sum(dist_no) > 0)  dist_no  <- dist_no  / sum(dist_no)
  if (sum(dist_noA_yesR) > 0) dist_noA_yesR <- dist_noA_yesR / sum(dist_noA_yesR)
  if (sum(dist_yesA_noR) > 0) dist_yesA_noR <- dist_yesA_noR / sum(dist_yesA_noR)
  
  PMF_sen      <- data.frame(LRO = 0:(length(dist_sen)-1), pmf = dist_sen, model = "Senescence")
  PMF_no       <- data.frame(LRO = 0:(length(dist_no)-1),  pmf = dist_no,  model = "No-senescence")
  PMF_noA_yesR <- data.frame(LRO = 0:(length(dist_noA_yesR)-1), pmf = dist_noA_yesR, model = "No-actuarial/Yes-reproductive")
  PMF_yesA_noR <- data.frame(LRO = 0:(length(dist_yesA_noR)-1), pmf = dist_yesA_noR, model = "Yes-actuarial/No-reproductive")
  
  PMF_all <- rbind(PMF_sen, PMF_no, PMF_noA_yesR, PMF_yesA_noR)
  if (!include_zero) PMF_all <- subset(PMF_all, LRO > 0)
  
  xmax_all <- max(PMF_all$LRO, na.rm = TRUE)#limit the range of x-axis
  
  PMF_all_vis <- subset(PMF_all, LRO <= xmax_all & is.finite(pmf))
  PMF_all_vis$model <- factor(PMF_all_vis$model, levels = c("Senescence", "No-senescence", "No-actuarial/Yes-reproductive", "Yes-actuarial/No-reproductive"))
  
  pLRO4 <- ggplot(PMF_all_vis, aes(x = LRO, y = pmf, color = model, fill = model)) +
    geom_line(linewidth = 1.2, alpha = 0.7) +
    geom_area(alpha = 0.10, position = "identity") +
    scale_color_manual(values = c("Senescence" = "#1F78B4", "No-senescence" = "#E31A1C", "No-actuarial/Yes-reproductive" = "#33A02C", "Yes-actuarial/No-reproductive" = "#6A3D9A")) +
    scale_fill_manual(values = c("Senescence" = "#1F78B4", "No-senescence" = "#E31A1C", "No-actuarial/Yes-reproductive" = "#33A02C", "Yes-actuarial/No-reproductive" = "#6A3D9A")) +
    labs(x = "LRO (T)", y = "Probability mass (PMF)", title = "LRO distribution (4 models, analytical, pre-breeding)") +
    coord_cartesian(xlim = c(0, xmax_all)) +
    theme_classic(base_size = 14) +
    theme(legend.title = element_blank(), legend.position = "top")
  
  pLRO4
}