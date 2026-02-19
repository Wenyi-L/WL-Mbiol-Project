## =========================================================
## senescence_plots.R
## Content: Pure plotting functions for survival, fecundity, LS, LRO
## Note: Requires library(ggplot2) to be loaded globally.
## =========================================================

## ---------------------------------------------------------
## plot_maturity_detection
## ---------------------------------------------------------
plot_maturity_detection <- function(ages = NULL, fx = NULL, dat = NULL,
                                    use_smoothed_fx = TRUE, span = 0.5,
                                    nls_control = list(), plot_file = NULL,
                                    title = NULL) {
  if (is.null(dat)) {
    if (is.null(ages) || is.null(fx)) stop("Provide either dat (with x and fx) or both ages and fx.")
    dat_in <- data.frame(x = ages, fx = fx)
  } else {
    if (!("x" %in% names(dat)) || !("fx" %in% names(dat))) stop("dat must contain columns 'x' and 'fx'.")
    dat_in <- data.frame(x = as.numeric(dat$x), fx = as.numeric(dat$fx))
  }
  
  # compute fit + get maturity age with fit data
  res <- detect_maturity_age_logistic50(dat_in$x, dat_in$fx,
                                        use_smoothed_fx = use_smoothed_fx,
                                        span = span,
                                        nls_control = nls_control,
                                        return_fit = TRUE)
  ma <- res$maturity_age
  fit_data <- res$fit_data
  fit_obj <- res$fit
  
  # prepare plotting data: points (age vs cumulative prob)
  df_plot <- data.frame(age = fit_data$age, cum_prob = fit_data$p, fx = fit_data$fx)
  
  # build smooth curve for plotting:
  if (!is.null(fit_obj) && inherits(fit_obj, "nls")) {
    age_seq <- seq(min(df_plot$age), max(df_plot$age), length.out = 300)
    coef_est <- coef(fit_obj)
    x50 <- as.numeric(coef_est["x50"]); s <- as.numeric(coef_est["s"])
    logistic_fun <- function(x) 1 / (1 + exp(-(x - x50)/s))
    df_curve <- data.frame(age = age_seq, cum_prob = logistic_fun(age_seq))
  } else {
    age_seq <- seq(min(df_plot$age), max(df_plot$age), length.out = 300)
    df_curve <- data.frame(age = age_seq, cum_prob = stats::approx(df_plot$age, df_plot$cum_prob, xout = age_seq)$y)
  }
  
  # Ensure ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it (install.packages('ggplot2')).")
  }
  
  p <- ggplot() +
    geom_line(data = df_curve, aes(x = age, y = cum_prob), linewidth = 1) +
    geom_point(data = df_plot, aes(x = age, y = cum_prob), size = 2.5) +
    geom_vline(xintercept = ma, linetype = "dashed") +
    annotate("text", x = ma, y = 0.05, label = sprintf("50%% at age = %.3f", as.numeric(ma)),
             hjust = 0, vjust = 0, size = 3.5) +
    labs(x = "Age", y = "Cumulative reproduction probability",
         title = if (is.null(title)) "Maturity detection (logistic50)" else title) +
    coord_cartesian(ylim = c(0, 1)) + 
    theme_bw(base_size = 13)
  
  
  if (!is.null(plot_file)) {
    ggsave(plot_file, plot = p, width = 7, height = 5)
  }
  p
}

## ---------------------------------------------------------
## Survival & fecundity plotting
## (Keeps original colors but adds linetypes for visibility)
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
    scale_linetype_manual(values = c("Senescence" = "solid", "No-senescence" = "dashed")) +
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
    scale_linetype_manual(values = c("Senescence" = "solid", "No-senescence" = "dashed"))
  p
}

## ---------------------------------------------------------
## Helper for adaptive global cutoff
## ---------------------------------------------------------
.choose_global_cutoff_age <- function(ages_all,
                                      pmf_global,
                                      q_max  = 0.999,
                                      q_min  = 0.95,
                                      step   = 0.001,
                                      factor = 2) {
  pmf_global <- as.numeric(pmf_global)
  pmf_global[!is.finite(pmf_global)] <- 0
  pmf_global <- pmf_global / sum(pmf_global)
  cdf <- cumsum(pmf_global)
  .get_age_at_q <- function(q) {
    idx <- which(cdf >= q)[1]
    if (is.na(idx)) {
      return(max(ages_all, na.rm = TRUE))
    } else {
      return(ages_all[idx])
    }
  }
  age_90 <- .get_age_at_q(0.90)
  if (!is.finite(age_90)) return(max(ages_all, na.rm = TRUE))
  qs <- seq(q_max, q_min, by = -step)
  chosen_age <- NA_real_
  for (q in qs) {
    age_q <- .get_age_at_q(q)
    if (is.finite(age_q) && age_q <= age_90 * factor) {
      chosen_age <- age_q
      break
    }
  }
  if (!is.finite(chosen_age)) chosen_age <- .get_age_at_q(q_min)
  ages_finite <- sort(unique(ages_all[is.finite(ages_all)]))
  if (length(ages_finite) >= 2) {
    min_cutoff <- ages_finite[2]
  } else {
    min_cutoff <- ages_finite[1]
  }
  if (!is.finite(chosen_age)) chosen_age <- max(ages_finite)
  if (chosen_age < min_cutoff) chosen_age <- min_cutoff
  chosen_age
}

## ---------------------------------------------------------
## Lifespan distributions
## ---------------------------------------------------------
plot_lifespan_distributions <- function(U_sen, U_no, mix_sen, mix_no, qcut = 0.999) {
  LS_exact_sen <- exact_lifespan_pmf(U_sen, mix_sen)
  LS_exact_no  <- exact_lifespan_pmf(U_no,  mix_no)
  ages_all <- sort(unique(c(LS_exact_sen$lifespan, LS_exact_no$lifespan)))
  pmf_sen_all <- rep(0, length(ages_all))
  pmf_no_all  <- rep(0, length(ages_all))
  idx_sen <- match(LS_exact_sen$lifespan, ages_all)
  idx_no  <- match(LS_exact_no$lifespan,  ages_all)
  pmf_sen_all[idx_sen] <- LS_exact_sen$pmf
  pmf_no_all[idx_no]   <- LS_exact_no$pmf
  pmf_global <- pmf_sen_all + pmf_no_all
  pmf_global <- pmf_global / sum(pmf_global)
  xmax_exact <- .choose_global_cutoff_age(ages_all = ages_all, pmf_global = pmf_global, q_max = qcut, q_min = 0.95, step = 0.001, factor = 2)
  LS_exact_sen$model <- "Senescence"
  LS_exact_no$model  <- "No-senescence"
  LS_all <- rbind(LS_exact_sen, LS_exact_no)
  LS_all_vis <- subset(LS_all, lifespan <= xmax_exact & is.finite(pmf))
  
  # Ensure factor levels for consistency
  LS_all_vis$model <- factor(LS_all_vis$model, levels = c("Senescence", "No-senescence"))
  
  pLS <- ggplot(LS_all_vis, aes(x = lifespan, y = pmf, color = model, fill = model)) +
    geom_line(linewidth = 1.2, alpha = 0.6) +
    geom_area(alpha = 0.10, position = "identity") +
    scale_color_manual(values = c("Senescence" = "#1F78B4", "No-senescence" = "#E31A1C")) +
    scale_fill_manual(values = c("Senescence" = "#1F78B4", "No-senescence" = "#E31A1C")) +
    labs(x = "Lifespan (steps, min = 1)", y = "Probability mass (PMF)", title = "Lifespan distribution (Exact; tail adaptively truncated by global CDF)") +
    coord_cartesian(xlim = c(1, xmax_exact)) +
    theme_classic(base_size = 14) +
    theme(legend.title = element_blank(), legend.position = "top")
  pLS
}

## ---------------------------------------------------------
## LRO distributions (ANALYTICAL, post-breeding)
## ---------------------------------------------------------
# Single-model (2-model wrapper kept for backward compat)
plot_LRO_distributions <- function(U_sen, U_no,
                                   F_sen, F_no,
                                   mix_sen, mix_no,
                                   maxClutchSize = 20,
                                   maxLRO = 50,
                                   include_zero = FALSE,
                                   qcut = 0.999) {
  # U_* and F_* are age-structured post-breeding U and F (as in build_MPM_*).
  # mix_* are initial state distributions (vectors length = nrow(U)).
  # We'll use calcDistLROPostBreedingNoEnv since models are no-env here.
  
  k_sen <- nrow(U_sen)
  k_no  <- nrow(U_no)
  
  # c0: newborn state distribution -- assume newborns start in first age/class
  c0_sen <- rep(0, k_sen); c0_sen[1] <- 1
  c0_no  <- rep(0, k_no);  c0_no[1] <- 1
  
  dist_sen <- tryCatch({
    calcDistLROPostBreedingNoEnv(U_sen, F_sen, c0_sen, maxClutchSize, maxLRO)
  }, error = function(e) {
    stop("calcDistLROPostBreedingNoEnv failed for Senescence model: ", conditionMessage(e))
  })
  dist_no <- tryCatch({
    calcDistLROPostBreedingNoEnv(U_no, F_no, c0_no, maxClutchSize, maxLRO)
  }, error = function(e) {
    stop("calcDistLROPostBreedingNoEnv failed for No-senescence model: ", conditionMessage(e))
  })
  
  # dist_* are numeric vectors indexed from LRO = 0..maxLRO
  # Ensure they are normalized
  dist_sen <- as.numeric(dist_sen); dist_no <- as.numeric(dist_no)
  if (sum(dist_sen) > 0) dist_sen <- dist_sen / sum(dist_sen)
  if (sum(dist_no) > 0)  dist_no  <- dist_no  / sum(dist_no)
  
  df_sen <- data.frame(LRO = 0:(length(dist_sen)-1), pmf = dist_sen, model = "Senescence")
  df_no  <- data.frame(LRO = 0:(length(dist_no)-1),  pmf = dist_no,  model = "No-senescence")
  df_all <- rbind(df_sen, df_no)
  
  # optionally drop LRO == 0 from plot if requested
  if (!include_zero) df_all <- subset(df_all, LRO > 0)
  
  # adaptive global cutoff in LRO-space:
  T_all <- sort(unique(df_all$LRO))
  pmf_global <- tapply(df_all$pmf, df_all$LRO, sum)
  pmf_global <- as.numeric(pmf_global); names(pmf_global) <- T_all
  xmax_all <- .choose_global_cutoff_age(ages_all = as.numeric(names(pmf_global)), pmf_global = pmf_global, q_max = qcut, q_min = 0.95, step = 0.001, factor = 2)
  
  pLRO <- ggplot(df_all, aes(x = LRO, y = pmf, color = model, fill = model)) +
    geom_line(linewidth = 1.2, alpha = 0.7) +
    geom_area(alpha = 0.10, position = "identity") +
    scale_color_manual(values = c("Senescence" = "#1F78B4", "No-senescence" = "#E31A1C")) +
    scale_fill_manual(values = c("Senescence" = "#1F78B4", "No-senescence" = "#E31A1C")) +
    labs(x = "LRO (T)", y = "Probability mass (PMF)", title = "LRO distribution (Analytical, post-breeding)") +
    coord_cartesian(xlim = c(0, xmax_all)) +
    theme_classic(base_size = 14) +
    theme(legend.title = element_blank(), legend.position = "top")
  
  pLRO
}


# Four-model version: uses analytical distributions
plot_LRO_distributions_4 <- function(
    U_sen, U_no, U_noA_yesR, U_yesA_noR,
    F_sen, F_no, F_noA_yesR, F_yesA_noR,
    mix_sen, mix_no, mix_noA_yesR, mix_yesA_noR,
    maxClutchSize = 20,
    maxLRO = 50,
    include_zero = FALSE,
    qcut = 0.999
) {
  
  # c0 for each model (assume newborns start in class 1)
  c0_sen <- rep(0, nrow(U_sen)); c0_sen[1] <- 1
  c0_no  <- rep(0, nrow(U_no));  c0_no[1]  <- 1
  c0_noA <- rep(0, nrow(U_noA_yesR)); c0_noA[1] <- 1
  c0_yesA <- rep(0, nrow(U_yesA_noR)); c0_yesA[1] <- 1
  
  dist_sen      <- calcDistLROPostBreedingNoEnv(U_sen,      F_sen,      c0_sen,    maxClutchSize, maxLRO)
  dist_no       <- calcDistLROPostBreedingNoEnv(U_no,       F_no,       c0_no,     maxClutchSize, maxLRO)
  dist_noA_yesR <- calcDistLROPostBreedingNoEnv(U_noA_yesR, F_noA_yesR, c0_noA,    maxClutchSize, maxLRO)
  dist_yesA_noR <- calcDistLROPostBreedingNoEnv(U_yesA_noR, F_yesA_noR, c0_yesA,   maxClutchSize, maxLRO)
  
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
  
  T_all <- sort(unique(PMF_all$LRO))
  pmf_global <- tapply(PMF_all$pmf, PMF_all$LRO, sum)
  pmf_global <- as.numeric(pmf_global); names(pmf_global) <- T_all
  xmax_all <- .choose_global_cutoff_age(ages_all = as.numeric(names(pmf_global)), pmf_global = pmf_global, q_max = qcut, q_min = 0.95, step = 0.001, factor = 2)
  
  PMF_all_vis <- subset(PMF_all, LRO <= xmax_all & is.finite(pmf))
  PMF_all_vis$model <- factor(PMF_all_vis$model, levels = c("Senescence", "No-senescence", "No-actuarial/Yes-reproductive", "Yes-actuarial/No-reproductive"))
  
  pLRO4 <- ggplot(PMF_all_vis, aes(x = LRO, y = pmf, color = model, fill = model)) +
    geom_line(linewidth = 1.2, alpha = 0.7) +
    geom_area(alpha = 0.10, position = "identity") +
    scale_color_manual(values = c("Senescence" = "#1F78B4", "No-senescence" = "#E31A1C", "No-actuarial/Yes-reproductive" = "#33A02C", "Yes-actuarial/No-reproductive" = "#6A3D9A")) +
    scale_fill_manual(values = c("Senescence" = "#1F78B4", "No-senescence" = "#E31A1C", "No-actuarial/Yes-reproductive" = "#33A02C", "Yes-actuarial/No-reproductive" = "#6A3D9A")) +
    labs(x = "LRO (T)", y = "Probability mass (PMF)", title = "LRO distribution (4 models, analytical, post-breeding)") +
    coord_cartesian(xlim = c(0, xmax_all)) +
    theme_classic(base_size = 14) +
    theme(legend.title = element_blank(), legend.position = "top")
  
  pLRO4
}