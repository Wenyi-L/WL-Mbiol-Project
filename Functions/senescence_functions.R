## =========================================================
## senescence_functions.R (Final Compatible Version)
##
## FIXES:
## 1. Reverted argument names to 'sx_senescence'/'fx_senescence'
##    to match your existing run scripts.
## 2. Logic remains WEIGHTED MEAN (Scientific Correctness).
## 3. Includes fixed exact_lifespan_pmf & helpers.
## =========================================================

## ---------------------------------------------------------
## Helper 1: Calculate lx-weighted mean for adult rates
## ---------------------------------------------------------
calculate_weighted_adult_means <- function(ages, sx, fx, maturity_age) {
  k <- length(ages)
  lx <- numeric(k); lx[1] <- 1
  if (k > 1) { for (i in 1:(k-1)) lx[i+1] <- lx[i] * sx[i] }
  
  idx_adult <- which(ages >= maturity_age)
  if (length(idx_adult) == 0) idx_adult <- k 
  
  w <- lx[idx_adult]
  if (sum(w) == 0) w <- rep(1, length(w))
  
  sx_mean <- weighted.mean(sx[idx_adult], w, na.rm = TRUE)
  fx_mean <- weighted.mean(fx[idx_adult], w, na.rm = TRUE)
  
  return(list(sx_mean = sx_mean, fx_mean = fx_mean))
}

## ---------------------------------------------------------
## Helper 2: Exact Lifespan PMF (Fixed Column Names & Robust)
## ---------------------------------------------------------
exact_lifespan_pmf <- function(U, mix, max_age = 500, ...) {
  # Added ... to swallow 'max_t' if passed
  pmf <- numeric(max_age)
  current_vec <- mix
  
  for(t in 1:max_age) {
    mortality_vec <- 1 - colSums(U)
    mortality_vec[mortality_vec < 0] <- 0
    pmf[t] <- sum(current_vec * mortality_vec)
    current_vec <- U %*% current_vec
  }
  
  if(sum(pmf) < 1) pmf[max_age] <- pmf[max_age] + (1 - sum(pmf))
  
  # Returns 'lifespan' and 'pmf' to satisfy senescence_plots.R
  return(data.frame(lifespan = 1:max_age, pmf = pmf))
}

## ---------------------------------------------------------
## Helper 3: Iterative LRO Dist
## ---------------------------------------------------------
calcDistLRO_iterative <- function(U, F, c0_vector, maxClutchSize = 20, maxLRO = 50, Fdist = "Poisson", ...) {
  ages <- 1:nrow(U)
  k_ages <- length(ages)
  sx <- numeric(k_ages)
  for(i in 1:(k_ages-1)) sx[i] <- U[i+1, i]
  sx[k_ages] <- U[k_ages, k_ages]
  
  B_list <- list()
  for (i in 1:k_ages) {
    mean_off <- sum(F[, i])
    surv_prob <- if (i < k_ages) sum(U[,i]) else U[k_ages, k_ages]
    lambda_val <- if (surv_prob > 0) mean_off / surv_prob else 0
    probs <- dpois(0:maxClutchSize, lambda = lambda_val)
    B_list[[i]] <- probs / sum(probs)
  }
  
  current_dist <- numeric(maxLRO + 1)
  if(!is.null(c0_vector) && length(c0_vector) == k_ages) current_dist[1] <- 1 else current_dist[1] <- 1
  final_dead_dist <- numeric(maxLRO + 1)
  
  for (age in 1:(k_ages*3)) { 
    idx <- if(age <= k_ages) age else k_ages
    s_rate <- sx[idx]
    clutch_probs <- B_list[[idx]]
    dead_frac <- 1 - s_rate
    final_dead_dist <- final_dead_dist + (current_dist * dead_frac)
    
    if (s_rate <= 0 || sum(current_dist) < 1e-6) break
    
    next_gen <- numeric(maxLRO + 1)
    idx_c <- which(current_dist > 0)
    idx_b <- which(clutch_probs > 0)
    
    for (i in idx_c) {
      for (j in idx_b) {
        ni <- i + j - 1
        if (ni <= (maxLRO + 1)) next_gen[ni] <- next_gen[ni] + current_dist[i] * clutch_probs[j]
        else next_gen[maxLRO + 1] <- next_gen[maxLRO + 1] + current_dist[i] * clutch_probs[j]
      }
    }
    current_dist <- next_gen * s_rate
  }
  
  if (sum(current_dist) > 0) final_dead_dist <- final_dead_dist + current_dist
  if(sum(final_dead_dist) > 0) final_dead_dist <- final_dead_dist / sum(final_dead_dist)
  return(final_dead_dist)
}

calcDistLROPostBreedingNoEnv <- function(U, F, c0, maxClutchSize = 20, maxLRO = 50, Fdist = "Poisson") {
  calcDistLRO_iterative(U, F, c0, maxClutchSize, maxLRO, Fdist)
}

## ---------------------------------------------------------
## Helper 4: Detect Maturity (Logistic)
## ---------------------------------------------------------
detect_maturity_age_logistic <- function(ages, fx, min_fx = 1e-8, use_smoothed_fx = TRUE, 
                                         span = 0.5, prob_threshold = 0.5, 
                                         nls_control = list(), return_fit = FALSE) {
  ages <- as.numeric(ages); fx <- as.numeric(fx)
  if (length(ages) != length(fx)) stop("ages/fx mismatch")
  if (all(is.na(fx)) || max(fx, na.rm=TRUE) <= 0) return(min(ages, na.rm=TRUE))
  
  fx_clean <- fx; fx_clean[!is.finite(fx_clean)] <- 0; fx_clean[fx_clean < 0] <- 0
  sm_vals <- fx_clean
  if (use_smoothed_fx && length(unique(ages[is.finite(fx_clean)])) >= 3) {
    try({
      lo <- stats::loess(fx_clean ~ ages, span = span, degree = 2, control = loess.control(surface = "interpolate"))
      sm_vals <- as.numeric(stats::predict(lo, newdata = ages))
    }, silent=TRUE)
  }
  sm_vals[!is.finite(sm_vals)] <- 0; sm_vals[sm_vals < 0] <- 0
  
  cum_raw <- cumsum(sm_vals)
  total <- tail(cum_raw, 1)
  if (!is.finite(total) || total <= 0) {
    idx <- which(fx_clean > min_fx)
    return(if(length(idx)>0) ages[min(idx)] else min(ages, na.rm=TRUE))
  }
  
  cum_prob <- cum_raw / total
  df_fit <- data.frame(age = ages, p = cum_prob, fx = sm_vals)
  
  start_x50 <- sum(df_fit$age * df_fit$fx)/sum(df_fit$fx)
  start_s <- (max(df_fit$age) - min(df_fit$age)) / 6
  if (!is.finite(start_s) || start_s <= 0) start_s <- 1
  
  fit_nls <- NULL
  try({ fit_nls <- stats::nls(p ~ 1 / (1 + exp(-(age - x50)/s)), data = df_fit, 
                              start = list(x50 = start_x50, s = start_s), control = do.call(nls.control, nls_control)) }, silent=TRUE)
  if (is.null(fit_nls) && requireNamespace("minpack.lm", quietly=TRUE)) {
    try({ fit_nls <- minpack.lm::nlsLM(p ~ 1 / (1 + exp(-(age - x50)/s)), data = df_fit, start = list(x50 = start_x50, s = start_s)) }, silent=TRUE)
  }
  
  ma_est <- NULL
  if (!is.null(fit_nls)) {
    coef_est <- coef(fit_nls)
    p_safe <- min(max(prob_threshold, 1e-6), 1 - 1e-6)
    ma_est <- coef_est["x50"] - coef_est["s"] * log(1/p_safe - 1)
  }
  
  if (is.null(ma_est)) {
    if (any(df_fit$p >= prob_threshold)) {
      i <- which(df_fit$p >= prob_threshold)[1]
      if (i == 1) ma_est <- df_fit$age[1]
      else {
        frac <- (prob_threshold - df_fit$p[i-1]) / (df_fit$p[i] - df_fit$p[i-1])
        ma_est <- df_fit$age[i-1] + frac * (df_fit$age[i] - df_fit$age[i-1])
      }
    } else ma_est <- max(ages)
  }
  ma_est <- max(ma_est, min(ages, na.rm=TRUE))
  if (return_fit) return(list(maturity_age = as.numeric(ma_est), fit = fit_nls))
  return(as.numeric(ma_est))
}

## ---------------------------------------------------------
## Data Prep
## ---------------------------------------------------------
prepare_demography_data_from_df <- function(dat, input_type="auto", maturity_age="auto", estimate_tail=TRUE, 
                                            tail_k=2, maturity_method="logistic", maturity_prob=0.5, 
                                            auto_min_fx=1e-8, auto_span=0.5, plot_maturity=FALSE, plot_path=NULL, ...) {
  if(maturity_method[1] == "logistic50") maturity_method <- "logistic"
  maturity_method <- match.arg(maturity_method, c("logistic", "absolute"))
  if (!("x" %in% names(dat))) stop("Dataset must contain column: x")
  dat <- dat[order(dat$x), ]
  ages <- dat$x; k <- nrow(dat)
  
  has_col <- function(n) n %in% names(dat) && any(!is.na(dat[[n]]))
  Nx <- if(has_col("Nx")) as.numeric(dat$Nx) else rep(NA, k)
  lx_raw <- if(has_col("lx")) as.numeric(dat$lx) else rep(NA, k)
  qx <- if(has_col("qx")) as.numeric(dat$qx) else rep(NA, k)
  noff <- if(has_col("noffspring")) as.numeric(dat$noffspring) else rep(NA, k)
  fert_mx <- if(has_col("fert.mx")) as.numeric(dat$fert.mx) else rep(NA, k)
  
  if(has_col("qx")) sx <- 1 - qx
  else if(has_col("lx") && k>=2) {
    sx <- rep(NA, k); idx <- which(!is.na(lx_raw[-k]) & lx_raw[-k]>0)
    sx[idx] <- lx_raw[idx+1]/lx_raw[idx]; sx[lx_raw==0] <- 0
  } else if(has_col("Nx") && k>=2) {
    sx <- rep(NA, k); idx <- which(!is.na(Nx[-k]) & Nx[-k]>0)
    sx[idx] <- Nx[idx+1]/Nx[idx]
  } else stop("No survival info")
  
  if(is.na(tail(sx,1)) && estimate_tail) {
    s_omega <- if(sum(!is.na(sx))>=tail_k) mean(tail(sx[!is.na(sx)], tail_k)) else 0
    sx[length(sx)] <- s_omega
  }
  sx[!is.na(sx)] <- pmin(pmax(sx[!is.na(sx)], 0), 0.9999); sx[!is.finite(sx)] <- 0
  
  if(has_col("fert.mx")) fx <- fert_mx
  else if(has_col("noffspring") && has_col("Nx")) fx <- noff/Nx
  else if(has_col("noffspring")) fx <- noff
  else stop("No repro info")
  fx[!is.finite(fx)] <- 0
  
  if(is.character(maturity_age) && maturity_age=="auto") {
    if(maturity_method=="logistic") {
      res <- detect_maturity_age_logistic(ages, fx, min_fx=auto_min_fx, span=auto_span, prob_threshold=maturity_prob)
      maturity_age_eff <- res
    } else {
      idx <- which(fx > auto_min_fx)
      maturity_age_eff <- if(length(idx)>0) ages[min(idx)] else min(ages)
    }
  } else maturity_age_eff <- maturity_age
  
  dat$sx <- sx; dat$fx <- fx
  list(data=dat, ages=ages, sx=sx, fx=fx, maturity_age=maturity_age_eff)
}

## ---------------------------------------------------------
## Model Builders (RESTORED ARGUMENT NAMES)
## ---------------------------------------------------------
build_MPM_senescence <- function(ages, sx, fx) {
  k <- length(ages)
  U <- matrix(0, k, k); if(k>1) for(i in 1:(k-1)) U[i+1,i] <- sx[i]; U[k,k] <- sx[k]
  F <- matrix(0, k, k); F[1,] <- fx
  list(ages=ages, sx=sx, fx=fx, U=U, F=F, A=U+F)
}

# NOTE: Reverted to sx_senescence / fx_senescence to fix your error
build_MPM_no_senescence <- function(ages, sx_senescence, fx_senescence, maturity_age) {
  # Internally uses the weighted mean logic, but accepts old arg names
  w <- calculate_weighted_adult_means(ages, sx_senescence, fx_senescence, maturity_age)
  idx <- which(ages >= maturity_age)
  
  sx_no <- sx_senescence
  sx_no[idx] <- w$sx_mean
  
  fx_no <- fx_senescence
  fx_no[idx] <- w$fx_mean
  
  build_MPM_senescence(ages, sx_no, fx_no)
}

build_MPM_no_actuarial_yes_reproductive <- function(ages, sx_senescence, fx_senescence, maturity_age) {
  w <- calculate_weighted_adult_means(ages, sx_senescence, fx_senescence, maturity_age)
  idx <- which(ages >= maturity_age)
  
  sx_no <- sx_senescence
  sx_no[idx] <- w$sx_mean
  
  # Fecundity stays as is (senescence)
  build_MPM_senescence(ages, sx_no, fx_senescence)
}

build_MPM_yes_actuarial_no_reproductive <- function(ages, sx_senescence, fx_senescence, maturity_age) {
  w <- calculate_weighted_adult_means(ages, sx_senescence, fx_senescence, maturity_age)
  idx <- which(ages >= maturity_age)
  
  fx_no <- fx_senescence
  fx_no[idx] <- w$fx_mean
  
  # Survival stays as is (senescence)
  build_MPM_senescence(ages, sx_senescence, fx_no)
}

compute_summary_table <- function(U_sen, U_no, U_noA_yesR, U_yesA_noR, F_sen, F_no, F_noA_yesR, F_yesA_noR, mix_sen, mix_no, mix_noA_yesR, mix_yesA_noR, repro_var="Poisson") {
  if(!exists("mean_lifespan")) stop("Source LuckFunctions.R first!")
  
  .get_moments <- function(U, F, mix) {
    mL <- mean_lifespan(U, mix) 
    vL <- var_lifespan(U, mix)
    sL <- skew_lifespan(U, mix)
    mR <- mean_LRO(U, F, mix)
    vR <- var_LRO_mcr(U, F, repro_var, mix)
    sR <- skew_LRO(U, F, repro_var, mix)
    c(mL=as.numeric(mL), vL=as.numeric(vL), sL=as.numeric(sL), mR=as.numeric(mR), vR=as.numeric(vR), sR=as.numeric(sR))
  }
  
  r1 <- .get_moments(U_sen, F_sen, mix_sen)
  r2 <- .get_moments(U_no, F_no, mix_no)
  r3 <- .get_moments(U_noA_yesR, F_noA_yesR, mix_noA_yesR)
  r4 <- .get_moments(U_yesA_noR, F_yesA_noR, mix_yesA_noR)
  
  data.frame(
    model = c("Senescence", "No-senescence", "No-actuarial/Yes-reproductive", "Yes-actuarial/No-reproductive"),
    mean_lifespan = c(r1["mL"], r2["mL"], r3["mL"], r4["mL"]),
    var_lifespan = c(r1["vL"], r2["vL"], r3["vL"], r4["vL"]),
    skew_lifespan = c(r1["sL"], r2["sL"], r3["sL"], r4["sL"]),
    mean_LRO = c(r1["mR"], r2["mR"], r3["mR"], r4["mR"]),
    var_LRO = c(r1["vR"], r2["vR"], r3["vR"], r4["vR"]),
    skew_LRO = c(r1["sR"], r2["sR"], r3["sR"], r4["sR"])
  )
}