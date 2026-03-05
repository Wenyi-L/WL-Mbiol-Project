## =========================================================
## senescence_functions.R (Updated 03/03/2026)
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
calculate_weighted_adult_means <- function(ages, sx, fx, senescence_onset_age) {
  k <- length(ages)
  lx <- numeric(k); lx[1] <- 1
  if (k > 1) { for (i in 1:(k-1)) lx[i+1] <- lx[i] * sx[i] }
  
  idx_adult <- which(ages >= senescence_onset_age)
  if (length(idx_adult) == 0) idx_adult <- k 
  
  w <- lx[idx_adult]
  if (sum(w) == 0) w <- rep(1, length(w))
  
  sx_mean <- weighted.mean(sx[idx_adult], w, na.rm = TRUE)
  fx_mean <- weighted.mean(fx[idx_adult], w, na.rm = TRUE)
  
  return(list(sx_mean = sx_mean, fx_mean = fx_mean))
}

## ---------------------------------------------------------
## Helper 2: Exact Lifespan (Updated according to calcDistLRO)
## ---------------------------------------------------------
calcDistLifespan <- function(U, c0_vector, Fdist = "Poisson", ...) {
  ages <- 1:nrow(U)
  k_ages <- length(ages)
  sx <- numeric(k_ages)
  for(i in 1:(k_ages-1)) sx[i] <- U[i+1, i]
  sx[k_ages] <- U[k_ages, k_ages]
 
  p_alive <- 1
  dist_lifespan <- numeric()
  age <- 1
  
  while (p_alive > 0.0001) {
    idx <- if(age <= k_ages) age else k_ages
    s_rate <- sx[idx]
    dead_prop <- p_alive*(1-sx[idx])
    dist_lifespan[age] <- dead_prop
    p_alive <- p_alive * sx[idx]
    age = age + 1
    
  }
  
  if (sum(dist_lifespan)>0) final_dist_lifespan <- dist_lifespan
  if (sum(final_dist_lifespan)>0) final_dist_lifespan <-final_dist_lifespan/sum(final_dist_lifespan)
  return(final_dist_lifespan)
  
}

## ---------------------------------------------------------
## Helper 3: USE calcDistLRONoEnv FROM Robin's code for LRO distribution calculation
## ---------------------------------------------------------

## ---------------------------------------------------------
## Data Prep
## ---------------------------------------------------------
prepare_demography_data_from_df <- function(dat, input_type="auto", study_type = "Unknown",...) {
  if (!("x" %in% names(dat))) stop("Dataset must contain column: x")
  dat <- dat[order(dat$x), ]
  ages <- dat$x; k <- nrow(dat)
  
  has_col <- function(n) n %in% names(dat) && any(!is.na(dat[[n]]))
  has_name <- function(n) n %in% names(dat)
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
  } else stop("No survival info") #This step is different methods to calculate sx from different variable stypes in the original dataset
  
  if (is.na(tail(sx, 1))) {
    if (grepl("IBCohort", study_type, ignore.case = TRUE)) {
      sx[!is.finite(sx)] <- 0
      message("     [Info] Cohort data detected: setting tail NA to 0")
    } else if (grepl("LTPeriod|Modelled", study_type, ignore.case = TRUE)) {
      valid_sx <- sx[!is.na(sx)]
      if (length(valid_sx) > 0) {
        sx[is.na(sx)] <- tail(valid_sx, 1) 
        message("     [Info] Period/Modelled data detected: carrying forward last survival rate")
      } else {
        sx[is.na(sx)] <- 0 
      }
    } else {
      sx[!is.finite(sx)] <- 0 
    }
  }
  # convert else NAs to 0
  sx[!is.finite(sx)] <-0
  
  if(has_col("fert.mx") && has_name("fx")) fx <- fert_mx
  else if(has_col("fert.mx") && has_name("mx")) fx <-fert_mx*sx[1] #converts mx to fx
  else if(has_col("noffspring") && has_col("Nx") && has_name("fx")) fx <- noff/Nx
  else if(has_col("noffspring") && has_col("Nx") && has_name("mx")) fx <- sx[1]*noff/Nx
  else stop("No repro info")
  fx[!is.finite(fx)] <- 0
  
  ##Senescence onset age##
  senescence_onset_age <- ages[which.max(fx)]
  half_peak <- 0.5*max(fx, na.rm = TRUE)
  above_half_peak <- which(fx >= half_peak)
  early_onset <- ages[min(above_half_peak)]
  late_onset <- ages[max(above_half_peak)]
  
  dat$sx <- sx; dat$fx <- fx
  list(data=dat, 
       ages=ages, 
       sx=sx, 
       fx=fx,
       senescence_onset_age = senescence_onset_age,
       early_onset = early_onset,
       late_onset = late_onset)
}




## ---------------------------------------------------------
## Model Builders (Removed MIXDIST) 
## ---------------------------------------------------------
build_MPM_senescence <- function(ages, sx, fx) {
  k <- length(ages)
  U <- matrix(0, k, k); if(k>1) for(i in 1:(k-1)) U[i+1,i] <- sx[i]; U[k,k] <- sx[k]
  F <- matrix(0, k, k); F[1,] <- fx
  list(ages=ages, sx=sx, fx=fx, U=U, F=F, A=U+F)
}

# NOTE: Reverted to sx_senescence / fx_senescence to fix your error
build_MPM_no_senescence <- function(ages, sx_senescence, fx_senescence, senescence_onset_age) {
  # Weighted mean
  w <- calculate_weighted_adult_means(ages, sx_senescence, fx_senescence, senescence_onset_age)
  idx <- which(ages > senescence_onset_age)
  
  sx_no <- sx_senescence
  sx_no[idx] <- w$sx_mean
  
  fx_no <- fx_senescence
  fx_no[idx] <- w$fx_mean
  
  build_MPM_senescence(ages, sx_no, fx_no)
}

build_MPM_no_actuarial_yes_reproductive <- function(ages, sx_senescence, fx_senescence, senescence_onset_age) {
  w <- calculate_weighted_adult_means(ages, sx_senescence, fx_senescence, senescence_onset_age)
  idx <- which(ages > senescence_onset_age)
  
  sx_no <- sx_senescence
  sx_no[idx] <- w$sx_mean
  
  # Fecundity stays as is (senescence)
  build_MPM_senescence(ages, sx_no, fx_senescence)
}

build_MPM_yes_actuarial_no_reproductive <- function(ages, sx_senescence, fx_senescence, senescence_onset_age) {
  w <- calculate_weighted_adult_means(ages, sx_senescence, fx_senescence, senescence_onset_age)
  idx <- which(ages > senescence_onset_age)
  
  fx_no <- fx_senescence
  fx_no[idx] <- w$fx_mean
  
  # Survival stays as is (senescence)
  build_MPM_senescence(ages, sx_senescence, fx_no)
}

compute_summary_table <- function(U_sen, U_no, U_noA_yesR, U_yesA_noR, F_sen, F_no, F_noA_yesR, F_yesA_noR, repro_var="Poisson") {
  if(!exists("mean_lifespan")) stop("Source LuckFunctions.R first!")
  
  .get_moments <- function(U, F) {
    mL <- mean_lifespan(U, mixdist = NULL)[1]
    vL <- var_lifespan(U, mixdist = NULL)[1]
    sL <- skew_lifespan(U, mixdist = NULL)[1]
    mR <- mean_LRO(U, F, mixdist = NULL)[1]
    vR <- var_LRO_mcr(U, F, repro_var, mixdist = NULL)[1]
    sR <- skew_LRO(U, F, repro_var, mixdist = NULL)[1]
    c(mL=as.numeric(mL), vL=as.numeric(vL), sL=as.numeric(sL), mR=as.numeric(mR), vR=as.numeric(vR), sR=as.numeric(sR))
  }
  
  r1 <- .get_moments(U_sen, F_sen)
  r2 <- .get_moments(U_no, F_no)
  r3 <- .get_moments(U_noA_yesR, F_noA_yesR)
  r4 <- .get_moments(U_yesA_noR, F_yesA_noR)
  
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

compute_sensitivity_table <- function(U_sen, U_no_peak, U_no_early, U_no_late, F_sen, F_no_peak, F_no_early, F_no_late, repro_var="Poisson") {
  if(!exists("mean_lifespan")) stop("Source LuckFunctions.R first!")
  
  .get_moments <- function(U, F) {
    mL <- mean_lifespan(U, mixdist = NULL)[1]
    vL <- var_lifespan(U, mixdist = NULL)[1]
    sL <- skew_lifespan(U, mixdist = NULL)[1]
    mR <- mean_LRO(U, F, mixdist = NULL)[1]
    vR <- var_LRO_mcr(U, F, repro_var, mixdist = NULL)[1]
    sR <- skew_LRO(U, F, repro_var, mixdist = NULL)[1]
    c(mL=as.numeric(mL), vL=as.numeric(vL), sL=as.numeric(sL), mR=as.numeric(mR), vR=as.numeric(vR), sR=as.numeric(sR))
  }
  
  r1 <- .get_moments(U_sen, F_sen)
  r2 <- .get_moments(U_no_peak, F_no_peak)
  r3 <- .get_moments(U_no_early, F_no_early)
  r4 <- .get_moments(U_no_late, F_no_late)
  
  data.frame(
    model = c("Senescence", "No-senescence-peak", "No-senescence-early", "No-senescence-late"),
    mean_lifespan = c(r1["mL"], r2["mL"], r3["mL"], r4["mL"]),
    var_lifespan = c(r1["vL"], r2["vL"], r3["vL"], r4["vL"]),
    skew_lifespan = c(r1["sL"], r2["sL"], r3["sL"], r4["sL"]),
    mean_LRO = c(r1["mR"], r2["mR"], r3["mR"], r4["mR"]),
    var_LRO = c(r1["vR"], r2["vR"], r3["vR"], r4["vR"]),
    skew_LRO = c(r1["sR"], r2["sR"], r3["sR"], r4["sR"])
  )
}

