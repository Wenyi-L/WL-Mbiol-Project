## ==========================================================
## Function: plateau_after_maturity
## Purpose : Create a "counterfactual" version of fx_vec or sx_vec
##           where values remain constant after sexual maturity.
## Inputs  :
##   v              - numeric vector (fx_vec or sx_vec)
##   maturity_age   - the age (integer) at which sexual maturity is reached
##   is_survival    - logical; TRUE if v is survival probabilities (sx_vec)
## Returns :
##   numeric vector where all ages >= maturity_age are constant at the
##   maturity-age value.
## Notes   :
##   - Age is assumed to start at 0 (age 0 = newborn).
##   - Values before maturity are unchanged.
##   - Non-finite or negative values are replaced with 0.
##   - Survival values are truncated to [0, 1].
## ==========================================================

plateau_after_maturity <- function(v, maturity_age, is_survival = FALSE) {
  v <- as.numeric(v)
  A <- length(v)
  ages <- seq_len(A) - 1L  # assume age 0 corresponds to v[1]
  
  ## Clean input values
  if (is_survival) {
    v[!is.finite(v) | v < 0 | v > 1] <- 0
  } else {
    v[!is.finite(v) | v < 0] <- 0
  }
  
  ## Basic bounds check for maturity_age
  if (!is.finite(maturity_age) || maturity_age < 0) maturity_age <- 0L
  if (maturity_age > max(ages)) {
    # if maturity_age exceeds the defined range, return v unchanged
    return(v)
  }
  
  ## Identify the maturity-age index and its constant value
  idx_m <- which(ages == maturity_age)
  if (length(idx_m) == 0) idx_m <- A
  const_val <- v[idx_m]
  
  ## Optional fallback: if maturity-age value is 0 or NA, 
  ## use the most recent finite/nonzero value before maturity.
  # if (!is.finite(const_val) || const_val == 0) {
  #   prev_ok <- tail(which(is.finite(v) & v > 0 & ages <= maturity_age), 1)
  #   if (length(prev_ok) == 1) const_val <- v[prev_ok]
  # }
  
  ## Replace all post-maturity ages (>= maturity_age) with const_val
  v[ages >= maturity_age] <- const_val
  
  ## Ensure survival values remain within [0,1]
  if (is_survival) v <- pmin(pmax(v, 0), 1)
  
  return(v)
}
