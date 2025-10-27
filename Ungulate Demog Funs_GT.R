#####Derived from Rees et al.'s paper#####


## IBM to generate the data for the simple ungulate example

## Survival kernel: from age x -> x+1
P_x1x <- function(x1, x, sx_vec){
  ifelse(x1 == x + 0,  # indices will be handled in mk_K_age; keep form simple here
         NA, NA)      
}

## Fecundity kernel: all recruits enter age 0
F_x1x <- function(x1, x, fx_vec){
  ifelse(x1 == 0, fx_vec[x + 1], 0)
}

## Build discrete age-structured kernel (Leslie-type)
## ages: integer vector of ages, e.g. 0:A
## sx_vec: length = length(ages), sx_vec[i] = P(age=ages[i] -> ages[i]+1)
## fx_vec: length = length(ages), per-capita fecundity at age i
## retain_last: whether last age retains survivors in last class

mk_K_age <- function(ages, sx_vec, fx_vec, retain_last = FALSE){
  stopifnot(all(ages == seq(min(ages), max(ages))), min(ages) == 0)
  A <- length(ages)
  P <- matrix(0, A, A)
  F <- matrix(0, A, A)
  
  ## fecundity: all newborns to age 0 (row 1)
  F[1, ] <- fx_vec*0.5     #assume half of the offsprings is female
  
  ## survival: subdiagonal P[i+1, i] = s_i  (i index = age i-1)
  for(i in 1:(A-1)){
    P[i+1, i] <- sx_vec[i]
  }
  if (retain_last) {
    P[A, A] <- sx_vec[A]   #whether last age retains survivors in last class
  }
  
  K <- P + F
  list(K = K, P = P, F = F, ages = ages)
}
