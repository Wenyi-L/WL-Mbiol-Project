## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Individual-based model (IBM)
## Simplified version: directly uses externally provided fx_vec (breeding probability)
## Required external inputs:
##   sx_vec        : age-specific survival probabilities
##   fx_vec        : age-specific breeding probabilities
##   mean_chick    : SINGLE scalar mean clutch size for breeding individuals (global)
##   great_tit     : data.frame with columns x (age) and Nx (abundance) for stable age distribution
##   init.pop.size : initial population size
##   n.yrs         : number of years to simulate
## Optional inputs:
##   prop_female   : probability of being female at birth (default = 0.5)
##   max.pop.size  : population size ceiling (default = 8000)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## (0) Initialize stable age distribution
stable_age_dist <- great_tit$Nx / sum(great_tit$Nx)
ages <- great_tit$x
pop <- sample(ages, size = init.pop.size, replace = TRUE, prob = stable_age_dist)

## (1) Initialize parameters and storage
pop.size.t <- numeric(n.yrs)
prop_female <- 0.5

## (2) Initialize individual attributes
sex <- rbinom(n = length(pop), size = 1, prob = prop_female)   # 1=female, 0=male
id  <- seq_len(init.pop.size)                                  # numeric unique ID
name <- paste0("ind_", id)                                     # readable label
mother_id <- rep(NA_integer_, init.pop.size)                   # founders have no mothers
next_id <- init.pop.size + 1L

## (3) Create historical data frame
sim.data <- data.frame(
    Year       = integer(),
    ID         = integer(),
    Name       = character(),
    Mother_ID  = integer(),
    Age        = integer(),
    Sex        = integer(),
    Surv       = integer(),
    Offspring  = integer(),
    stringsAsFactors = FALSE
)

## (4) Main simulation loop
yr <- 1
while (yr <= n.yrs && length(pop) < 8000) {
    
    ## (A) Current population size
    pop.size <- length(pop)
    pop.size.t[yr] <- pop.size
    
    ## (B) Survival process
    A <- length(sx_vec)
    p_surv <- ifelse(pop < (A - 1), sx_vec[pop + 1], sx_vec[A])
    surv <- rbinom(n = pop.size, size = 1, prob = p_surv)
    pop_surv_next <- pop[surv == 1] + 1
    
    ## (C) Reproduction process (uses externally defined fx_vec)
    age_idx_rep <- ifelse(pop < (A - 1), pop + 1, A)
    idx_female_surv <- which(surv == 1 & sex == 1)
    
    ## (C1) Breeding probability (from fx_vec)
    p_breed <- if (length(idx_female_surv) > 0) {
        pmin(pmax(fx_vec[age_idx_rep[idx_female_surv]], 0), 1)
    } else numeric(0)
    
    ## (C2) Draw whether each surviving female breeds (Bernoulli)
    do_breed <- if (length(idx_female_surv) > 0) {
        rbinom(length(idx_female_surv), size = 1, prob = p_breed)
    } else integer(0)
    
    idx_breed <- if (length(idx_female_surv) > 0) idx_female_surv[which(do_breed == 1)] else integer(0)
    
    ## (C3) Draw clutch size for breeders (Poisson, zero-truncated) — USE GLOBAL mean_chick
    n_offspring_vec <- integer(length(idx_female_surv))
    if (length(idx_breed) > 0) {
        lambda_val <- max(mean_chick, 0)                              # scalar guard (non-negative)
        draws <- rpois(length(idx_breed), lambda = lambda_val)
        draws[draws < 1] <- 1L                                         # zero-truncation
        n_offspring_vec[match(idx_breed, idx_female_surv)] <- draws
    }
    
    ## (C4) Record offspring count for each female (0 for non-breeders)
    offspring_vec <- integer(pop.size)
    if (length(idx_female_surv) > 0) offspring_vec[idx_female_surv] <- n_offspring_vec
    
    ## (D) Generate recruits (age=0) with IDs and mother links
    n_recruits   <- sum(n_offspring_vec)
    pop_recruits <- rep(0L, n_recruits)
    sex_recruits <- if (n_recruits > 0) rbinom(n_recruits, size = 1, prob = prop_female) else integer(0)
    id_recruits  <- if (n_recruits > 0) seq.int(next_id, length.out = n_recruits) else integer(0)
    name_recruits <- character(n_recruits)
    mother_id_recruits <- integer(n_recruits)
    
    if (n_recruits > 0) {
        ptr <- 1L
        for (k in seq_along(idx_female_surv)) {
            i <- idx_female_surv[k]
            cnt <- n_offspring_vec[k]
            if (cnt > 0) {
                mom_id <- id[i]
                for (j in seq_len(cnt)) {
                    name_recruits[ptr]      <- paste0("chick_", yr, "_mom", mom_id, "_", j)
                    mother_id_recruits[ptr] <- mom_id
                    ptr <- ptr + 1L
                }
            }
        }
    }
    next_id <- next_id + n_recruits
    
    ## (E) Append this year's individuals to sim.data
    ## (E1) Existing individuals
    if (pop.size > 0) {
        sim.data <- rbind(
            sim.data,
            data.frame(
                Year       = rep(yr, pop.size),
                ID         = id,
                Name       = name,
                Mother_ID  = mother_id,
                Age        = pop,
                Sex        = sex,
                Surv       = surv,
                Offspring  = offspring_vec,
                stringsAsFactors = FALSE
            )
        )
    }
    
    ## (E2) Newborns
    if (n_recruits > 0) {
        sim.data <- rbind(
            sim.data,
            data.frame(
                Year       = rep(yr, n_recruits),
                ID         = id_recruits,
                Name       = name_recruits,
                Mother_ID  = mother_id_recruits,
                Age        = pop_recruits,
                Sex        = sex_recruits,
                Surv       = NA_integer_,
                Offspring  = 0L,
                stringsAsFactors = FALSE
            )
        )
    }
    
    ## (F) Update population vectors for next year
    pop       <- c(pop_recruits,       pop_surv_next)
    sex       <- c(sex_recruits,       sex[surv == 1])
    id        <- c(id_recruits,        id[surv == 1])
    name      <- c(name_recruits,      name[surv == 1])
    mother_id <- c(mother_id_recruits, mother_id[surv == 1])
    
    ## Advance year
    yr <- yr + 1
}

## (G) Final cleanup
pop.size.t <- pop.size.t[pop.size.t > 0]

## sim.data contains all individuals ever existed:
## Columns:
##   Year, ID, Name, Mother_ID, Age, Sex, Surv, Offspring
