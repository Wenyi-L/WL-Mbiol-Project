## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 -
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define age range based on survival vector
stable_age_dist <- great_tit$Nx / sum(great_tit$Nx)
ages <- great_tit$x


## initial size distribution

pop <- sample(ages, size = init.pop.size, replace = TRUE, prob = stable_age_dist)


## vectors to store pop size and mean size

pop.size.t <- numeric(n.yrs)

## Iterate the model using the "true" parameters and store data in a data.frame
## ages start with 0
# sx_vec=P(survive x→x+1)
# fx_vec为Nx/noffsprings
## only female reproduce!

#Randomly assign gender to individuals
#1=female, 0=male
prop_female <- 0.5
sex <- rbinom(n = length(pop), size = 1, prob = prop_female)

last_age <- last_sex <- last_surv <- last_offspring <- NULL

id <- seq_len(init.pop.size)
next_id <- init.pop.size + 1L

death_records <- data.frame(id = integer(),
                            Age_at_death = integer(),
                            Sex = integer(),
                            stringsAsFactors = FALSE)

yr <- 1
while (yr <= n.yrs && length(pop) < 8000) {
    
    pop.size <- length(pop)
    pop.size.t[yr] <- pop.size
    
    ## survival
    A <- length(sx_vec)
    p_surv <- ifelse(pop < (A-1), sx_vec[pop + 1], sx_vec[A])
    surv   <- rbinom(n = length(pop), size = 1, prob = p_surv)
    
    ## age+1 for survivors
    pop_surv_next <- pop[surv == 1] + 1
    
    ## reproduction
    age_idx_rep <- ifelse(pop < (A-1), pop + 1, A)
    idx_female_surv <- which(surv == 1 & sex == 1)
    lam <- fx_vec[ age_idx_rep[idx_female_surv] ]
    n_offspring <- rpois(length(idx_female_surv), lambda = pmax(lam, 0))
    
    ## for new recruits, age = 0
    n_recruits   <- sum(n_offspring)
    pop_recruits <- rep(0L, n_recruits)
    sex_recruits <- rbinom(n_recruits, size = 1, prob = prop_female)
    id_recruits  <- if (n_recruits > 0) seq.int(next_id, length.out = n_recruits) else integer(0)
    next_id <- next_id + n_recruits
    
    ## temporarily save current snapshot (like original sim.data)
    offspring_vec <- integer(pop.size)
    if (length(idx_female_surv) > 0) offspring_vec[idx_female_surv] <- n_offspring
    
    last_age       <- pop
    last_sex       <- sex
    last_surv      <- surv
    last_offspring <- offspring_vec
    
    ##save death record
    idx_dead <- which(surv == 0)
    if (length(idx_dead) > 0) {
        deaths_this_year <- data.frame(
            id            = id[idx_dead],
            Age_at_death  = pop[idx_dead],
            Sex           = sex[idx_dead]
        )
        death_records <- rbind(death_records, deaths_this_year)
    }
    
    ## next year
    pop <- c(pop_recruits, pop_surv_next)
    sex <- c(sex_recruits, sex[surv == 1])
    id  <- c(id_recruits,  id[surv == 1])
    
    yr <- yr + 1
}


## trim the population size vector to remove the zeros at the end
pop.size.t <- pop.size.t[pop.size.t>0]

#save sim.data
sim.data <- data.frame(
    Age       = last_age,
    Sex       = last_sex,
    Surv      = last_surv,
    Offspring = last_offspring
)

