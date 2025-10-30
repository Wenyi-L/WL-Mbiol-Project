## ----------------------------------------------------------
## Compute realized LRO per mother
## Definition: number of offspring that survived to age 1
## Includes non-breeding females (LRO = 0)
## ----------------------------------------------------------
compute_LRO_realized <- function(sim.data) {
  # All adult females in the dataset
  all_females <- unique(sim.data$ID[sim.data$Sex == 1])
  
  # Records of individuals that have a known mother (i.e., are offspring)
  offspring_records <- sim.data[!is.na(sim.data$Mother_ID), ]
  
  # Identify offspring who survived to at least age 1
  survived_offspring <- offspring_records %>%
    group_by(ID) %>%
    summarise(max_age = max(Age, na.rm = TRUE), .groups = "drop") %>%
    filter(max_age >= 1)
  
  # Count how many of those belong to each mother
  LRO_df <- offspring_records %>%
    filter(ID %in% survived_offspring$ID) %>%
    group_by(Mother_ID) %>%
    summarise(LRO = n(), .groups = "drop")
  
  # Mothers who produced offspring but none survived → LRO = 0
  all_mothers_with_chicks <- unique(offspring_records$Mother_ID)
  zero_LRO_moms <- setdiff(all_mothers_with_chicks, LRO_df$Mother_ID)
  if (length(zero_LRO_moms) > 0) {
    LRO_df <- bind_rows(LRO_df, data.frame(Mother_ID = zero_LRO_moms, LRO = 0))
  }
  
  # Females who never reproduced → also assign LRO = 0
  zero_repro_females <- setdiff(all_females, all_mothers_with_chicks)
  if (length(zero_repro_females) > 0) {
    LRO_df <- bind_rows(LRO_df, data.frame(Mother_ID = zero_repro_females, LRO = 0))
  }
  
  # Return tidy dataframe (one row per female)
  LRO_df
}
