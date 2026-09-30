### 01. MULTINOMIAL MODELS ####

#### Reference levels ####

# Use mature natural-reef corals at Koh Tao as the predictor reference group.
condition_data <- condition_data %>%
  mutate(
    condition = relevel(condition, ref = "Healthy"),
    reef_type = relevel(factor(reef_type), ref = "Natural"),
    location = relevel(factor(location), ref = "Koh Tao")
  )

matched_location_data <- condition_data %>%
  filter(bleaching_period == "Middle") %>%
  droplevels()

kt_data <- condition_data %>%
  filter(location == "Koh Tao") %>%
  droplevels()

#### Random-effects structure ####

# CPCE points within quadrats are not independent. Models therefore include a
# quadrat-level random intercept. A diagonal category covariance estimates
# separate quadrat-level variance for each condition-vs-Healthy comparison
# without estimating correlations among those random effects.

random_formula <- ~ 1 | quadrat_id
category_covariance <- "diagonal"

#### M1: matched location comparison ####

# Compare AR recruits and mature NR corals between Koh Tao and Rayong using only
# the temporally matched Middle bleaching period. The reef_type × location
# interaction tests whether the AR-NR condition difference differs between locations.

model_1 <- mclogit::mblogit(
  condition ~ reef_type * location,
  random = random_formula,
  data = matched_location_data,
  estimator = "ML",
  catCov = category_covariance
)

summary(model_1)

##### Model 1 diagnostics ####

model_1_check <- tibble(
  converged = isTRUE(model_1$converged),
  finite_coefficients = all(is.finite(coef(model_1))),
  finite_vcov = all(is.finite(vcov(model_1))),
  residual_deviance = model_1$deviance
)

model_1_check

model_1$VarCov

#### M2: Koh Tao temporal comparison ####

# Test whether condition trajectories through the bleaching event differ between
# AR recruits and mature NR corals, while allowing those trajectories to differ
# between Mango and Tanote.

model_2 <- mclogit::mblogit(
  condition ~ reef_type * bleaching_period * site,
  random = random_formula,
  data = kt_data,
  estimator = "ML",
  catCov = category_covariance
)

summary(model_2)

##### Model 2 diagnostics ####

model_2_check <- tibble(
  converged = isTRUE(model_2$converged),
  finite_coefficients = all(is.finite(coef(model_2))),
  finite_vcov = all(is.finite(vcov(model_2))),
  residual_deviance = model_2$deviance
)

model_2_check



### 02. GENUS-SPECIFIC NONPARAMETRIC TESTS ####
#### Quadrat-level condition composition ####

q3_condition_quadrat <- condition_data %>%
  filter(bleaching_period == "Middle", !is.na(genus)) %>%
  count(location, site, genus, reef_type, quadrat_id, health, name = "points") %>%
  group_by(location, site, genus, reef_type, quadrat_id) %>%
  complete(health = condition_codes, fill = list(points = 0)) %>%
  mutate(proportion = points / sum(points)) %>%
  ungroup() %>%
  select(-points) %>%
  pivot_wider(names_from = health, values_from = proportion, values_fill = 0)


#### Genus support ####

q3_genera <- q3_condition_quadrat %>%
  distinct(location, genus, reef_type) %>%
  count(location, genus, name = "reef_types") %>%
  filter(reef_types == 2) %>%
  select(location, genus)


#### PERMANOVA ####

set.seed(42)
n_perm <- 99999

q3_permanova_test <- function(location_i, genus_i) {
  
  genus_data <- q3_condition_quadrat %>%
    filter(location == location_i, genus == genus_i) %>%
    droplevels()
  
  condition_matrix <- genus_data %>%
    select(H, PBL, FBL, PRK, FRK)
  
  if (n_distinct(genus_data$site) > 1) {
    test <- vegan::adonis2(
      condition_matrix ~ reef_type,
      data = genus_data,
      permutations = n_perm,
      method = "bray",
      strata = genus_data$site
    )
  } else {
    test <- vegan::adonis2(
      condition_matrix ~ reef_type,
      data = genus_data,
      permutations = n_perm,
      method = "bray"
    )
  }
  
  as.data.frame(test) %>%
    rownames_to_column("term") %>%
    slice(1) %>%
    transmute(
      location = location_i,
      genus = genus_i,
      df = Df,
      r2 = R2,
      pseudo_f = F,
      p = `Pr(>F)`
    )
}

q3_results <- map2_dfr(
  q3_genera$location,
  q3_genera$genus,
  q3_permanova_test
) %>%
  group_by(location) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  ungroup()

q3_results %>%
  arrange(location, p_adj) %>%
  print(n = Inf)


#### PERMANOVA dispersion checks ####

q3_dispersion_test <- function(location_i, genus_i) {
  
  genus_data <- q3_condition_quadrat %>%
    filter(location == location_i, genus == genus_i) %>%
    droplevels()
  
  condition_dist <- vegan::vegdist(
    genus_data %>% select(H, PBL, FBL, PRK, FRK),
    method = "bray"
  )
  
  dispersion <- vegan::betadisper(
    condition_dist,
    genus_data$reef_type,
    add = "lingoes"
  )
  
  test <- vegan::permutest(
    dispersion,
    permutations = n_perm
  )
  
  tibble(
    location = location_i,
    genus = genus_i,
    dispersion_f = test$tab[1, "F"],
    dispersion_p = test$tab[1, "Pr(>F)"]
  )
}

q3_dispersion <- map2_dfr(
  q3_genera$location,
  q3_genera$genus,
  q3_dispersion_test
) %>%
  group_by(location) %>%
  mutate(dispersion_p_adj = p.adjust(dispersion_p, method = "BH")) %>%
  ungroup()

q3_dispersion %>%
  arrange(location, dispersion_p_adj) %>%
  print(n = Inf)

q3_results <- q3_results %>%
  left_join(q3_dispersion, by = c("location", "genus"))

q3_results %>%
  arrange(location, p_adj) %>%
  print(n = Inf)

#### b: Quadrat-level susceptibility ####

genus_quadrat <- condition_data %>%
  filter(!is.na(genus)) %>%
  group_by(location, site, genus, reef_type, bleaching_period, quadrat_id) %>%
  summarise(
    coral_points = n(),
    affected_prop = mean(affected),
    .groups = "drop"
  )

# check 


#### b-Q3: AR vs NR within genus and location ####

n_perm <- 99999 
set.seed(42)


permute_reef_difference <- function(df, n_perm) {
  reef_type <- as.character(df$reef_type)
  affected_prop <- df$affected_prop
  
  n_artificial <- sum(reef_type == "Artificial")
  n_natural <- sum(reef_type == "Natural")
  
  artificial_mean <- mean(affected_prop[reef_type == "Artificial"])
  natural_mean <- mean(affected_prop[reef_type == "Natural"])
  observed <- artificial_mean - natural_mean
  
  if (n_distinct(affected_prop) == 1) {
    return(tibble(
      artificial_quadrats = n_artificial,
      natural_quadrats = n_natural,
      artificial_mean = artificial_mean,
      natural_mean = natural_mean,
      difference = observed,
      p = 1
    ))
  }
  
  strata <- split(df, df$site)
  y <- lapply(strata, \(x) x$affected_prop)
  n_artificial_site <- vapply(
    strata,
    \(x) sum(x$reef_type == "Artificial"),
    integer(1)
  )
  
  total_affected <- sum(affected_prop)
  
  permuted <- replicate(n_perm, {
    artificial_sum <- sum(vapply(
      seq_along(y),
      \(i) sum(sample(y[[i]], n_artificial_site[i])),
      numeric(1)
    ))
    
    artificial_sum / n_artificial -
      (total_affected - artificial_sum) / n_natural
  })
  
  tibble(
    artificial_quadrats = n_artificial,
    natural_quadrats = n_natural,
    artificial_mean = artificial_mean,
    natural_mean = natural_mean,
    difference = observed,
    p = (sum(abs(permuted) >= abs(observed)) + 1) / (n_perm + 1)
  )
}

q3_tests <- genus_quadrat %>%
  filter(bleaching_period == "Middle") %>%
  group_by(location, genus) %>%
  filter(n_distinct(reef_type) == 2) %>%
  group_modify(~ permute_reef_difference(.x, n_perm)) %>%
  ungroup() %>%
  group_by(location) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  ungroup()

q3_tests %>%
  arrange(location, p_adj) %>%
  print(n = Inf)

### 03. Q4 GENUS-SPECIFIC CONDITION TRAJECTORIES ####

#### Quadrat-level five-state composition ####

q4_condition_quadrat <- condition_data %>%
  filter(location == "Koh Tao", !is.na(genus)) %>%
  count(
    site, genus, reef_type, bleaching_period, quadrat_id, health,
    name = "points"
  ) %>%
  group_by(site, genus, reef_type, bleaching_period, quadrat_id) %>%
  complete(health = condition_codes, fill = list(points = 0)) %>%
  mutate(proportion = points / sum(points)) %>%
  ungroup() %>%
  select(-points) %>%
  pivot_wider(
    names_from = health,
    values_from = proportion,
    values_fill = 0
  )


#### Genus support ####

# The interaction requires both reef types to be represented in multiple periods.
q4_support <- q4_condition_quadrat %>%
  count(genus, bleaching_period, reef_type, name = "quadrats") %>%
  group_by(genus, bleaching_period) %>%
  summarise(
    both_reef_types = n_distinct(reef_type) == 2,
    min_quadrats = min(quadrats),
    .groups = "drop"
  ) %>%
  group_by(genus) %>%
  summarise(
    periods = n_distinct(bleaching_period),
    periods_both_reef_types = sum(both_reef_types),
    min_quadrats = min(min_quadrats),
    .groups = "drop"
  ) %>%
  arrange(desc(periods_both_reef_types), desc(min_quadrats))

q4_support %>% print(n = Inf)

q4_genera <- q4_support %>%
  filter(periods_both_reef_types >= 2) %>%
  pull(genus)


#### Five-state trajectories ####

# Means are calculated across quadrats so CPCE points remain subsamples rather
# than being treated as independent biological replicates.
q4_trajectories <- q4_condition_quadrat %>%
  filter(genus %in% q4_genera) %>%
  group_by(genus, reef_type, bleaching_period) %>%
  summarise(
    quadrats = n(),
    across(c(H, PBL, FBL, PRK, FRK), mean),
    .groups = "drop"
  )

q4_trajectories %>%
  arrange(genus, reef_type, bleaching_period) %>%
  print(n = Inf)


#### PERMANOVA ####

set.seed(42)
n_perm <- 99999

# Restrict each genus to bleaching periods where both reef types are represented.
q4_permanova_data <- q4_condition_quadrat %>%
  filter(genus %in% q4_genera) %>%
  group_by(genus, bleaching_period) %>%
  filter(n_distinct(reef_type) == 2) %>%
  ungroup()

q4_permanova <- set_names(q4_genera) %>%
  map(function(g) {
    
    genus_data <- q4_permanova_data %>%
      filter(genus == g) %>%
      droplevels()
    
    condition_matrix <- genus_data %>%
      select(H, PBL, FBL, PRK, FRK)
    
    vegan::adonis2(
      condition_matrix ~ site + reef_type * bleaching_period,
      data = genus_data,
      permutations = n_perm,
      method = "bray",
      strata = genus_data$site,
      by = "margin"
    )
  })

q4_permanova_results <- imap_dfr(
  q4_permanova,
  ~ as.data.frame(.x) %>%
    rownames_to_column("term") %>%
    filter(term %in% c(
      "reef_type",
      "bleaching_period",
      "reef_type:bleaching_period"
    )) %>%
    transmute(
      genus = .y,
      term,
      df = Df,
      r2 = R2,
      pseudo_f = F,
      p = `Pr(>F)`
    )
) %>%
  group_by(term) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  ungroup()

q4_permanova_results %>%
  arrange(term, p_adj) %>%
  print(n = Inf)

#### PERMANOVA dispersion checks ####

q4_dispersion <- set_names(q4_genera) %>%
  map(function(g) {
    
    genus_data <- q4_permanova_data %>%
      filter(genus == g) %>%
      droplevels()
    
    condition_dist <- vegan::vegdist(
      genus_data %>% select(H, PBL, FBL, PRK, FRK),
      method = "bray"
    )
    
    group <- interaction(
      genus_data$reef_type,
      genus_data$bleaching_period,
      drop = TRUE
    )
    
    dispersion <- vegan::betadisper(condition_dist, group)
    
    list(
      dispersion = dispersion,
      test = vegan::permutest(dispersion, permutations = n_perm)
    )
  })

