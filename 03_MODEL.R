### 03. MULTINOMIAL MODELS ####

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

#### Model 1: matched location comparison ####

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

#### Model 2: Koh Tao temporal comparison ####

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

#### Model 3: genus-specific matched location comparison ####

model_3_genera <- c(
  "Porites",
  "Platygyra",
  "Pocillopora",
  "Favites",
  "Goniastrea",
  "Dipsastrea",
  "Lobophyllia"
)

genus_location_quadrat <- matched_location_data %>%
  filter(genus %in% model_3_genera) %>%
  group_by(genus, location, reef_type, quadrat_id) %>%
  summarise(
    affected_points = sum(affected == 1),
    healthy_points = sum(affected == 0),
    .groups = "drop"
  )

model_3 <- map(
  model_3_genera,
  function(g) {
    
    genus_data <- genus_location_quadrat %>%
      filter(genus == g) %>%
      droplevels()
    
    glmmTMB(
      cbind(affected_points, healthy_points) ~ reef_type * location,
      family = betabinomial(link = "logit"),
      data = genus_data
    )
  }
)

names(model_3) <- model_3_genera

##### Model 3 diagnostics ####
# check fits 
model_3_checks <- imap_dfr(
  model_3,
  ~ tibble(
    genus = .y,
    convergence_code = .x$fit$convergence,
    pd_hessian = .x$sdr$pdHess,
    finite_coefficients = all(is.finite(fixef(.x)$cond)),
    max_abs_beta = max(abs(fixef(.x)$cond))
  )
)

model_3_checks
# inspect dispersion 
imap_dfr(
  model_3,
  ~ tibble(
    genus = .y,
    dispersion = sigma(.x)
  )
)

# print summaries 
imap(
  model_3,
  ~ {
    cat("\n\n###", .y, "###\n")
    print(summary(.x))
  }
)

#### Koh Tao site consistency ####

kt_middle_quadrat <- matched_location_data %>%
  filter(location == "Koh Tao") %>%
  group_by(site, reef_type, quadrat_id) %>%
  summarise(
    affected_points = sum(affected == 1),
    healthy_points = sum(affected == 0),
    .groups = "drop"
  )

kt_site_check <- glmmTMB(
  cbind(affected_points, healthy_points) ~ reef_type * site,
  family = betabinomial(link = "logit"),
  data = kt_middle_quadrat
)

summary(kt_site_check)
#### Model 4: genus-specific Koh Tao temporal susceptibility ####

model_4_genera <- candidate_genera

genus_temporal_quadrat <- kt_data %>%
  filter(genus %in% model_4_genera) %>%
  group_by(genus, site, reef_type, bleaching_period, quadrat_id) %>%
  summarise(
    affected_points = sum(affected == 1),
    healthy_points = sum(affected == 0),
    .groups = "drop"
  )

model_4 <- map(
  model_4_genera,
  function(g) {
    
    genus_data <- genus_temporal_quadrat %>%
      filter(genus == g) %>%
      droplevels()
    
    glmmTMB(
      cbind(affected_points, healthy_points) ~
        reef_type * bleaching_period + site,
      family = betabinomial(link = "logit"),
      data = genus_data
    )
  }
)

names(model_4) <- model_4_genera


##### Model 4 diagnostics ####

model_4_checks <- imap_dfr(
  model_4,
  ~ tibble(
    genus = .y,
    convergence_code = .x$fit$convergence,
    pd_hessian = .x$sdr$pdHess,
    finite_coefficients = all(is.finite(fixef(.x)$cond)),
    max_abs_beta = max(abs(fixef(.x)$cond)),
    dispersion = sigma(.x)
  )
)

model_4_checks
