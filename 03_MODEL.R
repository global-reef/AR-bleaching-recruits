###  Q1 MATCHED LOCATION COMPARISON ####

#### Data and reference levels ####

# Healthy, mature natural-reef corals at Koh Tao define the reference condition.
condition_data <- condition_data %>%
  mutate(
    condition = relevel(condition, ref = "Healthy"),
    reef_type = relevel(factor(reef_type), ref = "Natural"),
    location = relevel(factor(location), ref = "Koh Tao")
  )

q1_data <- condition_data %>%
  filter(bleaching_period == "Middle") %>%
  droplevels()

#### Multinomial mixed model ####

# CPCE points are subsamples within quadrats; quadrat_id accounts for clustering.
q1_model <- mclogit::mblogit(
  condition ~ reef_type * location,
  random = ~ 1 | quadrat_id,
  data = q1_data,
  estimator = "ML",
  catCov = "diagonal"
)

q1_diagnostics <- tibble(
  converged = isTRUE(q1_model$converged),
  finite_coefficients = all(is.finite(coef(q1_model))),
  finite_vcov = all(is.finite(vcov(q1_model))),
  residual_deviance = q1_model$deviance
)

saveRDS(q1_model, file.path(fits_dir, "q1_multinomial_model.rds"))
write_csv(q1_diagnostics, file.path(stats_dir, "q1_model_diagnostics.csv"))

summary(q1_model)
q1_diagnostics
q1_model$VarCov


### Q2 KOH TAO TEMPORAL COMPARISON ####

#### Data ####

q2_data <- condition_data %>%
  filter(location == "Koh Tao") %>%
  droplevels()

#### Multinomial mixed model ####

# The three-way interaction tests whether AR-NR temporal patterns differ by site.
q2_model <- mclogit::mblogit(
  condition ~ reef_type * bleaching_period * site,
  random = ~ 1 | quadrat_id,
  data = q2_data,
  estimator = "ML",
  catCov = "diagonal"
)

q2_diagnostics <- tibble(
  converged = isTRUE(q2_model$converged),
  finite_coefficients = all(is.finite(coef(q2_model))),
  finite_vcov = all(is.finite(vcov(q2_model))),
  residual_deviance = q2_model$deviance
)

saveRDS(q2_model, file.path(fits_dir, "q2_multinomial_model.rds"))
write_csv(q2_diagnostics, file.path(stats_dir, "q2_model_diagnostics.csv"))

summary(q2_model)
q2_diagnostics
q2_model$VarCov


###  Q3a GENUS-SPECIFIC MATCHED-PEAK CONDITION ####

#### Quadrat-level condition composition ####

condition_cols <- c("H", "PBL", "FBL", "PRK", "FRK")
n_perm <- 99999
set.seed(42)

q3_condition_quadrat <- condition_data %>%
  filter(bleaching_period == "Middle", !is.na(genus)) %>%
  count(location, site, genus, reef_type, quadrat_id, health, name = "points") %>%
  group_by(location, site, genus, reef_type, quadrat_id) %>%
  complete(health = condition_codes, fill = list(points = 0)) %>%
  mutate(
    coral_points = sum(points),
    proportion = points / coral_points
  ) %>%
  ungroup() %>%
  select(-points) %>%
  pivot_wider(names_from = health, values_from = proportion, values_fill = 0)

#### Genus support ####

q3_support <- q3_condition_quadrat %>%
  group_by(location, genus, reef_type) %>%
  summarise(
    quadrats = n(),
    min_coral_points = min(coral_points),
    median_coral_points = median(coral_points),
    max_coral_points = max(coral_points),
    .groups = "drop"
  )

q3_genera <- q3_condition_quadrat %>%
  distinct(location, genus, reef_type) %>%
  count(location, genus, name = "reef_types") %>%
  filter(reef_types == 2) %>%
  select(location, genus)


#### PERMANOVA ####

q3_permanova_test <- function(location_i, genus_i) {
  
  genus_data <- q3_condition_quadrat %>%
    filter(location == location_i, genus == genus_i) %>%
    droplevels()
  
  condition_matrix <- genus_data %>%
    select(all_of(condition_cols))
  
  if (n_distinct(genus_data$site, na.rm = TRUE) > 1) {
    test <- vegan::adonis2(
      condition_matrix ~ site + reef_type,
      data = genus_data,
      permutations = n_perm,
      method = "bray",
      strata = genus_data$site,
      by = "margin"
    )
  } else {
    test <- vegan::adonis2(
      condition_matrix ~ reef_type,
      data = genus_data,
      permutations = n_perm,
      method = "bray",
      by = "margin"
    )
  }
  
  as.data.frame(test) %>%
    rownames_to_column("term") %>%
    filter(term == "reef_type") %>%
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

#### Dispersion diagnostics ####

q3_dispersion_test <- function(location_i, genus_i) {
  
  genus_data <- q3_condition_quadrat %>%
    filter(location == location_i, genus == genus_i) %>%
    droplevels()
  
  dispersion <- genus_data %>%
    select(all_of(condition_cols)) %>%
    vegan::vegdist(method = "bray") %>%
    vegan::betadisper(genus_data$reef_type, add = "lingoes")
  
  test <- vegan::permutest(dispersion, permutations = n_perm)
  
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

q3_results <- q3_results %>%
  left_join(q3_dispersion, by = c("location", "genus"))

write_csv(q3_support, file.path(stats_dir, "q3_genus_support.csv"))
write_csv(q3_results, file.path(stats_dir, "q3_permanova_results.csv"))


q3_results %>%
  arrange(location, p_adj) %>%
  print(n = Inf)



### Q3b GENUS-SPECIFIC SITE DIFFERENCES ####

#### Genus selection ####

q3b_top_genera <- q3_condition_quadrat %>%
  group_by(genus) %>%
  summarise(coral_points = sum(coral_points), .groups = "drop") %>%
  arrange(desc(coral_points)) %>%
  slice_head(n = 10) %>%
  pull(genus)

q3_condition_quadrat %>%
  filter(genus %in% q3b_top_genera) %>%
  group_by(genus, reef_type, site) %>%
  summarise(
    quadrats = n(),
    coral_points = sum(coral_points),
    .groups = "drop"
  ) %>%
  arrange(match(genus, q3b_top_genera), reef_type, site) %>%
  print(n = Inf)

##### Genus support ####

q3b_support <- q3_condition_quadrat %>%
  filter(genus %in% q3b_top_genera) %>%
  count(genus, reef_type, site, name = "quadrats") %>%
  filter(quadrats >= 5)

q3b_genera <- q3b_support %>%
  count(genus, reef_type, name = "sites") %>%
  filter(sites >= 2) %>%
  select(genus, reef_type)

q3b_support %>%
  semi_join(q3b_genera, by = c("genus", "reef_type")) %>%
  arrange(match(genus, q3b_top_genera), reef_type, site) %>%
  print(n = Inf)

#### PERMANOVA ####

q3b_permanova_test <- function(genus_i, reef_type_i) {
  
  genus_data <- q3_condition_quadrat %>%
    filter(genus == genus_i, reef_type == reef_type_i) %>%
    semi_join(
      q3b_support %>%
        filter(genus == genus_i, reef_type == reef_type_i),
      by = c("genus", "reef_type", "site")
    ) %>%
    droplevels()
  
  condition_matrix <- genus_data %>%
    select(all_of(condition_cols))
  
  test <- vegan::adonis2(
    condition_matrix ~ site,
    data = genus_data,
    permutations = n_perm,
    method = "bray",
    by = "margin"
  )
  
  as.data.frame(test) %>%
    rownames_to_column("term") %>%
    filter(term == "site") %>%
    transmute(
      genus = genus_i,
      reef_type = reef_type_i,
      df = Df,
      r2 = R2,
      pseudo_f = F,
      p = `Pr(>F)`
    )
}

q3b_results <- map2_dfr(
  q3b_genera$genus,
  q3b_genera$reef_type,
  q3b_permanova_test
) %>%
  group_by(reef_type) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  ungroup()


#### Dispersion diagnostics ####

q3b_dispersion_test <- function(genus_i, reef_type_i) {
  
  genus_data <- q3_condition_quadrat %>%
    filter(genus == genus_i, reef_type == reef_type_i) %>%
    semi_join(
      q3b_support %>%
        filter(genus == genus_i, reef_type == reef_type_i),
      by = c("genus", "reef_type", "site")
    ) %>%
    droplevels()
  
  dispersion <- genus_data %>%
    select(all_of(condition_cols)) %>%
    vegan::vegdist(method = "bray") %>%
    vegan::betadisper(genus_data$site, add = "lingoes")
  
  test <- vegan::permutest(dispersion, permutations = n_perm)
  
  tibble(
    genus = genus_i,
    reef_type = reef_type_i,
    dispersion_f = test$tab[1, "F"],
    dispersion_p = test$tab[1, "Pr(>F)"]
  )
}

q3b_dispersion <- map2_dfr(
  q3b_genera$genus,
  q3b_genera$reef_type,
  q3b_dispersion_test
) %>%
  group_by(reef_type) %>%
  mutate(dispersion_p_adj = p.adjust(dispersion_p, method = "BH")) %>%
  ungroup()

q3b_results <- q3b_results %>%
  left_join(q3b_dispersion, by = c("genus", "reef_type"))

q3b_results


#### Pairwise site follow-ups ####

q3b_pairwise <- q3b_results %>%
  filter(df > 1, p_adj < 0.05) %>%
  select(genus, reef_type) %>%
  pmap_dfr(function(genus, reef_type) {
    
    genus_data <- q3_condition_quadrat %>%
      filter(genus == !!genus, reef_type == !!reef_type) %>%
      semi_join(
        q3b_support %>%
          filter(genus == !!genus, reef_type == !!reef_type),
        by = c("genus", "reef_type", "site")
      ) %>%
      droplevels()
    
    combn(unique(genus_data$site), 2, simplify = FALSE) %>%
      map_dfr(function(site_pair) {
        
        pair_data <- genus_data %>%
          filter(site %in% site_pair) %>%
          droplevels()
        
        test <- vegan::adonis2(
          pair_data %>% select(all_of(condition_cols)) ~ site,
          data = pair_data,
          permutations = n_perm,
          method = "bray",
          by = "margin"
        )
        
        as.data.frame(test) %>%
          rownames_to_column("term") %>%
          filter(term == "site") %>%
          transmute(
            genus = genus,
            reef_type = reef_type,
            site_1 = site_pair[1],
            site_2 = site_pair[2],
            df = Df,
            r2 = R2,
            pseudo_f = F,
            p = `Pr(>F)`
          )
      })
  }) %>%
  group_by(genus, reef_type) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  ungroup()

q3b_pairwise %>%
  arrange(genus, reef_type, p_adj) %>%
  print(n = Inf)

#### Pairwise dispersion diagnostics ####

q3b_pairwise_dispersion <- q3b_pairwise %>%
  select(genus, reef_type, site_1, site_2) %>%
  pmap_dfr(function(genus, reef_type, site_1, site_2) {
    
    pair_data <- q3_condition_quadrat %>%
      filter(
        genus == !!genus,
        reef_type == !!reef_type,
        site %in% c(site_1, site_2)
      ) %>%
      semi_join(
        q3b_support %>%
          filter(
            genus == !!genus,
            reef_type == !!reef_type,
            site %in% c(site_1, site_2)
          ),
        by = c("genus", "reef_type", "site")
      ) %>%
      droplevels()
    
    dispersion <- pair_data %>%
      select(all_of(condition_cols)) %>%
      vegan::vegdist(method = "bray") %>%
      vegan::betadisper(pair_data$site, add = "lingoes")
    
    test <- vegan::permutest(dispersion, permutations = n_perm)
    
    tibble(
      genus = genus,
      reef_type = reef_type,
      site_1 = site_1,
      site_2 = site_2,
      dispersion_f = test$tab[1, "F"],
      dispersion_p = test$tab[1, "Pr(>F)"]
    )
  }) %>%
  group_by(genus, reef_type) %>%
  mutate(dispersion_p_adj = p.adjust(dispersion_p, method = "BH")) %>%
  ungroup()

q3b_pairwise %>%
  left_join(
    q3b_pairwise_dispersion,
    by = c("genus", "reef_type", "site_1", "site_2")
  ) %>%
  arrange(genus, reef_type, p_adj) %>%
  print(n = Inf)

q3b_pairwise_results <- q3b_pairwise %>%
  left_join(
    q3b_pairwise_dispersion,
    by = c("genus", "reef_type", "site_1", "site_2")
  )

write_csv(q3b_support, file.path(stats_dir, "q3b_genus_site_support.csv"))
write_csv(q3b_results, file.path(stats_dir, "q3b_permanova_results.csv"))
write_csv(q3b_pairwise_results, file.path(stats_dir, "q3b_pairwise_results.csv"))

### Q4 GENUS-SPECIFIC TEMPORAL CONDITION ####

#### Quadrat-level condition composition ####

q4_condition_quadrat <- condition_data %>%
  filter(location == "Koh Tao", !is.na(genus)) %>%
  count(site, genus, reef_type, bleaching_period, quadrat_id, health, name = "points") %>%
  group_by(site, genus, reef_type, bleaching_period, quadrat_id) %>%
  complete(health = condition_codes, fill = list(points = 0)) %>%
  mutate(
    coral_points = sum(points),
    proportion = points / coral_points
  ) %>%
  ungroup() %>%
  select(-points) %>%
  pivot_wider(names_from = health, values_from = proportion, values_fill = 0)

#### Genus support ####

# The interaction requires both reef types in at least two bleaching periods.
q4_support <- q4_condition_quadrat %>%
  group_by(genus, bleaching_period, reef_type) %>%
  summarise(
    quadrats = n(),
    min_coral_points = min(coral_points),
    median_coral_points = median(coral_points),
    max_coral_points = max(coral_points),
    .groups = "drop"
  ) %>%
  group_by(genus, bleaching_period) %>%
  mutate(both_reef_types = n_distinct(reef_type) == 2) %>%
  ungroup()

q4_genera <- q4_support %>%
  distinct(genus, bleaching_period, both_reef_types) %>%
  group_by(genus) %>%
  summarise(periods_both_reef_types = sum(both_reef_types), .groups = "drop") %>%
  filter(periods_both_reef_types >= 2) %>%
  pull(genus)

#### Five-state trajectories ####

# Quadrat means describe population-level composition, not tracked colony transitions.
q4_trajectories <- q4_condition_quadrat %>%
  filter(genus %in% q4_genera) %>%
  group_by(genus, reef_type, bleaching_period) %>%
  summarise(
    quadrats = n(),
    across(all_of(condition_cols), mean),
    .groups = "drop"
  )

#### PERMANOVA ####

# Restrict each genus to periods represented by both reef types.
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
    
    condition_matrix <- genus_data %>% select(all_of(condition_cols))
    
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
    filter(term %in% c("reef_type", "bleaching_period", "reef_type:bleaching_period")) %>%
    transmute(
      genus = .y,
      term,
      df = Df,
      r2 = R2,
      pseudo_f = F,
      p = `Pr(>F)`
    )
)

# The prespecified Q4 inferential family is the AR × bleaching-period interaction across genera.
q4_interaction_results <- q4_permanova_results %>%
  filter(term == "reef_type:bleaching_period") %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(p_adj)

#### Dispersion diagnostics ####

q4_dispersion <- set_names(q4_genera) %>%
  map(function(g) {
    
    genus_data <- q4_permanova_data %>%
      filter(genus == g) %>%
      droplevels()
    
    group <- interaction(
      genus_data$reef_type,
      genus_data$bleaching_period,
      drop = TRUE
    )
    
    dispersion <- genus_data %>%
      select(all_of(condition_cols)) %>%
      vegan::vegdist(method = "bray") %>%
      vegan::betadisper(group, add = "lingoes")
    
    list(
      dispersion = dispersion,
      test = vegan::permutest(dispersion, permutations = n_perm)
    )
  })

q4_dispersion_results <- imap_dfr(
  q4_dispersion,
  ~ tibble(
    genus = .y,
    dispersion_f = .x$test$tab[1, "F"],
    dispersion_p = .x$test$tab[1, "Pr(>F)"]
  )
) %>%
  mutate(dispersion_p_adj = p.adjust(dispersion_p, method = "BH"))

saveRDS(q4_permanova, file.path(fits_dir, "q4_permanova_models.rds"))
saveRDS(q4_dispersion, file.path(fits_dir, "q4_dispersion_models.rds"))
write_csv(q4_support, file.path(stats_dir, "q4_genus_support.csv"))
write_csv(q4_trajectories, file.path(stats_dir, "q4_condition_trajectories.csv"))
write_csv(q4_permanova_results, file.path(stats_dir, "q4_permanova_all_terms.csv"))
write_csv(q4_interaction_results, file.path(stats_dir, "q4_interaction_results.csv"))
write_csv(q4_dispersion_results, file.path(stats_dir, "q4_dispersion_results.csv"))

q4_interaction_results
q4_dispersion_results

print("analysis done! fits and stats saved to the dated analysis folders")

summary(q1_model)
summary(q2_model)

q1_diagnostics
q2_diagnostics
