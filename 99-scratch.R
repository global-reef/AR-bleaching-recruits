## scratch 

### 04. GENUS-SPECIFIC NONPARAMETRIC TESTS ####

#### Quadrat-level susceptibility ####

genus_quadrat <- condition_data %>%
  filter(!is.na(genus)) %>%
  group_by(location, site, genus, reef_type, bleaching_period, quadrat_id) %>%
  summarise(
    coral_points = n(),
    affected_prop = mean(affected),
    .groups = "drop"
  )

# check 
#### Q3 permutation-test support ####

genus_quadrat %>%
  filter(bleaching_period == "Middle") %>%
  group_by(location, genus, reef_type) %>%
  summarise(
    quadrats = n(),
    unique_props = n_distinct(affected_prop),
    min_prop = min(affected_prop),
    max_prop = max(affected_prop),
    .groups = "drop"
  ) %>%
  arrange(location, genus, reef_type) %>%
  print(n = Inf)

#### Q3: AR vs NR within genus and location ####
#### Q3: AR vs NR within genus and location ####

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

### 05. Q4 GENUS-SPECIFIC CONDITION TRAJECTORIES ####

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

# Bray-Curtis compares the complete five-state condition composition.
# Permutations are restricted within site because Mango and Tanote are sampling
# blocks rather than the biological comparison of interest.
q4_permanova <- set_names(q4_genera) %>%
  map(function(g) {
    
    genus_data <- q4_condition_quadrat %>%
      filter(genus == g) %>%
      droplevels()
    
    condition_matrix <- genus_data %>%
      select(H, PBL, FBL, PRK, FRK)
    
    vegan::adonis2(
      condition_matrix ~ reef_type * bleaching_period,
      data = genus_data,
      permutations = n_perm,
      method = "bray",
      strata = genus_data$site,
      by = "terms"
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