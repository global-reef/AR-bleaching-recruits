### 02. EXPLORATORY DATA ANALYSIS ####

#### Analysis dataset ####

# Restrict the primary EDA to coral points that can enter the five-state condition model.
condition_data <- cpce_clean %>%
  filter(analysis_eligible)

#### Overall dataset structure ####

# Establish the scale of the analysis before examining biological patterns.
condition_data %>%
  summarise(
    coral_points = n(),
    quadrats = n_distinct(quadrat_id),
    sampling_units = n_distinct(sampling_unit),
    sites = n_distinct(site),
    genera = n_distinct(genus, na.rm = TRUE)
  )

condition_data %>%
  count(location, site, reef_type, bleaching_period, name = "coral_points") %>%
  arrange(location, site, reef_type, bleaching_period)

#### Sampling balance ####

# CPCE points are subsamples within quadrats, not independent biological replicates.
# Quadrat and sampling-unit replication therefore matter more than raw point counts.
condition_data %>%
  group_by(location, site, reef_type, bleaching_period) %>%
  summarise(
    sampling_units = n_distinct(sampling_unit),
    quadrats = n_distinct(quadrat_id),
    coral_points = n(),
    .groups = "drop"
  )

condition_data %>%
  count(location, site, reef_type, bleaching_period, quadrat_id, name = "coral_points") %>%
  group_by(location, site, reef_type, bleaching_period) %>%
  summarise(
    quadrats = n(),
    median_coral_points = median(coral_points),
    q1 = quantile(coral_points, 0.25),
    q3 = quantile(coral_points, 0.75),
    min_coral_points = min(coral_points),
    max_coral_points = max(coral_points),
    .groups = "drop"
  )

#### Sampling hierarchy ####

# Points are nested within quadrats, and quadrats are nested within blocks/transects.
# Quantify both levels before choosing the multinomial random-effects structure.
sampling_hierarchy <- condition_data %>%
  distinct(location, site, reef_type, bleaching_period, sampling_unit, quadrat_id) %>%
  count(location, site, reef_type, bleaching_period, sampling_unit, name = "quadrats_unit") %>%
  group_by(location, site, reef_type, bleaching_period) %>%
  summarise(
    sampling_units = n(),
    quadrats = sum(quadrats_unit),
    median_quadrats_unit = median(quadrats_unit),
    min_quadrats_unit = min(quadrats_unit),
    max_quadrats_unit = max(quadrats_unit),
    .groups = "drop"
  )

sampling_hierarchy

# Each quadrat must belong to exactly one sampling unit.
stopifnot(
  condition_data %>%
    distinct(quadrat_id, sampling_unit) %>%
    count(quadrat_id) %>%
    summarise(ok = all(n == 1)) %>%
    pull(ok)
)

#### Missingness and exclusions ####

# Document how many coral observations are excluded from the primary five-state response.
cpce_clean %>%
  filter(is_coral) %>%
  mutate(
    response_status = case_when(
      analysis_eligible ~ "Five-state response",
      is.na(health) ~ "Missing health",
      health == "DIS" ~ "Diseased",
      health == "OG" ~ "Overgrown",
      health == "UNKN" ~ "Unknown health",
      TRUE ~ "Other excluded"
    )
  ) %>%
  count(response_status, sort = TRUE) %>%
  mutate(percent = 100 * n / sum(n))

cpce_clean %>%
  filter(is_coral, !analysis_eligible) %>%
  count(location, site, reef_type, bleaching_period, health, sort = TRUE)

#### Five-state response distribution ####

# These are descriptive raw point proportions only. They do not account for
# clustering, unequal coral cover or differences in genus composition.
condition_data %>%
  count(condition, name = "n") %>%
  mutate(percent = 100 * n / sum(n))

condition_data %>%
  count(location, reef_type, condition, name = "n") %>%
  group_by(location, reef_type) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

condition_data %>%
  count(site, reef_type, bleaching_period, condition, name = "n") %>%
  group_by(site, reef_type, bleaching_period) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

#### Quadrat-level condition structure ####

# Check whether quadrats contain mixtures of condition states or are commonly
# dominated by a single state, which can influence multinomial model stability.
condition_data %>%
  count(location, site, reef_type, bleaching_period, quadrat_id, condition, name = "n") %>%
  group_by(location, site, reef_type, bleaching_period, quadrat_id) %>%
  summarise(
    coral_points = sum(n),
    states_present = n_distinct(condition),
    dominant_state_share = max(n) / sum(n),
    .groups = "drop"
  ) %>%
  group_by(location, site, reef_type, bleaching_period) %>%
  summarise(
    quadrats = n(),
    median_states = median(states_present),
    one_state_quadrats = sum(states_present == 1),
    pct_one_state = 100 * mean(states_present == 1),
    median_dominant_share = median(dominant_state_share),
    .groups = "drop"
  )

#### Temporally matched location comparison ####

# Rayong was sampled only at Middle, so all Koh Tao-Rayong comparisons must use
# Koh Tao Middle observations rather than pooling the full Koh Tao time series.
matched_location_data <- condition_data %>%
  filter(bleaching_period == "Middle")

matched_location_data %>%
  group_by(location, reef_type) %>%
  summarise(
    coral_points = n(),
    quadrats = n_distinct(quadrat_id),
    sampling_units = n_distinct(sampling_unit),
    .groups = "drop"
  )

matched_location_data %>%
  count(location, reef_type, condition, name = "n") %>%
  group_by(location, reef_type) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

#### Koh Tao temporal balance ####

kt_data <- condition_data %>%
  filter(location == "Koh Tao")

# Confirm representation of both reef types across all four event stages at both sites.
kt_data %>%
  group_by(site, reef_type, bleaching_period) %>%
  summarise(
    coral_points = n(),
    quadrats = n_distinct(quadrat_id),
    sampling_units = n_distinct(sampling_unit),
    .groups = "drop"
  )

# Raw temporal condition proportions are descriptive, not model estimates.
kt_data %>%
  count(site, reef_type, bleaching_period, condition, name = "n") %>%
  group_by(site, reef_type, bleaching_period) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

#### Genus representation ####

# Total points can overstate replication because multiple points may occur within
# the same quadrat. Assess both point counts and independent quadrat occurrence.
condition_data %>%
  filter(!is.na(genus)) %>%
  group_by(genus) %>%
  summarise(
    coral_points = n(),
    quadrats = n_distinct(quadrat_id),
    locations = n_distinct(location),
    reef_types = n_distinct(reef_type),
    periods = n_distinct(bleaching_period),
    .groups = "drop"
  ) %>%
  arrange(desc(quadrats), desc(coral_points)) %>%
  print(n = Inf)

#### Genus replication ####

# Build genus × reef type × period replication within each real site.
# Zeros are retained because absence from one side of an interaction affects estimability.
genus_replication <- condition_data %>%
  filter(!is.na(genus)) %>%
  group_by(location, site, genus, reef_type, bleaching_period) %>%
  summarise(
    coral_points = n(),
    quadrats = n_distinct(quadrat_id),
    .groups = "drop"
  ) %>%
  group_by(location, site, genus) %>%
  complete(
    reef_type = c("Artificial", "Natural"),
    bleaching_period = levels(condition_data$bleaching_period),
    fill = list(coral_points = 0, quadrats = 0)
  ) %>%
  ungroup()

#### Genus AR-NR overlap ####

# Within-genus AR-NR contrasts require the genus to occur in both reef types.
# Genera restricted to one reef type primarily describe community composition.
genus_overlap <- condition_data %>%
  filter(!is.na(genus)) %>%
  group_by(location, genus) %>%
  summarise(
    artificial_points = sum(reef_type == "Artificial"),
    natural_points = sum(reef_type == "Natural"),
    artificial_quadrats = n_distinct(quadrat_id[reef_type == "Artificial"]),
    natural_quadrats = n_distinct(quadrat_id[reef_type == "Natural"]),
    both_reef_types = artificial_points > 0 & natural_points > 0,
    .groups = "drop"
  ) %>%
  arrange(location, desc(both_reef_types), desc(artificial_quadrats + natural_quadrats))

genus_overlap %>%
  print(n = Inf)

#### Genus condition-state sparsity ####

# Assess representation of the five condition states within observed
# genus × reef type × period combinations. Whole-cell absences are already
# documented separately in genus_replication.

genus_state_summary <- condition_data %>%
  filter(!is.na(genus)) %>%
  count(location, site, genus, reef_type, bleaching_period, condition, name = "coral_points") %>%
  group_by(location, site, genus, reef_type, bleaching_period) %>%
  summarise(
    coral_points = sum(coral_points),
    states_present = n(),
    zero_states = 5L - states_present,
    min_nonzero_state_points = min(coral_points),
    .groups = "drop"
  )

#### Candidate genus temporal coverage ####

candidate_genera <- c(
  "Porites", "Platygyra", "Acropora", "Pocillopora", "Favites",
  "Goniastrea", "Dipsastrea", "Lobophyllia", "Montipora"
)

candidate_genus_replication <- genus_replication %>%
  filter(location == "Koh Tao", genus %in% candidate_genera) %>%
  group_by(genus, reef_type, bleaching_period) %>%
  summarise(
    coral_points = sum(coral_points),
    quadrats = sum(quadrats),
    sites_present = sum(coral_points > 0),
    .groups = "drop"
  )

candidate_genus_replication %>%
  arrange(genus, reef_type, bleaching_period) %>%
  print(n = Inf)

candidate_genus_states <- condition_data %>%
  filter(location == "Koh Tao", genus %in% candidate_genera) %>%
  count(genus, reef_type, bleaching_period, condition, name = "state_points") %>%
  group_by(genus, reef_type, bleaching_period) %>%
  summarise(
    coral_points = sum(state_points),
    states_present = n(),
    zero_states = 5L - states_present,
    min_nonzero_state_points = min(state_points),
    .groups = "drop"
  )

candidate_genus_states %>%
  arrange(genus, reef_type, bleaching_period) %>%
  print(n = Inf)

#### Genus binary-response support ####
# this will really tell us what the genus-based models can support 

genus_binary_support <- condition_data %>%
  filter(!is.na(genus)) %>%
  group_by(location, genus, reef_type, bleaching_period) %>%
  summarise(
    coral_points = n(),
    quadrats = n_distinct(quadrat_id),
    healthy = sum(affected == 0),
    affected = sum(affected == 1),
    .groups = "drop"
  )

genus_binary_support %>%
  arrange(location, genus, reef_type, bleaching_period) %>%
  print(n = Inf)

#### Potential confounding ####

# Reef type is structurally confounded with life stage:
# Artificial = juvenile recruits; Natural = mature reef corals.
stopifnot(
  all(condition_data$life_stage[condition_data$reef_type == "Artificial"] == "Juvenile"),
  all(condition_data$life_stage[condition_data$reef_type == "Natural"] == "Mature")
)

# Overall AR-NR differences may also partly reflect differing genus composition.
# Within-genus models address composition but cannot separate reef type from life stage.
condition_data %>%
  filter(!is.na(genus)) %>%
  count(location, reef_type, genus, name = "n") %>%
  group_by(location, reef_type) %>%
  mutate(proportion = n / sum(n)) %>%
  arrange(location, reef_type, desc(proportion)) %>%
  ungroup() %>%
  print(n = Inf)



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

#### EDA outputs ####

# Save only tables needed for model estimability decisions and audit trail.
write_csv(genus_replication, file.path(eda_dir, "genus_replication.csv"))
write_csv(genus_state_summary, file.path(eda_dir, "genus_state_replication.csv"))

sampling_hierarchy
print(candidate_genus_replication, n=Inf)
print(candidate_genus_states, n=Inf)



#### Total survey effort ####

cpce_long %>%
  mutate(
    sampling_unit = case_when(
      reef_type == "Artificial" ~ paste(location, site, date, "block", block, sep = "_"),
      reef_type == "Natural" ~ paste(location, site, date, "transect", transect, sep = "_")
    ),
    quadrat_id = paste(sampling_unit, "quadrat", quadrat, sep = "_")
  ) %>%
  summarise(
    cpce_points = n(),
    quadrats = n_distinct(quadrat_id),
    sampling_units = n_distinct(sampling_unit)
  )

#### Coral and analytical effort ####

cpce_clean %>%
  summarise(
    coral_points = sum(is_coral),
    coral_quadrats = n_distinct(quadrat_id[is_coral]),
    eligible_points = sum(analysis_eligible),
    eligible_quadrats = n_distinct(quadrat_id[analysis_eligible])
  )

cpce_long %>%
  mutate(
    sampling_unit = case_when(
      reef_type == "Artificial" ~ paste(location, site, date, "block", block, sep = "_"),
      reef_type == "Natural" ~ paste(location, site, date, "transect", transect, sep = "_")
    ),
    quadrat_id = paste(sampling_unit, "quadrat", quadrat, sep = "_")
  ) %>%
  group_by(reef_type) %>%
  summarise(
    cpce_points = n(),
    quadrats = n_distinct(quadrat_id),
    sampling_units = n_distinct(sampling_unit),
    .groups = "drop"
  )
