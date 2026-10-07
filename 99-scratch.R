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




###### testing plots for q3b #### 
### 05. F3b GENUS-SPECIFIC SITE CONDITION ####

#### Shared labels ####

reef_labels <- c(
  Natural = "Natural reef",
  Artificial = "Artificial reef"
)

#### Plot data ####

q3b_summary <- read_csv(file.path(stats_dir, "q3b_results_summary.csv"), show_col_types = FALSE)

q3b_plot_data <- read_csv(file.path(stats_dir, "q3b_condition_composition.csv"), show_col_types = FALSE) %>%
  left_join(
    q3b_summary %>%
      select(genus, reef_type, p_adj, significant, dispersion_issue),
    by = c("genus", "reef_type")
  ) %>%
  pivot_longer(
    cols = all_of(condition_order),
    names_to = "condition",
    values_to = "proportion"
  ) %>%
  mutate(
    site = factor(site, levels = c("Mango", "Tanote", "Rayong")),
    reef_type = factor(reef_type, levels = c("Natural", "Artificial")),
    condition = factor(condition, levels = condition_order),
    panel = case_when(
      significant & dispersion_issue ~ paste0(
        "<i>", genus, "</i><br>",
        reef_labels[as.character(reef_type)], " †"
      ),
      significant ~ paste0(
        "<i>", genus, "</i><br>",
        reef_labels[as.character(reef_type)], " *"
      ),
      TRUE ~ paste0(
        "<i>", genus, "</i><br>",
        reef_labels[as.character(reef_type)]
      )
    )
  )

#### Panel order ####

panel_order <- c(
  "<i>Platygyra</i><br>Artificial reef *",
  "<i>Porites</i><br>Artificial reef",
  "<i>Dipsastrea</i><br>Artificial reef †",
  "<i>Leptastrea</i><br>Artificial reef *",
  "<i>Diploastrea</i><br>Natural reef",
  "<i>Platygyra</i><br>Natural reef *",
  "<i>Porites</i><br>Natural reef *",
  "<i>Pavona</i><br>Natural reef *",
  "<i>Lobophyllia</i><br>Natural reef *",
  "<i>Favites</i><br>Natural reef"
)

q3b_plot_data <- q3b_plot_data %>%
  mutate(panel = factor(panel, levels = panel_order))

#### Horizontal stacked bars ####

F3b <- ggplot(q3b_plot_data, aes(x = site, y = proportion, fill = condition)) +
  geom_col(width = 0.72) +
  facet_wrap(~ panel, ncol = 5) +
  coord_flip() +
  scale_y_continuous(
    labels = scales::label_percent(),
    expand = c(0, 0),
    limits = c(0, 1)
  ) +
  scale_fill_manual(
    values = condition_palette,
    breaks = condition_order,
    labels = condition_labels
  ) +
  labs(
    x = NULL,
    y = "Condition composition",
    fill = "Condition"
  ) +
  theme_clean +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    legend.position = "bottom",
    strip.text = ggtext::element_markdown(),
    panel.spacing = unit(0.7, "lines"),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.4)
  )

F3b
save_figure(F3b, "F3b_genus_site_condition", figure_width, half_page_height)
F3b
