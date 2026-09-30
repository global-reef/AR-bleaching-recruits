### 01. Q4 CONDITION TRAJECTORIES ####


condition_palette <- c(
  H = "#2A2B59",
  PBL = "#3B4F7A",
  FBL = "#5E6FA3",
  PRK = "#8A8BC3",
  FRK = "#C5B4D9"
)

#### Condition palette ####

condition_palette <- c(
  H = "#E8D8B5",
  PBL = "#B9C7B5",
  FBL = "#7FA9A8",
  PRK = "#477F91",
  FRK = "#1F4E6B"
)
#### Plot data ####

q4_plot_data <- q4_trajectories %>%
  pivot_longer(
    cols = c(H, PBL, FBL, PRK, FRK),
    names_to = "condition",
    values_to = "proportion"
  ) %>%
  mutate(
    period_num = as.numeric(bleaching_period),
    condition = factor(condition, levels = c("FRK", "PRK", "FBL", "PBL", "H"))
  )


#### Condition rivers ####

ggplot(
  q4_plot_data,
  aes(period_num, proportion, fill = condition)
) +
  geom_area(position = "stack") +
  facet_grid(genus ~ reef_type) +
  scale_x_continuous(
    breaks = 1:4,
    labels = c("Start", "Middle", "Late", "Post"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    labels = scales::label_percent(),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = condition_palette,
    breaks = c("H", "PBL", "FBL", "PRK", "FRK"),
    labels = c(
      H = "Healthy",
      PBL = "Partially bleached",
      FBL = "Fully bleached",
      PRK = "Partially recently killed",
      FRK = "Fully recently killed"
    )
  ) +
  labs(
    x = "Bleaching period",
    y = "Condition composition",
    fill = "Condition"
  ) +
  theme_clean +
  theme(
    strip.text.y = element_text(face = "italic"),
    panel.spacing.y = unit(0.4, "lines")
  )

### alluvials 
### 01. Q4 CONDITION TRAJECTORIES ####

#### Alluvial trajectory ####

plot_genus <- "Porites"

q4_alluvial <- q4_trajectories %>%
  filter(genus == plot_genus) %>%
  pivot_longer(
    cols = c(H, PBL, FBL, PRK, FRK),
    names_to = "condition",
    values_to = "proportion"
  ) %>%
  complete(
    reef_type,
    bleaching_period,
    condition = c("H", "PBL", "FBL", "PRK", "FRK"),
    fill = list(proportion = 0)
  ) %>%
  mutate(
    bleaching_period = factor(bleaching_period, levels = c("Start", "Middle", "Late", "Post")),
    condition = factor(condition, levels = c("H", "PBL", "FBL", "PRK", "FRK"))
  )

ggplot(
  q4_alluvial,
  aes(
    x = bleaching_period,
    stratum = condition,
    alluvium = condition,
    y = proportion,
    fill = condition
  )
) +
  geom_flow(width = 0.18, alpha = 0.8) +
  geom_stratum(width = 0.18, colour = "white") +
  facet_wrap(~ reef_type, nrow = 1) +
  scale_fill_manual(
    values = condition_palette,
    labels = c(
      H = "Healthy",
      PBL = "Partially bleached",
      FBL = "Fully bleached",
      PRK = "Partially recently killed",
      FRK = "Fully recently killed"
    )
  ) +
  scale_y_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
  labs(
    x = "Bleaching period",
    y = "Mean quadrat-level condition proportion",
    fill = "Condition",
    title = plot_genus
  ) +
  theme_clean


# more than just porites 
### 01. Q4 CONDITION TRAJECTORIES ####

#### Focal genera ####

plot_genera <- c(
  "Porites", "Platygyra", "Acropora",
  "Pocillopora", "Goniastrea", "Lobophyllia",
  "Montipora", "Dipsastrea", "Favites"
)

q4_alluvial <- q4_trajectories %>%
  filter(genus %in% plot_genera) %>%
  pivot_longer(
    cols = c(H, PBL, FBL, PRK, FRK),
    names_to = "condition",
    values_to = "proportion"
  ) %>%
  complete(
    genus,
    reef_type,
    bleaching_period,
    condition = c("H", "PBL", "FBL", "PRK", "FRK"),
    fill = list(proportion = 0)
  ) %>%
  mutate(
    genus = factor(genus, levels = plot_genera),
    reef_type = factor(reef_type, levels = c("Natural", "Artificial")),
    bleaching_period = factor(bleaching_period, levels = c("Start", "Middle", "Late", "Post")),
    condition = factor(condition, levels = c("H", "PBL", "FBL", "PRK", "FRK"))
  )
# sampe size labels 
#### Sample-size labels ####

q4_missing <- q4_n_labels %>% filter(quadrats == 0)

ggplot(
  q4_alluvial,
  aes(
    x = bleaching_period,
    stratum = condition,
    alluvium = condition,
    y = proportion,
    fill = condition
  )
) +
  geom_flow(width = 0.15, alpha = 0.85) +
  geom_stratum(width = 0.15, colour = "white", linewidth = 0.3) +
  geom_text(data = q4_missing, aes(x = bleaching_period, y = 0.5, label = "n = 0"), 
            inherit.aes = FALSE, size = 2.7, colour = "grey45", angle = 90) + 
  facet_grid(genus ~ reef_type) +
  scale_fill_manual(
    values = condition_palette,
    breaks = c("H", "PBL", "FBL", "PRK", "FRK"),
    labels = c(
      H = "Healthy",
      PBL = "Partially bleached",
      FBL = "Fully bleached",
      PRK = "Partially recently killed",
      FRK = "Fully recently killed"
    )
  ) +
  scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
  labs(
    x = "Bleaching period",
    y = "Condition",
    fill = "Condition"
  ) +
  theme_clean +
  theme(
    legend.position = "bottom",
    strip.text.y = element_text(face = "italic"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
