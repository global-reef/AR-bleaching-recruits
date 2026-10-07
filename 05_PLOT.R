### 01. FIGURE SETUP ####

#### Shared settings ####
condition_order <- c("H", "PBL", "FBL", "PRK", "FRK")

condition_labels <- c(
  H = "Healthy",
  PBL = "Partially bleached",
  FBL = "Fully bleached",
  PRK = "Partially recently killed",
  FRK = "Fully recently killed"
)

condition_recode <- c(
  "Healthy" = "H",
  "Partially bleached" = "PBL",
  "Fully bleached" = "FBL",
  "Partially recently killed" = "PRK",
  "Fully recently killed" = "FRK"
)
reef_labels <- c(
  Natural = "Natural reef",
  Artificial = "Artificial reef"
)

theme_clean <- theme_minimal(base_family = "sans") # fix for pdfs 

save_figure <- function(plot, name, width, height) {
  ggsave(
    file.path(plots_dir, paste0(name, ".tiff")),
    plot,
    width = width, height = height, units = "in",
    dpi = 600, compression = "lzw"
  )
  
  ggsave(
    file.path(plots_dir, paste0(name, ".png")),
    plot,
    width = width, height = height, units = "in",
    dpi = 600
  )
}


figure_width <- 7.2
half_page_height <- 4.5
full_page_height <- 8

### 02. F1 MATCHED-MIDDLE CONDITION COMPOSITION ####

#### Plot data ####

q1_plot_data <- read_csv(file.path(stats_dir, "q1_estimated_probabilities.csv"), show_col_types = FALSE) %>%
  mutate(
    location = factor(location, levels = c("Koh Tao", "Rayong")),
    reef_type = factor(reef_type, levels = c("Natural", "Artificial")),
    condition = recode(as.character(condition), !!!condition_recode),
    condition = factor(condition, levels = condition_order)
  )
#### Predicted composition ####

F1 <- ggplot(q1_plot_data, aes(reef_type, probability, fill = condition)) +
  geom_col(width = 0.72) +
  facet_wrap(~ location, nrow = 1) +
  scale_x_discrete(labels = reef_labels) +
  scale_y_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
  scale_fill_manual(values = condition_palette, breaks = condition_order, labels = condition_labels) +
  labs(x = NULL, y = "Predicted condition probability", fill = "Condition") +
  theme_clean +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom")
F1



### 03. F2 KOH TAO TEMPORAL CONDITION TRAJECTORIES ####

#### Plot data ####
q2_plot_data <- read_csv(file.path(stats_dir, "q2_estimated_probabilities.csv"), show_col_types = FALSE) %>%
  mutate(
    site = factor(site, levels = c("Mango", "Tanote")),
    reef_type = factor(reef_type, levels = c("Natural", "Artificial")),
    bleaching_period = factor(bleaching_period, levels = c("Start", "Middle", "Late", "Post")),
    condition = recode(as.character(condition), !!!condition_recode),
    condition = factor(condition, levels = condition_order)
  )

#### Predicted trajectories ####
F2 <- ggplot(
  q2_plot_data,
  aes(
    x = bleaching_period,
    stratum = condition,
    alluvium = interaction(condition, site, reef_type),
    y = probability,
    fill = condition
  )
) +
  geom_flow(width = 0.16, alpha = 0.85) +
  geom_stratum(width = 0.16, colour = "white", linewidth = 0.3) +
  facet_grid(site ~ reef_type, labeller = labeller(reef_type = reef_labels)) +
  scale_y_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
  scale_fill_manual(values = condition_palette, breaks = condition_order, labels = condition_labels) +
  labs(x = "Bleaching period", y = "Predicted condition probability", fill = "Condition") +
  theme_clean +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(0.5, "lines")
  )
F2



### 04. F3 GENUS-SPECIFIC MATCHED-MIDDLE CONDITION ####

#### Plot data ####

q3_summary <- read_csv(file.path(stats_dir, "q3_results_summary.csv"), show_col_types = FALSE)

q3_plot_genera <- q3_summary %>%
  filter(location == "Rayong", significant) %>%
  arrange(p_adj) %>%
  pull(genus)

q3_plot_data <- condition_data %>%
  filter(location == "Rayong", bleaching_period == "Middle", genus %in% q3_plot_genera) %>%
  count(genus, reef_type, quadrat_id, health, name = "points") %>%
  group_by(genus, reef_type, quadrat_id) %>%
  complete(health = condition_order, fill = list(points = 0)) %>%
  mutate(proportion = points / sum(points)) %>%
  ungroup() %>%
  group_by(genus, reef_type, health) %>%
  summarise(proportion = mean(proportion), .groups = "drop") %>%
  mutate(
    genus = factor(genus, levels = q3_plot_genera),
    reef_type = factor(reef_type, levels = c("Natural", "Artificial")),
    health = factor(health, levels = condition_order)
  )

#### Condition composition ####

F3 <- ggplot(q3_plot_data, aes(reef_type, proportion, fill = health)) +
  geom_col(width = 0.72) +
  facet_wrap(~ genus, nrow = 1) +
  scale_x_discrete(labels = reef_labels) +
  scale_y_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
  scale_fill_manual(values = condition_palette, breaks = condition_order, labels = condition_labels) +
  labs(x = NULL, y = "Condition composition", fill = "Condition") +
  theme_clean +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "italic"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

F3

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

### 06. F4 GENUS-SPECIFIC TEMPORAL TRAJECTORIES ####

#### Plot data ####

plot_genera <- c(
  "Porites", "Platygyra", "Acropora",
  "Pocillopora", "Goniastrea", "Lobophyllia",
  "Montipora", "Dipsastrea", "Favites"
)

q4_raw <- read_csv(file.path(stats_dir, "q4_condition_trajectories.csv"), show_col_types = FALSE) %>%
  filter(genus %in% plot_genera) %>%
  mutate(
    genus = factor(genus, levels = plot_genera),
    reef_type = factor(reef_type, levels = c("Natural", "Artificial")),
    bleaching_period = factor(bleaching_period, levels = c("Start", "Middle", "Late", "Post"))
  )

q4_missing <- expand_grid(
  genus = factor(plot_genera, levels = plot_genera),
  reef_type = factor(c("Natural", "Artificial"), levels = c("Natural", "Artificial")),
  bleaching_period = factor(c("Start", "Middle", "Late", "Post"), levels = c("Start", "Middle", "Late", "Post"))
) %>%
  anti_join(q4_raw %>% distinct(genus, reef_type, bleaching_period), by = c("genus", "reef_type", "bleaching_period"))

q4_plot_data <- q4_raw %>%
  pivot_longer(cols = all_of(condition_order), names_to = "condition", values_to = "proportion") %>%
  complete(genus, reef_type, bleaching_period, condition = condition_order, fill = list(proportion = 0)) %>%
  mutate(condition = factor(condition, levels = condition_order))

#### Condition trajectories ####

F4 <- ggplot(
  q4_plot_data,
  aes(x = bleaching_period, stratum = condition, alluvium = condition, y = proportion, fill = condition)
) +
  ggalluvial::geom_flow(width = 0.15, alpha = 0.85) +
  ggalluvial::geom_stratum(width = 0.15, colour = "white", linewidth = 0.3) +
  geom_text(
    data = q4_missing,
    aes(x = bleaching_period, y = 0.5, label = "No data"),
    inherit.aes = FALSE, size = 2.4, colour = "grey45", angle = 90
  ) +
  facet_grid(genus ~ reef_type, labeller = labeller(reef_type = reef_labels)) +
  scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
  scale_fill_manual(values = condition_palette, breaks = condition_order, labels = condition_labels) +
  labs(x = "Bleaching period", y = "Condition composition", fill = "Condition") +
  theme_clean +
  theme(
    legend.position = "bottom",
    strip.text.y = element_text(face = "italic"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.y = grid::unit(0.4, "lines")
  )

save_figure(F1, "F1_matched_middle_condition", figure_width, half_page_height)
save_figure(F2, "F2_koh_tao_temporal_condition", figure_width, half_page_height)
save_figure(F3, "F3_rayong_genus_condition", figure_width, half_page_height)
save_figure(F4, "F4_genus_temporal_condition", figure_width, full_page_height)

print("plots done! publication figures saved to the dated plots folder")
F1
F2
F3
F4


