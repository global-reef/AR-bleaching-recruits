### 01. Q1 MATCHED LOCATION RESULTS ####

#### Model results ####

q1_model <- readRDS(file.path(fits_dir, "q1_multinomial_model.rds"))

q1_emm_latent <- emmeans(q1_model, ~ reef_type * location | condition, mode = "latent")
q1_interaction <- contrast(q1_emm_latent, interaction = c("pairwise", "pairwise"), by = "condition")

q1_emm <- emmeans(q1_model, ~ reef_type * location | condition, mode = "prob")

q1_posthoc <- contrast(q1_emm, method = "revpairwise", by = c("location", "condition"), adjust = "none") %>%
  as.data.frame() %>%
  as_tibble() %>%
  group_by(location) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

#### Detailed results ####

q1_interaction_table <- q1_interaction %>%
  as.data.frame() %>%
  as_tibble() %>%
  transmute(condition, estimate, SE, z_ratio = z.ratio, p = p.value)

q1_emm_table <- q1_emm %>%
  as.data.frame() %>%
  as_tibble() %>%
  transmute(location, reef_type, condition, probability = prob, SE, lower_CI = asymp.LCL, upper_CI = asymp.UCL)

q1_posthoc_table <- q1_posthoc %>%
  mutate(lower_CI = estimate - qnorm(0.975) * SE, upper_CI = estimate + qnorm(0.975) * SE) %>%
  transmute(location, condition, contrast, estimate, lower_CI, upper_CI, SE, z_ratio = z.ratio, p = p.value, p_adj)

write_csv(q1_interaction_table, file.path(stats_dir, "q1_interaction.csv"))
write_csv(q1_emm_table, file.path(stats_dir, "q1_estimated_probabilities.csv"))
write_csv(q1_posthoc_table, file.path(stats_dir, "q1_posthoc.csv"))

# post hoc between koh tao and rayong 

#### Q1 AR location contrasts ####

q1_ar_location <- emmeans(
  q1_model,
  ~ location | condition * reef_type,
  mode = "prob"
) %>%
  contrast(
    method = "revpairwise",
    by = c("condition", "reef_type")
  ) %>%
  summary(infer = TRUE, adjust = "none") %>%
  as.data.frame() %>%
  filter(reef_type == "Artificial") %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH")
  )

q1_ar_location
write_csv(q1_ar_location, file.path(stats_dir, "q1_ar_location.csv"))

### 02. Q2 KOH TAO TEMPORAL RESULTS ####

#### Model results ####

q2_model <- readRDS(file.path(fits_dir, "q2_multinomial_model.rds"))

q2_emm_latent <- emmeans(q2_model, ~ reef_type * bleaching_period * site | condition, mode = "latent")
q2_interaction <- contrast(q2_emm_latent, interaction = c("pairwise", "pairwise", "pairwise"), by = "condition")

q2_interaction_test <- test(q2_interaction, joint = TRUE, by = "condition") %>%
  as.data.frame() %>%
  as_tibble()

q2_emm <- emmeans(q2_model, ~ reef_type * bleaching_period * site | condition, mode = "prob")

q2_posthoc <- contrast(
  q2_emm, method = "revpairwise",
  by = c("site", "bleaching_period", "condition"),
  adjust = "none"
) %>%
  as.data.frame() %>%
  as_tibble() %>%
  group_by(site) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

#### Detailed results ####

q2_interaction_summary <- q2_interaction_test %>%
  transmute(condition, df = df1, chisq = Chisq, p = p.value)

q2_interaction_table <- q2_interaction %>%
  as.data.frame() %>%
  as_tibble() %>%
  transmute(
    condition,
    reef_type_contrast = reef_type_pairwise,
    period_contrast = bleaching_period_pairwise,
    site_contrast = site_pairwise,
    estimate, SE, z_ratio = z.ratio, p = p.value
  )

q2_emm_table <- q2_emm %>%
  as.data.frame() %>%
  as_tibble() %>%
  transmute(site, bleaching_period, reef_type, condition, probability = prob, SE, lower_CI = asymp.LCL, upper_CI = asymp.UCL)

q2_posthoc_table <- q2_posthoc %>%
  mutate(lower_CI = estimate - qnorm(0.975) * SE, upper_CI = estimate + qnorm(0.975) * SE) %>%
  transmute(site, bleaching_period, condition, contrast, estimate, lower_CI, upper_CI, SE, z_ratio = z.ratio, p = p.value, p_adj)

write_csv(q2_interaction_summary, file.path(stats_dir, "q2_interaction_summary.csv"))
write_csv(q2_interaction_table, file.path(stats_dir, "q2_interaction_components.csv"))
write_csv(q2_emm_table, file.path(stats_dir, "q2_estimated_probabilities.csv"))
write_csv(q2_posthoc_table, file.path(stats_dir, "q2_posthoc.csv"))


### 03. Q3 GENUS-SPECIFIC MATCHED-MIDDLE RESULTS ####

#### Results summary ####

q3_results <- read_csv(file.path(stats_dir, "q3_permanova_results.csv"), show_col_types = FALSE)

q3_summary <- q3_results %>%
  mutate(
    significant = p_adj < 0.05,
    dispersion_issue = dispersion_p_adj < 0.05,
    interpretation = case_when(
      significant & !dispersion_issue ~ "Significant",
      significant & dispersion_issue ~ "Significant; dispersion differs",
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(location, p_adj)

write_csv(q3_summary, file.path(stats_dir, "q3_results_summary.csv"))


### 04. Q4 GENUS-SPECIFIC TEMPORAL RESULTS ####

#### Results summary ####

q4_results <- read_csv(file.path(stats_dir, "q4_interaction_results.csv"), show_col_types = FALSE)
q4_dispersion <- read_csv(file.path(stats_dir, "q4_dispersion_results.csv"), show_col_types = FALSE)

q4_summary <- q4_results %>%
  left_join(q4_dispersion, by = "genus") %>%
  mutate(
    significant = p_adj < 0.05,
    dispersion_issue = dispersion_p_adj < 0.05,
    interpretation = case_when(
      significant & !dispersion_issue ~ "Significant",
      significant & dispersion_issue ~ "Significant; dispersion differs",
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(p_adj)

write_csv(q4_summary, file.path(stats_dir, "q4_results_summary.csv"))


### 05. PUBLICATION TABLES ####

#### T1: Genus-specific AR-NR condition differences ####

T1_genus_condition <- q3_summary %>%
  mutate(
    location = recode(location, "Koh Tao" = "kt", "Rayong" = "rayong"),
    dispersion = if_else(dispersion_issue, "Yes", "No")
  ) %>%
  select(genus, location, r2, p_adj, dispersion) %>%
  pivot_wider(
    names_from = location,
    values_from = c(r2, p_adj, dispersion),
    names_glue = "{location}_{.value}"
  ) %>%
  arrange(genus)

write_csv(T1_genus_condition, file.path(tables_dir, "T1_genus_condition_comparison.csv"))


#### S1: Q1-Q2 AR-NR post-hoc contrasts ####

q1_probabilities <- q1_emm_table %>%
  select(location, reef_type, condition, probability) %>%
  pivot_wider(names_from = reef_type, values_from = probability, names_glue = "{reef_type}_probability")

q2_probabilities <- q2_emm_table %>%
  select(site, bleaching_period, reef_type, condition, probability) %>%
  pivot_wider(names_from = reef_type, values_from = probability, names_glue = "{reef_type}_probability")

S1_posthoc <- bind_rows(
  q1_posthoc_table %>%
    left_join(q1_probabilities, by = c("location", "condition")) %>%
    transmute(
      question = "Q1", location, site = NA_character_, bleaching_period = "Middle", condition,
      natural_probability = Natural_probability,
      artificial_probability = Artificial_probability,
      ar_nr_difference = estimate, lower_CI, upper_CI, p_adj
    ),
  q2_posthoc_table %>%
    left_join(q2_probabilities, by = c("site", "bleaching_period", "condition")) %>%
    transmute(
      question = "Q2", location = "Koh Tao", site, bleaching_period, condition,
      natural_probability = Natural_probability,
      artificial_probability = Artificial_probability,
      ar_nr_difference = estimate, lower_CI, upper_CI, p_adj
    )
)

write_csv(S1_posthoc, file.path(tables_dir, "S1_primary_model_posthocs.csv"))


#### S2: Genus-specific temporal interactions ####

S2_genus_temporal <- q4_summary %>%
  transmute(
    genus, df, r2, pseudo_f, p_adj,
    dispersion_p_adj,
    dispersion_issue,
    interpretation
  )

write_csv(S2_genus_temporal, file.path(tables_dir, "S2_genus_temporal_interactions.csv"))


### 06. MANUSCRIPT OUTPUT ####

#### Results for text ####

cat("\nQ1: reef type × location interaction by condition\n")
print(q1_interaction_table, n = Inf)

cat("\nQ1: significant AR-NR post-hoc contrasts\n")
q1_posthoc_table %>% filter(p_adj < 0.05) %>% print(n = Inf)

cat("\nQ2: reef type × bleaching period × site interaction\n")
print(q2_interaction_summary, n = Inf)

cat("\nQ2: significant AR-NR post-hoc contrasts\n")
q2_posthoc_table %>% filter(p_adj < 0.05) %>% print(n = Inf)

cat("\nQ3: significant genus-specific AR-NR differences\n")
q3_summary %>% filter(significant) %>% print(n = Inf)

cat("\nQ4: genus-specific temporal interactions\n")
q4_summary %>% select(genus, r2, pseudo_f, p, p_adj, dispersion_p_adj) %>% print(n = Inf)

print("results done! detailed outputs saved to stats; publication tables saved to tables")