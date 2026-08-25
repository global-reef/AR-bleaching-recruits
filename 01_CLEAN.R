### 01. DATA CLEANING ####

#### Load data ####

cpce_clean <- read_csv(file.path(data_clean_dir, "cpce_long.csv"), show_col_types = FALSE)

#### Cleaning definitions ####

# TBA: confirm these codes against the original CPCE datasheets / field protocol.
# They are excluded temporarily rather than assigned an uncertain meaning.
substrate_tba <- c("COEL", "AV", "BA", "MM", "OAB")
health_tba <- c("N/S", "ORK", "POR", "PVL")

condition_codes <- c("H", "PBL", "FBL", "PRK", "FRK")

#### Standardise substrate and health codes ####

cpce_clean <- cpce_clean %>%
  mutate(
    substrate = str_to_upper(str_trim(substrate)),
    substrate = recode(
      substrate,
      "ACRP" = "ACRO",
      "FAV" = "FAVI",
      "PORI" = "POR",
      "PALTY" = "PLATY",
      "UNKWN" = "UNKN",
      "MSIC" = "MISC",
      "DIPS" = "DIPSA",
      "LETA" = "LEPTA",
      "FUNGI" = "FUNGII",
      "FUNGI?" = "FUNGII",
      "LITH0" = "LITH",
      "NA" = NA_character_,
      "-" = NA_character_
    ),
    health = str_to_upper(str_trim(health)),
    health = recode(
      health,
      "FB;" = "FBL",
      "PVL" = "PBL",
      "N/A" = NA_character_,
      "NA" = NA_character_,
      "HEALTH" = NA_character_,
      "-" = NA_character_
    )
  )

#### Audit and exclude unresolved codes ####

# Keep the unresolved codes visible in the workflow so their temporary exclusion
# is revisited once their meanings are confirmed.
cpce_clean %>%
  filter(substrate %in% substrate_tba | health %in% health_tba) %>%
  count(substrate, health, sort = TRUE) %>%
  print(n = Inf)

cpce_clean <- cpce_clean %>%
  filter(
    !substrate %in% substrate_tba,
    !health %in% health_tba
  )

#### Coral taxonomy ####

genus_lookup <- c(
  ACRO = "Acropora",
  ASTREO = "Astreopora",
  COELA = "Coelastrea",
  CTEN = "Ctenactis",
  CYCLO = "Cycloseris",
  CYPH = "Cyphastrea",
  DANA = "Danafungia",
  DIPLO = "Diploastrea",
  DIPSA = "Dipsastrea",
  DUNC = "Duncanopsammia",
  ECHINOPO = "Echinopora",
  FAVI = "Favites",
  FUNG = "Fungia",
  GALA = "Galaxea",
  GARD = "Gardineroseris",
  GONIA = "Goniastrea",
  GONIO = "Goniopora",
  HYD = "Hydnophora",
  LEPTA = "Leptastrea",
  LEPTOS = "Leptoseris",
  LITH = "Lithophyllon",
  LOBO = "Lobophyllia",
  MERU = "Merulina",
  MONTI = "Montipora",
  OULA = "Oulastrea",
  OULO = "Oulophyllia",
  PACHY = "Pachyseris",
  PARAGONI = "Paragoniastrea",
  PAVO = "Pavona",
  PECT = "Pectinia",
  PLATY = "Platygyra",
  PLESI = "Plesiastrea",
  PLEUR = "Pleuractis",
  POCIL = "Pocillopora",
  PODA = "Podabacia",
  POLY = "Polyphyllia",
  POR = "Porites",
  PSAM = "Psammocora",
  TURB = "Turbinaria"
)

coral_codes <- c(names(genus_lookup), "FUNGII")

cpce_clean <- cpce_clean %>%
  mutate(
    is_coral = substrate %in% coral_codes,
    genus = unname(genus_lookup[substrate]),
    taxon = case_when(
      !is.na(genus) ~ genus,
      substrate == "FUNGII" ~ "Fungiidae",
      TRUE ~ NA_character_
    )
  )

#### Bleaching period ####

# Periods represent equivalent stages of the 2024 bleaching event rather than identical calendar dates. 
# Rayong was sampled once at the Middle stage.
cpce_clean <- cpce_clean %>%
  mutate(
    bleaching_period = case_when(
      site == "Tanote" & date == as.Date("2024-04-21") ~ "Start",
      site == "Tanote" & date == as.Date("2024-05-30") ~ "Middle",
      site == "Tanote" & date == as.Date("2024-07-31") ~ "Late",
      site == "Tanote" & date == as.Date("2024-08-18") ~ "Post",
      
      site == "Mango" & date == as.Date("2024-05-09") ~ "Start",
      site == "Mango" & date %in% as.Date(c("2024-06-11", "2024-06-13")) ~ "Middle",
      site == "Mango" & date == as.Date("2024-08-03") ~ "Late",
      site == "Mango" & date == as.Date("2024-08-27") ~ "Post",
      
      location == "Rayong" & date == as.Date("2024-05-16") ~ "Middle",
      TRUE ~ NA_character_
    ),
    bleaching_period = factor(
      bleaching_period,
      levels = c("Start", "Middle", "Late", "Post")
    )
  )

#### Sampling identifiers ####

# Exact blocks, transects and quadrats were not repeatedly sampled through time.
# IDs therefore include date and identify independent sampling clusters within
# each survey event rather than longitudinally repeated physical units.
cpce_clean <- cpce_clean %>%
  mutate(
    sampling_unit = case_when(
      reef_type == "Artificial" ~ paste(location, site, date, "block", block, sep = "_"),
      reef_type == "Natural" ~ paste(location, site, date, "transect", transect, sep = "_")
    ),
    quadrat_id = paste(sampling_unit, "quadrat", quadrat, sep = "_"),
    life_stage = factor(
      if_else(reef_type == "Artificial", "Juvenile", "Mature"),
      levels = c("Mature", "Juvenile")
    )
  )

#### Coral-condition responses ####

# The primary response is nominal, not ordinal. PRK/FRK represent recent
# mortality and are not assumed to form a linear severity sequence with bleaching.
cpce_clean <- cpce_clean %>%
  mutate(
    condition = case_when(
      is_coral & health == "H" ~ "Healthy",
      is_coral & health == "PBL" ~ "Partially bleached",
      is_coral & health == "FBL" ~ "Fully bleached",
      is_coral & health == "PRK" ~ "Partially recently killed",
      is_coral & health == "FRK" ~ "Fully recently killed",
      TRUE ~ NA_character_
    ),
    condition = factor(
      condition,
      levels = c(
        "Healthy",
        "Partially bleached",
        "Fully bleached",
        "Partially recently killed",
        "Fully recently killed"
      )
    ),
    analysis_eligible = is_coral & health %in% condition_codes,
    
    # Locked binomial fallback A: healthy versus any bleaching/recent mortality.
    affected = case_when(
      is_coral & health == "H" ~ 0L,
      is_coral & health %in% c("PBL", "FBL", "PRK", "FRK") ~ 1L,
      TRUE ~ NA_integer_
    ),
    
    # Locked binomial fallback B: among affected coral, bleaching versus recent mortality.
    recent_mortality = case_when(
      is_coral & health %in% c("PBL", "FBL") ~ 0L,
      is_coral & health %in% c("PRK", "FRK") ~ 1L,
      TRUE ~ NA_integer_
    )
  )

#### Final checks  ####

# Every scheduled survey must map to one bleaching period.
stopifnot(!any(is.na(cpce_clean$bleaching_period)))

# Unresolved codes must not leak into the current analysis-ready dataset.
stopifnot(
  !any(cpce_clean$substrate %in% substrate_tba),
  !any(cpce_clean$health %in% health_tba)
)

# Primary-response eligibility must correspond exactly to the five retained states.
stopifnot(
  all(cpce_clean$analysis_eligible == (cpce_clean$is_coral & cpce_clean$health %in% condition_codes)),
  all(is.na(cpce_clean$condition) == !cpce_clean$analysis_eligible)
)

# Quadrat IDs must map to only one sampling hierarchy combination.
stopifnot(
  cpce_clean %>%
    distinct(location, site, reef_type, date, block, transect, quadrat, quadrat_id) %>%
    count(quadrat_id) %>%
    summarise(ok = all(n == 1)) %>%
    pull(ok)
)

cpce_clean %>%
  summarise(
    points = n(),
    coral_points = sum(is_coral),
    analysis_points = sum(analysis_eligible),
    genus_identified = sum(analysis_eligible & !is.na(genus)),
    family_only = sum(analysis_eligible & taxon == "Fungiidae", na.rm = TRUE),
    missing_health_coral = sum(is_coral & is.na(health))
  ) %>%
  print()

cpce_clean %>%
  filter(is_coral) %>%
  count(location, site, reef_type, bleaching_period, health, sort = TRUE) %>%
  print(n = Inf)

#### Save analysis-ready dataset ####

write_rds(cpce_clean, file.path(data_clean_dir, "cpce_clean.rds"))
