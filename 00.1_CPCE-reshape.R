### 00.1. CPCE RESHAPE ####

#### Load setup ####

source("00_SETUP.R")

#### File metadata ####

cpce_files <- tribble(
  ~file,          ~location, ~site,     ~reef_type,
  "KMN_AR.xlsx",  "Rayong",  "Rayong",  "Artificial",
  "KMN_NR.xlsx",  "Rayong",  "Rayong",  "Natural",
  "Mango_AR.xlsx","Koh Tao", "Mango",   "Artificial",
  "Mango_NR.xlsx","Koh Tao", "Mango",   "Natural",
  "Tanote_AR.xlsx","Koh Tao","Tanote",  "Artificial",
  "Tanote_NR.xlsx","Koh Tao","Tanote",  "Natural"
) %>%
  mutate(file = file.path(data_raw_dir, file))

#### File metadata ####
cpce_sheets <- tribble(
  ~file_name,       ~sheet,           ~date,
  "KMN_AR.xlsx",    "Dataset",        "2024-05-16",
  "KMN_NR.xlsx",    "16 May",         "2024-05-16",
  
  "Mango_AR.xlsx",  "09.05",          "2024-05-09",
  "Mango_AR.xlsx",  "13.06",          "2024-06-13",
  "Mango_AR.xlsx",  "03.08",          "2024-08-03",
  "Mango_AR.xlsx",  "27.08",          "2024-08-27",
  
  "Mango_NR.xlsx",  "Mango - 09.05",  "2024-05-09",
  "Mango_NR.xlsx",  "Mango - 11.06",  "2024-06-11",
  "Mango_NR.xlsx",  "Mango - 03.08",  "2024-08-03",
  "Mango_NR.xlsx",  "Mango - 27.08",  "2024-08-27",
  
  "Tanote_AR.xlsx", "21.04",           "2024-04-21",
  "Tanote_AR.xlsx", "30.05",           "2024-05-30",
  "Tanote_AR.xlsx", "31.07",           "2024-07-31",
  "Tanote_AR.xlsx", "18.08",           "2024-08-18",
  
  "Tanote_NR.xlsx", "21.04",           "2024-04-21",
  "Tanote_NR.xlsx", "30.05",           "2024-05-30",
  "Tanote_NR.xlsx", "31.07",           "2024-07-31",
  "Tanote_NR.xlsx", "18.08",           "2024-08-18"
) %>%
  mutate(date = as.Date(date))



#### Reshape CPCE files ####
# write one function for AR sheets and one for NR sheets, because their sampling hierarchy is different. Then bind them into one long dataset with common columns.

reshape_ar_sheet <- function(file, sheet, location, site, date) {
  
  raw <- read_excel(file, sheet = sheet, col_names = FALSE)
  
  raw %>%
    # Block labels occur in column 1 immediately before each set of 40 CPCE points.
    # cumsum() carries that block number down through its associated point rows.
    mutate(
      block = cumsum(coalesce(
        str_detect(as.character(...1), regex("^\\s*block", ignore_case = TRUE)),
        FALSE
      )),
      point = suppressWarnings(as.integer(...1))
    ) %>%
    
    # Retain only actual CPCE observations. This removes block/header/blank rows.
    filter(point %in% 1:40) %>%
    
    # AR surveys contain 10 quadrats per block, with paired substrate/health
    # columns for each quadrat. Only columns 2:21 therefore contain survey data.
    select(block, point, all_of(paste0("...", 2:21))) %>%
    pivot_longer(
      cols = starts_with("..."),
      names_to = "column",
      values_to = "value"
    ) %>%
    mutate(
      column = as.integer(str_remove(column, "^\\.\\.\\.")),
      quadrat = ((column - 2) %/% 2) + 1,
      field = if_else((column - 2) %% 2 == 0, "substrate", "health")
    ) %>%
    select(block, quadrat, point, field, value) %>%
    pivot_wider(names_from = field, values_from = value) %>%
    mutate(
      location = location,
      site = site,
      reef_type = "Artificial",
      date = as.Date(date),
      quadrat = as.integer(quadrat),
      health = na_if(health, "NA"),
      .before = 1
    )
}

# test one sheet 
test_ar <- reshape_ar_sheet(
  file = file.path(data_raw_dir, "Tanote_AR.xlsx"),
  sheet = "21.04",
  location = "Koh Tao",
  site = "Tanote",
  date = "2024-04-21"
)

glimpse(test_ar)
count(test_ar, block, quadrat)

# test across all AR 
ar_test <- cpce_sheets %>%
  filter(str_detect(file_name, "_AR\\.xlsx$")) %>%
  left_join(
    cpce_files %>% mutate(file_name = basename(file)),
    by = "file_name"
  ) %>%
  mutate(
    data = pmap(
      list(file, sheet, location, site, date),
      reshape_ar_sheet
    )
  )
# show whether any AR sheet has an unexpected number of blocks/quadrats before we bind everything together.
ar_test %>%
  transmute(file_name, sheet, rows = map_int(data, nrow))
# check 
ar_test %>%
  transmute(
    file_name,
    sheet,
    blocks = map_int(data, ~ n_distinct(.x$block)),
    quadrats = map_int(data, ~ n_distinct(paste(.x$block, .x$quadrat))),
    rows = map_int(data, nrow)
  )



reshape_nr_sheet <- function(file, sheet, location, site, date) {
  
  raw <- suppressMessages(read_excel(file, sheet = sheet, col_names = FALSE))
  
  # Locate rows containing quadrat headers (Q1-Q40). This avoids relying on
  # fixed column positions, which differ slightly among the source workbooks.
  header_rows <- which(
    apply(raw, 1, \(x) any(str_detect(as.character(x), "^Q[0-9]+$"), na.rm = TRUE))
  )
  
  map_dfr(header_rows, \(header_row) {
    
    quadrat_cols <- which(
      str_detect(as.character(raw[header_row, ]), "^Q[0-9]+$")
    )
    
    map_dfr(quadrat_cols, \(quadrat_col) {
      
      quadrat <- as.integer(
        str_remove(as.character(raw[[quadrat_col]][header_row]), "^Q")
      )
      
      tibble(
        transect = 1L,
        quadrat = quadrat,
        point = 1:40,
        substrate = as.character(raw[[quadrat_col]][(header_row + 2):(header_row + 41)]),
        health = as.character(raw[[quadrat_col + 1]][(header_row + 2):(header_row + 41)])
      )
    })
  }) %>%
    mutate(
      location = location,
      site = site,
      reef_type = "Natural",
      date = as.Date(date),
      health = na_if(health, "NA"),
      .before = 1
    ) %>%
    arrange(transect, quadrat, point)
}

# test tanote 
test_nr <- reshape_nr_sheet(
  file = file.path(data_raw_dir, "Tanote_NR.xlsx"),
  sheet = "30.05",
  location = "Koh Tao",
  site = "Tanote",
  date = "2024-05-30"
)

glimpse(test_nr)
count(test_nr, transect, quadrat)
nrow(test_nr)

# run across all nr 
nr_test <- cpce_sheets %>%
  filter(str_detect(file_name, "_NR\\.xlsx$")) %>%
  left_join(
    cpce_files %>% mutate(file_name = basename(file)),
    by = "file_name"
  ) %>%
  mutate(
    data = pmap(
      list(file, sheet, location, site, date),
      reshape_nr_sheet
    )
  )

# validate across all nr 
nr_test %>%
  transmute(
    file_name,
    sheet,
    quadrats = map_int(data, ~ n_distinct(.x$quadrat)),
    rows = map_int(data, nrow)
  )

# two only have 38 or 39 
nr_test %>%
  transmute(
    file_name,
    sheet,
    missing_quadrats = map_chr(
      data,
      ~ paste(setdiff(1:40, unique(.x$quadrat)), collapse = ", ")
    ),
    quadrats = map_int(data, ~ n_distinct(.x$quadrat)),
    rows = map_int(data, nrow)
  )

#### Bind reshaped data ####

cpce_long <- bind_rows(
  ar_test %>% select(data) %>% pull() %>% bind_rows(),
  nr_test %>% select(data) %>% pull() %>% bind_rows()
) %>%
  arrange(location, site, reef_type, date, block, transect, quadrat, point)

# check 
glimpse(cpce_long)

cpce_long %>%
  count(location, site, reef_type, date)

cpce_long %>%
  summarise(
    rows = n(),
    missing_substrate = sum(is.na(substrate)),
    missing_health = sum(is.na(health))
  )

#### Save long dataset ####

write_csv(cpce_long, file.path(data_clean_dir, "cpce_long.csv"))


# sampling effort 
cpce_long %>%
  group_by(location, site, reef_type, date) %>%
  summarise(
    blocks = n_distinct(block, na.rm = TRUE),
    transects = n_distinct(transect, na.rm = TRUE),
    quadrats = n_distinct(if_else(
      reef_type == "Artificial",
      paste(block, quadrat),
      paste(transect, quadrat)
    )),
    points = n(),
    .groups = "drop"
  )

